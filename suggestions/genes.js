#!/usr/bin/env node
var http = require("q-io/http")
  , request = require("sync-request")
  , fs = require("fs")
  , collections = require('gramene-mongodb-config');

var genesURL = process.argv[2];

// ---------------------------------------------------------------- scaling notes
// This script stalled for 20 hours on release 70 without emitting a single byte, and would have
// OOMed before finishing. It is not slow in proportion to the data: genes70 has 3x the documents
// of the sorghum build but 9.8x the distinct _terms values (20,983,024 against 2,135,970), because
// 267 genomes contribute their own id namespaces plus 224,197 tree ids across five tree sets.
// sorghum's whole stage 65 took 6m50s.
//
// Four things made that fatal, all fixed below:
//   1. find().toArray() materialised all 16,127,051 gene docs at once. Now streamed.
//   2. uniqueId / originalCase / term_freq / taxa_lut were plain objects. At tens of millions of
//      keys V8 drops those into dictionary mode, where `for...in` also materialises an enumeration
//      cache of every key. They are Maps now, which stay fast at this size and iterate without
//      building a key array.
//   3. the per-taxon loop called term.toUpperCase() TWICE per term, allocating two strings for
//      every one of ~21M terms per taxon pass. Hoisted to one call.
//   4. every emitted doc cost two console.log() syscalls; at 40M+ docs that alone is ~80M writes.
//      Output is buffered and flushed in 8 MB blocks.
// Also removed two structures that were written on every one of the 16.1M docs and never read
// (otherfeaturesId, treeRootNodeTaxonId) -- the latter was the only reason to load the genetree
// collection at all, so that query is gone too.
//
// NB _term_N ids may differ from a previous run: Map iterates in insertion order where `for...in`
// reorders integer-like keys numerically. The ids are synthetic and the core is recreated fresh
// each run, so nothing downstream pins them.

var OUT_MAX = 8 << 20, out_buf = [], out_len = 0;
function out(s) {
  out_buf.push(s); out_len += s.length;
  if (out_len >= OUT_MAX) { out_flush(); }
}
function out_flush() {
  if (out_len) { fs.writeSync(1, out_buf.join('')); out_buf = []; out_len = 0; }
}
var t_start = Date.now();
function progress(msg) {
  console.error('[' + Math.round((Date.now() - t_start) / 1000) + 's] ' + msg);
}

function getTaxa(query) {
  var url = genesURL + '/query?q='+query+'&rows=0&facet=true&facet.field=taxon_id&facet.limit=-1&facet.mincount=1&json.nl=map';
  var res = request('GET',url);
  var taxa = { ids:[], counts:[]};
  var taxTally = JSON.parse(res.getBody()).facet_counts.facet_fields.taxon_id;
  for (var t in taxTally) {
    taxa.ids.push(+t);
    taxa.counts.push(taxTally[t]);
  }
  return taxa;
}

// read all of the unique ids from the genes collection
// to avoid suggesting them as non-unique terms when they are mentioned as an xref or something in another gene
collections.genes.mongoCollection().then(async function(collection) {
  try {
    // STREAMED, not toArray(): 16,127,051 docs will not fit in heap as live JS objects.
    var uniqueId = new Map();      // UPPER(_id) -> taxon_id
    var uniqueTaxa = new Map();    // taxon_id   -> 1
    var originalCase = new Map();  // UPPER(_id) -> original-case _id
    // alternate ids (Sobic./Sb.../GRMZM/pan.../Vitvi.../LOC_Os…) get their own
    // resolvable suggestions that query the genes-core alt_id field, so an alt id and
    // its primary id both land on the same gene. UPPERCASE key -> { display: original
    // case alt id, taxa: {taxon_id: gene count} }.
    var altId = new Map();         // UPPER(alt) -> {display, taxa: Map(taxon_id -> count)}
    progress('streaming genes ...');
    var nseen = 0;
    var cursor = collection.find({}, {_id:1, taxon_id:1, db_type:1, alt_id:1});
    for await (const d of cursor) {
      var UP = d._id.toUpperCase();
      uniqueTaxa.set(d.taxon_id, 1);
      uniqueId.set(UP, d.taxon_id);
      originalCase.set(UP, d._id);
      if (d.alt_id) {
        for (var ai = 0; ai < d.alt_id.length; ai++) {
          var a = d.alt_id[ai];
          if (!a) { continue; }
          var A = a.toUpperCase();
          var rec = altId.get(A);
          if (!rec) { rec = { display: a, taxa: new Map() }; altId.set(A, rec); }
          rec.taxa.set(d.taxon_id, (rec.taxa.get(d.taxon_id) || 0) + 1);
        }
      }
      if (++nseen % 2000000 === 0) {
        progress('  ' + nseen + ' genes; ' + uniqueId.size + ' ids, ' + altId.size +
                 ' alt ids, heap ' + Math.round(process.memoryUsage().heapUsed / 1048576) + ' MB');
      }
    }
    progress('genes streamed: ' + nseen + ' docs, ' + uniqueId.size + ' ids, ' + altId.size +
             ' alt ids, ' + uniqueTaxa.size + ' taxa, heap ' +
             Math.round(process.memoryUsage().heapUsed / 1048576) + ' MB');

    // The genetree collection used to be loaded here purely to build treeRootNodeTaxonId, which
    // nothing ever read. Dropped along with it.
    collections.closeMongoDatabase();
    {

        // non-unique terms get their own doc with a relevant fq field
        // build a lookup table with the non-unique ids from these
        var url = genesURL + '/query?rows=0&facet=true&facet.limit=-1&json.nl=map&facet.pivot.mincount=2&facet.pivot=_terms,taxon_id';
        console.error(url);
        http.read(url).then(function(data) {
          var term_freq = new Map();   // term -> total gene count
          var taxa_lut = new Map();    // term -> Map(taxon_id -> count)
          progress('parsing the _terms,taxon_id pivot ...');
          JSON.parse(data).facet_counts.facet_pivot['_terms,taxon_id'].forEach(function(d) {
            var term = d.value; // .toUpperCase(); // don't do this it leads to problems
            var UP = term.toUpperCase();   // was computed twice per term
            if (!uniqueId.has(UP) && !altId.has(UP)) {
              if (term_freq.has(term)) {
                term_freq.set(term, term_freq.get(term) + d.count);
              }
              else {
                term_freq.set(term, d.count);
                taxa_lut.set(term, new Map()); // term may appear only once in some species, but we don't know which ones yet
              }

              if (d.pivot) { // term appears at least 2 times in some species
                var tl = taxa_lut.get(term);
                d.pivot.forEach(function(p) {
                  tl.set(p.value, (tl.get(p.value) || 0) + p.count);
                });
              }
            }
          });
          progress('pivot parsed: ' + term_freq.size + ' non-unique terms, heap ' +
                   Math.round(process.memoryUsage().heapUsed / 1048576) + ' MB');
          var n=0;
          // Find all the terms in each species. THIS is where 20 hours went on release 70: one
          // pass per taxon over that taxon's entire _terms facet, and for every term it called
          // toUpperCase() twice and did two hasOwnProperty lookups against 16.1M-key dictionary
          // objects. Now one toUpperCase and Map.has(). Logged per taxon so a stall is visible.
          var taxaList = Array.from(uniqueTaxa.keys());
          progress('scanning _terms per taxon (' + taxaList.length + ' taxa)');
          taxaList.forEach(function(taxon_id, ti) {
            var url = genesURL + '/query?rows=0&facet=true&facet.field=_terms&facet.limit=-1&json.nl=map&facet.mincount=1&q=taxon_id:'+taxon_id;
            var res = request('GET',url);
            var term_tally = JSON.parse(res.getBody()).facet_counts.facet_fields._terms;
            var nt = 0;
            for (var t in term_tally) {
              var term = t; // .toUpperCase(); // don't do this it leads to problems
              var UP = term.toUpperCase();
              nt++;
              if (!uniqueId.has(UP) && !altId.has(UP)) {
                var tl = taxa_lut.get(term);
                if (tl) {
                  tl.set(taxon_id, (tl.get(taxon_id) || 0) + term_tally[t]);
                }
                else if (/[A-Za-z]/.test(term)) { // this is a unique term that contains a word character
                  var solr = {
                    category    : 'Genes: symbol/name',
                    subcategory : 'term',
                    id          : '_term_'+ ++n,
                    display_name: term,
                    name        : term,
                    fq_field    : '_terms',
                    fq_value    : term,
                    taxon_id    : [+taxon_id],
                    taxon_freq  : [1],
                    num_genes   : 1,
                    relevance   : 1
                  };
                  if (n===1) out('[\n');
                  else out(',\n');
                  out(JSON.stringify(solr) + '\n');
                }
              }
            }
            progress('  taxon ' + (ti + 1) + '/' + taxaList.length + ' (' + taxon_id + '): ' +
                     nt + ' terms, ' + n + ' docs so far, heap ' +
                     Math.round(process.memoryUsage().heapUsed / 1048576) + ' MB');
          });

          progress('finished building taxa_lut. n=' + n + '; emitting ' + term_freq.size +
                   ' non-unique term docs');
          for (const [term, tf] of term_freq) {
            var taxa = {ids:[],counts:[]};
            var tlm = taxa_lut.get(term);
            if (tlm) {
              for (const [taxon_id, cnt] of tlm) {
                taxa.ids.push(+taxon_id);
                taxa.counts.push(cnt);
              }
            }
            var solr = {
              category    : 'Genes: symbol/name',
              subcategory : 'term',
              id          : '_term_'+ ++n,
              display_name: term,
              name        : term,
              fq_field    : '_terms',
              fq_value    : term,
              taxon_id    : taxa.ids,
              taxon_freq  : taxa.counts,
              num_genes   : tf,
              relevance   : tf < 10 ? 1.1 : 1
            };
            if (n===1) out('[\n');
            else out(',\n');
            out(JSON.stringify(solr) + '\n');
          }
      
          // gene tree suggestion
          // var url = genesURL + '/query?rows=0&facet=true&facet.limit=-1&json.nl=map&facet.mincount=1&facet.field=gene_tree';
          var url = genesURL
          + '/query?rows=0&facet=true&facet.limit=-1&facet.mincount=1&json.nl=map'
          + '&facet.pivot=gene_tree,taxon_id';
          console.error(url);
          http.read(url).then(function(data) {
            JSON.parse(data).facet_counts.facet_pivot['gene_tree,taxon_id'].forEach(function(d) {
              var taxa = {ids:[],counts:[]};
              d.pivot.forEach(function(p) {
                taxa.ids.push(p.value);
                taxa.counts.push(p.count);
              });
              var solr = {
                category    : 'Gene trees',
                id          : '_term_'+ ++n,
                taxon_id    : taxa.ids,
                taxon_freq  : taxa.counts,
                display_name: d.value,
                name        : d.value,
                fq_field    : 'gene_tree',
                fq_value    : d.value,
                num_genes   : d.count,
                relevance   : 1
              }
              out(',\n');
              out(JSON.stringify(solr) + '\n');
            });

            // biotype suggestion
            var url = genesURL + '/query?rows=0&facet=true&facet.limit=-1&json.nl=map&facet.mincount=1&facet.field=biotype';
            console.error(url);
            http.read(url).then(function(data) {
              var biotype_count = JSON.parse(data).facet_counts.facet_fields.biotype;
              for (var biotype in biotype_count) {
                var taxa = getTaxa('biotype:'+biotype);
                var solr = {
                  category    : 'Biotypes',
                  id          : '_term_'+ ++n,
                  display_name: biotype,
                  name        : biotype,
                  fq_field    : 'biotype',
                  fq_value    : biotype,
                  taxon_id    : taxa.ids,
                  taxon_freq  : taxa.counts,
                  num_genes   : biotype_count[biotype],
                  relevance   : 1
                };
                out(',\n');
                out(JSON.stringify(solr) + '\n');
              }
      
              // output unique IDs
              progress('output uniqueIds (' + uniqueId.size + ')');
              for (const [uid, uid_taxon] of uniqueId) {
                out(',\n');
                out(JSON.stringify({
                  category : 'Genes: primary id',
                  subcategory : 'id',
                  fq_field : 'id',
                  fq_value : originalCase.get(uid),
                  id       : originalCase.get(uid),
                  display_name : originalCase.get(uid),
                  num_genes : 1,
                  relevance : 1.2,
                  taxon_id : uid_taxon,
                  taxon_freq : [1]
                }) + '\n');
              }

              // output alternate ids — each is matched on `name` and resolves by
              // querying the genes-core alt_id field, so an alt id and its primary id
              // land on the same gene. Shown under the "Alternate IDs" subcategory.
              progress('output alt_ids (' + altId.size + ')');
              for (const [A, at] of altId) {
                var ataxa = { ids: [], counts: [] };
                var atotal = 0;
                for (const [tid, cnt] of at.taxa) {
                  ataxa.ids.push(+tid);
                  ataxa.counts.push(cnt);
                  atotal += cnt;
                }
                out(',\n');
                out(JSON.stringify({
                  category : 'Genes: alternate id',
                  subcategory : 'Alternate IDs',
                  fq_field : 'alt_id',
                  fq_value : at.display,
                  id       : 'alt_' + at.display,
                  display_name : at.display,
                  name     : at.display,
                  num_genes : atotal,
                  relevance : 1.2,
                  taxon_id : ataxa.ids,
                  taxon_freq : ataxa.counts
                }) + '\n');
              }
              out(']\n');
              out_flush();
              progress('done: ' + n + ' suggestion docs written');
            });   // biotype
          });     // gene_tree pivot
        });       // _terms pivot
    }             // (block that replaced the removed genetree callbacks)
  } catch (e) {
    // The original swallowed nothing but also surfaced nothing: an error inside the streaming
    // loop would have rejected silently and left a 0-byte genes.json, which is exactly the state
    // the 20-hour run was mistaken for. Fail loudly instead.
    out_flush();
    console.error('genes.js FAILED: ' + (e && e.stack || e));
    process.exit(1);
  }
});

