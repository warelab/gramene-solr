#!/usr/bin/env node
// genes.js — gene-level suggestion docs (ids, alt ids, free-text terms, gene trees, biotypes).
//
// usage: node genes.js <genes-core-url-IGNORED> [--terms <sorted terms file>]
//        SUGGEST_TERMS_TSV=<file>   same as --terms
//
// WHY THIS NO LONGER TALKS TO SOLR
// It used to recover the term/taxon tallies by faceting the genes core: a 588s
// _terms,taxon_id pivot plus one unbounded facet per taxon (1083s for 267), to get back a field
// mongo2solr.join.js had itself derived from mongo. 1,743s of round trip for data we already had,
// and it was the only part that ever failed -- on 2026-10-18 the gene_tree pivot's q-io promise
// never settled, the event loop drained, and the process exited 0 leaving a 4.49 GB genes.json
// with no closing ']' and no error, because q-io runs on Q and Q keeps rejections to itself.
//
// Now: mongo2solr.join.js writes the (term, taxon_id) pairs as it builds each solr doc
// (SOLR_TERMS_TSV), that file is sorted by term, and this script aggregates it in a STREAMING
// pass. Nothing large is held in memory for terms, there are no Solr queries at all, and the
// terms cannot drift from what the genes core actually indexes because one code path produces
// both. Cross-checked against the facets it replaces: 20,983,024 distinct terms and 2,526,566
// with count>=2, both exactly matching Solr's own numbers.
//
// The output docs, categories and schema are unchanged. The suggestions core and API are untouched.
'use strict';

var fs = require('fs');
var readline = require('readline');
var collections = require('gramene-mongodb-config');

var argv = process.argv.slice(2);
function opt(flag) { var i = argv.indexOf(flag); return i === -1 ? null : argv[i + 1]; }
var TERMS = opt('--terms') || process.env.SUGGEST_TERMS_TSV ||
            '/scratch/olson/release70_solr/terms.sorted.tsv';

// ---------------------------------------------------------------- buffered stdout
var OUT_MAX = 8 << 20, out_buf = [], out_len = 0;
function out(s) { out_buf.push(s); out_len += s.length; if (out_len >= OUT_MAX) { out_flush(); } }
function out_flush() { if (out_len) { fs.writeSync(1, out_buf.join('')); out_buf = []; out_len = 0; } }

var t_start = Date.now();
function progress(msg) { console.error('[' + Math.round((Date.now() - t_start) / 1000) + 's] ' + msg); }

// A truncated genes.json must never look like success: the stage only checks for a closing ']'.
var finished = false;
process.on('exit', function (code) {
  try { out_flush(); } catch (e) { /* nothing useful to do while exiting */ }
  if (!finished) {
    console.error('genes.js EXITED WITHOUT FINISHING (code ' + code + '): output is truncated, ' +
                  'no closing "]" was written.');
    if (code === 0) { process.exitCode = 1; }
  }
});

var n = 0;
function emit(doc) {
  out(n === 0 ? '[\n' : ',\n');
  n++;
  out(JSON.stringify(doc) + '\n');
}

(async function main() {
  if (!fs.existsSync(TERMS)) {
    throw new Error('terms file not found: ' + TERMS +
                    '\n  Produce it with: SOLR_TERMS_TSV=<out> node ../genes/mongo2solr.join.js >/dev/null' +
                    '\n  then: LC_ALL=C sort -t$\'\\t\' -k1,1 <out> > ' + TERMS);
  }

  // ---------------------------------------------------------------- 1. one mongo pass
  // Everything the remaining categories need: the id/alt-id suggestions, the biotype tallies and
  // the gene-tree tallies. All were previously separate Solr facet queries.
  var genesCol = await collections.genes.mongoCollection();
  var uniqueId = new Map();     // UPPER(_id) -> taxon_id
  var originalCase = new Map(); // UPPER(_id) -> original _id
  var altId = new Map();        // UPPER(alt) -> {display, taxa: Map(taxon -> count)}
  var biotype = new Map();      // biotype    -> Map(taxon -> count)
  var geneTree = new Map();     // tree id    -> Map(taxon -> count)
  var seen = 0;

  progress('streaming genes from mongo ...');
  var cursor = genesCol.find({}, {
    _id: 1, taxon_id: 1, db_type: 1, alt_id: 1, biotype: 1, 'homology.gene_tree.id': 1
  });
  for await (const d of cursor) {
    var UP = d._id.toUpperCase();
    uniqueId.set(UP, d.taxon_id);
    originalCase.set(UP, d._id);

    if (d.alt_id) {
      for (var ai = 0; ai < d.alt_id.length; ai++) {
        var a = d.alt_id[ai];
        if (typeof a !== 'string' || !a.length) { continue; }
        var A = a.toUpperCase();
        var rec = altId.get(A);
        if (!rec) { rec = { display: a, taxa: new Map() }; altId.set(A, rec); }
        rec.taxa.set(d.taxon_id, (rec.taxa.get(d.taxon_id) || 0) + 1);
      }
    }
    if (d.biotype) { bump(biotype, d.biotype, d.taxon_id); }
    var gt = d.homology && d.homology.gene_tree && d.homology.gene_tree.id;
    if (gt) { bump(geneTree, gt, d.taxon_id); }

    if (++seen % 4000000 === 0) {
      progress('  ' + seen + ' genes, heap ' + heapMB() + ' MB');
    }
  }
  progress('genes streamed: ' + seen + ' docs; ' + uniqueId.size + ' ids, ' + altId.size +
           ' alt ids, ' + biotype.size + ' biotypes, ' + geneTree.size + ' gene trees, heap ' +
           heapMB() + ' MB');
  collections.closeMongoDatabase();

  // ---------------------------------------------------------------- 2. terms, streamed
  // One group per distinct term, so nothing accumulates. Terms that ARE a gene id or an alt id are
  // skipped here exactly as before -- they get their own id/alt-id suggestions below.
  progress('aggregating ' + TERMS + ' ...');
  var rows = 0, terms_seen = 0, singles = 0, merged = 0, skipped_id = 0, skipped_nonword = 0, bad_decode = 0;
  var curKey = null;      // JSON-encoded term, used only for grouping/comparison
  var curTerm = null;     // the decoded term, which is what goes in the doc
  var curTaxa = null, curTotal = 0;
  var hasLetter = /[A-Za-z]/;

  function flushTerm() {
    if (curKey === null) { return; }
    if (curTerm === null) { bad_decode++; return; }   // unparseable row, counted not silently dropped
    terms_seen++;
    var UP = curTerm.toUpperCase();
    if (uniqueId.has(UP) || altId.has(UP)) { skipped_id++; return; }
    if (curTotal === 1) {
      // Previously emitted inside the per-taxon facet loop; same shape, same guard.
      if (!hasLetter.test(curTerm)) { skipped_nonword++; return; }
      var only = curTaxa.keys().next().value;
      emit({
        category: 'Genes: symbol/name', subcategory: 'term', id: '_term_' + (n + 1),
        display_name: curTerm, name: curTerm, fq_field: '_terms', fq_value: curTerm,
        taxon_id: [+only], taxon_freq: [1], num_genes: 1, relevance: 1
      });
      singles++;
    } else {
      var ids = [], counts = [];
      for (const [tax, cnt] of curTaxa) { ids.push(+tax); counts.push(cnt); }
      emit({
        category: 'Genes: symbol/name', subcategory: 'term', id: '_term_' + (n + 1),
        display_name: curTerm, name: curTerm, fq_field: '_terms', fq_value: curTerm,
        taxon_id: ids, taxon_freq: counts, num_genes: curTotal,
        relevance: curTotal < 10 ? 1.1 : 1
      });
      merged++;
    }
  }

  var rl = readline.createInterface({ input: fs.createReadStream(TERMS), crlfDelay: Infinity });
  for await (const line of rl) {
    if (!line) { continue; }
    var tab = line.lastIndexOf('\t');
    if (tab < 0) { continue; }
    var encoded = line.slice(0, tab);
    var taxon = +line.slice(tab + 1);
    if (encoded !== curKey) {
      flushTerm();
      curKey = encoded;
      curTaxa = new Map();
      curTotal = 0;
      // decode once per distinct term, not once per row
      try { curTerm = JSON.parse(encoded); } catch (e) { curTerm = null; }
    }
    curTaxa.set(taxon, (curTaxa.get(taxon) || 0) + 1);
    curTotal++;
    if (++rows % 10000000 === 0) {
      progress('  ' + rows + ' rows, ' + terms_seen + ' terms, ' + n + ' docs, heap ' + heapMB() + ' MB');
    }
  }
  flushTerm();
  progress('terms done: ' + rows + ' rows, ' + terms_seen + ' distinct terms -> ' +
           singles + ' singleton + ' + merged + ' merged docs (' + skipped_id +
           ' were ids/alt ids, ' + skipped_nonword + ' had no word character, ' +
           bad_decode + ' undecodable)');
  if (bad_decode) { throw new Error(bad_decode + ' term rows could not be JSON-decoded'); }

  // ---------------------------------------------------------------- 3. gene trees
  for (const [tree, taxa] of geneTree) {
    var gids = [], gcounts = [], gtot = 0;
    for (const [tax, cnt] of taxa) { gids.push(tax); gcounts.push(cnt); gtot += cnt; }
    emit({
      category: 'Gene trees', id: '_term_' + (n + 1), taxon_id: gids, taxon_freq: gcounts,
      display_name: tree, name: tree, fq_field: 'gene_tree', fq_value: tree,
      num_genes: gtot, relevance: 1
    });
  }
  progress('gene tree docs: ' + geneTree.size);

  // ---------------------------------------------------------------- 4. biotypes
  for (const [bt, taxa] of biotype) {
    var bids = [], bcounts = [], btot = 0;
    for (const [tax, cnt] of taxa) { bids.push(+tax); bcounts.push(cnt); btot += cnt; }
    emit({
      category: 'Biotypes', id: '_term_' + (n + 1), display_name: bt, name: bt,
      fq_field: 'biotype', fq_value: bt, taxon_id: bids, taxon_freq: bcounts,
      num_genes: btot, relevance: 1
    });
  }
  progress('biotype docs: ' + biotype.size);

  // ---------------------------------------------------------------- 5. primary ids
  progress('emitting ' + uniqueId.size + ' primary id docs');
  for (const [uid, uid_taxon] of uniqueId) {
    var orig = originalCase.get(uid);
    emit({
      category: 'Genes: primary id', subcategory: 'id', fq_field: 'id', fq_value: orig,
      id: orig, display_name: orig, num_genes: 1, relevance: 1.2,
      taxon_id: uid_taxon, taxon_freq: [1]
    });
  }

  // ---------------------------------------------------------------- 6. alternate ids
  progress('emitting ' + altId.size + ' alternate id docs');
  for (const [A, at] of altId) {
    var aids = [], acounts = [], atotal = 0;
    for (const [tid, cnt] of at.taxa) { aids.push(+tid); acounts.push(cnt); atotal += cnt; }
    emit({
      category: 'Genes: alternate id', subcategory: 'Alternate IDs', fq_field: 'alt_id',
      fq_value: at.display, id: 'alt_' + at.display, display_name: at.display, name: at.display,
      num_genes: atotal, relevance: 1.2, taxon_id: aids, taxon_freq: acounts
    });
  }

  out(']\n');
  out_flush();
  finished = true;
  progress('done: ' + n + ' suggestion docs');
})().catch(function (e) {
  out_flush();
  console.error('genes.js FAILED: ' + (e && e.stack || e));
  process.exit(1);
});

function heapMB() { return Math.round(process.memoryUsage().heapUsed / 1048576); }
function bump(map, key, taxon) {
  var m = map.get(key);
  if (!m) { m = new Map(); map.set(key, m); }
  m.set(taxon, (m.get(taxon) || 0) + 1);
}
