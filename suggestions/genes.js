#!/usr/bin/env node
var http = require("q-io/http")
  , request = require("sync-request")
  , collections = require('gramene-mongodb-config');

var genesURL = process.argv[2];

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
collections.genes.mongoCollection().then(function(collection) {
  collection.find({},{_id:1, taxon_id:1}).toArray(function(err,docs) {
    if (err) throw err;
    var uniqueId = {};
    var uniqueTaxa = {};
    docs.forEach(function(d) {
      uniqueTaxa[d.taxon_id] = 1;
      uniqueId[d._id.toUpperCase()] = d.taxon_id;
    });

    collections.genetrees.mongoCollection().then(function(collection) {
      collection.find({},{_id:1,taxon_id:1}).toArray(function(err,docs) {
        if (err) throw err;
        collections.closeMongoDatabase();
        var treeRootNodeTaxonId = {};
        docs.forEach(function(d) {
          treeRootNodeTaxonId[d._id] = d.taxon_id;
        });

        // non-unique terms get their own doc with a relevant fq field
        // build a lookup table with the non-unique ids from these
        var url = genesURL + '/query?rows=0&facet=true&facet.limit=-1&json.nl=map&facet.pivot.mincount=2&facet.pivot=_terms,taxon_id';
        console.error(url);
        http.read(url).then(function(data) {
          var term_freq = {};
          var taxa_lut = {};
          JSON.parse(data).facet_counts.facet_pivot['_terms,taxon_id'].forEach(function(d) {
            var term = d.value.toUpperCase();
            if (!uniqueId.hasOwnProperty(term)) {
              if (term_freq.hasOwnProperty(term)) {
                term_freq[term] += d.count;
              }
              else {
                term_freq[term] = d.count;
                taxa_lut[term] = {}; // term may appear only once in some species, but we don't know which ones yet
              }
              
              if (d.pivot) { // term appears at least 2 times in some species
                d.pivot.forEach(function(p) {
                  if (taxa_lut[term].hasOwnProperty(p.value)) {
                    taxa_lut[term][p.value] += p.count;
                  }
                  else {
                    taxa_lut[term][p.value] = p.count;
                  }
                });
              }
            }
          });

          // find all the terms in each species
          Object.keys(uniqueTaxa).forEach(function(taxon_id) {
            var url = genesURL + '/query?rows=0&facet=true&facet.field=_terms&facet.limit=-1&json.nl=map&facet.mincount=1&q=taxon_id:'+taxon_id;
            console.error(url);
            var res = request('GET',url);
            var term_tally = JSON.parse(res.getBody()).facet_counts.facet_fields._terms;
            for (var t in term_tally) {
              var term = t.toUpperCase();
              if (term_freq.hasOwnProperty(term) && !taxa_lut[term].hasOwnProperty(taxon_id)) { // this is a non-unique term that appears only once in this genome
                taxa_lut[term][taxon_id] = term_tally[t];
              }
            }
          });
          
          console.error("finished building taxa_lut");
          var n=0;
          for (var term in term_freq) {
            var tf = term_freq[term];
            var taxa = {ids:[],counts:[]};
            for(var taxon_id in taxa_lut[term]) {
              taxa.ids.push(+taxon_id);
              taxa.counts.push(taxa_lut[term][taxon_id]);
            }
            var solr = {
              category    : 'Gene',
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
            if (n===1) console.log('[');
            else console.log(',');
            console.log(JSON.stringify(solr));
          }
      
          // gene tree suggestion
          var url = genesURL + '/query?rows=0&facet=true&facet.limit=-1&json.nl=map&facet.mincount=1&facet.field=gene_tree';
          console.error(url);
          http.read(url).then(function(data) {
            var tree_size = JSON.parse(data).facet_counts.facet_fields.gene_tree;
            for (var tree in tree_size) {
              var solr = {
                category    : 'Gene tree',
                id          : '_term_'+ ++n,
                taxon_id    : [1, treeRootNodeTaxonId[tree]],
                taxon_freq  : [tree_size[tree],tree_size[tree]],
                display_name: tree,
                name        : tree,
                fq_field    : 'gene_tree',
                fq_value    : tree,
                num_genes   : tree_size[tree],
                relevance   : 1
              };
              console.log(',');
              console.log(JSON.stringify(solr));
            }

            // biotype suggestion
            var url = genesURL + '/query?rows=0&facet=true&facet.limit=-1&json.nl=map&facet.mincount=1&facet.field=biotype';
            console.error(url);
            http.read(url).then(function(data) {
              var biotype_count = JSON.parse(data).facet_counts.facet_fields.biotype;
              for (var biotype in biotype_count) {
                var taxa = getTaxa('biotype:'+biotype);
                var solr = {
                  category    : 'Biotype',
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
                console.log(',');
                console.log(JSON.stringify(solr));
              }
      
              // read the mongo genes docs and generate suggestions with unique values only
              require('readline').createInterface(
                {
                  input: process.stdin,
                  terminal: false
                }
              ).on('line', function(line) { // one JSON object per line
                var mongo = JSON.parse(line);

                console.log(',');
                console.log(JSON.stringify({
                  category : 'Gene',
                  fq_field : 'id',
                  fq_value : mongo._id,
                  id       : mongo._id,
                  display_name : mongo._id,
                  num_genes : 1,
                  relevance : 1.2,
                  taxon_id : [mongo.taxon_id],
                  taxon_freq : [1]
                }));

                // add uniquely identifying xrefs
                if (mongo.xrefs) {
                  var xref_h = {};
                  mongo.xrefs.forEach(function(xref) {
                    xref.ids.forEach(function(id) {
                      xref_h[id]=xref.db;
                    });
                  });
                  Object.keys(xref_h).filter(function(xr) {
                    return !term_freq[xr.toUpperCase()];
                  }).forEach(function(xr) {
                    console.log(',');
                    console.log(JSON.stringify({
                      category : 'Gene',
                      fq_field : 'id',
                      fq_value : mongo._id,
                      id       : '_term_'+ ++n,
                      xref     : xr,
                      display_name : xr,
                      num_genes : 1,
                      relevance : 1.2,
                      taxon_id : [mongo.taxon_id],
                      taxon_freq : [1]
                    }));
                  });
                }

                // add uniquely identifying synonyms
                if (mongo.hasOwnProperty('synonyms')) {
                  mongo.synonyms.filter(function(syn) {
                    return !term_freq[syn.toUpperCase()];
                  }).forEach(function(syn) {
                    console.log(',');
                    console.log(JSON.stringify({
                      category : 'Gene',
                      fq_field : 'id',
                      fq_value : mongo._id,
                      id       : '_term_'+ ++n,
                      synonym  : syn,
                      display_name : syn,
                      num_genes : 1,
                      relevance : 1.1,
                      taxon_id : [mongo.taxon_id],
                      taxon_freq : [1]
                    }));
                  });
                }

              }).on('close', function() {
                console.log(']');
              });
            });
          });
        });
      });
    });
  });
});

