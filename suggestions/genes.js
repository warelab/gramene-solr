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
  collection.find({},{_id:1, taxon_id:1, db_type:1, synonyms:1}).toArray(function(err,docs) {
    if (err) throw err;
    var uniqueId = {};
    var otherfeaturesId = {};
    var uniqueTaxa = {};
    var originalCase = {};
    var synOf = {};
    docs.forEach(function(d) {
      if (d.db_type != 'core') {
        otherfeaturesId[d._id.toUpperCase()] = d.taxon_id;
      }
      uniqueTaxa[d.taxon_id] = 1;
      uniqueId[d._id.toUpperCase()] = d.taxon_id;
      originalCase[d._id.toUpperCase()] = d._id;
      if (d.synonyms) {
        d.synonyms.forEach(function(syn) {
          uniqueId[syn.toUpperCase()] = d.taxon_id;
          originalCase[syn.toUpperCase()] = syn;
          synOf[syn.toUpperCase()] = d._id;
        });
      }
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
            var term = d.value; // .toUpperCase(); // don't do this it leads to problems
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
          var n=0;
          // find all the terms in each species
          Object.keys(uniqueTaxa).forEach(function(taxon_id) {
            var url = genesURL + '/query?rows=0&facet=true&facet.field=_terms&facet.limit=-1&json.nl=map&facet.mincount=1&q=taxon_id:'+taxon_id;
            console.error(url);
            var res = request('GET',url);
            var term_tally = JSON.parse(res.getBody()).facet_counts.facet_fields._terms;
            for (var t in term_tally) {
              var term = t; // .toUpperCase(); // don't do this it leads to problems
              if (!uniqueId.hasOwnProperty(term)) {
                if (term_freq.hasOwnProperty(term)) {
                  if (taxa_lut[term].hasOwnProperty(taxon_id)) {
                    taxa_lut[term][taxon_id] += term_tally[t];
                  }
                  else {
                    taxa_lut[term][taxon_id] = term_tally[t];
                  }
                }
                else if (/[A-Za-z]/.test(term)) { // this is a unique term that contains a word character
                  var solr = {
                    category    : 'Genes',
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
                  if (n===1) console.log('[');
                  else console.log(',');
                  console.log(JSON.stringify(solr));
                }
              }
            }
          });
          
          console.error("finished building taxa_lut. n=",n);
          for (var term in term_freq) {
            var tf = term_freq[term];
            var taxa = {ids:[],counts:[]};
            for(var taxon_id in taxa_lut[term]) {
              taxa.ids.push(+taxon_id);
              taxa.counts.push(taxa_lut[term][taxon_id]);
            }
            var solr = {
              category    : 'Genes',
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
            if (n===1) console.log('[');
            else console.log(',');
            console.log(JSON.stringify(solr));
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
              console.log(',');
              console.log(JSON.stringify(solr));
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
                console.log(',');
                console.log(JSON.stringify(solr));
              }
      
              // output unique IDs
              console.error('output uniqueIds');
              for (var uid in uniqueId) {
                console.log(',');
                console.log(JSON.stringify({
                  category : 'Genes',
                  subcategory : 'id',
                  fq_field : 'id',
                  fq_value : synOf[uid] || originalCase[uid],
                  id       : originalCase[uid],
                  display_name : originalCase[uid],
                  num_genes : 1,
                  relevance : 1.2,
                  taxon_id : uniqueId[uid],
                  taxon_freq : [1]
                }));
              }
              console.log(']');
            });
          });
        });
      });
    });
  });
});

