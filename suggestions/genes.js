#!/usr/bin/env node
var http = require("q-io/http")
  , collections = require('gramene-mongodb-config');

var genesURL = process.argv[2];

// read all of the unique ids from the genes collection
// to avoid suggesting them as non-unique terms when they are mentioned as an xref or something in another gene
collections.genes.mongoCollection().then(function(collection) {
  collection.find({},{_id:1}).toArray(function(err,docs) {
    if (err) throw err;
    collections.closeMongoDatabase();
    var uniqueId = {};
    docs.forEach(function(d) {
      uniqueId[d._id.toUpperCase()] = 1;
    });

    // non-unique terms get their own doc with a relevant fq field
    // build a lookup table with the non-unique ids from these
    var url = genesURL + '/query?rows=0&facet=true&facet.limit=-1&json.nl=map&facet.mincount=2&facet.field=_terms';
    console.error(url);
    http.read(url).then(function(data) {
      var term_freq = JSON.parse(data).facet_counts.facet_fields._terms;
      var n=0;
      for (var term in term_freq) {
        var ucTerm = term.toUpperCase();
        term_freq[ucTerm] = term_freq[term];
        if (!uniqueId.hasOwnProperty(ucTerm)) {
          var tf = term_freq[ucTerm];
          var solr = {
            category    : 'Gene',
            id          : '_term_'+ ++n,
            display_name: term,
            name        : term,
            fq_field    : '_terms',
            fq_value    : term,
            num_genes   : tf,
            relevance   : tf > 100 ? 1 : tf/100
          };
          if (n===1) console.log('[');
          else console.log(',');
          console.log(JSON.stringify(solr));
        }
      }
      
      // gene tree suggestion
      var url = genesURL + '/query?rows=0&facet=true&facet.limit=-1&json.nl=map&facet.mincount=1&facet.field=grm_gene_tree';
      console.error(url);
      http.read(url).then(function(data) {
        var tree_size = JSON.parse(data).facet_counts.facet_fields.grm_gene_tree;
        for (var tree in tree_size) {
          var solr = {
            category    : 'Gene tree',
            id          : '_term_'+ ++n,
            display_name: tree,
            name        : tree,
            fq_field    : 'grm_gene_tree',
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
            var solr = {
              category    : 'Biotype',
              id          : '_term_'+ ++n,
              display_name: biotype,
              name        : biotype,
              fq_field    : 'biotype',
              fq_value    : biotype,
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
              relevance : 0,
              taxon_id : mongo.taxon_id
            }));

            // add uniquely identifying xrefs
            var xref_h = {};
            for (var db in mongo.xrefs) {
              if (!collections.hasOwnProperty(db)) {
                mongo.xrefs[db].forEach(function(xr) {
                  xref_h[xr]=db;
                });
              }
            }
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
                relevance : 0,
                taxon_id : mongo.taxon_id
              }));
            });

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
                  relevance : 0,
                  taxon_id : mongo.taxon_id
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
