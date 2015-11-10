#!/usr/bin/env node
var http = require("q-io/http")
  , collections = require('gramene-mongodb-config');

collections.taxonomy.mongoCollection().then(function(collection) {
  collection.find().toArray(function(err,docs) {
    if (err) throw err;
    collections.closeMongoDatabase();
    var taxonomy={};
    docs.forEach(function(doc) {
      taxonomy[doc._id] = doc;
    });
  
    var genesURL = process.argv[2];

    // non-unique terms get their own doc with a relevant fq field
    // build a lookup table with the non-unique ids from these
    var url = genesURL + '/query?rows=0&facet=true&facet.limit=-1&json.nl=map&facet.mincount=2&facet.field=_terms';
    console.error(url);
    http.read(url).then(function(data) {
      var term_freq = JSON.parse(data).facet_counts.facet_fields._terms;
      var n=0;
      for (var term in term_freq) {
        var solr = {
          category    : 'genes',
          id          : '_term_'+ ++n,
          display_name: term,
          name        : term,
          fq_field    : '_terms',
          fq_value    : term,
          _genes      : term_freq[term],
          relevance   : term_freq[term] > 100 ? 1 : term_freq[term]/100
        };
        if (n===1) console.log('[');
        else console.log(',');
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
        if (!mongo.description) {
          mongo.description = 'unknown';
        }
  
        var solr = {
          category     : 'genes',
          fq_field     : 'id',
          fq_value     : mongo._id,
          id           : mongo._id,
          description  : mongo.description.replace(/\s+\[Source:.*/,''), // strip off the [Source:...]
          display_name : mongo._id,
          _genes       : 1,
          relevance    : 0
        };
        if (mongo.name !== mongo._id && !term_freq.hasOwnProperty(mongo.name)) {
          solr.name = mongo.name;
          solr.relevance=1;
          solr.display_name += ' ' + mongo.name;
        }
        // append species name to displayName
        solr.display_name += ' ('+taxonomy[mongo.taxon_id].name+')';
        var xref_h = {};
        for (var db in mongo.xrefs) {
          if (!mongo.ancestors.hasOwnProperty(db)) { // aux cores
            mongo.xrefs[db].forEach(function(xr) {
              xref_h[xr]=db;
            });
          }
        }
        solr.xref = Object.keys(xref_h).filter(function(xr) {
          return !term_freq[xr];
        });

        if (mongo.hasOwnProperty('synonyms')) {
          solr.synonym = mongo.synonyms.filter(function(syn) {
            return !term_freq[syn];
          });
        }
  
        console.log(',');
        console.log(JSON.stringify(solr));
      }).on('close', function() {
        console.log(']');
      });
    });
  });
});