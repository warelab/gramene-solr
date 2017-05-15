#!/usr/bin/env node
var http = require("q-io/http")
  , _ = require('lodash')
  , fs = require('fs')
  , collections = require('gramene-mongodb-config');

var genesURL = process.argv[2];

function getURL(fq_field,fq_value,category,label) {
  var fq = fq_field + ':' + fq_value;
  var filters = {};
  filters[fq] = {
    category: category,
    display_name: label,
    fq: fq,
    exclude: false
  };
  var hash = encodeURI(JSON.stringify({
    filters: filters,
    taxa: {}
  }));
  return 'http://search.gramene.org/#' + hash;
}

var taxLUT;
collections.taxonomy.mongoCollection().then(function(tax) {
  tax.find().toArray(function(err,docs) {
    if (err) throw err;
    collections.closeMongoDatabase();
    taxLUT = _.keyBy(docs,'_id');


    // get the genes docs from solr
    var url = genesURL
    + '/query?q=taxonomy__ancestors:147368&rows=1000000&fl=id,name,biotype,description,taxon_id,annotations';
    return http.read(url).then(function(databuffer) {
      var data = JSON.parse(databuffer);
      var solrDocs = data.response.docs.map(function(doc) {
        var solr = {
          entry_type: 'Gene',
          database_name: 'Gramene',
          db_id: doc.id,
          species: taxLUT[doc.taxon_id].name,
          description: [doc.name, doc.biotype],
          url: getURL('id',doc.id,'Gene',doc.name)
        };
        if (doc.description !== 'unknown') {
          solr.description.push(doc.description);
        }
        if (doc.annotations) {
          solr.description.push(doc.annotations);
        }
        solr.description = _.uniq(_.flattenDeep(solr.description));
        return solr;
      });
      fs.writeFile('genes.json',JSON.stringify(solrDocs,null,'  '), function(err) {
        if (err) throw(err);
        console.error('wrote to json file');
      });
    })
  })
})
