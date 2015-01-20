#!/usr/bin/env node
var request = require('request');
var settings = require('./config/settings.json');
var gramene = require('./gramene-lib.js');

var url = settings.hostname + ':' + settings.port + '/solr/taxonomy/select?q=*:*&rows=1000&wt=json&indent=true';
request.get(url,function(err,response, body) {
  var r = JSON.parse(response.body);

  var tax_tree = gramene.parseTaxa(r.response.docs);

  console.log(JSON.stringify(gramene.speciesTree(tax_tree),null,' '));
  
  console.log(gramene.rootTaxon(tax_tree, [3702]));
  console.log(gramene.rootTaxon(tax_tree, [3702,81972]));
  console.log(gramene.rootTaxon(tax_tree, [3702,29760]));
  console.log(gramene.rootTaxon(tax_tree, [39947,4528,4538]));
  console.log(gramene.rootTaxon(tax_tree, [39947,4528,4538,4577]));
});
