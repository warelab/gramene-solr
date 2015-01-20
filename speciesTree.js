#!/usr/bin/env node
var request = require('request');
var settings = require('./config/settings.json');

var url = settings.hostname + ':' + settings.port + '/solr/taxonomy/select?q=*:*&rows=1000&wt=json&indent=true';
request.get(url,function(err,response, body) {
  var r = JSON.parse(response.body);
  var taxes = r.response.docs;
  // console.log(taxes);
  var tax_node={};
  for (var i=0;i<taxes.length; i++) {
    var tax = taxes[i];
    tax_node[tax.id] = {
      name : tax.name_s,
      children : []
    };
    if (tax.hasOwnProperty('is_a_is')) {
      tax_node[tax.id].parent = tax.is_a_is[0];
    }
  }
  // console.log(JSON.stringify(tax_node));
  for(var id in tax_node) {
    if (tax_node[id].hasOwnProperty('parent')) {
      var p = tax_node[id].parent;
      tax_node[p].children.push(id);
    }
  }
  // don't write nodes with one child
  function bf (tree,id) {
    if (tree[id].children.length === 0) { // leaf node
      return {id:+id,name:tree[id].name};
    }
    var children = [];
    for(var child in tree[id].children) {
      children.push(bf(tree,tree[id].children[child]));
    }
    if (tree[id].children.length > 1) { // branching happens here
      return {
        id: +id,
        name: tree[id].name,
        children: children
      };
    }
    if (tree[id].children.length === 1) {
      return children[0];
    }
  }
  var ctree = bf(tax_node,"1",0);
  console.log(JSON.stringify(ctree,null,'  '));
});
