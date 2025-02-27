#!/usr/bin/env node

var genesURL = process.argv[2];
var collections = require("gramene-mongodb-config");
var fs = require('fs');

collections.germplasm.mongoCollection().then(function(coll) {
  coll.find().toArray(function(err,germplasm) {
    if (err) throw(err);
    console.error(`got ${germplasm.length} accessions from mongo`);
    // collections.closeMongoDatabase();
    fetchData('EMS', germplasm);
    fetchData('NAT', germplasm);
  })
});


async function fetchData(dtype, germplasm) {
  var url = genesURL + '/query?rows=0&facet=true&facet.limit=-1&facet.mincount=1&json.nl=map'
  + '&facet.pivot=VEP__merged__'+dtype+'__attr_ss,taxon_id';
  console.error(url);
  try {
    const response = await fetch(url);
    if (!response.ok) {
      throw new Error('Network response was not ok');
    }
    const data = await response.json();
    // console.error('got data',data.facet_counts.facet_pivot);
    var numGenes = {};
    var taxa = {};
    data.facet_counts.facet_pivot[`VEP__merged__${dtype}__attr_ss,taxon_id`].forEach(function(d) {
      numGenes[d.value] = d.count;
      taxa[d.value] = {ids:[],counts:[]};
      d.pivot.forEach(function(p) {
        taxa[d.value].ids.push(p.value);
        taxa[d.value].counts.push(p.count);
      })
    });
    var solrDocs = germplasm.filter(g => numGenes.hasOwnProperty(g.ens_id)).map(g => {
      var solr = {
        category: 'Germplasm PTVs',
        subcategory: "Natural variation",
        id: g.ens_id,
        display_name: g.pub_id,
        name: g.pub_id,
        fq_field: `VEP__merged__${dtype}__attr_ss`,
        fq_value: g.ens_id,
        num_genes: numGenes[g.ens_id],
        relevance: 1.0,
        taxon_id: taxa[g.ens_id].ids,
        taxon_freq: taxa[g.ens_id].counts
      }
      if (dtype === "EMS") {
        solr.subcategory = g.pop_id === "1" ? "EMS (Purdue)" : "EMS (Lubbock)";
        if (g.pop_id === "2") {
          solr.display_name = `${g.ens_id} (${g.pub_id})`
        }
      }
      return solr
    })
    fs.writeFile(`${dtype}_germplasm.json`, JSON.stringify(solrDocs, null, '  '), function(err) {
      if (err) throw new Error(err);
      console.error(`${dtype} written to json`);
    })
  } catch (error) {
    console.error('error fetching data:', error);
  }
}