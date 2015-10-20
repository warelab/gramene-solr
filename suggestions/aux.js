#!/usr/bin/env node

var http = require("q-io/http");
var key = process.argv[2];
var genesURL = process.argv[3];

// do a facet query on the genes core on the  key__ancestors field
var url = genesURL
  + '/query?rows=0&facet=true&facet.limit=-1&json.nl=map'
  + '&facet.field=' + key + '__ancestors';

http.read(url).then(function(databuffer) {
  var data = JSON.parse(databuffer);

  var assoc = data.facet_counts.facet_fields[key+'__ancestors'];
  
  // setup reader
  var n=0;
  require('readline').createInterface(
    {
      input: process.stdin,
      terminal: false
    }
  ).on('line', function(line) { // one JSON object per line
    var mongo = JSON.parse(line);
    var solr = mongo2solr[key](mongo,assoc);
    if (n===0) console.log('[');
    else console.log(',');
    console.log(JSON.stringify(solr));
    n++;
  }).on('close', function() {
    console.log(']');
  });
});

var optionalFields = ['comment','xref','synonym'];
var mongo2solr = {
  GO: function(doc,assoc) {
    var solr = {
      category: doc.namespace, // biological_process, molecular_function, cellular_component
      int_id: doc._id,
      id: doc.id,
      name: doc.name,
      def: doc.def,
      fqField: 'GO__ancestors',
      _genes: assoc.hasOwnProperty(doc._id) ? assoc[doc._id] : 0,
      relevance: assoc.hasOwnProperty(doc._id) ?
        doc.ancestors.length/100 // more weight to more specific terms
      : 0 // demote GO terms with no genes associated with them
    };
    optionalFields.forEach(function(f) {
      if (doc.hasOwnProperty(f)) {
        solr[f] = doc[f];
      }
    });
    return solr;
  },
  PO: function(doc,assoc) {
    var solr = {
      category: doc.namespace, // plant_anatomy plant_structural_developmental_stage
      int_id: doc._id,
      id: doc.id,
      name: doc.name,
      def: doc.def,
      fqField: 'PO__ancestors',
      _genes: assoc.hasOwnProperty(doc._id) ? assoc[doc._id] : 0,
      relevance: assoc.hasOwnProperty(doc._id) ?
        doc.ancestors.length/100 // more weight to more specific terms
      : 0 // demote GO terms with no genes associated with them
    };
    optionalFields.forEach(function(f) {
      if (doc.hasOwnProperty(f)) {
        solr[f] = doc[f];
      }
    });
    return solr;
  },
  taxonomy: function(doc,assoc) {
    var solr = {
      category: doc.namespace, // ncbi_taxonomy
      int_id: doc._id,
      id: doc.id,
      name: doc.name,
      fqField: 'taxonomy__ancestors',
      _genes: assoc.hasOwnProperty(doc._id) ? assoc[doc._id] : 0,
      relevance: assoc.hasOwnProperty(doc._id) ?
        doc.ancestors.length/100
      : 0
    };
    if (doc._id === 3702) { // hard coded boost for arabidopsis thaliana (over lyrata subsp, lyrata)
      solr.relevance *= 1.2;
    }
    if (doc.hasOwnProperty('synonym')) {
      solr.synonym = doc.synonym;
    }
    return solr;
  },
  domains: function(doc,assoc) {
    var solr = {
      category: doc.type, // Active_site Binding_site Conserved_site Domain Family PTM Repeat
      int_id: doc._id,
      id: doc.id,
      name: doc.name,
      description: doc.description,
      abstract: doc.abstract,
      xref: [],
      fqField: 'domains__ancestors',
      _genes: assoc.hasOwnProperty(doc._id) ? assoc[doc._id] : 0,
      relevance: assoc.hasOwnProperty(doc._id) ?
        doc.ancestors.length/100
      : 0
    };
    for (var f in doc) {
      if (!(solr.hasOwnProperty(f) || f === 'ancestors' || f === 'type')) {
        if (Array.isArray(doc[f])) {
          Array.prototype.push.apply(solr.xref,doc[f]);
        }
        else {
          solr.xref.push(doc[f]);
        }
      }
    }
    return solr;
  },
  pathways: function(doc,assoc) {
    var solr = {
      category: doc.type, // Reaction or Pathway
      int_id: doc._id,
      id: doc.id,
      name: doc.name,
      synonym: doc.synonyms,
      fqField: 'pathways__ancestors',
      _genes: assoc.hasOwnProperty(doc._id) ? assoc[doc._id] : 0,
      relevance: assoc.hasOwnProperty(doc._id) ?
        doc.ancestors.length/100
      : 0
    };
    return solr;
  }
};
