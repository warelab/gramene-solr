#!/usr/bin/env node

var http = require('q-io/http')
  , fs = require('fs')
  , collections = require('gramene-mongodb-config')
  , _ = require('lodash')
  , Q = require('q');

var genesURL = process.argv[2];

var optionalFields = ['comment','xref','synonym'];
var mongo2solr = {
  GO: function(doc,genes,specificity) {
    var categoryLabel = {
      molecular_function : 'GO function',
      biological_process : 'GO process',
      cellular_component : 'GO component'
    };
    var solr = {
      category: categoryLabel[doc.namespace],
      int_id: doc._id,
      id: doc.id,
      displayName: doc.name,
      name: doc.name,
      def: doc.def,
      fqField: 'GO__ancestors',
      _genes: genes,
      relevance: genes
        ? 1/Math.sqrt(specificity) // prioritize more specific terms
        : -0.5 // penalize suggestions without genes
    };
    optionalFields.forEach(function(f) {
      if (doc.hasOwnProperty(f)) {
        solr[f] = doc[f];
      }
    });
    return solr;
  },
  PO: function(doc,genes,specificity) {
    var solr = {
      category: 'Plant ontology', //doc.namespace, // plant_anatomy plant_structural_developmental_stage
      int_id: doc._id,
      id: doc.id,
      displayName: doc.name,
      name: doc.name,
      def: doc.def,
      fqField: 'PO__ancestors',
      _genes: genes,
      relevance: genes
        ? 1/Math.sqrt(specificity) // prioritize more specific terms
        : -0.5 // penalize suggestions without genes
    };
    optionalFields.forEach(function(f) {
      if (doc.hasOwnProperty(f)) {
        solr[f] = doc[f];
      }
    });
    return solr;
  },
  taxonomy: function(doc,genes,specificity) {
    function getRank(doc) {
      if (doc.hasOwnProperty('property_value')) {
        var rank = doc.property_value.match(/has_rank NCBITaxon:(.*)/);
        if (rank.length===2) {
          return ' (' + rank[1] + ')';
        }
      }
      return '';
    }
    var solr = {
      category: 'Taxonomy', //doc.namespace, // ncbi_taxonomy
      int_id: doc._id,
      id: doc.id,
      displayName: doc.name + getRank(doc),
      name: doc.name,
      fqField: 'taxonomy__ancestors',
      _genes: genes,
      relevance: genes
        ? 1/Math.sqrt(specificity) // prioritize more specific terms
        : -0.5 // penalize suggestions without genes
    };
    if (doc._id === 3702) { // hard coded boost for arabidopsis thaliana (over lyrata subsp, lyrata)
      solr.relevance *= 1.2;
    }
    if (doc.hasOwnProperty('synonym')) {
      solr.synonym = doc.synonym;
    }
    return solr;
  },
  domains: function(doc,genes,specificity) {
    var solr = {
      category: 'InterPro', //doc.type, // Active_site Binding_site Conserved_site Domain Family PTM Repeat
      int_id: doc._id,
      id: doc.id,
      displayName: doc.name + ' (' + doc.type + ')',
      name: doc.name,
      description: doc.description,
      abstract: doc.abstract,
      xref: [],
      fqField: 'domains__ancestors',
      _genes: genes,
      relevance: genes
        ? 1/Math.sqrt(specificity) // prioritize more specific terms
        : -0.5 // penalize suggestions without genes
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
  pathways: function(doc,genes,specificity) {
    var solr = {
      category: 'Plant Reactome',
      int_id: doc._id,
      id: doc.id,
      displayName: doc.name,
      name: doc.name,
      synonym: doc.synonyms,
      fqField: 'pathways__ancestors',
      _genes: genes,
      relevance: genes
        ? 1/Math.sqrt(specificity) // prioritize more specific terms
        : -0.5 // penalize suggestions without genes
    };
    return solr;
  }
};

var promises = _.map(mongo2solr, function(f,key) {
  console.error(key,'started');
  // do a facet query on the genes core on the  key__ancestors field
  var url = genesURL
    + '/query?rows=0&facet=true&facet.limit=-1&json.nl=map'
    + '&facet.field=' + key + '__ancestors';
  return http.read(url).then(function(databuffer) {
    console.error(key,'got facets');
    var data = JSON.parse(databuffer);
    var deferred = Q.defer();

    var assoc = data.facet_counts.facet_fields[key+'__ancestors'];
    // instead of reading from stdin, get the docs from mongodb
    collections[key].mongoCollection().then(function(collection) {
      collection.find().toArray(function(err,docs) {
        if (err) deferred.reject(err);
        console.error('got '+docs.length + ' docs from collection '+key);
        // leaf nodes get the maximum boost
        // so for each node we need to know how many nodes are below it (descendants)
        // and how many are above it (ancestors)
        // count frequency of doc ids in the ancestors fields.
        var termSpecificity = {};
        docs.forEach(function(doc) {
          doc.ancestors.forEach(function(id) {
            if (! termSpecificity.hasOwnProperty(id)) {
              termSpecificity[id]=0;
            }
            termSpecificity[id]++;
          });
        });
        var solrDocs = docs.map(function(doc) {
          return f(doc,
            assoc.hasOwnProperty(doc._id) ? assoc[doc._id] : 0,
            termSpecificity[doc._id]);
        });
        fs.writeFile(key + '.json',JSON.stringify(solrDocs,null,'  '), function(err) {
          if (err) deferred.reject(err);
          console.error(key,'wrote to json file');
          deferred.resolve(true);
        });
      });
    });

    return deferred.promise;
  });
});

Q.all(promises).then(function(arrayOfTrues) {
  collections.closeMongoDatabase();
});
