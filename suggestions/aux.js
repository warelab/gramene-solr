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
      display_name: doc.name,
      name: doc.name,
      def: doc.def,
      fq_field: 'GO__ancestors',
      fq_value: doc._id,
      num_genes: genes,
      relevance: 1.1 - 0.1/Math.sqrt(specificity) // prioritize less specific GO terms
    };
    if (doc.subset && doc.subset.indexOf('goslim_plant') !== -1) {
      solr.relevance+=0.1; // give goslim_plant terms a nice boost
    }
    optionalFields.forEach(function(f) {
      if (doc.hasOwnProperty(f)) {
        solr[f] = doc[f];
      }
    });
    return solr;
  },
  PO: function(doc,genes,specificity) {
    var categoryLabel = {
      plant_anatomy : 'Plant anatomy',
      plant_structure_development_stage: 'Plant structural/developmental stage'
    };
    var solr = {
      category: categoryLabel[doc.namespace], // 'Plant ontology',
      int_id: doc._id,
      id: doc.id,
      display_name: doc.name,
      name: doc.name,
      def: doc.def,
      fq_field: 'PO__ancestors',
      fq_value: doc._id,
      num_genes: genes,
      relevance: 1.1 - 0.1/Math.sqrt(specificity) // prioritize less specific terms
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
          return rank[1];
        }
      }
      return '';
    }
    var solr = {
      category: 'Taxonomy',
      subcategory: getRank(doc),
      int_id: doc._id,
      id: doc.id,
      display_name: doc.name,
      name: doc.name,
      fq_field: 'taxonomy__ancestors',
      fq_value: doc._id,
      num_genes: genes,
      relevance: 1.0 + 0.1/Math.sqrt(specificity) // prioritize more specific terms
    };
    if (doc._id === 3702) { // hard coded boost for arabidopsis thaliana (over lyrata subsp, lyrata)
      solr.relevance *= 1.1;
    }
    if (doc.hasOwnProperty('synonym')) {
      solr.synonym = doc.synonym;
    }
    return solr;
  },
  domains: function(doc,genes,specificity) {
    var solr = {
      category: 'InterPro',
      subcategory: doc.type,
      int_id: doc._id,
      id: doc.id,
      display_name: doc.name,
      name: doc.name,
      description: doc.description,
      abstract: doc.abstract,
      xref: [],
      fq_field: 'domains__ancestors',
      fq_value: doc._id,
      num_genes: genes,
      relevance: 1.05 + 0.1/Math.sqrt(specificity) // prioritize more specific terms
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
      display_name: doc.name,
      name: doc.name,
      synonym: doc.synonyms,
      fq_field: 'pathways__ancestors',
      fq_value: doc._id,
      num_genes: genes,
      relevance: 1.2 - 0.1/Math.sqrt(specificity) // prioritize 
    };
    return solr;
  }
};

var promises = _.map(mongo2solr, function(f,key) {
  console.error(key,'started');
  // do a facet query on the genes core on the  key__ancestors field
  var url = genesURL
    + '/query?rows=0&facet=true&facet.limit=-1&facet.mincount=1&json.nl=map'
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
        var solrDocs = _.filter(docs, function(doc) {
          return assoc.hasOwnProperty(doc._id);
        }).map(function(doc) {
          return f(doc, assoc[doc._id], termSpecificity[doc._id]);
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
