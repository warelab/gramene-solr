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
      molecular_function : 'Gene Ontology: molecular function',
      biological_process : 'Gene Ontology: biological process',
      cellular_component : 'Gene Ontology: cellular component'
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
    if (doc.subset && doc.subset.indexOf('goslim_plant') === -1) {
      solr.relevance-=0.05; // penalize non-goslim_plant terms
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
      plant_anatomy : 'Plant Ontology: anatomy',
      plant_structure_development_stage: 'Plant Ontology: structural/developmental stage'
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
  TO: function(doc,genes,specificity) {
    var solr = {
      category: 'Trait ontology',
      int_id: doc._id,
      id: doc.id,
      display_name: doc.name,
      name: doc.name,
      def: doc.def,
      fq_field: 'TO__ancestors',
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
  QTL_TO: function(doc,genes,specificity) {
    var solr = {
      category: 'QTLs',
      int_id: doc._id,
      id: doc.id,
      display_name: doc.name,
      name: doc.name,
      def: doc.def,
      fq_field: 'QTL_TO__ancestors',
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
    const rank = doc.rank === "no rank" ? '' : doc.rank;
    var solr = {
      category: 'Taxonomy',
      subcategory: rank,
      int_id: doc._id,
      id: doc.id,
      display_name: `${doc.name} (${rank})`,
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
      category: `InterPro: ${doc.type}`,
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
      category: `Plant Reactome: ${doc.type}`,
      int_id: doc._id,
      id: 'R-OSA-'+doc._id,
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
  + '&facet.pivot=' + key + '__ancestors,taxon_id';
  console.error(url);
  return http.read(url).then(function(databuffer) {
    console.error(key,'got facets');
    var data = JSON.parse(databuffer);
    var deferred = Q.defer();

    var numGenes = {};
    var taxa = {};//{ids:[],counts:[]};
    data.facet_counts.facet_pivot[key+'__ancestors,taxon_id'].forEach(function(d) {
      numGenes[d.value] = d.count;
      taxa[d.value] = {ids:[],counts:[]};
      d.pivot.forEach(function(p) {
        taxa[d.value].ids.push(p.value);
        taxa[d.value].counts.push(p.count);
      });
    });
    let fkey = key;
    if (key === "QTL_TO") {
      // generate suggestions from the qtls collection that do range queries
      collections.qtls.mongoCollection().then(function(qtls_coll) {
        qtls_coll.find().toArray(function(err, qtls) {
          if (err) deferred.reject(err);
          if (!qtls) {
            console.error('qtls not found');
            deferred.reject('qtls not found');
          }
          console.error('got '+qtls.length + ' qtls from mongo');
          var solrDocs = qtls.map(function(qtl) {
            const l = qtl.location;
            let syns = [];
            const id_ver = qtl._id.split('.');
            if (id_ver.length === 2) {
              syns.push(id_ver[0]);
              const digits = /[0-9]+$/;
              const pre_chr = id_ver[0].replace(digits,'');
              if (pre_chr !== id_ver[0]) {
                syns.push(pre_chr)
              }
            }
            var solr = {
              category: 'QTLs',//qtl.source.replaceAll('_',' '),
              id: qtl._id,
              display_name: qtl._id,
              name: qtl.description,
              synonym: syns,
              fq_field: 'location',
              fq_value: `(map:${l.map} AND region:${l.region} AND start:[${l.start} TO ${l.end}])`,
              num_genes: 1,
              relevance: 1.0,
              taxon_id: [4558,39947],
              taxon_freq: [1,1]
            };
            return solr
          })
          
          fs.writeFile('qtls.json', JSON.stringify(solrDocs, null, '  '), function(err) {
            if (err) deferred.reject(err);
            console.error('qtls written to json');
            // deferred.resolve(true);
          })
        })
      })
      // change the key to TO so we get QTL suggestions
      key = "TO";
    }
    // instead of reading from stdin, get the docs from mongodb
    collections[key].mongoCollection().then(function(collection) {
      collection.find().toArray(function(err,docs) {
        if (err) deferred.reject(err);
        if (!docs) {
          console.error('docs is not defined '+ key);
          deferred.reject('docs undefined for collection '+ key);
        }
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
          return numGenes.hasOwnProperty(doc._id);
        }).map(function(doc) {
          return f(doc, numGenes[doc._id], termSpecificity[doc._id]);
        }).map(function(doc) {
          doc.taxon_id = taxa[doc.int_id].ids;
          doc.taxon_freq = taxa[doc.int_id].counts;
          return doc;
        });
        
        fs.writeFile(fkey + '.json',JSON.stringify(solrDocs,null,'  '), function(err) {
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
