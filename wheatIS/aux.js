#!/usr/bin/env node

var http = require('q-io/http')
  , fs = require('fs')
  , collections = require('gramene-mongodb-config')
  , _ = require('lodash')
, Q = require('q');

var genesURL = process.argv[2];

var optionalFields = ['comment','xref','synonym'];

function getURL(fq_field,fq_value,category,label) {
  var fq = fq_field + ':' + fq_value;
  var filters = {};
  filters[fq] = {
    category: category,
    display_name: label,
    fq: fq,
    exclude: false
  };
  filters['taxonomy__ancestors:147368'] = {
    category: 'Taxonomy',
    display_name: 'Pooideae',
    fq: 'taxonomy__ancestors:147368',
    exclude: false
  };
  var hash = encodeURI(JSON.stringify({
    filters: filters,
    taxa: {}
  }));
  return 'http://search.gramene.org/#' + hash;
}

var mongo2solr = {
  GO: function(doc,genes,specificity) {
    var categoryLabel = {
      molecular_function : 'GO function',
      biological_process : 'GO process',
      cellular_component : 'GO component'
    };
    var solr = {
      entry_type: categoryLabel[doc.namespace],
      database_name: 'Gramene',
      id: doc.id,
      int_id: doc._id,
      db_id: doc.id,
      description: [doc.name,doc.def],
      url: getURL('GO__ancestors',doc._id, categoryLabel[doc.namespace], doc.name)
    };
    optionalFields.forEach(function(f) {
      if (doc.hasOwnProperty(f)) {
        solr.description.push(doc[f]);
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
      entry_type: categoryLabel[doc.namespace], // 'Plant ontology',
      database_name: 'Gramene',
      id: doc.id,
      int_id: doc._id,
      db_id: doc.id,
      description: [doc.name,doc.def],
      url: getURL('PO__ancestors',doc._id, categoryLabel[doc.namespace], doc.name)
    };
    optionalFields.forEach(function(f) {
      if (doc.hasOwnProperty(f)) {
        solr.description.push(doc[f]);
      }
    });
    return solr;
  },
  taxonomy: function(doc,genes,specificity) {
    var solr = {
      entry_type: 'Taxonomy',
      database_name: 'Gramene',
      id: doc.id,
      int_id: doc._id,
      db_id: doc.id,
      description: [doc.name],
      url: getURL('taxonomy__ancestors',doc._id, 'Taxonomy', doc.name)
    };
    if (doc.hasOwnProperty('synonym')) {
      solr.description.push(doc.synonym);
    }
    if (doc.hasOwnProperty('property_value')) {
      var rank = doc.property_value.match(/has_rank NCBITaxon:(.*)/);
      if (rank.length===2) {
        solr.description.push(rank[1]);
      }
    }
    return solr;
  },
  domains: function(doc,genes,specificity) {
    var solr = {
      entry_type: 'InterPro ' + doc.type,
      database_name: 'Gramene',
      id: doc.id,
      int_id: doc._id,
      db_id: doc.id,
      description: [doc.name,doc.description,doc.abstract],
      url: getURL('domains__ancestors',doc._id, 'InterPro', doc.name)
    };
    return solr;
  },
  pathways: function(doc,genes,specificity) {
    var solr = {
      database_name: 'Gramene',
      entry_type: doc.type,
      id: doc.id,
      int_id: doc._id,
      db_id: doc.id,
      description: [doc.name,doc.synonyms],
      url: getURL('pathways__ancestors',doc._id, 'Plant Reactome', doc.name)
    };
    return solr;
  }
};

var promises = _.map(mongo2solr, function(f,key) {
  console.error(key,'started');
  // do a facet query on the genes core on the  key__ancestors field
  var url = genesURL
  + '/query?q=taxonomy__ancestors:147368&rows=0&facet=true&facet.limit=-1&facet.mincount=1&json.nl=map'
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
          doc.species = taxa[doc.int_id].ids.map(function(tid) {
            return taxLUT[tid].name;
          });
          doc.description = _.uniq(_.flattenDeep(doc.description));
          delete doc.int_id;
          delete doc.id;
          return doc;
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
var taxLUT;
collections.taxonomy.mongoCollection().then(function(tax) {
  tax.find().toArray(function(err,docs) {
    if (err) throw err;
    taxLUT = _.keyBy(docs,'_id');
    Q.all(promises).then(function(arrayOfTrues) {
      collections.closeMongoDatabase();
    });
  })
})
