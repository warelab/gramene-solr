#!/usr/bin/env node
var collections = require('gramene-mongodb-config');

collections.genes.mongoCollection().then(function(collection) {
  var cursor = collection.find().sort([{'taxon_id':1},{'location.region':1},{'location.start':1}]);
  var n=0;
  var current_taxon=0;
  var gene_offset = 0;
  cursor.each(function(err,mongo) {
    if (err) throw err;
    if (mongo == null) {
      console.log(']');
      collections.closeMongoDatabase();
    }
    else {
      var location = mongo.location;
      if (!mongo.description) {
        mongo.description = 'unknown';
      }

      var solr = {
        // required data for list view
        id          : mongo._id,
        name        : mongo.name,
        description : mongo.description.replace(/\s+\[Source:.*/,''), // strip off the [Source:...]
        taxon_id    : mongo.taxon_id,

        // fields useful for genomic interval queries
        map         : location.map,
        region      : location.region,
        start       : location.start,
        end         : location.end,
        strand      : location.strand,

        // additional field(s) for query/faceting
        biotype : mongo.biotype,
        synonyms: mongo.synonyms,
        
        // so we know if it's coming from core or otherfeatures
        db_type : mongo.db_type,
        system_name : mongo.system_name,
        
        capabilities : ['location']
      };

      if (current_taxon !== mongo.taxon_id) {
        gene_offset = 0;
        current_taxon = mongo.taxon_id;
      }
      else {
        gene_offset++;
      }

      solr.gene_idx = gene_offset;

      // representative homolog (for display purposes)
      if (mongo.hasOwnProperty('representative') && mongo.representative.id !== solr.id) {
        solr.rep_id = mongo.representative.id;
        solr.rep_taxon_id = mongo.representative.taxon_id;
        if (mongo.representative.hasOwnProperty('name')) {
          solr.rep_name = mongo.representative.name;
        }
        if (mongo.representative.hasOwnProperty('description')) {
          solr.rep_desc = mongo.representative.description.replace(/\s+\[Source:.*/,'');
        }
      }

      // facet counting on bin fields drives the taxagenomic distribution
      for (var field in mongo.bins) {
        solr[field + '__bin'] = mongo.bins[field];
      }

      // homology fields
      if (mongo.hasOwnProperty('homology')) {
        solr.capabilities.push('homology');
        solr.grm_gene_tree_root_taxon_id = mongo.grm_gene_tree_root_taxon_id;
        solr.epl_gene_tree = mongo.epl_gene_tree;
        solr.grm_gene_tree = mongo.grm_gene_tree;
        if (mongo.hasOwnProperty('epl_sibling_trees')) {
          solr.epl_sibling_trees = mongo.epl_sibling_trees;
        }
        for (htype in mongo.homology) {
          solr['homology__'+htype] = mongo.homology[htype];
        }
      }

      // protein annotation fields
      if (mongo.hasOwnProperty("canonical_translation")) {
        if (!!mongo.canonical_translation.domain_roots) {
          solr.domain_roots = mongo.canonical_translation.domain_roots;
        }
        ['avg_res_weight','charge','iso_point','length','molecular_weight'].forEach(function(fname) {
          solr['protein__'+fname] = mongo.canonical_translation[fname];
        });
      }

      // transcript properties
      if (mongo.hasOwnProperty('canonical_transcript')) {
        solr.transcript__length = mongo.canonical_transcript.length;
        solr.transcript__exons = mongo.canonical_transcript.exons.length;
      }

      // convert xrefs:{db:[list]} to db__xrefs:[list]
      var ancestorFields = [];
      var hasXrefs = false;
      for (var db in mongo.xrefs) {
        if (collections.hasOwnProperty(db)) { // except for these
          ancestorFields.push(db);
        }
        else {
          solr[db + '__xrefs'] = mongo.xrefs[db];
          hasXrefs=true;
        }
      }
      if (hasXrefs) {
        solr.capabilities.push('xrefs');
      }
      ancestorFields.forEach(function(f) {
        var solrField = f + '__ancestors';
        solr.capabilities.push(f);
        solr[solrField] = mongo.xrefs[f];
        if (mongo.ancestors.hasOwnProperty(f)) {
          mongo.ancestors[f].forEach(function(r) {
            solr[solrField].push(r);
          });
        }
      });
      // special case for the ancestors of the grm_gene_tree_root_taxon_id
      if (mongo.ancestors.hasOwnProperty('gene_family')) {
        solr['gene_family__ancestors'] = mongo.ancestors.gene_family;
      }

      if (n===0) console.log('[');
      else console.log(',');
      console.log(JSON.stringify(solr));
      n++;
    }
  });
});
