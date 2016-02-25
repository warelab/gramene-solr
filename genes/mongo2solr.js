#!/usr/bin/env node
var collections = require('gramene-mongodb-config');

function get_rep(c) {
  var rep = {
    id: c.id,
    taxon_id: c.taxon_id
  };
  if (c.hasOwnProperty('name')) {
    rep.name = c.name;
  }
  if (c.hasOwnProperty('description')) {
    rep.description = c.description.replace(/\s+\[Source:.*/,'');
  }
  return rep;
}

collections.genes.mongoCollection().then(function(collection) {
  var cursor = collection.find().sort([{'db_type':1},{'taxon_id':1},{'gene_idx':1}]);
  var n=0;
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
        gene_idx    : mongo.gene_idx,

        // additional field(s) for query/faceting
        biotype : mongo.biotype,
        synonyms: mongo.synonyms,
        
        // so we know if it's coming from core or otherfeatures
        db_type : mongo.db_type,
        system_name : mongo.system_name,
        
        capabilities : ['location']
      };

      // facet counting on bin fields drives the taxagenomic distribution
      for (var field in mongo.bins) {
        solr[field + '__bin'] = mongo.bins[field];
      }

      // homology fields
      if (mongo.homology) {
        solr.capabilities.push('homology');
        // representative homolog (for display purposes)
        if (mongo.homology.gene_tree) {
          solr.gene_tree_root_taxon_id = mongo.homology.gene_tree.root_taxon_id;
          solr.gene_tree = mongo.homology.gene_tree.id;
          if (mongo.homology.gene_tree.representative) {
            var mhgr = mongo.homology.gene_tree.representative;
            if (mhgr.closest) {
              var rep = get_rep(mhgr.closest);
              for (var f in rep) {
                solr['closest_rep_'+f] = rep[f];
              }
            }

            if (mhgr.hasOwnProperty('model')) {
              var rep = get_rep(mhgr.model);
              for (var f in rep) {
                solr['model_rep_'+f] = rep[f];
              }
            }
          }
        }

        
        solr.homology__all_orthologs = [mongo._id];
        solr.homology__all_homeologs = [mongo._id];
        for (htype in mongo.homology) {
          if (Array.isArray(mongo.homology[htype])) {
            solr['homology__'+htype] = mongo.homology[htype];
            if (htype.match(/ortholog/)) {
              mongo.homology[htype].forEach(function(o) {
                solr.homology__all_orthologs.push(o);
              });
            }
            if (htype.match(/homeolog/)) {
              mongo.homology[htype].forEach(function(h) {
                solr.homology__all_homeologs.push(h);
              });
            }
          }
        }
        if (solr.homology__all_homeologs.length === 1) {
          delete solr.homology__all_homeologs;
        }
        if (solr.homology__all_orthologs.length === 1) {
          delete solr.homology__all_orthologs;
        }
        if (solr.hasOwnProperty('homology__within_species_paralog')) {
          solr.homology__within_species_paralog.push(mongo._id);
        }
        if (solr.hasOwnProperty('homology__gene_split')) {
          solr.homology__gene_split.push(mongo._id);
        }
      }

      // check if canonical transcript's translation has a domain architecture
      var ct = mongo.gene_structure.canonical_transcript;
      if (!(ct && mongo.gene_structure.transcripts[ct])) {
        console.error("no canonical transcript in gene",mongo._id);
        process.exit(2);
      }
      if (mongo.gene_structure.transcripts[ct].translation) {
        var tl = mongo.gene_structure.transcripts[ct].translation;
        if (tl.features.domain) {
          solr.domain_roots = tl.features.domain.roots;
        }
        solr.protein__length = tl.length;
      }

      // canonical transcript properties
      solr.transcript__length = mongo.gene_structure.transcripts[ct].length;
      solr.transcript__exons = mongo.gene_structure.transcripts[ct].exons.length;

      solr.transcript__count = Object.keys(mongo.gene_structure.transcripts).length;

      // convert xrefs:{db:[list]} to db__xrefs:[list]
      var hasXrefs = false;
      for (var db in mongo.xrefs) {
        solr[db + '__xrefs'] = mongo.xrefs[db];
        hasXrefs=true;
      }
      if (hasXrefs) {
        solr.capabilities.push('xrefs');
      }

      // add ancestors fields from the annotations section
      for (var f in mongo.annotations) {
        if (mongo.annotations[f].ancestors) {
          solr.capabilities.push(f);
          solrField = f + '__ancestors';
          solr[solrField] = mongo.annotations[f].ancestors;
        }
      }

      if (n===0) console.log('[');
      else console.log(',');
      console.log(JSON.stringify(solr));
      n++;
    }
  });
});
