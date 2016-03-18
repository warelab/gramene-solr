#!/usr/bin/env node
var collections = require('gramene-mongodb-config');
var _ = require('lodash');

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

Array.prototype.swap = function (x,y) {
  var b = this[x];
  this[x] = this[y];
  this[y] = b;
  return this;
}

collections.genes.mongoCollection().then(function(collection) {
  collection
  .find({taxon_id:45770004,db_type:'otherfeatures'},{_id:1,location:1,gene_idx:1})
  .toArray(function(err,projectedGenes) {
    if (err) throw err;
    console.error("we have",projectedGenes.length,"projectedGenes");
    // make sure they are sorted here because mongo wasn't having it
    for(var i=0;i<projectedGenes.length;i++) {
      if (i !== projectedGenes[i].gene_idx) projectedGenes.swap(i,projectedGenes[i].gene_idx);
    }
    
    // make a lookup table of v3 genes that were projected to v4
    var projected = {};
    var v3v4 = {};
    for(var i=0;i<projectedGenes.length;i++) {
      var doc = projectedGenes[i];
      doc.id = doc._id.replace(/_projected/,'');
      projected[doc.id] = 1;
      v3v4[doc.location.region] = v3v4[doc.location.region] || [];
      v3v4[doc.location.region].push(doc)
    }
    
    var doc_i = {};
    for (var region in v3v4) {
      doc_i[region] = 0;
    }
    var cursor = collection.find().sort([{'species_idx':1},{'db_type':1},{'gene_idx':1}]);
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

        if (mongo.system_name === "zea_mays" && !projected[mongo._id]) {
          solr.description += " -- not projected to v4";
        }

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
        var tIdx = _.keyBy(mongo.gene_structure.transcripts,'id');
        if (!(ct && tIdx[ct])) {
          console.error("no canonical transcript in gene",mongo._id);
          process.exit(2);
        }
        if (tIdx[ct].translation) {
          var tl = tIdx[ct].translation;
          if (tl.features && tl.features.domain) {
            solr.domain_roots = tl.features.domain.roots;
          }
          solr.protein__length = tl.length;
        }

        // canonical transcript properties
        solr.transcript__length = tIdx[ct].length;
        solr.transcript__exons = tIdx[ct].exons.length;

        solr.transcript__count = mongo.gene_structure.transcripts.length;

        // convert xrefs
        if (mongo.xrefs) {
          mongo.xrefs.forEach(function(xref) {
            solr[xref.db + '__xrefs'] = xref.ids;
          });
          solr.capabilities.push('xrefs');
        }

        // add ancestors fields from the annotations section
        for (var f in mongo.annotations) {
          if (mongo.annotations[f]) {
            solr.capabilities.push(f);
            solrField = f + '__ancestors';
            if (mongo.annotations[f].ancestors) {
              solr[solrField] = mongo.annotations[f].ancestors;
            }
            else {
              solr[solrField] = [];
            }
            if (mongo.annotations[f].entries) {
              mongo.annotations[f].entries.forEach(function(e) {
                if (e._id) {
                  solr[solrField].push(e._id);
                }
                else {
                  solr[solrField].push(parseInt(e.id.match(/\d+/)[0]));
                }
              })
            }
          }
        }

        // check for overlap with v3 projected genes
        if (mongo.system_name === "zea_mays4m" && mongo.db_type === "core") {
          // rewind until projected gene preceeds current gene
          function findOverlaps(i,v3,mongo,solr) {
            while (i > 0 && i < v3.length && v3[i].location.end > mongo.location.start) {
              i--;
            }
          // skip docs that end before this gene starts
            while (i < v3.length && v3[i].location.end < mongo.location.start) {
              i++;
            }
            // gather overlapping genes and add ids to synonyms
            while (i < v3.length && v3[i].location.start < mongo.location.end) {
              solr.synonyms = solr.synonyms || [];
              solr.synonyms.push(v3[i].id);
              i++;
            }
          };

          if (v3v4[mongo.location.region]) {
            findOverlaps(doc_i[mongo.location.region],v3v4[mongo.location.region],mongo,solr);
          }
        }
        if (n===0) console.log('[');
        else console.log(',');
        console.log(JSON.stringify(solr));
        n++;
      }
    });
  });
  
});
