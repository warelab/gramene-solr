#!/usr/bin/env node
//
// NOT THE MAPPER THE BUILD RUNS. 60_solr_genes.sh runs mongo2solr.join.js, which joins the
// expression collection server-side instead of loading it into node memory (this one OOMed at
// release-70 scale). This file is kept only because older recipes referenced it.
//
// Edit mongo2solr.join.js, not this. Commit 49679f3 added the cultivar tree-id block HERE by
// mistake; the ids therefore never reached solr, every gate passed, and the omission survived two
// full rebuilds before anyone noticed the genes core exposed gene_tree only. The per-set supertree
// block added alongside this comment is likewise only in mongo2solr.join.js.
var collections = require('gramene-mongodb-config');
var _ = require('lodash');

function get_rep(c) {
  var rep = {
    id: c.id,
    taxon_id: c.taxon_id,
    identity: c.percent_identity
  };
  if (c.hasOwnProperty('name')) {
    rep.name = c.name;
  }
  if (c.hasOwnProperty('description')) {
    rep.description = c.description.replace(/\s+\[Source:.*/,'');
  }
  return rep;
}

collections.expression.mongoCollection().then(function(atlas) {
  atlas.find().toArray(function(err,docs) {
    console.error("got docs from atlas");
    var expressionData = _.keyBy(docs, '_id');
    collections.genes.mongoCollection().then(function(collection) {
      var cursor = collection.find({},{sort:{'species_idx':1,'db_type':1,'gene_idx':1}});
      var n=0;
      var p=0;
      var terminator = {}; // key is lowercase version, value is original term
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
            gene_idx    : n,
            compara_idx : p,

            // additional field(s) for query/faceting
            biotype : mongo.biotype,
            synonyms: mongo.synonyms || [],
            annotations: [],
        
            // so we know if it's coming from core or otherfeatures
            db_type : mongo.db_type,
            system_name : mongo.system_name,
        
            capabilities : ['location'],
          };

          function generateNeighborhood(start,n) {
            var arr = Array.apply(null, Array(n));
            return arr.map(function (x, i) { return start + i });
          }
          solr.neighbors_5 = generateNeighborhood(n-5, 2*5 + 1);
          solr.neighbors_10 = generateNeighborhood(n-10, 2*10 + 1);
          solr.neighbors_20 = generateNeighborhood(n-20, 2*20 + 1);
          solr.compara_neighbors_5 = generateNeighborhood(p-5, 2*5 + 1);
          solr.compara_neighbors_10 = generateNeighborhood(p-10, 2*10 + 1);
          solr.compara_neighbors_20 = generateNeighborhood(p-20, 2*20 + 1);

          // solr.description.split(/\s+/).forEach(function(w) {
          //   if (w.match(/^[a-z].*[0-9]$/i)) {
          //     solr.synonyms.push(w);
          //   }
          // });
          if (mongo.atlas_id) {
            solr.synonyms.push(mongo.atlas_id);
            solr.atlas_id = mongo.atlas_id;
          }

          // uniqify synonyms
          var uniq = {};
          var lcName = solr.name.toLowerCase();
          var lcId = solr.id.toLowerCase();
          uniq[lcName]=solr.name;
          uniq[lcId] = solr.id;
          solr.synonyms.forEach(function(syn) {
            var lc = syn.toLowerCase();
            if (!uniq.hasOwnProperty(lc)) {
              uniq[lc]=syn;
            }
          });
          delete uniq[lcName];
          delete uniq[lcId];
          if (Object.keys(uniq).length > 0) {
            solr.synonyms = Object.keys(uniq).map(function(k){return uniq[k]});
          }
          else {
            delete solr.synonyms;
          }
          // put the id and name back in
          uniq[lcName] = solr.name;
          uniq[lcId] = solr.id;
      
          if (mongo.summary) {
            solr.summary = mongo.summary;
          }

          // facet counting on bin fields drives the taxagenomic distribution
          for (var field in mongo.bins) {
            solr[field + '__bin'] = mongo.bins[field];
          }

          // check if this gene has expression data
          if (expressionData[mongo._id]) {
            solr.capabilities.push('expression');
            // make some __expr fields
            _.forEach(expressionData[mongo._id], function(samples, exp_id) {
              if (exp_id !== '_id') {
                let solrKey = exp_id.replace(/-/g,'_');
                _.forEach(samples, function(sample) {
                  if (sample.value) {
                    solr[`${solrKey}_${sample.group}__expr`] = sample.value;
                  }
                  if (sample.l2fc) {
                    solr[`${solrKey}_${sample.group}_l2fc_attr_f`] = sample.l2fc;
                  }
                  if (sample.p_value) {
                    solr[`${solrKey}_${sample.group}_pval_attr_f`] = sample.p_value;
                  }
                });
              }
            });
          }

          // homology fields
          if (mongo.homology) {
            solr.capabilities.push('homology');
            p++;
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
                  if (solr.closest_rep_name) uniq[solr.closest_rep_name.toLowerCase()] = solr.closest_rep_name;
                }

                if (mhgr.hasOwnProperty('model')) {
                  var rep = get_rep(mhgr.model);
                  for (var f in rep) {
                    solr['model_rep_'+f] = rep[f];
                  }
                  if (solr.model_rep_name) uniq[solr.model_rep_name.toLowerCase()] = solr.model_rep_name;
                }
              }
            }
        
            // The crop cultivar tree sets. Field names match the mongo homology keys exactly so
            // the two sides cannot drift. A gene may carry several at once -- a wheat cultivar
            // gene is in both the default plants tree and the wheat_cultivars tree -- so these are
            // additive, and an absent field simply means the gene is not in that set.
            ['barley', 'oat', 'wheat', 'rice'].forEach(function (set) {
              if (mongo.homology[set]) {
                solr[set] = mongo.homology[set].id;
                solr[set + '_root_taxon_id'] = mongo.homology[set].root_taxon_id;
              }
            });

            // Left as it was: pan-homology trees are not loaded, so homology.pan_tree is never
            // written. Kept rather than removed because it predates this change and was already
            // unpopulated, so deleting it would be churn for no gain.
            if (mongo.homology.pan_tree) {
              solr.pan_tree = mongo.homology.pan_tree.id;
              solr.pan_tree_root_taxon_id = mongo.homology.pan_tree.root_taxon_id;
            }

            if (mongo.homology.homologous_genes) {
              solr.homology__all_orthologs = [mongo._id];
              solr.homology__all_homeologs = [mongo._id];
              for (htype in mongo.homology.homologous_genes) {
                if (Array.isArray(mongo.homology.homologous_genes[htype])) {
                  solr['homology__'+htype] = mongo.homology.homologous_genes[htype];
                  if (htype.match(/ortholog/)) {
                    mongo.homology.homologous_genes[htype].forEach(function(o) {
                      solr.homology__all_orthologs.push(o);
                    });
                  }
                  if (htype.match(/homeolog/)) {
                    mongo.homology.homologous_genes[htype].forEach(function(h) {
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
            let pubs=[];
            let have_xrefs=false;
            mongo.xrefs.forEach(function(xref) {
              if (xref.db === 'PUBMED') {
                pubs.push(xref.ids[0])
              }
              else if (!mongo.annotations.hasOwnProperty(xref.db)) {
                solr[xref.db + '__xrefs'] = xref.ids;
                have_xrefs=true;
                xref.ids.forEach(function(id) {
                  if (_.isString(id)) {
                    var lc = id.toLowerCase();
                    if (!uniq.hasOwnProperty(lc)) {
                      uniq[lc]=id;
                    }
                  }
                })
              }
            });
            if (pubs.length > 0) {
              solr.capabilities.push('pubs');
              solr['PUBMED__xrefs'] = pubs;
            }
            if (have_xrefs) {
              solr.capabilities.push('xrefs');
            }
          }
          // get we don't want id in _terms
          delete uniq[lcId];
          solr._terms = Object.keys(uniq).map(function(k){return uniq[k]});
      
          // add ancestors fields from the annotations section
          // and text of annotations (except taxonomy)
          for (var f in mongo.annotations) {
            if (mongo.annotations[f] && (mongo.annotations[f].ancestors || mongo.annotations[f].entries)) {
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
                  if (f !== "taxonomy") {
                    if (e.id) solr.annotations.push(e.id);
                    if (e.name) solr.annotations.push(e.name);
                    if (e.description) solr.annotations.push(e.description);
                    if (e.def) solr.annotations.push(e.def);
                  }
                });
              }
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
});
