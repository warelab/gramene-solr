#!/usr/bin/env node
var collections = require('gramene-mongodb-config');
var _ = require('lodash');
var fs = require('fs');

// OPTIONAL SIDE-OUTPUT: term -> taxon_id pairs, for suggestions/genes.js.
//
// genes.js used to recover _terms by faceting the genes SOLR core -- a 588s pivot plus one
// unbounded facet per taxon (1083s for 267), to get back a field THIS script derived from mongo in
// the first place. That round trip was also the only thing that ever failed: on 2026-10-18 the
// gene_tree pivot's q-io promise never settled, the event loop drained, and the process exited 0
// leaving a 4.49 GB genes.json with no closing ']' and no error -- Q keeps rejections to itself.
//
// Emitting the pairs HERE rather than re-deriving them in genes.js is deliberate. _terms is not a
// simple union: it folds in name, id, synonyms, alt_id, description (unless "unknown") and the ids
// of every xref whose db is absent from mongo.annotations, then removes the gene's own id -- and
// name/synonyms/description each have their own upstream fixups. Re-implementing that elsewhere
// would drift, and suggestions would start offering terms the genes core cannot match.
//
// Off unless SOLR_TERMS_TSV is set, so the normal stage-60 run is byte-for-byte unchanged.
// Format: JSON.stringify(term) TAB taxon_id  -- the term is JSON-encoded because descriptions and
// xref ids can contain tabs, quotes and non-ASCII, and because identical terms must encode
// identically for a later sort to group them.
var TERMS_TSV = process.env.SOLR_TERMS_TSV || '';
var terms_fd = TERMS_TSV ? fs.openSync(TERMS_TSV, 'w') : -1;
var terms_buf = [], terms_len = 0, terms_rows = 0;
function terms_emit(termList, taxon_id) {
  if (terms_fd < 0) { return; }
  for (var i = 0; i < termList.length; i++) {
    var t = termList[i];
    if (typeof t !== 'string' || !t.length) { continue; }
    terms_buf.push(JSON.stringify(t)); terms_buf.push('\t'); terms_buf.push(String(taxon_id)); terms_buf.push('\n');
    terms_len += 1; terms_rows++;
    if (terms_len >= 200000) { fs.writeSync(terms_fd, terms_buf.join('')); terms_buf = []; terms_len = 0; }
  }
}
function terms_close() {
  if (terms_fd < 0) { return; }
  if (terms_len) { fs.writeSync(terms_fd, terms_buf.join('')); terms_buf = []; terms_len = 0; }
  fs.closeSync(terms_fd);
  console.error('SOLR_TERMS_TSV: wrote ' + terms_rows + ' term/taxon rows to ' + TERMS_TSV);
}
process.on('exit', terms_close);

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

collections.genes.mongoCollection().then(function(genesCol) {
  const cursor = genesCol.aggregate([
    { $sort: { species_idx: 1, db_type: 1, gene_idx: 1 } },

    // Join genes._id to expression._id
    {
      $lookup: {
        from: "expression",
        localField: "_id",
        foreignField: "_id",
        as: "expression"
      }
    },

    // Convert array -> single doc (or leave null if missing)
    { $set: { expression: { $first: "$expression" } } },

    // Optional: if you don't want the _id inside the embedded expression
    // { $set: { "expression._id": "$$REMOVE" } },
  ], { allowDiskUse: true });

  var n=0;
  var p_padding = 100;
  var p= 0 - p_padding;
  var terminator = {}; // key is lowercase version, value is original term
  var current_region = "not likely";
  
  cursor.forEach(
    function(mongo) {
      var location = mongo.location;
      if (location.region !== current_region) {
        p += p_padding;
        current_region = location.region;
      }
      if (!mongo.description) {
        mongo.description = 'unknown';
      }

      var solr = {
        // required data for list view
        id          : mongo._id,
        name        : mongo.name || mongo._id,
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
        saved_search : [] // initialize empty list. TODO: at release time, update this field based on previous core
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

      // alternate ids (Sobic./Sb.../GRMZM/AC/pan.../Vitvi.../LOC_Os…) get their own
      // indexed + facetable `alt_id` field (powers the "Alternate IDs" suggestion's exact
      // fq=alt_id resolution) and are also folded into `_terms`. Free-text/API search
      // reaches them via the schema copyField alt_id -> text (mirroring synonyms), since
      // they no longer live in gene.synonyms.
      if (Array.isArray(mongo.alt_id) && mongo.alt_id.length) {
        var altUniq = {};
        mongo.alt_id.forEach(function(a) {
          if (typeof a === 'string' && a.length) {
            var alc = a.toLowerCase();
            if (!altUniq.hasOwnProperty(alc)) altUniq[alc] = a;
            if (!uniq.hasOwnProperty(alc)) uniq[alc] = a; // include in _terms
          }
        });
        var altVals = Object.keys(altUniq).map(function(k){return altUniq[k];});
        if (altVals.length) solr.alt_id = altVals;
      }


      // add the description
      if (solr.description !== "unknown") {
        uniq[solr.description.toLowerCase()] = solr.description;
      }
  
      if (mongo.summary) {
        solr.summary = mongo.summary;
      }

      // facet counting on bin fields drives the taxagenomic distribution
      for (var field in mongo.bins) {
        solr[field + '__bin'] = mongo.bins[field];
      }
      if (mongo && mongo.expression) {
        const exprDoc = mongo.expression;

        // ----- your merge code skeleton -----
        const expressedInStudy = new Set();
        for (const [exp_id, samples] of Object.entries(exprDoc)) {
          if (exp_id === "_id") continue;
          const solrKey = exp_id.replace(/-/g, "_");

          for (const sample of samples) {
            if (sample.value != null) {
              solr[`${solrKey}_${sample.group}__expr`] = sample.value;
              if (sample.value >= 0.5) expressedInStudy.add(exp_id);
            }
            if (sample.l2fc != null) {
              solr[`${solrKey}_${sample.group}_l2fc_attr_f`] = sample.l2fc;
            }
            if (sample.p_value != null) {
              solr[`${solrKey}_${sample.group}_pval_attr_f`] = sample.p_value;
              if (sample.p_value < 0.05 && Math.abs(sample.l2fc ?? 0) >= 1.0) {
                expressedInStudy.add(exp_id);
              }
            }
          }
        }
        if (expressedInStudy.size > 0) {
          solr.capabilities ??= [];
          solr.capabilities.push("expression");
          solr.expressed_in_gxa_attr_ss = Array.from(expressedInStudy);
        }
        // ------------------------------------
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
              // if (solr.closest_rep_name) uniq[solr.closest_rep_name.toLowerCase()] = solr.closest_rep_name;
            }

            if (mhgr.hasOwnProperty('model')) {
              var rep = get_rep(mhgr.model);
              for (var f in rep) {
                solr['model_rep_'+f] = rep[f];
              }
              // if (solr.model_rep_name) uniq[solr.model_rep_name.toLowerCase()] = solr.model_rep_name;
            }
          }
        }
    
        // The crop cultivar tree sets. Field names match the mongo homology keys exactly so the
        // two sides cannot drift, and the genes schema already declares all eight as stored, so
        // they reach the API with no schema change. A gene may carry several at once -- a wheat
        // cultivar gene is in both the default plants tree and the wheat_cultivars tree -- so
        // these are additive, and an absent field simply means the gene is not in that set. That
        // is what lets a client offer the user a choice of which tree to load.
        //
        // This block was added to mongo2solr.js in 49679f3, but stage 60 runs THIS file, so until
        // now the cultivar tree ids never reached solr and the API exposed gene_tree only.
        var inCultivarTree = false;
        ['barley', 'oat', 'wheat', 'rice'].forEach(function (set) {
          if (mongo.homology[set]) {
            solr[set] = mongo.homology[set].id;
            solr[set + '_root_taxon_id'] = mongo.homology[set].root_taxon_id;
            inCultivarTree = true;
          }
        });

        // Which tree sets is this gene in? capabilities is what tells a client what to offer for
        // a gene -- the 'homology' value pushed above is what enables its Homology tab -- and the
        // tab cannot offer a choice of tree without knowing which sets exist for this gene. The
        // <set> id fields are not in the search response's default fl, so say it in capabilities:
        // one value per set, plus a cultivar_tree roll-up for the plain "is this gene in any
        // cultivar set tree" question, which is otherwise four checks on the client.
        //
        // The default set is included so a chooser can read all five from one place. familyRoot
        // already happens to mark the same genes, but only as an xref side-effect (ancestor_adder
        // adds an NCBITaxon xref for the tree's root taxon); its name says nothing about tree
        // membership and it would silently stop tracking it if that xref ever changed.
        if (mongo.homology.gene_tree) { solr.capabilities.push('gene_tree'); }
        ['barley', 'oat', 'wheat', 'rice'].forEach(function (set) {
          if (mongo.homology[set]) { solr.capabilities.push(set + '_tree'); }
        });
        if (inCultivarTree) { solr.capabilities.push('cultivar_tree'); }

        // Supertrees, per set. decorate writes homology.<set>.supertree (search/supertrees.js),
        // so this reads the same shape the tree ids use rather than the single flat
        // homology.supertree that release 69 produced -- a gene can be in a default supertree AND
        // a cultivar one, and one field cannot carry both.
        //
        // The default set keeps r69's field name, supertree_attr_s, so anything already reading it
        // keeps working; the cultivar sets get <set>_supertree_attr_s. No schema change is needed:
        // the genes managed-schema's *_attr_s dynamicField covers all five.
        ['gene_tree', 'barley', 'oat', 'wheat', 'rice'].forEach(function (set) {
          var h = mongo.homology[set];
          if (h && h.supertree) {
            solr[(set === 'gene_tree' ? 'supertree' : set + '_supertree') + '_attr_s'] = h.supertree;
          }
        });
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
      solr.canonical_transcript__attr_s = ct;
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
          else if (xref.db === 'PanOryza') {
            solr['PanOryza__xrefs'] = xref.ids;
            solr.capabilities.push('PanOryza');
          }
          else if (xref.db === 'SB_PanGene') {
            solr['SB_PanGene__xrefs'] = xref.ids;
            solr.capabilities.push('SB_PanGene');
          }
          else if (xref.db === 'Grassius') {
            solr['Grassius__xrefs'] = xref.ids;
            solr.capabilities.push('Grassius');
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
      // side channel for suggestions/genes.js -- same array, same moment, so it cannot drift
      terms_emit(solr._terms, solr.taxon_id);
  
      // add ancestors fields from the annotations section
      // and text of annotations (except taxonomy)
      for (var f in mongo.annotations) {
        if (mongo.annotations[f] && (mongo.annotations[f].ancestors || mongo.annotations[f].entries)) {
          solr.capabilities.push(f);
          solrField = f + '__ancestors';
          if (mongo.annotations[f].ancestors) {
            solr[solrField] = mongo.annotations[f].ancestors.filter(a => typeof a === 'number');
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
    },
    function(err) {
      if (err) throw err;
      console.log(']');
      collections.closeMongoDatabase();
    }
  );
});
