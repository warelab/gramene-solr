#!/usr/bin/env node

const collections = require('gramene-mongodb-config');

const genesURL = process.argv[2];

async function getTaxa(query) {
  const url = `${genesURL}/query?q=${encodeURIComponent(query)}&rows=0&facet=true&facet.field=taxon_id&facet.limit=-1&facet.mincount=1&json.nl=map`;
  console.error("getTaxa", query, url)
  const res = await fetch(url);
  const json = await res.json();
  const taxTally = json.facet_counts.facet_fields.taxon_id;
  return {
    ids: Object.keys(taxTally).map(Number),
    counts: Object.values(taxTally),
  };
}

(async () => {
  const geneCol = await collections.genes.mongoCollection();
  console.error("fetch synonyms");

  const cursor = geneCol.find(
    { taxon_id: { $exists: true } }
  );
  cursor.project({ _id: 1, taxon_id: 1, db_type: 1, synonyms: 1 });
  const docs = await cursor.toArray();
  
  console.error("got synonyms from", docs.length, "genes docs");

  const uniqueId = new Map();
  const otherfeaturesId = new Map();
  const uniqueTaxa = new Set();
  const originalCase = new Map();
  const synOf = new Map();

  for (const d of docs) {
    const upperId = d._id.toUpperCase();
    if (d.db_type !== 'core') {
      otherfeaturesId.set(upperId, d.taxon_id);
    }
    uniqueTaxa.add(d.taxon_id);
    uniqueId.set(upperId, d.taxon_id);
    originalCase.set(upperId, d._id);
    if (d.synonyms) {
      d.synonyms.filter(syn => syn.length >= 10).forEach(syn => {
        const upperSyn = syn.toUpperCase();
        uniqueId.set(upperSyn, d.taxon_id);
        originalCase.set(upperSyn, syn);
        synOf.set(upperSyn, d._id);
      });
    }
  }

  console.error("loaded synonyms into LUT");
  const taxaArray = [...uniqueTaxa];
  console.error("uniqueTaxa: ", taxaArray);


  await collections.closeMongoDatabase();

  const url = `${genesURL}/query?rows=0&facet=true&facet.limit=-1&json.nl=map&facet.pivot.mincount=2&facet.pivot=_terms,taxon_id`;
  console.error(url);
  const term_freq = new Map();
  const taxa_lut = new Map();
  try {
    const res = await fetch(url);
    console.error("got res");
    const data = await res.json();
    console.error("parsed")

    data.facet_counts.facet_pivot['_terms,taxon_id'].forEach(d => {
      const term = d.value; // don't .toUpperCase() here
      if (!uniqueId.has(term.toUpperCase())) {
        if (term_freq.has(term)) {
          term_freq.set(term, term_freq.get(term) + d.count);
        } else {
          term_freq.set(term, d.count);
          taxa_lut.set(term, new Map());
        }

        if (d.pivot) {
          d.pivot.forEach(p => {
            const tlut = taxa_lut.get(term);
            if (tlut.has(p.value)) {
              tlut.set(p.value, tlut.get(p.value) + p.count);
            } else {
              tlut.set(p.value, p.count);
            }
          });
        }
      }
    });
    console.error("loaded term_freq and taxa_lut maps");
  } catch (err) {
    console.error('Error fetching Solr pivot data:', err)
  }
  let n = 0;
  // find all the terms in each species
  for( let taxon_id of uniqueTaxa) {
    var surl = genesURL + '/query?rows=0&facet=true&facet.field=_terms&facet.limit=-1&json.nl=map&facet.mincount=1&q=taxon_id:' + taxon_id;
    console.error(surl);
    try {
      const sres = await fetch(surl);
      const sdata = await sres.json();

      var term_tally = sdata.facet_counts.facet_fields._terms;
      console.error("got term_tally for taxon_id",taxon_id)
      for (var t in term_tally) {
        var term = t; // .toUpperCase(); // don't do this it leads to problems
        if (!uniqueId.has(term.toUpperCase())) {
          if (term_freq.has(term)) {
            const tlut = taxa_lut.get(term);
            if (tlut.has(taxon_id)) {
              tlut.set(taxon_id, tlut.get(taxon_id) + term_tally[t]);
            } else {
              tlut.set(taxon_id, term_tally[t]);
            }
          } else if (/[A-Za-z]/.test(term)) { // this is a unique term that contains a word character
            var solr = {
              category: 'Genes',
              subcategory: 'term',
              id: '_term_' + ++n,
              display_name: term,
              name: term,
              fq_field: '_terms',
              fq_value: term,
              taxon_id: [+taxon_id],
              taxon_freq: [1],
              num_genes: 1,
              relevance: 1
            };
            if (n === 1) console.log('[');
            else console.log(',');
            console.log(JSON.stringify(solr));
          }
        }
      }
    } catch (err) {
      console.error("Error fetching pivot data for taxon_id",taxon_id, err)
    }
  };

  console.error("finished building taxa_lut. n=", n);
  
  // output suggestions for Genes
  for (var [term, tf] of term_freq) {
    var taxa = {
      ids: [],
      counts: []
    };
    for (var [taxon_id,count] of taxa_lut.get(term)) {
      taxa.ids.push(+taxon_id);
      taxa.counts.push(count);
    }
    var solr = {
      category: 'Genes',
      subcategory: 'term',
      id: '_term_' + ++n,
      display_name: term,
      name: term,
      fq_field: '_terms',
      fq_value: term,
      taxon_id: taxa.ids,
      taxon_freq: taxa.counts,
      num_genes: tf,
      relevance: tf < 10 ? 1.1 : 1
    };
    if (n === 1) console.log('[');
    else console.log(',');
    console.log(JSON.stringify(solr));
  }

  // Gene trees
  const gturl = genesURL +
            '/query?rows=0&facet=true&facet.limit=-1&facet.mincount=1&json.nl=map' +
            '&facet.pivot=gene_tree,taxon_id&q=gene_tree:*';
  console.error(gturl);
  try {
    const gtres = await fetch(gturl);
    const gtdata = await gtres.json();
    gtdata.facet_counts.facet_pivot['gene_tree,taxon_id'].forEach(function(d) {
      var taxa = {
        ids: [],
        counts: []
      };
      d.pivot.forEach(function(p) {
        taxa.ids.push(p.value);
        taxa.counts.push(p.count);
      });
      var solr = {
        category: 'Gene trees',
        id: '_term_' + ++n,
        taxon_id: taxa.ids,
        taxon_freq: taxa.counts,
        display_name: d.value,
        name: d.value,
        fq_field: 'gene_tree',
        fq_value: d.value,
        num_genes: d.count,
        relevance: 1
      }
      console.log(',');
      console.log(JSON.stringify(solr));
    })
  } catch (err) {
    console.error("Error fetching pivot data for gene trees", err)
  }
  console.error("finished gene tree suggestions")
  
  // SB PanGene xrefs
  const pgurl = genesURL + '/query?rows=0&facet=true&facet.limit=-1&json.nl=map&facet.mincount=1&facet.pivot=SB_PanGene__xrefs,taxon_id&q=SB_PanGene__xrefs:*';
  console.error(pgurl);
  try {
    const pgres = await fetch(pgurl);
    const pgdata = await pgres.json();
    if (pgdata) {
      pgdata.facet_counts.facet_pivot['SB_PanGene__xrefs,taxon_id'].forEach(function(d) {
        var taxa = {
          ids: [],
          counts: []
        };
        d.pivot.forEach(function(p) {
          taxa.ids.push(p.value);
          taxa.counts.push(p.count);
        });
        var solr = {
          category: 'Sorghum pan genes',
          id: '_term_' + ++n,
          taxon_id: taxa.ids,
          taxon_freq: taxa.counts,
          display_name: d.value,
          name: d.value,
          fq_field: 'SB_PanGene__xrefs',
          fq_value: d.value,
          num_genes: d.count,
          relevance: 1
        }
        console.log(',');
        console.log(JSON.stringify(solr));
      })
    }
  } catch (err) {
    console.error("Error fetching pivot data for SB_PanGene__xrefs", err);
  }
  
  // biotype
  const biotypeURL = genesURL + '/query?rows=0&facet=true&facet.limit=-1&json.nl=map&facet.mincount=1&facet.field=biotype&q=biotype:*';
  console.error(biotypeURL);
  try {
    const biotypeRes = await fetch(biotypeURL);
    const biotypeData = await biotypeRes.json();
    var biotype_count = biotypeData.facet_counts.facet_fields.biotype;
    for (var biotype in biotype_count) {
      var taxa = await getTaxa('biotype:' + biotype);
      console.error(biotype, taxa);
      var solr = {
        category: 'Biotypes',
        id: '_term_' + ++n,
        display_name: biotype,
        name: biotype,
        fq_field: 'biotype',
        fq_value: biotype,
        taxon_id: taxa.ids,
        taxon_freq: taxa.counts,
        num_genes: biotype_count[biotype],
        relevance: 1
      };
      console.log(',');
      console.log(JSON.stringify(solr));
    }
  } catch (err) {
    console.error("Error fetching biotype counts", err);
  }
  
  // grassius
  const grassiusurl = genesURL + '/query?rows=0&facet=true&facet.limit=-1&json.nl=map&facet.mincount=1&facet.pivot=Grassius__xrefs,taxon_id&q=Grassius__xrefs:*';
  console.error(grassiusurl);
  try {
    const grassres = await fetch(grassiusurl);
    const grassdata = await grassres.json();
    grassdata.facet_counts.facet_pivot['Grassius__xrefs,taxon_id'].forEach(function(d) {
      var taxa = {
        ids: [],
        counts: []
      };
      d.pivot.forEach(function(p) {
        taxa.ids.push(p.value);
        taxa.counts.push(p.count);
      });
      var solr = {
        category: 'Grassius',
        id: '_term_' + ++n,
        taxon_id: taxa.ids,
        taxon_freq: taxa.counts,
        display_name: d.value,
        name: d.value,
        fq_field: 'Grassius__xrefs',
        fq_value: d.value,
        num_genes: d.count,
        relevance: 1
      }
      console.log(',');
      console.log(JSON.stringify(solr));
    })
  } catch (err) {
    console.error("Error fetching pivot data for SB_PanGene__xrefs", err);
  }
  
  // output unique IDs
  console.error('output uniqueIds');
  for (let [uid, val] of uniqueId) {
    console.log(',');
    const ouid = originalCase.get(uid);
    console.log(JSON.stringify({
      category: 'Genes',
      subcategory: 'id',
      fq_field: 'id',
      fq_value: synOf.get(uid) || ouid,
      id: ouid,
      display_name: ouid,
      num_genes: 1,
      relevance: 1.2,
      taxon_id: val,
      taxon_freq: [1]
    }));
  }
  console.log(']');
  
})();
