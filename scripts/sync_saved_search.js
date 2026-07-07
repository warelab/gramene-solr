#!/usr/bin/env node
'use strict';
// sync_saved_search.js — propagate user gene-list `saved_search` hashes from an OLD gramene genes
// core to a NEW one. Fulfills the "// TODO: at release time, update this field based on previous
// core" left at gramene-solr/genes/mongo2solr.join.js:84.
//
// A user gene list tags its member genes with a murmur3-x86-32 `hash` in the multi-valued
// `saved_search` field of the genes core (written by /gene_lists/validate at save time). The core is
// rebuilt per release, so those tags do not carry over; the list metadata lives in the shared
// userData1.genelists collection. This tool copies the still-valid tags forward: for every gene in
// the OLD core it keeps only the saved_search hashes that STILL exist in userData1.genelists (a
// force-deleted / 30-day-purged list is gone from that collection, so its hash is dropped) and
// atomic-`add-distinct`s them onto the same gene in the NEW core. Additive (never removes) and
// idempotent (add-distinct). Trashed-but-restorable lists are kept (their hashes are still present).
//
//   config via env: MONGO_URL, SOURCE, TARGET, READ_BATCH, EXIST_BATCH, POST_BATCH
//   flags:          --dry-run   (compute + report, no writes)
//                   --skip-exist (skip the existence check — trust that every id exists in the new core)
//
//   node scripts/sync_saved_search.js --dry-run
//   node scripts/sync_saved_search.js
//   SOURCE=https://data.sorghumbase.org/sorghum_v10/search node scripts/sync_saved_search.js
//
// Uses Node's global fetch (v18+) + the mongodb driver (already a gramene-solr dependency); no new deps.

const MONGO_URL = process.env.MONGO_URL || 'mongodb://localhost:27017/userData1';
const rawSource = process.env.SOURCE   || 'http://localhost:8983/solr/sorghum_genes10';
const rawTarget = process.env.TARGET   || 'http://localhost:8983/solr/sorghum_genes11';

const argv = process.argv.slice(2);
const DRY_RUN    = argv.includes('--dry-run')    || process.env.DRY_RUN === '1';
const SKIP_EXIST = argv.includes('--skip-exist') || process.env.SKIP_EXIST_FILTER === '1';

const READ_BATCH  = parseInt(process.env.READ_BATCH  || '1000', 10);
const EXIST_BATCH = parseInt(process.env.EXIST_BATCH || '5000', 10);
const POST_BATCH  = parseInt(process.env.POST_BATCH  || '5000', 10);
const SEP = '|';   // {!terms} separator; gene ids never contain '|'

// OLD core read URL: append /select unless it already ends in /select or /search (so the public
// https://.../sorghum_v10/search endpoint works as-is). NEW core is a direct solr core (needs /update).
const SOURCE_SELECT = /\/(select|search)$/.test(rawSource) ? rawSource : rawSource.replace(/\/+$/, '') + '/select';
const TARGET_BASE   = rawTarget.replace(/\/+$/, '');
const TARGET_SELECT = TARGET_BASE + '/select';
const TARGET_UPDATE = TARGET_BASE + '/update';

async function getJson(url, params) {
  const r = await fetch(url + '?' + new URLSearchParams(params).toString());
  if (!r.ok) throw new Error('GET ' + url + ' failed ' + r.status + ' ' + (await r.text()).slice(0, 300));
  return r.json();
}
async function postSelect(url, params) {
  const r = await fetch(url, { method: 'POST', headers: { 'Content-Type': 'application/x-www-form-urlencoded' },
                               body: new URLSearchParams(params).toString() });
  if (!r.ok) throw new Error('POST select ' + url + ' failed ' + r.status + ' ' + (await r.text()).slice(0, 300));
  return r.json();
}
async function postUpdate(body, commit) {
  const r = await fetch(TARGET_UPDATE + (commit ? '?commit=true' : ''),
                        { method: 'POST', headers: { 'Content-Type': 'application/json' },
                          body: typeof body === 'string' ? body : JSON.stringify(body) });
  if (!r.ok) throw new Error('update failed ' + r.status + ' ' + (await r.text()).slice(0, 300));
  return r.json();
}

// 1. hashes still present in userData1.genelists (active + trashed/restorable lists)
async function loadHashSet() {
  const { MongoClient } = require('mongodb');
  const client = new MongoClient(MONGO_URL);
  await client.connect();
  try {
    const hashes = await client.db().collection('genelists').distinct('hash');
    return new Set(hashes.map(Number).filter(h => !Number.isNaN(h)));
  } finally { await client.close(); }
}

// 2. scan the OLD core, keeping per gene only the saved_search hashes that are in the map
async function scanSource(hashSet) {
  const contributing = new Map();   // id -> [in-map hashes]
  let start = 0, numFound = Infinity, scanned = 0;
  while (start < numFound) {
    const j = await getJson(SOURCE_SELECT,
      { q: 'saved_search:*', fl: 'id,saved_search', rows: String(READ_BATCH), start: String(start), wt: 'json' });
    numFound = j.response.numFound;
    const docs = j.response.docs || [];
    if (!docs.length) break;
    for (const d of docs) {
      scanned++;
      const vals = Array.isArray(d.saved_search) ? d.saved_search : (d.saved_search != null ? [d.saved_search] : []);
      const kept = [...new Set(vals.map(Number).filter(h => !Number.isNaN(h) && hashSet.has(h)))];
      if (kept.length) contributing.set(d.id, kept);
    }
    start += docs.length;
    console.error(`  scanned ${Math.min(start, numFound)}/${numFound}`);
  }
  return { contributing, scanned };
}

// 3. which of these ids actually exist in the NEW core (an atomic update on a missing id would
//    create an orphan partial doc). {!terms} in the POST body scales past URL length limits.
async function existingInTarget(ids) {
  if (SKIP_EXIST) return new Set(ids);
  const present = new Set();
  for (let i = 0; i < ids.length; i += EXIST_BATCH) {
    const batch = ids.slice(i, i + EXIST_BATCH);
    const j = await postSelect(TARGET_SELECT,
      { q: '{!terms f=id separator="' + SEP + '"}' + batch.join(SEP), fl: 'id', rows: String(batch.length), wt: 'json' });
    for (const d of j.response.docs) present.add(d.id);
  }
  return present;
}

(async () => {
  console.error(`MONGO_URL = ${MONGO_URL}`);
  console.error(`SOURCE    = ${SOURCE_SELECT}`);
  console.error(`TARGET    = ${TARGET_UPDATE}${DRY_RUN ? '   (DRY RUN — no writes)' : ''}`);

  const hashSet = await loadHashSet();
  console.error(`hash map: ${hashSet.size} distinct saved-list hashes in genelists`);
  if (!hashSet.size) { console.error('empty hash map — nothing to sync'); process.exit(0); }

  const { contributing, scanned } = await scanSource(hashSet);
  console.error(`old core: scanned ${scanned} saved_search docs; ${contributing.size} carry an in-map hash`);
  if (!contributing.size) { console.error('no genes to update'); process.exit(0); }

  const ids = [...contributing.keys()];
  const present = await existingInTarget(ids);
  console.error(`new core: ${present.size}/${ids.length} ids present (${ids.length - present.size} absent — skipped)`);

  const docs = [];
  for (const [id, hashes] of contributing) if (present.has(id)) docs.push({ id, saved_search: { 'add-distinct': hashes } });
  console.error(`atomic add-distinct updates to apply: ${docs.length}`);

  if (DRY_RUN) { console.error('DRY RUN — no writes made'); process.exit(0); }

  for (let i = 0; i < docs.length; i += POST_BATCH) {
    await postUpdate(docs.slice(i, i + POST_BATCH), false);
    console.error(`  update ${Math.min(i + POST_BATCH, docs.length)}/${docs.length}`);
  }
  await postUpdate('{"commit":{}}', false);
  console.error(`DONE — ${docs.length} genes updated + committed on the new core`);
})().catch(e => { console.error('FAILED:', e && e.message || e); process.exit(1); });
