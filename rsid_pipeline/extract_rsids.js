#!/usr/bin/env node
/**
 * extract_rsids.js — rsIDs overlapping a genome's gene models, for the rsid__attr_ss attribute.
 *
 *   node extract_rsids.js --genome sorghum_bicolor --vcf /path/sorghum_bicolor.vcf.gz > rows.tsv
 *
 * Emits attribute-table rows (no header; the driver writes that):
 *   <gene id>\trsid\t<rsID,rsID,...>
 *
 * A variant counts for a gene when it falls inside the gene model, or within the UPSTREAM flank
 * (strand-aware: below start on +, above end on -).
 * Intronic variants further than INTRON_MAX from a canonical exon are dropped: a deep intronic
 * variant says little about the gene and only adds noise to the field.
 *
 * Coordinates go through gramene-gene-positions, never through local arithmetic.
 * gene_structure.exons are gene-RELATIVE offsets in TRANSCRIPTION order, so on the minus strand
 * genomic coordinates DECREASE as gene coordinates increase. remap() knows that; hand-rolled
 * arithmetic silently gets every minus-strand gene wrong.
 *
 * Genotypes are never materialised: bcftools emits CHROM/POS/ID only. The same rice callset on this
 * host is 194 GB with genotypes and 67 KB without, and none of it is relevant here.
 */
'use strict';
const { spawn } = require('child_process');
const readline = require('readline');
const positions = require('gramene-gene-positions');
const collections = require('gramene-mongodb-config');

function arg(name, dflt) {
  const i = process.argv.indexOf('--' + name);
  return i !== -1 && process.argv[i + 1] ? process.argv[i + 1] : dflt;
}
const GENOME     = arg('genome');
const VCF        = arg('vcf');
// Flank is UPSTREAM-ONLY by default and therefore strand-dependent: on the + strand it extends
// below location.start, on the - strand above location.end. A symmetric window is not the same
// thing and would silently include downstream variants on half the genes.
const FLANK_UP   = +arg('flank-up', process.env.RSID_FLANK_UP || 200);
const FLANK_DOWN = +arg('flank-down', process.env.RSID_FLANK_DOWN || 0);
const INTRON_MAX = +arg('intron-max', process.env.RSID_INTRON_MAX_DIST || 10);
const MAX_PER_GENE = +arg('max-per-gene', process.env.RSID_MAX_PER_GENE || 5000);
const MIN_CHROM_FRAC = +arg('min-chrom-frac', process.env.RSID_MIN_CHROM_FRAC || 0.5);
const CHROM_MAP_FILE = arg('chrom-map', process.env.RSID_CHROM_MAP || '');
const CLASS_OUT = arg('emit-classes', process.env.RSID_EMIT_CLASSES || '');   // rsid\tgene\tclass
// MODE:
//   gene (default) - window is the gene model +/- flank; keep exonic, near-intronic, flank
//   cds            - window is the CDS +/- CDS_UP/CDS_DOWN; keep only coding-region variants,
//                    splice-region variants inside the coding region, and the flanks. UTR variants
//                    are dropped, which is the point: some models here carry a 50 kb 5'UTR around a
//                    360 bp CDS, and the gene-model window counts all of it.
const MODE     = arg('mode', process.env.RSID_MODE || 'gene');
const CDS_UP   = +arg('cds-up', process.env.RSID_CDS_UP || 1000);     // bp before the start codon
const CDS_DOWN = +arg('cds-down', process.env.RSID_CDS_DOWN || 500);  // bp after the stop codon
const SPLICE   = +arg('splice', process.env.RSID_SPLICE || 10);       // nt at each intron end
const NC_UP    = +arg('nc-up', process.env.RSID_NC_UP || 200);        // non-coding: bp before transcript start
const SPLICE_PTV = +arg('splice-ptv', process.env.RSID_SPLICE_PTV || 2); // nt each side counting as a splice-site PTV
const FASTA_DIR  = arg('fasta-dir', process.env.RSID_CDS_FASTA_DIR || '/scratch/olson/fasta');
const NO_CONSEQ  = process.argv.includes('--no-consequence');
const MAX_REF_MISMATCH = +arg('max-ref-mismatch', process.env.RSID_MAX_REF_MISMATCH || 0.001);

if (!GENOME || !VCF) {
  console.error('usage: extract_rsids.js --genome <system_name> --vcf <file.vcf[.gz]>');
  process.exit(2);
}
const die = m => { console.error('FATAL [' + GENOME + ']: ' + m); process.exit(1); };
const log = m => console.error('  ' + m);

// optional CHROM normalisation, e.g. "chr1<TAB>1" per line. Explicit map, never guessing.
const chromMap = {};
if (CHROM_MAP_FILE) {
  require('fs').readFileSync(CHROM_MAP_FILE, 'utf8').split('\n').forEach(l => {
    const [a, b] = l.split('\t');
    if (a && b) chromMap[a.trim()] = b.trim();
  });
  log('loaded ' + Object.keys(chromMap).length + ' CHROM normalisations');
}
const normChrom = c => (chromMap[c] !== undefined ? chromMap[c] : c);

// ── interval index: per region, intervals sorted by start with a running max-end, so containment
// works by binary search and does NOT assume the VCF is coordinate-sorted.
function buildIndex(genes) {
  const byRegion = {};
  genes.forEach(g => {
    const r = String(g.location.region);
    const plus = g.location.strand !== -1;
    // In cds mode the precise window is computed in pass 2 (it needs gene_structure); here we cast
    // a superset net -- the CDS lies inside the gene, so gene +/- max(CDS_UP,CDS_DOWN) contains it.
    const pad = (MODE === 'cds') ? Math.max(CDS_UP, CDS_DOWN) : 0;
    const up = (MODE === 'cds') ? pad : (plus ? FLANK_UP : FLANK_DOWN);
    const dn = (MODE === 'cds') ? pad : (plus ? FLANK_DOWN : FLANK_UP);
    (byRegion[r] = byRegion[r] || []).push({
      id: g._id,
      lo: Math.max(1, g.location.start - up),
      hi: g.location.end + dn,
      gs: g.location.start, ge: g.location.end
    });
  });
  Object.keys(byRegion).forEach(r => {
    const a = byRegion[r];
    a.sort((x, y) => x.lo - y.lo);
    let m = 0;
    a.forEach(iv => { m = Math.max(m, iv.hi); iv.maxEnd = m; });
  });
  return byRegion;
}
function hits(list, pos) {
  if (!list) return [];
  let lo = 0, hi = list.length - 1, last = -1;      // last interval with lo <= pos
  while (lo <= hi) { const m = (lo + hi) >> 1; if (list[m].lo <= pos) { last = m; lo = m + 1; } else hi = m - 1; }
  const out = [];
  for (let i = last; i >= 0; i--) {
    if (list[i].maxEnd < pos) break;                // nothing further back can reach pos
    if (list[i].lo <= pos && pos <= list[i].hi) out.push(list[i]);
  }
  return out;
}

// ── canonical exon intervals in GENE coordinates, and the distance from a gene position to the
// nearest one. Gene space is a linear reversal of genomic space on the minus strand, so distances
// measured here equal genomic distances — no strand handling needed for the distance itself.
function canonicalExons(gene) {
  const gs = gene.gene_structure;
  if (!gs || !gs.canonical_transcript || !gs.transcripts || !gs.exons) return null;
  const t = gs.transcripts.find(x => x.id === gs.canonical_transcript);
  if (!t || !t.exons || !t.exons.length) return null;
  const byId = {};
  gs.exons.forEach(e => { byId[e.id] = e; });
  const ex = t.exons.map(id => byId[id]).filter(Boolean).map(e => ({ start: e.start, end: e.end }));
  return ex.length ? ex.sort((a, b) => a.start - b.start) : null;
}
function distToExon(exons, gpos) {
  let best = Infinity;
  for (const e of exons) {
    if (gpos >= e.start && gpos <= e.end) return 0;
    best = Math.min(best, gpos < e.start ? e.start - gpos : gpos - e.end);
  }
  return best;
}

// canonical transcript record (needs .cds, in TRANSCRIPT coordinates -- transcriptToProtein in
// gramene-gene-positions compares a transcript position directly against cds.start/cds.end)
function canonicalTranscript(gene) {
  const gs = gene.gene_structure;
  if (!gs || !gs.canonical_transcript || !gs.transcripts) return null;
  return gs.transcripts.find(x => x.id === gs.canonical_transcript) || null;
}
// Everything in the cds-mode window that depends only on the GENE, resolved once per gene instead
// of once per variant. Previously classifyCds() rebuilt the exon array and called remap() twice for
// the CDS bounds on every single variant -- with ~50M variants per genome that dominated the run
// (295 s vs 74 s for gene mode). remap() also deep-picks the gene and rebuilds its internal exon
// keyBy on each call, so the waste compounded.
function cdsContext(gene, exons) {
  if (!exons) return { kind: 'nostruct' };
  const gs = gene.gene_structure;
  const tid = gs.canonical_transcript;
  const t = canonicalTranscript(gene);
  const plus = gene.location.strand !== -1;
  if (!t || !t.cds) {
    const a = positions.remap(gene, exons[0].start, 'gene', 'genome');
    const b = positions.remap(gene, exons[exons.length - 1].end, 'gene', 'genome');
    if (a === -1 || b === -1) return { kind: 'nostruct' };
    const txLo = Math.min(a, b), txHi = Math.max(a, b);
    return { kind: 'noncoding', exons, tid, plus, txLo, txHi,
             winLo: plus ? txLo - NC_UP : txLo,
             winHi: plus ? txHi : txHi + NC_UP };
  }
  const gStart = positions.remap(gene, t.cds.start, 'transcript', 'genome', tid); // first coding base
  const gEnd   = positions.remap(gene, t.cds.end,   'transcript', 'genome', tid); // last coding base
  if (gStart === -1 || gEnd === -1) return { kind: 'nostruct' };
  const cdsLo = Math.min(gStart, gEnd), cdsHi = Math.max(gStart, gEnd);
  const seq = CDS_SEQS ? CDS_SEQS.get(tid) : null;
  if (CDS_SEQS && !seq) noSeq++;
  if (seq && seq.length !== (t.cds.end - t.cds.start + 1)) lenBad++;
  return { kind: 'coding', exons, tid, cds: t.cds, plus, cdsLo, cdsHi, seq,
           winLo: plus ? cdsLo - CDS_UP   : cdsLo - CDS_DOWN,
           winHi: plus ? cdsHi + CDS_DOWN : cdsHi + CDS_UP };
}

// ── genetic code ─────────────────────────────────────────────────────────────────────
const CODON = (() => {
  const b = 'TCAG', aa = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG';
  const t = {}; let i = 0;
  for (const a of b) for (const c of b) for (const d of b) t[a + c + d] = aa[i++];
  return t;
})();
const translate = c => CODON[c] || 'X';
const COMP = { A: 'T', C: 'G', G: 'C', T: 'A', N: 'N', a: 't', c: 'g', g: 'c', t: 'a', n: 'n' };
const revcomp = sq => { let o = ''; for (let i = sq.length - 1; i >= 0; i--) o += (COMP[sq[i]] || 'N'); return o; };

// CDS sequences for one genome, keyed by transcript id. Loaded once: the largest of these files is
// 20 MB gzipped (~65 MB of sequence), which is far cheaper than spawning samtools per gene.
function loadCdsFasta(systemName) {
  const glob = require('child_process')
    .execSync('ls ' + FASTA_DIR + '/' + systemName + '/cds/*.fa.gz 2>/dev/null || true')
    .toString().trim().split('\n').filter(Boolean);
  if (!glob.length) return null;
  const raw = require('child_process').execSync('zcat ' + glob[0], { maxBuffer: 1 << 30 }).toString();
  const seqs = new Map();
  let id = null, buf = [];
  for (const line of raw.split('\n')) {
    if (line.startsWith('>')) {
      if (id) seqs.set(id, buf.join(''));
      id = line.slice(1).split(/\s+/)[0]; buf = [];
    } else if (line) buf.push(line.trim());
  }
  if (id) seqs.set(id, buf.join(''));
  return seqs;
}

// consequence tallies
const CONSEQ = {};
const bump = k => { CONSEQ[k] = (CONSEQ[k] || 0) + 1; };
let refOk = 0, refBad = 0, noSeq = 0, lenBad = 0;
// codon-position tally per consequence. A uniform 1/3 split is NOT expected: synonymous changes
// concentrate at position 3 (wobble) and missense at positions 1-2. A rotated peak, or a flat
// synonymous profile, would mean the CDS offset is off by one.
const FRAME = {};
const frameBump = (cls, p) => { (FRAME[cls] = FRAME[cls] || [0,0,0])[p]++; };
let CDS_SEQS = null;                       // Map(transcript id -> CDS sequence), one genome's worth

const PTV_SET = new Set(['splice', 'frameshift', 'stop_gained', 'stop_lost', 'start_lost']);
const PAV_SET = new Set(['missense', 'inframe_indel']);   // PTVs are added to PAV separately

// Effect of the ALT allele on the canonical transcript. Returns a consequence label or null when
// the variant is not in the coding region (flanks, UTR exons, deep introns) or cannot be judged.
//
// The ALT allele is given on the genomic PLUS strand; the CDS sequence is in transcript orientation.
// On a minus-strand gene the allele must therefore be complemented before it is substituted into a
// codon. Getting this wrong still yields a plausible-looking amino acid change, so it is checked by
// verifying the REF allele against the CDS base on every SNV.
function consequenceOf(ctx, gene, pos, entry) {
  const { exons, tid, cds, cdsLo, cdsHi, plus, seq } = ctx;
  if (!seq) { noSeq++; return null; }
  const ref = entry.ref, alt = entry.alt;
  if (!ref || !alt || alt === '.' || alt === '*') return null;
  const isSnv = (ref.length === 1 && alt.length === 1);
  // VCF indels are left-anchored: REF=TA/ALT=T deletes the base at pos+1, so the first base the
  // variant actually changes is one to the right of POS.
  const effPos = isSnv ? pos : pos + 1;
  if (effPos < cdsLo || effPos > cdsHi) return null;

  const gp = positions.remap(gene, effPos, 'genome', 'gene');
  if (gp === -1) return null;
  const tp = positions.remap(gene, gp, 'gene', 'transcript', tid);

  if (tp === -1) {                                  // intronic inside the coding span
    return (distToExon(exons, gp) <= SPLICE_PTV) ? 'splice' : null;
  }
  if (tp < cds.start || tp > cds.end) return null;  // exonic but UTR
  const off = tp - cds.start;                       // 0-based offset into the CDS

  if (!isSnv) {
    return ((alt.length - ref.length) % 3 !== 0) ? 'frameshift' : 'inframe_indel';
  }

  const base = seq[off];
  if (base === undefined) { noSeq++; return null; }
  const wantRef = plus ? ref.toUpperCase() : (COMP[ref.toUpperCase()] || 'N');
  if (base.toUpperCase() !== wantRef) { refBad++; return null; }
  refOk++;

  const ci = Math.floor(off / 3), cp = off % 3;
  const refCodon = seq.substr(ci * 3, 3).toUpperCase();
  if (refCodon.length < 3) return null;
  const altBase = plus ? alt.toUpperCase() : (COMP[alt.toUpperCase()] || 'N');
  const altCodon = refCodon.slice(0, cp) + altBase + refCodon.slice(cp + 1);
  const aRef = translate(refCodon), aAlt = translate(altCodon);

  let call;
  if (ci === 0 && aRef === 'M' && aAlt !== 'M') call = 'start_lost';
  else if (aRef !== '*' && aAlt === '*') call = 'stop_gained';
  else if (aRef === '*' && aAlt !== '*') call = 'stop_lost';
  else if (aRef !== aAlt) call = 'missense';
  else call = 'synonymous';
  frameBump(call, cp);
  frameBump('ALL', cp);
  return call;
}

// module scope: classifyCds() below is top-level and needs these
const classFh = CLASS_OUT ? require('fs').createWriteStream(CLASS_OUT) : null;
const emit = (rsids, gid, cls) => {
  if (classFh) rsids.forEach(r => classFh.write(r + '\t' + gid + '\t' + cls + '\n'));
};

let cCds = 0, cSplice = 0, cUp = 0, cDown = 0, cUtr = 0, cDeep = 0, cOut = 0, cNoCds = 0;
let nExon = 0, nSplice = 0, nUp = 0, nDeep = 0, nOut = 0;

// Non-coding fallback. Transcript bounds come from the canonical exon list in GENE coordinates,
// converted with remap so the upstream side follows the strand (below start on +, above end on -).
function classifyNonCoding(ctx, gene, pos, rsids, gid) {
  const { exons, tid, txLo, txHi, winLo, winHi } = ctx;

  if (pos < winLo || pos > winHi) { nOut += rsids.size; emit(rsids, gid, 'nc_out'); return false; }
  if (pos < txLo || pos > txHi) { nUp += rsids.size; emit(rsids, gid, 'nc_upstream'); return true; }
  const gpos = positions.remap(gene, pos, 'genome', 'gene');
  const tpos = (gpos === -1) ? -1 : positions.remap(gene, gpos, 'gene', 'transcript', tid);
  if (tpos !== -1) { nExon += rsids.size; emit(rsids, gid, 'nc_exon'); return true; }
  const d = distToExon(exons, gpos);
  if (d <= SPLICE) { nSplice += rsids.size; emit(rsids, gid, 'nc_splice'); return true; }
  nDeep += rsids.size; emit(rsids, gid, 'nc_deepintron'); return false;
}

// cds mode. Window is 1 kb before the start codon to 500 bp after the stop codon, in GENOMIC space
// and strand-aware: on the minus strand the start codon sits at the HIGHER genomic coordinate, so
// "upstream" extends upward. Inside the coding span we keep coding-sequence variants and variants
// within SPLICE nt of an intron end; UTR and deep-intronic are dropped.
function classifyCds(ctx, gene, pos, rsids, gid) {
  if (ctx.kind === 'nostruct') { cNoCds += rsids.size; emit(rsids, gid, 'nostruct'); return false; }
  // Non-coding gene (no CDS on the canonical transcript): there is no start/stop codon to anchor a
  // window to, so fall back to the transcript itself -- all exonic, splice-region intronic, and
  // NC_UP bp upstream of the transcript start. Without this these genes contribute nothing at all,
  // which silently loses every ncRNA/pseudogene from the index.
  if (ctx.kind === 'noncoding') return classifyNonCoding(ctx, gene, pos, rsids, gid);
  const { exons, tid, cds: t_cds, cdsLo, cdsHi, winLo, winHi, plus } = ctx;
  const t = { cds: t_cds };

  if (pos < winLo || pos > winHi) { cOut += rsids.size; emit(rsids, gid, 'outside'); return false; }
  if (pos < cdsLo || pos > cdsHi) {                      // in a flank
    const upstream = plus ? (pos < cdsLo) : (pos > cdsHi);
    if (upstream) { cUp += rsids.size; emit(rsids, gid, 'upstream'); }
    else          { cDown += rsids.size; emit(rsids, gid, 'downstream'); }
    return true;
  }
  // inside the coding span
  const gpos = positions.remap(gene, pos, 'genome', 'gene');
  const tpos = (gpos === -1) ? -1 : positions.remap(gene, gpos, 'gene', 'transcript', tid);
  if (tpos !== -1) {
    if (tpos >= t.cds.start && tpos <= t.cds.end) { cCds += rsids.size; emit(rsids, gid, 'cds'); return true; }
    cUtr += rsids.size; emit(rsids, gid, 'utr'); return false;   // exonic but not coding
  }
  const d = distToExon(exons, gpos);                     // intronic: distance to nearest exon edge
  if (d <= SPLICE) { cSplice += rsids.size; emit(rsids, gid, 'splice'); return true; }
  cDeep += rsids.size; emit(rsids, gid, 'deepintron'); return false;
}

(async () => {
  const genesCol = await collections.genes.mongoCollection();

  // pass 1 needs coordinates only — gene_structure is bulky and only genes that actually get a hit
  // need it, so it is fetched in pass 2 for that subset.
  // NB .project(), not find(query,{projection:...}): gramene-mongodb-config pins mongodb 2.2.29,
  // where the second find() argument IS the field selector and a {projection:...} object silently
  // returns _id only.
  const genes = await genesCol.find({ system_name: GENOME })
    .project({ _id: 1, 'location.region': 1, 'location.start': 1, 'location.end': 1,
                'location.strand': 1 }).toArray();
  if (!genes.length) die('no genes in mongo for system_name "' + GENOME + '" — is that a real genome?');
  if (MODE === 'cds' && !NO_CONSEQ) {
    CDS_SEQS = loadCdsFasta(GENOME);
    if (!CDS_SEQS) die('no CDS fasta under ' + FASTA_DIR + '/' + GENOME + '/cds/ — needed for ' +
                       'consequence calling; pass --no-consequence to skip PTV/PAV');
    log('loaded ' + CDS_SEQS.size + ' CDS sequences from ' + FASTA_DIR + '/' + GENOME + '/cds');
  }
  const index = buildIndex(genes);
  const regions = new Set(Object.keys(index));
  log(genes.length + ' genes across ' + regions.size + ' regions; mode=' + MODE +
      (MODE === 'cds' ? '  cds_up=' + CDS_UP + ' cds_down=' + CDS_DOWN + ' splice=' + SPLICE + ' nc_up=' + NC_UP
                      : '  flank up=' + FLANK_UP + ' down=' + FLANK_DOWN + ' intron_max=' + INTRON_MAX));

  // ── stream variants ────────────────────────────────────────────────────────────────
  const bc = spawn('bcftools', ['query', '-f', '%CHROM\t%POS\t%ID\t%REF\t%ALT\n', VCF], { stdio: ['ignore', 'pipe', 'pipe'] });
  let bcErr = '';
  bc.stderr.on('data', d => { bcErr += d; });

  const candidates = new Map();          // gene id -> Map(pos -> {ids:Set, ref, alt})
  const vcfChroms = new Set();
  const matchedChroms = new Set();
  let nVar = 0, nWithId = 0, nHitPairs = 0;

  const rl = readline.createInterface({ input: bc.stdout, terminal: false });
  for await (const line of rl) {
    if (!line) continue;
    const f = line.split('\t');
    if (f.length < 3) continue;
    nVar++;
    const chrom = normChrom(f[0]);
    vcfChroms.add(chrom);
    const ids = f[2].split(';').map(s => s.trim()).filter(s => s && s !== '.');
    if (!ids.length) continue;
    nWithId++;
    for (const id of ids) {
      if (id.includes(',')) die('rsID "' + id + '" contains a comma; add_attributes.pl splits multi-valued cells on comma, so this would corrupt the field');
    }
    const pos = +f[1];
    const hs = hits(index[chrom], pos);
    if (hs.length) matchedChroms.add(chrom);
    for (const h of hs) {
      nHitPairs++;
      let m = candidates.get(h.id);
      if (!m) { m = new Map(); candidates.set(h.id, m); }
      let e = m.get(pos);
      if (!e) { e = { ids: new Set(), ref: f[3] || '', alt: (f[4] || '').split(',')[0] }; m.set(pos, e); }
      ids.forEach(i => e.ids.add(i));
    }
  }
  await new Promise((res, rej) => bc.on('close', c => c === 0 ? res() : rej(new Error('bcftools exit ' + c + ': ' + bcErr.slice(0, 300)))))
    .catch(e => die(e.message));

  // ── chromosome-naming guard: the classic silent failure is a table that builds fine with no rows
  const overlap = [...vcfChroms].filter(c => regions.has(c));
  log('VCF: ' + nVar + ' variants, ' + nWithId + ' with an ID; CHROMs ' + vcfChroms.size +
      ', matching gene regions ' + overlap.length);
  if (!overlap.length) {
    die('no VCF CHROM matches any gene region for ' + GENOME + '.\n' +
        '        VCF CHROMs : ' + [...vcfChroms].slice(0, 10).join(', ') + '\n' +
        '        gene regions: ' + [...regions].slice(0, 10).join(', ') + '\n' +
        '        Supply --chrom-map to translate them; refusing to emit an empty table.');
  }
  // Measure coverage in GENES, not regions. A genome here has 110 regions (10 chromosomes plus
  // ~100 scaffolds) while a callset normally covers only the chromosomes, so a region-count ratio
  // reads 9% on a perfectly good VCF. Genes are what we are actually annotating.
  const matchedGenes = genes.filter(g => overlap.includes(String(g.location.region))).length;
  const frac = matchedGenes / genes.length;
  log('gene coverage of matched regions: ' + matchedGenes + '/' + genes.length +
      ' (' + (frac * 100).toFixed(1) + '%)');
  if (frac < MIN_CHROM_FRAC) {
    die('VCF CHROMs cover only ' + (frac * 100).toFixed(1) + '% of this genome\'s genes (' +
        matchedGenes + '/' + genes.length + '), below --min-chrom-frac ' + MIN_CHROM_FRAC + '.\n' +
        '        matched regions: ' + overlap.slice(0, 10).join(', ') + '\n' +
        '        Looks like a naming mismatch or the wrong VCF for this genome.');
  }

  // ── pass 2: classify, needing gene_structure only for genes that got a hit ──────────
  const ids = [...candidates.keys()];
  const struct = new Map();
  const CH = 2000;
  for (let i = 0; i < ids.length; i += CH) {
    const docs = await genesCol.find({ _id: { $in: ids.slice(i, i + CH) } })
      .project({ _id: 1, location: 1, gene_structure: 1 }).toArray();
    docs.forEach(d => struct.set(d._id, d));
  }

  let kExon = 0, kIntron = 0, dIntron = 0, kFlank = 0, noStruct = 0, capped = 0;
  const out = [];
  for (const [gid, posMap] of candidates) {
    const gene = struct.get(gid);
    if (!gene) { noStruct++; continue; }
    const exons = canonicalExons(gene);
    const ctx = (MODE === 'cds') ? cdsContext(gene, exons) : null;   // ONCE per gene, not per variant
    const keep = new Set(), ptv = new Set(), pav = new Set();
    for (const [pos, entry] of posMap) {
      const rsids = entry.ids;
      let take;
      if (MODE === 'cds') {
        take = classifyCds(ctx, gene, pos, rsids, gid);
        if (take) rsids.forEach(r => keep.add(r));
        if (ctx && ctx.kind === 'coding' && !NO_CONSEQ) {
          const c = consequenceOf(ctx, gene, pos, entry);
          if (c) {
            bump(c);
            if (PTV_SET.has(c)) rsids.forEach(r => { ptv.add(r); pav.add(r); });
            else if (PAV_SET.has(c)) rsids.forEach(r => pav.add(r));
          }
        }
        continue;
      }
      const inGene = pos >= gene.location.start && pos <= gene.location.end;
      if (!inGene) { take = true; kFlank += rsids.size; emit(rsids, gid, 'flank'); }
      else if (!exons) { take = true; noStruct++; emit(rsids, gid, 'nostruct'); }
      else {
        const gpos = positions.remap(gene, pos, 'genome', 'gene');   // strand-aware
        if (gpos === -1) { take = true; kFlank += rsids.size; emit(rsids, gid, 'flank'); }
        else {
          const d = distToExon(exons, gpos);
          if (d === 0) { take = true; kExon += rsids.size; emit(rsids, gid, 'exon'); }
          else if (d <= INTRON_MAX) { take = true; kIntron += rsids.size; emit(rsids, gid, 'nearintron'); }
          else { take = false; dIntron += rsids.size; emit(rsids, gid, 'deepintron'); }
        }
      }
      if (take) rsids.forEach(r => keep.add(r));
    }
    if (!keep.size) continue;
    let list = [...keep].sort();
    if (list.length > MAX_PER_GENE) { capped++; list = list.slice(0, MAX_PER_GENE); }
    // PTV/PAV are strict subsets and far smaller, so they are not capped. Emitting them empty is
    // fine: add_attributes.pl skips empty cells, so the field is simply absent on that gene.
    const sub = st => [...st].sort().join(',');
    out.push(gid + '\trsid\t' + list.join(',') + '\t' + sub(ptv) + '\t' + sub(pav));
  }
  out.sort();
  out.forEach(l => process.stdout.write(l + '\n'));
  if (classFh) await new Promise(r => classFh.end(r));

  if (MODE === 'cds') {
    log('coding   kept: cds ' + cCds + ', splice(<=' + SPLICE + 'nt) ' + cSplice +
        ', upstream(' + CDS_UP + ') ' + cUp + ', downstream(' + CDS_DOWN + ') ' + cDown +
        '  |  dropped: utr ' + cUtr + ', deep intronic ' + cDeep + ', outside ' + cOut);
    log('noncoding kept: exonic ' + nExon + ', splice ' + nSplice + ', upstream(' + NC_UP + ') ' + nUp +
        '  |  dropped: deep intronic ' + nDeep + ', outside ' + nOut +
        (cNoCds ? '  |  no usable structure: ' + cNoCds : ''));
  } else {
    log('kept: exonic ' + kExon + ', intronic<=' + INTRON_MAX + 'nt ' + kIntron + ', flank ' + kFlank +
        '  |  dropped: deep intronic ' + dIntron);
  }
  if (MODE === 'cds' && !NO_CONSEQ) {
    const order = ['synonymous','missense','inframe_indel','frameshift','stop_gained','stop_lost','start_lost','splice'];
    const tot = order.reduce((a, k) => a + (CONSEQ[k] || 0), 0);
    log('consequences: ' + order.map(k => k + ' ' + (CONSEQ[k] || 0)).join(', '));
    const ptvN = ['frameshift','stop_gained','stop_lost','start_lost','splice'].reduce((a,k)=>a+(CONSEQ[k]||0),0);
    const mis = CONSEQ.missense || 0, syn = CONSEQ.synonymous || 0;
    log('  PTV ' + ptvN + ' (' + (tot ? (100*ptvN/tot).toFixed(2) : 0) + '% of coding), ' +
        'missense:synonymous = ' + (syn ? (mis/syn).toFixed(2) : 'n/a'));
    // The REF allele must match the CDS base. A systematic mismatch means the FASTA and the gene
    // models are out of sync for this genome, and EVERY consequence call would be quietly wrong.
    const checked = refOk + refBad;
    if (checked) {
      const rate = refBad / checked;
      log('  REF vs CDS sequence: ' + refOk + '/' + checked + ' agree (' + (100*rate).toFixed(4) + '% mismatch)');
      if (rate > MAX_REF_MISMATCH) {
        die('REF allele disagrees with the CDS sequence for ' + (100*rate).toFixed(2) + '% of coding SNVs ' +
            '(' + refBad + '/' + checked + '), above --max-ref-mismatch ' + MAX_REF_MISMATCH + '. The CDS ' +
            'fasta and the gene models are out of sync for ' + GENOME + ' — consequence calls cannot be trusted.');
      }
    } else if (CDS_SEQS) {
      die('no coding SNV was checkable against the CDS sequence — the offset logic is not reaching the sequence');
    }
    log('  codon position of coding SNVs (pos1/pos2/pos3):');
    ['ALL','synonymous','missense','stop_gained'].forEach(k => {
      const f = FRAME[k]; if (!f) return;
      const n = f[0]+f[1]+f[2]; if (!n) return;
      log('    ' + k.padEnd(12) + f.map(v => (100*v/n).toFixed(1) + '%').join(' / ') +
          '   (n=' + n + ')');
    });
    if (noSeq)  log('  WARNING: ' + noSeq + ' lookup(s) had no CDS sequence for the canonical transcript');
    if (lenBad) log('  WARNING: ' + lenBad + ' transcript(s) whose fasta length != cds.end-cds.start+1');
  }
  log('genes with rsIDs: ' + out.length + ' (of ' + candidates.size + ' with an overlapping variant)');
  if (noStruct) log('WARNING: ' + noStruct + ' gene(s) had no usable canonical transcript — variants kept unfiltered');
  if (capped)   log('WARNING: ' + capped + ' gene(s) exceeded --max-per-gene ' + MAX_PER_GENE + ' and were truncated');
  await collections.closeMongoDatabase();
})().catch(e => { console.error('FATAL [' + GENOME + ']: ' + (e && e.stack || e)); process.exit(1); });
