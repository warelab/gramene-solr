#!/usr/bin/env bash
# build_rsid_suggestions.sh — one suggestion doc per distinct rsID, straight from the
# per-genome extractor output in ${RSID_WORK_DIR}.
#
# WHY NOT FACETED, unlike every other suggestion category here: aux.js builds TF/expression/
# ontology categories with facet.pivot=<field>,taxon_id&facet.limit=-1 against the genes core.
# rsid__attr_ss carries ~10M distinct terms over ~460M values; that pivot would materialize the
# whole term space inside the 8 GB heap this Solr shares across 90 cores (it OOMed 2026-09-12).
# The files below hold the same information and cost nothing but disk.
#
# SHAPE: 104 per-genome sorts that each fit in RAM, then ONE streaming k-way merge. A single
# 452M-line sort would need ~13 GB of temp and a multi-pass external merge; `sort -m` over
# pre-sorted inputs needs no temp files at all — provided --batch-size exceeds the file count
# (sort's default merge batch is 16, which would silently reintroduce multi-pass spilling).
#
# usage: build_rsid_suggestions.sh [output.jsonl]
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
. "${HERE}/../../build/config.sh" >/dev/null

OUT="${1:-${RSID_SUGG_JSON:-/scratch/olson/rsid_projection/rsid_suggestions.jsonl}}"
TMP="${RSID_SUGG_TMP:-/scratch/olson/rsid_projection/sugg_tmp}"   # /scratch has TB; NEVER /tmp (20 GB root)
AGG="${TMP}/agg"
JOBS="${RSID_SUGG_JOBS:-8}"
mkdir -p "${AGG}" "$(dirname "${OUT}")"

# ── system_name -> taxon_id ───────────────────────────────────────────────────
# The genes core is the authority: it indexes PER-GENOME taxon ids (sorghum_bicolor is 4558006,
# not 4558), and a suggestion's taxon_id array must match what the genes core is filtered on.
# One pivot, ~4 s, one taxon per genome. Cached so reruns and the refresh path don't re-query.
MAP="${TMP}/system_name_taxon.tsv"
if [ ! -s "${MAP}" ]; then
  echo "building system_name -> taxon_id map from ${SOLR_GENES_CORE}" >&2
  curl -s "${SOLR_URL}/${SOLR_GENES_CORE}/query?rows=0&facet=true&facet.limit=-1&facet.mincount=1&json.nl=map&facet.pivot=system_name,taxon_id" \
  | "${NODE_BIN}" -e 'let s="";process.stdin.on("data",d=>s+=d).on("end",()=>{
      let p;
      try { p = JSON.parse(s).facet_counts.facet_pivot["system_name,taxon_id"] || []; }
      catch (e) { console.error("could not parse pivot response: " + e.message); process.exit(1); }
      let n = 0;
      for (const d of p) {
        const t = d.pivot || [];
        if (t.length !== 1) { console.error("ambiguous taxon for " + d.value + " (" + t.length + ")"); continue; }
        console.log(d.value + "\t" + t[0].value); n++;
      }
      if (!n) { console.error("empty system_name pivot"); process.exit(1); }
    })' > "${MAP}.tmp"
  mv "${MAP}.tmp" "${MAP}"
fi
[ -s "${MAP}" ] || { echo "could not build ${MAP}" >&2; exit 1; }
echo "taxon map: $(wc -l < "${MAP}") genomes" >&2

# ── phase A: per genome -> sorted "rsid \t taxon \t genes \t maxflag" ─────────
# maxflag 2=PTV, 1=PAV, 0=plain. PTV ⊆ PAV ⊆ rsid, so the max over the genes a variant hits in
# this genome is its strongest consequence here.
phaseA() {
  local tsv="$1" g t
  g="$(basename "${tsv}" .tsv)"
  case " ${RSID_SKIP_GENOMES} " in
    *" ${g} "*) echo "skip ${g} (RSID_SKIP_GENOMES)" >&2; return 0 ;;
  esac
  # An extractor run that ABORTS still leaves its partial .tsv behind; the .ok marker is written
  # only on success, so it — not the presence of the file — is what says the data is complete.
  # sorghum_pi154844 is exactly this: 26,351 rows from a run that failed the REF guard at 74.52%
  # mismatch, the only one of 104 .tsv files with no marker. It is also in RSID_SKIP_GENOMES, so
  # this is belt and braces: the next genome to fail won't be, and its bad gene->rsID assignments
  # would otherwise join the output silently.
  if [ ! -e "${RSID_WORK_DIR}/${g}.ok" ]; then
    echo "skip ${g} — no .ok marker (incomplete/aborted extraction)" >&2
    return 0
  fi
  t="$(awk -F'\t' -v g="${g}" '$1==g{print $2; exit}' "${MAP}")"
  [ -n "${t}" ] || { echo "no taxon_id for ${g} — not a maps.system_name?" >&2; return 0; }
  awk -F'\t' -v OFS='\t' '
    { delete ptv; delete pav
      if ($4 != "") { n = split($4, b, ","); for (i=1;i<=n;i++) if (b[i] != "") ptv[b[i]] = 1 }
      if ($5 != "") { n = split($5, c, ","); for (i=1;i<=n;i++) if (c[i] != "") pav[c[i]] = 1 }
      if ($3 == "") next
      n = split($3, a, ",")
      for (i=1;i<=n;i++) { r = a[i]; if (r == "") continue
        print r, (r in ptv ? 2 : (r in pav ? 1 : 0)) } }' "${tsv}" \
  | LC_ALL=C sort -t$'\t' -k1,1 -S 512M -T "${TMP}" \
  | awk -F'\t' -v OFS='\t' -v T="${t}" '
      { if ($1 != p) { if (p != "") print p, T, c, m; p = $1; c = 0; m = 0 }
        c++; if ($2+0 > m) m = $2+0 }
      END { if (p != "") print p, T, c, m }' > "${AGG}/${g}.agg.tmp"
  mv "${AGG}/${g}.agg.tmp" "${AGG}/${g}.agg"
  echo "  ${g}: $(wc -l < "${AGG}/${g}.agg") distinct rsIDs (taxon ${t})" >&2
}

echo "phase A: aggregating ${JOBS} genomes at a time" >&2
i=0
for tsv in "${RSID_WORK_DIR}"/*.tsv; do
  [ -e "${tsv}" ] || { echo "no .tsv files in ${RSID_WORK_DIR}" >&2; exit 1; }
  phaseA "${tsv}" &
  i=$((i+1))
  [ $((i % JOBS)) -eq 0 ] && wait
done
wait
# Floor guards against a silently partial run (a half-populated work dir would produce a
# plausible-looking file with millions of rsIDs missing). Lower it only to test the pipeline.
n_agg=$(ls "${AGG}"/*.agg 2>/dev/null | wc -l)
[ "${n_agg}" -ge "${RSID_MIN_GENOMES:-100}" ] \
  || { echo "only ${n_agg} per-genome aggregates (need ${RSID_MIN_GENOMES:-100}) — aborting" >&2; exit 1; }
echo "phase A done: ${n_agg} genomes" >&2

# ── phase B: k-way merge -> one JSONL doc per distinct rsID ───────────────────
# LC_ALL=C and -t$'\t' -k1,1 MUST match phase A exactly. A locale or key mismatch makes the merge
# output unordered, which splits one rsID across several docs — same id, last write wins, counts
# silently wrong.
echo "phase B: merging ${n_agg} sorted files" >&2
LC_ALL=C sort -m --batch-size=128 -t$'\t' -k1,1 -S 256M -T "${TMP}" "${AGG}"/*.agg \
| awk -F'\t' '
  BEGIN { BAD = 0 }
  # NB: scat, not sub — sub() is an awk built-in and cannot be used as a variable name.
  function flush(   i, ids, frq, scat, rel) {
    if (p == "") return
    if (p !~ /^rs[0-9]+$/) { BAD++; return }
    ids = ""; frq = ""
    for (i=1;i<=nt;i++) { ids = ids (i>1 ? "," : "") tid[i]; frq = frq (i>1 ? "," : "") tcnt[i] }
    if      (maxf == 2) { scat = "PTV";  rel = "1.15" }
    else if (maxf == 1) { scat = "PAV";  rel = "1.10" }
    else                { scat = "rsID"; rel = "1.05" }
    printf "{\"category\":\"Variants: rsID\",\"subcategory\":\"%s\",\"fq_field\":\"rsid__attr_ss\",\"fq_value\":\"%s\",\"id\":\"rsid_%s\",\"display_name\":\"%s\",\"name\":\"%s\",\"num_genes\":%d,\"relevance\":%s,\"taxon_id\":[%s],\"taxon_freq\":[%s]}\n", \
      scat, p, p, p, p, tot, rel, ids, frq
  }
  { if ($1 != p) { flush(); p = $1; nt = 0; tot = 0; maxf = 0 }
    nt++; tid[nt] = $2; tcnt[nt] = $3; tot += $3; if ($4+0 > maxf) maxf = $4+0 }
  END { flush(); if (BAD) { printf "REJECTED %d malformed rsIDs\n", BAD > "/dev/stderr"; exit 1 } }' \
> "${OUT}.tmp"
mv "${OUT}.tmp" "${OUT}"
echo "wrote $(wc -l < "${OUT}") rsID suggestions to ${OUT}" >&2
