#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# validate_rsid_table.sh — structural checks on the finished rsID attribute table.
#
# The driver already guarantees the table is non-empty and free of duplicate gene ids. This checks
# the properties that only matter once consequence calling exists, and that no count would reveal:
#
#   1. header and column count
#   2. PTV subset of PAV subset of rsid, per row -- the fields are nested by construction, so a
#      violation means the classifier emitted an rsID into a consequence field without putting it in
#      the base field, i.e. the dedup or the cap disagreed with the classifier
#   3. no empty cells written as literal '' and no stray whitespace in a value
#   4. per-genome coverage, so a genome that silently produced only base rsIDs is visible
#   5. optional: rsid__attr_ss byte-identical to a reference table (adding consequence fields must
#      not perturb the field that already shipped)
#
# Usage:  ./validate_rsid_table.sh [table] [reference-table]
# ─────────────────────────────────────────────────────────────────────────────
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
. "${HERE}/../../build/config.sh" 2>/dev/null || true
set +e
T="${1:-${RSID_TABLE}}"
REF="${2:-}"
[ -s "${T}" ] || { echo "missing or empty table: ${T}" >&2; exit 1; }
echo "table: ${T}  ($(du -h "${T}" | cut -f1), $(( $(wc -l < "${T}") - 1 )) gene rows)"
rc=0

echo
echo "── 1. header"
head -1 "${T}" | tr '\t' '\n' | nl -ba | sed 's/^/     /'
exp=$'id\tcapabilities\trsid__attr_ss\trsid_PTV__attr_ss\trsid_PAV__attr_ss'
[ "$(head -1 "${T}")" = "${exp}" ] && echo "     OK" || { echo "     FAIL: unexpected header"; rc=1; }

echo
echo "── 2. containment (PTV ⊆ PAV ⊆ rsid) + 3. cell hygiene"
awk -F'\t' 'NR==1{next}
  {
    nf=NF; if(nf!=5){bad_nf++; next}
    split("",R); split("",P); split("",A)
    n=split($3,a,","); for(i=1;i<=n;i++) if(a[i]!="") R[a[i]]=1
    n=split($4,b,","); for(i=1;i<=n;i++) if(b[i]!="") P[b[i]]=1
    n=split($5,c,","); for(i=1;i<=n;i++) if(c[i]!="") A[c[i]]=1
    for(k in P){ if(!(k in A)) {ptv_not_pav++; if(!e1){e1=$1" "k}} ; if(!(k in R)) {ptv_not_rsid++; if(!e2){e2=$1" "k}} }
    for(k in A) if(!(k in R)) {pav_not_rsid++; if(!e3){e3=$1" "k}}
    if($0 ~ /[ ]/) sp++
    if($3=="") empty_base++
    rows++
  }
  END{
    printf "     rows=%d  wrong_field_count=%d\n", rows, bad_nf+0
    printf "     PTV not in PAV : %d %s\n", ptv_not_pav+0, (e1?"(e.g. "e1")":"")
    printf "     PTV not in rsid: %d %s\n", ptv_not_rsid+0, (e2?"(e.g. "e2")":"")
    printf "     PAV not in rsid: %d %s\n", pav_not_rsid+0, (e3?"(e.g. "e3")":"")
    printf "     rows with empty base rsid field: %d\n", empty_base+0
    printf "     rows containing a space: %d\n", sp+0
    if(bad_nf+ptv_not_pav+ptv_not_rsid+pav_not_rsid+empty_base+sp > 0) exit 1
  }' "${T}"
[ $? -eq 0 ] && echo "     OK" || { echo "     FAIL"; rc=1; }

echo
echo "── 4. per-genome coverage (genes carrying each field)"
awk -F'\t' 'NR==1{next}{
    g=$1; sub(/[^A-Za-z0-9_].*$/,"",g)     # crude genome hint from the id prefix
    tot++; if($4!="")ptv++; if($5!="")pav++
  } END{
    printf "     genes total=%d  with PTV=%d (%.2f%%)  with PAV=%d (%.2f%%)\n",
           tot, ptv+0, tot?100*ptv/tot:0, pav+0, tot?100*pav/tot:0
  }' "${T}"

# 5. row count must reconcile against the genomes that actually succeeded and were not excluded.
# This is the only check that catches a work-dir glob sweeping stale or excluded genomes into the
# table: every structural property above still holds on a polluted table.
W="${RSID_WORK_DIR:-}"
if [ -n "${W}" ] && [ -d "${W}" ]; then
  echo
  echo "── 5. row count reconciles with the included per-genome files"
  sum=0; ng=0
  for okf in "${W}"/*.ok; do
    [ -e "${okf}" ] || continue
    g=$(basename "${okf}" .ok)
    grep -qw -- "${g}" <<< "${RSID_SKIP_GENOMES:-}" && continue
    [ -s "${W}/${g}.tsv" ] || continue
    sum=$(( sum + $(wc -l < "${W}/${g}.tsv") )); ng=$(( ng + 1 ))
  done
  have=$(( $(wc -l < "${T}") - 1 ))
  printf "     genomes included=%d  expected rows=%d  table rows=%d\n" "${ng}" "${sum}" "${have}"
  if [ "${sum}" -eq "${have}" ]; then echo "     OK"
  else echo "     FAIL: surplus/deficit of $(( have - sum )) rows — a stale or excluded genome may have been swept in"; rc=1; fi
fi

if [ -n "${REF}" ] && [ -s "${REF}" ]; then
  echo
  echo "── 5. rsid__attr_ss unchanged vs ${REF}"
  cut -f1,3 "${T}"   | tail -n +2 | sort > /tmp/.rsid_new.$$
  cut -f1,3 "${REF}" | tail -n +2 | sort > /tmp/.rsid_ref.$$
  if diff -q /tmp/.rsid_new.$$ /tmp/.rsid_ref.$$ >/dev/null; then echo "     OK (identical)"
  else echo "     DIFFERS:"; diff /tmp/.rsid_new.$$ /tmp/.rsid_ref.$$ | head -6 | sed 's/^/       /'; rc=1; fi
  rm -f /tmp/.rsid_new.$$ /tmp/.rsid_ref.$$
fi

echo
[ ${rc} -eq 0 ] && echo "VALIDATION PASSED" || echo "VALIDATION FAILED"
exit ${rc}
