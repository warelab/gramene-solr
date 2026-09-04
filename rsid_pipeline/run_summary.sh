#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# run_summary.sh — per-genome QC for a finished (or in-flight) rsID extraction run.
#
# Reads the per-genome logs in the work dir and reports the three numbers that reveal a bad genome:
#
#   REF%     agreement between the VCF REF allele and the CDS base. Anything below 99.9% means the
#            VCF is not on the assembly our gene models use -- see check_vcf_assembly.sh and the
#            "wrong assembly" section of README.md. Real genomes score 100.00%.
#   PTV%     protein-truncating variants as a share of coding variants. Observed 6.4-10.6% across
#            the pan-genome. 0% or ~50% is a frame bug, not a finding.
#   mis:syn  missense:synonymous ratio, ~1.0-1.5 here. A wild departure means the codon frame is off.
#
# Parse note: read the extractor's own "PTV N (X% of coding)" line. Do NOT re-derive PTV by summing
# the per-class counts with a loose regex -- the codon-position block further down the log also
# contains "stop_gained <pct>", and a naive match overwrites the real count and silently understates
# PTV by ~12%.
#
# Usage:  ./run_summary.sh [work_dir]
# ─────────────────────────────────────────────────────────────────────────────
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
. "${HERE}/../../build/config.sh" 2>/dev/null || true
set +e
W="${1:-${RSID_WORK_DIR:-/scratch/olson/rsid_projection/work}}"
[ -d "${W}" ] || { echo "no work dir: ${W}" >&2; exit 1; }

printf '%-32s %9s %8s %10s %8s %8s\n' GENOME ROWS REF% PTV PTV% MIS:SYN
for f in "${W}"/*.log; do
  [ -e "${f}" ] || continue
  g=$(basename "${f}" .log)
  rows=$([ -s "${W}/${g}.tsv" ] && wc -l < "${W}/${g}.tsv" || echo 0)
  done_ok=$([ -f "${W}/${g}.ok" ] && echo yes || echo NO)
  skipped=$(grep -qw -- "${g}" <<< "${RSID_SKIP_GENOMES:-}" && echo 1 || echo 0)
  awk -v g="${g}" -v rows="${rows}" -v ok="${done_ok}" -v skipped="${skipped}" '
    /REF vs CDS sequence:/ { if (match($0,/\(([0-9.]+)% mismatch\)/,m)) refmm=m[1]+0 }
    /PTV [0-9]+ \([0-9.]+% of coding\)/ {
      if (match($0,/PTV ([0-9]+) \(([0-9.]+)%/,p)) { ptv=p[1]; ptvp=p[2]+0 }
      if (match($0,/missense:synonymous = ([0-9.]+)/,r)) ms=r[1]+0 }
    /REF allele disagrees/ { failed=1 }
    END{
      flag=""
      if (skipped=="1")        flag=" <-- excluded (RSID_SKIP_GENOMES), not in the table"
      else if (failed || refmm>0.1) flag=" <-- WRONG ASSEMBLY"
      else if (ok=="NO")       flag=" <-- incomplete"
      else if (ptvp<3||ptvp>20) flag=" <-- PTV out of range"
      else if (ms<0.7||ms>2.0)  flag=" <-- mis:syn out of range"
      printf "%-32s %9d %8s %10s %8s %8s%s\n", g, rows,
             (refmm==""?"-":sprintf("%.2f",100-refmm)), (ptv==""?"-":ptv),
             (ptvp==""?"-":sprintf("%.2f",ptvp)), (ms==""?"-":sprintf("%.2f",ms)), flag
    }' "${f}"
done | sort > /tmp/.rsid_sum.$$

cat /tmp/.rsid_sum.$$
awk '!/<--/ && $5!="-" {n++; s+=$5; m+=$6; if(n==1||$5<mn)mn=$5; if($5>mx)mx=$5}
     /WRONG ASSEMBLY/{bad++} /incomplete/{inc++} /out of range/{oor++} /excluded \(RSID_SKIP/{exc++}
     END{ printf "\n  genomes with results = %d\n", n+0
          if(n) printf "  PTV%%    mean=%.2f  range %.2f-%.2f   (flagged genomes excluded)\n  mis:syn mean=%.2f\n", s/n, mn, mx, m/n
          printf "  wrong assembly=%d  incomplete=%d  out of range=%d  excluded=%d\n", bad+0, inc+0, oor+0, exc+0 }' /tmp/.rsid_sum.$$
bad=$(grep -c 'WRONG ASSEMBLY\|out of range' /tmp/.rsid_sum.$$)   # excluded rows are labelled differently and do not count; rm -f /tmp/.rsid_sum.$$
[ "${bad}" -eq 0 ]
