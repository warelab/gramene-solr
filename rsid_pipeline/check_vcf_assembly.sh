#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# check_vcf_assembly.sh — do the rsID VCFs address the assemblies our gene models use?
#
# WHY THIS EXISTS. The rsIDs are projected from the anchor onto each pan-genome assembly. If a VCF
# was projected onto a *different build* of that accession than the one in our cores, every position
# still looks valid -- it is in range, it lands inside genes, it produces a plausible number of rows
# -- but it addresses the wrong sequence. Nothing downstream can notice: gene->rsID assignments are
# wrong and consequence calls are noise. This is the pipeline's usual failure shape (see
# build/docs/04-troubleshooting.md): well-formed, right size, wrong data.
#
# THE TEST. Compare the VCF REF allele against the genomic base at that coordinate. A VCF on the
# right assembly agrees ~100% (measured: 1500/1500 on every concordant genome). A VCF on the wrong
# one agrees ~25% -- chance, i.e. one base in four. There is no middle ground, so the 90% threshold
# below is not a tuning knob.
#
# Usage:  ./check_vcf_assembly.sh [genome ...]     (default: every VCF in RSID_VCF_DIR)
# Exit 1 if any genome is discordant, so it can gate a run.
# ─────────────────────────────────────────────────────────────────────────────
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
. "${HERE}/../../build/config.sh" 2>/dev/null || true
# config.sh turns on `set -e`; this script deliberately keeps going past a discordant genome so it
# can report ALL of them in one pass, and signals the verdict through its own exit code instead.
set +e
VCF_DIR="${RSID_VCF_DIR:-/scratch/olson/rsid_projection/vcf}"
FASTA_DIR="${RSID_CDS_FASTA_DIR:-/scratch/olson/fasta}"
N="${CHECK_N:-1500}"          # SNVs sampled per genome
MIN_PCT="${CHECK_MIN_PCT:-90}"

command -v samtools >/dev/null || { echo "samtools not found" >&2; exit 1; }

if [ "$#" -gt 0 ]; then GENOMES=("$@"); else
  mapfile -t GENOMES < <(ls "${VCF_DIR}"/*.vcf.gz 2>/dev/null | xargs -n1 basename | sed 's/\.vcf\.gz$//' | sort)
fi
[ "${#GENOMES[@]}" -gt 0 ] || { echo "no VCFs in ${VCF_DIR}" >&2; exit 1; }

one() {
  local g="$1" fa sf snv chrom hi seq clen
  fa=$(ls "${FASTA_DIR}/${g}"/dna/*dna.toplevel.fa.gz 2>/dev/null | head -1)
  [ -n "${fa}" ] || { printf '%-32s SKIP  no genomic fasta under %s/%s/dna\n' "${g}" "${FASTA_DIR}" "${g}"; return 0; }
  [ -f "${fa}.fai" ] || samtools faidx "${fa}" 2>/dev/null || { printf '%-32s ERROR faidx failed\n' "${g}"; return 2; }
  # first N biallelic SNVs on whichever sequence the VCF starts with; only the head of the file is read
  chrom=$(zcat "${VCF_DIR}/${g}.vcf.gz" 2>/dev/null | awk -F'\t' '!/^#/{print $1; exit}')
  [ -n "${chrom}" ] || { printf '%-32s ERROR no records\n' "${g}"; return 2; }
  snv=$(zcat "${VCF_DIR}/${g}.vcf.gz" 2>/dev/null | awk -F'\t' -v c="${chrom}" -v n="${N}" \
        '!/^#/ && $1==c && length($4)==1 && length($5)==1 {print $2"\t"toupper($4); k++} k>=n{exit}')
  [ -n "${snv}" ] || { printf '%-32s ERROR no SNVs on %s\n' "${g}" "${chrom}"; return 2; }
  hi=$(awk '{if($1>m)m=$1}END{print m}' <<< "${snv}")
  # the sequence goes through a file: an 80 Mb argv blows E2BIG
  sf=$(mktemp) || return 2
  samtools faidx "${fa}" "${chrom}:1-${hi}" 2>/dev/null | tail -n +2 | tr -d '\n' | tr 'a-z' 'A-Z' > "${sf}"
  [ -s "${sf}" ] || { printf '%-32s ERROR %s not in fasta\n' "${g}" "${chrom}"; rm -f "${sf}"; return 2; }
  clen=$(awk -v c="${chrom}" '$1==c{print $2}' "${fa}.fai")
  awk -v sf="${sf}" -v g="${g}" -v c="${chrom}" -v clen="${clen}" -v hi="${hi}" -v min="${MIN_PCT}" '
    BEGIN{ getline s < sf }
    { if(substr(s,$1,1)==$2) ok++; tot++ }
    END{ pct = tot ? 100*ok/tot : 0
         printf "%-32s REF_match=%6.2f%% (%d/%d)  %s_len=%s  sampled_to=%s%s\n",
                g, pct, ok, tot, c, clen, hi, (pct<min ? "   <-- WRONG ASSEMBLY" : "")
         exit (pct<min ? 3 : 0) }' <<< "${snv}"
  local rc=$?; rm -f "${sf}"; return ${rc}
}
export -f one; export VCF_DIR FASTA_DIR N MIN_PCT

RES=$(mktemp)
printf '%s\n' "${GENOMES[@]}" | xargs -P "${CHECK_JOBS:-6}" -I{} bash -c 'one "$@"' _ {} | sort > "${RES}"
cat "${RES}"
awk '/WRONG ASSEMBLY/{bad++} /REF_match/{n++} /SKIP|ERROR/{o++}
  END{ printf "\n  screened=%d  concordant=%d  WRONG_ASSEMBLY=%d  other=%d\n", n+0, n-bad, bad+0, o+0
       if (bad>0) printf "  -> exclude these via RSID_SKIP_GENOMES (see README)\n" }' "${RES}"
bad=$(grep -c 'WRONG ASSEMBLY' "${RES}"); rm -f "${RES}"
[ "${bad}" -eq 0 ]
