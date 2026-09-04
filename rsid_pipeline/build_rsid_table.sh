#!/usr/bin/env bash
# build_rsid_table.sh — rsID attribute table for the solr gene docs.
#
# Run out of band when new/updated VCFs arrive; this is not a build stage. Reruns are idempotent.
#
#   bash build_rsid_table.sh                       # all VCFs in RSID_VCF_DIR
#   bash build_rsid_table.sh sorghum_bicolor ...   # only these genomes
#
# Input : ${RSID_VCF_DIR}/<system_name>.vcf.gz     (one per target genome)
# Output: ${RSID_TABLE}  ->  id \t capabilities \t rsid__attr_ss
#
# The rsIDs are projected from an anchor genome onto each pan-genome assembly, so the same rsID
# appears in several genomes' VCFs at their own coordinates -- searching one rsID then finds the
# corresponding gene in every genome that carries it.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"

# Defer to build/config.sh when it is reachable: it is the single place that knows paths for this
# release, and duplicating its defaults here is how the two silently drift apart. Falls back to the
# values below so the pipeline still runs standalone from a bare checkout.
_CFG="${HERE}/../../build/config.sh"
if [ -z "${RSID_VCF_DIR:-}" ] && [ -r "${_CFG}" ]; then
  # shellcheck disable=SC1090
  . "${_CFG}" >/dev/null 2>&1 || true
fi

VCF_DIR="${RSID_VCF_DIR:-/scratch/olson/rsid_projection/vcf}"
OUT="${RSID_TABLE:-${HERE}/../genes/rsid_attrib_table_for_solr.txt}"
FLANK_UP="${RSID_FLANK_UP:-200}"
FLANK_DOWN="${RSID_FLANK_DOWN:-0}"
INTRON_MAX="${RSID_INTRON_MAX_DIST:-10}"
MAX_PER_GENE="${RSID_MAX_PER_GENE:-5000}"
JOBS="${RSID_JOBS:-4}"
MODE="${RSID_MODE:-cds}"
CDS_UP="${RSID_CDS_UP:-1000}"
CDS_DOWN="${RSID_CDS_DOWN:-500}"
SPLICE="${RSID_SPLICE:-10}"
SPLICE_PTV="${RSID_SPLICE_PTV:-2}"
NC_UP="${RSID_NC_UP:-200}"
FASTA_DIR="${RSID_CDS_FASTA_DIR:-/scratch/olson/fasta}"
CHROM_MAP="${RSID_CHROM_MAP:-}"
MONGO_DB_NAME="${MONGO_DB:-}"

[ -d "${VCF_DIR}" ] || { echo "missing VCF dir: ${VCF_DIR} (set RSID_VCF_DIR)" >&2; exit 1; }
command -v bcftools >/dev/null || { echo "bcftools not found" >&2; exit 1; }

# which genomes?
declare -a GENOMES=()
if [ "$#" -gt 0 ]; then
  GENOMES=("$@")
else
  shopt -s nullglob
  for f in "${VCF_DIR}"/*.vcf.gz "${VCF_DIR}"/*.vcf; do
    b="$(basename "$f")"; GENOMES+=("${b%%.vcf*}")
  done
  shopt -u nullglob
fi
[ "${#GENOMES[@]}" -gt 0 ] || { echo "no VCFs found in ${VCF_DIR}" >&2; exit 1; }

# every stem must be a real genome -- a typo must not silently produce a smaller table
KNOWN="$(mongosh --quiet ${MONGO_DB_NAME:+"${MONGO_DB_NAME}"} --eval \
  'db.maps.distinct("system_name").join("\n")' 2>/dev/null || true)"
# A VCF whose stem is not a genome in this release is skipped with a warning rather than aborting:
# the projection set can legitimately contain accessions this release does not carry (e.g.
# sorghum_pi533792, which has no maps/genes/core entry here). Skipping is reported loudly so it
# cannot pass unnoticed.
if [ -n "${KNOWN}" ]; then
  declare -a KEEP=()
  for g in "${GENOMES[@]}"; do
    if grep -qxF "${g}" <<< "${KNOWN}"; then KEEP+=("${g}")
    else echo "WARNING: skipping '${g}' — not a maps.system_name in this release" >&2; fi
  done
  GENOMES=("${KEEP[@]}")
  [ "${#GENOMES[@]}" -gt 0 ] || { echo "FATAL: no VCF matched a genome in this release" >&2; exit 1; }
fi

# Genomes whose VCF is on a DIFFERENT ASSEMBLY than this release's gene models. Their variant
# coordinates do not address our sequence at all, so both the consequence calls and the gene->rsID
# assignments would be wrong -- silently, since a variant still lands inside *some* gene. Detected by
# check_vcf_assembly.sh (VCF REF vs the genomic base: concordant genomes score ~100%, these ~24%,
# which is chance). Excluding them is a deliberate, visible decision: fix the projection upstream and
# empty RSID_SKIP_GENOMES rather than letting a bad genome through.
if [ -n "${RSID_SKIP_GENOMES:-}" ]; then
  declare -a KEEP2=()
  for g in "${GENOMES[@]}"; do
    if grep -qw -- "${g}" <<< "${RSID_SKIP_GENOMES//,/ }"; then
      echo "WARNING: EXCLUDING '${g}' — RSID_SKIP_GENOMES (VCF is on a different assembly; see rsid_pipeline/README.md)" >&2
    else KEEP2+=("${g}"); fi
  done
  GENOMES=("${KEEP2[@]}")
  [ "${#GENOMES[@]}" -gt 0 ] || { echo "FATAL: every genome was excluded" >&2; exit 1; }
fi

# Durable work dir: the run writes ~4 GB of per-genome TSVs and the .ok markers are what make it
# resumable, so this must not live on /tmp (20 GB root fs here -- a run already died of ENOSPC
# there, losing 62 completed genomes). RSID_TMPDIR still overrides for one-off runs.
TMP="${RSID_TMPDIR:-${RSID_WORK_DIR:-${VCF_DIR}/../work}}"
mkdir -p "${TMP}" || { echo "FATAL: cannot create work dir ${TMP}" >&2; exit 1; }
avail_kb="$(df -Pk "${TMP}" | awk 'NR==2{print $4}')"
[ "${avail_kb:-0}" -ge 8000000 ] || { echo "FATAL: only $((avail_kb/1024)) MB free at ${TMP}; a full run needs ~4 GB (set RSID_WORK_DIR to a bigger filesystem)" >&2; exit 1; }
echo "rsID table: ${#GENOMES[@]} genome(s), mode=${MODE} cds_up=${CDS_UP} cds_down=${CDS_DOWN} splice=${SPLICE}/${SPLICE_PTV} nc_up=${NC_UP} cap=${MAX_PER_GENE} jobs=${JOBS}" >&2

one_genome() {
  local g="$1" vcf=""
  for cand in "${VCF_DIR}/${g}.vcf.gz" "${VCF_DIR}/${g}.vcf"; do
    [ -s "${cand}" ] && { vcf="${cand}"; break; }
  done
  [ -n "${vcf}" ] || { echo "FATAL: no VCF for ${g} in ${VCF_DIR}" >&2; return 1; }
  # resumable: a completed genome is not redone (each run is ~75s, 100+ genomes)
  if [ -s "${TMP}/${g}.tsv" ] && [ -f "${TMP}/${g}.ok" ]; then
    echo "── ${g}  (cached, $(wc -l < "${TMP}/${g}.tsv") rows)" >&2; return 0
  fi
  echo "── ${g}  ($(basename "${vcf}"))" >&2
  node --max-old-space-size=24576 "${HERE}/extract_rsids.js" --genome "${g}" --vcf "${vcf}" \
       --mode "${MODE}" --flank-up "${FLANK_UP}" --flank-down "${FLANK_DOWN}" \
       --intron-max "${INTRON_MAX}" --max-per-gene "${MAX_PER_GENE}" \
       --cds-up "${CDS_UP}" --cds-down "${CDS_DOWN}" --splice "${SPLICE}" \
       --splice-ptv "${SPLICE_PTV}" --nc-up "${NC_UP}" --fasta-dir "${FASTA_DIR}" \
       ${CHROM_MAP:+--chrom-map "${CHROM_MAP}"} > "${TMP}/${g}.tsv" 2> "${TMP}/${g}.log" \
    || { echo "FATAL: ${g} failed — see ${TMP}/${g}.log" >&2; tail -3 "${TMP}/${g}.log" >&2
         rm -f "${TMP}/${g}.tsv"   # no .ok, and no partial file for a later run to sweep up
         return 1; }
  touch "${TMP}/${g}.ok"
  echo "   ${g}: $(wc -l < "${TMP}/${g}.tsv") rows" >&2
  return 0
}
export -f one_genome
export VCF_DIR HERE TMP FLANK_UP FLANK_DOWN INTRON_MAX MAX_PER_GENE CHROM_MAP
export MODE CDS_UP CDS_DOWN SPLICE SPLICE_PTV NC_UP FASTA_DIR

printf '%s\n' "${GENOMES[@]}" | xargs -P "${JOBS}" -I{} bash -c 'one_genome "$@"' _ {} \
  || { echo "FATAL: one or more genomes failed; nothing written" >&2; exit 1; }

# every requested genome must have produced a file — a silent gap would ship a partial table
for g in "${GENOMES[@]}"; do
  [ -f "${TMP}/${g}.ok" ] || { echo "FATAL: ${g} produced no output" >&2; exit 1; }
done

# one table, header first. Gene ids are unique across genomes, so a plain concatenation is correct.
# Concatenate the SELECTED genomes by name -- never a glob over the work dir. The work dir is
# durable and resumable, so it also holds .tsv files for genomes excluded by RSID_SKIP_GENOMES
# (partial output from the run in which they failed) and for genomes from earlier runs with other
# parameters. A `cat "${TMP}"/*.tsv` silently swept those into the table: it added 58,181 rows from
# the two wrong-assembly genomes that had just been deliberately excluded, and nothing downstream
# would have noticed.
{ printf 'id\tcapabilities\trsid__attr_ss\trsid_PTV__attr_ss\trsid_PAV__attr_ss\n'
  for g in "${GENOMES[@]}"; do cat "${TMP}/${g}.tsv"; done; } > "${OUT}.part"
rows=$(( $(wc -l < "${OUT}.part") - 1 ))
[ "${rows}" -gt 0 ] || { echo "FATAL: produced 0 rows — refusing to write an empty table" >&2; exit 1; }

# ids must be unique: a duplicate would mean two genomes claim the same gene id, and
# add_attributes.pl would keep only one of them
dups=$(tail -n +2 "${OUT}.part" | cut -f1 | sort | uniq -d | head -5)
[ -z "${dups}" ] || { echo "FATAL: duplicate gene ids across genomes: ${dups}" >&2; exit 1; }

mv "${OUT}.part" "${OUT}"
echo "wrote ${OUT} (${rows} rows)" >&2
