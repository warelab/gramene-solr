#!/usr/bin/env bash
# build_maker_tables.sh — MAKER attribute table for the sorghum v11 build.
#
# Adds the 9 pan-genome accessions that arrived with MAKER GFFs to the 28 genomes already
# covered by the release-10 table, and writes a merged 37-genome table for v11.
#
# Run once, out of band; this is not a build stage. Reruns are idempotent.
#
#   bash build_maker_tables.sh
#
# NB /scratch/olson/28_sorg_maker_AED is READ-ONLY here: the v10 release serves its
# maker_attrib_table_for_solr.txt through a symlink, so it must keep what it has.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
GFF_DIR="/scratch/olson/sorghum_v11_maker_gff"
OLD_DIR="/scratch/olson/28_sorg_maker_AED"
V11_GENES="/usr/local/gramene/subsites/sorghum/v11/gramene-solr/genes"
OUT_TABLE="${V11_GENES}/maker_attrib_table_for_solr.txt"

MUNGE="${OLD_DIR}/munge_MAKER.pl"
OLD_TABLE="${OLD_DIR}/maker_attrib_table_for_solr.txt"

[ -d "${GFF_DIR}" ]  || { echo "missing ${GFF_DIR}" >&2; exit 1; }
[ -x "${MUNGE}" ]    || { echo "missing ${MUNGE}" >&2; exit 1; }
[ -s "${OLD_TABLE}" ]|| { echo "missing ${OLD_TABLE}" >&2; exit 1; }

cd "${HERE}"
rm -f ./*.AED_QI.scores geneID.canonical.txt maker_attrib_table_9new.txt
mkdir -p per_genome

# Each genome is munged SEPARATELY. Transcript ids are NOT unique across these pan-gene sets —
# 3381.casb001g000110.635.1 occurs in five of the nine accessions — so a single shared
# geneID.canonical.txt collapses them (391,601 canonical transcripts -> 114,114 distinct gene
# ids) and emits rows attributing one genome's scores to another's genes. Scoping the lookup to
# one genome at a time keeps each transcript id unambiguous.
echo "== extracting and munging per genome, from ${GFF_DIR}"
n_genomes=0
for gff in "${GFF_DIR}"/*.pan-gene.coding.gff; do
  genome="$(basename "${gff}" .pan-gene.coding.gff)"
  idmap="$(ls "${GFF_DIR}"/assign_stable_ids_*_"${genome}".out 2>/dev/null | head -1)"
  if [ -z "${idmap}" ]; then
    echo "  !! no stable-id file for ${genome} — skipping" >&2
    continue
  fi
  python3 "${HERE}/extract_maker_scores.py" "${gff}" "${idmap}" "${genome}" "${HERE}" \
    > "per_genome/${genome}.canonical.txt"

  # munge_MAKER.pl reads ./geneID.canonical.txt from its cwd, keeps only canonical transcripts
  # and blanks -1 QI fields. Reused as-is so the -1 rule and the 12-column header live in one
  # place; run in a scratch dir holding just this genome's lookup.
  rm -rf .munge && mkdir -p .munge
  cp "per_genome/${genome}.canonical.txt" .munge/geneID.canonical.txt
  ( cd .munge && perl "${MUNGE}" < "${HERE}/${genome}.AED_QI.scores" ) \
    > "per_genome/${genome}.attrib.txt"
  cat "per_genome/${genome}.canonical.txt" >> geneID.canonical.txt
  n_genomes=$((n_genomes+1))
done
rm -rf .munge
echo "   genomes processed: ${n_genomes}"
echo "   canonical transcripts: $(wc -l < geneID.canonical.txt)"

echo "== combining the per-genome tables"
{
  head -1 "$(ls per_genome/*.attrib.txt | head -1)"        # header once
  for f in per_genome/*.attrib.txt; do tail -n +2 "$f"; done
} > maker_attrib_table_9new.txt
echo "   new rows: $(( $(wc -l < maker_attrib_table_9new.txt) - 1 ))"
echo "   distinct gene ids: $(tail -n +2 maker_attrib_table_9new.txt | cut -f1 | sort -u | wc -l)"

echo "== merging with the existing 28-genome table"
mkdir -p "${V11_GENES}"
{
  cat "${OLD_TABLE}"                       # header + the 28 existing genomes
  tail -n +2 maker_attrib_table_9new.txt   # new rows, header dropped
} > "${OUT_TABLE}.tmp"
mv -f "${OUT_TABLE}.tmp" "${OUT_TABLE}"

echo "   ${OUT_TABLE}"
echo "   total rows: $(( $(wc -l < "${OUT_TABLE}") - 1 ))  (was $(( $(wc -l < "${OLD_TABLE}") - 1 )))"
echo "   ${OLD_DIR} untouched: $(ls -la "${OLD_TABLE}" | awk '{print $6,$7,$8}')"
