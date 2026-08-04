# MAKER attribute pipeline

Turns MAKER pan-gene GFFs into the `maker_attrib_table_for_solr.txt` that `60_solr_genes`
merges into the gene docs via `genes/add_attributes.pl`, giving genes the `MAKER` capability
and the `MAKER__AED__attr_f` / `MAKER__QI1..9__attr_*` fields.

Run out of band when new MAKER annotations arrive — this is not a build stage.

```bash
bash build_maker_tables.sh
```

Inputs (paths at the top of the driver):

| | |
|---|---|
| `/scratch/olson/sorghum_v11_maker_gff` | `<genome>.pan-gene.coding.gff` + `assign_stable_ids_*_<genome>.out` |
| `/scratch/olson/28_sorg_maker_AED` | the release-10 table (28 genomes) and `munge_MAKER.pl`, **read-only** |

Output: `genes/maker_attrib_table_for_solr.txt` — the release-10 genomes plus the new ones.
`config.sh` points `MAKER_TABLE` at it. The v10 copy is left alone: that release serves it
through a symlink.

`munge_MAKER.pl` is reused rather than reimplemented, so the "`-1` means missing, emit an empty
cell" rule and the exact 12-column header stay in one place.

## Canonical transcript

Longest summed CDS, ties broken by longest summed exon length, then by transcript id so repeat
runs are byte-identical. Note this is *not* always the longest transcript — a gene may well have
a longer isoform with a shorter CDS, and the CDS wins.

## Two things that will bite you

**Transcript ids are not unique across genomes.** These are pan-gene sets, so the same source id
(`3381.casb001g000110.635.1`) appears in several accessions — five of the nine, in the v11 batch.
`munge_MAKER.pl` keys its lookup on the transcript id, so a *shared* `geneID.canonical.txt`
collapses them: 391,601 canonical transcripts became 114,114 distinct gene ids, and rows were
emitted attributing one genome's scores to another's genes. Each genome is therefore munged
separately, against a lookup holding only its own transcripts. Do not "optimise" this by
concatenating the score files first.

**Gene ids join to the stable-id map only partially.** The GFF and the stable-id assignment do
not always pick the same source annotation as a gene's representative id. For `sorghum_pi655993`
only 14,387 of 46,077 genes join by id — yet every id that *does* join has byte-identical
coordinates, and each map's column 3 correctly names its own accession, so the files are
correctly paired. The extractor therefore joins by id first (authoritative) and falls back to
`(chr, start, end)`, which recovers ~46,000. Coordinate keys are near-unique; the 15-23 ambiguous
ones per file are dropped rather than guessed. The per-genome log line reports the split:

```
sorghum_pi655993  mRNA=110772 genes=46077 scores=110772 canonical=46056 (by_id=14387 by_coord=31669) unmapped=21
```

## Expected coverage

The GFFs are coding-only, so a genome yields attributes for ~89-90% of its indexed genes. The
v11 batch added 435,854 rows over 9 genomes, taking `capabilities:MAKER` from 1,015,322 to
1,451,176 and coverage from 28 to 37 genomes.

## Verifying

* every row has 12 columns; canonical count == row count == distinct gene ids
* table ids are a subset of the solr ids for those taxa (`taxon_id:<t>` per genome)
* the pre-existing rows are unchanged — `head -<n> new | diff - old`
* `-1` survives only inside gene ids such as `S369-1.001G000100`, never in a value column
* after `make refresh-attributes`, `capabilities:MAKER` rises by the new row count and the
  already-covered genomes hold steady
