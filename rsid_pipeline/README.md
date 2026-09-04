# rsID attribute pipeline

Builds `genes/rsid_attrib_table_for_solr.txt`, which gives gene docs the `rsid` capability and the
multi-valued `rsid__attr_ss` field, so users can **search for genes by variant rsID**.

Run out of band when new or updated VCFs arrive — this is not a build stage.

```bash
bash build_rsid_table.sh                      # every VCF in RSID_VCF_DIR
bash build_rsid_table.sh sorghum_bicolor      # just these genomes
```

Then merge it into the index, either as part of a full build (`60_solr_genes` picks the table up
automatically) or in place on the live core with no downtime:

```bash
make refresh-rsid                             # = 62_attr_atomic.sh rsid
```

## Inputs

| | default | |
|---|---|---|
| `RSID_VCF_DIR` | `/scratch/olson/rsid_projection/vcf` | one VCF per genome, named **`<system_name>.vcf.gz`** |
| `RSID_FLANK_UP` | 200 | bp **upstream** of the gene that still counts (strand-aware) |
| `RSID_FLANK_DOWN` | 0 | bp downstream; 0 means downstream variants are ignored |
| `RSID_INTRON_MAX_DIST` | 10 | intronic variants further than this from a canonical exon are dropped |
| `RSID_MAX_PER_GENE` | 5000 | per-gene ceiling |
| `RSID_MODE` | `cds` | `cds` anchors the window on the coding sequence; `gene` on the gene model |
| `RSID_CDS_UP` / `RSID_CDS_DOWN` | 1000 / 500 | bp before the start codon / after the stop (strand-aware) |
| `RSID_SPLICE` / `RSID_SPLICE_PTV` | 10 / 2 | nt each side of an intron: kept in rsid / counting as a PTV |
| `RSID_NC_UP` | 200 | non-coding genes: bp before the transcript start |
| `RSID_CDS_FASTA_DIR` | `/scratch/olson/fasta` | CDS sequence for consequence calling |
| `RSID_CHROM_MAP` | — | optional TSV `<vcf CHROM>\t<gene region>` when naming differs |
| `RSID_JOBS` | 4 | genomes extracted in parallel |

The driver sources `build/config.sh` when it can reach it, so those defaults live in one place.

A filename stem that is not a `maps.system_name` is **skipped with a warning**, not fatal: the
projection set legitimately contains accessions this release does not carry (`sorghum_pi533792`).
Per-genome runs are resumable — a genome with a `.ok` marker in the temp dir is not redone.

rsIDs are projected from an anchor genome onto each pan-genome assembly, so the same rsID appears
in several genomes' VCFs at their own coordinates. Searching one rsID therefore finds the
corresponding gene in every genome that carries it — which is the point.

## Which variants count

1. Inside the gene model, **exonic** in the canonical transcript → kept.
2. Inside the gene model, **intronic within `RSID_INTRON_MAX_DIST`** (10 nt) of a canonical exon → kept.
3. Inside the gene model, **intronic further than that** from any canonical exon → **dropped**.
   A deep intronic variant says little about the gene and only adds noise to the search field.
4. **Upstream** of the gene model, within `RSID_FLANK_UP` (200 bp) → kept. The window is
   strand-aware: below `location.start` on `+`, above `location.end` on `-`. Downstream variants
   are excluded by default.

### Why these numbers

Measured on `sorghum_bicolor` (49.9 M variants). A 2 kb symmetric flank with a 50 nt intron buffer
produced **10.6 M** rsID values per genome — about 1.09 billion across the pan-genome, ~175x the
largest attribute table this pipeline has ever merged. 200 bp upstream-only with a 10 nt buffer
gives **3.0 M** per genome (−72%), because the flank component alone falls 95%.

The per-gene cap is 5000 because only ~105 genes exceed 1000 and the largest genuinely holds 3,772.
A truncated gene is **unfindable by its dropped rsIDs**, so the cap is there to stop pathological
bloat, not to trim ordinary long genes.

"Exonic" is judged against the **canonical transcript only**. A variant intronic in the canonical
but exonic in a minor isoform is dropped unless it falls within 50 nt.

## Two things that will bite you

**Coordinates must go through `gramene-gene-positions`.** `gene_structure.exons` are gene-relative
offsets in **transcription order, not genomic order**, so on the minus strand genomic coordinates
*decrease* as gene coordinates increase. `remap()` handles it; hand-rolled arithmetic silently gets
every minus-strand gene wrong. Verified on `SORBI_3008G061701` (strand −1): a variant 50 nt into a
2,063 nt intron is kept and one at 51 nt is dropped, in both directions from the flanking exons.

**The upstream flank flips with strand.** A symmetric window is not the same thing: on a `-` strand
gene "upstream" is *above* `location.end`. Verified with synthetic variants 100 bp either side of a
`+` and a `-` gene — upstream kept, downstream dropped, on both.

**Chromosome naming.** Gene `location.region` is a bare name (`1`…`10`). A VCF using `Chr1` matches
nothing, and the failure mode is a table that builds perfectly with zero rows. The extractor
compares the VCF's CHROMs against the genome's regions and **aborts** unless they cover at least
`--min-chrom-frac` (default 0.5) of that genome's **genes** — genes, not regions, because a genome
here has ~10 chromosomes plus ~100 scaffolds and a normal callset covers only the chromosomes
(that ratio reads 9% on a perfectly good VCF, which is exactly the false positive it produced on
the first run). Use `--chrom-map` to translate; the pipeline never guesses.

## Genomes whose VCF is on the wrong assembly

Two of the 104 projected VCFs — **`sorghum_353`** and **`sorghum_pi154844`** — do not address the
assemblies this release's gene models are built on. They are excluded via `RSID_SKIP_GENOMES` in
`build/config.sh`.

This is the pipeline's classic failure shape: nothing about those files looks wrong. The
chromosome names match (`1`–`10`), the coordinate ranges are plausible, the variant counts are
normal, and every position lands inside some gene. Only the *sequence* disagrees.

**How it was caught.** The consequence caller compares the VCF REF allele against the CDS base it
computes an offset for. `sorghum_353` disagreed on **74.36%** of coding SNVs — very close to the
75% you get from a random base. That guard exists for exactly this and is the only reason the
problem surfaced at all; without consequence calling these two genomes would have shipped confident,
wrong gene→rsID assignments.

**Confirmed independently**, against the genome rather than the CDS, which rules out any error in
the CDS offset arithmetic:

* VCF REF vs the genomic base: **24.67%** and **23.13%** agreement, where all 102 other genomes
  score exactly **100.00%**. The split is bimodal — there is no borderline case.
* `sorghum_353`'s variants extend **7.5 Mb past the end of chr10** (VCF max 64,725,190 against an
  assembly length of 57,189,457), and 3–13% past the end on other chromosomes. No coordinate offset
  or liftover slip can put a variant beyond the end of the sequence; this is a different assembly.

**The screen is cheap — run it before a build, not after:**

```bash
./check_vcf_assembly.sh                      # all VCFs, ~2 min, exit 1 if any are discordant
./check_vcf_assembly.sh sorghum_353          # one genome
```

It samples 1,500 SNVs per genome and compares REF to the genomic base. The 90% threshold is not a
tuning knob: real genomes score 100%, wrong ones score ~25%.

**To fix rather than exclude**, the projection for these accessions has to be redone against the
assemblies we serve. Then clear `RSID_SKIP_GENOMES` and re-run — `check_vcf_assembly.sh` should
report 104 concordant.

Separately, `sorghum_pi533792` has a VCF but is **not a genome in this release**; the driver already
skips unknown stems with a warning, and it needs no entry in the skip list.

## Performance

Genotypes are never materialised — `bcftools query -f '%CHROM\t%POS\t%ID\n'` only. The same rice
callset on this host is **194 GB** with genotypes and **67 KB** without.

Overlap is a per-region interval index with binary search, not `bedtools` (not installed here) and
not a `tabix` call per gene (infeasible across 5.4M genes). `gene_structure` is fetched only for
genes that actually have an overlapping variant, since it is far bulkier than coordinates.

## Consequence fields

Two subset fields are derived from the alternate allele's effect on the **canonical transcript**:

| field | contents |
|---|---|
| `rsid__attr_ss` | every rsID in the window (below) |
| `rsid_PAV__attr_ss` | protein-altering: all PTVs, plus missense and in-frame indels |
| `rsid_PTV__attr_ss` | protein-truncating: stop gained, stop lost, start lost, frameshift, splice site |

`PTV ⊆ PAV ⊆ rsid`, asserted for every gene. All three share the single `rsid` capability token.

| condition | consequence | fields |
|---|---|---|
| intronic, ≤`RSID_SPLICE_PTV` (2) nt from an exon edge, intron inside the CDS | splice | PTV, PAV |
| indel in CDS, length change not divisible by 3 | frameshift | PTV, PAV |
| indel in CDS, length change divisible by 3 | in-frame indel | PAV |
| SNV, alt codon is a stop and ref is not | stop gained | PTV, PAV |
| SNV, ref codon is a stop and alt is not | stop lost | PTV, PAV |
| SNV in codon 1 disrupting `ATG` | start lost | PTV, PAV |
| SNV, different amino acid | missense | PAV |
| SNV, same amino acid | synonymous | neither |

Non-coding genes (no CDS on the canonical transcript) get `rsid__attr_ss` only — there is no protein
to truncate. Sequence comes from `RSID_CDS_FASTA_DIR/<system_name>/cds/*.fa.gz`, loaded once per
genome (largest is 20 MB gzipped; do **not** spawn `samtools faidx` per gene).

### The trap: the ALT allele is on the plus strand

The CDS sequence is in transcript orientation, but VCF REF/ALT are always on the genomic plus
strand. On a minus-strand gene the allele **must be complemented** before it is substituted into a
codon. Get this wrong and every call still looks like a plausible amino-acid change — nothing
crashes, nothing looks odd, the biology is just wrong.

Two guards make that detectable rather than silent:

* **REF must match the CDS base.** Every coding SNV is checked; the run **aborts** above
  `--max-ref-mismatch` (0.1%). On `sorghum_bicolor` this reads **1,037,137/1,037,137 agree,
  0.0000% mismatch** — a whole-genome check, not a sample.
* **Codon-position tally per consequence**, reported each run. The genetic code makes the expected
  shape unmistakable:

```
  ALL          27.2% / 23.7% / 49.0%    (n=1,037,137)
  synonymous    4.3% /  0.1% / 95.6%
  missense     45.9% / 43.8% / 10.2%
  stop_gained  58.0% / 22.5% / 19.5%
```

Synonymous at **0.1% second-position** is the sharpest signal: in the standard code a second-position
change is almost always non-synonymous. An off-by-one in the CDS offset would rotate this pattern and
destroy the complementarity between synonymous and missense. A flat-ish overall split is *not*
expected and would itself indicate a bug.

## Why `string` and not a numeric field

rsIDs strip trivially to integers (`rs6453096324` -> `6453096324`), so a `plong` field looks like an
obvious win. It was measured, on two identical cores holding 200k docs x 95 real rsIDs (19M values):

| | `string` | `plong` |
|---|---|---|
| stored fields (`.fdt`) | 121.9 MB | 84.2 MB (-31%) |
| indexed structure | 74.4 MB (`.tim`/`.tip`) | 77.3 MB (`.kdd`) — **a wash** |
| docValues (`.dvd`) | 107.8 MB | 95.8 MB |
| **index dir total** | **291 MB** | **247 MB (-15%)** |
| index time | 105 s | 29 s |
| exact-match QTime | 0 ms | 0 ms |

The 12-bytes-to-8 saving shows up only in the *stored* section. On the *indexed* side there is none:
the string term dictionary is FST prefix-compressed and every id shares the `rs` prefix plus long
numeric runs, so it ends up slightly smaller than the BKD points tree. Query latency is identical —
this field only ever serves exact-match lookups, which both formats answer sub-millisecond.

15% of one field's index does not justify the cost: Solr here has **no `plong` dynamic field**, so
adding one forces a full reindex (~2.7 h), `add_attributes.pl` would need its
`^\w+_attr_[sif]s?$` guard widened, and ids would be permanently constrained to `rs\d+`.

Note also that `*_attr_is` is **`pint`** (32-bit) and is not an option regardless: 80% of these
rsIDs exceed int32 (max observed 5,453,096,324).

## Debugging

`--emit-classes <file>` writes one `rsid<TAB>gene<TAB>class` row per (gene, variant) pair, where
class is `exon` / `nearintron` / `flank` / `deepintron` / `nostruct`. That is how the flank-vs-exon
trade-offs above were measured — e.g. establishing that dropping a flank mapping when the rsID is
exonic in another gene only removes 0.2% of kept pairs at a 200 bp window, so it is not worth the
extra logic. Expect ~13M rows for one genome; it is a diagnostic, not part of a normal run.

## Verifying

Three scripts, in the order you would run them:

```bash
./check_vcf_assembly.sh                       # BEFORE a run: are the VCFs on our assemblies? (~2 min)
./run_summary.sh  /scratch/olson/rsid_projection/work_conseq   # per-genome QC: REF%, PTV%, mis:syn
./validate_rsid_table.sh                      # structural checks on the finished table
```

Each exits non-zero on a problem, so they can gate a build. What they look for:

| script | catches |
| --- | --- |
| `check_vcf_assembly.sh` | a VCF projected onto a different assembly (see above) |
| `run_summary.sh` | a genome with an anomalous PTV rate or missense:synonymous ratio, an unfinished genome, and REF disagreement. Genomes named in `RSID_SKIP_GENOMES` are reported as expected exclusions rather than failures |
| `validate_rsid_table.sh` | header drift, `PTV ⊄ PAV ⊄ rsid`, whitespace in a value, and — given a reference table — any change to `rsid__attr_ss` |

**Reconcile the row count.** The table's row count must equal the sum of the per-genome `.tsv`
files for the genomes actually included. This is the one check that catches a work-dir glob sweeping
in excluded or stale genomes — every structural check still passes on a polluted table:

```bash
for f in "$RSID_WORK_DIR"/*.tsv; do ... ; done   # sum only the included genomes
```

Reference values measured across the pan-genome: **REF agreement 100.00%** on every concordant
genome, **PTV 6.4–10.6%** of coding variants, **missense:synonymous 1.04–1.52**. A genome on the
wrong assembly lands far outside all three at once (REF ~25%, PTV ~35%, mis:syn ~2.3), which is
what makes these cheap checks decisive rather than suggestive.

Non-coding genes carry `rsid__attr_ss` but neither consequence field — 1,089/1,089 ncRNA genes on
the anchor, by construction.

Then, the older manual checks:

* per-genome log line reports kept exonic / kept intronic / kept flank / **dropped deep intronic**;
  a dropped count of 0 on a dense VCF means the filter is not doing anything
* cross-check a few genes against `tabix <vcf> <region>:<start-2000>-<end+2000>` with
  `--intron-max 100000` (filter off) — the sets should match exactly
* `--intron-max 100000` must yield a strict superset; `--flank 0` a strict subset
* after merging, `capabilities:rsid` in the genes core should equal the table's row count
