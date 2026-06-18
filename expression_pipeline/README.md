# Expression Attribute Pipeline

Builds per-gene **expression attributes** for any Gramene/SorghumBase-style genome and
emits them as MongoDB-ready JSONL. For each gene it summarises where it is expressed
(per-clade organ panel, with `not/low/medium/high/very_high` level calls), how broadly
(tau specificity → `tissue_specific` / `tissue_enhanced` / `broadly_expressed`), and which
stresses it responds to (`activated_by` / `repressed_by`). Designed for display on Compara
gene trees, so attributes are harmonised across species via Plant Ontology.

**Pure Python standard library — no third-party dependencies.**

## Files

| File | Role |
|---|---|
| `manifest.json` | **The only per-site config.** API base + the list of genomes (`system_name`, `taxon`, `clade`, `short`, `name`). |
| `po_map.py` | Curation surface: universal PO→organ map, label/substring fallbacks, per-clade organ panels (`CLADE_PANELS`), optional `CLADE_OF_TAXON`. |
| `expr_lib.py` | Shared helper library: API fetch, expr-field parsing, tau, level binning, stress-contrast logic. |
| `discover.py` | Taxon-agnostic discovery: finds a genome's experiments from gene expr-fields, caches assay metadata (`<taxon>.assays_cache.json`) + an organ/stress inventory (`inv_<taxon>.json`). |
| `refbuild.py` | Per-species quantile reference distributions (`ref_<taxon>.json`). |
| `assemble_config.py` | Builds `expression_panel.json` from the manifest + references. Holds the pipeline constants (level/specificity thresholds, DE cutoffs, **stress vocabulary**). |
| `score.py` | The scorer: one gene doc → one attribute record. |
| `batch.py` | Batch runner: paginates a genome, streams `_id`-keyed JSONL, resumable. |
| `build.py` | End-to-end driver: discover → references → config → score, for every genome in the manifest. |

Dependency graph: `po_map` and `expr_lib` are leaves; `refbuild`, `assemble_config`, `score`
import `po_map`; `score` + `batch` import `expr_lib`; `batch` imports `score`; `build` shells
out to the rest.

## Build

```
python3 build.py                    # full build for all genomes in manifest.json
python3 build.py --only 4558 39947  # just these taxa
python3 build.py --rediscover       # ignore cached assays and re-discover
```

Outputs one `<short>_attributes.jsonl` per genome. Load into one MongoDB collection keyed
on `_id` (the gene id is globally unique across genomes):

```
mongoimport --collection expression_attributes --file sorghum_attributes.jsonl
```

Steps are idempotent: existing `*.assays_cache.json` / `ref_*.json` are reused, and
`batch.py` resumes from its JSONL line count, so an interrupted build continues where it
left off. To run a single stage manually:

```
python3 discover.py sorghum_bicolor 4558
python3 refbuild.py 4558
python3 assemble_config.py
python3 batch.py --taxon 4558
python3 score.py SORBI_3001G241000 --tsv     # ad-hoc single-gene lookup
```

(`refbuild.py` / `batch.py` run to completion by default; both accept a `--budget` seconds
limit and resume on re-invocation, for sandboxed environments with per-call time caps.)

## Output schema (one JSON object per JSONL line)

```json
{
  "_id": "SORBI_3001G241000", "taxon_id": 4558, "species": "sorghum", "clade": "grass",
  "name": "...", "description": "...",
  "attributes": {
    "specific_to": ["embryo"], "enhanced_in": [], "broadly_expressed": false,
    "ubiquitously_detected": false, "not_expressed_in_panel": false,
    "highly_expressed_in": [], "lowly_expressed_in": ["pericarp"],
    "activated_by": [], "repressed_by": ["salinity"]
  },
  "breadth": "tissue_specific", "tau": 0.968,
  "n_organs_with_data": 13, "n_organs_detected": 3, "max_tpm": 16.0,
  "peak_organs": ["embryo"], "low_confidence_peak_organs": [],
  "per_organ": { "...": {"tpm": ..., "level": "...", "n_samples": ...} },
  "stress": { "salinity": [ {"tissue": "leaf", "log2fc": -6.8, "p_value": ..., "direction": "repressed", ...} ] }
}
```

## Adding a genome / a new site

1. **New site** (different data API): edit `manifest.json` — set `api_base` and list the genomes.
2. **New genome**: add a manifest entry. If it's from a clade already in `po_map.CLADE_PANELS`
   (`grass`, `eudicot`, `bryophyte`, `algae`), nothing else is needed.
3. **New clade**: add an organ panel to `po_map.CLADE_PANELS` and use that clade name in the manifest.
4. **New organs / pathogens**: extend `po_map` (`PO_TO_ORGAN`, `SUBSTR_FALLBACK`) and the
   `STRESS_CLASSIFIER` in `assemble_config.py`. Run `discover.py` and read `inv_<taxon>.json`
   to see which organism-part terms and stress factors a genome actually uses.

Notes: genomes without baseline organ data (e.g. unicellular algae) get stress + level
attributes only (`breadth = no_baseline`); genomes whose baseline covers < 5 organs get
levels + stress but no specificity call (`breadth = limited_panel`). Both are honest, expected
states. Level bands are per-species quantiles, so "high" means the same rank within each
species; tau is computed only over the organs a species actually samples, so sparser panels
aren't biased toward false specificity.
