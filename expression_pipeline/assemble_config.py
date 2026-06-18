#!/usr/bin/env python3
"""
assemble_config.py — build expression_panel.json from manifest.json + the per-species
reference distributions (ref_<taxon>.json). Manifest-driven: api base, genome list, and
each genome's clade/short/name come from manifest.json; the organ panels come from po_map.py.

    python3 assemble_config.py            # writes expression_panel.json

Pipeline constants (level bands, specificity thresholds, DE cutoffs, the stress vocabulary)
live here — they are the same across sites; only manifest.json changes per site.
"""
import json, os

HERE = os.path.dirname(os.path.abspath(__file__))


def ref_block(taxon):
    f = os.path.join(HERE, f"ref_{taxon}.json")
    if not os.path.exists(f):
        return None
    r = json.load(open(f)); g = r["global"]
    if g.get("p25") is None:                       # no expressed baseline organ data
        return None
    return {"n_sampled": r["n_sampled"],
            "global_breakpoints_tpm": {"p25": g["p25"], "p75": g["p75"], "p90": g["p90"]},
            "per_organ_breakpoints_tpm": {o: {"p25": r["per_organ"][o]["p25"],
                "p75": r["per_organ"][o]["p75"], "p90": r["per_organ"][o]["p90"]}
                for o in r["per_organ"] if r["per_organ"][o]["n_expressed"] > 0}}


# --- stress vocabulary (curation surface; extend as new experiments appear) ---
BASELINE_CONTROL_LABELS = ["none", "mock", "water", "control", "buffer", "dimethyl sulfoxide",
                           "dmso", "normal watering", "wild type", "ambient", "untreated",
                           "0 hour", "treated directly"]

STRESS_CLASSIFIER = [
 {"category": "drought", "factor_type": "environmental stress", "match": ["drought", "water deficit", "dehydration"]},
 {"category": "salinity", "factor_type": "environmental stress", "match": ["saline", "salt"], "intensity_tiers": True},
 {"category": "cold", "factor_type": "environmental stress", "match": ["cold", "frost", "10 degree"]},
 {"category": "heat", "factor_type": "environmental stress", "match": ["warm/hot", "heat"]},
 {"category": "submergence", "factor_type": "environmental stress", "match": ["submergence", "flooding", "waterlog", "anaerobic", "anoxi"]},
 {"category": "nitrogen", "factor_type": "environmental stress", "match": ["nitrogen deprivation"]},
 {"category": "phosphate", "factor_type": "environmental stress", "match": ["phosphate"]},
 {"category": "ozone", "factor_type": "environmental stress", "match": ["ozone"]},
 {"category": "wounding", "factor_type": "environmental stress", "match": ["wound"]},
 {"category": "bacterial_pathogen", "factor_type": "environmental stress", "match": ["pseudomonas", "xanthomonas"]},
 {"category": "herbivore", "factor_type": "environmental stress", "match": ["tetranychus", "oligonychus"]},
 {"category": "salinity", "factor_type": "compound", "match": ["sodium chloride"]},
 {"category": "osmotic", "factor_type": "compound", "match": ["mannitol", "polyethylene glycol"]},
 {"category": "phosphate", "factor_type": "compound", "match": ["phosphate"]},
 {"category": "ABA", "factor_type": "compound", "match": ["abscisic acid"]},
 {"category": "jasmonate", "factor_type": "compound", "match": ["jasmon"]},
 {"category": "salicylic_acid", "factor_type": "compound", "match": ["salicylic", "acibenzolar"]},
 {"category": "auxin", "factor_type": "compound", "match": ["indole-3-acetic", "auxin"]},
 {"category": "cadmium", "factor_type": "compound", "match": ["cadmium"]},
 {"category": "elicitor", "factor_type": "compound", "match": ["flg22", "mamp", "flagellin"]},
 {"category": "ER_stress", "factor_type": "compound", "match": ["tunicamycin", "dithiothreitol"]},
 {"category": "ozone", "factor_type": "compound", "match": ["ozone"]},
 {"category": "salinity", "factor_type": "growth condition", "match": ["sodium chloride"]},
 {"category": "drought", "factor_type": "growth condition", "match": ["drought"]},
 {"category": "wounding", "factor_type": "stimulus", "match": ["injury", "wound"]},
 {"category": "submergence", "factor_type": "stimulus", "match": ["anaerobic", "anoxi"]},
 {"category": "fungal_pathogen", "factor_type": "infect", "match": ["fusarium", "blumeria", "puccinia", "botrytis",
    "neofusicoccum", "pyrenochaeta", "phaeomoniella", "exserohilum", "erysiphe", "magnaporthe", "verticillium",
    "ustilago", "pyricularia"]},
 {"category": "oomycete", "factor_type": "infect", "match": ["phytophthora", "plasmopara", "pythium"]},
 {"category": "bacterial_pathogen", "factor_type": "infect", "match": ["xanthomonas", "pseudomonas", "xylella", "ralstonia", "erwinia"]},
 {"category": "nematode", "factor_type": "infect", "match": ["globodera", "meloidogyne", "heterodera", "nematode"]},
 {"category": "virus", "factor_type": "infect", "match": ["virus"]},
 {"category": "herbivore", "factor_type": "infect", "match": ["tetranychus", "oligonychus", "aphid", "sesamia"]},
]


def main():
    import po_map as M
    import mongo_source
    manifest = json.load(open(os.path.join(HERE, "manifest.json")))
    genomes = mongo_source.discover_genomes(manifest.get("genomes"))   # from mongo (manifest overrides)
    species = {}
    for g in genomes:
        t = g["taxon"]
        species[str(t)] = {"name": g["name"], "short": g["short"],
                           "system_name": g["system_name"], "taxon_id": t,
                           "clade": g["clade"], "reference": ref_block(t)}
    cfg = {
        "schema_version": "1.0",
        "description": "Expression attribute scheme. Per-clade organ panels via Plant Ontology; "
                       "per-species quantile level bands; shared specificity + stress vocabulary.",
        "source": {"type": "mongodb", "db": mongo_source.MONGO_DB,
                   "collections": ["expression", "assays", "experiments", "genes"]},
        "clade_panels": M.CLADE_PANELS,
        "expression_level": {"method": "per-species genome-wide quantile bins on per-(gene,organ) median TPM",
            "detection_floor_tpm": 1.0,
            "bands": [{"label": "not_expressed", "rule": "<1"}, {"label": "low", "rule": "1..P25"},
                      {"label": "medium", "rule": "P25..P75"}, {"label": "high", "rule": "P75..P90"},
                      {"label": "very_high", "rule": ">=P90"}],
            "use_per_organ_breakpoints": False},
        "specificity": {"method": "tau over organs the species samples (Yanai 2005)", "min_organs_required": 5,
            "tau_tissue_specific": 0.85, "tau_broadly_expressed": 0.50, "broad_min_frac_detected": 0.70,
            "ubiquitous_min_frac_detected": 0.90},
        "differential": {"de_cutoffs": {"abs_log2fc_min": 1.0, "p_value_max": 0.05},
            "direction": "log2FC oriented so positive = rises with stress; other side must be an unstressed baseline.",
            "baseline_control_labels": BASELINE_CONTROL_LABELS,
            "stress_classifier": STRESS_CLASSIFIER},
        "species": species,
    }
    out = os.path.join(HERE, "expression_panel.json")
    json.dump(cfg, open(out, "w"), indent=1)
    nbase = sum(1 for s in species.values() if s["reference"])
    print(f"wrote expression_panel.json: {len(species)} genomes, {nbase} with references")


if __name__ == "__main__":
    main()
