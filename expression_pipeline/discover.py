#!/usr/bin/env python3
"""
discover.py — per-species discovery + inventory, read DIRECTLY from MongoDB.

  python3 discover.py <system_name> <taxon>

Reads all assays of the experiments this taxon participates in (mongo assays + experiments
collections), caches them to <taxon>.assays_cache.json (the format score.py / refbuild.py
expect), and writes inv_<taxon>.json with the organism-part PO histogram (baseline) and the
stress-factor histogram (differential). The scoring/curation logic is unchanged; only the
data source moved from the HTTP API to mongo.
"""
import json, sys
from collections import Counter, defaultdict
import mongo_source


def main():
    system_name, taxon = sys.argv[1], int(sys.argv[2])

    assays = mongo_source.load_assays(system_name)
    experiments = mongo_source.experiments_for_genome(system_name)
    json.dump(assays, open(f"{taxon}.assays_cache.json", "w"))

    # inventory
    organ_hist = Counter()                 # (po_int,label) -> n  (baseline assays)
    stress_hist = defaultdict(Counter)     # factor_type -> label -> n (differential assays)
    base_exps = {e for e, d in experiments.items() if d.get("type") == "Baseline"}
    diff_exps = {e for e, d in experiments.items() if d.get("type") == "Differential"}
    for aid, a in assays.items():
        exp = a.get("experiment")
        op = None
        for src in ("factor", "characteristic"):
            for f in a.get(src, []) or []:
                if f.get("type") == "organism part":
                    op = f
                    break
            if op:
                break
        if exp in base_exps:
            organ_hist[(op.get("int_id") if op else None, op.get("label") if op else None)] += 1
        if exp in diff_exps:
            for f in a.get("factor", []) or []:
                if f.get("type") in ("organism part",):
                    continue
                stress_hist[f["type"]][f.get("label")] += 1

    inv = {"system_name": system_name, "taxon": taxon,
           "n_experiments": len(experiments),
           "n_baseline": len(base_exps), "n_differential": len(diff_exps),
           "n_assays": len(assays),
           "organs": [{"po": po, "label": lab, "n": n} for (po, lab), n in organ_hist.most_common()],
           "stress_factors": {ft: dict(c.most_common()) for ft, c in stress_hist.items()}}
    json.dump(inv, open(f"inv_{taxon}.json", "w"), indent=1)
    print(f"{system_name} taxon={taxon}: exps={len(experiments)} (B={len(base_exps)} D={len(diff_exps)}) "
          f"assays={len(assays)} organs={len(organ_hist)}", flush=True)


if __name__ == "__main__":
    main()
