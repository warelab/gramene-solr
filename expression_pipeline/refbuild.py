#!/usr/bin/env python3
"""Per-species reference distributions on the harmonized organ panel, read DIRECTLY from
MongoDB (streams the genome's expression-bearing genes via mongo_source.iter_genome).
   python3 refbuild.py <taxon_id>  ->  ref_<taxon>.json
"""
import json, sys, os, statistics
import po_map as M
import mongo_source


def fof(a, t):
    for src in ("factor", "characteristic"):
        for f in a.get(src, []) or []:
            if f.get("type") == t:
                return f
    return None


def parse_baseline_field(field):
    if not field.endswith("__expr"):
        return None
    parts = field[:-len("__expr")].split("_")
    g = parts[-1]
    if not (g.startswith("g") and g[1:].isdigit()):
        return None
    return "-".join(parts[:-1]), g


def pct(data, p):
    if not data:
        return None
    s = sorted(data); k = (len(s) - 1) * p / 100
    f = int(k); c = min(f + 1, len(s) - 1)
    return s[f] + (s[c] - s[f]) * (k - f)


def main():
    taxon = int(sys.argv[1])
    FLOOR = 1.0
    system_name = mongo_source.system_name_for_taxon(taxon)

    # sample_to_organ for this species, from the discovery assay cache.
    cache = f"{taxon}.assays_cache.json"
    assays = json.load(open(cache)) if os.path.exists(cache) else {}
    s2o = {}
    for aid, a in assays.items():
        op = fof(a, "organism part")
        organ = M.organ_for(op.get("int_id") if op else None, op["label"] if op else None)
        if organ:
            s2o[aid] = organ
    print(f"taxon {taxon} ({system_name}): {len(s2o)} baseline sample->organ from cache", flush=True)

    if not s2o or not system_name:   # no baseline organ data -> empty reference (stress-only species)
        json.dump({"taxon_id": taxon, "n_sampled": 0, "floor_tpm": FLOOR,
                   "global": {"p25": None, "p50": None, "p75": None, "p90": None, "n": 0},
                   "per_organ": {}}, open(f"ref_{taxon}.json", "w"), indent=2)
        print(f"  no baseline organs; wrote empty ref_{taxon}.json", flush=True)
        return

    # stream every expression-bearing gene of the genome; accumulate per-organ medians.
    organ_vals = {o: [] for o in M.ALL_ORGANS}
    n = 0
    for doc in mongo_source.iter_genome(system_name):
        n += 1
        per = {}
        for field, val in doc.items():
            if not field.endswith("__expr") or val is None:
                continue
            pb = parse_baseline_field(field)
            if not pb:
                continue
            organ = s2o.get(f"{pb[0]}.{pb[1]}")
            if organ:
                per.setdefault(organ, []).append(float(val))
        for organ, vals in per.items():
            if organ in organ_vals:
                organ_vals[organ].append(statistics.median(vals))

    glob = [v for o in M.ALL_ORGANS for v in organ_vals[o] if v >= FLOOR]
    ref = {"taxon_id": taxon, "n_sampled": n, "floor_tpm": FLOOR,
           "global": {"p25": pct(glob, 25), "p50": pct(glob, 50),
                      "p75": pct(glob, 75), "p90": pct(glob, 90), "n": len(glob)},
           "per_organ": {}}
    for o in M.ALL_ORGANS:
        v = organ_vals[o]; ex = [x for x in v if x >= FLOOR]
        ref["per_organ"][o] = {"n_total": len(v), "n_expressed": len(ex),
            "p25": pct(ex, 25), "p75": pct(ex, 75), "p90": pct(ex, 90)}
    json.dump(ref, open(f"ref_{taxon}.json", "w"), indent=2)
    g = ref["global"]
    print(f"  global P25={g['p25']} P50={g['p50']} P75={g['p75']} P90={g['p90']} (n={g['n']})", flush=True)
    print(f"  wrote ref_{taxon}.json (genes streamed={n})", flush=True)


if __name__ == "__main__":
    main()
