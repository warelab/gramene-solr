#!/usr/bin/env python3
"""
all_species_attributes.py — clade-aware expression attributes for every Gramene genome
with expression data. Organs resolve from Plant Ontology via po_map.py to per-clade panels
(grass / eudicot; moss & algae are stress-only). Level bands are per-species quantiles;
specificity (tau over sampled organs) and the stress vocabulary are shared.

    python3 all_species_attributes.py <geneID> [<geneID> ...] [--tsv] [--out FILE]

Reuses the validated helpers from grass_expression_attributes; config = all_species_expression_panel.json.
Assay metadata is read from the discovery caches (<taxon>.assays_cache.json).
"""
import argparse, json, os, sys, statistics
from collections import defaultdict
import po_map as M
import expr_lib as G

HERE = os.path.dirname(os.path.abspath(__file__))
CFG_PATH = os.path.join(HERE, "expression_panel.json")


def interpret_contrast(cf, tf, l2fc, cfg):
    """Like the grass version, but a side carrying a baseline keyword is never treated as a
    stress match (handles e.g. 'control ozone treatment' vs 'ozone treatment')."""
    base = [b.lower() for b in cfg["differential"]["baseline_control_labels"]]

    def is_base(label):
        return label is None or any(b in label for b in base)

    out = []
    for rule in cfg["differential"]["stress_classifier"]:
        ft, cat = rule["factor_type"], rule["category"]
        match = [m.lower() for m in rule["match"]]
        cl = (cf.get(ft) or "").lower() or None
        tl = (tf.get(ft) or "").lower() or None
        c_match = cl is not None and not is_base(cl) and any(m in cl for m in match)
        t_match = tl is not None and not is_base(tl) and any(m in tl for m in match)
        if t_match and not c_match:
            if not is_base(cl):
                continue
            stressed, oriented = tl, l2fc
        elif c_match and not t_match:
            if not is_base(tl):
                continue
            stressed, oriented = cl, -l2fc
        elif c_match and t_match:
            cv = G.parse_intensity_value(cl) or 0
            tv = G.parse_intensity_value(tl) or 0
            if tv == cv:
                continue
            stressed, oriented = (tl, l2fc) if tv > cv else (cl, -l2fc)
        else:
            continue
        out.append((cat, stressed, oriented))
    return out


def setup_species(cfg, taxon):
    """(organs, sample_to_organ, assays, intensity_map, breakpoints) from the discovery cache."""
    cache = os.path.join(HERE, f"{taxon}.assays_cache.json")
    assays = json.load(open(cache)) if os.path.exists(cache) else {}
    s2o = {}
    for aid, a in assays.items():
        op = G.factor_of(a, "organism part")
        organ = M.organ_for(op.get("int_id") if op else None, op["label"] if op else None)
        if organ:
            s2o[aid] = organ
    imap = G.build_intensity_map(assays, cfg)
    clade = cfg["species"].get(str(taxon), {}).get("clade") or M.CLADE_OF_TAXON.get(taxon, "other")
    organs = M.CLADE_PANELS.get(clade, [])
    ref = cfg["species"].get(str(taxon), {}).get("reference")
    bp = ref["global_breakpoints_tpm"] if ref else None
    return organs, s2o, assays, imap, bp, clade


def score_gene(gene_id, doc, taxon, organs, s2o, assays, cfg, imap, bp, clade):
    floor = cfg["expression_level"]["detection_floor_tpm"]
    spec = cfg["specificity"]

    organ_samples = defaultdict(list)
    for field, val in doc.items():
        if val is None:
            continue
        pb = G.parse_baseline_field(field)
        if pb:
            organ = s2o.get(f"{pb[0]}.{pb[1]}")
            if organ and organ in organs:
                organ_samples[organ].append(float(val))

    organ_median = {o: statistics.median(v) for o, v in organ_samples.items()}
    sampled = [o for o in organs if o in organ_median]
    vec = [organ_median[o] for o in sampled]
    n_with_data = len(sampled)

    per_organ = {}
    for o in organs:
        if o in organ_median:
            lvl = G.level_call(organ_median[o], bp, floor) if bp else "no_ref"
            per_organ[o] = {"tpm": round(organ_median[o], 2), "level": lvl,
                            "n_samples": len(organ_samples[o])}
        else:
            per_organ[o] = {"tpm": None, "level": "no_data", "n_samples": 0}

    tau = G.compute_tau(vec)
    detected = [o for o in sampled if organ_median[o] >= floor]
    frac_detected = len(detected) / n_with_data if n_with_data else 0
    mx = max(vec) if vec else 0

    if n_with_data < spec["min_organs_required"] or mx < floor:
        breadth = ("insufficient_data" if (mx < floor and n_with_data >= 1)
                   else ("no_baseline" if n_with_data == 0 else "limited_panel"))
        peak = []
    elif tau is not None and tau >= spec["tau_tissue_specific"]:
        breadth = "tissue_specific"; peak = [o for o in organs if organ_median.get(o, 0) >= 0.5 * mx]
    elif tau is not None and tau < spec["tau_broadly_expressed"] and frac_detected >= spec["broad_min_frac_detected"]:
        breadth = "broadly_expressed"; peak = []
    else:
        breadth = "tissue_enhanced"
        dv = [organ_median[o] for o in detected] or [0]; pmed = statistics.median(dv)
        peak = [o for o in organs if organ_median.get(o, 0) >= 0.8 * mx and organ_median.get(o, 0) >= 2 * pmed]

    single_peak = [o for o in peak if len(organ_samples.get(o, [])) <= 1]
    ubiquitous = n_with_data > 0 and frac_detected >= spec["ubiquitous_min_frac_detected"]

    # differential / stress
    cut = cfg["differential"]["de_cutoffs"]
    contrasts = defaultdict(dict)
    for field, val in doc.items():
        if val is None:
            continue
        pdf = G.parse_diff_field(field)
        if pdf:
            exp, gc, gt, metric = pdf
            contrasts[(exp, gc, gt)][metric] = float(val)

    records = []; genotype_skipped = 0
    for (exp, gc, gt), m in contrasts.items():
        l2fc, p = m.get("l2fc"), m.get("p")
        if l2fc is None or p is None:
            continue
        ca, ta = assays.get(f"{exp}.{gc}"), assays.get(f"{exp}.{gt}")
        cf, tf = G.all_factors(ca), G.all_factors(ta)
        interp = interpret_contrast(cf, tf, l2fc, cfg)
        if not interp:
            genotype_skipped += 1
            continue
        for cat, stressed, oriented in interp:
            if abs(oriented) < cut["abs_log2fc_min"] or p > cut["p_value_max"]:
                continue
            direction = "activated" if oriented > 0 else "repressed"
            sa = ca if (cf.get(G._rule_factor_type(cfg, cat)) or "").lower() == stressed else ta
            tissue = (G.factor_of(sa, "organism part") or {}).get("label")
            records.append({"cat": cat, "tier_label": G.stress_label(exp, cat, stressed, imap),
                            "direction": direction,
                            "detail": {"tissue": tissue, "log2fc": round(oriented, 2), "p_value": p,
                                       "direction": direction, "experiment": exp,
                                       "contrast": f"{gc}->{gt}", "stressed_condition": stressed}})
    cat_dirs = defaultdict(set)
    for r in records:
        cat_dirs[r["cat"]].add(r["direction"])
    stress_detail = defaultdict(list); activated, repressed = set(), set()
    for r in records:
        label = r["tier_label"] if len(cat_dirs[r["cat"]]) > 1 else r["cat"]
        (activated if r["direction"] == "activated" else repressed).add(label)
        stress_detail[label].append(r["detail"])

    highly = [o for o in organs if per_organ[o]["level"] in ("high", "very_high")]
    lowly = [o for o in organs if per_organ[o]["level"] == "low"]
    attributes = {
        "specific_to": peak if breadth == "tissue_specific" else [],
        "enhanced_in": peak if breadth == "tissue_enhanced" else [],
        "broadly_expressed": breadth == "broadly_expressed",
        "ubiquitously_detected": ubiquitous,
        "not_expressed_in_panel": breadth == "insufficient_data",
        "highly_expressed_in": highly, "lowly_expressed_in": lowly,
        "activated_by": sorted(activated), "repressed_by": sorted(repressed),
    }
    return {
        "_id": gene_id, "taxon_id": taxon,
        "species": cfg["species"].get(str(taxon), {}).get("short", str(taxon)), "clade": clade,
        "name": doc.get("name"), "description": doc.get("description"),
        "attributes": attributes, "breadth": breadth,
        "tau": round(tau, 3) if tau is not None else None,
        "n_organs_with_data": n_with_data, "n_organs_detected": len(detected),
        "frac_organs_detected": round(frac_detected, 2), "ubiquitously_detected": ubiquitous,
        "max_tpm": round(mx, 2), "peak_organs": peak, "low_confidence_peak_organs": single_peak,
        "per_organ": per_organ, "stress": {k: v for k, v in stress_detail.items()},
        "n_genotype_comparisons_skipped": genotype_skipped,
    }


def tsv_row(gid, r):
    a = r["attributes"]
    if a["specific_to"]:
        loc = "specific:" + ",".join(a["specific_to"])
    elif a["broadly_expressed"]:
        loc = "broad"
    elif a["enhanced_in"]:
        loc = "enhanced:" + ",".join(a["enhanced_in"])
    elif a["not_expressed_in_panel"]:
        loc = "not_expressed"
    elif r["breadth"] == "no_baseline":
        loc = "no_baseline"
    else:
        loc = "ubiquitous" if a["ubiquitously_detected"] else "-"
    return (f"{gid}\t{r['species']}\t{r['clade']}\t{r['breadth']}\t{r['tau']}\t{r['max_tpm']}\t{loc}\t"
            f"{','.join(a['highly_expressed_in'])}\t{','.join(a['activated_by'])}\t{','.join(a['repressed_by'])}")


TSV_HEADER = "gene_id\tspecies\tclade\tbreadth\ttau\tmax_tpm\tlocation\thigh_in\tactivated_by\trepressed_by\n"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("genes", nargs="*")
    ap.add_argument("--config", default=CFG_PATH)
    ap.add_argument("--tsv", action="store_true")
    ap.add_argument("--out")
    args = ap.parse_args()
    cfg = json.load(open(args.config))
    docs = G.fetch_expression_docs(cfg, args.genes)

    by_tax = defaultdict(list)
    for gid in args.genes:
        d = docs.get(gid)
        if d is not None and d.get("taxon_id") is not None:
            by_tax[int(d["taxon_id"]) // 1000].append(gid)
    setups = {t: setup_species(cfg, t) for t in by_tax if str(t) in cfg["species"]}

    results = {}
    for gid in args.genes:
        d = docs.get(gid)
        if d is None:
            results[gid] = {"_id": gid, "error": "no expression data"}
            continue
        t = int(d["taxon_id"]) // 1000
        if t not in setups:
            results[gid] = {"_id": gid, "taxon_id": t, "error": "species not configured"}
            continue
        organs, s2o, assays, imap, bp, clade = setups[t]
        results[gid] = score_gene(gid, d, t, organs, s2o, assays, cfg, imap, bp, clade)

    if args.out:
        json.dump(results, open(args.out, "w"), indent=2); print(f"wrote {args.out}", file=sys.stderr)
    if args.tsv or not args.out:
        sys.stdout.write(TSV_HEADER)
        for gid in args.genes:
            r = results[gid]
            print(f"{gid}\t\t\tERROR\t\t\t{r['error']}" if "error" in r else tsv_row(gid, r))
    if not args.out:
        print("\n--- detail ---"); print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
