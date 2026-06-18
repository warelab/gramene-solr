#!/usr/bin/env python3
"""
grass_expression_attributes.py — cross-species (sorghum / maize / rice) gene expression
attributes on a harmonized organ panel, so attributes are comparable on a Compara gene tree.

Organs are derived from Plant Ontology ids (species-agnostic); expression-level bands are
per-species quantiles; specificity (tau) and stress vocabulary are shared. Config lives in
grass_expression_panel.json.

    python3 grass_expression_attributes.py SORBI_3001G241000 Zm00001eb... Os01g... --tsv
    python3 grass_expression_attributes.py --file genes.txt --out attrs.json

Only the Python standard library is required.
"""
import argparse, json, os, re, sys, time, urllib.request, urllib.parse, statistics
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
CFG_PATH = os.path.join(HERE, "grass_expression_panel.json")

SEVERITY_WORDS = {"mild": 1, "light": 1, "low": 1, "weak": 1, "moderate": 2, "medium": 2,
                  "intermediate": 2, "high": 3, "strong": 3, "severe": 3, "extreme": 3}


# ----------------------------------------------------------------------------- API
def api_get(base, path, params, retries=4, timeout=120):
    url = base + path + ("?" + urllib.parse.urlencode(params, doseq=True) if params else "")
    last = None
    for _ in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as r:
                return json.load(r)
        except Exception as e:  # noqa: BLE001
            last = e; time.sleep(2)
    raise RuntimeError(f"GET {url} failed: {last}")


def fetch_expression_docs(cfg, gene_ids):
    """Gene docs (flat expr fields + name/description/taxon_id) read DIRECTLY from MongoDB
    (the expression + genes collections), not the HTTP API. Same return shape as before."""
    import mongo_source
    return mongo_source.fetch_docs(gene_ids)


# ------------------------------------------------------------------- field parsing
def parse_baseline_field(field):
    if not field.endswith("__expr"):
        return None
    parts = field[:-len("__expr")].split("_")
    g = parts[-1]
    if not (g.startswith("g") and g[1:].isdigit()):
        return None
    return "-".join(parts[:-1]), g


def parse_diff_field(field):
    if field.endswith("_l2fc_attr_f"):
        metric, stem = "l2fc", field[:-len("_l2fc_attr_f")]
    elif field.endswith("_pval_attr_f"):
        metric, stem = "p", field[:-len("_pval_attr_f")]
    else:
        return None
    parts = stem.split("_")
    if len(parts) < 3:
        return None
    gc, gt = parts[-2], parts[-1]
    if not (gc.startswith("g") and gc[1:].isdigit() and gt.startswith("g") and gt[1:].isdigit()):
        return None
    return "-".join(parts[:-2]), gc, gt, metric


def factor_of(assay, typ):
    if not assay:
        return None
    for src in ("factor", "characteristic"):
        for f in assay.get(src, []) or []:
            if f.get("type") == typ:
                return f
    return None


def all_factors(assay):
    out = {}
    if assay:
        for f in assay.get("factor", []) or []:
            out[f["type"]] = f.get("label")
    return out


# ---------------------------------------------------------------- organ resolution
def make_organ_for(cfg):
    op = cfg["organ_panel"]
    po2 = {int(k): v for k, v in op["po_to_organ"].items()}
    fb = {k.lower(): v for k, v in op["label_fallbacks"].items()}
    ex_po = set(op["exclude_po"]); ex_lab = {l.lower() for l in op["exclude_labels"]}

    def f(po_int, label):
        if po_int in ex_po:
            return None
        if label and label.lower() in ex_lab:
            return None
        if po_int in po2:
            return po2[po_int]
        if label and label.lower() in fb:
            return fb[label.lower()]
        return None
    return f


# ------------------------------------------------------------------- level + tau
def level_call(tpm, bp, floor):
    if tpm < floor:
        return "not_expressed"
    if tpm < bp["p25"]:
        return "low"
    if tpm < bp["p75"]:
        return "medium"
    if tpm < bp["p90"]:
        return "high"
    return "very_high"


def compute_tau(values):
    if not values:
        return None
    mx = max(values)
    if mx <= 0 or len(values) < 2:
        return None
    return sum(1 - (v / mx) for v in values) / (len(values) - 1)


# -------------------------------------------------------------------- stress logic
def parse_intensity_value(label_lower):
    m = re.search(r"(\d+(?:\.\d+)?)\s*(percent|%|molar|micromolar|millimolar|mm|um|µm)", label_lower)
    if m:
        return float(m.group(1))
    for w, v in SEVERITY_WORDS.items():
        if re.search(r"\b" + re.escape(w) + r"\b", label_lower):
            return float(v)
    return None


def build_intensity_map(assays, cfg):
    base = [b.lower() for b in cfg["differential"]["baseline_control_labels"]]
    rules = [r for r in cfg["differential"]["stress_classifier"] if r.get("intensity_tiers")]
    acc = defaultdict(dict)
    for a in assays.values():
        exp = a.get("experiment")
        facs = {f["type"]: (f.get("label") or "") for f in a.get("factor", []) or []}
        for r in rules:
            lab = facs.get(r["factor_type"])
            if not lab:
                continue
            ll = lab.lower()
            if any(b in ll for b in base) or not any(m in ll for m in r["match"]):
                continue
            acc[(exp, r["category"])][ll] = parse_intensity_value(ll)
    tiers = {}
    for (exp, cat), labmap in acc.items():
        uniq = sorted(labmap, key=lambda x: (labmap[x] is None, labmap[x] if labmap[x] is not None else 0))
        if len(uniq) <= 1:
            for l in uniq:
                tiers[(exp, cat, l)] = ""
        elif len(uniq) == 2:
            tiers[(exp, cat, uniq[0])] = "low"; tiers[(exp, cat, uniq[1])] = "high"
        else:
            for i, l in enumerate(uniq):
                frac = i / (len(uniq) - 1)
                tiers[(exp, cat, l)] = "low" if frac < 0.34 else ("medium" if frac < 0.67 else "high")
    return tiers


def interpret_contrast(cf, tf, l2fc, cfg):
    base = [b.lower() for b in cfg["differential"]["baseline_control_labels"]]

    def is_base(label):
        return label is None or any(b in label for b in base)

    out = []
    for rule in cfg["differential"]["stress_classifier"]:
        ft, cat = rule["factor_type"], rule["category"]
        match = [m.lower() for m in rule["match"]]
        cl = (cf.get(ft) or "").lower() or None
        tl = (tf.get(ft) or "").lower() or None
        c_match = cl is not None and any(m in cl for m in match)
        t_match = tl is not None and any(m in tl for m in match)
        if t_match and not c_match:
            if not is_base(cl):
                continue
            stressed, oriented = tl, l2fc
        elif c_match and not t_match:
            if not is_base(tl):
                continue
            stressed, oriented = cl, -l2fc
        elif c_match and t_match:
            cv = parse_intensity_value(cl) or 0
            tv = parse_intensity_value(tl) or 0
            if tv == cv:
                continue
            stressed, oriented = (tl, l2fc) if tv > cv else (cl, -l2fc)
        else:
            continue
        out.append((cat, stressed, oriented))
    return out


def stress_label(exp, cat, stressed_lower, intensity_map):
    tier = intensity_map.get((exp, cat, stressed_lower), "")
    return f"{tier}_{cat}" if tier else cat


def _rule_factor_type(cfg, category):
    for r in cfg["differential"]["stress_classifier"]:
        if r["category"] == category:
            return r["factor_type"]
    return None


# ------------------------------------------------------------------- species setup
def resolve_taxon(doc, cfg):
    """Gene docs carry a composite taxon_id (e.g. 4558001 = NCBI taxon * 1000 + genome idx).
    Return the NCBI taxon present in the config, or None."""
    raw = doc.get("taxon_id")
    if raw is None:
        return None
    for cand in (int(raw), int(raw) // 1000):
        if str(cand) in cfg["species"]:
            return cand
    return None


def setup_species(cfg, taxon, assays_override=None):
    """Return (sample_to_organ, assays, intensity_map, breakpoints) for a species."""
    base = cfg["api"]["base"]
    if assays_override is not None:
        assays = assays_override
    else:
        exps = api_get(base, cfg["api"]["experiments_endpoint"], {"taxon_id": taxon})
        assays = {}
        for e in exps:
            for a in api_get(base, cfg["api"]["assays_endpoint"], {"experiment": e["_id"]}):
                assays[a["_id"]] = a
    organ_for = make_organ_for(cfg)
    s2o = {}
    for aid, a in assays.items():
        op = factor_of(a, "organism part")
        s2o[aid] = organ_for(op.get("int_id") if op else None, op["label"] if op else None)
    imap = build_intensity_map(assays, cfg)
    bp = cfg["species"][str(taxon)]["reference"]["global_breakpoints_tpm"]
    return s2o, assays, imap, bp


# --------------------------------------------------------------------- scoring
def score_gene(gene_id, doc, taxon, s2o, assays, cfg, intensity_map, breakpoints):
    organs = cfg["organ_panel"]["organs"]
    floor = cfg["expression_level"]["detection_floor_tpm"]
    spec = cfg["specificity"]

    organ_samples = defaultdict(list)
    for field, val in doc.items():
        if val is None:
            continue
        pb = parse_baseline_field(field)
        if not pb:
            continue
        organ = s2o.get(f"{pb[0]}.{pb[1]}")
        if organ:
            organ_samples[organ].append(float(val))

    organ_median = {o: statistics.median(v) for o, v in organ_samples.items()}
    # tau / breadth are computed ONLY over organs this species actually samples, so that
    # species with sparser panels are not biased toward false tissue-specificity by
    # structural zeros. Unsampled organs are 'no_data', distinct from 'not_expressed'.
    sampled = [o for o in organs if o in organ_median]
    vec = [organ_median[o] for o in sampled]
    n_with_data = len(sampled)

    per_organ = {}
    for o in organs:
        if o in organ_median:
            per_organ[o] = {"tpm": round(organ_median[o], 2),
                            "level": level_call(organ_median[o], breakpoints, floor),
                            "n_samples": len(organ_samples[o])}
        else:
            per_organ[o] = {"tpm": None, "level": "no_data", "n_samples": 0}

    tau = compute_tau(vec)
    detected = [o for o in sampled if organ_median[o] >= floor]
    frac_detected = len(detected) / n_with_data if n_with_data else 0
    mx = max(vec) if vec else 0

    if n_with_data < spec["min_organs_required"] or mx < floor:
        breadth = "insufficient_data" if mx < floor else "limited_panel"
        peak = []
    elif tau is not None and tau >= spec["tau_tissue_specific"]:
        breadth = "tissue_specific"
        peak = [o for o in organs if organ_median.get(o, 0) >= 0.5 * mx]
    elif tau is not None and tau < spec["tau_broadly_expressed"] and frac_detected >= spec["broad_min_frac_detected"]:
        breadth = "broadly_expressed"
        peak = []
    else:
        breadth = "tissue_enhanced"
        dv = [organ_median[o] for o in detected] or [0]
        pmed = statistics.median(dv)
        peak = [o for o in organs if organ_median.get(o, 0) >= 0.8 * mx and organ_median.get(o, 0) >= 2 * pmed]

    single_sample_peak = [o for o in peak if len(organ_samples.get(o, [])) <= 1]
    ubiquitous = frac_detected >= spec["ubiquitous_min_frac_detected"]

    # differential
    cut = cfg["differential"]["de_cutoffs"]
    contrasts = defaultdict(dict)
    for field, val in doc.items():
        if val is None:
            continue
        pdf = parse_diff_field(field)
        if pdf:
            exp, gc, gt, metric = pdf
            contrasts[(exp, gc, gt)][metric] = float(val)

    records = []
    genotype_comparisons = 0
    for (exp, gc, gt), m in contrasts.items():
        l2fc, p = m.get("l2fc"), m.get("p")
        if l2fc is None or p is None:
            continue
        ca, ta = assays.get(f"{exp}.{gc}"), assays.get(f"{exp}.{gt}")
        cf, tf = all_factors(ca), all_factors(ta)
        interp = interpret_contrast(cf, tf, l2fc, cfg)
        if not interp:
            genotype_comparisons += 1
            continue
        for cat, stressed, oriented in interp:
            if abs(oriented) < cut["abs_log2fc_min"] or p > cut["p_value_max"]:
                continue
            direction = "activated" if oriented > 0 else "repressed"
            sa = ca if (cf.get(_rule_factor_type(cfg, cat)) or "").lower() == stressed else ta
            tissue = (factor_of(sa, "organism part") or {}).get("label")
            records.append({"cat": cat, "tier_label": stress_label(exp, cat, stressed, intensity_map),
                            "direction": direction,
                            "detail": {"tissue": tissue, "log2fc": round(oriented, 2), "p_value": p,
                                       "direction": direction, "experiment": exp,
                                       "contrast": f"{gc}->{gt}", "stressed_condition": stressed}})

    cat_dirs = defaultdict(set)
    for r in records:
        cat_dirs[r["cat"]].add(r["direction"])
    stress_detail = defaultdict(list)
    activated, repressed = set(), set()
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
        "highly_expressed_in": highly,
        "lowly_expressed_in": lowly,
        "activated_by": sorted(activated),
        "repressed_by": sorted(repressed),
    }
    return {
        "gene_id": gene_id, "taxon_id": taxon,
        "species": cfg["species"][str(taxon)]["short"],
        "name": doc.get("name"), "description": doc.get("description"),
        "attributes": attributes, "breadth": breadth,
        "tau": round(tau, 3) if tau is not None else None,
        "n_organs_with_data": n_with_data, "n_organs_detected": len(detected),
        "frac_organs_detected": round(frac_detected, 2),
        "ubiquitously_detected": ubiquitous, "max_tpm": round(mx, 2),
        "peak_organs": peak, "low_confidence_peak_organs": single_sample_peak,
        "per_organ": per_organ, "stress": {k: v for k, v in stress_detail.items()},
        "n_genotype_comparisons_skipped": genotype_comparisons,
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
    else:
        loc = "ubiquitous" if a["ubiquitously_detected"] else "-"
    return (f"{gid}\t{r['species']}\t{r['breadth']}\t{r['tau']}\t{r['max_tpm']}\t{loc}\t"
            f"{','.join(a['highly_expressed_in'])}\t{','.join(a['activated_by'])}\t{','.join(a['repressed_by'])}")


TSV_HEADER = "gene_id\tspecies\tbreadth\ttau\tmax_tpm\tlocation\thigh_in\tactivated_by\trepressed_by\n"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("genes", nargs="*")
    ap.add_argument("--file")
    ap.add_argument("--config", default=CFG_PATH)
    ap.add_argument("--out")
    ap.add_argument("--tsv", action="store_true")
    args = ap.parse_args()

    cfg = json.load(open(args.config))
    gene_ids = list(args.genes)
    if args.file:
        with open(args.file) as fh:
            gene_ids += [ln.strip() for ln in fh if ln.strip() and not ln.startswith("#")]
    if not gene_ids:
        ap.error("no gene IDs given")

    docs = fetch_expression_docs(cfg, gene_ids)
    # group by taxon, set up each species once
    by_tax = defaultdict(list)
    for gid in gene_ids:
        d = docs.get(gid)
        if d is None:
            continue
        by_tax[resolve_taxon(d, cfg)].append(gid)
    setups = {}
    for tx in by_tax:
        if tx is not None:
            setups[tx] = setup_species(cfg, tx)

    results = {}
    for gid in gene_ids:
        d = docs.get(gid)
        if d is None:
            results[gid] = {"gene_id": gid, "error": "no expression data / gene not found"}
            continue
        tx = resolve_taxon(d, cfg)
        if tx not in setups:
            results[gid] = {"gene_id": gid, "taxon_id": tx, "error": "species not in grass panel"}
            continue
        s2o, assays, imap, bp = setups[tx]
        results[gid] = score_gene(gid, d, tx, s2o, assays, cfg, imap, bp)

    if args.out:
        json.dump(results, open(args.out, "w"), indent=2)
        print(f"wrote {args.out}", file=sys.stderr)
    if args.tsv or not args.out:
        sys.stdout.write(TSV_HEADER)
        for gid in gene_ids:
            r = results[gid]
            print(f"{gid}\t\tERROR\t\t\t\t\t\t{r['error']}" if "error" in r else tsv_row(gid, r))
    if not args.out:
        print("\n--- full detail (JSON) ---")
        print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
