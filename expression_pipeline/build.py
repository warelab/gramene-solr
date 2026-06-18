#!/usr/bin/env python3
"""
build.py — end-to-end driver. For every genome in manifest.json:
  1. discover.py   -> <taxon>.assays_cache.json + inv_<taxon>.json   (experiments/organs/stresses)
  2. refbuild.py   -> ref_<taxon>.json                                (per-species quantile references)
  3. assemble_config.py -> expression_panel.json                      (once, after all references)
  4. batch.py      -> <short>_attributes.jsonl                        (one _id-keyed JSON object per line)

Each step is idempotent / resumable: existing caches and references are skipped, and
batch.py resumes from its JSONL line count. Load the resulting <short>_attributes.jsonl
into MongoDB keyed on _id (gene id), e.g.:

    mongoimport --collection expression_attributes --file sorghum_attributes.jsonl

Usage:
    python3 build.py                      # full build for all genomes in the manifest
    python3 build.py --only 4558 39947    # restrict to specific taxa
    python3 build.py --rediscover         # force re-discovery (ignore cached assays)
    python3 build.py --rerefs             # force reference recomputation

Pure Python standard library; no third-party dependencies.
"""
import argparse, json, os, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))


def run(script, *args):
    cmd = [sys.executable, os.path.join(HERE, script), *map(str, args)]
    print(">>", " ".join(cmd), flush=True)
    subprocess.check_call(cmd, cwd=HERE)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", nargs="*", type=int, help="restrict to these taxa")
    ap.add_argument("--rediscover", action="store_true")
    ap.add_argument("--rerefs", action="store_true")
    args = ap.parse_args()

    import mongo_source
    manifest = json.load(open(os.path.join(HERE, "manifest.json")))
    # The genome list comes from mongo — every genome that actually has expression data in
    # the build's `expression` collection — so coverage tracks the database, not a hand-kept
    # list. Manifest entries act as curation overrides (nice short/name + clade pin).
    genomes = mongo_source.discover_genomes(manifest.get("genomes"))
    if args.only:
        genomes = [g for g in genomes if g["taxon"] in args.only]
    print(">> genomes with expression in mongo: " +
          ", ".join(f"{g['short']}({g['taxon']},{g['clade']})" for g in genomes), flush=True)

    for g in genomes:
        cache = os.path.join(HERE, f"{g['taxon']}.assays_cache.json")
        if args.rediscover or not os.path.exists(cache):
            run("discover.py", g["system_name"], g["taxon"])

    for g in genomes:
        ref = os.path.join(HERE, f"ref_{g['taxon']}.json")
        if args.rerefs or not os.path.exists(ref):
            run("refbuild.py", g["taxon"])

    run("assemble_config.py")

    for g in genomes:
        run("batch.py", "--taxon", g["taxon"])

    print("\nDONE. Per-genome outputs: <short>_attributes.jsonl (load into MongoDB on _id).")


if __name__ == "__main__":
    main()
