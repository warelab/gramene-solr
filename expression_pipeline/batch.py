#!/usr/bin/env python3
"""
all_species_batch.py — score every expression gene of one genome, streaming to
<short>_attributes.{jsonl,tsv}. Resumable + time-budgeted (invoke repeatedly).

  python3 all_species_batch.py --taxon 4565            # wheat
  python3 all_species_batch.py --taxon 4565 --assemble # build .json from .jsonl
"""
import argparse, json, os, time
import score as A
import mongo_source


def assemble(prefix):
    with open(prefix + ".jsonl") as src, open(prefix + ".json", "w") as out:
        out.write("{\n"); first = True
        for line in src:
            line = line.strip()
            if not line:
                continue
            rec = json.loads(line)
            if not first:
                out.write(",\n")
            first = False
            out.write(json.dumps(rec["_id"]) + ": " + json.dumps(rec))
        out.write("\n}\n")
    print("assembled", prefix + ".json")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--taxon", type=int, required=True)
    ap.add_argument("--config", default=A.CFG_PATH)
    ap.add_argument("--rows", type=int, default=400)
    ap.add_argument("--budget", type=float, default=1e9)
    ap.add_argument("--assemble", action="store_true")
    args = ap.parse_args()

    cfg = json.load(open(args.config))
    sp = cfg["species"][str(args.taxon)]
    short, system_name = sp["short"], sp["system_name"]
    prefix = f"{short}_attributes"
    if args.assemble:
        assemble(prefix); return

    organs, s2o, assays, imap, bp, clade = A.setup_species(cfg, args.taxon)

    jpath, tpath = prefix + ".jsonl", prefix + ".tsv"
    done = sum(1 for _ in open(jpath)) if os.path.exists(jpath) else 0
    if not os.path.exists(tpath) or os.path.getsize(tpath) == 0:
        open(tpath, "w").write(A.TSV_HEADER)

    jf, tf = open(jpath, "a"), open(tpath, "a")
    t0 = time.time()
    # stream every expression-bearing gene of the genome directly from mongo, resuming
    # past the `done` already written (budget lets a sandboxed run stop + resume).
    for doc in mongo_source.iter_genome(system_name, skip=done):
        if time.time() - t0 > args.budget:
            break
        rec = A.score_gene(doc["id"], doc, args.taxon, organs, s2o, assays, cfg, imap, bp, clade)
        jf.write(json.dumps(rec) + "\n")
        tf.write(A.tsv_row(doc["id"], rec) + "\n")
        done += 1
        if done % 2000 == 0:
            jf.flush(); tf.flush()
            print(f"[{time.strftime('%H:%M:%S')}] {short} done={done}", flush=True)
    jf.flush(); tf.flush(); jf.close(); tf.close()
    print(f"COMPLETE: {short} {done}", flush=True)


if __name__ == "__main__":
    main()
