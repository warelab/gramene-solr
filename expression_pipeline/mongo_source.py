#!/usr/bin/env python3
"""
mongo_source.py — data access for the expression-attribute pipeline, reading DIRECTLY
from the build's MongoDB collections instead of the HTTP data API. Pure standard library:
it shells out to `mongosh` (no pymongo dependency), preserving the pipeline's no-deps design.

Connection comes from the environment (the build's config.sh exports these):
    MONGO_URI   default mongodb://localhost:27017
    MONGO_DB    default sorghum11

Collections consumed:
  expression : keyed by gene _id; each non-_id key is a GXA experiment accession ->
               [{group,value} baseline | {group:"gc_gt", l2fc, p_value} differential]
  assays     : _id="<exp>.<group>", experiment, group, taxon_id (NCBI base),
               factor[]/characteristic[] with {type,label,int_id}
  experiments: _id, type ("Baseline"/"Differential"), taxon_id
  genes      : _id, system_name, taxon_id (NCBI*1000+idx), name, description

The scorer consumes a flat "doc" (the same shape the solr /search API returned):
  {id, name, description, taxon_id,
   "<exp>_<grp>__expr": tpm,
   "<exp>_<gc>_<gt>_l2fc_attr_f": l2fc, "<exp>_<gc>_<gt>_pval_attr_f": p}
flatten_expr() builds those fields from a mongo expression document (mongo2solr does the
exact same flattening when it writes the solr core, so the two stay consistent).
"""
import json, os, subprocess, tempfile

MONGO_URI = os.environ.get("MONGO_URI", "mongodb://localhost:27017")
MONGO_DB = os.environ.get("MONGO_DB", "sorghum11")
# where to spool mongosh output (see _stream); defaults next to this file (on /usr/local,
# ample) rather than /tmp which may be a small root fs.
SPOOL_DIR = os.environ.get("MONGO_SPOOL_DIR") or os.path.dirname(os.path.abspath(__file__))


def _target():
    return MONGO_URI.rstrip("/") + "/" + MONGO_DB


def _stream(js):
    """Run a mongosh --eval snippet that print()s one JSON object per line; yield parsed objects.

    mongosh's stdout is SPOOLED to a temp file first, then parsed. This is essential: if we read
    the pipe lazily while doing slow per-doc work (scoring), the OS pipe back-pressures mongosh,
    its server-side cursor goes idle and gets killed, and the stream is silently truncated — worst
    for genomes with large joined docs (e.g. maize: ~34 KB/doc x 44k genes). Spooling lets mongosh
    run to completion at disk speed regardless of consumer speed. We also CHECK the exit status and
    raise on failure instead of swallowing it (the old code sent stderr to /dev/null)."""
    fd, path = tempfile.mkstemp(suffix=".ndjson", dir=SPOOL_DIR)
    try:
        with os.fdopen(fd, "wb") as out:
            r = subprocess.run(["mongosh", "--quiet", _target(), "--eval", js],
                               stdout=out, stderr=subprocess.PIPE)
        if r.returncode != 0:
            err = (r.stderr or b"").decode("utf-8", "replace").strip()
            raise RuntimeError("mongosh failed (rc=%d): %s" % (r.returncode, err[-2000:]))
        with open(path, "r", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                line = line.strip()
                if not line or line[0] not in "{[":
                    continue
                try:
                    yield json.loads(line)
                except ValueError:
                    continue
    finally:
        try:
            os.unlink(path)
        except OSError:
            pass


# ------------------------------------------------------------------- expr flattening
def flatten_expr(expr):
    """mongo expression doc -> flat solr-style field dict that the scorer parses."""
    out = {}
    if not expr:
        return out
    for exp_id, samples in expr.items():
        if exp_id == "_id" or not isinstance(samples, list):
            continue
        key = exp_id.replace("-", "_")
        for s in samples:
            grp = s.get("group")
            if grp is None:
                continue
            if s.get("value") is not None:
                out["%s_%s__expr" % (key, grp)] = s["value"]
            if s.get("l2fc") is not None:
                out["%s_%s_l2fc_attr_f" % (key, grp)] = s["l2fc"]
            if s.get("p_value") is not None:
                out["%s_%s_pval_attr_f" % (key, grp)] = s["p_value"]
    return out


def _doc_from(g):
    """{_id,name,description,taxon_id,expr} (genes ⋈ expression) -> scorer doc."""
    d = {"id": g["_id"], "name": g.get("name"),
         "description": g.get("description"), "taxon_id": g.get("taxon_id")}
    d.update(flatten_expr(g.get("expr")))
    return d


# --------------------------------------------------------------------- gene docs
def fetch_docs(gene_ids):
    """{gene_id: scorer_doc} for a list of gene ids (ad-hoc lookups / score.py)."""
    ids = [g for g in gene_ids if g]
    if not ids:
        return {}
    js = ("var ids=%s;"
          "db.genes.find({_id:{$in:ids}},{name:1,description:1,taxon_id:1}).forEach(function(g){"
          "  g.expr=db.expression.findOne({_id:g._id});"
          "  print(JSON.stringify(g));"
          "});") % json.dumps(ids)
    return {g["_id"]: _doc_from(g) for g in _stream(js)}


def iter_genome(system_name, skip=0):
    """Stream scorer docs for every expression-bearing gene of a genome, sorted by _id.
    Joins genes(system_name) ⋈ expression server-side (like mongo2solr.join.js). `skip`
    supports resume. Yields one flat doc per gene."""
    js = ("db.genes.aggregate(["
          "  {$match:{system_name:%s}},"
          "  {$lookup:{from:'expression',localField:'_id',foreignField:'_id',as:'expr'}},"
          "  {$match:{'expr.0':{$exists:true}}},"
          "  {$sort:{_id:1}},"
          "  {$skip:%d},"
          "  {$project:{name:1,description:1,taxon_id:1,expr:{$first:'$expr'}}}"
          "],{allowDiskUse:true}).forEach(function(g){print(JSON.stringify(g));});"
          ) % (json.dumps(system_name), int(skip))
    for g in _stream(js):
        yield _doc_from(g)


# --------------------------------------------------------------------- assays / experiments
def load_assays(taxon):
    """{assay_id: assay} for ALL assays of the experiments this taxon participates in
    (taxon = NCBI base id). Same shape as the old discovery cache."""
    js = ("var tx=%d;"
          "var exps=db.assays.distinct('experiment',{taxon_id:tx});"
          "db.assays.find({experiment:{$in:exps}}).forEach(function(a){print(JSON.stringify(a));});"
          ) % int(taxon)
    return {a["_id"]: a for a in _stream(js)}


def experiments_for_taxon(taxon):
    """{exp_id: experiment} for the experiments this taxon participates in."""
    js = ("var tx=%d;"
          "var exps=db.assays.distinct('experiment',{taxon_id:tx});"
          "db.experiments.find({_id:{$in:exps}}).forEach(function(e){print(JSON.stringify(e));});"
          ) % int(taxon)
    return {e["_id"]: e for e in _stream(js)}


# --------------------------------------------------------------------- genome discovery
_GENOMES = None


def discover_genomes(manifest_overrides=None):
    """Genomes that actually have expression data in mongo (so coverage tracks the build,
    not a hand-maintained list): [{system_name,taxon,clade,short,name,n_genes}], sorted by size.
    taxon = NCBI base id (taxon_id // 1000). clade comes from the manifest override if present,
    else po_map.CLADE_OF_TAXON, else 'other' (stress-only)."""
    global _GENOMES
    if _GENOMES is not None and manifest_overrides is None:
        return _GENOMES
    import po_map as M
    ov = {g["system_name"]: g for g in (manifest_overrides or [])}
    js = ("db.expression.aggregate(["
          "  {$lookup:{from:'genes',localField:'_id',foreignField:'_id',as:'g'}},"
          "  {$unwind:'$g'},"
          "  {$group:{_id:{sn:'$g.system_name',tx:'$g.taxon_id'},n:{$sum:1}}}"
          "],{allowDiskUse:true}).forEach(function(d){"
          "  print(JSON.stringify({system_name:d._id.sn,taxon_id:d._id.tx,n:d.n}));});")
    genomes = []
    for d in _stream(js):
        sn = d["system_name"]
        base = int(d["taxon_id"]) // 1000
        o = ov.get(sn, {})
        genomes.append({
            "system_name": sn, "taxon": base,
            "clade": o.get("clade") or M.CLADE_OF_TAXON.get(base, "other"),
            "short": o.get("short", sn), "name": o.get("name", sn),
            "n_genes": d["n"]})
    genomes.sort(key=lambda x: -x["n_genes"])
    if manifest_overrides is None:
        _GENOMES = genomes
    return genomes


def system_name_for_taxon(taxon):
    for g in discover_genomes():
        if g["taxon"] == int(taxon):
            return g["system_name"]
    return None


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == "genomes":
        for g in discover_genomes():
            print("%-32s taxon=%-7d clade=%-10s genes_with_expr=%d"
                  % (g["system_name"], g["taxon"], g["clade"], g["n_genes"]))
    else:
        print("usage: mongo_source.py genomes", file=sys.stderr)
