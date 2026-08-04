#!/usr/bin/env python3
"""
extract_maker_scores.py <gff> <stable_ids.out> <genome> <outdir>

Pulls MAKER quality scores out of one pan-gene coding GFF and emits the two artifacts the
release-10 MAKER pipeline consumes:

  <outdir>/<genome>.AED_QI.scores   transcript_id \t _AED=<v> \t _QI=<9 pipe-separated>
                                    every mRNA, matching the shape of the existing
                                    28_sorg_maker_AED/*.AED_QI.scores files
  stdout                            <new stable gene id> \t <canonical transcript id>
                                    appended by the driver into geneID.canonical.txt

Why the id mapping is mandatory: these GFFs are pan-gene sets, so their gene ids come from
whichever source annotation each gene was drawn from (SbPI513676…, IS12661…, Sobic.…,
3381.casb…) and do NOT match the assembly they belong to. assign_stable_ids_*.out maps
old gene id (col2) -> new stable id (col3), and the new id is what the solr index uses.

Canonical transcript = longest summed CDS, ties broken by longest summed exon length, then
by transcript id so repeat runs are byte-identical.
"""
import sys, os, re


def attrs(col9):
    """GFF3 column 9 -> dict. Values are left as-is (no URL-unescaping needed here)."""
    d = {}
    for field in col9.rstrip().split(';'):
        if not field:
            continue
        k, _, v = field.partition('=')
        d[k] = v
    return d


def main():
    if len(sys.argv) != 5:
        sys.exit(__doc__)
    gff, idmap_path, genome, outdir = sys.argv[1:5]

    # old gene id -> new stable gene id, plus a (chr,start,end) -> new id fallback.
    #
    # The id join is authoritative but incomplete: these are pan-gene sets, and the GFF and the
    # stable-id assignment do not always pick the same source annotation as a gene's representative
    # id. For sorghum_pi655993 only 14,387 of 46,077 genes join by id, yet every one that does has
    # byte-identical coordinates, and joining on coordinates instead recovers 46,061. So: try the id
    # first, fall back to position. Coordinate keys are near-unique (15-23 collisions per file); the
    # ambiguous ones are dropped rather than guessed.
    old2new = {}
    coord_hits = {}
    with open(idmap_path) as fh:
        for line in fh:
            c = line.rstrip('\n').split('\t')
            if len(c) >= 6 and c[1] and c[2]:
                old2new[c[1]] = c[2]
                coord_hits.setdefault((c[3], c[4], c[5]), []).append(c[2])
    coord2new = {k: v[0] for k, v in coord_hits.items() if len(v) == 1}

    mrna = {}          # tid -> {"gene": old gene id, "aed": str, "qi": str}
    exon_len = {}      # tid -> summed exon length
    cds_len = {}       # tid -> summed CDS length
    gene_coord_hits = {}   # (chr,start,end) -> [gene ids]  (to drop ambiguous GFF coords)
    gene_coord = {}        # gene id -> (chr,start,end)

    with open(gff) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            c = line.rstrip('\n').split('\t')
            if len(c) < 9:
                continue
            kind = c[2]
            if kind == 'gene':
                a = attrs(c[8])
                gid = a.get('ID')
                if gid:
                    key = (c[0], c[3], c[4])
                    gene_coord[gid] = key
                    gene_coord_hits.setdefault(key, []).append(gid)
            elif kind == 'mRNA':
                a = attrs(c[8])
                tid = a.get('ID')
                if not tid:
                    continue
                mrna[tid] = {'gene': a.get('Parent', ''),
                             'aed': a.get('_AED', ''),
                             'qi': a.get('_QI', '')}
            elif kind in ('exon', 'CDS'):
                a = attrs(c[8])
                parent = a.get('Parent', '')
                if not parent:
                    continue
                try:
                    length = int(c[4]) - int(c[3]) + 1
                except ValueError:
                    continue
                # a feature may legitimately list several parents
                for p in parent.split(','):
                    tgt = cds_len if kind == 'CDS' else exon_len
                    tgt[p] = tgt.get(p, 0) + length

    # ---- scores file: every mRNA, in the 3-column release-10 shape --------------------
    os.makedirs(outdir, exist_ok=True)
    scores_path = os.path.join(outdir, '%s.AED_QI.scores' % genome)
    n_scores = 0
    with open(scores_path, 'w') as out:
        for tid in sorted(mrna):
            m = mrna[tid]
            if not m['aed'] or not m['qi']:
                continue
            out.write('%s\t_AED=%s\t_QI=%s\n' % (tid, m['aed'], m['qi']))
            n_scores += 1

    # ---- canonical transcript per gene ------------------------------------------------
    by_gene = {}
    for tid, m in mrna.items():
        by_gene.setdefault(m['gene'], []).append(tid)

    n_canon = n_unmapped = n_by_id = n_by_coord = 0
    for gene_old in sorted(by_gene):
        new_id = old2new.get(gene_old)
        if new_id:
            n_by_id += 1
        else:
            key = gene_coord.get(gene_old)
            # only fall back when the position is unambiguous on BOTH sides
            if key and len(gene_coord_hits.get(key, [])) == 1:
                new_id = coord2new.get(key)
                if new_id:
                    n_by_coord += 1
        if not new_id:
            n_unmapped += 1
            continue
        # longest CDS, then longest transcript, then id (deterministic)
        best = sorted(by_gene[gene_old],
                      key=lambda t: (-cds_len.get(t, 0), -exon_len.get(t, 0), t))[0]
        if not mrna[best]['aed'] or not mrna[best]['qi']:
            continue
        sys.stdout.write('%s\t%s\n' % (new_id, best))
        n_canon += 1

    sys.stderr.write(
        '  %-26s mRNA=%-7d genes=%-6d scores=%-7d canonical=%-6d '
        '(by_id=%d by_coord=%d) unmapped=%d\n'
        % (genome, len(mrna), len(by_gene), n_scores, n_canon,
           n_by_id, n_by_coord, n_unmapped))


if __name__ == '__main__':
    main()
