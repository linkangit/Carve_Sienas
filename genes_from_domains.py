#!/usr/bin/env python3
"""
genes_from_domains.py
=====================
List the genes that FLANK a qualifying siena (>= --min-siena-len bp), organised
into gene CLUSTERS, with optional fusion of nearby clusters.

Pipeline
--------
1. Keep domains with  log2FoldChange > --min-log2fc.
2. Carve each surviving domain into intergenic sienas (interval subtraction
   against gene bodies), exactly as carve_sienas.py does.
3. Keep sienas with length >= --min-siena-len ("valid" sienas).
4. Report only the genes that FLANK a valid siena (its left_gene / right_gene,
   excluding domain edges). Genes boxed in by sub-threshold sienas are dropped.
5. CLUSTER: the flanking genes of a domain are strung into one ordered cluster.
   With --fuse-gap N, contributing domains within N bp of each other (same
   chromosome, single-linkage chaining) are fused so their genes share one
   cluster. Without --fuse-gap, each domain is its own cluster.
       >=2 genes in the (fused) cluster -> multi_gene, labelled cluster_1, ...
       1  gene                          -> single_gene (no cluster label)

Columns
-------
domain_id    domain_1, domain_2, ... -- the parent epic domain of each gene
             (assigned to every contributing domain, in genomic order).
domain_class single_gene / multi_gene, following the (fused) CLUSTER size --
             NOT the raw count of genes overlapping the domain.
cluster_id   cluster_1, cluster_2, ... for clusters of >=2 genes.
cluster_size / cluster_order   size of the cluster and the gene's 5'->3' rank in it.
n_domains_in_cluster           how many epic domains the cluster spans (>1 = fused).

Gene bodies are built as in carve_sienas.py: feature-union span per gene_id
(--feature, default exon), optional --strip-suffix to collapse isoforms, strand
from GTF column 7.

Usage
-----
python3 genes_from_domains.py \
    --gtf annotation.gtf --domains diff_domains.csv \
    --feature exon --min-log2fc 1 --min-siena-len 1000 --fuse-gap 10000 \
    --out gene_clusters.tsv
"""
import argparse, re, sys
from collections import defaultdict, Counter
import numpy as np
import pandas as pd

GENE_ID_RE = re.compile(r'gene_id "([^"]+)"')


def build_gene_bodies(gtf_path, feature="exon", strip_suffix=None):
    span, strands = {}, {}
    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != feature:
                continue
            m = GENE_ID_RE.search(f[8])
            if not m:
                continue
            gid = m.group(1)
            if strip_suffix:
                gid = gid.split(strip_suffix)[0]
            chrom, s, e, st = f[0], int(f[3]), int(f[4]), f[6]
            if gid not in span:
                span[gid] = [chrom, s, e]; strands[gid] = {st}
            else:
                span[gid][1] = min(span[gid][1], s)
                span[gid][2] = max(span[gid][2], e)
                strands[gid].add(st)
    if not span:
        sys.exit(f"ERROR: no '{feature}' records with a gene_id found in {gtf_path}")
    strand_by_id, body_by_id = {}, {}
    for gid, sset in strands.items():
        sset = {x for x in sset if x in ("+", "-")}
        strand_by_id[gid] = sset.pop() if len(sset) == 1 else "."
    for gid, (chrom, s, e) in span.items():
        body_by_id[gid] = (chrom, s, e)
    by_chrom = defaultdict(list)
    for gid, (chrom, s, e) in span.items():
        by_chrom[chrom].append((s, e, gid))
    genes = {}
    for chrom, vals in by_chrom.items():
        vals.sort()
        genes[chrom] = (np.array([v[0] for v in vals]),
                        np.array([v[1] for v in vals]),
                        np.array([v[2] for v in vals], dtype=object))
    return genes, strand_by_id, body_by_id, len(span)


def carve_domain(chrom, d_start, d_end, genes):
    """Return (sienas, n_overlap). sienas = (start, end, left_gene, right_gene)."""
    if chrom not in genes:
        return [(d_start, d_end, "domain_start", "domain_end")], 0
    g_start, g_end, g_id = genes[chrom]
    mask = (g_start <= d_end) & (g_end >= d_start)
    n_overlap = int(mask.sum())
    if n_overlap == 0:
        return [(d_start, d_end, "domain_start", "domain_end")], 0
    cs = np.maximum(g_start[mask], d_start)
    ce = np.minimum(g_end[mask], d_end)
    gid = g_id[mask]
    order = np.argsort(cs)
    cs, ce, gid = cs[order], ce[order], gid[order]
    blocks = []
    for s, e, gi in zip(cs, ce, gid):
        if blocks and s <= blocks[-1][1] + 1:
            blocks[-1][1] = max(blocks[-1][1], e); blocks[-1][3] = gi
        else:
            blocks.append([s, e, gi, gi])
    sienas, cursor, prev = [], d_start, "domain_start"
    for bs, be, lid, rid in blocks:
        if bs > cursor:
            sienas.append((cursor, bs - 1, prev, lid))
        cursor = max(cursor, be + 1); prev = rid
    if cursor <= d_end:
        sienas.append((cursor, d_end, prev, "domain_end"))
    return sienas, n_overlap


def main():
    ap = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--domains", required=True)
    ap.add_argument("--feature", default="exon")
    ap.add_argument("--strip-suffix", default=None)
    ap.add_argument("--sep", default=None)
    ap.add_argument("--domain-chrom", default="#Chromosome")
    ap.add_argument("--domain-start", default="Start")
    ap.add_argument("--domain-end", default="End")
    ap.add_argument("--min-log2fc", type=float, default=1.0,
                    help="Keep domains with log2FoldChange STRICTLY GREATER than this")
    ap.add_argument("--min-siena-len", type=int, default=1000,
                    help="Keep sienas with length >= this (bp); a gene is reported "
                         "only if such a siena flanks it")
    ap.add_argument("--fuse-gap", type=int, default=None,
                    help="Fuse gene clusters from DIFFERENT epic domains lying "
                         "within this many bp of each other (same chromosome, "
                         "single-linkage chaining) into one cluster. Default: no "
                         "fusion (one cluster per domain).")
    ap.add_argument("--out", required=True, help="Output TSV (one row per gene)")
    args = ap.parse_args()

    genes, strand_by_id, body_by_id, n_units = build_gene_bodies(
        args.gtf, feature=args.feature, strip_suffix=args.strip_suffix)
    print(f"[gtf] gene bodies: {n_units:,}", file=sys.stderr)

    sep = args.sep
    if sep is not None:
        sep = {"tab": "\t", "comma": ",", "\\t": "\t"}.get(sep, sep)
    dom = pd.read_csv(args.domains, sep=sep, engine="python")
    if dom.shape[1] == 1:
        sys.exit("ERROR: domain table parsed to one column -- wrong delimiter; pass --sep.")
    dom = dom.rename(columns={args.domain_chrom: "Chrom",
                              args.domain_start: "domain_start",
                              args.domain_end: "domain_end"})
    for c in ["Chrom", "domain_start", "domain_end", "log2FoldChange"]:
        if c not in dom.columns:
            sys.exit(f"ERROR: column '{c}' not found. Present: {list(dom.columns)}")
    print(f"[domains] total: {len(dom):,}", file=sys.stderr)

    n0 = len(dom)
    dom = dom[dom.log2FoldChange > args.min_log2fc].copy()
    print(f"[filter] log2FC > {args.min_log2fc}: {n0:,} -> {len(dom):,} domains", file=sys.stderr)
    if len(dom) == 0:
        sys.exit("ERROR: no domains pass the log2FC filter.")

    dom = dom.sort_values(["Chrom", "domain_start", "domain_end"]).reset_index(drop=True)

    # ---- pass 1: per-domain flanking genes of valid (>= min_siena_len) sienas ----
    contrib = []
    n_genefree = n_no_qual = 0
    for r in dom.itertuples(index=False):
        c, ds, de = r.Chrom, int(r.domain_start), int(r.domain_end)
        sien, n_overlap = carve_domain(c, ds, de, genes)
        if n_overlap == 0:
            n_genefree += 1
            continue
        flank = {}
        for (a, b, lg, rg) in sien:
            if (b - a + 1) < args.min_siena_len:
                continue
            for g in (lg, rg):
                if g not in ("domain_start", "domain_end"):
                    flank[g] = flank.get(g, 0) + 1
        if not flank:
            n_no_qual += 1
            continue
        contrib.append({"chrom": c, "ds": ds, "de": de,
                        "log2fc": r.log2FoldChange, "genes": flank})

    contrib.sort(key=lambda d: (d["chrom"], d["ds"], d["de"]))
    for i, d in enumerate(contrib, start=1):
        d["domain_id"] = f"domain_{i}"

    # ---- pass 2: fuse contributing domains within --fuse-gap into groups ----
    # (single-linkage chaining; gap measured between consecutive contributing
    #  domains on the same chromosome). Without --fuse-gap, each domain is its
    #  own group.
    groups = []
    for d in contrib:
        if (groups and groups[-1]["chrom"] == d["chrom"]
                and args.fuse_gap is not None
                and (d["ds"] - groups[-1]["max_end"] - 1) <= args.fuse_gap):
            groups[-1]["domains"].append(d)
            groups[-1]["max_end"] = max(groups[-1]["max_end"], d["de"])
        else:
            groups.append({"chrom": d["chrom"], "max_end": d["de"], "domains": [d]})

    # ---- pass 3: classify (fused) clusters and emit one row per gene ----
    rows = []
    cluster_counter = 0
    n_single = n_multi = n_fused = 0
    cluster_sizes = []
    for grp in groups:
        gene_count, gene_domain = {}, {}
        for d in grp["domains"]:
            for g, cnt in d["genes"].items():
                gene_count[g] = gene_count.get(g, 0) + cnt
                gene_domain.setdefault(g, d)      # leftmost domain that has the gene
        genes_ordered = sorted(gene_count, key=lambda g: body_by_id[g][1])
        csize = len(genes_ordered)
        ndoms = len(grp["domains"])
        if ndoms > 1:
            n_fused += 1
        if csize >= 2:
            cluster_counter += 1
            clust_id, dom_class = f"cluster_{cluster_counter}", "multi_gene"
            cluster_sizes.append(csize)
            n_multi += 1
        else:
            clust_id, dom_class = "", "single_gene"
            n_single += 1
        for pos, g in enumerate(genes_ordered, start=1):
            d = gene_domain[g]
            _, gs, ge = body_by_id[g]
            rows.append({
                "domain_id": d["domain_id"],
                "domain_class": dom_class,
                "cluster_id": clust_id,
                "cluster_size": csize,
                "cluster_order": pos if clust_id else "",
                "n_domains_in_cluster": ndoms,
                "Chrom": d["chrom"], "domain_start": d["ds"], "domain_end": d["de"],
                "log2FoldChange": d["log2fc"],
                "gene": g, "gene_strand": strand_by_id.get(g, "."),
                "gene_start": gs, "gene_end": ge,
                "n_flanking_qual_sienas": gene_count[g],
            })

    out = pd.DataFrame(rows).sort_values(["Chrom", "gene_start"])
    out.to_csv(args.out, sep="\t", index=False)

    print(f"[result] single_gene (1 flanking gene)    : {n_single:,}", file=sys.stderr)
    print(f"[result] multi_gene clusters (>=2 genes)  : {n_multi:,}", file=sys.stderr)
    if args.fuse_gap is not None:
        print(f"[result] clusters fused across >1 domain  : {n_fused:,} "
              f"(--fuse-gap {args.fuse_gap:,})", file=sys.stderr)
    if cluster_sizes:
        dist = ", ".join(f"{k}-gene x{v}" for k, v in sorted(Counter(cluster_sizes).items()))
        print(f"[result]   cluster size distribution     : {dist}", file=sys.stderr)
    print(f"[result] domains w/ genes but no >= {args.min_siena_len} bp siena: {n_no_qual:,}", file=sys.stderr)
    print(f"[result] gene_free domains (dropped)      : {n_genefree:,}", file=sys.stderr)
    print(f"[result] genes written: {len(out):,} ({out.gene.nunique():,} unique) -> {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
