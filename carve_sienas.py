#!/usr/bin/env python3
"""
carve_sienas.py
===============
Define SIENAs (Stimulus-Induced ENhancer-locus Annotations) directly from a set
of differential ChIP domains and a gene annotation -- no ChIPseeker required.

A SIENA is an intergenic sub-interval of an induced domain: take each domain,
remove every gene-body span it overlaps, and keep the leftover non-genic pieces.
A domain with N internal genes yields up to N+1 sienas; a domain lying entirely
inside a gene body yields none. Boundaries use the NEAR EDGE of each gene body
(strand-agnostic) and inclusive 1-bp gaps (siena ends at gene_start-1, the next
resumes at gene_end+1).

Each siena inherits its parent domain's statistics (log2FoldChange, FDR, ...),
so the "stimulus-induced" status of the domain carries through. Each siena is
also labelled with a `domain_class`:
    gene_free    parent domain spans 0 genes (fully intergenic domain)
    single_gene  parent domain spans exactly 1 gene
    multi_gene   parent domain spans 2+ genes
single_gene and multi_gene sienas are, by construction, adjacent to a gene body
that lies inside an induced domain (i.e. "near an induced gene body").

Inputs
------
--gtf       Gene annotation (GTF). Gene bodies are reconstructed as the
            exon-union span per gene_id (min exon start -> max exon end,
            introns included). Works with annotations that contain only
            `exon`/`CDS` records (e.g. Liftoff output) -- no `gene` lines needed.
--domains   Differential domain table. Delimiter is auto-detected by default
            (comma or tab); override with --sep. Must contain chromosome/start/
            end plus any stat columns you want carried onto the sienas.

Outputs
-------
--out-csv      Full siena table (CSV), incl. domain_class.
--out-bed      Optional siena BED6 (0-based, half-open) for browsers/bedtools.
--out-genebed  Optional gene-body BED6: full spans of the genes inside
               single_gene + multi_gene domains that yield a qualifying siena.

Example
-------
python3 carve_sienas.py \
    --gtf      Benthi.gtf \
    --domains  H3K9ac_Benthi_2h_vs_ctrl.csv \
    --domain-chrom '#Chromosome' --domain-start Start --domain-end End \
    --carry log2FoldChange FDR PValue Score ChIPCount InputCount \
    --min-log2fc 0.5 --min-len 500 \
    --out-csv      sienas_classified.csv \
    --out-bed      sienas.bed \
    --out-genebed  gene_bodies.bed
"""

import argparse
import re
import sys
from collections import defaultdict

import numpy as np
import pandas as pd

GENE_ID_RE = re.compile(r'gene_id "([^"]+)"')


def build_gene_bodies(gtf_path, feature="exon"):
    """Return {chrom: (starts[], ends[], ids[])} sorted by start.

    Gene body = exon-union span per gene_id: minimum feature start to maximum
    feature end. Introns are included; strand is ignored (near-edge convention).
    """
    span = {}  # gene_id -> [chrom, start, end]
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
            gid, chrom, start, end = m.group(1), f[0], int(f[3]), int(f[4])
            if gid not in span:
                span[gid] = [chrom, start, end]
            else:
                span[gid][1] = min(span[gid][1], start)
                span[gid][2] = max(span[gid][2], end)

    if not span:
        sys.exit(f"ERROR: no '{feature}' records with a gene_id found in {gtf_path}")

    by_chrom = defaultdict(list)
    for gid, (chrom, start, end) in span.items():
        by_chrom[chrom].append((start, end, gid))

    genes = {}
    for chrom, vals in by_chrom.items():
        vals.sort()  # by start
        genes[chrom] = (
            np.array([v[0] for v in vals]),
            np.array([v[1] for v in vals]),
            np.array([v[2] for v in vals], dtype=object),
        )
    return genes, len(span)


def carve_domain(chrom, d_start, d_end, genes):
    """Carve one domain into intergenic sienas.

    Returns (sienas, gene_bodies):
      sienas      -- list of (siena_start, siena_end, left_gene, right_gene,
                     n_genes); left/right_gene is a gene_id or
                     'domain_start'/'domain_end'.
      gene_bodies -- list of (gene_start, gene_end, gene_id) for every gene
                     overlapping the domain (FULL span, not clipped).
    """
    if chrom not in genes:
        # No annotation on this contig: the whole domain is one siena.
        return [(d_start, d_end, "domain_start", "domain_end", 0)], []

    g_start, g_end, g_id = genes[chrom]
    # Genes overlapping the domain: gene_start <= d_end AND gene_end >= d_start
    mask = (g_start <= d_end) & (g_end >= d_start)
    n_genes = int(mask.sum())
    gene_bodies = list(zip(g_start[mask].tolist(),
                           g_end[mask].tolist(),
                           g_id[mask].tolist()))
    if n_genes == 0:
        return [(d_start, d_end, "domain_start", "domain_end", 0)], []

    # Clip overlapping gene bodies to the domain, sort by clipped start.
    cs = np.maximum(g_start[mask], d_start)
    ce = np.minimum(g_end[mask], d_end)
    gid = g_id[mask]
    order = np.argsort(cs)
    cs, ce, gid = cs[order], ce[order], gid[order]

    # Merge touching/overlapping clipped gene intervals into blocks, tracking the
    # gene_id at each block's left and right edge (for boundary labelling).
    blocks = []  # [start, end, left_id, right_id]
    for s, e, gi in zip(cs, ce, gid):
        if blocks and s <= blocks[-1][1] + 1:
            blocks[-1][1] = max(blocks[-1][1], e)
            blocks[-1][3] = gi
        else:
            blocks.append([s, e, gi, gi])

    # Subtract gene blocks from [d_start, d_end]; the gaps are the sienas.
    sienas = []
    cursor = d_start
    prev_gene = "domain_start"
    for b_start, b_end, left_id, right_id in blocks:
        if b_start > cursor:
            sienas.append((cursor, b_start - 1, prev_gene, left_id, n_genes))
        cursor = max(cursor, b_end + 1)
        prev_gene = right_id
    if cursor <= d_end:
        sienas.append((cursor, d_end, prev_gene, "domain_end", n_genes))

    return sienas, gene_bodies


def main():
    ap = argparse.ArgumentParser(
        description="Carve induced ChIP domains into intergenic sienas.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument("--gtf", required=True, help="Gene annotation GTF")
    ap.add_argument("--domains", required=True, help="Differential domain table (CSV or TSV)")
    ap.add_argument("--feature", default="exon",
                    help="GTF feature used to build gene bodies")
    ap.add_argument("--sep", default=None,
                    help="Delimiter for the domain table. Default: auto-detect "
                         "(comma or tab). Use ',' or $'\\t' to force one.")
    ap.add_argument("--domain-chrom", default="#Chromosome",
                    help="Chromosome column name in the domain table")
    ap.add_argument("--domain-start", default="Start",
                    help="Start column name in the domain table")
    ap.add_argument("--domain-end", default="End",
                    help="End column name in the domain table")
    ap.add_argument("--carry", nargs="*",
                    default=["log2FoldChange", "FDR", "PValue", "Score",
                             "ChIPCount", "InputCount"],
                    help="Domain columns to inherit onto each siena")
    # --- domain-level significance / effect-size thresholds (applied BEFORE
    #     carving; sienas inherit the parent domain's stats, so these define
    #     which domains are "differential" enough to contribute sienas) ---
    ap.add_argument("--min-log2fc", type=float, default=None,
                    help="Keep domains with log2FoldChange >= this (signed; "
                         "e.g. 1.0 for >=2-fold gain)")
    ap.add_argument("--min-abs-log2fc", type=float, default=None,
                    help="Keep domains with |log2FoldChange| >= this "
                         "(use for two-sided gains AND losses)")
    ap.add_argument("--max-fdr", type=float, default=None,
                    help="Keep domains with FDR <= this (e.g. 0.05)")
    ap.add_argument("--max-pvalue", type=float, default=None,
                    help="Keep domains with PValue <= this")
    ap.add_argument("--min-score", type=float, default=None,
                    help="Keep domains with Score >= this")
    ap.add_argument("--min-len", type=int, default=0,
                    help="Drop sienas shorter than this many bp (0 = keep all)")
    ap.add_argument("--require-genic-domain", action="store_true",
                    help="Keep only single_gene + multi_gene sienas (those next "
                         "to a gene body inside an induced domain); drops the "
                         "gene_free class from the siena outputs")
    ap.add_argument("--out-csv", required=True, help="Output siena table (CSV)")
    ap.add_argument("--out-bed", default=None,
                    help="Optional BED6 output (0-based, half-open)")
    ap.add_argument("--out-genebed", default=None,
                    help="Optional gene-body BED6 for single_gene + multi_gene "
                         "domains that yield >=1 qualifying siena "
                         "(0-based, half-open; full gene spans)")
    args = ap.parse_args()

    genes, n_genes_total = build_gene_bodies(args.gtf, feature=args.feature)
    print(f"[gtf] genes reconstructed: {n_genes_total:,}", file=sys.stderr)

    # Read the domain table. sep=None + engine="python" auto-detects the
    # delimiter (comma or tab), so a raw epic2 TSV works without reformatting.
    # An explicit --sep overrides the sniffer; accept friendly aliases and an
    # escaped tab ('\t') because a bare literal tab is easily lost by the shell.
    sep = args.sep
    if sep is not None:
        sep = {"tab": "\t", "comma": ",", "\\t": "\t", "\\\\t": "\t"}.get(sep, sep)
    dom = pd.read_csv(args.domains, sep=sep, engine="python")
    if dom.shape[1] == 1:
        sys.exit(
            "ERROR: the domain table parsed into a single column -- the "
            "delimiter was guessed wrong.\n"
            f"  Parsed column: {list(dom.columns)[0]!r}\n"
            "Pass the right delimiter explicitly, e.g. --sep $'\\t' for a "
            "tab-separated epic2 table or --sep ',' for CSV."
        )

    dom = dom.rename(columns={args.domain_chrom: "Chrom",
                              args.domain_start: "domain_start",
                              args.domain_end: "domain_end"})

    # Fail loudly if the chromosome/coordinate columns did not resolve, instead
    # of crashing later inside the carve loop on a missing attribute.
    need = ["Chrom", "domain_start", "domain_end"]
    miss_cols = [c for c in need if c not in dom.columns]
    if miss_cols:
        sys.exit(
            f"ERROR: required column(s) {miss_cols} not found after renaming.\n"
            f"  Columns present: {list(dom.columns)}\n"
            "Check --domain-chrom / --domain-start / --domain-end (and --sep)."
        )

    carry = [c for c in args.carry if c in dom.columns]
    missing = [c for c in args.carry if c not in dom.columns]
    if missing:
        print(f"[warn] carry columns not found, skipping: {missing}", file=sys.stderr)
    print(f"[domains] induced domains: {len(dom):,}", file=sys.stderr)

    # Guard: warn (or stop) if domain chromosomes don't match the GTF, the most
    # common silent cause of "everything is gene_free".
    gtf_chroms = set(genes.keys())
    dom_chroms = set(dom["Chrom"].astype(str).unique())
    overlap = gtf_chroms & dom_chroms
    if not overlap:
        sys.exit(
            "ERROR: no chromosome names are shared between the GTF and the "
            "domain table -- every domain would be 'gene_free'.\n"
            f"  GTF examples   : {sorted(gtf_chroms)[:5]}\n"
            f"  domain examples: {sorted(dom_chroms)[:5]}\n"
            "Rename one side so the chromosome names match."
        )
    missing_chroms = dom_chroms - gtf_chroms
    if missing_chroms:
        print(f"[warn] {len(missing_chroms)} domain chrom(s) absent from the GTF "
              f"(their domains become gene_free): {sorted(missing_chroms)[:10]}",
              file=sys.stderr)

    # Apply domain-level thresholds before carving.
    def _apply(df, col, op, val, label):
        if val is None:
            return df
        if col not in df.columns:
            print(f"[warn] threshold '{label}' set but column '{col}' missing; "
                  f"skipping", file=sys.stderr)
            return df
        n0 = len(df)
        df = df[op(df[col], val)].copy()
        print(f"[filter] {label}: {n0:,} -> {len(df):,} domains", file=sys.stderr)
        return df

    import operator as _op
    dom = _apply(dom, "log2FoldChange", _op.ge, args.min_log2fc, f"log2FC >= {args.min_log2fc}")
    if args.min_abs_log2fc is not None and "log2FoldChange" in dom.columns:
        n0 = len(dom)
        dom = dom[dom["log2FoldChange"].abs() >= args.min_abs_log2fc].copy()
        print(f"[filter] |log2FC| >= {args.min_abs_log2fc}: {n0:,} -> {len(dom):,} domains",
              file=sys.stderr)
    dom = _apply(dom, "FDR", _op.le, args.max_fdr, f"FDR <= {args.max_fdr}")
    dom = _apply(dom, "PValue", _op.le, args.max_pvalue, f"PValue <= {args.max_pvalue}")
    dom = _apply(dom, "Score", _op.ge, args.min_score, f"Score >= {args.min_score}")
    if len(dom) == 0:
        sys.exit("ERROR: no domains pass the thresholds; nothing to carve.")

    rows = []
    gene_rows = []   # (Chrom, gene_start, gene_end, gene_id, domain_start, domain_end, n_genes)
    domains_with_siena = 0
    for r in dom.itertuples(index=False):
        c = getattr(r, "Chrom")
        ds = int(getattr(r, "domain_start"))
        de = int(getattr(r, "domain_end"))
        sien, gbodies = carve_domain(c, ds, de, genes)
        if sien:
            domains_with_siena += 1
        n_in_dom = len(sien)
        for k, (a, b, lg, rg, ng) in enumerate(sien, start=1):
            rows.append((c, ds, de, a, b, lg, rg, ng, k, n_in_dom))
        for gs, ge, gid in gbodies:
            gene_rows.append((c, gs, ge, gid, ds, de, len(gbodies)))

    cols = ["Chrom", "domain_start", "domain_end", "siena_start", "siena_end",
            "left_gene", "right_gene", "n_genes_in_domain",
            "siena_idx_in_domain", "n_sienas_in_domain"]
    S = pd.DataFrame(rows, columns=cols)
    S["siena_len"] = S.siena_end - S.siena_start + 1

    # domain_class from gene count: 0 -> gene_free, 1 -> single_gene, >=2 -> multi_gene
    def _dom_class(n):
        return "gene_free" if n == 0 else ("single_gene" if n == 1 else "multi_gene")
    S["domain_class"] = S.n_genes_in_domain.map(_dom_class)

    # Inherit parent-domain statistics.
    if carry:
        S = S.merge(dom[["Chrom", "domain_start", "domain_end", *carry]],
                    on=["Chrom", "domain_start", "domain_end"], how="left")

    # Optional length filter.
    n_before = len(S)
    if args.min_len > 0:
        S = S[S.siena_len >= args.min_len].copy()

    # Optional: keep only sienas next to an induced gene body (single/multi gene).
    if args.require_genic_domain:
        n_pre = len(S)
        S = S[S.domain_class.isin(["single_gene", "multi_gene"])].copy()
        print(f"[filter] require-genic-domain: {n_pre:,} -> {len(S):,} sienas "
              f"(dropped gene_free)", file=sys.stderr)

    # Stable IDs after sorting by genomic position.
    S = S.sort_values(["Chrom", "siena_start", "siena_end"]).reset_index(drop=True)
    S.insert(0, "siena_id", [f"siena_{i:05d}" for i in range(1, len(S) + 1)])

    S.to_csv(args.out_csv, index=False)

    if args.out_bed:
        # BED is 0-based, half-open: start-1, end stays as-is.
        bed = pd.DataFrame({
            "chrom": S.Chrom,
            "start": S.siena_start - 1,
            "end": S.siena_end,
            "name": S.siena_id,
            "score": (S["log2FoldChange"] if "log2FoldChange" in S.columns else 0),
            "strand": ".",
        })
        bed.to_csv(args.out_bed, sep="\t", header=False, index=False)

    if args.out_genebed:
        # Gene bodies belonging to single_gene + multi_gene domains that still
        # have >=1 siena after all thresholds. Keyed on the surviving sienas in S.
        keep_dom = S.loc[S.domain_class.isin(["single_gene", "multi_gene"]),
                         ["Chrom", "domain_start", "domain_end"]].drop_duplicates()
        GB = pd.DataFrame(gene_rows, columns=["Chrom", "gene_start", "gene_end",
                                              "gene_id", "domain_start",
                                              "domain_end", "n_genes_in_domain"])
        GB = GB.merge(keep_dom, on=["Chrom", "domain_start", "domain_end"], how="inner")
        GB = GB.drop_duplicates(["Chrom", "gene_start", "gene_end", "gene_id"]) \
               .sort_values(["Chrom", "gene_start", "gene_end"])
        gene_bed = pd.DataFrame({
            "chrom": GB.Chrom,
            "start": GB.gene_start - 1,      # 1-based inclusive -> 0-based half-open
            "end": GB.gene_end,
            "name": GB.gene_id,
            "score": 0,
            "strand": ".",
        })
        gene_bed.to_csv(args.out_genebed, sep="\t", header=False, index=False)
        print(f"[result] gene bodies in genebed     : {len(gene_bed):,} "
              f"(single+multi gene domains with a qualifying siena)", file=sys.stderr)

    # Console summary.
    print(f"[result] domains yielding >=1 siena : {domains_with_siena:,}", file=sys.stderr)
    print(f"[result] domains yielding 0 sienas  : {len(dom) - domains_with_siena:,} "
          f"(wholly inside a gene body)", file=sys.stderr)
    print(f"[result] sienas before length filter: {n_before:,}", file=sys.stderr)
    print(f"[result] sienas written             : {len(S):,}"
          + (f"  (>= {args.min_len} bp)" if args.min_len else ""), file=sys.stderr)
    print(f"[result] wrote {args.out_csv}"
          + (f" and {args.out_bed}" if args.out_bed else ""), file=sys.stderr)


if __name__ == "__main__":
    main()
