#!/usr/bin/env python3
"""
carve_sienas.py
===============
Define SIENAs (Stimulus-Induced ENhancer-locus Annotations) directly from a set
of differential ChIP domains and a gene annotation -- no ChIPseeker required.

A SIENA is an intergenic sub-interval of an induced domain: take each domain,
remove every gene-body span it overlaps, and keep the leftover non-genic pieces.
A domain with N internal genes yields N-1 to N+1 sienas (depending on whether
each domain edge lands in intergenic space); a domain lying entirely inside a
gene body yields none. Boundaries use the NEAR EDGE of each gene body
(strand-agnostic) and inclusive 1-bp gaps (siena ends at gene_start-1, the next
resumes at gene_end+1).

Each siena inherits its parent domain's statistics (log2FoldChange, FDR, ...),
so the "stimulus-induced" status of the domain carries through.

Two labels are attached to every siena:

`domain_class` -- how many genes the PARENT DOMAIN spans:
    gene_free    parent domain spans 0 genes (fully intergenic domain)
    single_gene  parent domain spans exactly 1 gene
    multi_gene   parent domain spans 2+ genes

`siena_class` -- refines that by the siena's PROMOTER SIDE (strand-aware), so you
can see at a glance which sienas are 5'/promoter intervals and which genic ones
are left out of --out-promoter-bed:
    gene_free            no flanking gene (== domain_class gene_free)
    single_gene_5prime   genic flank AND the siena is the 5'/upstream (promoter)
    multi_gene_5prime    interval of >=1 flanking gene: right gene is '+' OR left
                         gene is '-'. These are EXACTLY the sienas that enter
                         --out-promoter-bed.
    single_gene_3prime   genic flank(s) but the siena is the promoter of NO gene
    multi_gene_3prime    (downstream-only, convergent, or strand-unresolved).
                         These are the qualifying genic sienas LEFT OUT of
                         --out-promoter-bed.

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
--out-csv      Full siena table (CSV), incl. domain_class, siena_class,
               left_strand, right_strand.
--out-bed      Optional siena BED6 (0-based, half-open) for browsers/bedtools.
--out-genebed  Optional gene-body BED6: full spans of the genes inside
               single_gene + multi_gene domains that yield a qualifying siena
               (strand-aware; column 6 carries the gene strand).
--out-promoter-bed  Optional STRAND-AWARE BED6 with one 5'-flanking (promoter)
               siena per locus -- the siena_class *_5prime set -- for deepTools
               `computeMatrix scale-regions`.

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


def build_gene_bodies(gtf_path, feature="exon", strip_suffix=None):
    """Return (genes, strand_by_id, n_units, n_mixed_strand).

    genes          -- {chrom: (starts[], ends[], ids[])} sorted by start.
    strand_by_id   -- {id: '+'|'-'|'.'}; '.' if a locus's records disagree on
                      strand (then it is excluded from strand-aware outputs).
    n_units        -- number of bodies built (loci if strip_suffix, else genes).
    n_mixed_strand -- number of ids dropped to '.' for strand disagreement.

    Gene body = feature-union span per id (min start -> max end, introns
    included). The id is `gene_id` from the GTF, optionally truncated at the
    first occurrence of `strip_suffix` to collapse isoforms onto one locus
    (e.g. strip_suffix='-mRNA' maps 'Nb01g01059-mRNA' -> 'Nb01g01059'). Strand
    is read from GTF column 7.
    """
    span = {}      # id -> [chrom, start, end]
    strands = {}   # id -> set of strands seen
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
                gid = gid.split(strip_suffix)[0]   # locus id
            chrom, start, end, strand = f[0], int(f[3]), int(f[4]), f[6]
            if gid not in span:
                span[gid] = [chrom, start, end]
                strands[gid] = {strand}
            else:
                span[gid][1] = min(span[gid][1], start)
                span[gid][2] = max(span[gid][2], end)
                strands[gid].add(strand)

    if not span:
        sys.exit(f"ERROR: no '{feature}' records with a gene_id found in {gtf_path}")

    strand_by_id = {}
    n_mixed = 0
    for gid, sset in strands.items():
        sset = {s for s in sset if s in ("+", "-")}
        if len(sset) == 1:
            strand_by_id[gid] = sset.pop()
        else:
            strand_by_id[gid] = "."     # none, or conflicting
            if len(sset) > 1:
                n_mixed += 1

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
    return genes, strand_by_id, len(span), n_mixed


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
    ap.add_argument("--strip-suffix", default=None,
                    help="Collapse isoforms to loci: truncate each gene_id at "
                         "the first occurrence of this string to form the locus "
                         "id (e.g. '-mRNA' maps Nb01g01059-mRNA -> Nb01g01059). "
                         "Required for one-promoter-per-locus output.")
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
    ap.add_argument("--out-flanking-genebed", default=None,
                    help="Optional gene-body BED6 of ONLY the genes directly "
                         "flanking a siena that passed all thresholds (i.e. a "
                         "left_gene/right_gene of the final siena set). Stricter "
                         "than --out-genebed, which keeps every gene in a "
                         "qualifying domain. Strand-aware; one row per gene.")
    ap.add_argument("--out-downstream-bed", default=None,
                    help="Optional STRAND-AWARE BED6 with exactly one 3'-flanking "
                         "(downstream) siena per gene -- the mirror of "
                         "--out-promoter-bed. + gene -> the siena to its right; "
                         "- gene -> the siena to its left. Column 6 carries the "
                         "gene strand. A convergent siena (+ gene left, - gene "
                         "right) is the 3' region of both and is emitted once per "
                         "gene. Use with --strip-suffix for per-locus.")
    ap.add_argument("--out-promoter-bed", default=None,
                    help="Optional STRAND-AWARE BED6 with exactly one 5'-flanking "
                         "(promoter) siena per locus, for deepTools "
                         "`computeMatrix scale-regions`. + locus -> its left "
                         "siena; - locus -> its right siena. Column 6 carries the "
                         "locus strand. Use with --strip-suffix for per-locus.")
    args = ap.parse_args()

    genes, strand_by_id, n_genes_total, n_mixed = build_gene_bodies(
        args.gtf, feature=args.feature, strip_suffix=args.strip_suffix)
    unit = "loci" if args.strip_suffix else "genes"
    print(f"[gtf] {unit} reconstructed: {n_genes_total:,}"
          + (f" (isoforms collapsed on '{args.strip_suffix}')" if args.strip_suffix else ""),
          file=sys.stderr)
    if n_mixed:
        print(f"[warn] {n_mixed:,} {unit} had conflicting strands -> strand '.', "
              f"excluded from --out-promoter-bed", file=sys.stderr)


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

    # Strand of each flank (gene/locus strand; '.' at a domain edge or if the
    # GTF strand was missing/conflicting).
    S["left_strand"] = S.left_gene.map(lambda g: strand_by_id.get(g, "."))
    S["right_strand"] = S.right_gene.map(lambda g: strand_by_id.get(g, "."))

    # siena_class: refine domain_class by the siena's promoter side.
    #   *_5prime  = siena is the 5'/upstream (promoter) interval of >=1 flanking
    #               gene  (right gene '+'  OR  left gene '-')  -> enters
    #               --out-promoter-bed.
    #   *_3prime  = genic flank(s) but promoter of none (downstream-only,
    #               convergent, or strand unresolved) -> left OUT of promoter bed.
    is_promoter = (S.right_strand == "+") | (S.left_strand == "-")
    prefix = np.where(S.n_genes_in_domain == 1, "single_gene", "multi_gene")
    S["siena_class"] = np.where(
        S.n_genes_in_domain == 0, "gene_free",
        np.where(is_promoter, prefix + "_5prime", prefix + "_3prime"))

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
            "strand": GB.gene_id.map(lambda g: strand_by_id.get(g, ".")),
        })
        gene_bed.to_csv(args.out_genebed, sep="\t", header=False, index=False)
        print(f"[result] gene bodies in genebed     : {len(gene_bed):,} "
              f"(single+multi gene domains with a qualifying siena)", file=sys.stderr)

    if args.out_flanking_genebed:
        # ONLY the genes directly flanking a siena that survived all thresholds:
        # the left_gene / right_gene of the final siena table S. Stricter than
        # --out-genebed (which keeps every gene in a qualifying domain). Full
        # gene spans are looked up from gene_rows (unclipped), strand-aware.
        body_by_id = {}
        for (c, gs, ge, gid, ds, de, ng) in gene_rows:
            body_by_id[gid] = (c, gs, ge)   # full span; consistent per gene id

        flank_ids = {g for g in pd.concat([S.left_gene, S.right_gene],
                                          ignore_index=True).unique()
                     if g not in ("domain_start", "domain_end")}
        fb_rows = []
        for gid in flank_ids:
            if gid in body_by_id:
                c, gs, ge = body_by_id[gid]
                fb_rows.append((c, gs - 1, ge, gid, 0,
                                strand_by_id.get(gid, ".")))
        fb = (pd.DataFrame(fb_rows,
                           columns=["chrom", "start", "end", "name", "score", "strand"])
                .drop_duplicates()
                .sort_values(["chrom", "start", "end"]))
        fb.to_csv(args.out_flanking_genebed, sep="\t", header=False, index=False)
        print(f"[result] flanking gene bodies       : {len(fb):,} "
              f"(genes directly bordering a surviving siena)", file=sys.stderr)

    if args.out_promoter_bed:
        # One 5'-flanking (promoter) siena per locus, strand-aware.
        #   + locus: promoter = siena whose right_gene == locus (interval on its
        #            left, i.e. upstream of a + gene).
        #   - locus: promoter = siena whose left_gene == locus (interval on its
        #            right, i.e. upstream of a - gene).
        # A single intergenic siena flanked by a - gene on its left and a + gene
        # on its right is the (shared) promoter of BOTH -> emitted once per locus
        # (bidirectional promoter). Loci whose 5' side is genic, or whose
        # promoter siena was removed by --min-len, have no promoter and drop out.
        has_l2fc = "log2FoldChange" in S.columns

        plus = S[S.right_gene.map(lambda g: strand_by_id.get(g, ".")) == "+"].copy()
        plus["locus"] = plus.right_gene
        plus["pstrand"] = "+"

        minus = S[S.left_gene.map(lambda g: strand_by_id.get(g, ".")) == "-"].copy()
        minus["locus"] = minus.left_gene
        minus["pstrand"] = "-"

        P = pd.concat([plus, minus], ignore_index=True)
        if len(P) == 0:
            print("[warn] --out-promoter-bed: no promoter sienas found "
                  "(did you pass --strip-suffix, and does the GTF carry strand?)",
                  file=sys.stderr)
        # Exactly one per locus: if a locus has several candidates (e.g. it
        # borders two induced domains), keep the strongest |log2FC|.
        P["_rank"] = P["log2FoldChange"].abs() if has_l2fc else 0.0
        P = (P.sort_values(["locus", "_rank"], ascending=[True, False])
               .drop_duplicates("locus", keep="first"))

        prom = pd.DataFrame({
            "chrom": P.Chrom,
            "start": P.siena_start - 1,       # 1-based inclusive -> 0-based half-open
            "end": P.siena_end,
            "name": P.locus,
            "score": (P["log2FoldChange"] if has_l2fc else 0),
            "strand": P.pstrand,
        }).sort_values(["chrom", "start", "end"])
        prom.to_csv(args.out_promoter_bed, sep="\t", header=False, index=False)

        # Report how many strand-bearing loci in induced domains got no promoter.
        loci_in_play = {g for g in pd.concat([S.left_gene, S.right_gene]).unique()
                        if strand_by_id.get(g, ".") in ("+", "-")}
        dropped = len(loci_in_play) - len(prom)
        npos = int((prom.strand == "+").sum())
        nneg = int((prom.strand == "-").sum())
        print(f"[result] promoter sienas (1/locus)  : {len(prom):,} "
              f"(+ {npos:,} / - {nneg:,})", file=sys.stderr)
        print(f"[result] loci w/ no 5' siena dropped: {dropped:,} "
              f"(5' side genic, or promoter siena below --min-len)", file=sys.stderr)

    if args.out_downstream_bed:
        # One 3'-flanking (downstream) siena per gene, strand-aware -- the mirror
        # of --out-promoter-bed.
        #   + gene: downstream = siena whose left_gene == gene (interval on its
        #           right, i.e. 3' of a + gene).
        #   - gene: downstream = siena whose right_gene == gene (interval on its
        #           left, i.e. 3' of a - gene).
        # A siena flanked by a + gene on its left and a - gene on its right is the
        # 3' region of BOTH (convergent) -> emitted once per gene. Genes whose 3'
        # side is genic, or whose 3' siena was removed by --min-len, drop out.
        has_l2fc = "log2FoldChange" in S.columns

        dplus = S[S.left_gene.map(lambda g: strand_by_id.get(g, ".")) == "+"].copy()
        dplus["locus"] = dplus.left_gene
        dplus["pstrand"] = "+"

        dminus = S[S.right_gene.map(lambda g: strand_by_id.get(g, ".")) == "-"].copy()
        dminus["locus"] = dminus.right_gene
        dminus["pstrand"] = "-"

        D = pd.concat([dplus, dminus], ignore_index=True)
        if len(D) == 0:
            print("[warn] --out-downstream-bed: no downstream sienas found "
                  "(does the GTF carry strand?)", file=sys.stderr)
        # Exactly one per gene: if a gene borders two domains, keep strongest |log2FC|.
        D["_rank"] = D["log2FoldChange"].abs() if has_l2fc else 0.0
        D = (D.sort_values(["locus", "_rank"], ascending=[True, False])
               .drop_duplicates("locus", keep="first"))

        down = pd.DataFrame({
            "chrom": D.Chrom,
            "start": D.siena_start - 1,       # 1-based inclusive -> 0-based half-open
            "end": D.siena_end,
            "name": D.locus,
            "score": (D["log2FoldChange"] if has_l2fc else 0),
            "strand": D.pstrand,
        }).sort_values(["chrom", "start", "end"])
        down.to_csv(args.out_downstream_bed, sep="\t", header=False, index=False)

        genic = {g for g in pd.concat([S.left_gene, S.right_gene]).unique()
                 if strand_by_id.get(g, ".") in ("+", "-")}
        dropped_d = len(genic) - len(down)
        npos = int((down.strand == "+").sum())
        nneg = int((down.strand == "-").sum())
        print(f"[result] downstream sienas (1/gene) : {len(down):,} "
              f"(+ {npos:,} / - {nneg:,})", file=sys.stderr)
        print(f"[result] genes w/ no 3' siena       : {dropped_d:,} "
              f"(3' side genic, or 3' siena below --min-len)", file=sys.stderr)

    # Console summary.
    print(f"[result] siena_class breakdown:", file=sys.stderr)
    for k in ["gene_free", "single_gene_5prime", "single_gene_3prime",
              "multi_gene_5prime", "multi_gene_3prime"]:
        n = int((S.siena_class == k).sum())
        if n:
            tail = "  -> --out-promoter-bed" if k.endswith("_5prime") else (
                   "  (left out of promoter bed)" if k.endswith("_3prime") else "")
            print(f"    {k:<20} {n:,}{tail}", file=sys.stderr)
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
