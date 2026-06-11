# 🧬 SIENA caller

**S**timulus-**I**nduced **EN**hancer-locus **A**nnotation — carve the intergenic
parts of differential ChIP domains into candidate regulatory loci, straight from
a domain table and a gene annotation.

> No peak-annotation package (ChIPseeker, etc.) required. Pure `pandas` + `numpy`.
> Works on **any genome, organism, mark, or sample** — anything that yields a
> table of differential domains plus a GTF.

A SIENA is the **intergenic portion of an induced ChIP domain** — what's left
after you subtract the gene bodies. Each siena inherits its parent domain's
statistics and is labelled by how it sits relative to the surrounding genes, so
you can immediately separate 5′/promoter intervals from 3′/downstream ones, and
emit strand-aware tracks ready for deepTools.

```
domain      |===========================================|
genes              [ gene 1 ]          [ gene 2 ]
sienas      |=====|           |========|          |=====|
            siena 1            siena 2             siena 3
```

---

## Contents

- [Quick start](#quick-start)
- [Install](#install)
- [How carving works](#how-carving-works)
- [The two labels: `domain_class` and `siena_class`](#the-two-labels)
- [Inputs](#inputs)
- [All flags](#all-flags)
- [Outputs](#outputs)
- [Which track for deepTools?](#which-track-for-deeptools)
- [Adapting to your GTF](#adapting-to-your-gtf)
- [Troubleshooting](#troubleshooting)
- [Notes & caveats](#notes--caveats)
- [Citation](#citation)

---

## Quick start

```bash
pip install pandas numpy

python3 carve_sienas.py \
  --gtf      annotation.gtf \
  --domains  diff_domains.csv \
  --min-log2fc 1.0 \
  --min-len    1000 \
  --out-csv            sienas_classified.csv \
  --out-bed            sienas.bed \
  --out-genebed        gene_bodies.bed \
  --out-flanking-genebed flanking_gene_bodies.bed \
  --out-promoter-bed   promoter_sienas.bed \
  --out-downstream-bed downstream_sienas.bed
```

Every output except `--out-csv` is optional — ask only for the tracks you need.
Two flags adapt the tool to how your GTF is built (see
[Adapting to your GTF](#adapting-to-your-gtf)):

- **`--feature transcript`** if `gene_id` is on `transcript` lines, not `exon`
  lines. Default is `exon`.
- **`--strip-suffix=<str>`** to collapse alternative transcripts onto one locus
  (e.g. `--strip-suffix=.t`). ⚠️ Use the `=` form — a leading-dash value is
  misread by argparse.

The domain table's delimiter (comma **or** tab) is auto-detected, so a raw caller
table (e.g. an epic2 TSV) works untouched.

---

## Install

```bash
pip install -r requirements.txt   # pandas, numpy   (Python 3.8+)
```

---

## How carving works

A stimulus-induced signal often appears as a broad domain of an activating mark
that grows on induction, spanning gene bodies *and* the intergenic space between
them. SIENA isolates the intergenic space by **interval subtraction**: take a
domain, remove every gene-body span it overlaps, keep the leftover non-genic
pieces.

**Conventions**

- A domain with **N** internal genes yields **N−1 to N+1** sienas. A domain
  wholly inside a gene body yields none.
- Boundaries use the **nearest edge** of each flanking gene body (strand-agnostic)
  with inclusive 1-bp gaps: a siena ends at `gene_start − 1`; the next resumes at
  `gene_end + 1`.
- A **gene body** is the feature-union span per gene/locus (min feature start →
  max feature end, introns included).

### Case A — domain begins and ends in intergenic space

```
domain      |===========================================|
genes              [ gene 1 ]          [ gene 2 ]
sienas      |=====|           |========|          |=====|
            siena 1            siena 2             siena 3
```

### Case B — a domain edge falls *inside* a gene body

If the domain **starts inside** a gene (or **ends inside** one), the genic stretch
is **not** emitted as a siena — carving resumes at that gene's far edge. No genic
sequence leaks into a siena:

```
domain            |=====================================|
genes      [ ---- gene 1 ---- ]        [ gene 2 ]
sienas                          |======|          |=====|
                                siena 1            siena 2
```

Worked example — domain `11000–25000` opens inside gene 1 (`9000–13000`) and ends
intergenic at `25000`, gene 2 at `18000–20000`:

| siena | coordinates | bounded by | note |
|-------|-------------|------------|------|
| — | (11000–13000) | — | inside gene 1 → **no siena** |
| 1 | `13001–17999` | gene 1 → gene 2 | starts at gene 1's far edge |
| 2 | `20001–25000` | gene 2 → domain end | runs to the intergenic domain end |

The rule is symmetric for a domain ending inside a gene.

---

## The two labels

Every siena carries two labels so you can slice the output without re-deriving
anything.

### `domain_class` — what the parent **domain** spans

| value | genes in parent domain |
|-------|------------------------|
| `gene_free` | 0 — fully intergenic domain |
| `single_gene` | 1 |
| `multi_gene` | ≥ 2 |

### `siena_class` — the siena's **promoter side** (strand-aware)

| value | meaning |
|-------|---------|
| `gene_free` | no flanking gene |
| `single_gene_5prime` · `multi_gene_5prime` | siena is the **5′ / promoter** interval of ≥ 1 flanking gene — right gene is `+` **or** left gene is `−` |
| `single_gene_3prime` · `multi_gene_3prime` | genic flank(s) but promoter of **none** — downstream-only, convergent (`+` left / `−` right), or strand-unresolved |

The CSV also reports `left_strand` / `right_strand`, so every call is auditable.

> **Divergent / convergent sienas serve two genes.** A siena with a `−` gene on
> its left and a `+` gene on its right is the shared 5′ promoter of *both*; one
> with `+` on the left and `−` on the right is the 3′ region of both. The per-gene
> tracks (`--out-promoter-bed`, `--out-downstream-bed`) emit such a siena once per
> gene, each with that gene's strand.

---

## Inputs

### Domain table — `--domains`

One row per domain. Defaults expect epic2-style headers (override as needed):

| field | default column | flag |
|-------|----------------|------|
| chromosome | `#Chromosome` | `--domain-chrom` |
| start | `Start` | `--domain-start` |
| end | `End` | `--domain-end` |
| stats to inherit | `log2FoldChange FDR PValue Score ChIPCount InputCount` | `--carry` |

Coordinates are 1-based inclusive. Delimiter (comma/tab) is auto-detected; force
it with `--sep` (`,`, `'\t'`, `comma`, or `tab`). Any differential-domain caller
works once the chrom/start/end columns are mapped.

### Gene annotation — `--gtf`

Gene bodies = feature-union span per `gene_id`, strand from column 7.
**The `gene_id "…"` attribute must be on whatever feature `--feature` selects.**
See [Adapting to your GTF](#adapting-to-your-gtf).

---

## All flags

<details>
<summary><b>Parsing & annotation</b></summary>

| flag | default | purpose |
|------|---------|---------|
| `--gtf` | *(required)* | gene annotation GTF |
| `--domains` | *(required)* | differential domain table (CSV or TSV) |
| `--feature` | `exon` | feature to build bodies from; **must carry `gene_id`** |
| `--strip-suffix` | none | collapse isoforms to loci (`--strip-suffix=.t`) |
| `--sep` | auto | domain delimiter (`,`, `'\t'`, `comma`, `tab`) |
| `--domain-chrom` / `--domain-start` / `--domain-end` | `#Chromosome` / `Start` / `End` | domain column names |
| `--carry` | log2FoldChange FDR PValue Score ChIPCount InputCount | columns inherited onto each siena |

</details>

<details>
<summary><b>Thresholds</b> (applied to parent domains, before carving)</summary>

| flag | keeps | typical |
|------|-------|---------|
| `--min-log2fc` | `log2FoldChange ≥ value` | `1.0` (≥ 2-fold gain) |
| `--min-abs-log2fc` | `|log2FoldChange| ≥ value` | two-sided (gains + losses) |
| `--max-fdr` | `FDR ≤ value` | `0.05` / `0.01` |
| `--max-pvalue` | `PValue ≤ value` | raw-p alternative |
| `--min-score` | `Score ≥ value` | caller signal floor |
| `--min-len` | siena length `≥ value` bp | **acts on sienas only**, never gene bodies |
| `--require-genic-domain` | `single_gene` + `multi_gene` only | drops `gene_free` sienas |

Stats are inherited, so you can carve once with no threshold and filter the CSV
afterwards to explore cutoffs without rerunning.

</details>

<details>
<summary><b>Outputs</b></summary>

| flag | writes |
|------|--------|
| `--out-csv` | *(required)* full per-siena table |
| `--out-bed` | all siena intervals (BED6, **strand-agnostic**) |
| `--out-genebed` | every gene in a qualifying domain (BED6, strand-aware) |
| `--out-flanking-genebed` | only genes bordering a surviving siena (BED6, strand-aware) |
| `--out-promoter-bed` | one 5′ siena per gene (BED6, strand-resolved) |
| `--out-downstream-bed` | one 3′ siena per gene (BED6, strand-resolved) |

</details>

### Built-in guards

The tool fails loudly instead of producing misleading output:

| condition | behaviour |
|-----------|-----------|
| domain table parsed to 1 column | **error** — wrong delimiter |
| `Chrom`/`Start`/`End` unresolved | **error** — check `--domain-*` / `--sep` |
| no chromosome names shared by GTF & domains | **error** — rename one side |
| some domain chroms absent from GTF | warning — those become `gene_free` |
| no `--feature` records with `gene_id` | **error** — wrong feature / non-GTF attrs |
| a locus's records disagree on strand | warning — strand `.`, dropped from strand-resolved tracks |

---

## Outputs

### `--out-csv` — one row per siena

Columns: `siena_id`, `Chrom`, `domain_start/end`, `siena_start/end` (1-based
inclusive), `left_gene` / `right_gene`, `left_strand` / `right_strand`,
`n_genes_in_domain`, `siena_idx_in_domain` / `n_sienas_in_domain`, `siena_len`,
`domain_class`, `siena_class`, + every `--carry` column.

```bash
# class breakdown at a glance
cut -d, -f16 sienas_classified.csv | tail -n +2 | sort | uniq -c
```
*(`siena_class` is column 16 in the default order — confirm with
`head -1 sienas_classified.csv` before scripting against the index.)*

### `--out-bed` — all siena intervals

BED6, 0-based half-open, `siena_id` in name, `log2FoldChange` in score.
**Column 6 is `.`** — a siena is intergenic and has no intrinsic strand (it can
flank genes of opposite orientation), so this track is strand-agnostic by design.

### `--out-genebed` — induced gene bodies *(strand-aware)*

Full gene-body spans inside `single_gene` + `multi_gene` domains that yield ≥ 1
qualifying siena, deduplicated. Column 6 carries the gene strand.

### `--out-flanking-genebed` — genes flanking a siena *(strand-aware)*

A **stricter** version of `--out-genebed`: only genes that are a
`left_gene`/`right_gene` of a siena that **survived all thresholds**. Interior
genes whose flanking sienas were all filtered out are excluded. One row per gene.

| output | example domain `G1 G2 G3`, short gaps around `G2` |
|--------|-----------------------------------------------------|
| `--out-genebed` | `G1, G2, G3` (all genes in the qualifying domain) |
| `--out-flanking-genebed` | `G1, G3` (only genes bordering a surviving siena) |

### `--out-promoter-bed` — one 5′ siena per gene *(strand-resolved)*

For each gene, the siena on its **5′/upstream** side (left siena for `+`, right
siena for `−`), named by the gene, strand = the gene's strand. Divergent
promoters are emitted once per gene. Feed to `computeMatrix scale-regions`.

### `--out-downstream-bed` — one 3′ siena per gene *(strand-resolved)*

The mirror of the promoter track: for each gene, the siena on its **3′/downstream**
side (right siena for `+`, left siena for `−`), named by the gene, strand = the
gene's strand. Convergent sienas are emitted once per gene.

> `--out-promoter-bed` and `--out-downstream-bed` are name-matched per gene, so a
> gene's promoter row and downstream row share a key — handy for paired analyses.

---

## Which track for deepTools?

| goal | anchor on | mode |
|------|-----------|------|
| Gene-focused metagene | `gene_bodies.bed` | `scale-regions` / `reference-point` |
| Gene-focused, sienas-only | `flanking_gene_bodies.bed` | `scale-regions` / `reference-point` |
| Promoter (5′) metagene | `promoter_sienas.bed` | `scale-regions` |
| Downstream (3′) metagene | `downstream_sienas.bed` | `scale-regions` |

All strand-aware BEDs let deepTools flip `−` rows so everything aligns 5′→3′. The
siena BED and gene-body BED have **different row counts by design** — expected for
paired tracks, not an error.

> For a single **continuous** 5′→gene→3′ metagene, give deepTools the gene bodies
> with flanking windows in one `scale-regions` call (`--beforeRegionStartLength` /
> `--afterRegionStartLength`) rather than stitching three separate plots — each
> separate `scale-regions` run rescales its panel independently, so joins won't be
> smooth.

---

## Adapting to your GTF

The only thing that varies between annotations is where `gene_id` lives and how
isoforms are named. Set `--feature` and `--strip-suffix` accordingly:

| your GTF | `--feature` | `--strip-suffix` |
|----------|-------------|------------------|
| `gene_id` on `exon` lines (Ensembl/ITAG-style) | `exon` (default) | only if isoform id splits the locus |
| `gene` / `transcript` lines present, all carry `gene_id` | `gene` (full body to gene end) | usually none |
| `gene_id` only on `transcript` lines | `transcript` | as needed |
| only `exon`/`CDS` records (e.g. Liftoff) | `exon` (or `CDS`) | as needed |
| GFF3 (`ID=…;Parent=…`, no `gene_id`) | — | convert first: `gffread in.gff3 -T -o out.gtf` |

**Isoforms:** if `gene_id` is already the bare locus (e.g. `AT1G01010`, with the
`.1/.2` in `transcript_id`), you need **no** `--strip-suffix` — isoforms collapse
automatically. If the locus id itself carries an isoform suffix (e.g.
`GENE.1`, `GENE.2`), use `--strip-suffix=.`. Inspect with:

```bash
grep -m1 -P '\texon\t' annotation.gtf      # see which attrs the exon lines carry
```

**Chromosome names must match** between the GTF and the domain table (`Chr1` vs
`chr1` vs `1` vs RefSeq `NC_003070.9`). Quick check:

```bash
comm -12 <(grep -v '^#' annotation.gtf | cut -f1 | sort -u) \
         <(tail -n +2 diff_domains.csv | cut -f1 | sort -u)
```

---

## Troubleshooting

<details>
<summary><b>"All my sienas are <code>gene_free</code>" (but IGV shows genes there)</b></summary>

1. **Wrong delimiter.** Many callers emit tab-separated tables even when named
   `.csv`. Check `head -2 file | cat -A` (tabs show as `^I`); auto-detect usually
   handles it, else pass `--sep tab`.
2. **Chromosome-name mismatch** (`Chr1` vs `chr1` vs `1`). No shared names →
   error; partial → warning.
3. **`gene_id` not on the chosen `--feature`** *(most common)*. If `exon` rows
   carry only `transcript_id`, use `--feature transcript`.
4. **GFF3 mislabelled `.gtf`** (`ID=…;Parent=…`, no `gene_id`). Convert with
   `gffread in.gff3 -T -o out.gtf`.

After fixing, delete the stale `*_classified.csv` and rerun. The stderr summary
should show a real `genes reconstructed` count and a real `siena_class breakdown`.

</details>

---

## Notes & caveats

- **Domain-level fold-change is a proxy for induction.** A siena/gene is "induced"
  because its domain's aggregate log2FoldChange passed threshold — not because the
  mark was measured over that interval. To strengthen the claim, recount reads
  over the siena coordinates directly and test treatment vs control there.
- **Siena count ≠ gene-body count, by construction** — N−1…N+1 sienas vs N bodies
  per domain; `--min-len` filters sienas only; gene BEDs are deduplicated.
- **No strand logic in carving.** Boundaries are gene-body edges; strand enters
  only in the labels and the strand-aware BEDs.
- **Reproducibility.** All behaviour is flag-driven; the stderr summary logs gene
  count, domain count, every threshold, and the `siena_class` breakdown — a logged
  command fully reproduces a result.

---

## Citation

If you use this in a publication, please cite this repository together with the
differential-domain caller (e.g. epic2) and the gene annotation you supplied.
