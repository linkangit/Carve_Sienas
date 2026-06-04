# SIENA caller — intergenic enhancer loci from differential ChIP domains

**SIENA** = *Stimulus-Induced ENhancer-locus Annotation*: an intergenic
sub-interval of a ChIP domain that gains a histone mark on stimulation.
`carve_sienas.py` defines sienas directly from a set of differential ChIP
domains plus a gene annotation — **no peak-annotation package (e.g. ChIPseeker)
required**.

The worked examples use a tomato H3K9ac dataset (+JA vs. control), and the tool
has also been run on *Nicotiana benthamiana* (see §5.3 and §9 for the
annotation-format note that case surfaced). Nothing in the tool is species- or
mark-specific.

---

## Contents

| File | Purpose |
|------|---------|
| `carve_sienas.py` | The tool |
| `diagnose_gtf.py` | Diagnostic: reports what the GTF parser actually sees (use when results look wrong) |
| `run_example.sh` | One-command reproduction of the worked set |
| `requirements.txt` | Python dependencies |

---

## 1. The idea

A stimulus-induced enhancer signal often appears as a broad domain of an
activating mark (here H3K9ac) that grows on induction. But a domain called by a
broad-peak caller usually spans gene bodies as well as the intergenic space
between them. The regulatory element of interest is the **intergenic portion**
of that domain — the part not inside a transcribed gene.

A siena is obtained by **interval subtraction**: take an induced domain, remove
every gene-body span it overlaps, and keep the leftover intergenic pieces. A
domain with *N* internal genes yields **N−1 to N+1 sienas** (see §2), or none if
it lies entirely inside a gene body.

Two conventions make this exact and reproducible:

- **Near-edge boundaries.** A siena is bounded by the nearest edge of each
  flanking gene body, regardless of transcription direction (strand-agnostic).
- **Gene body = feature-union span.** Per gene, the body runs from the minimum
  feature start to the maximum feature end (introns included). The feature used
  is set by `--feature` (default `exon`) and **must carry the `gene_id`
  attribute** — see §5.2, which is the single most common source of empty
  results.

Inclusive 1-bp gaps are used throughout: a siena ends at `gene_start − 1` and
the next resumes at `gene_end + 1`.

---

## 2. How a domain is carved

**Case A — domain begins in intergenic space.** The first siena runs from the
domain start to the near edge of the first gene; each subsequent siena spans the
gap between consecutive genes; the last runs from the final gene to the domain
end.

```
domain        |=================================================|
genes                 [ gene 1 ]            [ gene 2 ]
sienas        |======|          |==========|          |=========|
              siena 1            siena 2                siena 3
```

**Case B — domain begins inside a gene body.** No siena starts at the domain
start (it is genic); the first siena begins at that gene's far edge.

```
domain              |===========================================|
genes        [ ---- gene 1 ---- ]         [ gene 2 ]
sienas                           |========|          |==========|
                                  siena 1             siena 2
```

**Siena count vs. gene count.** An N-gene domain produces N−1, N, or N+1 sienas
depending on its edges: each domain edge that falls in intergenic space adds a
flanking siena, each edge inside a gene body does not. Both edges genic → N−1
sienas; both intergenic → N+1; mixed → N. This is why a siena track and a
gene-body track for the same domains have different row counts by construction
(see §9).

Real rows from the tomato run:

| Case | Domain | Gene(s) in domain | Sienas produced |
|------|--------|-------------------|-----------------|
| A | `ch01:370,200–374,599` | Solyc01g005000.3 (+, 371,778–374,003) | `370,200–371,777` and `374,004–374,599` |
| B | `ch01:401,000–404,999` | Solyc01g005020.3 (−, 389,305–402,708); Solyc01g005030.4 (+, 404,753–412,977) | `402,709–404,752` |

In Case B the domain opens inside the minus-strand gene Solyc01g005020, so the
siena begins at that gene's near (right) edge — the near-edge rule in action.

---

## 3. Domain classes

Every siena is labelled by how many genes its parent domain spans:

| `domain_class` | Genes in parent domain | Meaning |
|----------------|------------------------|---------|
| `gene_free` | 0 | Fully intergenic domain; the whole domain is one siena. Not adjacent to any induced gene body. |
| `single_gene` | 1 | Domain spans one gene. |
| `multi_gene` | ≥2 | Domain spans two or more genes. |

`single_gene` and `multi_gene` sienas are, by construction, adjacent to a gene
body lying **inside an induced domain** — i.e. "near an induced gene body."
`gene_free` sienas are not. The `--require-genic-domain` flag (§6) keeps only
the single/multi classes when that adjacency is part of your definition.

> **If you build from `transcript` and `gene_id == transcript_id`** (some
> annotations, including the Nb GTF in §5.3), each isoform counts as its own
> "gene," so a single multi-isoform gene can be labelled `multi_gene`. Carving
> stays correct (overlapping isoforms merge into one genic block), but the
> *counts* and `left_gene`/`right_gene` labels are per-transcript. Collapse
> isoforms first if true gene-level classes matter.

---

## 4. Install

```bash
pip install -r requirements.txt   # pandas, numpy
```

Python 3.8+.

---

## 5. Inputs

### 5.1 Differential domain table (`--domains`)

One row per domain, with chromosome / start / end plus any statistics to carry
onto the sienas. Defaults expect epic2-style headers; override if yours differ:

| What | Default column | Flag |
|------|----------------|------|
| Chromosome | `#Chromosome` | `--domain-chrom` |
| Start | `Start` | `--domain-start` |
| End | `End` | `--domain-end` |
| Stats to inherit | `log2FoldChange FDR PValue Score ChIPCount InputCount` | `--carry` |

Coordinates are treated as 1-based inclusive (as epic2 emits).

**Delimiter is auto-detected.** The table may be comma- *or* tab-separated; raw
epic2 output (which is tab-separated even when named `.csv`) works untouched.
Override with `--sep` if needed — it accepts `,`, the escaped tab `'\t'`, or the
aliases `comma` / `tab`. (A bare literal tab, `--sep $'\t'`, is easily dropped by
the shell; prefer `--sep '\t'` or `--sep tab`.)

If the table parses into a single column, the tool stops with a clear error
naming the delimiter it guessed — that means the separator was wrong, not that
the data is bad.

### 5.2 Gene annotation (`--gtf`)

Gene bodies are reconstructed as the **feature-union span per `gene_id`**
(minimum feature start → maximum feature end, introns included). Works with
annotations that have no explicit `gene` lines (e.g. Liftoff output).

**The `gene_id "..."` attribute must be present on the feature you build from**
(set by `--feature`, default `exon`). This is the critical requirement:

- Tomato/ITAG-style GTF — `gene_id` is on the `exon` lines → use the default
  (`--feature exon`).
- Some GTFs put `gene_id` only on `transcript` lines while `exon` lines carry
  just `transcript_id`. With the default `--feature exon` the parser finds no
  `gene_id`, builds (almost) no gene bodies, and **every domain comes back
  `gene_free`** even though genes clearly sit at those loci in IGV. Build from
  the feature that carries `gene_id` instead — usually `--feature transcript`.

Run `diagnose_gtf.py` (§10) if unsure which feature carries `gene_id`.

### 5.3 *N. benthamiana* example (the `--feature transcript` case)

The Nb GTF used here has `gene_id` on `transcript` rows only; its `exon` rows
carry `transcript_id` alone. The fix is one flag — `--feature transcript`:

```bash
python3 carve_sienas.py \
  --gtf          Benthi.gtf \
  --domains      H3K9ac_Benthi_2h_vs_ctrl.csv \
  --feature      transcript \
  --domain-chrom '#Chromosome' --domain-start Start --domain-end End \
  --carry        log2FoldChange FDR PValue Score ChIPCount InputCount \
  --min-log2fc   0.5 \
  --min-len      500 \
  --out-csv      sienas_classified.csv \
  --out-bed      sienas.bed \
  --out-genebed  gene_bodies.bed
```

`transcript` is the right choice here because it carries `gene_id` and its
coordinates already span the full transcript (UTRs included).

### 5.4 Input validation (built-in guards)

Before carving, the tool fails loudly rather than producing misleading output:

| Guard | Condition | Behaviour |
|-------|-----------|-----------|
| Single-column | domain table parsed to 1 column | **Error** — delimiter wrong; pass `--sep` |
| Missing columns | `Chrom`/`Start`/`End` unresolved after rename | **Error** — check `--domain-chrom/-start/-end` and `--sep` |
| No shared chromosomes | GTF and domains share no chrom names | **Error** — every domain would be `gene_free`; rename one side |
| Partial chromosome match | some domain chroms absent from the GTF | **Warning** — those domains become `gene_free` |
| No genes | no `--feature` records with `gene_id` in the GTF | **Error** — wrong `--feature` or non-GTF attribute format |

---

## 6. Usage

Minimal:

```bash
python3 carve_sienas.py --gtf annotation.gtf --domains diff_domains.csv \
  --out-csv sienas.csv
```

Full worked set (also in `run_example.sh`):

```bash
python3 carve_sienas.py \
  --gtf          SLM_r2_0-ITAG4_0.gtf \
  --domains      H3K9ac_tomato_2h_vs_ctrl.csv \
  --domain-chrom '#Chromosome' --domain-start Start --domain-end End \
  --carry        log2FoldChange FDR PValue Score ChIPCount InputCount \
  --min-log2fc   1.0 \
  --min-len      1000 \
  --require-genic-domain \
  --out-csv      sienas_classified.csv \
  --out-bed      sienas.bed \
  --out-genebed  gene_bodies.bed
```

### Parsing flags

| Flag | Default | Purpose |
|------|---------|---------|
| `--feature` | `exon` | GTF feature to build gene bodies from; **must carry `gene_id`** (use `transcript` if `gene_id` is not on exons) |
| `--sep` | auto-detect | Domain-table delimiter; `,`, `'\t'`, `comma`, or `tab` |
| `--domain-chrom` / `--domain-start` / `--domain-end` | `#Chromosome` / `Start` / `End` | Column names in the domain table |
| `--carry` | log2FoldChange FDR PValue Score ChIPCount InputCount | Domain columns inherited onto each siena |

### Thresholds

Significance / effect-size thresholds apply to the **parent domains** and act
**before** carving — sienas inherit their domain's statistics, so a
"differential" cutoff filters which domains may contribute sienas. Each
threshold logs how many domains it removed.

| Flag | Keeps | Typical use |
|------|-------|-------------|
| `--min-log2fc` | `log2FoldChange >= value` | `1.0` for ≥2-fold gains |
| `--min-abs-log2fc` | `|log2FoldChange| >= value` | two-sided (gains **and** losses) |
| `--max-fdr` | `FDR <= value` | `0.05` standard, `0.01` strict |
| `--max-pvalue` | `PValue <= value` | raw-p alternative |
| `--min-score` | `Score >= value` | caller signal floor |
| `--min-len` | siena length `>= value` bp | drop short slivers (**acts on sienas only**, never on gene bodies) |
| `--require-genic-domain` | single_gene + multi_gene only | enforce "near an induced gene body" |

**Choosing values.** No universal cutoff exists. Use FDR as the significance
lever (`0.05` sensitive, `0.01` strict). Treat the log2FC floor as a
coverage-vs-stringency trade-off: histone marks shift more modestly than mRNA,
so a hard 2-fold cut discards real, modest induction — many people pair a looser
effect-size floor (~1.4–1.5×) with a strict FDR. Because statistics are
inherited, you can carve once with no threshold and filter the output CSV later
to explore cutoffs without rerunning.

---

## 7. Outputs

### `--out-csv` — one row per siena

| Column | Meaning |
|--------|---------|
| `siena_id` | Stable ID assigned after sorting by position |
| `Chrom` | Chromosome |
| `domain_start`, `domain_end` | Parent-domain coordinates |
| `siena_start`, `siena_end` | Siena coordinates (1-based inclusive) |
| `left_gene`, `right_gene` | Bounding gene ID per side, or `domain_start` / `domain_end` at a domain edge |
| `n_genes_in_domain` | Genes the parent domain overlaps |
| `siena_idx_in_domain`, `n_sienas_in_domain` | Index of this siena and total from its domain |
| `siena_len` | Length in bp |
| `domain_class` | `gene_free` / `single_gene` / `multi_gene` |
| *(carried)* | Each `--carry` column, inherited from the parent domain |

The `left_gene` / `right_gene` pair gives an immediate candidate-target view for
siena-to-gene linking.

### `--out-bed` — siena intervals

BED6, **0-based half-open** (a siena at 361,800–366,799 inclusive becomes
`361799  366799`), `siena_id` in the name field, `log2FoldChange` in the score
field for browser shading.

### `--out-genebed` — induced gene bodies

BED6 of the **full gene-body spans** belonging to `single_gene` + `multi_gene`
domains that yield at least one qualifying siena — i.e. the induced gene bodies
that the sienas sit next to. Deduplicated, so a gene shared by adjacent domains
appears once. Pairs with `--out-bed` as a genic track alongside the intergenic
sienas.

> The siena BED and the gene-body BED have **different row counts by design**
> (see §2 and §9). This is expected for paired metagene tracks, not an error.

---

## 8. Worked result (tomato H3K9ac, +JA vs ctrl)

With `--min-log2fc 1.0 --min-len 1000 --require-genic-domain`:

- 17,797 induced domains → 1,899 pass log2FC ≥ 1.0
- Sienas before the genic-domain filter: 1,922 (`gene_free` 1,094, `single_gene`
  431, `multi_gene` 397)
- **Final sienas: 828** (single + multi gene)
- **Gene bodies (`gene_bodies.bed`): 839** induced gene bodies adjacent to those
  sienas

---

## 9. Notes and caveats

- **Siena count ≠ gene-body count, by construction.** For an N-gene domain the
  siena count is N−1 to N+1 (domain-edge geometry, §2), while the gene-body count
  is N. On top of that: `--min-len` drops short *sienas* but never gene bodies;
  `--min-log2fc` filters whole domains so it scales both together;
  `gene_bodies.bed` is deduplicated; and if you build from `transcript` with
  `gene_id == transcript_id`, isoforms inflate the body count without each adding
  a siena. Net: the two tracks legitimately differ in size — a "gap" in paired
  metagene plots needs no reconciliation.
- **Domain-level fold-change is a proxy for gene-body induction.** A gene body is
  treated as "induced" because it lies inside a domain whose aggregate
  log2FoldChange passed threshold — not because H3K9ac was measured over that
  gene body specifically (per-gene counts are not used).
- **`gene_id` must be on the `--feature` you select.** If results are all
  `gene_free` but IGV shows genes there, the feature/`gene_id` pairing is almost
  always the cause (§5.2, §11). Run `diagnose_gtf.py`.
- **Feature-union gene bodies.** A domain extending past a gene's last annotated
  exon into a long 3′-UTR-like tail treats that tail as siena. Provide a
  UTR-aware annotation (or build from `transcript`, which includes UTRs) if this
  matters.
- **Compact genomes.** "Intergenic" means "not inside a gene body," not "far
  from any gene." The `left_gene` / `right_gene` columns let you check proximity
  per siena.
- **No strand logic.** Boundaries are gene-body edges; the tool never branches on
  strand. Intentional, so every intergenic gap is kept.
- **Reproducibility.** All behaviour is driven by flags; the stderr run summary
  records gene count, domain count, and every threshold applied, so a logged
  command fully reproduces a result.

---

## 10. Diagnostics (`diagnose_gtf.py`)

When results look wrong, this mirrors the tool's own parser and prints exactly
what it sees:

```bash
python3 diagnose_gtf.py annotation.gtf diff_domains.csv
```

It reports the GTF's column-3 feature types and counts, the column-1 chromosome
names, sample column-9 attribute strings (so you can see whether `gene_id "..."`
is present and on which feature), how many gene bodies the tool would build, and
whether your first few domains overlap any of them. The number of genes it would
build is the single most decisive value: if it is tiny, your `--feature` /
`gene_id` pairing is wrong.

---

## 11. Troubleshooting: "all my sienas are `gene_free`"

Work down this list; each step has a matching guard or diagnostic.

1. **Wrong delimiter.** Raw epic2 tables are tab-separated even when named
   `.csv`. Auto-detection handles this; if you forced `--sep` incorrectly the
   tool errors with the column it parsed. Check with `head -2 file | cat -A`
   (tabs show as `^I`).
2. **Chromosome-name mismatch.** The domains and GTF must use the same names
   (`Chr01` vs `chr1` vs `Niben…`). No shared names → error; partial overlap →
   warning. Confirm with
   `comm -12 <(cut -f1 ann.gtf | grep -v '^#' | sort -u) <(cut -f1 dom | sort -u)`.
3. **`gene_id` not on the `--feature` you chose (most common).** If `exon` rows
   carry only `transcript_id`, default `--feature exon` finds no `gene_id` and
   builds almost no genes. Switch to `--feature transcript` (or whichever feature
   carries `gene_id`). Confirm with `diagnose_gtf.py` or
   `grep -m1 -P '\texon\t' ann.gtf` to inspect the attribute column.
4. **GFF3 mislabelled `.gtf`.** Attributes like `ID=...;Parent=...` (no
   `gene_id "..."`) are GFF3. Convert to true GTF, e.g.
   `gffread in.gff3 -T -o out.gtf`.

After fixing, delete the stale `*_classified.csv` so you don't confuse it with
the new run, then check the stderr summary: `[gtf] genes reconstructed:` should
be in the tens of thousands and `single_gene` / `multi_gene` classes should
appear.

---

## 12. Citation

If you use this in a publication, please cite this repository together with the
domain caller (e.g. epic2) and gene annotation you supplied. Issues and pull
requests welcome.
