# SIENA caller — intergenic enhancer loci from differential ChIP domains

**SIENA** = *Stimulus-Induced ENhancer-locus Annotation*: an intergenic
sub-interval of a ChIP domain that gains a histone mark on stimulation.
`carve_sienas.py` defines sienas directly from a set of differential ChIP
domains plus a gene annotation — **no peak-annotation package (e.g. ChIPseeker)
required**.

The worked examples use a tomato H3K9ac dataset (+JA vs. control), but nothing
in the tool is species- or mark-specific.

---

## Contents

| File | Purpose |
|------|---------|
| `carve_sienas.py` | The tool |
| `run_example.sh` | One-command reproduction of the worked set |
| `requirements.txt` | Python dependencies |
| `LICENSE` | MIT (add your name) |

---

## 1. The idea

A stimulus-induced enhancer signal often appears as a broad domain of an
activating mark (here H3K9ac) that grows on induction. But a domain called by a
broad-peak caller usually spans gene bodies as well as the intergenic space
between them. The regulatory element of interest is the **intergenic portion**
of that domain — the part not inside a transcribed gene.

A siena is obtained by **interval subtraction**: take an induced domain, remove
every gene-body span it overlaps, and keep the leftover intergenic pieces. A
domain with *N* internal genes yields up to *N*+1 sienas; a domain lying
entirely inside a gene body yields none.

Two conventions make this exact and reproducible:

- **Near-edge boundaries.** A siena is bounded by the nearest edge of each
  flanking gene body, regardless of transcription direction (strand-agnostic).
- **Gene body = exon-union span.** Per gene, the body runs from the minimum exon
  start to the maximum exon end (introns included). This works with annotations
  containing only `exon`/`CDS` records and no explicit `gene` lines (e.g.
  Liftoff output).

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
`gene_free` sienas are not. The `--require-genic-domain` flag (Section 6) keeps
only the single/multi classes when that adjacency is part of your definition.

---

## 4. Install

```bash
pip install -r requirements.txt   # pandas, numpy
```

Python 3.8+.

---

## 5. Inputs

**Differential domain table (`--domains`, CSV).** One row per domain, with
chromosome / start / end plus any statistics to carry onto the sienas. Defaults
expect epic2-style headers; override if yours differ:

| What | Default column | Flag |
|------|----------------|------|
| Chromosome | `#Chromosome` | `--domain-chrom` |
| Start | `Start` | `--domain-start` |
| End | `End` | `--domain-end` |
| Stats to inherit | `log2FoldChange FDR PValue Score ChIPCount InputCount` | `--carry` |

Coordinates are treated as 1-based inclusive (as epic2 emits).

**Gene annotation (`--gtf`, GTF).** Gene bodies are reconstructed from `exon`
records by default (`--feature` to change). The `gene_id` attribute is required;
`gene`/`transcript`/UTR lines are not.

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
| `--min-len` | siena length `>= value` bp | drop short slivers (acts on sienas) |
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

- **Domain-level fold-change is a proxy for gene-body induction.** A gene body is
  treated as "induced" because it lies inside a domain whose aggregate
  log2FoldChange passed threshold — not because H3K9ac was measured over that
  gene body specifically (per-gene counts are not used).
- **Exon-union gene bodies.** A domain extending past a gene's last annotated
  exon into a long 3′-UTR-like tail treats that tail as siena. Provide a
  UTR-aware annotation if this matters.
- **Compact genomes.** "Intergenic" means "not inside a gene body," not "far
  from any gene." The `left_gene` / `right_gene` columns let you check proximity
  per siena.
- **No strand logic.** Boundaries are gene-body edges; the tool never branches on
  strand. Intentional, so every intergenic gap is kept.
- **Reproducibility.** All behaviour is driven by flags; the stderr run summary
  records gene count, domain count, and every threshold applied, so a logged
  command fully reproduces a result.

---

## 10. Citation

If you use this in a publication, please cite this repository together with the
domain caller (e.g. epic2) and gene annotation you supplied. Issues and pull
requests welcome.
