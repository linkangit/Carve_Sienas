# Calling SIENAs from differential H3K9ac domains

**SIENA** = *Stimulus-Induced Enhancer Acetylation*: an intergenic
sub-interval of a ChIP domain that gains a histone mark on stimulation. This
repository contains `carve_sienas.py`, a small, dependency-light tool that
defines sienas directly from a set of differential ChIP domains and a gene
annotation — **no peak-annotation package (e.g. ChIPseeker) required**.

This README walks through the idea, the inputs, the command line, and how to
read the output. The worked examples use a tomato H3K9ac dataset (+JA vs.
control), but nothing in the tool is species- or mark-specific.

---

## 1. The idea

A stimulus-induced enhancer signal often shows up as a broad domain of an
activating mark (here H3K9ac) that grows on induction. But a domain called by a
broad-peak caller frequently spans gene bodies as well as the intergenic space
between genes. The regulatory element we care about is the **intergenic
portion** of that domain — the part that is not inside a transcribed gene.

So a siena is obtained by **interval subtraction**: take an induced domain,
remove every gene-body span it overlaps, and keep the leftover intergenic
pieces. A domain with *N* internal genes yields up to *N*+1 sienas; a domain
lying entirely inside a gene body yields none.

Two conventions make this exact and reproducible:

- **Near-edge boundaries.** A siena is bounded by the nearest edge of each
  flanking gene body, regardless of transcription direction. (The boundary is a
  TSS or a TES depending on the gene's strand, but the tool does not branch on
  strand — it simply uses the gene-body edge.)
- **Gene body = exon-union span.** For each gene, the gene body is taken as the
  minimum exon start to the maximum exon end (introns included). This works even
  with annotations that contain only `exon`/`CDS` records and no explicit
  `gene` lines, such as Liftoff output.

Inclusive 1-bp gaps are used throughout: a siena ends at `gene_start − 1` and
the next one resumes at `gene_end + 1`.

---

## 2. How a domain is carved

Consider a single induced domain and the genes that fall within it.

**Case A — the domain begins in intergenic space.** The first siena runs from
the domain start to the near edge of the first gene; each subsequent siena runs
from where one gene ends to where the next begins; the last runs from the final
gene to the domain end.

```
domain        |=================================================|
genes                 [ gene 1 ]            [ gene 2 ]
sienas        |======|          |==========|          |=========|
              siena 1            siena 2                siena 3
```

**Case B — the domain begins inside a gene body.** No siena starts at the
domain start (it is genic); the first siena begins at that gene's far edge.

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

## 3. Requirements

- Python 3.8+
- `pandas`
- `numpy`

```bash
pip install pandas numpy
```

---

## 4. Inputs

**Differential domain table (`--domains`, CSV).** One row per domain, with
chromosome / start / end columns plus any statistics you want carried onto the
sienas. The defaults expect epic2-style headers; override them if yours differ:

| What | Default column | Flag to change it |
|------|----------------|-------------------|
| Chromosome | `#Chromosome` | `--domain-chrom` |
| Start | `Start` | `--domain-start` |
| End | `End` | `--domain-end` |
| Stats to inherit | `log2FoldChange FDR PValue Score ChIPCount InputCount` | `--carry` |

Coordinates are treated as 1-based inclusive (as epic2 emits).

**Gene annotation (`--gtf`, GTF).** Gene bodies are reconstructed from `exon`
records by default (change with `--feature`). The `gene_id` attribute in column
9 is required; `gene`/`transcript`/UTR lines are not.

> The domain table is assumed to already contain the domains you consider
> induced (e.g. positive fold-change, significant FDR). The thresholds in
> Section 6 let you tighten that further, but the tool does not impose any
> direction by default.

---

## 5. Basic usage

```bash
python3 carve_sienas.py \
  --gtf      annotation.gtf \
  --domains  diff_domains.csv \
  --out-csv  sienas.csv
```

With non-default column names and a BED export for genome browsers / bedtools:

```bash
python3 carve_sienas.py \
  --gtf      SLM_r2_0-ITAG4_0.gtf \
  --domains  H3K9ac_tomato_2h_vs_ctrl.csv \
  --domain-chrom '#Chromosome' --domain-start Start --domain-end End \
  --carry    log2FoldChange FDR PValue Score ChIPCount InputCount \
  --out-csv  sienas_from_domains.csv \
  --out-bed  sienas_from_domains.bed
```

The tool prints a run summary to `stderr`, e.g.:

```
[gtf] genes reconstructed: 33,823
[domains] induced domains: 17,797
[result] domains yielding >=1 siena : 17,139
[result] domains yielding 0 sienas  : 658 (wholly inside a gene body)
[result] sienas written             : 28,915
```

---

## 6. Thresholds

Significance and effect-size thresholds apply to the **parent domains** and are
applied **before** carving — each siena inherits its domain's statistics, so a
"differential" cutoff is really a filter on which domains are allowed to
contribute sienas. Every threshold logs how many domains it removed.

| Flag | Keeps domains where | Typical use |
|------|---------------------|-------------|
| `--min-log2fc` | `log2FoldChange >= value` | `1.0` for ≥2-fold gains |
| `--min-abs-log2fc` | `|log2FoldChange| >= value` | two-sided (gains **and** losses) |
| `--max-fdr` | `FDR <= value` | `0.05` standard, `0.01` strict |
| `--max-pvalue` | `PValue <= value` | raw-p alternative to FDR |
| `--min-score` | `Score >= value` | caller-specific signal floor |
| `--min-len` | siena length `>= value` (bp) | drop short slivers, e.g. `250` |

`--min-len` is the one threshold that acts on the sienas themselves rather than
the domains.

Example — a stringent shortlist:

```bash
python3 carve_sienas.py \
  --gtf annotation.gtf --domains diff_domains.csv \
  --min-log2fc 1.0 --max-fdr 0.01 --min-len 250 \
  --out-csv sienas_strict.csv
```

### Choosing values

There is no universal cutoff. A practical approach:

- Use **FDR** as the significance lever — `0.05` for a sensitive catalogue,
  `0.01` for a high-confidence core.
- Treat the **log2FC** floor as a coverage-vs-stringency trade-off. Histone
  marks often shift more modestly than mRNA, so a hard 2-fold cut can discard
  real, modest induction. A looser effect-size floor (~1.4–1.5-fold) paired with
  a strict FDR is a common compromise.
- Because statistics are inherited, you can carve once with no threshold and
  then filter the output CSV in R/pandas to explore cutoffs without rerunning.
  If the siena set is stable across a range of thresholds, your conclusions are
  not threshold-sensitive; if it collapses, the signal is concentrated in a few
  strong domains.

---

## 7. Output columns

`--out-csv` is one row per siena:

| Column | Meaning |
|--------|---------|
| `siena_id` | Stable ID assigned after sorting by genomic position |
| `Chrom` | Chromosome |
| `domain_start`, `domain_end` | Coordinates of the parent domain |
| `siena_start`, `siena_end` | Coordinates of the siena (1-based inclusive) |
| `left_gene`, `right_gene` | Gene ID bounding each side, or `domain_start` / `domain_end` if bounded by the domain edge |
| `n_genes_in_domain` | How many genes the parent domain overlaps |
| `siena_idx_in_domain` | Index of this siena within its parent domain |
| `n_sienas_in_domain` | Total sienas produced by the parent domain |
| `siena_len` | Siena length in bp |
| *(carried)* | Each `--carry` column, inherited from the parent domain |

The `left_gene` / `right_gene` pair gives an immediate candidate-target view:
a siena flanked by two genes lists both, which is the natural starting point for
siena-to-gene linking.

`--out-bed` (optional) is BED6, **0-based half-open** (a siena at 361,800–366,799
inclusive becomes `361799  366799`), with `siena_id` in the name field and
`log2FoldChange` in the score field for quick browser shading.

---

## 8. Notes and caveats

- **Exon-union gene bodies.** A domain extending past a gene's last annotated
  exon into a long 3′-UTR-like tail will treat that tail as siena, because the
  gene body ends at the last exon. Provide a UTR-aware annotation if this
  matters for your locus.
- **Compact genomes.** "Intergenic" means "not inside a gene body," not
  "far from any gene." On a gene-dense genome many sienas sit only a few kb from
  the nearest TSS and may regulate a flanking gene rather than act over long
  range. The `left_gene` / `right_gene` columns let you check this per siena.
- **No strand logic.** Boundaries are gene-body edges; the tool never branches
  on strand. This is intentional (every intergenic gap is kept) and was a
  deliberate design choice for enhancer discovery.
- **Reproducibility.** All behaviour is driven by command-line flags; the run
  summary on `stderr` records the gene count, domain count, and every threshold
  applied, so a logged command fully reproduces a result.

---

## 9. Citation / contact

If you use this in a publication, please cite the repository and the underlying
domain caller (e.g. epic2) and annotation you supplied. Issues and pull requests
welcome.
