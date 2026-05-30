#!/usr/bin/env bash
# Reproduce the worked tomato H3K9ac (+JA vs ctrl) siena set.
# Edit the paths/thresholds to match your data.
set -euo pipefail

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
