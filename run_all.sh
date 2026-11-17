#!/usr/bin/env bash
# Download inputs (once), then run the four steps in order.
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p data results figures
XENA_GDC=https://gdc-hub.s3.us-east-1.amazonaws.com/download
XENA_PANCAN=https://tcga-pancan-atlas-hub.s3.us-east-1.amazonaws.com/download
MSIGDB=https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.2.Hs
fetch() { [ -s "data/$(basename "$1")" ] || curl -sSf -o "data/$(basename "$1")" "$1"; }
fetch $XENA_GDC/TCGA-UCEC.star_counts.tsv.gz
fetch $XENA_GDC/gencode.v36.annotation.gtf.gene.probemap
fetch $XENA_PANCAN/TCGASubtype.20170308.tsv.gz
fetch $XENA_PANCAN/Survival_SupplementalTable_S1_20171025_xena_sp
for s in h.all c2.cp.reactome c2.cp.kegg_medicus; do fetch $MSIGDB/$s.v2023.2.Hs.symbols.gmt; done
for step in scripts/0*.py; do echo "== $step"; python "$step"; done
