#!/usr/bin/env bash
# Runs scripts/deseq2_analysis.R in the DESeq2 container used by nf-core/rnaseq 3.21.0 (DESeq2 1.28).
# Expects results/all/kallisto/ (from run_nfcore.sh) and inputs/sample_metadata.csv.
set -euo pipefail
cd "$(dirname "$0")"

mkdir -p results/deseq2
python3 scripts/make_gene_annotation.py inputs/GMI1000_chr_megaplasmid_cleaned.gff results/deseq2/gene_annotation.tsv

IMG=quay.io/biocontainers/mulled-v2-8849acf39a43cdd6c839a369a74c0adc823e2f91:ab110436faf952a33575c64dd74615a84011450b-0
docker run --rm --user "$(id -u):$(id -g)" -v "$PWD":"$PWD" -w "$PWD" "$IMG" \
    Rscript scripts/deseq2_analysis.R

# the container's PNG device crashes, so plots are written as PDF; convert if poppler is available
if command -v pdftoppm >/dev/null; then
    for f in results/deseq2/plots/*.pdf; do pdftoppm -png -r 150 -singlefile "$f" "${f%.pdf}"; done
fi
