#!/usr/bin/env bash
# nf-core/rnaseq 3.21.0, pseudo-alignment only (kallisto), on reads that are already trimmed.
#
# usage: ./run_nfcore.sh <samplesheet.csv> <outdir>
# Needs Java 17+, Nextflow 25.04.x and Docker, e.g. from a conda env:
#   mamba create -n nf25 -c conda-forge -c bioconda openjdk=17 nextflow=25.04.3
#   CONDA_ENV=nf25 ./run_nfcore.sh inputs/samplesheet.csv results/all
set -euo pipefail
cd "$(dirname "$0")"

INPUT=${1:?usage: run_nfcore.sh <samplesheet.csv> <outdir>}
OUTDIR=${2:?usage: run_nfcore.sh <samplesheet.csv> <outdir>}

if [[ -n "${CONDA_ENV:-}" ]]; then
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "$CONDA_ENV"
fi

nextflow run nf-core/rnaseq -r 3.21.0 -profile docker -c conf/local.config \
    --input "$INPUT" \
    --outdir "$OUTDIR" \
    --fasta inputs/GMI1000_chr_megaplasmid.fasta \
    --gff inputs/GMI1000_chr_megaplasmid_cleaned.gff \
    --skip_alignment --skip_trimming --skip_gtf_filter \
    --pseudo_aligner kallisto \
    -work-dir work \
    -resume
