# Ralstonia solanacearum GMI1000 RNA-seq: kallisto + DESeq2

Code for a bacterial RNA-seq analysis of *Ralstonia solanacearum* GMI1000 with a 2 × 2 design
(wild type vs a mutant, at two time points).

This repository contains code, configuration and templates only. No sequencing data or results are included.

## Workflow

1. **Reference.** `scripts/build_reference.sh` downloads the chromosome (`AL646052.1`) and the
   megaplasmid (`AL646053.1`) from NCBI. `scripts/clean_ncbi_gff.py` then rewrites the annotation so
   that every gene has a gene → transcript → exon structure, which RSEM/kallisto require. It also drops
   repeat regions, mobile elements and parentless `misc_RNA` features, and makes duplicated gene IDs unique.
   Including the megaplasmid matters: without it, roughly a third of the genes are missing and the
   pseudo-alignment rate drops considerably.
2. **Quantification.** `run_nfcore.sh` runs [nf-core/rnaseq](https://nf-co.re/rnaseq/3.21.0) 3.21.0
   with `--pseudo_aligner kallisto`, `--skip_alignment`, `--skip_trimming` (the reads are already trimmed)
   and `--skip_gtf_filter`. Strandedness is inferred automatically (`strandedness = auto`).
3. **Differential expression.** `run_deseq2.sh` runs `scripts/deseq2_analysis.R` inside the DESeq2
   container that nf-core/rnaseq uses:
   - rRNA genes are removed; genes need ≥ 10 counts in ≥ 5 samples
   - genotype and time levels are read from `inputs/sample_metadata.csv`; the reference genotype is `WT`,
     and the first time level in sort order is the reference time
   - pairwise contrasts use the cell-means model `~ condition`
   - the interaction test uses `~ genotype + time + genotype:time`
   - genes are called significant at padj < 0.05 and |log2FC| ≥ 1

## Requirements

- Docker (runnable without sudo)
- Java 17+ and Nextflow 25.04.x, for example:
  ```bash
  mamba create -n nf25 -c conda-forge -c bioconda openjdk=17 nextflow=25.04.3
  ```
- Python 3, `curl`, and optionally `pdftoppm` (poppler) to convert the plots to PNG

## Usage

```bash
# 1. reference -> inputs/
./scripts/build_reference.sh

# 2. samplesheet (see templates/samplesheet.csv), then quantification
CONDA_ENV=nf25 ./run_nfcore.sh inputs/samplesheet.csv results/all

# 3. sample design (see templates/sample_metadata.csv) -> inputs/sample_metadata.csv, then DESeq2
./run_deseq2.sh
```

Adjust `conf/local.config` to your machine's CPUs and memory. nf-core/rnaseq 3.21.0 uses `resourceLimits`;
`--max_cpus` and `--max_memory` no longer exist.

## Outputs (`results/deseq2/`)

| File | Content |
|---|---|
| `DE_<contrast>_all.tsv` / `_significant.tsv` | log2FC (MLE and shrunk), padj, group means, gene annotation |
| `DE_summary.tsv` | number of significant genes per contrast |
| `normalized_counts.tsv`, `vst_counts.tsv` | normalised and variance-stabilised counts |
| `plots/` | volcano plots, MA plots, PCA, sample-distance heatmap |

## Interpretation notes

- Identical multi-copy genes, such as IS-element transposases, receive equal shares of ambiguous reads
  from kallisto. They therefore show identical statistics and should be read as one signal.
- Genes directly downstream of a deletion can change because of the deletion construct (a polar effect)
  rather than because of the deleted gene's function.

## Software

nf-core/rnaseq 3.21.0 · Nextflow 25.04.3 · kallisto 0.51.1 · DESeq2 1.28.0
