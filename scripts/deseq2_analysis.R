#!/usr/bin/env Rscript
# DESeq2 analysis for a 2 x 2 design (genotype x time) on nf-core/rnaseq kallisto gene counts.
# inputs/sample_metadata.csv needs columns sample, well, group, genotype, time:
#   genotype: reference level "WT" plus exactly one other level
#   time:     exactly two levels; the first in sort order is the reference (e.g. T1/T2)
# Input: nf-core/rnaseq 3.21.0 kallisto gene counts (results/all), sample metadata, gene annotation.
# Run inside the nf-core DESeq2 container from the project root (see run_deseq2.sh).

suppressPackageStartupMessages({
    library(DESeq2)
    library(ggplot2)
    library(pheatmap)
    library(RColorBrewer)
})

counts_file <- "results/all/kallisto/kallisto.merged.gene_counts.tsv"
meta_file   <- "inputs/sample_metadata.csv"
annot_file  <- "results/deseq2/gene_annotation.tsv"
outdir      <- "results/deseq2"
padj_cut    <- 0.05
lfc_cut     <- 1          # |log2FC| >= 1 for the "significant" lists

dir.create(file.path(outdir, "plots"), showWarnings = FALSE, recursive = TRUE)

## ---- data -------------------------------------------------------------------
cts   <- read.delim(counts_file, check.names = FALSE)
annot <- read.delim(annot_file, check.names = FALSE, quote = "")
meta  <- read.csv(meta_file, stringsAsFactors = FALSE)

rownames(cts) <- cts$gene_id
cts <- as.matrix(cts[, meta$sample])
cts <- round(cts)
storage.mode(cts) <- "integer"

# rRNA removed: depletion efficiency varies between libraries and is not biology
is_rrna <- grepl("rRNA", rownames(cts))
cat("Removing", sum(is_rrna), "rRNA genes\n")
cts <- cts[!is_rrna, ]

geno_ref <- "WT"
geno_alt <- setdiff(unique(meta$genotype), geno_ref)
times    <- sort(unique(meta$time))
stopifnot(geno_ref %in% meta$genotype, length(geno_alt) == 1, length(times) == 2)
t1 <- times[1]; t2 <- times[2]
cnd <- function(g, t) paste(g, t, sep = "_")

meta$genotype  <- factor(meta$genotype, levels = c(geno_ref, geno_alt))
meta$time      <- factor(meta$time, levels = times)
meta$group     <- factor(meta$group)
meta$condition <- factor(cnd(meta$genotype, meta$time),
                         levels = c(cnd(geno_ref, t1), cnd(geno_alt, t1), cnd(geno_ref, t2), cnd(geno_alt, t2)))
rownames(meta) <- meta$sample
stopifnot(identical(colnames(cts), rownames(meta)))

# pre-filter: >= 10 counts in >= 5 samples (smallest group size)
keep <- rowSums(cts >= 10) >= 5
cat("Genes kept after filtering:", sum(keep), "of", nrow(cts), "\n")
cts <- cts[keep, ]

## ---- models -----------------------------------------------------------------
# 1) cell-means model for the pairwise contrasts
dds <- DESeqDataSetFromMatrix(cts, meta, design = ~ condition)
dds <- DESeq(dds, quiet = TRUE)

# 2) full factorial model for the interaction test (same fit, different parameterisation)
dds_int <- DESeqDataSetFromMatrix(cts, meta, design = ~ genotype + time + genotype:time)
dds_int <- DESeq(dds_int, quiet = TRUE)
cat("Interaction model coefficients:", paste(resultsNames(dds_int), collapse = ", "), "\n")

## ---- contrasts --------------------------------------------------------------
contrasts <- list(
    c("condition", cnd(geno_alt, t1), cnd(geno_ref, t1)),   # genotype effect at time 1
    c("condition", cnd(geno_alt, t2), cnd(geno_ref, t2)),   # genotype effect at time 2
    c("condition", cnd(geno_ref, t2), cnd(geno_ref, t1)),   # time effect in reference genotype
    c("condition", cnd(geno_alt, t2), cnd(geno_alt, t1))    # time effect in the other genotype
)
names(contrasts) <- c(paste0(geno_alt, "_vs_", geno_ref, "_", t1),
                      paste0(geno_alt, "_vs_", geno_ref, "_", t2),
                      paste0(t2, "_vs_", t1, "_", geno_ref),
                      paste0(t2, "_vs_", t1, "_", geno_alt))

norm <- counts(dds, normalized = TRUE)
grp_mean <- sapply(levels(meta$condition), function(g)
    rowMeans(norm[, meta$condition == g, drop = FALSE]))
colnames(grp_mean) <- paste0("mean_norm_", colnames(grp_mean))

annotate <- function(res, shrunk) {
    df <- data.frame(gene_id = rownames(res),
                     baseMean = res$baseMean,
                     log2FoldChange = res$log2FoldChange,
                     log2FoldChange_shrunk = shrunk$log2FoldChange,
                     lfcSE = res$lfcSE,
                     stat = res$stat,
                     pvalue = res$pvalue,
                     padj = res$padj,
                     check.names = FALSE)
    df <- merge(annot, df, by = "gene_id", all.y = TRUE, sort = FALSE)
    df <- cbind(df, grp_mean[df$gene_id, , drop = FALSE])
    df[order(df$padj, -abs(df$log2FoldChange), na.last = TRUE), ]
}

volcano <- function(df, title, file) {
    df$sig <- ifelse(!is.na(df$padj) & df$padj < padj_cut & abs(df$log2FoldChange) >= lfc_cut,
                     ifelse(df$log2FoldChange > 0, "up", "down"), "ns")
    df$label <- ifelse(df$gene_name != "", paste0(df$locus_tag, " (", df$gene_name, ")"), df$locus_tag)
    top <- head(df[df$sig != "ns", ], 15)
    df$mlog10 <- -log10(pmax(df$padj, 1e-300))
    top$mlog10 <- -log10(pmax(top$padj, 1e-300))
    n_up <- sum(df$sig == "up"); n_dn <- sum(df$sig == "down")
    p <- ggplot(df, aes(log2FoldChange, mlog10, colour = sig)) +
        geom_point(size = 1, alpha = 0.7) +
        scale_colour_manual(values = c(down = "#2166ac", ns = "grey70", up = "#b2182b")) +
        geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed", colour = "grey40") +
        geom_hline(yintercept = -log10(padj_cut), linetype = "dashed", colour = "grey40") +
        geom_text(data = top, aes(label = label), size = 2.5, vjust = -0.6, colour = "black",
                  check_overlap = TRUE) +
        labs(title = title, subtitle = sprintf("padj < %.2f & |log2FC| >= %g: %d up, %d down",
                                               padj_cut, lfc_cut, n_up, n_dn),
             x = "log2 fold change", y = "-log10 adjusted p-value", colour = NULL) +
        theme_bw()
    ggsave(file, p, width = 7, height = 6)
}

summary_rows <- list()
for (nm in names(contrasts)) {
    res <- results(dds, contrast = contrasts[[nm]], alpha = padj_cut)
    shr <- lfcShrink(dds, contrast = contrasts[[nm]], res = res, type = "normal", quiet = TRUE)
    df  <- annotate(res, shr)
    write.table(df, file.path(outdir, paste0("DE_", nm, "_all.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
    sig <- df[!is.na(df$padj) & df$padj < padj_cut & abs(df$log2FoldChange) >= lfc_cut, ]
    write.table(sig, file.path(outdir, paste0("DE_", nm, "_significant.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
    volcano(df, nm, file.path(outdir, "plots", paste0("volcano_", nm, ".pdf")))
    pdf(file.path(outdir, "plots", paste0("MA_", nm, ".pdf")), width = 7, height = 5.5)
    plotMA(shr, main = nm, alpha = padj_cut, ylim = c(-6, 6)); dev.off()
    summary_rows[[nm]] <- data.frame(contrast = nm,
        padj05_any_lfc = sum(df$padj < padj_cut, na.rm = TRUE),
        up = sum(sig$log2FoldChange > 0), down = sum(sig$log2FoldChange < 0))
}

# interaction: does the genotype effect differ between the two time points?
int_name <- tail(resultsNames(dds_int), 1)
res <- results(dds_int, name = int_name, alpha = padj_cut)
# type="normal" shrinkage is not available for interaction terms in DESeq2 1.28 (apeglm/ashr not installed)
shr <- data.frame(log2FoldChange = rep(NA_real_, nrow(res)))
df  <- annotate(res, shr)
write.table(df, file.path(outdir, "DE_interaction_genotype_x_time_all.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
sig <- df[!is.na(df$padj) & df$padj < padj_cut & abs(df$log2FoldChange) >= lfc_cut, ]
write.table(sig, file.path(outdir, "DE_interaction_genotype_x_time_significant.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
volcano(df, "interaction genotype x time", file.path(outdir, "plots", "volcano_interaction_genotype_x_time.pdf"))
summary_rows[["interaction"]] <- data.frame(contrast = "interaction_genotype_x_time",
    padj05_any_lfc = sum(df$padj < padj_cut, na.rm = TRUE),
    up = sum(sig$log2FoldChange > 0), down = sum(sig$log2FoldChange < 0))

summ <- do.call(rbind, summary_rows)
write.table(summ, file.path(outdir, "DE_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(summ, row.names = FALSE)

## ---- sample-level QC --------------------------------------------------------
vsd <- vst(dds, blind = TRUE)
pca <- plotPCA(vsd, intgroup = c("genotype", "time"), returnData = TRUE)
pv  <- round(100 * attr(pca, "percentVar"))
pca$well <- meta[rownames(pca), "well"]
p <- ggplot(pca, aes(PC1, PC2, colour = genotype, shape = time)) +
    geom_point(size = 3) +
    geom_text(aes(label = well), size = 2.5, vjust = -1, show.legend = FALSE) +
    labs(x = paste0("PC1: ", pv[1], "% variance"), y = paste0("PC2: ", pv[2], "% variance"),
         title = "PCA (VST, 500 most variable genes)") +
    theme_bw()
ggsave(file.path(outdir, "plots", "PCA.pdf"), p, width = 7, height = 5.5)

d <- as.matrix(dist(t(assay(vsd))))
ann <- data.frame(genotype = meta$genotype, time = meta$time, row.names = meta$sample)
rownames(d) <- colnames(d) <- meta$sample
pdf(file.path(outdir, "plots", "sample_distances.pdf"), width = 9, height = 8)
pheatmap(d, annotation_row = ann, annotation_col = ann,
         labels_row = meta$well, labels_col = meta$well,
         col = colorRampPalette(rev(brewer.pal(9, "Blues")))(100))
dev.off()

nc <- merge(annot, data.frame(gene_id = rownames(norm), norm, check.names = FALSE), by = "gene_id", sort = FALSE)
write.table(nc, file.path(outdir, "normalized_counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(gene_id = rownames(vsd), assay(vsd), check.names = FALSE),
            file.path(outdir, "vst_counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(dds, file.path(outdir, "dds_condition.rds"))
saveRDS(dds_int, file.path(outdir, "dds_interaction.rds"))

writeLines(capture.output(sessionInfo()), file.path(outdir, "R_sessionInfo.txt"))
cat("Done.\n")
