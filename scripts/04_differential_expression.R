# ============================================================
# 04_differential_expression.R
# Differential Expression Analysis of TCGA-LUAD RNA-seq
#
# Project:
# From Transcriptional Dysregulation to Therapeutic Target
# Prioritization in Lung Adenocarcinoma
#
# Design:
#   Paired tumor-normal analysis
#   ~ patient + condition
#
# Input:
#   data/processed/TCGA_LUAD_raw_counts.rds
#   data/processed/paired_metadata.rds
#
# Outputs:
#   data/processed/TCGA_LUAD_dds_DESeq2.rds
#   data/processed/TCGA_LUAD_DESeq2_raw_results.rds
#   data/processed/TCGA_LUAD_DESeq2_shrunken_results.rds
#
# ============================================================


# ------------------------------------------------------------
# 1. Setup
# ------------------------------------------------------------

project_root <- getwd()

processed_dir <- file.path(
  project_root,
  "data",
  "processed"
)

figures_dir <- file.path(
  project_root,
  "figures"
)

dir.create(
  processed_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  figures_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


required_packages <- c(
  "DESeq2",
  "apeglm",
  "ggplot2"
)

for (pkg in required_packages) {
  
  if (!requireNamespace(pkg, quietly = TRUE)) {
    
    stop(
      "Required package not installed: ",
      pkg,
      "\nPlease install it before running this script."
    )
    
  }
}

library(DESeq2)
library(apeglm)
library(ggplot2)


# ------------------------------------------------------------
# 2. Load count matrix and metadata
# ------------------------------------------------------------

counts <- readRDS(
  file.path(
    processed_dir,
    "TCGA_LUAD_raw_counts.rds"
  )
)

metadata <- readRDS(
  file.path(
    processed_dir,
    "paired_metadata.rds"
  )
)


cat("Count matrix dimensions:\n")
print(dim(counts))

cat("\nMetadata dimensions:\n")
print(dim(metadata))


# ------------------------------------------------------------
# 3. Verify sample alignment
# ------------------------------------------------------------

if (!all(colnames(counts) %in% metadata$sample)) {
  
  missing_samples <- setdiff(
    colnames(counts),
    metadata$sample
  )
  
  stop(
    "Samples missing from metadata:\n",
    paste(
      missing_samples,
      collapse = ", "
    )
  )
}


metadata <- metadata[
  match(
    colnames(counts),
    metadata$sample
  ),
  ,
  drop = FALSE
]


if (!all(metadata$sample == colnames(counts))) {
  
  stop(
    "Metadata and count matrix are not aligned."
  )
  
}


cat(
  "\nSample alignment check: PASSED\n"
)


# ------------------------------------------------------------
# 4. Verify paired design
# ------------------------------------------------------------

metadata$condition <- factor(
  metadata$condition,
  levels = c(
    "Normal",
    "Tumor"
  )
)

metadata$patient <- factor(
  metadata$patient
)


pair_table <- table(
  metadata$patient,
  metadata$condition
)


complete_pairs <- sum(
  pair_table[, "Normal"] == 1 &
    pair_table[, "Tumor"] == 1
)


cat(
  "\nPaired design:\n"
)

cat(
  "Patients:",
  nrow(pair_table),
  "\n"
)

cat(
  "Complete tumor-normal pairs:",
  complete_pairs,
  "\n"
)


if (
  complete_pairs != nrow(pair_table)
) {
  
  stop(
    "Incomplete or duplicated tumor-normal pairs detected."
  )
  
}


# ------------------------------------------------------------
# 5. Construct DESeq2 object
#
# The patient term accounts for the matched nature of the
# tumor-normal samples.
#
# condition_Tumor_vs_Normal is the biological contrast of
# interest.
# ------------------------------------------------------------

dds <- DESeqDataSetFromMatrix(
  
  countData = round(
    as.matrix(counts)
  ),
  
  colData = metadata,
  
  design = ~ patient + condition
  
)


cat(
  "\nDESeqDataSet constructed.\n"
)

cat(
  "Genes:",
  nrow(dds),
  "\n"
)

cat(
  "Samples:",
  ncol(dds),
  "\n"
)


# ------------------------------------------------------------
# 6. Independent filtering
#
# DESeq2 performs independent filtering of low-count genes
# during results() by default.
#
# We do not impose an arbitrary CPM/TPM threshold here.
# ------------------------------------------------------------


# ------------------------------------------------------------
# 7. Run DESeq2
# ------------------------------------------------------------

cat(
  "\nRunning DESeq2...\n"
)

dds <- DESeq(
  dds
)


cat(
  "DESeq2 completed successfully.\n"
)


# ------------------------------------------------------------
# 8. Inspect available coefficients
# ------------------------------------------------------------

cat(
  "\nAvailable model coefficients:\n"
)

print(
  resultsNames(dds)
)


# Expected coefficient:
#
# condition_Tumor_vs_Normal
#


if (
  !"condition_Tumor_vs_Normal" %in%
  resultsNames(dds)
) {
  
  stop(
    "Expected coefficient 'condition_Tumor_vs_Normal' was not found."
  )
  
}


# ------------------------------------------------------------
# 9. Extract raw DESeq2 results
#
# Contrast:
#   Tumor vs Normal
#
# log2FoldChange > 0 = higher expression in Tumor
# log2FoldChange < 0 = lower expression in Tumor
#
# padj = Benjamini-Hochberg adjusted p-value
# ------------------------------------------------------------

res_raw <- results(
  dds,
  name = "condition_Tumor_vs_Normal",
  alpha = 0.05
)


res_raw_df <- as.data.frame(
  res_raw
)


res_raw_df$gene <- rownames(
  res_raw_df
)


# Reorder columns

res_raw_df <- res_raw_df[
  ,
  c(
    "gene",
    "baseMean",
    "log2FoldChange",
    "lfcSE",
    "pvalue",
    "padj"
  )
]


# Save raw results

saveRDS(
  res_raw_df,
  file.path(
    processed_dir,
    "TCGA_LUAD_DESeq2_raw_results.rds"
  )
)


write.csv(
  res_raw_df,
  file.path(
    processed_dir,
    "TCGA_LUAD_DESeq2_raw_results.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 10. Summarize raw differential expression
# ------------------------------------------------------------

raw_significant <- res_raw_df[
  !is.na(res_raw_df$padj) &
    res_raw_df$padj < 0.05,
  ,
  drop = FALSE
]


raw_up <- raw_significant[
  raw_significant$log2FoldChange > 0,
  ,
  drop = FALSE
]


raw_down <- raw_significant[
  raw_significant$log2FoldChange < 0,
  ,
  drop = FALSE
]


cat(
  "\nRaw DESeq2 results:\n"
)

cat(
  "Genes tested:",
  sum(!is.na(res_raw_df$pvalue)),
  "\n"
)

cat(
  "FDR < 0.05:",
  nrow(raw_significant),
  "\n"
)

cat(
  "Upregulated:",
  nrow(raw_up),
  "\n"
)

cat(
  "Downregulated:",
  nrow(raw_down),
  "\n"
)


# ------------------------------------------------------------
# 11. Shrink log2 fold changes using apeglm
#
# Shrinkage improves interpretability of effect sizes,
# particularly for genes with limited information.
#
# Significance testing remains based on the DESeq2 model;
# shrinkage is used to obtain more stable effect-size
# estimates.
# ------------------------------------------------------------

cat(
  "\nApplying apeglm log2 fold-change shrinkage...\n"
)


res_shrunken <- lfcShrink(
  dds,
  coef = "condition_Tumor_vs_Normal",
  type = "apeglm"
)


res_shrunken_df <- as.data.frame(
  res_shrunken
)


res_shrunken_df$gene <- rownames(
  res_shrunken_df
)


# Keep the same useful columns

res_shrunken_df <- res_shrunken_df[
  ,
  c(
    "gene",
    "baseMean",
    "log2FoldChange",
    "lfcSE",
    "pvalue",
    "padj"
  )
]


# Save shrunken results

saveRDS(
  res_shrunken_df,
  file.path(
    processed_dir,
    "TCGA_LUAD_DESeq2_shrunken_results.rds"
  )
)


write.csv(
  res_shrunken_df,
  file.path(
    processed_dir,
    "TCGA_LUAD_DESeq2_shrunken_results.csv"
  ),
  row.names = FALSE
)


cat(
  "apeglm shrinkage completed successfully.\n"
)


# ------------------------------------------------------------
# 12. Define practical differential-expression set
#
# Criteria:
#
#   padj < 0.05
#   |shrunken log2FC| >= 1
#
# This corresponds to an FDR-controlled result with at
# least approximately a two-fold expression difference.
#
# The threshold is used as a transparent practical
# prioritization criterion, not as a claim that genes below
# it are biologically irrelevant.
# ------------------------------------------------------------

res_shrunken_df$DE_category <- "Not significant"


res_shrunken_df$DE_category[
  !is.na(res_shrunken_df$padj) &
    res_shrunken_df$padj < 0.05 &
    res_shrunken_df$log2FoldChange >= 1
] <- "Upregulated"


res_shrunken_df$DE_category[
  !is.na(res_shrunken_df$padj) &
    res_shrunken_df$padj < 0.05 &
    res_shrunken_df$log2FoldChange <= -1
] <- "Downregulated"


practical_de <- res_shrunken_df[
  res_shrunken_df$DE_category != "Not significant",
  ,
  drop = FALSE
]


practical_up <- practical_de[
  practical_de$DE_category == "Upregulated",
  ,
  drop = FALSE
]


practical_down <- practical_de[
  practical_de$DE_category == "Downregulated",
  ,
  drop = FALSE
]


cat(
  "\nPractical DE set:\n"
)

cat(
  "FDR < 0.05 and |shrunken log2FC| >= 1:",
  nrow(practical_de),
  "\n"
)

cat(
  "Upregulated:",
  nrow(practical_up),
  "\n"
)

cat(
  "Downregulated:",
  nrow(practical_down),
  "\n"
)


# ------------------------------------------------------------
# 13. Save DESeq2 object
# ------------------------------------------------------------

saveRDS(
  dds,
  file.path(
    processed_dir,
    "TCGA_LUAD_dds_DESeq2.rds"
  )
)


# ------------------------------------------------------------
# 14. Volcano plot
# ------------------------------------------------------------

volcano_df <- res_shrunken_df

volcano_df$neg_log10_padj <- -log10(
  pmax(
    volcano_df$padj,
    .Machine$double.xmin,
    na.rm = FALSE
  )
)


volcano_df$plot_category <- "Not significant"


volcano_df$plot_category[
  volcano_df$DE_category == "Upregulated"
] <- "Upregulated"


volcano_df$plot_category[
  volcano_df$DE_category == "Downregulated"
] <- "Downregulated"


volcano_df$plot_category <- factor(
  volcano_df$plot_category,
  levels = c(
    "Downregulated",
    "Not significant",
    "Upregulated"
  )
)


p_volcano <- ggplot(
  volcano_df,
  aes(
    x = log2FoldChange,
    y = neg_log10_padj,
    color = plot_category
  )
) +
  
  geom_point(
    alpha = 0.55,
    size = 1.2,
    na.rm = TRUE
  ) +
  
  geom_vline(
    xintercept = c(
      -1,
      1
    ),
    linetype = "dashed"
  ) +
  
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed"
  ) +
  
  scale_color_manual(
    values = c(
      "Downregulated" = "#377EB8",
      "Not significant" = "grey75",
      "Upregulated" = "#E41A1C"
    )
  ) +
  
  labs(
    title = "TCGA-LUAD differential expression",
    x = "Shrunken log2 fold change",
    y = expression(-log[10]("adjusted p-value")),
    color = NULL
  ) +
  
  theme_bw()


ggsave(
  file.path(
    figures_dir,
    "TCGA_LUAD_volcano.pdf"
  ),
  p_volcano,
  width = 8,
  height = 6
)


# ------------------------------------------------------------
# 15. Final report
# ------------------------------------------------------------

cat(
  "\n============================================================\n"
)

cat(
  "DIFFERENTIAL EXPRESSION ANALYSIS COMPLETED\n"
)

cat(
  "============================================================\n"
)

cat(
  "Samples:",
  ncol(dds),
  "\n"
)

cat(
  "Genes:",
  nrow(dds),
  "\n"
)

cat(
  "Patients:",
  nrow(pair_table),
  "\n"
)

cat(
  "Complete pairs:",
  complete_pairs,
  "\n"
)

cat(
  "Raw FDR < 0.05:",
  nrow(raw_significant),
  "\n"
)

cat(
  "Upregulated:",
  nrow(raw_up),
  "\n"
)

cat(
  "Downregulated:",
  nrow(raw_down),
  "\n"
)

cat(
  "Practical DE set:",
  nrow(practical_de),
  "\n"
)

cat(
  "Practical upregulated:",
  nrow(practical_up),
  "\n"
)

cat(
  "Practical downregulated:",
  nrow(practical_down),
  "\n"
)

cat(
  "\nOutputs saved successfully.\n"
)

cat(
  "============================================================\n"
)