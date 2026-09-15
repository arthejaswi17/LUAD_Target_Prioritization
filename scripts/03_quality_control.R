# ============================================================
# 03_quality_control.R
# Quality Control and Exploratory Analysis of TCGA-LUAD RNA-seq
#
# Project:
# From Transcriptional Dysregulation to Therapeutic Target
# Prioritization in Lung Adenocarcinoma
#
# Input:
#   data/processed/TCGA_LUAD_raw_counts.rds
#   data/processed/paired_metadata.rds
#
# Outputs:
#   data/processed/TCGA_LUAD_vsd.rds
#   data/processed/TCGA_LUAD_qc_summary.rds
#
# Figures:
#   figures/TCGA_LUAD_library_size.pdf
#   figures/TCGA_LUAD_detected_genes.pdf
#   figures/TCGA_LUAD_library_vs_detected.pdf
#   figures/TCGA_LUAD_PCA.pdf
#   figures/TCGA_LUAD_sample_correlation.pdf
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


# Required packages

required_packages <- c(
  "DESeq2",
  "ggplot2",
  "pheatmap",
  "scales"
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
library(ggplot2)
library(pheatmap)
library(scales)


# ------------------------------------------------------------
# 2. Load data
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
# 3. Basic integrity checks
# ------------------------------------------------------------

if (is.null(rownames(counts))) {
  stop("Count matrix has no gene identifiers.")
}

if (is.null(colnames(counts))) {
  stop("Count matrix has no sample identifiers.")
}

if (anyDuplicated(rownames(counts)) > 0) {
  stop("Duplicated gene identifiers detected.")
}

if (anyDuplicated(colnames(counts)) > 0) {
  stop("Duplicated sample identifiers detected.")
}


# Verify that all count-matrix samples occur in metadata

if (!all(colnames(counts) %in% metadata$sample)) {
  
  missing_samples <- setdiff(
    colnames(counts),
    metadata$sample
  )
  
  stop(
    "Some count-matrix samples are missing from metadata:\n",
    paste(missing_samples, collapse = ", ")
  )
}


# Align metadata to count-matrix column order

metadata <- metadata[
  match(
    colnames(counts),
    metadata$sample
  ),
  ,
  drop = FALSE
]


# Confirm exact alignment

if (!all(metadata$sample == colnames(counts))) {
  stop(
    "Metadata could not be aligned exactly to count-matrix columns."
  )
}

cat("\nSample alignment check: PASSED\n")


# ------------------------------------------------------------
# 4. Count matrix validity
# ------------------------------------------------------------

missing_count_values <- sum(
  is.na(counts)
)

negative_count_values <- sum(
  counts < 0
)

non_integer_values <- sum(
  counts != round(counts)
)


cat("\nCount matrix validity:\n")

cat(
  "Missing values:",
  missing_count_values,
  "\n"
)

cat(
  "Negative values:",
  negative_count_values,
  "\n"
)

cat(
  "Non-integer values:",
  non_integer_values,
  "\n"
)


if (
  missing_count_values > 0 ||
  negative_count_values > 0 ||
  non_integer_values > 0
) {
  
  stop(
    "Count matrix contains invalid count values."
  )
  
}


# ------------------------------------------------------------
# 5. Library size
# ------------------------------------------------------------

library_size <- colSums(counts)

library_size_summary <- summary(
  library_size
)

cat("\nLibrary size summary:\n")
print(library_size_summary)


library_size_df <- data.frame(
  
  sample_id = names(library_size),
  
  library_size = as.numeric(
    library_size
  ),
  
  condition = metadata$condition,
  
  patient_id = metadata$patient
  
)


# Library size plot

p_library <- ggplot(
  library_size_df,
  aes(
    x = reorder(
      sample_id,
      library_size
    ),
    y = library_size,
    fill = condition
  )
) +
  
  geom_col() +
  
  scale_y_continuous(
    labels = scales::comma
  ) +
  
  labs(
    title = "TCGA-LUAD RNA-seq library sizes",
    x = "Sample",
    y = "Total mapped reads"
  ) +
  
  theme_bw() +
  
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )


ggsave(
  file.path(
    figures_dir,
    "TCGA_LUAD_library_size.pdf"
  ),
  p_library,
  width = 8,
  height = 5
)


# ------------------------------------------------------------
# 6. Number of detected genes
#
# A gene is considered detected when its raw count > 0.
#
# This is a descriptive QC metric and is NOT used as an
# arbitrary sample-filtering criterion.
# ------------------------------------------------------------

genes_detected <- colSums(
  counts > 0
)

genes_detected_summary <- summary(
  genes_detected
)


cat("\nDetected genes per sample:\n")
print(genes_detected_summary)


detected_df <- data.frame(
  
  sample_id = names(
    genes_detected
  ),
  
  genes_detected = as.numeric(
    genes_detected
  ),
  
  condition = metadata$condition,
  
  patient_id = metadata$patient
  
)


# Detected genes plot

p_detected <- ggplot(
  detected_df,
  aes(
    x = reorder(
      sample_id,
      genes_detected
    ),
    y = genes_detected,
    fill = condition
  )
) +
  
  geom_col() +
  
  scale_y_continuous(
    labels = scales::comma
  ) +
  
  labs(
    title = "Genes detected per TCGA-LUAD sample",
    x = "Sample",
    y = "Genes with count > 0"
  ) +
  
  theme_bw() +
  
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )


ggsave(
  file.path(
    figures_dir,
    "TCGA_LUAD_detected_genes.pdf"
  ),
  p_detected,
  width = 8,
  height = 5
)


# ------------------------------------------------------------
# 7. Library size versus detected genes
# ------------------------------------------------------------

library_detected_df <- data.frame(
  
  sample_id = colnames(counts),
  
  library_size = as.numeric(
    library_size
  ),
  
  genes_detected = as.numeric(
    genes_detected
  ),
  
  condition = metadata$condition
  
)


p_library_detected <- ggplot(
  library_detected_df,
  aes(
    x = library_size,
    y = genes_detected,
    color = condition
  )
) +
  
  geom_point(
    size = 2.5,
    alpha = 0.8
  ) +
  
  scale_x_continuous(
    labels = scales::comma
  ) +
  
  scale_y_continuous(
    labels = scales::comma
  ) +
  
  labs(
    title = "Library size versus genes detected",
    x = "Library size",
    y = "Genes detected"
  ) +
  
  theme_bw()


ggsave(
  file.path(
    figures_dir,
    "TCGA_LUAD_library_vs_detected.pdf"
  ),
  p_library_detected,
  width = 7,
  height = 5
)


# ------------------------------------------------------------
# 8. Construct DESeq2 object for exploratory transformation
#
# Primary study design:
#
#   ~ patient + condition
#
# This accounts for the matched tumor-normal structure.
#
# IMPORTANT:
# The VST generated here is used for exploratory QC only.
# It is NOT used as input for differential expression testing.
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


dds_qc <- DESeqDataSetFromMatrix(
  
  countData = round(
    as.matrix(counts)
  ),
  
  colData = metadata,
  
  design = ~ patient + condition
  
)


cat(
  "\nDESeqDataSet created successfully.\n"
)

cat(
  "Dimensions:",
  nrow(dds_qc),
  "genes x",
  ncol(dds_qc),
  "samples\n"
)


# ------------------------------------------------------------
# 9. Variance-stabilizing transformation
# ------------------------------------------------------------

vsd <- vst(
  dds_qc,
  blind = TRUE
)


saveRDS(
  vsd,
  file.path(
    processed_dir,
    "TCGA_LUAD_vsd.rds"
  )
)


cat(
  "\nVST transformation completed.\n"
)


# ------------------------------------------------------------
# 10. Principal component analysis
# ------------------------------------------------------------

pca_data <- plotPCA(
  vsd,
  intgroup = "condition",
  returnData = TRUE
)


percent_var <- round(
  100 *
    attr(
      pca_data,
      "percentVar"
    ),
  1
)


pca_plot <- ggplot(
  pca_data,
  aes(
    x = PC1,
    y = PC2,
    color = condition
  )
) +
  
  geom_point(
    size = 3,
    alpha = 0.85
  ) +
  
  labs(
    title = "TCGA-LUAD principal component analysis",
    
    x = paste0(
      "PC1 (",
      percent_var[1],
      "%)"
    ),
    
    y = paste0(
      "PC2 (",
      percent_var[2],
      "%)"
    )
  ) +
  
  theme_bw()


ggsave(
  file.path(
    figures_dir,
    "TCGA_LUAD_PCA.pdf"
  ),
  pca_plot,
  width = 7,
  height = 5.5
)


cat(
  "\nPCA variance explained:\n"
)

cat(
  "PC1:",
  percent_var[1],
  "%\n"
)

cat(
  "PC2:",
  percent_var[2],
  "%\n"
)


# ------------------------------------------------------------
# 11. Sample-sample correlation
#
# Pearson correlation is calculated on VST-transformed
# expression values.
# ------------------------------------------------------------

vsd_matrix <- assay(
  vsd
)


sample_cor <- cor(
  vsd_matrix,
  method = "pearson"
)


cor_values <- sample_cor[
  upper.tri(sample_cor)
]


cat(
  "\nSample correlation summary:\n"
)

print(
  summary(cor_values)
)


# Correlation heatmap

pdf(
  file.path(
    figures_dir,
    "TCGA_LUAD_sample_correlation.pdf"
  ),
  width = 9,
  height = 9
)


pheatmap(
  sample_cor,
  
  show_rownames = FALSE,
  
  show_colnames = FALSE,
  
  clustering_distance_rows = "correlation",
  
  clustering_distance_cols = "correlation",
  
  main = "TCGA-LUAD sample correlation"
)


dev.off()


# ------------------------------------------------------------
# 12. Tumor-normal pair integrity
# ------------------------------------------------------------

pair_table <- table(
  metadata$patient,
  metadata$condition
)


cat(
  "\nPaired design check:\n"
)

cat(
  "Patients:",
  nrow(pair_table),
  "\n"
)


if (
  !all(
    c(
      "Normal",
      "Tumor"
    ) %in%
    colnames(pair_table)
  )
) {
  
  stop(
    "Both Normal and Tumor conditions are required."
  )
  
}


complete_pairs <- sum(
  
  pair_table[, "Normal"] == 1 &
    
    pair_table[, "Tumor"] == 1
  
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
    "Not every patient has exactly one Normal and one Tumor sample."
  )
  
}


cat(
  "Pair integrity check: PASSED\n"
)


# ------------------------------------------------------------
# 13. Save QC summary
# ------------------------------------------------------------

qc_summary <- list(
  
  count_matrix_dimensions = dim(
    counts
  ),
  
  missing_values = missing_count_values,
  
  negative_values = negative_count_values,
  
  non_integer_values = non_integer_values,
  
  library_size = library_size,
  
  library_size_summary = library_size_summary,
  
  genes_detected = genes_detected,
  
  genes_detected_summary = genes_detected_summary,
  
  pca_percent_variance = percent_var,
  
  sample_correlation_summary = summary(
    cor_values
  ),
  
  complete_pairs = complete_pairs,
  
  total_patients = nrow(
    pair_table
  ),
  
  session_info = sessionInfo()
  
)


saveRDS(
  qc_summary,
  file.path(
    processed_dir,
    "TCGA_LUAD_qc_summary.rds"
  )
)


# ------------------------------------------------------------
# 14. Final report
# ------------------------------------------------------------

cat(
  "\n============================================================\n"
)

cat(
  "QC ANALYSIS COMPLETED\n"
)

cat(
  "============================================================\n"
)


cat(
  "Samples:",
  ncol(counts),
  "\n"
)

cat(
  "Genes:",
  nrow(counts),
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
  "Library size range:",
  format(
    min(library_size),
    big.mark = ","
  ),
  "to",
  format(
    max(library_size),
    big.mark = ","
  ),
  "\n"
)


cat(
  "Detected genes range:",
  format(
    min(genes_detected),
    big.mark = ","
  ),
  "to",
  format(
    max(genes_detected),
    big.mark = ","
  ),
  "\n"
)


cat(
  "PCA: PC1 =",
  percent_var[1],
  "%; PC2 =",
  percent_var[2],
  "%\n"
)


cat(
  "\nQC outputs saved successfully.\n"
)

cat(
  "============================================================\n"
)