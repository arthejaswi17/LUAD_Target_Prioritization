# ============================================================
# Script 11: Independent Validation in GSE32863
# Project:
# From Transcriptional Dysregulation to Therapeutic Target
# Prioritization in Lung Adenocarcinoma
#
# Purpose:
# Validate the 12 prioritized candidate genes in an independent
# paired LUAD tumor-normal cohort using limma.
#
# Validation dataset:
# GSE32863
# Platform:
# GPL6884 - Illumina HumanWG-6 v3.0
#
# Statistical design:
# Paired tumor vs adjacent normal tissue
# Model: ~ patient + condition
#
# Output:
# - Representative probe mapping
# - Paired limma results
# - TCGA vs GSE32863 concordance
# - Validation classification
# - Main validation plot
# ============================================================


# ============================================================
# 0. PACKAGES
# ============================================================

required_packages <- c(
  "GEOquery",
  "limma",
  "Biobase",
  "ggplot2"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      paste0(
        "Package '", pkg,
        "' is required but not installed."
      )
    )
  }
}

library(GEOquery)
library(limma)
library(Biobase)
library(ggplot2)


# ============================================================
# 1. DIRECTORIES
# ============================================================

project_root <- getwd()

processed_dir <- file.path(
  project_root,
  "data",
  "processed"
)

figure_dir <- file.path(
  project_root,
  "figures"
)

dir.create(
  processed_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

message("Project root: ", project_root)
message("Processed directory: ", processed_dir)

# ============================================================
# 2. LOAD TCGA SHORTLIST
# ============================================================
#
# The 12 genes were defined BEFORE looking at GSE32863.
# They must not be changed based on validation results.
# ============================================================

candidate_genes <- c(
  "AURKA",
  "BRCA1",
  "CDCA5",
  "CDK1",
  "CHEK1",
  "DNA2",
  "GMNN",
  "PLK1",
  "PRKDC",
  "SGO1",
  "SKA1",
  "TOP2A"
)

candidate_entrez <- c(
  AURKA = 6790,
  BRCA1 = 672,
  CDCA5 = 113130,
  CDK1 = 983,
  CHEK1 = 1111,
  DNA2 = 1763,
  GMNN = 51053,
  PLK1 = 5347,
  PRKDC = 5591,
  SGO1 = 151648,
  SKA1 = 220134,
  TOP2A = 7153
)

candidate_ensembl <- c(
  AURKA = "ENSG00000087586",
  BRCA1 = "ENSG00000012048",
  CDCA5 = "ENSG00000146670",
  CDK1  = "ENSG00000170312",
  CHEK1 = "ENSG00000149554",
  DNA2  = "ENSG00000138346",
  GMNN  = "ENSG00000112312",
  PLK1  = "ENSG00000166851",
  PRKDC = "ENSG00000253729",
  SGO1  = "ENSG00000129810",
  SKA1  = "ENSG00000154839",
  TOP2A = "ENSG00000131747"
)


# ============================================================
# 3. LOAD GSE32863
# ============================================================

message("Loading GSE32863...")

gse32863_samples <- getGEO(
  "GSE32863",
  GSEMatrix = TRUE,
  getGPL = FALSE
)

if (length(gse32863_samples) != 1) {
  stop("Expected exactly one ExpressionSet for GSE32863.")
}

eset <- gse32863_samples[[1]]

message(
  "Expression matrix dimensions: ",
  paste(dim(exprs(eset)), collapse = " x ")
)


# ============================================================
# 4. CONSTRUCT SAMPLE METADATA
# ============================================================

pheno <- pData(eset)

# Extract paired sample code from titles.
#
# Example:
# "... 05L12_N (expression)"
# "... 05L12_T (expression)"
#
# The patient/sample code is therefore 05L12.

pheno$sample_code <- sub(
  ".* ([^ ]+)_([NT]) \\(expression\\)$",
  "\\1",
  pheno$title
)

# Determine tissue condition.
pheno$condition <- ifelse(
  grepl("Lung adenocarcinoma", pheno$title),
  "Tumor",
  "Normal"
)

pheno$condition <- factor(
  pheno$condition,
  levels = c("Normal", "Tumor")
)

# Convert patient/sample identifier to factor.
pheno$sample_code <- factor(pheno$sample_code)

message("\nCondition counts:")
print(table(pheno$condition))


# ============================================================
# 5. IDENTIFY GENUINELY PAIRED SAMPLES
# ============================================================

pair_counts <- table(pheno$sample_code)

# A genuine pair must occur exactly twice.
paired_codes <- names(
  pair_counts[pair_counts == 2]
)

paired_pheno <- pheno[
  as.character(pheno$sample_code) %in% paired_codes,
  ,
  drop = FALSE
]

# Re-sort deterministically.
paired_pheno <- paired_pheno[
  order(
    paired_pheno$sample_code,
    paired_pheno$condition
  ),
  ,
  drop = FALSE
]

# ------------------------------------------------------------
# Explicit pairing verification
# ------------------------------------------------------------

pair_table <- table(
  as.character(paired_pheno$sample_code),
  as.character(paired_pheno$condition)
)

message("\nPaired cohort:")
print(table(paired_pheno$condition))

message(
  "Number of paired patients: ",
  length(unique(paired_pheno$sample_code))
)

message(
  "Number of paired samples: ",
  nrow(paired_pheno)
)

# Verify that every retained patient has exactly
# one Normal and one Tumor sample.
if (
  !all(
    rowSums(pair_table) == 2
  )
) {
  stop("Some patients do not have exactly two samples.")
}

if (
  !all(
    pair_table[, "Normal"] == 1
  )
) {
  stop("Some patients do not have exactly one Normal sample.")
}

if (
  !all(
    pair_table[, "Tumor"] == 1
  )
) {
  stop("Some patients do not have exactly one Tumor sample.")
}

message("Pairing check: PASSED")

# ============================================================
# 6. SUBSET EXPRESSION MATRIX TO PAIRED COHORT
# ============================================================

expr_paired <- exprs(eset)[
  ,
  rownames(paired_pheno),
  drop = FALSE
]

if (!all(
  colnames(expr_paired) == rownames(paired_pheno)
)) {
  stop(
    "Expression matrix columns do not match paired metadata."
  )
}

message(
  "\nPaired expression matrix dimensions: ",
  paste(dim(expr_paired), collapse = " x ")
)


# ============================================================
# 7. LOAD GPL6884 ANNOTATION
# ============================================================

message("\nLoading GPL6884 annotation...")

gpl6884 <- getGEO(
  "GPL6884",
  AnnotGPL = TRUE
)

gpl_table <- Table(gpl6884)

message("\nGPL6884 annotation columns:")
print(colnames(gpl_table))

# ------------------------------------------------------------
# GPL6884 annotation fields
# ------------------------------------------------------------
#
# Current GEO annotation uses:
#   ID          = probe ID
#   Gene symbol = gene symbol
#   Gene ID     = Entrez Gene ID
#
# ------------------------------------------------------------

if (!"ID" %in% colnames(gpl_table)) {
  stop("GPL6884 annotation does not contain an 'ID' column.")
}

if (!"Gene symbol" %in% colnames(gpl_table)) {
  stop("GPL6884 annotation does not contain 'Gene symbol'.")
}

if (!"Gene ID" %in% colnames(gpl_table)) {
  stop("GPL6884 annotation does not contain 'Gene ID'.")
}

# Construct standardized annotation table.
#
# We rename the GEO-specific columns to stable internal names
# so the rest of the script does not depend on GEO naming.

probe_map <- data.frame(
  ID = as.character(gpl_table[["ID"]]),
  Symbol = as.character(gpl_table[["Gene symbol"]]),
  Entrez_Gene_ID = as.character(gpl_table[["Gene ID"]]),
  stringsAsFactors = FALSE
)

# Clean Entrez IDs.
probe_map$Entrez_Gene_ID <- trimws(
  probe_map$Entrez_Gene_ID
)

# Keep probes with valid Entrez IDs.
probe_map <- probe_map[
  !is.na(probe_map$Entrez_Gene_ID) &
    probe_map$Entrez_Gene_ID != "",
  ,
  drop = FALSE
]

message(
  "\nNumber of annotated probes with Entrez IDs: ",
  nrow(probe_map)
)


# ============================================================
# 8. MAP THE 12 CANDIDATE GENES TO GPL6884 PROBES
# ============================================================

candidate_entrez_df <- data.frame(
  candidate = names(candidate_entrez),
  entrez = as.character(candidate_entrez),
  stringsAsFactors = FALSE
)

candidate_probe_map <- merge(
  candidate_entrez_df,
  probe_map,
  by.x = "entrez",
  by.y = "Entrez_Gene_ID",
  all.x = TRUE
)

message("\nCandidate probe mapping:")
print(
  candidate_probe_map[
    ,
    c(
      "candidate",
      "entrez",
      "ID",
      "Symbol"
    )
  ]
)


# ============================================================
# 9. CALCULATE MEAN EXPRESSION FOR EACH CANDIDATE PROBE
# ============================================================
#
# IMPORTANT:
#
# Probe selection is performed WITHOUT looking at tumor-normal
# differences.
#
# For genes represented by multiple probes, the probe with the
# highest mean expression across ALL 114 paired samples is chosen.
#
# This rule was defined before validation statistics were examined.
# ============================================================

probe_stats <- candidate_probe_map[
  !is.na(candidate_probe_map$ID),
  ,
  drop = FALSE
]

probe_stats$mean_expression <- sapply(
  probe_stats$ID,
  function(probe_id) {
    mean(
      expr_paired[probe_id, ],
      na.rm = TRUE
    )
  }
)

# Sort within gene by decreasing mean expression.
probe_stats <- probe_stats[
  order(
    probe_stats$entrez,
    -probe_stats$mean_expression
  ),
  ,
  drop = FALSE
]

# Select exactly one representative probe per candidate.
selected_probes <- probe_stats[
  !duplicated(probe_stats$entrez),
  ,
  drop = FALSE
]

# Restore candidate order.
selected_probes <- selected_probes[
  match(
    candidate_genes,
    selected_probes$candidate
  ),
  ,
  drop = FALSE
]

rownames(selected_probes) <- selected_probes$candidate

message("\nSelected representative probes:")
print(
  selected_probes[
    ,
    c(
      "candidate",
      "ID",
      "Symbol",
      "entrez",
      "mean_expression"
    )
  ]
)

# Ensure all candidates have a representative probe.
if (any(is.na(selected_probes$ID))) {
  missing_genes <- selected_probes$candidate[
    is.na(selected_probes$ID)
  ]
  
  stop(
    paste(
      "No GPL6884 probe found for:",
      paste(missing_genes, collapse = ", ")
    )
  )
}


# ============================================================
# 10. BUILD 12-GENE VALIDATION MATRIX
# ============================================================

validation_expr <- expr_paired[
  selected_probes$ID,
  ,
  drop = FALSE
]

rownames(validation_expr) <- selected_probes$candidate

if (!all(
  rownames(validation_expr) == candidate_genes
)) {
  stop("Candidate gene order is incorrect.")
}

if (!all(
  colnames(validation_expr) == rownames(paired_pheno)
)) {
  stop("Validation expression and metadata are not aligned.")
}

message(
  "\nValidation matrix dimensions: ",
  paste(dim(validation_expr), collapse = " x ")
)


# ============================================================
# 11. PAIRED LIMMA MODEL
# ============================================================
#
# Model:
#
# expression ~ patient + condition
#
# The patient/sample_code terms control for patient-specific
# baseline expression differences.
#
# The condition coefficient represents:
#
# Tumor - Normal
# ============================================================

paired_pheno$sample_code <- factor(
  paired_pheno$sample_code
)

paired_pheno$condition <- factor(
  paired_pheno$condition,
  levels = c("Normal", "Tumor")
)

design <- model.matrix(
  ~ sample_code + condition,
  data = paired_pheno
)

message("\nDesign matrix dimensions:")
print(dim(design))

fit <- lmFit(
  validation_expr,
  design
)

fit <- eBayes(fit)

# Extract the Tumor vs Normal coefficient.
validation_results <- topTable(
  fit,
  coef = "conditionTumor",
  number = Inf,
  sort.by = "none"
)

validation_results$gene <- rownames(
  validation_results
)

validation_results <- validation_results[
  match(
    candidate_genes,
    validation_results$gene
  ),
  ,
  drop = FALSE
]

# Explicitly calculate BH FDR across the 12-gene validation panel.
validation_results$FDR <- p.adjust(
  validation_results$P.Value,
  method = "BH"
)

# Rename logFC to make direction explicit.
validation_results$log2FC_Tumor_vs_Normal <-
  validation_results$logFC

# Add representative probe information.
validation_results$probe <- selected_probes[
  validation_results$gene,
  "ID"
]

validation_results$GPL6884_symbol <- selected_probes[
  validation_results$gene,
  "Symbol"
]

validation_results$entrez <- selected_probes[
  validation_results$gene,
  "entrez"
]

# Reorder columns.
validation_results <- validation_results[
  ,
  c(
    "gene",
    "probe",
    "GPL6884_symbol",
    "entrez",
    "log2FC_Tumor_vs_Normal",
    "AveExpr",
    "t",
    "P.Value",
    "FDR",
    "B"
  )
]


# ============================================================
# 12. COMPARE VALIDATION DIRECTION WITH TCGA
# ============================================================
#
# Load the TCGA shrunken DESeq2 results.
# TCGA results use Ensembl gene IDs as row names.
# ============================================================

tcga_file <- file.path(
  processed_dir,
  "TCGA_LUAD_DESeq2_shrunken_results.rds"
)

if (!file.exists(tcga_file)) {
  stop(
    paste(
      "Could not find TCGA shrunken results:",
      tcga_file
    )
  )
}

tcga_results <- readRDS(tcga_file)

# ------------------------------------------------------------
# TCGA candidate gene mapping
# ------------------------------------------------------------

candidate_ensembl <- c(
  AURKA = "ENSG00000087586",
  BRCA1 = "ENSG00000012048",
  CDCA5 = "ENSG00000146670",
  CDK1  = "ENSG00000170312",
  CHEK1 = "ENSG00000149554",
  DNA2  = "ENSG00000138346",
  GMNN  = "ENSG00000112312",
  PLK1  = "ENSG00000166851",
  PRKDC = "ENSG00000253729",
  SGO1  = "ENSG00000129810",
  SKA1  = "ENSG00000154839",
  TOP2A = "ENSG00000131747"
)

# ------------------------------------------------------------
# Identify the shrunken log2FC column
# ------------------------------------------------------------

lfc_candidates <- c(
  "log2FoldChange",
  "log2FC",
  "lfc"
)

lfc_col <- lfc_candidates[
  lfc_candidates %in% colnames(tcga_results)
][1]

if (is.na(lfc_col)) {
  stop(
    "Could not identify TCGA log2FC column."
  )
}

# ------------------------------------------------------------
# Match candidate genes using Ensembl IDs
# ------------------------------------------------------------

tcga_ensembl <- rownames(tcga_results)

tcga_idx <- match(
  candidate_ensembl,
  tcga_ensembl
)

# Check that all candidates were found
if (any(is.na(tcga_idx))) {
  
  missing_genes <- names(candidate_ensembl)[
    is.na(tcga_idx)
  ]
  
  stop(
    paste(
      "The following candidate genes were not found in TCGA:",
      paste(missing_genes, collapse = ", ")
    )
  )
}

# ------------------------------------------------------------
# Extract TCGA shrunken log2FC
# ------------------------------------------------------------

tcga_log2fc <- tcga_results[
  tcga_idx,
  lfc_col
]

# ------------------------------------------------------------
# Build TCGA lookup table
# ------------------------------------------------------------

tcga_lookup <- data.frame(
  gene = names(candidate_ensembl),
  ensembl_id = unname(candidate_ensembl),
  TCGA_shrunken_log2FC = as.numeric(tcga_log2fc),
  stringsAsFactors = FALSE
)

print(tcga_lookup)

# ------------------------------------------------------------
# Merge TCGA results with validation results
# ------------------------------------------------------------

validation_results <- merge(
  validation_results,
  tcga_lookup,
  by = "gene",
  all.x = TRUE,
  sort = FALSE
)

# Restore original candidate-gene order
validation_results <- validation_results[
  match(
    candidate_genes,
    validation_results$gene
  ),
  ,
  drop = FALSE
]


# ============================================================
# 13. VALIDATION CLASSIFICATION
# ============================================================

validation_results$TCGA_direction <- ifelse(
  validation_results$TCGA_shrunken_log2FC > 0,
  "Up",
  "Down"
)

validation_results$Validation_direction <- ifelse(
  validation_results$log2FC_Tumor_vs_Normal > 0,
  "Up",
  "Down"
)

validation_results$direction_concordant <-
  validation_results$TCGA_direction ==
  validation_results$Validation_direction

# Primary validation definition:
#
# Significant in GSE32863 AND same direction as TCGA.
#
# FDR threshold:
# BH-adjusted FDR < 0.05

validation_results$validation_class <- ifelse(
  validation_results$direction_concordant &
    validation_results$FDR < 0.05,
  "Validated",
  ifelse(
    validation_results$direction_concordant,
    "Directionally consistent",
    "Discordant"
  )
)


# ============================================================
# 14. DISPLAY FINAL VALIDATION TABLE
# ============================================================

final_validation_table <- validation_results[
  ,
  c(
    "gene",
    "probe",
    "GPL6884_symbol",
    "TCGA_shrunken_log2FC",
    "log2FC_Tumor_vs_Normal",
    "P.Value",
    "FDR",
    "direction_concordant",
    "validation_class"
  )
]

message("\n============================================")
message("FINAL GSE32863 VALIDATION RESULTS")
message("============================================")

print(
  final_validation_table,
  row.names = FALSE
)


# ============================================================
# 15. SUMMARY COUNTS
# ============================================================

n_validated <- sum(
  validation_results$validation_class ==
    "Validated",
  na.rm = TRUE
)

n_consistent <- sum(
  validation_results$validation_class ==
    "Directionally consistent",
  na.rm = TRUE
)

n_discordant <- sum(
  validation_results$validation_class ==
    "Discordant",
  na.rm = TRUE
)

n_concordant <- sum(
  validation_results$direction_concordant,
  na.rm = TRUE
)

message("\nValidation summary:")
message(
  "Directionally concordant: ",
  n_concordant,
  "/",
  length(candidate_genes)
)

message(
  "Validated (same direction + FDR < 0.05): ",
  n_validated
)

message(
  "Directionally consistent but FDR >= 0.05: ",
  n_consistent
)

message(
  "Discordant: ",
  n_discordant
)


# ============================================================
# 16. CORRELATION OF TCGA AND VALIDATION EFFECT SIZES
# ============================================================
#
# This is descriptive, not a causal/statistical proof.
# ============================================================

complete_effects <- validation_results[
  complete.cases(
    validation_results$TCGA_shrunken_log2FC,
    validation_results$log2FC_Tumor_vs_Normal
  ),
  ,
  drop = FALSE
]

effect_correlation <- cor.test(
  complete_effects$TCGA_shrunken_log2FC,
  complete_effects$log2FC_Tumor_vs_Normal,
  method = "spearman",
  exact = FALSE
)

message("\nTCGA vs GSE32863 effect-size correlation:")
print(effect_correlation)


# ============================================================
# 17. SAVE RESULTS
# ============================================================

saveRDS(
  selected_probes,
  file.path(
    processed_dir,
    "GSE32863_selected_probes.rds"
  )
)

saveRDS(
  validation_expr,
  file.path(
    processed_dir,
    "GSE32863_validation_expression_12genes.rds"
  )
)

saveRDS(
  paired_pheno,
  file.path(
    processed_dir,
    "GSE32863_paired_metadata.rds"
  )
)

saveRDS(
  validation_results,
  file.path(
    processed_dir,
    "GSE32863_validation_results.rds"
  )
)

write.csv(
  final_validation_table,
  file.path(
    processed_dir,
    "GSE32863_validation_results.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 18. VALIDATION EFFECT-SIZE PLOT
# ============================================================

plot_data <- validation_results

plot_data$gene <- factor(
  plot_data$gene,
  levels = candidate_genes
)

p <- ggplot(
  plot_data,
  aes(
    x = TCGA_shrunken_log2FC,
    y = log2FC_Tumor_vs_Normal,
    label = gene
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  geom_point(
    size = 3
  ) +
  geom_text(
    vjust = -0.7,
    size = 3
  ) +
  labs(
    title = "Independent Validation of Prioritized LUAD Targets",
    subtitle = "TCGA-LUAD vs paired GSE32863 tumor-normal effect sizes",
    x = "TCGA shrunken log2FC (Tumor vs Normal)",
    y = "GSE32863 limma log2FC (Tumor vs Normal)"
  ) +
  theme_bw()

ggsave(
  filename = file.path(
    figure_dir,
    "Figure_GSE32863_TCGA_effect_size_validation.png"
  ),
  plot = p,
  width = 8,
  height = 7,
  dpi = 300
)


# ============================================================
# 19. SESSION INFORMATION
# ============================================================

session_file <- file.path(
  processed_dir,
  "GSE32863_validation_sessionInfo.txt"
)

sink(session_file)

cat(
  "Script 11: GSE32863 Independent Validation\n\n"
)

cat(
  "Number of paired patients: ",
  length(unique(paired_pheno$sample_code)),
  "\n\n"
)

cat(
  "Number of candidate genes: ",
  length(candidate_genes),
  "\n\n"
)

cat(
  "Validated candidates: ",
  n_validated,
  "\n"
)

cat(
  "Directionally consistent candidates: ",
  n_consistent,
  "\n"
)

cat(
  "Discordant candidates: ",
  n_discordant,
  "\n\n"
)

print(sessionInfo())

sink()


# ============================================================
# 20. FINAL MESSAGE
# ============================================================

message("\n============================================")
message("SCRIPT 11 COMPLETED")
message("============================================")

message(
  "Results saved to: ",
  processed_dir
)

message(
  "Validation figure saved to: ",
  file.path(
    figure_dir,
    "Figure_GSE32863_TCGA_effect_size_validation.png"
  )
)

message("\nNext step:")
message(
  "Inspect the validation results before proceeding to "
)
message(
  "final target-evidence integration."
)