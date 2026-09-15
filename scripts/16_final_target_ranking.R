# ============================================================
# SCRIPT 16: FINAL TARGET RANKING AND BIOLOGICAL SELECTION
# ============================================================

cat("\n")
cat("============================================================\n")
cat("SCRIPT 16: FINAL TARGET RANKING AND BIOLOGICAL SELECTION\n")
cat("============================================================\n")

# ------------------------------------------------------------
# 0. PACKAGE SETUP
# ------------------------------------------------------------

required_packages <- c(
  "dplyr",
  "ggplot2",
  "readr"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      paste0(
        "Required package '", pkg,
        "' is not installed. Please install it before running Script 16."
      )
    )
  }
}

library(dplyr)
library(ggplot2)
library(readr)

# ------------------------------------------------------------
# 1. DIRECTORIES
# ------------------------------------------------------------

project_dir <- getwd()

processed_dir <- file.path(
  project_dir,
  "data",
  "processed"
)

figures_dir <- file.path(
  project_dir,
  "figures"
)

if (!dir.exists(processed_dir)) {
  dir.create(processed_dir, recursive = TRUE)
}

if (!dir.exists(figures_dir)) {
  dir.create(figures_dir, recursive = TRUE)
}

cat("Project directory:", project_dir, "\n")
cat("Processed directory:", processed_dir, "\n")
cat("Figures directory:", figures_dir, "\n")

# ------------------------------------------------------------
# 2. SAFE FILE LOADING FUNCTION
# ------------------------------------------------------------

load_rds_safe <- function(path, description) {
  
  if (!file.exists(path)) {
    stop(
      paste0(
        "\nMissing required file:\n",
        path,
        "\n\nRequired for: ",
        description,
        "\n"
      )
    )
  }
  
  obj <- readRDS(path)
  
  cat(
    description,
    "loaded:",
    if (is.data.frame(obj)) nrow(obj) else "object",
    "\n"
  )
  
  return(obj)
}

# ------------------------------------------------------------
# 3. LOAD SCRIPT 15 RESULTS
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("LOADING SCRIPT 15 RESULTS\n")
cat("============================================================\n")

integration_path <- file.path(
  processed_dir,
  "final_biological_evidence_integration.rds"
)

summary_path <- file.path(
  processed_dir,
  "final_target_evidence_summary.rds"
)

integration <- load_rds_safe(
  integration_path,
  "Integrated biological evidence"
)

target_summary <- load_rds_safe(
  summary_path,
  "Final target evidence summary"
)

# ------------------------------------------------------------
# 4. VALIDATE INPUT
# ------------------------------------------------------------

cat("\nValidating Script 15 results...\n")

if (!is.data.frame(integration)) {
  stop("final_biological_evidence_integration.rds is not a data.frame.")
}

if (!is.data.frame(target_summary)) {
  stop("final_target_evidence_summary.rds is not a data.frame.")
}

cat("Integration rows:", nrow(integration), "\n")
cat("Integration columns:", ncol(integration), "\n")

cat("Summary rows:", nrow(target_summary), "\n")
cat("Summary columns:", ncol(target_summary), "\n")

# ------------------------------------------------------------
# 5. IDENTIFY SYMBOL COLUMN
# ------------------------------------------------------------

find_symbol_column <- function(df) {
  
  candidates <- c(
    "SYMBOL",
    "gene_symbol",
    "symbol",
    "gene"
  )
  
  found <- candidates[candidates %in% names(df)]
  
  if (length(found) == 0) {
    stop(
      "Could not identify gene-symbol column in Script 15 output."
    )
  }
  
  found[1]
}

symbol_col <- find_symbol_column(integration)

cat("Gene symbol column:", symbol_col, "\n")

# ------------------------------------------------------------
# 6. STANDARDIZE SYMBOL COLUMN
# ------------------------------------------------------------

integration$FINAL_SYMBOL <- as.character(
  integration[[symbol_col]]
)

integration$FINAL_SYMBOL <- trimws(
  integration$FINAL_SYMBOL
)

integration <- integration[
  !is.na(integration$FINAL_SYMBOL) &
    integration$FINAL_SYMBOL != "",
]

integration <- integration[
  !duplicated(integration$FINAL_SYMBOL),
]

cat(
  "Unique final targets:",
  nrow(integration),
  "\n"
)

# ------------------------------------------------------------
# 7. EXPECTED FINAL TARGETS
# ------------------------------------------------------------

expected_targets <- c(
  "TOP2A",
  "SGO1",
  "PLK1",
  "AURKA",
  "CDCA5",
  "CDK1",
  "CHEK1",
  "SKA1",
  "GMNN",
  "BRCA1",
  "DNA2",
  "PRKDC"
)

cat("\n")
cat("Expected final target set:\n")
cat(paste(expected_targets, collapse = ", "), "\n")

missing_expected <- setdiff(
  expected_targets,
  integration$FINAL_SYMBOL
)

extra_targets <- setdiff(
  integration$FINAL_SYMBOL,
  expected_targets
)

if (length(missing_expected) > 0) {
  
  stop(
    paste0(
      "Expected targets missing from Script 15 output: ",
      paste(missing_expected, collapse = ", ")
    )
  )
}

if (length(extra_targets) > 0) {
  
  cat(
    "Additional targets detected:",
    paste(extra_targets, collapse = ", "),
    "\n"
  )
}

# ------------------------------------------------------------
# 8. RESTRICT TO FINAL 12
# ------------------------------------------------------------

final_targets <- integration[
  integration$FINAL_SYMBOL %in% expected_targets,
]

# Preserve the 12-target set
final_targets <- final_targets[
  match(
    expected_targets,
    final_targets$FINAL_SYMBOL
  ),
]

rownames(final_targets) <- NULL

cat(
  "Final targets retained:",
  nrow(final_targets),
  "\n"
)

if (nrow(final_targets) != 12) {
  stop("Final target count is not 12.")
}

# ------------------------------------------------------------
# 9. HELPER FOR SAFE NUMERIC EXTRACTION
# ------------------------------------------------------------

get_numeric <- function(df, column, default = 0) {
  
  if (column %in% names(df)) {
    
    x <- suppressWarnings(
      as.numeric(df[[column]])
    )
    
    x[is.na(x)] <- default
    
    return(x)
    
  } else {
    
    return(
      rep(default, nrow(df))
    )
  }
}

get_logical <- function(df, column, default = FALSE) {
  
  if (column %in% names(df)) {
    
    x <- as.logical(df[[column]])
    
    x[is.na(x)] <- default
    
    return(x)
    
  } else {
    
    return(
      rep(default, nrow(df))
    )
  }
}

# ------------------------------------------------------------
# 10. EXTRACT CORE EVIDENCE VARIABLES
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("EXTRACTING CORE EVIDENCE DOMAINS\n")
cat("============================================================\n")

final_targets$TCGA_log2FC <- get_numeric(
  final_targets,
  "TCGA_log2FC"
)

final_targets$TCGA_FDR <- get_numeric(
  final_targets,
  "TCGA_FDR",
  default = 1
)

final_targets$GSE32863_log2FC <- get_numeric(
  final_targets,
  "GSE32863_log2FC"
)

final_targets$GSE32863_FDR <- get_numeric(
  final_targets,
  "GSE32863_FDR",
  default = 1
)

final_targets$degree <- get_numeric(
  final_targets,
  "degree"
)

final_targets$n_pathways <- get_numeric(
  final_targets,
  "n_pathways"
)

final_targets$n_models <- get_numeric(
  final_targets,
  "n_models"
)

final_targets$n_models_non_luad <- get_numeric(
  final_targets,
  "n_models_non_luad"
)

final_targets$LUAD_selectivity_difference <- get_numeric(
  final_targets,
  "LUAD_selectivity_difference"
)

final_targets$FDR <- get_numeric(
  final_targets,
  "FDR",
  default = 1
)

final_targets$n_hallmark_pathways <- get_numeric(
  final_targets,
  "n_hallmark_pathways"
)

final_targets$n_significant_hallmark_support <- get_numeric(
  final_targets,
  "n_significant_hallmark_support"
)

final_targets$core_evidence_domains <- get_numeric(
  final_targets,
  "core_evidence_domains"
)

final_targets$integrated_biological_score <- get_numeric(
  final_targets,
  "integrated_biological_score"
)

final_targets$independent_validation <- get_logical(
  final_targets,
  "independent_validation"
)

final_targets$strong_LUAD_dependency <- get_logical(
  final_targets,
  "strong_LUAD_dependency"
)

final_targets$LUAD_selective_FDR05 <- get_logical(
  final_targets,
  "LUAD_selective_FDR05"
)

final_targets$transcriptomic_evidence <- get_logical(
  final_targets,
  "transcriptomic_evidence"
)

# ------------------------------------------------------------
# 11. EVIDENCE FLAGS
# ------------------------------------------------------------

cat("\nCalculating evidence flags...\n")

final_targets$strong_TCGA_DE <- (
  abs(final_targets$TCGA_log2FC) >= 2 &
    final_targets$TCGA_FDR < 0.05
)

final_targets$very_strong_TCGA_DE <- (
  abs(final_targets$TCGA_log2FC) >= 3 &
    final_targets$TCGA_FDR < 0.05
)

final_targets$strong_independent_validation <- (
  final_targets$independent_validation &
    final_targets$GSE32863_FDR < 0.05
)

final_targets$network_supported <- (
  final_targets$degree >= 30
)

final_targets$hallmark_supported <- (
  final_targets$n_hallmark_pathways > 0
)

final_targets$significant_hallmark_supported <- (
  final_targets$n_significant_hallmark_support > 0
)

final_targets$depmap_supported <- (
  final_targets$strong_LUAD_dependency
)

final_targets$selectivity_supported <- (
  final_targets$LUAD_selective_FDR05
)

# ------------------------------------------------------------
# 12. EVIDENCE DOMAIN COUNT
# ------------------------------------------------------------

final_targets$evidence_domain_count <- rowSums(
  cbind(
    final_targets$strong_TCGA_DE,
    final_targets$strong_independent_validation,
    final_targets$depmap_supported,
    final_targets$selectivity_supported,
    final_targets$network_supported,
    final_targets$hallmark_supported
  )
)

# ------------------------------------------------------------
# 13. FINAL SELECTION SCORE
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("CALCULATING FINAL SELECTION SCORE\n")
cat("============================================================\n")

# Start with the integrated biological score from Script 15.
final_targets$final_selection_score <- (
  final_targets$integrated_biological_score
)

# Add only modest bonuses for independent convergence.
# This does NOT replace the Script 15 evidence score.

final_targets$final_selection_score <- (
  final_targets$final_selection_score +
    
    ifelse(
      final_targets$strong_TCGA_DE,
      0.50,
      0
    ) +
    
    ifelse(
      final_targets$strong_independent_validation,
      0.75,
      0
    ) +
    
    ifelse(
      final_targets$depmap_supported,
      0.75,
      0
    ) +
    
    ifelse(
      final_targets$selectivity_supported,
      0.75,
      0
    ) +
    
    ifelse(
      final_targets$significant_hallmark_supported,
      0.50,
      0
    )
)

# ------------------------------------------------------------
# 14. FINAL RANKING
# ------------------------------------------------------------

final_targets <- final_targets[
  order(
    -final_targets$final_selection_score,
    -final_targets$core_evidence_domains,
    -final_targets$integrated_biological_score,
    -abs(final_targets$TCGA_log2FC)
  ),
]

rownames(final_targets) <- NULL

final_targets$final_rank <- seq_len(
  nrow(final_targets)
)

# ------------------------------------------------------------
# 15. FINAL PRIORITY CATEGORIES
# ------------------------------------------------------------

final_targets$final_priority <- ifelse(
  final_targets$final_rank <= 5,
  "Primary candidate",
  ifelse(
    final_targets$final_rank <= 8,
    "Strong secondary candidate",
    "Supporting candidate"
  )
)

# ------------------------------------------------------------
# 16. FINAL SELECTION STATUS
# ------------------------------------------------------------

final_targets$selection_status <- ifelse(
  final_targets$final_rank <= 5,
  "Selected for primary biological interpretation",
  ifelse(
    final_targets$final_rank <= 8,
    "Retained as strong supporting candidate",
    "Retained as secondary candidate"
  )
)

# ------------------------------------------------------------
# 17. PRINT FINAL RANKING
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("FINAL TARGET RANKING\n")
cat("============================================================\n")

display_columns <- c(
  "final_rank",
  "FINAL_SYMBOL",
  "final_selection_score",
  "TCGA_log2FC",
  "TCGA_FDR",
  "GSE32863_log2FC",
  "GSE32863_FDR",
  "strong_LUAD_dependency",
  "LUAD_selective_FDR05",
  "degree",
  "n_hallmark_pathways",
  "core_evidence_domains",
  "final_priority"
)

display_columns <- display_columns[
  display_columns %in% names(final_targets)
]

print(
  as.data.frame(
    final_targets[, display_columns, drop = FALSE]
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 18. PRIMARY TARGETS
# ------------------------------------------------------------

primary_targets <- final_targets[
  final_targets$final_rank <= 5,
]

cat("\n")
cat("============================================================\n")
cat("PRIMARY TARGETS\n")
cat("============================================================\n")

cat(
  paste(
    primary_targets$FINAL_SYMBOL,
    collapse = ", "
  ),
  "\n"
)

# ------------------------------------------------------------
# 19. STRONG SECONDARY TARGETS
# ------------------------------------------------------------

strong_secondary <- final_targets[
  final_targets$final_rank >= 6 &
    final_targets$final_rank <= 8,
]

cat("\n")
cat("============================================================\n")
cat("STRONG SECONDARY TARGETS\n")
cat("============================================================\n")

cat(
  paste(
    strong_secondary$FINAL_SYMBOL,
    collapse = ", "
  ),
  "\n"
)

# ------------------------------------------------------------
# 20. SUPPORTING TARGETS
# ------------------------------------------------------------

supporting_targets <- final_targets[
  final_targets$final_rank >= 9,
]

cat("\n")
cat("============================================================\n")
cat("SUPPORTING TARGETS\n")
cat("============================================================\n")

cat(
  paste(
    supporting_targets$FINAL_SYMBOL,
    collapse = ", "
  ),
  "\n"
)

# ------------------------------------------------------------
# 21. EVIDENCE SUMMARY
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("FINAL EVIDENCE SUMMARY\n")
cat("============================================================\n")

cat(
  "Total final targets:",
  nrow(final_targets),
  "\n"
)

cat(
  "Strong TCGA DE:",
  sum(final_targets$strong_TCGA_DE),
  "\n"
)

cat(
  "Very strong TCGA DE:",
  sum(final_targets$very_strong_TCGA_DE),
  "\n"
)

cat(
  "Independent validation:",
  sum(final_targets$strong_independent_validation),
  "\n"
)

cat(
  "Strong LUAD dependencies:",
  sum(final_targets$depmap_supported),
  "\n"
)

cat(
  "FDR-supported LUAD selectivity:",
  sum(final_targets$selectivity_supported),
  "\n"
)

cat(
  "Network-supported:",
  sum(final_targets$network_supported),
  "\n"
)

cat(
  "Hallmark-supported:",
  sum(final_targets$hallmark_supported),
  "\n"
)

cat(
  "Significant Hallmark-supported:",
  sum(final_targets$significant_hallmark_supported),
  "\n"
)

# ------------------------------------------------------------
# 22. CREATE PUBLICATION-STYLE FINAL TABLE
# ------------------------------------------------------------

final_priority_table <- final_targets[
  ,
  c(
    "final_rank",
    "FINAL_SYMBOL",
    "final_priority",
    "selection_status",
    "final_selection_score",
    "integrated_biological_score",
    "core_evidence_domains",
    "TCGA_log2FC",
    "TCGA_FDR",
    "GSE32863_log2FC",
    "GSE32863_FDR",
    "strong_LUAD_dependency",
    "LUAD_selective_FDR05",
    "LUAD_selectivity_difference",
    "degree",
    "n_pathways",
    "n_hallmark_pathways",
    "n_significant_hallmark_support",
    "n_models",
    "n_models_non_luad"
  )
]

# ------------------------------------------------------------
# 23. FINAL SHORTLIST TABLE
# ------------------------------------------------------------

primary_target_table <- final_priority_table[
  final_priority_table$final_rank <= 5,
]

# ------------------------------------------------------------
# 24. SAVE RDS OUTPUTS
# ------------------------------------------------------------

cat("\n")
cat("Saving final ranking outputs...\n")

saveRDS(
  final_targets,
  file.path(
    processed_dir,
    "Script16_final_target_ranking.rds"
  )
)

saveRDS(
  final_priority_table,
  file.path(
    processed_dir,
    "Script16_final_priority_table.rds"
  )
)

saveRDS(
  primary_target_table,
  file.path(
    processed_dir,
    "Script16_primary_target_shortlist.rds"
  )
)

# ------------------------------------------------------------
# 25. SAVE CSV OUTPUTS
# ------------------------------------------------------------

write.csv(
  final_priority_table,
  file.path(
    processed_dir,
    "Script16_final_priority_table.csv"
  ),
  row.names = FALSE
)

write.csv(
  primary_target_table,
  file.path(
    processed_dir,
    "Script16_primary_target_shortlist.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 26. FINAL SCORE FIGURE
# ------------------------------------------------------------

cat("\nCreating final target ranking figure...\n")

plot_data <- final_targets

plot_data$FINAL_SYMBOL <- factor(
  plot_data$FINAL_SYMBOL,
  levels = rev(
    plot_data$FINAL_SYMBOL
  )
)

p1 <- ggplot(
  plot_data,
  aes(
    x = final_selection_score,
    y = FINAL_SYMBOL
  )
) +
  geom_col() +
  labs(
    title = "Final Target Selection Ranking",
    x = "Integrated final selection score",
    y = "Target gene"
  ) +
  theme_minimal(base_size = 12)

ggsave(
  file.path(
    figures_dir,
    "Script16_final_target_ranking.png"
  ),
  p1,
  width = 9,
  height = 7,
  dpi = 300
)

ggsave(
  file.path(
    figures_dir,
    "Script16_final_target_ranking.pdf"
  ),
  p1,
  width = 9,
  height = 7
)

# ------------------------------------------------------------
# 27. EVIDENCE DOMAIN FIGURE
# ------------------------------------------------------------

cat("Creating evidence-domain figure...\n")

domain_plot <- data.frame(
  Evidence_Domain = c(
    "Strong TCGA DE",
    "Independent validation",
    "Strong LUAD dependency",
    "LUAD selectivity",
    "Network support",
    "Hallmark support"
  ),
  Number_of_Targets = c(
    sum(final_targets$strong_TCGA_DE),
    sum(final_targets$strong_independent_validation),
    sum(final_targets$depmap_supported),
    sum(final_targets$selectivity_supported),
    sum(final_targets$network_supported),
    sum(final_targets$hallmark_supported)
  )
)

domain_plot$Evidence_Domain <- factor(
  domain_plot$Evidence_Domain,
  levels = rev(domain_plot$Evidence_Domain)
)

p2 <- ggplot(
  domain_plot,
  aes(
    x = Number_of_Targets,
    y = Evidence_Domain
  )
) +
  geom_col() +
  labs(
    title = "Evidence Convergence Across Final Targets",
    x = "Number of final targets",
    y = "Evidence domain"
  ) +
  theme_minimal(base_size = 12)

ggsave(
  file.path(
    figures_dir,
    "Script16_evidence_domain_summary.png"
  ),
  p2,
  width = 9,
  height = 6,
  dpi = 300
)

ggsave(
  file.path(
    figures_dir,
    "Script16_evidence_domain_summary.pdf"
  ),
  p2,
  width = 9,
  height = 6
)

# ------------------------------------------------------------
# 28. FINAL QUALITY CHECKS
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("FINAL QUALITY CHECKS\n")
cat("============================================================\n")

cat(
  "Final target count:",
  nrow(final_targets),
  "\n"
)

cat(
  "Unique target count:",
  length(unique(final_targets$FINAL_SYMBOL)),
  "\n"
)

cat(
  "Missing final scores:",
  sum(
    is.na(
      final_targets$final_selection_score
    )
  ),
  "\n"
)

cat(
  "Missing target symbols:",
  sum(
    is.na(
      final_targets$FINAL_SYMBOL
    ) |
      final_targets$FINAL_SYMBOL == ""
  ),
  "\n"
)

if (
  nrow(final_targets) == 12 &&
  length(unique(final_targets$FINAL_SYMBOL)) == 12 &&
  all(!is.na(final_targets$final_selection_score)) &&
  all(
    !is.na(final_targets$FINAL_SYMBOL) &
    final_targets$FINAL_SYMBOL != ""
  )
) {
  
  cat("Target count: PASSED\n")
  cat("Target uniqueness: PASSED\n")
  cat("Score completeness: PASSED\n")
  cat("Symbol completeness: PASSED\n")
  
} else {
  
  stop(
    "One or more final quality checks failed."
  )
}

# ------------------------------------------------------------
# 29. SESSION INFORMATION
# ------------------------------------------------------------

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    processed_dir,
    "Script16_sessionInfo.txt"
  )
)

# ------------------------------------------------------------
# 30. FINAL OUTPUT
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("SCRIPT 16 COMPLETED SUCCESSFULLY\n")
cat("============================================================\n")

cat("\nFinal primary targets:\n")
cat(
  paste(
    primary_targets$FINAL_SYMBOL,
    collapse = ", "
  ),
  "\n"
)

cat("\nOutputs saved:\n")
cat("- Script16_final_target_ranking.rds\n")
cat("- Script16_final_priority_table.rds\n")
cat("- Script16_final_priority_table.csv\n")
cat("- Script16_primary_target_shortlist.rds\n")
cat("- Script16_primary_target_shortlist.csv\n")
cat("- figures/Script16_final_target_ranking.png\n")
cat("- figures/Script16_final_target_ranking.pdf\n")
cat("- figures/Script16_evidence_domain_summary.png\n")
cat("- figures/Script16_evidence_domain_summary.pdf\n")
cat("- Script16_sessionInfo.txt\n")

cat("\n")
cat("============================================================\n")
cat("PROJECT TARGET-SELECTION PIPELINE COMPLETED\n")
cat("============================================================\n")