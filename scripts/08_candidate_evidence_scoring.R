# SCRIPT 08 CANDIDATE EVIDENCE SCORING

# Integrate the strongest evidence currently available from differential expression and GSEA to create a transparent candidate evidence score.

# Inputs:
#   candidate_gene_annotation.rds
#   top_30_candidate_annotation.rds

# Outputs:
#   candidate_evidence_scoring.rds
#   candidate_evidence_scoring.csv
#   top_30_candidate_evidence.rds
#   top_30_candidate_evidence.csv
#   candidate_evidence_scoring_sessionInfo.txt


# 1. SETUP

cat("\n")
cat("SCRIPT 08: CANDIDATE EVIDENCE SCORING\n")

options(stringsAsFactors = FALSE)

processed_dir <- "data/processed"

if (!dir.exists(processed_dir)) {
  dir.create(
    processed_dir,
    recursive = TRUE
  )
}


# 2. LOAD PACKAGES

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})


# 3. LOAD INPUT FILES

cat("Loading Script 07 outputs...\n")

candidate_file <- file.path(
  processed_dir,
  "candidate_gene_annotation.rds"
)

top_candidate_file <- file.path(
  processed_dir,
  "top_30_candidate_annotation.rds"
)

if (!file.exists(candidate_file)) {
  stop(
    "Missing file: ",
    candidate_file,
    "\nRun Script 07 first."
  )
}

if (!file.exists(top_candidate_file)) {
  stop(
    "Missing file: ",
    top_candidate_file,
    "\nRun Script 07 first."
  )
}

candidate_table <- readRDS(candidate_file)
top_candidates <- readRDS(top_candidate_file)

cat(
  "Candidate genes loaded:",
  nrow(candidate_table),
  "\n"
)

cat(
  "Top candidates loaded:",
  nrow(top_candidates),
  "\n"
)


# 4. VALIDATE REQUIRED COLUMNS

required_columns <- c(
  "gene_symbol",
  "ensembl_gene",
  "baseMean",
  "log2FoldChange",
  "lfcSE",
  "pvalue",
  "padj",
  "DE_direction",
  "n_GSEA_pathways",
  "GSEA_pathways",
  "prioritization_score"
)

missing_columns <- setdiff(
  required_columns,
  colnames(candidate_table)
)

if (length(missing_columns) > 0) {
  
  stop(
    paste(
      "Missing required columns:",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  )
}

cat("Required columns: PASSED\n")


# 5. CREATE BASIC EVIDENCE VARIABLES

cat("\nCreating evidence variables...\n")

candidate_table <- candidate_table %>%
  mutate(
    
    # Absolute effect size
    abs_log2FC = abs(
      as.numeric(log2FoldChange)
    ),
    
    # Statistical significance
    neg_log10_padj = -log10(
      pmax(
        as.numeric(padj),
        .Machine$double.xmin
      )
    ),
    
    # GSEA support
    GSEA_supported = (
      !is.na(n_GSEA_pathways) &
        n_GSEA_pathways > 0
    ),
    
    # Strong GSEA support
    GSEA_strong = (
      !is.na(n_GSEA_pathways) &
        n_GSEA_pathways >= 2
    ),
    
    # Strong expression effect
    strong_DE = (
      !is.na(abs_log2FC) &
        abs_log2FC >= 2
    ),
    
    # Very strong expression effect
    very_strong_DE = (
      !is.na(abs_log2FC) &
        abs_log2FC >= 3
    )
  )


# 6. CREATE NORMALIZED COMPONENT SCORES

# Each component is scaled from 0 to 1.
# DE effect:
#   Larger absolute log2FC = stronger score
# Statistical evidence:
#   Smaller FDR = stronger score
# GSEA:
#   More leading-edge pathway support = stronger score
# This avoids letting one extremely small p-value dominate the entire prioritization.

cat("Calculating normalized evidence scores...\n")


# DE effect score

max_abs_fc <- max(
  candidate_table$abs_log2FC,
  na.rm = TRUE
)

candidate_table$DE_effect_score <- pmin(
  candidate_table$abs_log2FC / max_abs_fc,
  1
)


# Statistical significance score

max_neg_log10_padj <- max(
  candidate_table$neg_log10_padj[
    is.finite(
      candidate_table$neg_log10_padj
    )
  ],
  na.rm = TRUE
)

candidate_table$DE_significance_score <- pmin(
  candidate_table$neg_log10_padj /
    max_neg_log10_padj,
  1
)


# GSEA support score

max_gsea <- max(
  candidate_table$n_GSEA_pathways,
  na.rm = TRUE
)

if (
  is.na(max_gsea) ||
  max_gsea == 0
) {
  
  candidate_table$GSEA_support_score <- 0
  
} else {
  
  candidate_table$GSEA_support_score <- pmin(
    candidate_table$n_GSEA_pathways /
      max_gsea,
    1
  )
}


# 7. COMPOSITE EVIDENCE SCORE

# Weighting:
#   45% DE effect
#   30% statistical significance
#   25% GSEA pathway support

candidate_table <- candidate_table %>%
  mutate(
    
    expression_evidence_score =
      100 * (
        0.45 * DE_effect_score +
          0.30 * DE_significance_score +
          0.25 * GSEA_support_score
      )
    
  )


# 8. EVIDENCE CATEGORY

candidate_table <- candidate_table %>%
  mutate(
    
    evidence_category = case_when(
      
      expression_evidence_score >= 75 &
        GSEA_supported &
        strong_DE ~
        "High expression evidence",
      
      expression_evidence_score >= 50 &
        (
          GSEA_supported |
            strong_DE
        ) ~
        "Moderate expression evidence",
      
      TRUE ~
        "Lower expression evidence"
    )
  )


# 9. FINAL RANK

candidate_table <- candidate_table %>%
  arrange(
    desc(expression_evidence_score),
    desc(abs_log2FC),
    padj
  ) %>%
  mutate(
    expression_evidence_rank =
      row_number()
  )


# 10. SELECT FINAL OUTPUT COLUMNS

candidate_evidence <- candidate_table %>%
  dplyr::select(
    gene_symbol,
    ensembl_gene,
    baseMean,
    log2FoldChange,
    abs_log2FC,
    lfcSE,
    pvalue,
    padj,
    DE_direction,
    n_GSEA_pathways,
    GSEA_pathways,
    GSEA_supported,
    GSEA_strong,
    strong_DE,
    very_strong_DE,
    DE_effect_score,
    DE_significance_score,
    GSEA_support_score,
    expression_evidence_score,
    evidence_category,
    expression_evidence_rank,
    prioritization_score
  )


# 11. TOP 30

top_30_evidence <- candidate_evidence %>%
  slice_head(
    n = 30
  )


# 12. SUMMARY

cat("\n")
cat("CANDIDATE EVIDENCE SCORING SUMMARY\n")

cat(
  "Total candidates:",
  nrow(candidate_evidence),
  "\n"
)

cat(
  "GSEA-supported candidates:",
  sum(
    candidate_evidence$GSEA_supported,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Strong DE candidates (|log2FC| >= 2):",
  sum(
    candidate_evidence$strong_DE,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Very strong DE candidates (|log2FC| >= 3):",
  sum(
    candidate_evidence$very_strong_DE,
    na.rm = TRUE
  ),
  "\n"
)

cat("\nEvidence categories:\n")

print(
  table(
    candidate_evidence$evidence_category
  )
)


# 13. TOP 30 DISPLAY

cat("\n")
cat("Top 30 candidates by expression evidence:\n\n")

top_display <- top_30_evidence %>%
  dplyr::select(
    gene_symbol,
    log2FoldChange,
    padj,
    DE_direction,
    n_GSEA_pathways,
    expression_evidence_score,
    evidence_category
  ) %>%
  mutate(
    log2FoldChange = round(
      log2FoldChange,
      3
    ),
    padj = signif(
      padj,
      3
    ),
    expression_evidence_score = round(
      expression_evidence_score,
      2
    )
  ) %>%
  as.data.frame()

print(
  top_display,
  row.names = FALSE
)


# 14. SAVE OUTPUTS

saveRDS(
  candidate_evidence,
  file.path(
    processed_dir,
    "candidate_evidence_scoring.rds"
  )
)

readr::write_csv(
  candidate_evidence,
  file.path(
    processed_dir,
    "candidate_evidence_scoring.csv"
  )
)


saveRDS(
  top_30_evidence,
  file.path(
    processed_dir,
    "top_30_candidate_evidence.rds"
  )
)

readr::write_csv(
  top_30_evidence,
  file.path(
    processed_dir,
    "top_30_candidate_evidence.csv"
  )
)


# 15. SESSION INFO

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    processed_dir,
    "candidate_evidence_scoring_sessionInfo.txt"
  )
)


# 16. FINAL CHECKS

cat("\n")
cat("Running final checks...\n")

if (
  nrow(candidate_evidence) !=
  nrow(candidate_table)
) {
  
  stop(
    "Output row count does not match input."
  )
  
}

if (
  any(
    !is.finite(
      candidate_evidence$expression_evidence_score
    )
  )
) {
  
  stop(
    "Non-finite expression evidence scores detected."
  )
  
}

if (
  any(
    candidate_evidence$expression_evidence_score < 0 |
    candidate_evidence$expression_evidence_score > 100
  )
) {
  
  stop(
    "Evidence score outside expected 0-100 range."
  )
  
}

cat("Row count check: PASSED\n")
cat("Score validity check: PASSED\n")


# 17. COMPLETION

cat("\n")
cat("SCRIPT 08 COMPLETED SUCCESSFULLY\n")

cat("Outputs saved:\n")
cat("- candidate_evidence_scoring.rds\n")
cat("- candidate_evidence_scoring.csv\n")
cat("- top_30_candidate_evidence.rds\n")
cat("- top_30_candidate_evidence.csv\n")
cat("- candidate_evidence_scoring_sessionInfo.txt\n")

cat("\n")
cat("Next stage: ChEMBL pharmacological evidence (Script 09)\n")

cat("\n")