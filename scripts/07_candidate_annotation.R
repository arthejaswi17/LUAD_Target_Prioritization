# ============================================================
# SCRIPT 07: CANDIDATE GENE ANNOTATION
# ============================================================
#
# Purpose:
#   Prepare the candidate genes from Script 06 for downstream
#   biological and external evidence integration.
#
# Inputs:
#   candidate_gene_prioritization.rds
#   top_30_candidate_genes.rds
#
# Outputs:
#   candidate_gene_annotation.rds
#   candidate_gene_annotation.csv
#   top_30_candidate_annotation.rds
#   top_30_candidate_annotation.csv
#   candidate_gene_annotation_sessionInfo.txt
#
# ============================================================


# ============================================================
# 1. SETUP
# ============================================================

cat("\n")
cat("============================================================\n")
cat("SCRIPT 07: CANDIDATE GENE ANNOTATION\n")
cat("============================================================\n\n")

options(stringsAsFactors = FALSE)

# Directories
processed_dir <- "data/processed"
figures_dir <- "figures"

if (!dir.exists(processed_dir)) {
  dir.create(
    processed_dir,
    recursive = TRUE
  )
}

if (!dir.exists(figures_dir)) {
  dir.create(
    figures_dir,
    recursive = TRUE
  )
}


# ============================================================
# 2. LOAD REQUIRED PACKAGES
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})


# ============================================================
# 3. LOAD SCRIPT 06 OUTPUTS
# ============================================================

cat("Loading Script 06 candidate results...\n")

candidate_file <- file.path(
  processed_dir,
  "candidate_gene_prioritization.rds"
)

top_candidate_file <- file.path(
  processed_dir,
  "top_30_candidate_genes.rds"
)

if (!file.exists(candidate_file)) {
  stop(
    "Missing input file: ",
    candidate_file,
    "\nRun Script 06 first."
  )
}

if (!file.exists(top_candidate_file)) {
  stop(
    "Missing input file: ",
    top_candidate_file,
    "\nRun Script 06 first."
  )
}

candidate_table <- readRDS(candidate_file)

top_candidates <- readRDS(top_candidate_file)

cat(
  "Candidate rows loaded:",
  nrow(candidate_table),
  "\n"
)

cat(
  "Top candidate rows loaded:",
  nrow(top_candidates),
  "\n"
)


# ============================================================
# 4. VALIDATE INPUT STRUCTURE
# ============================================================

cat("\nValidating candidate table...\n")

required_columns <- c(
  "gene",
  "baseMean",
  "log2FoldChange",
  "lfcSE",
  "pvalue",
  "padj",
  "n_GSEA_pathways",
  "GSEA_pathways",
  "DE_direction",
  "abs_log2FC",
  "neg_log10_FDR",
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


# ============================================================
# 5. CLEAN ENSEMBL IDS
# ============================================================

cat("\nCleaning Ensembl gene IDs...\n")

candidate_table <- candidate_table %>%
  mutate(
    ensembl_gene = sub(
      "\\..*$",
      "",
      as.character(gene)
    )
  )

top_candidates <- top_candidates %>%
  mutate(
    ensembl_gene = sub(
      "\\..*$",
      "",
      as.character(gene)
    )
  )

# Remove accidental whitespace
candidate_table$ensembl_gene <- trimws(
  candidate_table$ensembl_gene
)

top_candidates$ensembl_gene <- trimws(
  top_candidates$ensembl_gene
)

cat(
  "Unique Ensembl IDs:",
  length(
    unique(candidate_table$ensembl_gene)
  ),
  "\n"
)


# ============================================================
# 6. CHECK ENSEMBL ID FORMAT
# ============================================================

valid_ensembl <- grepl(
  "^ENSG[0-9]+$",
  candidate_table$ensembl_gene
)

cat(
  "Valid Ensembl-format IDs:",
  sum(valid_ensembl),
  "\n"
)

cat(
  "Invalid/non-Ensembl IDs:",
  sum(!valid_ensembl),
  "\n"
)


# ============================================================
# 7. GENE SYMBOL ANNOTATION
# ============================================================
#
# We intentionally do NOT require org.Hs.eg.db here.
#
# The script attempts to use biomaRt if available.
# If biomaRt is unavailable or the connection fails,
# the analysis continues and preserves the Ensembl IDs.
#
# This prevents annotation from becoming a hard dependency.
# ============================================================

cat("\n")
cat("Attempting Ensembl → gene symbol annotation...\n")

annotation_table <- NULL

if (
  requireNamespace(
    "biomaRt",
    quietly = TRUE
  )
) {
  
  cat("biomaRt package detected.\n")
  cat("Connecting to Ensembl...\n")
  
  annotation_table <- tryCatch(
    
    {
      
      mart <- biomaRt::useEnsembl(
        biomart = "genes",
        dataset = "hsapiens_gene_ensembl"
      )
      
      unique_ids <- unique(
        candidate_table$ensembl_gene[
          valid_ensembl
        ]
      )
      
      if (length(unique_ids) == 0) {
        NULL
      } else {
        
        biomaRt::getBM(
          attributes = c(
            "ensembl_gene_id",
            "external_gene_name"
          ),
          filters = "ensembl_gene_id",
          values = unique_ids,
          mart = mart
        )
      }
      
    },
    
    error = function(e) {
      
      cat(
        "biomaRt annotation failed:\n",
        conditionMessage(e),
        "\n"
      )
      
      NULL
    }
  )
  
} else {
  
  cat(
    "biomaRt is not installed.\n"
  )
  
  cat(
    "Continuing without external annotation.\n"
  )
}


# ============================================================
# 8. MERGE GENE SYMBOLS
# ============================================================

if (
  !is.null(annotation_table) &&
  nrow(annotation_table) > 0
) {
  
  annotation_table <- annotation_table %>%
    dplyr::rename(
      ensembl_gene = ensembl_gene_id,
      gene_symbol = external_gene_name
    ) %>%
    dplyr::distinct(
      ensembl_gene,
      .keep_all = TRUE
    )
  
  candidate_table <- candidate_table %>%
    dplyr::left_join(
      annotation_table,
      by = "ensembl_gene"
    )
  
  top_candidates <- top_candidates %>%
    dplyr::left_join(
      annotation_table,
      by = "ensembl_gene"
    )
  
} else {
  
  candidate_table$gene_symbol <- NA_character_
  
  top_candidates$gene_symbol <- NA_character_
}


# ============================================================
# 9. FALLBACK FOR UNANNOTATED GENES
# ============================================================

candidate_table <- candidate_table %>%
  mutate(
    gene_symbol = ifelse(
      is.na(gene_symbol) |
        gene_symbol == "",
      ensembl_gene,
      gene_symbol
    )
  )

top_candidates <- top_candidates %>%
  mutate(
    gene_symbol = ifelse(
      is.na(gene_symbol) |
        gene_symbol == "",
      ensembl_gene,
      gene_symbol
    )
  )


# ============================================================
# 10. REORDER COLUMNS
# ============================================================

candidate_table <- candidate_table %>%
  dplyr::select(
    gene_symbol,
    ensembl_gene,
    everything()
  )

top_candidates <- top_candidates %>%
  dplyr::select(
    gene_symbol,
    ensembl_gene,
    everything()
  )


# ============================================================
# 11. ANNOTATION SUMMARY
# ============================================================

annotated_count <- sum(
  candidate_table$gene_symbol !=
    candidate_table$ensembl_gene
)

total_count <- nrow(candidate_table)

cat("\n")
cat("============================================================\n")
cat("ANNOTATION SUMMARY\n")
cat("============================================================\n")

cat(
  "Total candidate genes:",
  total_count,
  "\n"
)

cat(
  "Genes resolved to symbols:",
  annotated_count,
  "\n"
)

cat(
  "Genes without symbol annotation:",
  total_count - annotated_count,
  "\n"
)

cat(
  "Annotation rate:",
  round(
    100 * annotated_count / total_count,
    2
  ),
  "%\n"
)


# ============================================================
# 12. TOP CANDIDATE SUMMARY
# ============================================================

cat("\n")
cat("Top candidate genes:\n\n")

top_display <- top_candidates %>%
  dplyr::select(
    gene_symbol,
    ensembl_gene,
    log2FoldChange,
    padj,
    DE_direction,
    n_GSEA_pathways,
    GSEA_pathways,
    prioritization_score
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
    prioritization_score = round(
      prioritization_score,
      3
    )
  ) %>%
  as.data.frame()

print(
  top_display,
  row.names = FALSE
)


# ============================================================
# 13. SAVE COMPLETE ANNOTATED TABLE
# ============================================================

saveRDS(
  candidate_table,
  file.path(
    processed_dir,
    "candidate_gene_annotation.rds"
  )
)

readr::write_csv(
  candidate_table,
  file.path(
    processed_dir,
    "candidate_gene_annotation.csv"
  )
)


# ============================================================
# 14. SAVE TOP 30 ANNOTATED TABLE
# ============================================================

saveRDS(
  top_candidates,
  file.path(
    processed_dir,
    "top_30_candidate_annotation.rds"
  )
)

readr::write_csv(
  top_candidates,
  file.path(
    processed_dir,
    "top_30_candidate_annotation.csv"
  )
)


# ============================================================
# 15. BASIC QUALITY CHECKS
# ============================================================

cat("\n")
cat("Running final quality checks...\n")

if (
  nrow(candidate_table) !=
  nrow(readRDS(candidate_file))
) {
  
  stop(
    "Candidate row count changed unexpectedly."
  )
}

if (
  anyDuplicated(
    candidate_table$ensembl_gene
  ) > 0
) {
  
  cat(
    "Note: duplicate Ensembl IDs are present.\n"
  )
  
} else {
  
  cat(
    "Ensembl ID uniqueness: PASSED\n"
  )
}

if (
  any(
    is.na(
      candidate_table$log2FoldChange
    )
  )
) {
  
  cat(
    "Warning: NA log2FC values detected.\n"
  )
  
} else {
  
  cat(
    "log2FC completeness: PASSED\n"
  )
}


# ============================================================
# 16. SESSION INFORMATION
# ============================================================

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    processed_dir,
    "candidate_gene_annotation_sessionInfo.txt"
  )
)


# ============================================================
# 17. FINAL MESSAGE
# ============================================================

cat("\n")
cat("============================================================\n")
cat("SCRIPT 07 COMPLETED SUCCESSFULLY\n")
cat("============================================================\n\n")

cat("Outputs saved:\n")

cat(
  "- candidate_gene_annotation.rds\n"
)

cat(
  "- candidate_gene_annotation.csv\n"
)

cat(
  "- top_30_candidate_annotation.rds\n"
)

cat(
  "- top_30_candidate_annotation.csv\n"
)

cat(
  "- candidate_gene_annotation_sessionInfo.txt\n"
)

cat("\n")
cat("============================================================\n")