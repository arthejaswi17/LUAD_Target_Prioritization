# SCRIPT 15 FINAL BIOLOGICAL EVIDENCE INTEGRATION

# Inputs:
#   Script 12:
#     final_target_evidence_integrated.rds
#     final_target_priority_table.rds

#   Script 14:
#     final_target_functional_summary.rds
#     final_target_Hallmark_enrichment.rds
#     final_target_Hallmark_significant.rds
#     final_target_Hallmark_matrix.rds
#     final_target_module_summary.rds

# Outputs:
#   final_biological_evidence_integration.rds
#   final_biological_evidence_integration.csv
#   final_target_evidence_summary.csv
#   final_target_evidence_summary.rds
#   final_target_evidence_heatmap.png
#   final_target_evidence_heatmap.pdf
#   Script15_sessionInfo.txt


cat("\n")
cat("SCRIPT 15: FINAL BIOLOGICAL EVIDENCE INTEGRATION\n")


# SECTION 1: DIRECTORY SETUP

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

dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

cat("Project directory:", project_dir, "\n")
cat("Processed directory:", processed_dir, "\n")
cat("Figures directory:", figures_dir, "\n\n")


# SECTION 2: REQUIRED PACKAGES

required_packages <- c(
  "ggplot2",
  "pheatmap"
)

for (pkg in required_packages) {
  
  if (!requireNamespace(pkg, quietly = TRUE)) {
    
    stop(
      paste0(
        "Required package '",
        pkg,
        "' is not installed.\n",
        "Install it using install.packages('",
        pkg,
        "')"
      )
    )
    
  }
}

cat("Required packages detected.\n\n")


# SECTION 3: HELPER FUNCTIONS

safe_read_rds <- function(path, description) {
  
  if (!file.exists(path)) {
    
    stop(
      paste0(
        "\nMissing required file: ",
        path,
        "\nExpected input: ",
        description,
        "\n"
      )
    )
    
  }
  
  readRDS(path)
}


find_column <- function(df, candidates) {
  
  available <- colnames(df)
  
  exact_match <- candidates[
    candidates %in% available
  ]
  
  if (length(exact_match) > 0) {
    return(exact_match[1])
  }
  
  lower_available <- tolower(available)
  lower_candidates <- tolower(candidates)
  
  idx <- match(
    lower_candidates,
    lower_available
  )
  
  idx <- idx[!is.na(idx)]
  
  if (length(idx) > 0) {
    return(available[idx[1]])
  }
  
  return(NA_character_)
}


numeric_safe <- function(x) {
  
  if (is.numeric(x)) {
    return(x)
  }
  
  suppressWarnings(
    as.numeric(as.character(x))
  )
}


logical_safe <- function(x) {
  
  if (is.logical(x)) {
    return(x)
  }
  
  y <- tolower(trimws(as.character(x)))
  
  y[y %in% c("true", "t", "yes", "y", "1")] <- TRUE
  y[y %in% c("false", "f", "no", "n", "0")] <- FALSE
  
  result <- suppressWarnings(
    as.logical(y)
  )
  
  result
}


# SECTION 4: LOAD SCRIPT 12 FINAL TARGET TABLE

cat("Loading Script 12 final target table...\n")

priority_path <- file.path(
  processed_dir,
  "final_target_priority_table.rds"
)

priority <- safe_read_rds(
  priority_path,
  "Script 12 final target priority table"
)

priority <- as.data.frame(priority)

cat("Rows loaded:", nrow(priority), "\n")
cat("Columns loaded:", ncol(priority), "\n\n")


# SECTION 5: VALIDATE FINAL TARGET TABLE

symbol_col <- find_column(
  priority,
  c("SYMBOL", "gene_symbol", "symbol", "gene")
)

if (is.na(symbol_col)) {
  
  stop(
    "Could not identify the gene symbol column in the Script 12 table."
  )
  
}

priority$SYMBOL <- as.character(
  priority[[symbol_col]]
)

priority$SYMBOL <- trimws(
  priority$SYMBOL
)

priority <- priority[
  !is.na(priority$SYMBOL) &
    priority$SYMBOL != "",
  ,
  drop = FALSE
]

priority <- priority[
  !duplicated(priority$SYMBOL),
  ,
  drop = FALSE
]

cat("Final targets identified:", nrow(priority), "\n\n")

cat("Final target genes:\n")
cat(
  paste(priority$SYMBOL, collapse = ", "),
  "\n\n"
)


# SECTION 6: LOAD SCRIPT 14 FUNCTIONAL SUMMARY

cat("Loading Script 14 functional summary...\n")

functional_path <- file.path(
  processed_dir,
  "final_target_functional_summary.rds"
)

functional <- safe_read_rds(
  functional_path,
  "Script 14 functional summary"
)

functional <- as.data.frame(functional)

cat(
  "Functional summary rows:",
  nrow(functional),
  "\n"
)

cat(
  "Functional summary columns:",
  ncol(functional),
  "\n\n"
)


# SECTION 7: MERGE FUNCTIONAL INFORMATION

functional_symbol_col <- find_column(
  functional,
  c("SYMBOL", "gene_symbol", "symbol", "gene")
)

if (!is.na(functional_symbol_col)) {
  
  functional$SYMBOL <- as.character(
    functional[[functional_symbol_col]]
  )
  
  functional$SYMBOL <- trimws(
    functional$SYMBOL
  )
  
  functional <- functional[
    !duplicated(functional$SYMBOL),
    ,
    drop = FALSE
  ]
  
}


# SECTION 8: LOAD HALLMARK ENRICHMENT

cat("Loading Script 14 Hallmark enrichment...\n")

hallmark_path <- file.path(
  processed_dir,
  "final_target_Hallmark_enrichment.rds"
)

hallmark <- safe_read_rds(
  hallmark_path,
  "Script 14 Hallmark enrichment"
)

hallmark <- as.data.frame(hallmark)

cat(
  "Hallmark enrichment rows:",
  nrow(hallmark),
  "\n\n"
)


# SECTION 9: LOAD SIGNIFICANT HALLMARK PATHWAYS

cat("Loading significant Hallmark pathways...\n")

hallmark_sig_path <- file.path(
  processed_dir,
  "final_target_Hallmark_significant.rds"
)

hallmark_sig <- safe_read_rds(
  hallmark_sig_path,
  "Script 14 significant Hallmark pathways"
)

hallmark_sig <- as.data.frame(
  hallmark_sig
)

cat(
  "Significant Hallmark pathways:",
  nrow(hallmark_sig),
  "\n\n"
)


# SECTION 10: LOAD HALLMARK MATRIX

cat("Loading target × Hallmark matrix...\n")

matrix_path <- file.path(
  processed_dir,
  "final_target_Hallmark_matrix.rds"
)

hallmark_matrix <- safe_read_rds(
  matrix_path,
  "Script 14 target × Hallmark matrix"
)

hallmark_matrix <- as.data.frame(
  hallmark_matrix,
  check.names = FALSE
)

cat(
  "Hallmark matrix dimensions:",
  nrow(hallmark_matrix),
  "×",
  ncol(hallmark_matrix),
  "\n\n"
)


# SECTION 11: LOAD MODULE SUMMARY

cat("Loading module summary...\n")

module_path <- file.path(
  processed_dir,
  "final_target_module_summary.rds"
)

module_summary <- safe_read_rds(
  module_path,
  "Script 14 module summary"
)

module_summary <- as.data.frame(
  module_summary
)

cat(
  "Module summary rows:",
  nrow(module_summary),
  "\n\n"
)


# SECTION 12: CREATE EVIDENCE TABLE

cat("CREATING INTEGRATED EVIDENCE TABLE\n")

evidence <- data.frame(
  SYMBOL = priority$SYMBOL,
  stringsAsFactors = FALSE
)


# SECTION 13: TRANSCRIPTOMIC EVIDENCE

logfc_col <- find_column(
  priority,
  c(
    "TCGA_shrunken_log2FC",
    "log2FoldChange",
    "log2FC"
  )
)

fdr_col <- find_column(
  priority,
  c(
    "TCGA_FDR",
    "padj",
    "FDR"
  )
)

gse_logfc_col <- find_column(
  priority,
  c(
    "GSE32863_log2FC",
    "GSE_log2FC"
  )
)

gse_fdr_col <- find_column(
  priority,
  c(
    "GSE32863_FDR",
    "GSE_FDR"
  )
)

validated_col <- find_column(
  priority,
  c(
    "GSE32863_validated",
    "independent_validation"
  )
)


if (!is.na(logfc_col)) {
  
  evidence$TCGA_log2FC <- numeric_safe(
    priority[[logfc_col]]
  )
  
} else {
  
  evidence$TCGA_log2FC <- NA_real_
  
}


if (!is.na(fdr_col)) {
  
  evidence$TCGA_FDR <- numeric_safe(
    priority[[fdr_col]]
  )
  
} else {
  
  evidence$TCGA_FDR <- NA_real_
  
}


if (!is.na(gse_logfc_col)) {
  
  evidence$GSE32863_log2FC <- numeric_safe(
    priority[[gse_logfc_col]]
  )
  
} else {
  
  evidence$GSE32863_log2FC <- NA_real_
  
}


if (!is.na(gse_fdr_col)) {
  
  evidence$GSE32863_FDR <- numeric_safe(
    priority[[gse_fdr_col]]
  )
  
} else {
  
  evidence$GSE32863_FDR <- NA_real_
  
}


if (!is.na(validated_col)) {
  
  evidence$GSE32863_validated <- logical_safe(
    priority[[validated_col]]
  )
  
} else {
  
  evidence$GSE32863_validated <- FALSE
  
}


# SECTION 14: DE DIRECTION

evidence$TCGA_direction <- ifelse(
  is.na(evidence$TCGA_log2FC),
  "Unknown",
  ifelse(
    evidence$TCGA_log2FC > 0,
    "Upregulated",
    ifelse(
      evidence$TCGA_log2FC < 0,
      "Downregulated",
      "No change"
    )
  )
)


# SECTION 15: DE STRENGTH

evidence$strong_DE <- (
  !is.na(evidence$TCGA_log2FC) &
    abs(evidence$TCGA_log2FC) >= 2 &
    !is.na(evidence$TCGA_FDR) &
    evidence$TCGA_FDR < 0.05
)


evidence$very_strong_DE <- (
  !is.na(evidence$TCGA_log2FC) &
    abs(evidence$TCGA_log2FC) >= 3 &
    !is.na(evidence$TCGA_FDR) &
    evidence$TCGA_FDR < 0.05
)


# SECTION 16: INDEPENDENT VALIDATION

evidence$independent_validation <- (
  evidence$GSE32863_validated &
    !is.na(evidence$GSE32863_FDR) &
    evidence$GSE32863_FDR < 0.05
)


# SECTION 17: DEPMAP EVIDENCE

cat("Extracting DepMap evidence...\n")

depmap_cols <- c(
  "n_models",
  "n_models_non_luad",
  "median_gene_effect",
  "median_gene_effect_non_luad",
  "LUAD_selectivity_difference",
  "p_value",
  "FDR",
  "strong_LUAD_dependency",
  "LUAD_selective_FDR05"
)

for (col_name in depmap_cols) {
  
  if (col_name %in% colnames(priority)) {
    
    evidence[[col_name]] <- priority[[col_name]]
    
  } else {
    
    evidence[[col_name]] <- NA
    
  }
  
}


if ("strong_LUAD_dependency" %in% colnames(evidence)) {
  
  evidence$strong_LUAD_dependency <- logical_safe(
    evidence$strong_LUAD_dependency
  )
  
} else {
  
  evidence$strong_LUAD_dependency <- FALSE
  
}


if ("LUAD_selective_FDR05" %in% colnames(evidence)) {
  
  evidence$LUAD_selective_FDR05 <- logical_safe(
    evidence$LUAD_selective_FDR05
  )
  
} else {
  
  evidence$LUAD_selective_FDR05 <- FALSE
  
}


# SECTION 18: NETWORK EVIDENCE

cat("Extracting network evidence...\n")

network_cols <- c(
  "degree",
  "module",
  "n_pathways",
  "core_evidence_count",
  "priority_tier"
)

for (col_name in network_cols) {
  
  if (col_name %in% colnames(priority)) {
    
    evidence[[col_name]] <- priority[[col_name]]
    
  } else {
    
    evidence[[col_name]] <- NA
    
  }
  
}


evidence$degree <- numeric_safe(
  evidence$degree
)

evidence$n_pathways <- numeric_safe(
  evidence$n_pathways
)

evidence$core_evidence_count <- numeric_safe(
  evidence$core_evidence_count
)


# SECTION 19: HALLMARK PATHWAY COVERAGE

cat("Calculating Hallmark pathway coverage per target...\n")


hallmark_gene_col <- find_column(
  hallmark_matrix,
  c(
    "SYMBOL",
    "gene_symbol",
    "gene"
  )
)


if (!is.na(hallmark_gene_col)) {
  
  hallmark_matrix$SYMBOL <- as.character(
    hallmark_matrix[[hallmark_gene_col]]
  )
  
  hallmark_matrix$SYMBOL <- trimws(
    hallmark_matrix$SYMBOL
  )
  
}


hallmark_pathway_columns <- setdiff(
  colnames(hallmark_matrix),
  c(
    hallmark_gene_col,
    "SYMBOL",
    "gene_symbol",
    "gene"
  )
)


hallmark_counts <- data.frame(
  SYMBOL = evidence$SYMBOL,
  n_hallmark_pathways = 0,
  stringsAsFactors = FALSE
)


hallmark_names <- character(
  length(evidence$SYMBOL)
)


for (i in seq_len(nrow(evidence))) {
  
  gene <- evidence$SYMBOL[i]
  
  if (
    !is.na(hallmark_gene_col) &&
    gene %in% hallmark_matrix$SYMBOL
  ) {
    
    row_idx <- which(
      hallmark_matrix$SYMBOL == gene
    )[1]
    
    values <- hallmark_matrix[
      row_idx,
      hallmark_pathway_columns,
      drop = TRUE
    ]
    
    values_num <- suppressWarnings(
      as.numeric(as.character(values))
    )
    
    present <- sum(
      !is.na(values_num) &
        values_num > 0
    )
    
    hallmark_counts$n_hallmark_pathways[i] <- present
    
    if (present > 0) {
      
      pathway_names <- hallmark_pathway_columns[
        !is.na(values_num) &
          values_num > 0
      ]
      
      hallmark_names[i] <- paste(
        pathway_names,
        collapse = "; "
      )
      
    } else {
      
      hallmark_names[i] <- ""
      
    }
    
  } else {
    
    hallmark_counts$n_hallmark_pathways[i] <- 0
    hallmark_names[i] <- ""
    
  }
  
}


evidence$n_hallmark_pathways <-
  hallmark_counts$n_hallmark_pathways

evidence$hallmark_pathways <-
  hallmark_names


# SECTION 20: SIGNIFICANT HALLMARK SUPPORT

sig_pathway_col <- find_column(
  hallmark_sig,
  c(
    "gs_name",
    "pathway"
  )
)

sig_fdr_col <- find_column(
  hallmark_sig,
  c(
    "FDR",
    "padj"
  )
)

sig_target_col <- find_column(
  hallmark_sig,
  c(
    "target_genes",
    "genes"
  )
)


significant_hallmark_pathways <- character(0)

if (!is.na(sig_pathway_col)) {
  
  if (!is.na(sig_fdr_col)) {
    
    sig_fdr_values <- numeric_safe(
      hallmark_sig[[sig_fdr_col]]
    )
    
    significant_hallmark_pathways <-
      as.character(
        hallmark_sig[[sig_pathway_col]][
          !is.na(sig_fdr_values) &
            sig_fdr_values < 0.05
        ]
      )
    
  } else {
    
    significant_hallmark_pathways <-
      as.character(
        hallmark_sig[[sig_pathway_col]]
      )
    
  }
  
}


evidence$n_significant_hallmark_support <- 0

evidence$significant_hallmark_support <-
  ""


if (
  length(significant_hallmark_pathways) > 0 &&
  !is.na(hallmark_gene_col)
) {
  
  for (i in seq_len(nrow(evidence))) {
    
    gene <- evidence$SYMBOL[i]
    
    if (gene %in% hallmark_matrix$SYMBOL) {
      
      row_idx <- which(
        hallmark_matrix$SYMBOL == gene
      )[1]
      
      values <- hallmark_matrix[
        row_idx,
        hallmark_pathway_columns,
        drop = TRUE
      ]
      
      values_num <- suppressWarnings(
        as.numeric(as.character(values))
      )
      
      pathway_names <- hallmark_pathway_columns[
        !is.na(values_num) &
          values_num > 0
      ]
      
      supported <- intersect(
        pathway_names,
        significant_hallmark_pathways
      )
      
      evidence$n_significant_hallmark_support[i] <-
        length(supported)
      
      if (length(supported) > 0) {
        
        evidence$significant_hallmark_support[i] <-
          paste(
            supported,
            collapse = "; "
          )
        
      }
      
    }
    
  }
  
}


# SECTION 21: EVIDENCE DOMAIN SCORES

cat("\nCalculating evidence-domain scores...\n")


# Transcriptomic score

evidence$transcriptomic_score <- 0

evidence$transcriptomic_score[
  !is.na(evidence$TCGA_FDR) &
    evidence$TCGA_FDR < 0.05
] <- evidence$transcriptomic_score[
  !is.na(evidence$TCGA_FDR) &
    evidence$TCGA_FDR < 0.05
] + 1

evidence$transcriptomic_score[
  evidence$strong_DE
] <- evidence$transcriptomic_score[
  evidence$strong_DE
] + 1

evidence$transcriptomic_score[
  evidence$very_strong_DE
] <- evidence$transcriptomic_score[
  evidence$very_strong_DE
] + 1


# Independent validation score

evidence$validation_score <- ifelse(
  evidence$independent_validation,
  2,
  0
)


# DepMap score

evidence$depmap_score <- 0

evidence$depmap_score[
  evidence$strong_LUAD_dependency
] <- evidence$depmap_score[
  evidence$strong_LUAD_dependency
] + 1

evidence$depmap_score[
  evidence$LUAD_selective_FDR05
] <- evidence$depmap_score[
  evidence$LUAD_selective_FDR05
] + 2


# Network score

degree_values <- evidence$degree

degree_values[
  is.na(degree_values)
] <- 0


degree_rank <- rank(
  degree_values,
  ties.method = "average"
)

if (max(degree_rank) > 0) {
  
  evidence$network_score <-
    round(
      2 * degree_rank /
        max(degree_rank),
      3
    )
  
} else {
  
  evidence$network_score <- 0
  
}

# Hallmark score

evidence$hallmark_score <- 0

evidence$hallmark_score <-
  pmin(
    evidence$n_hallmark_pathways,
    4
  ) * 0.5

evidence$hallmark_score <-
  evidence$hallmark_score +
  pmin(
    evidence$n_significant_hallmark_support,
    2
  )


# ChEMBL / pharmacological evidence

chembl_candidates <- c(
  "chembl_evidence",
  "ChEMBL_evidence",
  "chembl_score",
  "ChEMBL_score",
  "n_chembl",
  "chembl_hits",
  "drug_evidence"
)

chembl_col <- find_column(
  priority,
  chembl_candidates
)

if (!is.na(chembl_col)) {
  
  evidence$ChEMBL_evidence <-
    priority[[chembl_col]]
  
  chembl_numeric <- numeric_safe(
    priority[[chembl_col]]
  )
  
  chembl_numeric[
    is.na(chembl_numeric)
  ] <- 0
  
  evidence$chembl_score <-
    ifelse(
      chembl_numeric > 0,
      2,
      0
    )
  
} else {
  
  evidence$ChEMBL_evidence <- NA
  evidence$chembl_score <- 0
  
}


# SECTION 22: OVERALL BIOLOGICAL EVIDENCE SCORE

cat("Calculating integrated biological evidence score...\n")


evidence$integrated_biological_score <-
  evidence$transcriptomic_score +
  evidence$validation_score +
  evidence$depmap_score +
  evidence$network_score +
  evidence$hallmark_score +
  evidence$chembl_score


# SECTION 23: CORE EVIDENCE COUNT

evidence$core_evidence_domains <- 0

evidence$core_evidence_domains <-
  evidence$core_evidence_domains +
  as.integer(
    !is.na(evidence$TCGA_FDR) &
      evidence$TCGA_FDR < 0.05
  )

evidence$core_evidence_domains <-
  evidence$core_evidence_domains +
  as.integer(
    evidence$independent_validation
  )

evidence$core_evidence_domains <-
  evidence$core_evidence_domains +
  as.integer(
    evidence$strong_LUAD_dependency
  )

evidence$core_evidence_domains <-
  evidence$core_evidence_domains +
  as.integer(
    evidence$LUAD_selective_FDR05
  )

evidence$core_evidence_domains <-
  evidence$core_evidence_domains +
  as.integer(
    evidence$n_hallmark_pathways > 0
  )

evidence$core_evidence_domains <-
  evidence$core_evidence_domains +
  as.integer(
    evidence$n_significant_hallmark_support > 0
  )

evidence$core_evidence_domains <-
  evidence$core_evidence_domains +
  as.integer(
    evidence$degree > 0
  )


# SECTION 24: BIOLOGICAL EVIDENCE CATEGORY

evidence$biological_evidence_category <- "Limited"


evidence$biological_evidence_category[
  evidence$core_evidence_domains >= 5
] <- "Strong"


evidence$biological_evidence_category[
  evidence$core_evidence_domains >= 6
] <- "Very strong"


# SECTION 25: RANK TARGETS

evidence <- evidence[
  order(
    -evidence$integrated_biological_score,
    -evidence$core_evidence_domains,
    -abs(evidence$TCGA_log2FC)
  ),
  ,
  drop = FALSE
]

evidence$final_biological_rank <-
  seq_len(nrow(evidence))


# SECTION 26: FINAL PRIORITY INTERPRETATION

evidence$interpretation <- "Secondary candidate"


evidence$interpretation[
  evidence$biological_evidence_category ==
    "Very strong"
] <- "High-confidence biological candidate"


evidence$interpretation[
  evidence$biological_evidence_category ==
    "Strong"
] <- "Strong biological candidate"


# SECTION 27: SUMMARY

cat("\n")
cat("FINAL BIOLOGICAL EVIDENCE SUMMARY\n")
cat(
  "Final targets evaluated:",
  nrow(evidence),
  "\n"
)

cat(
  "Strong DE candidates:",
  sum(evidence$strong_DE, na.rm = TRUE),
  "\n"
)

cat(
  "Very strong DE candidates:",
  sum(evidence$very_strong_DE, na.rm = TRUE),
  "\n"
)

cat(
  "Independently validated:",
  sum(evidence$independent_validation, na.rm = TRUE),
  "\n"
)

cat(
  "Strong LUAD dependencies:",
  sum(evidence$strong_LUAD_dependency, na.rm = TRUE),
  "\n"
)

cat(
  "FDR-supported LUAD selectivity:",
  sum(evidence$LUAD_selective_FDR05, na.rm = TRUE),
  "\n"
)

cat(
  "Targets represented in Hallmark:",
  sum(evidence$n_hallmark_pathways > 0, na.rm = TRUE),
  "\n"
)

cat(
  "Targets supported by significant Hallmark pathways:",
  sum(
    evidence$n_significant_hallmark_support > 0,
    na.rm = TRUE
  ),
  "\n\n"
)


# SECTION 28: PRINT FINAL TABLE

display_cols <- c(
  "final_biological_rank",
  "SYMBOL",
  "TCGA_log2FC",
  "TCGA_FDR",
  "GSE32863_log2FC",
  "GSE32863_FDR",
  "independent_validation",
  "strong_LUAD_dependency",
  "LUAD_selective_FDR05",
  "degree",
  "n_hallmark_pathways",
  "n_significant_hallmark_support",
  "core_evidence_domains",
  "integrated_biological_score",
  "biological_evidence_category",
  "interpretation"
)

display_cols <- display_cols[
  display_cols %in%
    colnames(evidence)
]


cat("FINAL BIOLOGICAL RANKING:\n\n")

print(
  evidence[
    ,
    display_cols,
    drop = FALSE
  ],
  row.names = FALSE
)


# SECTION 29: TOP CANDIDATES

cat("\n")
cat("TOP BIOLOGICAL CANDIDATES\n")


top_n <- min(
  5,
  nrow(evidence)
)


top_candidates <- evidence[
  seq_len(top_n),
  ,
  drop = FALSE
]


cat(
  paste(
    seq_len(top_n),
    top_candidates$SYMBOL,
    sep = ". ",
    collapse = "\n"
  ),
  "\n\n"
)


# SECTION 30: HEATMAP DATA

cat("Creating evidence heatmap...\n")


heatmap_data <- evidence[
  ,
  c(
    "SYMBOL",
    "transcriptomic_score",
    "validation_score",
    "depmap_score",
    "network_score",
    "hallmark_score",
    "chembl_score"
  ),
  drop = FALSE
]


rownames(heatmap_data) <-
  heatmap_data$SYMBOL


heatmap_data$SYMBOL <- NULL


heatmap_matrix <- as.matrix(
  heatmap_data
)


# SECTION 31: HEATMAP

png_path <- file.path(
  figures_dir,
  "Script15_final_biological_evidence_heatmap.png"
)

pdf_path <- file.path(
  figures_dir,
  "Script15_final_biological_evidence_heatmap.pdf"
)


png(
  png_path,
  width = 1800,
  height = 1200,
  res = 180
)

pheatmap::pheatmap(
  heatmap_matrix,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  border_color = NA,
  main = "Integrated Biological Evidence Across Final LUAD Targets",
  fontsize = 10,
  fontsize_row = 10,
  fontsize_col = 9
)

dev.off()


pdf(
  pdf_path,
  width = 11,
  height = 8
)

pheatmap::pheatmap(
  heatmap_matrix,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  border_color = NA,
  main = "Integrated Biological Evidence Across Final LUAD Targets",
  fontsize = 10,
  fontsize_row = 10,
  fontsize_col = 9
)

dev.off()


# SECTION 32: SAVE OUTPUT TABLES

cat("\nSaving integrated evidence outputs...\n")


rds_path <- file.path(
  processed_dir,
  "final_biological_evidence_integration.rds"
)

csv_path <- file.path(
  processed_dir,
  "final_biological_evidence_integration.csv"
)

summary_rds_path <- file.path(
  processed_dir,
  "final_target_evidence_summary.rds"
)

summary_csv_path <- file.path(
  processed_dir,
  "final_target_evidence_summary.csv"
)


saveRDS(
  evidence,
  rds_path
)


write.csv(
  evidence,
  csv_path,
  row.names = FALSE,
  na = ""
)


summary_table <- evidence[
  ,
  display_cols,
  drop = FALSE
]


saveRDS(
  summary_table,
  summary_rds_path
)


write.csv(
  summary_table,
  summary_csv_path,
  row.names = FALSE,
  na = ""
)


# SECTION 33: SAVE SESSION INFORMATION

session_path <- file.path(
  processed_dir,
  "Script15_sessionInfo.txt"
)


capture.output(
  sessionInfo(),
  file = session_path
)


# SECTION 34: FINAL QUALITY CHECKS

cat("\n")
cat("FINAL QUALITY CHECKS\n")


cat(
  "Final target count:",
  nrow(evidence),
  "\n"
)


cat(
  "Unique target count:",
  length(unique(evidence$SYMBOL)),
  "\n"
)


cat(
  "Missing integrated scores:",
  sum(
    is.na(
      evidence$integrated_biological_score
    )
  ),
  "\n"
)


cat(
  "Missing target symbols:",
  sum(
    is.na(evidence$SYMBOL) |
      evidence$SYMBOL == ""
  ),
  "\n"
)


if (
  nrow(evidence) == length(unique(evidence$SYMBOL))
) {
  
  cat(
    "Target uniqueness: PASSED\n"
  )
  
} else {
  
  cat(
    "Target uniqueness: FAILED\n"
  )
  
}


if (
  all(
    !is.na(
      evidence$integrated_biological_score
    )
  )
) {
  
  cat(
    "Integrated score completeness: PASSED\n"
  )
  
} else {
  
  cat(
    "Integrated score completeness: FAILED\n"
  )
  
}


if (
  file.exists(rds_path) &&
  file.exists(csv_path) &&
  file.exists(summary_rds_path) &&
  file.exists(summary_csv_path)
) {
  
  cat(
    "Core output files: PASSED\n"
  )
  
} else {
  
  cat(
    "Core output files: FAILED\n"
  )
  
}


# SECTION 35: FINAL MESSAGE

cat("\n")
cat("SCRIPT 15 COMPLETED SUCCESSFULLY\n")
cat("Outputs saved:\n")

cat(
  "- final_biological_evidence_integration.rds\n"
)

cat(
  "- final_biological_evidence_integration.csv\n"
)

cat(
  "- final_target_evidence_summary.rds\n"
)

cat(
  "- final_target_evidence_summary.csv\n"
)

cat(
  "- figures/Script15_final_biological_evidence_heatmap.png\n"
)

cat(
  "- figures/Script15_final_biological_evidence_heatmap.pdf\n"
)

cat(
  "- Script15_sessionInfo.txt\n"
)

cat("\n")
cat(
  "Next stage: Final target ranking and biological selection\n"
)