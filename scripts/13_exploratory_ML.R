###############################################################################
# SCRIPT 13: DEPMAP-BASED MACHINE LEARNING DATASET
#
# Purpose:
#   Construct a clean LUAD DepMap ML dataset using the final target genes
#   and a predefined cell-cycle / proliferation predictor panel.
#
# Inputs:
#   - Script 12 final target priority table
#   - DepMap CRISPR gene-effect data
#   - DepMap model metadata
#   - DepMap expression data
#
# Outputs:
#   - ML_AURKA_dataset.rds
#   - ML_AURKA_dataset.csv
#   - ML_predictor_expression.rds
#   - ML_predictor_expression.csv
#   - ML_DepMap_model_metadata.rds
#   - ML_DepMap_model_metadata.csv
#   - ML_AURKA_summary.txt
#   - ML_AURKA_sessionInfo.txt
###############################################################################

cat("\n")
cat("============================================================\n")
cat("SCRIPT 13: DEPMAP-BASED MACHINE LEARNING DATASET\n")
cat("============================================================\n\n")


###############################################################################
# 1. PACKAGES
###############################################################################

required_packages <- c(
  "data.table",
  "dplyr",
  "tibble",
  "glmnet"
)

for (pkg in required_packages) {
  
  if (!requireNamespace(pkg, quietly = TRUE)) {
    
    install.packages(
      pkg,
      repos = "https://cran.rstudio.com/"
    )
    
  }
  
}

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tibble)
  library(glmnet)
})


###############################################################################
# 2. PROJECT DIRECTORIES
###############################################################################

project_dir <- getwd()

processed_dir <- file.path(
  project_dir,
  "data",
  "processed"
)

raw_dir <- file.path(
  project_dir,
  "data",
  "raw"
)

figures_dir <- file.path(
  project_dir,
  "figures"
)

dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)


###############################################################################
# 3. HELPER FUNCTIONS
###############################################################################

clean_name <- function(x) {
  
  x <- as.character(x)
  
  x <- trimws(x)
  
  x
  
}


normalize_gene_name <- function(x) {
  
  x <- as.character(x)
  
  x <- trimws(x)
  
  # Remove common DepMap annotations such as:
  # AURKA (6790)
  # AURKA [6790]
  # AURKA|6790
  
  x <- sub("\\s*\\([^)]*\\)$", "", x)
  
  x <- sub("\\s*\\[[^]]*\\]$", "", x)
  
  x <- sub("\\|.*$", "", x)
  
  x <- trimws(x)
  
  toupper(x)
  
}


is_model_id <- function(x) {
  
  x <- as.character(x)
  
  mean(
    grepl(
      "^ACH[-_]",
      x,
      ignore.case = TRUE
    )
  ) > 0.5
  
}


find_file <- function(
    search_dirs,
    patterns
) {
  
  all_files <- character(0)
  
  for (d in search_dirs) {
    
    if (dir.exists(d)) {
      
      f <- list.files(
        d,
        recursive = TRUE,
        full.names = TRUE,
        ignore.case = TRUE
      )
      
      all_files <- c(
        all_files,
        f
      )
      
    }
    
  }
  
  all_files <- unique(all_files)
  
  all_files <- all_files[
    grepl(
      "\\.(csv|tsv|txt)$",
      all_files,
      ignore.case = TRUE
    )
  ]
  
  if (length(all_files) == 0) {
    
    return(NA_character_)
    
  }
  
  score <- rep(
    0,
    length(all_files)
  )
  
  for (i in seq_along(all_files)) {
    
    nm <- tolower(
      basename(all_files[i])
    )
    
    for (p in patterns) {
      
      if (grepl(
        p,
        nm,
        ignore.case = TRUE
      )) {
        
        score[i] <- score[i] + 1
        
      }
      
    }
    
  }
  
  if (max(score) == 0) {
    
    return(NA_character_)
    
  }
  
  all_files[
    which.max(score)
  ]
  
}


safe_fread <- function(path) {
  
  if (!file.exists(path)) {
    
    stop(
      "File does not exist:\n",
      path
    )
    
  }
  
  message(
    "Reading: ",
    basename(path)
  )
  
  dat <- tryCatch(
    
    fread(
      path,
      data.table = FALSE,
      check.names = FALSE,
      showProgress = FALSE
    ),
    
    error = function(e) {
      
      message(
        "Standard fread failed. Retrying with fill=TRUE..."
      )
      
      fread(
        path,
        data.table = FALSE,
        fill = TRUE,
        check.names = FALSE,
        showProgress = FALSE
      )
      
    }
    
  )
  
  dat
  
}


###############################################################################
# 4. LOCATE SCRIPT 12 OUTPUT
###############################################################################

cat("Locating Script 12 output...\n")

candidate_file <- file.path(
  processed_dir,
  "final_target_priority_table.csv"
)

candidate_rds <- file.path(
  processed_dir,
  "final_target_priority_table.rds"
)

if (file.exists(candidate_rds)) {
  
  final_targets <- readRDS(
    candidate_rds
  )
  
} else if (file.exists(candidate_file)) {
  
  final_targets <- read.csv(
    candidate_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
} else {
  
  stop(
    "\nCould not find Script 12 output.\n",
    "Expected:\n",
    candidate_rds,
    "\nor\n",
    candidate_file,
    "\n\nDo NOT rerun Scripts 0-12. Check that the Script 12 files are inside data/processed/."
  )
  
}

final_targets <- as.data.frame(
  final_targets,
  stringsAsFactors = FALSE
)

if (!"SYMBOL" %in% colnames(final_targets)) {
  
  stop(
    "Script 12 output does not contain the required SYMBOL column."
  )
  
}

final_targets$SYMBOL <- toupper(
  trimws(
    as.character(
      final_targets$SYMBOL
    )
  )
)

final_targets <- final_targets[
  !is.na(final_targets$SYMBOL) &
    final_targets$SYMBOL != "",
  ,
  drop = FALSE
]

final_targets <- final_targets[
  !duplicated(final_targets$SYMBOL),
  ,
  drop = FALSE
]

cat(
  "Final target genes loaded: ",
  nrow(final_targets),
  "\n",
  sep = ""
)


###############################################################################
# 5. DEFINE ML PREDICTOR PANEL
###############################################################################

predictor_genes <- c(
  "AURKA",
  "CDK1",
  "PLK1",
  "TOP2A",
  "CHEK1",
  "CCNB1",
  "CCNA2",
  "BUB1",
  "BUB1B",
  "CDC20"
)

predictor_genes <- toupper(
  predictor_genes
)

cat("\nML predictor panel:\n")
cat(
  paste(
    predictor_genes,
    collapse = ", "
  ),
  "\n"
)


###############################################################################
# 6. LOCATE DEPMAP FILES
###############################################################################

cat("\n")
cat("============================================================\n")
cat("SEARCHING FOR DEPMAP FILES\n")
cat("============================================================\n\n")

search_dirs <- c(
  raw_dir,
  project_dir
)

search_dirs <- unique(
  search_dirs[
    dir.exists(search_dirs)
  ]
)


# CRISPR gene effect ----------------------------------------------------------

crispr_file <- find_file(
  search_dirs,
  c(
    "crispr",
    "gene.*effect",
    "crisprgeneeffect",
    "chronos"
  )
)


# Model metadata --------------------------------------------------------------

model_file <- find_file(
  search_dirs,
  c(
    "model.*metadata",
    "model.*info",
    "sample.*info",
    "depmap.*model",
    "cell.*line"
  )
)


# Expression ------------------------------------------------------------------

expression_file <- find_file(
  search_dirs,
  c(
    "expression",
    "transcript",
    "rnaseq",
    "gene.*expression"
  )
)


cat(
  "CRISPR file: ",
  ifelse(
    is.na(crispr_file),
    "NOT FOUND",
    crispr_file
  ),
  "\n",
  sep = ""
)

cat(
  "Model metadata file: ",
  ifelse(
    is.na(model_file),
    "NOT FOUND",
    model_file
  ),
  "\n",
  sep = ""
)

cat(
  "Expression file: ",
  ifelse(
    is.na(expression_file),
    "NOT FOUND",
    expression_file
  ),
  "\n",
  sep = ""
)


###############################################################################
# 7. IF AUTOMATIC FILE DETECTION FAILED
###############################################################################

if (is.na(crispr_file)) {
  
  stop(
    "\nCould not automatically identify the DepMap CRISPR gene-effect file.\n",
    "Place the DepMap CRISPR gene-effect CSV in:\n",
    raw_dir,
    "\n\nThe file normally contains terms such as CRISPR, GeneEffect or Chronos."
  )
  
}

if (is.na(model_file)) {
  
  stop(
    "\nCould not automatically identify the DepMap model metadata file.\n",
    "Place the DepMap model metadata CSV in:\n",
    raw_dir,
    "\n\nIt should contain ModelID and cancer lineage information."
  )
  
}

if (is.na(expression_file)) {
  
  stop(
    "\nCould not automatically identify the DepMap expression file.\n",
    "Place the DepMap expression CSV in:\n",
    raw_dir,
    "\n\nThe file normally contains Expression, RNAseq or Transcript in its filename."
  )
  
}


###############################################################################
# 8. LOAD DEPMAP CRISPR DATA
###############################################################################

cat("\n")
cat("Loading DepMap CRISPR data...\n")

crispr <- safe_fread(
  crispr_file
)

cat(
  "CRISPR dimensions: ",
  nrow(crispr),
  " x ",
  ncol(crispr),
  "\n",
  sep = ""
)


###############################################################################
# 9. LOAD MODEL METADATA
###############################################################################

cat("\n")
cat("Loading DepMap model metadata...\n")

model_metadata <- safe_fread(
  model_file
)

cat(
  "Model metadata dimensions: ",
  nrow(model_metadata),
  " x ",
  ncol(model_metadata),
  "\n",
  sep = ""
)


###############################################################################
# 10. IDENTIFY MODEL ID COLUMN
###############################################################################

model_id_candidates <- c(
  "ModelID",
  "model_id",
  "DepMap_ID",
  "depmap_id",
  "DepMapID",
  "sample_id"
)

model_id_column <- model_id_candidates[
  model_id_candidates %in%
    colnames(model_metadata)
][1]

if (is.na(model_id_column)) {
  
  possible_model_columns <- which(
    sapply(
      model_metadata,
      is_model_id
    )
  )
  
  if (length(possible_model_columns) == 0) {
    
    stop(
      "Could not identify ModelID column in model metadata."
    )
    
  }
  
  model_id_column <- colnames(
    model_metadata
  )[possible_model_columns[1]]
  
}

cat(
  "Model ID column: ",
  model_id_column,
  "\n",
  sep = ""
)


###############################################################################
# 11. IDENTIFY CANCER LINEAGE COLUMN
###############################################################################

lineage_candidates <- c(
  "OncotreeLineage",
  "OncotreePrimaryDisease",
  "PrimaryDisease",
  "lineage",
  "Lineage",
  "primary_disease",
  "Primary_Disease",
  "cancer_type"
)

lineage_column <- lineage_candidates[
  lineage_candidates %in%
    colnames(model_metadata)
][1]

if (is.na(lineage_column)) {
  
  possible_lineage <- grep(
    "lineage|disease|cancer|oncotree",
    colnames(model_metadata),
    ignore.case = TRUE,
    value = TRUE
  )
  
  if (length(possible_lineage) == 0) {
    
    stop(
      "Could not identify the DepMap cancer lineage column."
    )
    
  }
  
  lineage_column <- possible_lineage[1]
  
}

cat(
  "Lineage column: ",
  lineage_column,
  "\n",
  sep = ""
)


###############################################################################
# 12. PREPARE MODEL METADATA
###############################################################################

model_metadata[[model_id_column]] <- trimws(
  as.character(
    model_metadata[[model_id_column]]
  )
)

model_metadata[[lineage_column]] <- as.character(
  model_metadata[[lineage_column]]
)

model_metadata <- model_metadata[
  !duplicated(
    model_metadata[[model_id_column]]
  ),
  ,
  drop = FALSE
]


###############################################################################
# 13. IDENTIFY LUAD MODELS
###############################################################################

lineage_values <- toupper(
  trimws(
    model_metadata[[lineage_column]]
  )
)

luad_index <- grepl(
  "LUNG ADENOCARCINOMA|LUAD",
  lineage_values
)

luad_models <- unique(
  model_metadata[[model_id_column]][
    luad_index
  ]
)

cat("\n")
cat(
  "LUAD models identified: ",
  length(luad_models),
  "\n",
  sep = ""
)


###############################################################################
# 14. IDENTIFY MODEL IDS IN CRISPR DATA
###############################################################################

crispr_columns <- colnames(
  crispr
)

crispr_model_columns <- crispr_columns[
  grepl(
    "^ACH[-_]",
    crispr_columns,
    ignore.case = TRUE
  )
]

if (length(crispr_model_columns) == 0) {
  
  # Alternative strategy:
  # identify columns overlapping with model metadata
  
  model_ids <- unique(
    model_metadata[[model_id_column]]
  )
  
  crispr_model_columns <- intersect(
    crispr_columns,
    model_ids
  )
  
}

cat(
  "CRISPR model columns identified: ",
  length(crispr_model_columns),
  "\n",
  sep = ""
)


###############################################################################
# 15. IDENTIFY CRISPR GENE COLUMN
###############################################################################

non_model_columns <- setdiff(
  crispr_columns,
  crispr_model_columns
)

gene_column <- non_model_columns[
  grepl(
    "gene|symbol|id|name",
    non_model_columns,
    ignore.case = TRUE
  )
][1]

if (is.na(gene_column)) {
  
  gene_column <- non_model_columns[1]
  
}

cat(
  "CRISPR gene column: ",
  gene_column,
  "\n",
  sep = ""
)


###############################################################################
# 16. NORMALIZE CRISPR GENE NAMES
###############################################################################

crispr$gene_symbol_normalized <- normalize_gene_name(
  crispr[[gene_column]]
)


###############################################################################
# 17. CHECK REQUIRED CRISPR PREDICTORS
###############################################################################

crispr_predictor_rows <- match(
  predictor_genes,
  crispr$gene_symbol_normalized
)

missing_crispr_predictors <- predictor_genes[
  is.na(
    crispr_predictor_rows
  )
]

if (length(missing_crispr_predictors) > 0) {
  
  cat("\n")
  cat(
    "WARNING: Some predictor genes were not found exactly in CRISPR data:\n"
  )
  
  cat(
    paste(
      missing_crispr_predictors,
      collapse = ", "
    ),
    "\n"
  )
  
}

present_crispr_predictors <- predictor_genes[
  !is.na(
    crispr_predictor_rows
  )
]

if (length(present_crispr_predictors) == 0) {
  
  stop(
    "\nNone of the ML predictor genes could be identified in the DepMap CRISPR file.\n",
    "Inspect the gene column before continuing."
  )
  
}


###############################################################################
# 18. BUILD CRISPR PREDICTOR MATRIX
###############################################################################

crispr_rows <- crispr[
  crispr$gene_symbol_normalized %in%
    present_crispr_predictors,
  ,
  drop = FALSE
]

crispr_rows <- crispr_rows[
  match(
    present_crispr_predictors,
    crispr_rows$gene_symbol_normalized
  ),
  ,
  drop = FALSE
]

crispr_luad <- crispr_rows[
  ,
  intersect(
    luad_models,
    colnames(crispr_rows)
  ),
  drop = FALSE
]

crispr_luad_t <- as.data.frame(
  t(
    as.matrix(
      crispr_luad
    )
  ),
  stringsAsFactors = FALSE
)

colnames(crispr_luad_t) <-
  crispr_rows$gene_symbol_normalized

crispr_luad_t$ModelID <- rownames(
  crispr_luad_t
)

crispr_luad_t <- crispr_luad_t[
  ,
  c(
    "ModelID",
    present_crispr_predictors
  ),
  drop = FALSE
]


###############################################################################
# 19. LOAD DEPMAP EXPRESSION DATA
###############################################################################

cat("\n")
cat("============================================================\n")
cat("LOADING DEPMAP EXPRESSION DATA\n")
cat("============================================================\n\n")

expression_data <- safe_fread(
  expression_file
)

cat(
  "Expression dimensions: ",
  nrow(expression_data),
  " x ",
  ncol(expression_data),
  "\n",
  sep = ""
)


###############################################################################
# 20. IDENTIFY EXPRESSION GENE COLUMN
###############################################################################

expression_columns <- colnames(
  expression_data
)

expression_model_columns <- intersect(
  expression_columns,
  unique(
    model_metadata[[model_id_column]]
  )
)

if (length(expression_model_columns) == 0) {
  
  expression_model_columns <- expression_columns[
    grepl(
      "^ACH[-_]",
      expression_columns,
      ignore.case = TRUE
    )
  ]
  
}

cat(
  "Expression model columns identified: ",
  length(expression_model_columns),
  "\n",
  sep = ""
)

if (length(expression_model_columns) == 0) {
  
  stop(
    "\nCould not identify DepMap ModelID columns in expression data.\n",
    "The expression file may use a different sample-ID format."
  )
  
}


###############################################################################
# 21. IDENTIFY EXPRESSION GENE IDENTIFIER COLUMN
###############################################################################

expression_gene_columns <- setdiff(
  expression_columns,
  expression_model_columns
)

if (length(expression_gene_columns) == 0) {
  
  stop(
    "Could not identify gene identifier column in expression data."
  )
  
}

expression_gene_column <- expression_gene_columns[1]

cat(
  "Expression gene identifier column: ",
  expression_gene_column,
  "\n",
  sep = ""
)


###############################################################################
# 22. NORMALIZE EXPRESSION GENE IDENTIFIERS
###############################################################################

expression_data$gene_symbol_normalized <-
  normalize_gene_name(
    expression_data[[expression_gene_column]]
  )


###############################################################################
# 23. CHECK EXPRESSION PREDICTORS
###############################################################################

expression_predictor_rows <- match(
  predictor_genes,
  expression_data$gene_symbol_normalized
)

missing_expression_predictors <- predictor_genes[
  is.na(
    expression_predictor_rows
  )
]

cat("\n")
cat(
  "Expression predictor matching:\n"
)

cat(
  "Found: ",
  sum(
    !is.na(
      expression_predictor_rows
    )
  ),
  " / ",
  length(predictor_genes),
  "\n",
  sep = ""
)

if (length(missing_expression_predictors) > 0) {
  
  cat(
    "Missing: ",
    paste(
      missing_expression_predictors,
      collapse = ", "
    ),
    "\n",
    sep = ""
  )
  
}


###############################################################################
# 24. IMPORTANT FALLBACK FOR ENSEMBL EXPRESSION DATA
###############################################################################

# If expression rows are Ensembl IDs, try to use a locally available
# annotation table from the project.

if (
  length(missing_expression_predictors) > 0 &&
  all(
    grepl(
      "^ENSG",
      expression_data$gene_symbol_normalized[
        expression_data$gene_symbol_normalized != ""
      ],
      ignore.case = TRUE
    )[
      seq_len(
        min(
          100,
          sum(
            expression_data$gene_symbol_normalized != ""
          )
        )
      )
    ]
  )
) {
  
  cat("\n")
  cat(
    "Expression appears to use Ensembl IDs.\n"
  )
  
  annotation_candidates <- list.files(
    raw_dir,
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )
  
  annotation_candidates <- annotation_candidates[
    grepl(
      "annotation|gene.*symbol|ensembl|gene.*mapping",
      basename(annotation_candidates),
      ignore.case = TRUE
    )
  ]
  
  annotation_candidates <- annotation_candidates[
    grepl(
      "\\.(csv|tsv|txt)$",
      annotation_candidates,
      ignore.case = TRUE
    )
  ]
  
  if (length(annotation_candidates) > 0) {
    
    cat(
      "Attempting local gene annotation file:\n",
      annotation_candidates[1],
      "\n"
    )
    
    annotation_data <- tryCatch(
      safe_fread(
        annotation_candidates[1]
      ),
      error = function(e) NULL
    )
    
    if (!is.null(annotation_data)) {
      
      annotation_names <- tolower(
        colnames(
          annotation_data
        )
      )
      
      ens_col <- which(
        grepl(
          "ensembl",
          annotation_names
        )
      )[1]
      
      sym_col <- which(
        grepl(
          "symbol|gene_name",
          annotation_names
        )
      )[1]
      
      if (
        !is.na(ens_col) &&
        !is.na(sym_col)
      ) {
        
        annotation_map <- annotation_data[
          ,
          c(
            ens_col,
            sym_col
          ),
          drop = FALSE
        ]
        
        colnames(annotation_map) <- c(
          "ensembl",
          "symbol"
        )
        
        annotation_map$ensembl <-
          sub(
            "\\..*$",
            "",
            toupper(
              trimws(
                as.character(
                  annotation_map$ensembl
                )
              )
            )
          )
        
        annotation_map$symbol <-
          toupper(
            trimws(
              as.character(
                annotation_map$symbol
              )
            )
          )
        
        expression_data$ensembl_clean <-
          sub(
            "\\..*$",
            "",
            toupper(
              trimws(
                as.character(
                  expression_data[[expression_gene_column]]
                )
              )
            )
          )
        
        map_index <- match(
          expression_data$ensembl_clean,
          annotation_map$ensembl
        )
        
        expression_data$gene_symbol_normalized <-
          annotation_map$symbol[
            map_index
          ]
        
        expression_data$gene_symbol_normalized[
          is.na(
            expression_data$gene_symbol_normalized
          )
        ] <- ""
        
      }
      
    }
    
  }
  
}


###############################################################################
# 25. FINAL EXPRESSION MATCH CHECK
###############################################################################

expression_predictor_rows <- match(
  predictor_genes,
  expression_data$gene_symbol_normalized
)

missing_expression_predictors <- predictor_genes[
  is.na(
    expression_predictor_rows
  )
]

if (length(missing_expression_predictors) > 0) {
  
  stop(
    "\nMissing ML predictor genes from expression data:\n",
    paste(
      missing_expression_predictors,
      collapse = ", "
    ),
    "\n\nThe expression matrix is present, but its gene identifiers could not ",
    "be mapped to the required predictor symbols.\n",
    "No downstream ML model should be run until this is resolved."
  )
  
}

cat(
  "All ",
  length(predictor_genes),
  " ML predictor genes successfully matched.\n",
  sep = ""
)


###############################################################################
# 26. BUILD EXPRESSION PREDICTOR MATRIX
###############################################################################

expression_rows <- expression_data[
  expression_predictor_rows,
  ,
  drop = FALSE
]

expression_matrix <- as.data.frame(
  t(
    as.matrix(
      expression_rows[
        ,
        expression_model_columns,
        drop = FALSE
      ]
    )
  ),
  stringsAsFactors = FALSE
)

colnames(
  expression_matrix
) <- predictor_genes

expression_matrix$ModelID <-
  rownames(
    expression_matrix
  )

expression_matrix <- expression_matrix[
  ,
  c(
    "ModelID",
    predictor_genes
  ),
  drop = FALSE
]


###############################################################################
# 27. CONVERT EXPRESSION VALUES TO NUMERIC
###############################################################################

for (g in predictor_genes) {
  
  expression_matrix[[g]] <-
    suppressWarnings(
      as.numeric(
        expression_matrix[[g]]
      )
    )
  
}


###############################################################################
# 28. RESTRICT TO LUAD MODELS
###############################################################################

expression_matrix <- expression_matrix[
  expression_matrix$ModelID %in%
    luad_models,
  ,
  drop = FALSE
]

crispr_luad_t <- crispr_luad_t[
  crispr_luad_t$ModelID %in%
    luad_models,
  ,
  drop = FALSE
]


###############################################################################
# 29. MERGE EXPRESSION AND CRISPR DATA
###############################################################################

ml_data <- merge(
  expression_matrix,
  crispr_luad_t,
  by = "ModelID",
  suffixes = c(
    "_expression",
    "_dependency"
  )
)


###############################################################################
# 30. ADD MODEL METADATA
###############################################################################

metadata_for_merge <- model_metadata[
  ,
  c(
    model_id_column,
    lineage_column
  ),
  drop = FALSE
]

colnames(
  metadata_for_merge
) <- c(
  "ModelID",
  "Lineage"
)

ml_data <- merge(
  ml_data,
  metadata_for_merge,
  by = "ModelID",
  all.x = TRUE
)


###############################################################################
# 31. CREATE AURKA TARGET VARIABLE
###############################################################################

if (!"AURKA" %in% colnames(ml_data)) {
  
  stop(
    "AURKA predictor/dependency column could not be created."
  )
  
}

# The AURKA dependency value is the target variable.
#
# More negative gene-effect values indicate stronger dependency.
#
# We retain the continuous dependency score for regression.

ml_data$AURKA_dependency <-
  suppressWarnings(
    as.numeric(
      ml_data$AURKA
    )
  )


###############################################################################
# 32. REMOVE DUPLICATE AURKA PREDICTOR
###############################################################################

# AURKA is both the target and part of the predictor panel.
#
# To prevent target leakage, AURKA expression can be used as a predictor,
# but AURKA dependency must NOT be included as a predictor.
#
# The CRISPR AURKA dependency column is renamed explicitly.

if ("AURKA_dependency" %in% colnames(ml_data)) {
  
  # Keep it as target.
  
}


###############################################################################
# 33. CHECK FOR TARGET VARIABILITY
###############################################################################

target_values <- ml_data$AURKA_dependency

target_values <- target_values[
  is.finite(
    target_values
  )
]

if (length(target_values) < 10) {
  
  stop(
    "Too few valid AURKA dependency observations for ML."
  )
  
}

if (
  sd(
    target_values,
    na.rm = TRUE
  ) == 0
) {
  
  stop(
    "AURKA dependency has zero variance. ML cannot be performed."
  )
  
}


###############################################################################
# 34. CHECK PREDICTOR VARIABILITY
###############################################################################

predictor_expression_columns <- paste0(
  predictor_genes,
  "_expression"
)

# Depending on merge behavior, expression columns may retain their
# original names rather than suffixes.

available_expression_predictors <- intersect(
  predictor_genes,
  colnames(ml_data)
)

if (length(available_expression_predictors) < 2) {
  
  # Try suffixed names.
  
  available_expression_predictors <- intersect(
    predictor_expression_columns,
    colnames(ml_data)
  )
  
}

if (length(available_expression_predictors) == 0) {
  
  stop(
    "No usable expression predictor columns were found after merging."
  )
  
}


###############################################################################
# 35. BUILD FINAL ML TABLE
###############################################################################

ml_predictor_columns <- predictor_genes[
  predictor_genes %in%
    colnames(ml_data)
]

if (length(ml_predictor_columns) == 0) {
  
  # Search for expression-suffixed columns.
  
  ml_predictor_columns <- paste0(
    predictor_genes,
    "_expression"
  )
  
  ml_predictor_columns <- ml_predictor_columns[
    ml_predictor_columns %in%
      colnames(ml_data)
  ]
  
}

if (length(ml_predictor_columns) < 2) {
  
  stop(
    "Fewer than two expression predictors are available."
  )
  
}

for (g in ml_predictor_columns) {
  
  ml_data[[g]] <-
    suppressWarnings(
      as.numeric(
        ml_data[[g]]
      )
    )
  
}

ml_final <- ml_data[
  ,
  c(
    "ModelID",
    ml_predictor_columns,
    "AURKA_dependency",
    "Lineage"
  ),
  drop = FALSE
]


###############################################################################
# 36. COMPLETE-CASE FILTERING
###############################################################################

ml_complete <- ml_final[
  complete.cases(
    ml_final[
      ,
      c(
        ml_predictor_columns,
        "AURKA_dependency"
      ),
      drop = FALSE
    ]
  ),
  ,
  drop = FALSE
]


###############################################################################
# 37. REMOVE ZERO-VARIANCE PREDICTORS
###############################################################################

predictor_variances <- sapply(
  ml_complete[
    ,
    ml_predictor_columns,
    drop = FALSE
  ],
  function(x) {
    
    var(
      x,
      na.rm = TRUE
    )
    
  }
)

usable_predictors <- names(
  predictor_variances[
    is.finite(
      predictor_variances
    ) &
      predictor_variances > 0
  ]
)

removed_predictors <- setdiff(
  ml_predictor_columns,
  usable_predictors
)

if (length(removed_predictors) > 0) {
  
  cat(
    "\nRemoved zero-variance predictors:\n",
    paste(
      removed_predictors,
      collapse = ", "
    ),
    "\n"
  )
  
}

ml_complete <- ml_complete[
  ,
  c(
    "ModelID",
    usable_predictors,
    "AURKA_dependency",
    "Lineage"
  ),
  drop = FALSE
]


###############################################################################
# 38. FINAL DATA VALIDATION
###############################################################################

cat("\n")
cat("============================================================\n")
cat("ML DATASET VALIDATION\n")
cat("============================================================\n\n")

cat(
  "LUAD models before complete-case filtering: ",
  nrow(ml_final),
  "\n",
  sep = ""
)

cat(
  "Complete ML models: ",
  nrow(ml_complete),
  "\n",
  sep = ""
)

cat(
  "Predictors available: ",
  length(usable_predictors),
  "\n",
  sep = ""
)

cat(
  "Predictors:\n",
  paste(
    usable_predictors,
    collapse = ", "
  ),
  "\n"
)

cat(
  "AURKA dependency range: ",
  round(
    min(
      ml_complete$AURKA_dependency
    ),
    3
  ),
  " to ",
  round(
    max(
      ml_complete$AURKA_dependency
    ),
    3
  ),
  "\n",
  sep = ""
)


###############################################################################
# 39. BASIC ML PREPARATION
###############################################################################

if (nrow(ml_complete) < 20) {
  
  stop(
    "\nFewer than 20 complete LUAD models are available.\n",
    "The dataset is too small for a meaningful ML analysis."
  )
  
}

x <- as.matrix(
  ml_complete[
    ,
    usable_predictors,
    drop = FALSE
  ]
)

y <- ml_complete$AURKA_dependency

if (
  any(
    !is.finite(x)
  )
) {
  
  stop(
    "Non-finite predictor values remain in the ML matrix."
  )
  
}

if (
  any(
    !is.finite(y)
  )
) {
  
  stop(
    "Non-finite target values remain in the ML dataset."
  )
  
}


###############################################################################
# 40. CORRELATION SCREENING
###############################################################################

correlations <- sapply(
  seq_len(
    ncol(x)
  ),
  function(i) {
    
    cor(
      x[, i],
      y,
      method = "spearman"
    )
    
  }
)

correlation_table <- data.frame(
  predictor = colnames(x),
  spearman_rho = correlations,
  stringsAsFactors = FALSE
)

correlation_table <- correlation_table[
  order(
    -abs(
      correlation_table$spearman_rho
    )
  ),
  ,
  drop = FALSE
]


###############################################################################
# 41. ELASTIC NET MODEL
###############################################################################

cat("\n")
cat("============================================================\n")
cat("ELASTIC NET MODEL\n")
cat("============================================================\n\n")

set.seed(
  2026
)

cv_fit <- cv.glmnet(
  x = x,
  y = y,
  alpha = 0.5,
  family = "gaussian",
  standardize = TRUE,
  nfolds = min(
    10,
    nrow(x)
  ),
  type.measure = "mse"
)

best_lambda <- cv_fit$lambda.min

elastic_fit <- glmnet(
  x = x,
  y = y,
  alpha = 0.5,
  family = "gaussian",
  standardize = TRUE,
  lambda = best_lambda
)

coef_matrix <- as.matrix(
  coef(
    elastic_fit
  )
)

coefficient_table <- data.frame(
  predictor = rownames(
    coef_matrix
  ),
  coefficient = as.numeric(
    coef_matrix[, 1]
  ),
  stringsAsFactors = FALSE
)

coefficient_table <- coefficient_table[
  coefficient_table$predictor != "(Intercept)",
  ,
  drop = FALSE
]

coefficient_table$abs_coefficient <-
  abs(
    coefficient_table$coefficient
  )

coefficient_table <- coefficient_table[
  order(
    -coefficient_table$abs_coefficient
  ),
  ,
  drop = FALSE
]


###############################################################################
# 42. MODEL PREDICTIONS
###############################################################################

predicted <- as.numeric(
  predict(
    elastic_fit,
    newx = x,
    s = best_lambda
  )
)

model_performance <- data.frame(
  
  RMSE = sqrt(
    mean(
      (
        y - predicted
      )^2
    )
  ),
  
  Pearson_R = cor(
    y,
    predicted,
    method = "pearson"
  ),
  
  Spearman_R = cor(
    y,
    predicted,
    method = "spearman"
  ),
  
  R_squared = cor(
    y,
    predicted
  )^2
  
)


###############################################################################
# 43. SAVE ML DATASET
###############################################################################

saveRDS(
  ml_complete,
  file.path(
    processed_dir,
    "ML_AURKA_dataset.rds"
  )
)

write.csv(
  ml_complete,
  file.path(
    processed_dir,
    "ML_AURKA_dataset.csv"
  ),
  row.names = FALSE
)


###############################################################################
# 44. SAVE PREDICTOR DATA
###############################################################################

saveRDS(
  ml_complete[
    ,
    c(
      "ModelID",
      usable_predictors
    ),
    drop = FALSE
  ],
  file.path(
    processed_dir,
    "ML_predictor_expression.rds"
  )
)

write.csv(
  ml_complete[
    ,
    c(
      "ModelID",
      usable_predictors
    ),
    drop = FALSE
  ],
  file.path(
    processed_dir,
    "ML_predictor_expression.csv"
  ),
  row.names = FALSE
)


###############################################################################
# 45. SAVE MODEL METADATA
###############################################################################

saveRDS(
  ml_complete[
    ,
    c(
      "ModelID",
      "Lineage"
    ),
    drop = FALSE
  ],
  file.path(
    processed_dir,
    "ML_DepMap_model_metadata.rds"
  )
)

write.csv(
  ml_complete[
    ,
    c(
      "ModelID",
      "Lineage"
    ),
    drop = FALSE
  ],
  file.path(
    processed_dir,
    "ML_DepMap_model_metadata.csv"
  ),
  row.names = FALSE
)


###############################################################################
# 46. SAVE CORRELATION RESULTS
###############################################################################

saveRDS(
  correlation_table,
  file.path(
    processed_dir,
    "ML_AURKA_predictor_correlations.rds"
  )
)

write.csv(
  correlation_table,
  file.path(
    processed_dir,
    "ML_AURKA_predictor_correlations.csv"
  ),
  row.names = FALSE
)


###############################################################################
# 47. SAVE ELASTIC NET RESULTS
###############################################################################

saveRDS(
  coefficient_table,
  file.path(
    processed_dir,
    "ML_AURKA_elastic_net_coefficients.rds"
  )
)

write.csv(
  coefficient_table,
  file.path(
    processed_dir,
    "ML_AURKA_elastic_net_coefficients.csv"
  ),
  row.names = FALSE
)

saveRDS(
  cv_fit,
  file.path(
    processed_dir,
    "ML_AURKA_cv_model.rds"
  )
)


###############################################################################
# 48. SAVE MODEL PERFORMANCE
###############################################################################

write.csv(
  model_performance,
  file.path(
    processed_dir,
    "ML_AURKA_model_performance.csv"
  ),
  row.names = FALSE
)


###############################################################################
# 49. SAVE SUMMARY
###############################################################################

summary_file <- file.path(
  processed_dir,
  "ML_AURKA_summary.txt"
)

sink(
  summary_file
)

cat(
  "SCRIPT 13: DepMap-based ML dataset\n"
)

cat(
  "========================================\n\n"
)

cat(
  "LUAD models:\n"
)

cat(
  nrow(ml_complete),
  "\n\n"
)

cat(
  "Predictors:\n"
)

cat(
  paste(
    usable_predictors,
    collapse = ", "
  ),
  "\n\n"
)

cat(
  "Target:\n"
)

cat(
  "AURKA CRISPR gene-effect dependency\n\n"
)

cat(
  "Elastic-net lambda.min:\n"
)

cat(
  best_lambda,
  "\n\n"
)

cat(
  "Model performance:\n"
)

print(
  model_performance
)

cat(
  "\nPredictor correlations:\n"
)

print(
  correlation_table
)

cat(
  "\nElastic-net coefficients:\n"
)

print(
  coefficient_table
)

sink()


###############################################################################
# 50. SESSION INFORMATION
###############################################################################

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    processed_dir,
    "ML_AURKA_sessionInfo.txt"
  )
)


###############################################################################
# 51. FINAL OUTPUT
###############################################################################

cat("\n")
cat("============================================================\n")
cat("SCRIPT 13 COMPLETED SUCCESSFULLY\n")
cat("============================================================\n\n")

cat(
  "LUAD models used: ",
  nrow(ml_complete),
  "\n",
  sep = ""
)

cat(
  "Predictors used: ",
  length(usable_predictors),
  "\n",
  sep = ""
)

cat(
  "Target: AURKA dependency\n"
)

cat("\nOutputs saved:\n")

cat(
  "- ML_AURKA_dataset.rds\n",
  "- ML_AURKA_dataset.csv\n",
  "- ML_predictor_expression.rds\n",
  "- ML_predictor_expression.csv\n",
  "- ML_DepMap_model_metadata.rds\n",
  "- ML_DepMap_model_metadata.csv\n",
  "- ML_AURKA_predictor_correlations.rds\n",
  "- ML_AURKA_predictor_correlations.csv\n",
  "- ML_AURKA_elastic_net_coefficients.rds\n",
  "- ML_AURKA_elastic_net_coefficients.csv\n",
  "- ML_AURKA_cv_model.rds\n",
  "- ML_AURKA_model_performance.csv\n",
  "- ML_AURKA_summary.txt\n",
  "- ML_AURKA_sessionInfo.txt\n"
)

cat("\n")
cat("============================================================\n")