# ============================================================
# SCRIPT 12
# FINAL TARGET EVIDENCE INTEGRATION
# ============================================================
#
# Purpose:
# Integrate:
#   1. TCGA-LUAD differential expression
#   2. Independent GSE32863 validation
#   3. STRING network/module evidence
#   4. DepMap LUAD dependency/selectivity
#   5. Open Targets tractability
#   6. ChEMBL pharmacological evidence
#
# Output:
#   - final_target_evidence_integrated.rds
#   - final_target_priority_table.rds
#   - final_target_priority_table.csv
#
# Important:
#   - DepMap Gene Effect: more negative = stronger dependency
#   - Strong dependency requires >=10 LUAD models
#   - Selectivity requires >=10 LUAD and >=10 non-LUAD models
#   - Statistical selectivity uses Wilcoxon rank-sum test
#   - Multiple testing correction uses Benjamini-Hochberg FDR
#
# ============================================================


# ============================================================
# 1. PACKAGES
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
})


# ============================================================
# 2. DIRECTORIES
# ============================================================

project_root <- getwd()

processed_dir <- file.path(
  project_root,
  "data",
  "processed"
)

raw_depmap_dir <- file.path(
  project_root,
  "data",
  "raw",
  "DepMap"
)


# ============================================================
# 3. CANDIDATE GENES
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


# ============================================================
# 4. LOAD NETWORK / MODULE RESULTS
# ============================================================

network_candidates <- readRDS(
  file.path(
    processed_dir,
    "network_candidates_modules.rds"
  )
)

network_degree <- readRDS(
  file.path(
    processed_dir,
    "STRING_network_degree.rds"
  )
)

module_membership <- readRDS(
  file.path(
    processed_dir,
    "STRING_module_membership.rds"
  )
)

module_enrichment_summary <- readRDS(
  file.path(
    processed_dir,
    "STRING_module_enrichment_summary.rds"
  )
)


# ============================================================
# 5. LOAD DRUG / TRACTABILITY EVIDENCE
# ============================================================

opentargets_summary <- readRDS(
  file.path(
    processed_dir,
    "OpenTargets_tractability_summary.rds"
  )
)

chembl_summary <- readRDS(
  file.path(
    processed_dir,
    "ChEMBL_evidence_summary.rds"
  )
)


# ============================================================
# 6. LOAD TCGA DESEQ2 RESULTS
# ============================================================

tcga_results <- readRDS(
  file.path(
    processed_dir,
    "TCGA_LUAD_DESeq2_shrunken_results.rds"
  )
)

tcga_df <- as.data.frame(tcga_results)

tcga_df$Ensembl <- rownames(tcga_df)


# ============================================================
# 7. CANDIDATE ENSEMBL IDS
# ============================================================

candidate_ensembl <- c(
  AURKA = "ENSG00000087586",
  BRCA1 = "ENSG00000012048",
  CDCA5 = "ENSG00000146670",
  CDK1 = "ENSG00000170312",
  CHEK1 = "ENSG00000149554",
  DNA2 = "ENSG00000138346",
  GMNN = "ENSG00000112312",
  PLK1 = "ENSG00000166851",
  PRKDC = "ENSG00000253729",
  SGO1 = "ENSG00000129810",
  SKA1 = "ENSG00000154839",
  TOP2A = "ENSG00000131747"
)


tcga_lookup <- tibble(
  SYMBOL = names(candidate_ensembl),
  Ensembl = unname(candidate_ensembl)
) %>%
  dplyr::left_join(
    tcga_df %>%
      dplyr::select(
        Ensembl,
        TCGA_shrunken_log2FC =
          log2FoldChange,
        TCGA_FDR =
          padj
      ),
    by = "Ensembl"
  )


# ============================================================
# 8. LOAD GSE32863 VALIDATION RESULTS
# ============================================================

validation <- readRDS(
  file.path(
    processed_dir,
    "GSE32863_validation_results.rds"
  )
)

validation_df <- validation %>%
  dplyr::transmute(
    SYMBOL =
      as.character(gene),
    
    GSE32863_log2FC =
      as.numeric(log2FC_Tumor_vs_Normal),
    
    GSE32863_FDR =
      as.numeric(FDR),
    
    Validation_direction =
      as.character(Validation_direction),
    
    direction_concordant =
      as.logical(direction_concordant),
    
    validation_class =
      as.character(validation_class)
  ) %>%
  dplyr::mutate(
    GSE32863_validated =
      validation_class == "Validated"
  )


cat(
  "\nGSE32863 validated candidates:",
  sum(
    validation_df$GSE32863_validated,
    na.rm = TRUE
  ),
  "/",
  nrow(validation_df),
  "\n"
)


# ============================================================
# 9. LOAD DEPMAP
# ============================================================

dep_effect <- readr::read_csv(
  file.path(
    raw_depmap_dir,
    "CRISPRGeneEffect.csv"
  ),
  show_col_types = FALSE
)

dep_model <- readr::read_csv(
  file.path(
    raw_depmap_dir,
    "Model.csv"
  ),
  show_col_types = FALSE
)


# ============================================================
# 10. IDENTIFY MODEL ID COLUMNS
# ============================================================

dep_model_id_col <- names(dep_model)[1]

dep_effect_id_col <- names(dep_effect)[1]

cat(
  "\nDepMap effect model ID column:",
  dep_effect_id_col,
  "\n"
)

cat(
  "DepMap model table ID column:",
  dep_model_id_col,
  "\n"
)


# ============================================================
# 11. IDENTIFY LUAD MODELS
# ============================================================

luad_model_ids <- dep_model %>%
  dplyr::filter(
    grepl(
      "lung",
      OncotreePrimaryDisease,
      ignore.case = TRUE
    ),
    grepl(
      "adenocarcinoma",
      OncotreeSubtype,
      ignore.case = TRUE
    )
  ) %>%
  dplyr::pull(
    dplyr::all_of(dep_model_id_col)
  ) %>%
  as.character()


# Keep only LUAD models represented in CRISPR data

luad_model_ids <- intersect(
  luad_model_ids,
  as.character(
    dep_effect[[dep_effect_id_col]]
  )
)


cat(
  "LUAD models with CRISPR data:",
  length(luad_model_ids),
  "\n"
)


# ============================================================
# 12. IDENTIFY NON-LUAD MODELS
# ============================================================

all_depmap_ids <- as.character(
  dep_effect[[dep_effect_id_col]]
)

non_luad_model_ids <- setdiff(
  all_depmap_ids,
  luad_model_ids
)


cat(
  "Non-LUAD models with CRISPR data:",
  length(non_luad_model_ids),
  "\n"
)


# ============================================================
# 13. DEPMAP CANDIDATE COLUMNS
# ============================================================

candidate_entrez <- c(
  AURKA = "6790",
  BRCA1 = "672",
  CDCA5 = "113130",
  CDK1 = "983",
  CHEK1 = "1111",
  DNA2 = "1763",
  GMNN = "51053",
  PLK1 = "5347",
  PRKDC = "5591",
  SGO1 = "151648",
  SKA1 = "220134",
  TOP2A = "7153"
)


candidate_columns <- paste0(
  candidate_genes,
  " (",
  candidate_entrez[candidate_genes],
  ")"
)


missing_columns <- setdiff(
  candidate_columns,
  names(dep_effect)
)


if (length(missing_columns) > 0) {
  
  stop(
    paste(
      "Missing DepMap columns:",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  )
  
}


cat(
  "All candidate DepMap columns found.\n"
)


# ============================================================
# 14. EXTRACT DEPMAP DATA
# ============================================================

dep_long <- dep_effect %>%
  
  dplyr::select(
    ModelID =
      dplyr::all_of(dep_effect_id_col),
    
    dplyr::all_of(candidate_columns)
  ) %>%
  
  tidyr::pivot_longer(
    cols = -ModelID,
    names_to = "DepMap_column",
    values_to = "GeneEffect"
  ) %>%
  
  dplyr::mutate(
    
    SYMBOL =
      sub(
        " \\([0-9]+\\)$",
        "",
        DepMap_column
      ),
    
    ModelID =
      as.character(ModelID),
    
    GeneEffect =
      as.numeric(GeneEffect)
  )


# ============================================================
# 15. LUAD DEPENDENCY SUMMARY
# ============================================================

luad_effect <- dep_long %>%
  dplyr::filter(
    ModelID %in% luad_model_ids
  )


luad_dependency_summary <- luad_effect %>%
  dplyr::group_by(SYMBOL) %>%
  dplyr::summarise(
    
    n_models =
      sum(!is.na(GeneEffect)),
    
    median_gene_effect =
      ifelse(
        n_models > 0,
        median(
          GeneEffect,
          na.rm = TRUE
        ),
        NA_real_
      ),
    
    mean_gene_effect =
      ifelse(
        n_models > 0,
        mean(
          GeneEffect,
          na.rm = TRUE
        ),
        NA_real_
      ),
    
    fraction_effect_below_m05 =
      ifelse(
        n_models > 0,
        mean(
          GeneEffect < -0.5,
          na.rm = TRUE
        ),
        NA_real_
      ),
    
    fraction_effect_below_m1 =
      ifelse(
        n_models > 0,
        mean(
          GeneEffect < -1,
          na.rm = TRUE
        ),
        NA_real_
      ),
    
    .groups = "drop"
  )


# ============================================================
# 16. NON-LUAD DEPENDENCY SUMMARY
# ============================================================

non_luad_effect <- dep_long %>%
  dplyr::filter(
    ModelID %in% non_luad_model_ids
  )


non_luad_dependency_summary <- non_luad_effect %>%
  dplyr::group_by(SYMBOL) %>%
  dplyr::summarise(
    
    n_models_non_luad =
      sum(!is.na(GeneEffect)),
    
    median_gene_effect_non_luad =
      ifelse(
        n_models_non_luad > 0,
        median(
          GeneEffect,
          na.rm = TRUE
        ),
        NA_real_
      ),
    
    .groups = "drop"
  )


# ============================================================
# 17. DEPMAP SELECTIVITY
# ============================================================

dependency_comparison <- luad_dependency_summary %>%
  
  dplyr::left_join(
    non_luad_dependency_summary,
    by = "SYMBOL"
  ) %>%
  
  dplyr::mutate(
    
    # --------------------------------------------------------
    # Positive difference means stronger dependency in LUAD
    # --------------------------------------------------------
    
    LUAD_selectivity_difference =
      median_gene_effect_non_luad -
      median_gene_effect,
    
    # --------------------------------------------------------
    # Strong dependency
    # --------------------------------------------------------
    
    strong_LUAD_dependency =
      n_models >= 10 &
      median_gene_effect <= -1,
    
    # --------------------------------------------------------
    # Directional selectivity
    # --------------------------------------------------------
    
    positive_LUAD_selectivity =
      n_models >= 10 &
      n_models_non_luad >= 10 &
      LUAD_selectivity_difference > 0
  )


# ============================================================
# 18. STATISTICAL LUAD VS NON-LUAD COMPARISON
# ============================================================

dependency_pvalues <- list()


for (gene in candidate_genes) {
  
  luad_values <- luad_effect %>%
    dplyr::filter(
      SYMBOL == gene
    ) %>%
    dplyr::pull(
      GeneEffect
    ) %>%
    na.omit()
  
  
  non_luad_values <- non_luad_effect %>%
    dplyr::filter(
      SYMBOL == gene
    ) %>%
    dplyr::pull(
      GeneEffect
    ) %>%
    na.omit()
  
  
  # ----------------------------------------------------------
  # Require sufficient observations in BOTH groups
  # ----------------------------------------------------------
  
  if (
    length(luad_values) >= 10 &
    length(non_luad_values) >= 10
  ) {
    
    wt <- wilcox.test(
      luad_values,
      non_luad_values,
      alternative = "two.sided",
      exact = FALSE
    )
    
    p_value <- wt$p.value
    
  } else {
    
    p_value <- NA_real_
    
  }
  
  
  dependency_pvalues[[gene]] <- tibble(
    SYMBOL = gene,
    p_value = p_value
  )
}


dependency_pvalues <- bind_rows(
  dependency_pvalues
)


# ============================================================
# 19. BENJAMINI-HOCHBERG FDR
# ============================================================

dependency_pvalues <- dependency_pvalues %>%
  
  dplyr::mutate(
    FDR =
      p.adjust(
        p_value,
        method = "BH"
      )
  )


dependency_comparison <- dependency_comparison %>%
  
  dplyr::left_join(
    dependency_pvalues,
    by = "SYMBOL"
  )


# ============================================================
# 20. DEPMAP SANITY CHECK
# ============================================================

cat(
  "\n============================================\n"
)

cat(
  "DEPMAP COVERAGE AND SELECTIVITY\n"
)

cat(
  "============================================\n\n"
)


print(
  dependency_comparison %>%
    dplyr::select(
      SYMBOL,
      n_models,
      n_models_non_luad,
      median_gene_effect,
      median_gene_effect_non_luad,
      LUAD_selectivity_difference,
      p_value,
      FDR,
      strong_LUAD_dependency,
      positive_LUAD_selectivity
    ),
  n = Inf
)


# ============================================================
# 21. PREPARE ONE-ROW-PER-GENE NETWORK DATA
# ============================================================

# network_candidates can contain multiple rows per gene.
# We reduce it to one row per candidate gene.
#
# Module enrichment remains stored separately because one
# module can have multiple enriched pathways.


network_clean <- network_candidates %>%
  
  dplyr::filter(
    SYMBOL %in% candidate_genes
  ) %>%
  
  dplyr::group_by(
    SYMBOL
  ) %>%
  
  dplyr::summarise(
    
    module =
      dplyr::first(
        as.character(module)
      ),
    
    n_pathways =
      dplyr::first(n_pathways),
    
    .groups = "drop"
  )


# ============================================================
# 22. PREPARE NETWORK DEGREE DATA
# ============================================================

network_degree_clean <- network_degree %>%
  
  dplyr::filter(
    SYMBOL %in% candidate_genes
  ) %>%
  
  dplyr::group_by(
    SYMBOL
  ) %>%
  
  dplyr::summarise(
    
    degree =
      dplyr::first(degree),
    
    .groups = "drop"
  )


# ============================================================
# 23. PREPARE OPEN TARGETS DATA
# ============================================================

opentargets_clean <- opentargets_summary %>%
  
  dplyr::filter(
    SYMBOL %in% candidate_genes
  ) %>%
  
  dplyr::distinct(
    SYMBOL,
    .keep_all = TRUE
  )


# ============================================================
# 24. PREPARE CHEMBL DATA
# ============================================================

chembl_clean <- chembl_summary %>%
  
  dplyr::filter(
    SYMBOL %in% candidate_genes
  ) %>%
  
  dplyr::distinct(
    SYMBOL,
    .keep_all = TRUE
  )


# ============================================================
# 25. MERGE ALL EVIDENCE
# ============================================================

# Start with exactly one row per candidate gene.
# Every joined table is reduced to one row per gene.
#
# module_enrichment_summary is intentionally NOT joined here
# because it is one-to-many. It remains available as a
# separate saved analysis result.


final_evidence <- tibble(
  SYMBOL = candidate_genes
) %>%
  
  # ----------------------------------------------------------
# TCGA
# ----------------------------------------------------------

dplyr::left_join(
  tcga_lookup %>%
    dplyr::distinct(
      SYMBOL,
      .keep_all = TRUE
    ),
  by = "SYMBOL"
) %>%
  
  # ----------------------------------------------------------
# GSE32863
# ----------------------------------------------------------

dplyr::left_join(
  validation_df %>%
    dplyr::distinct(
      SYMBOL,
      .keep_all = TRUE
    ),
  by = "SYMBOL"
) %>%
  
  # ----------------------------------------------------------
# Network
# ----------------------------------------------------------

dplyr::left_join(
  network_clean,
  by = "SYMBOL"
) %>%
  
  # ----------------------------------------------------------
# Network degree
# ----------------------------------------------------------

dplyr::left_join(
  network_degree_clean,
  by = "SYMBOL"
) %>%
  
  # ----------------------------------------------------------
# Open Targets
# ----------------------------------------------------------

dplyr::left_join(
  opentargets_clean,
  by = "SYMBOL"
) %>%
  
  # ----------------------------------------------------------
# ChEMBL
# ----------------------------------------------------------

dplyr::left_join(
  chembl_clean,
  by = "SYMBOL"
) %>%
  
  # ----------------------------------------------------------
# DepMap
# ----------------------------------------------------------

dplyr::left_join(
  dependency_comparison %>%
    dplyr::distinct(
      SYMBOL,
      .keep_all = TRUE
    ),
  by = "SYMBOL"
)


# ============================================================
# 26. CORE EVIDENCE FLAGS
# ============================================================

final_evidence <- final_evidence %>%
  
  dplyr::mutate(
    
    # --------------------------------------------------------
    # TCGA transcriptomic evidence
    # --------------------------------------------------------
    
    TCGA_significant =
      !is.na(TCGA_FDR) &
      TCGA_FDR < 0.05,
    
    TCGA_effect_size =
      !is.na(TCGA_shrunken_log2FC) &
      abs(TCGA_shrunken_log2FC) >= 1,
    
    transcriptomic_evidence =
      TCGA_significant &
      TCGA_effect_size,
    
    
    # --------------------------------------------------------
    # Independent validation
    # --------------------------------------------------------
    
    GSE32863_validated =
      validation_class == "Validated",
    
    independent_validation =
      GSE32863_validated,
    
    
    # --------------------------------------------------------
    # Functional dependency
    # --------------------------------------------------------
    
    functional_dependency =
      !is.na(median_gene_effect) &
      n_models >= 10 &
      median_gene_effect <= -1,
    
    
    # --------------------------------------------------------
    # Directional LUAD selectivity
    # --------------------------------------------------------
    
    LUAD_selective =
      !is.na(LUAD_selectivity_difference) &
      n_models >= 10 &
      n_models_non_luad >= 10 &
      LUAD_selectivity_difference > 0,
    
    
    # --------------------------------------------------------
    # Statistically supported LUAD selectivity
    # --------------------------------------------------------
    
    LUAD_selective_FDR05 =
      !is.na(FDR) &
      n_models >= 10 &
      n_models_non_luad >= 10 &
      LUAD_selectivity_difference > 0 &
      FDR < 0.05
  )


# ============================================================
# 27. CONVERGENT CORE EVIDENCE
# ============================================================

# Core evidence count deliberately uses directional
# LUAD selectivity rather than the stricter FDR flag.
#
# The FDR-based flag is retained separately and reported
# explicitly so that statistical strength is not hidden.


final_evidence <- final_evidence %>%
  
  dplyr::mutate(
    
    core_evidence_count =
      
      rowSums(
        cbind(
          transcriptomic_evidence,
          independent_validation,
          functional_dependency,
          LUAD_selective
        ),
        na.rm = TRUE
      ),
    
    
    priority_tier =
      
      dplyr::case_when(
        
        core_evidence_count >= 4 ~
          "Tier 1: Strong convergent evidence",
        
        core_evidence_count == 3 ~
          "Tier 2: Strong candidate",
        
        core_evidence_count == 2 ~
          "Tier 3: Moderate candidate",
        
        TRUE ~
          "Tier 4: Limited convergence"
      )
  )


# ============================================================
# 28. FINAL PRIORITY TABLE
# ============================================================

final_priority_table <- final_evidence %>%
  
  dplyr::select(
    
    SYMBOL,
    
    # --------------------------------------------------------
    # Transcriptomics
    # --------------------------------------------------------
    
    TCGA_shrunken_log2FC,
    TCGA_FDR,
    
    # --------------------------------------------------------
    # Independent validation
    # --------------------------------------------------------
    
    GSE32863_log2FC,
    GSE32863_FDR,
    GSE32863_validated,
    
    # --------------------------------------------------------
    # Network
    # --------------------------------------------------------
    
    degree,
    module,
    n_pathways,
    
    # --------------------------------------------------------
    # DepMap
    # --------------------------------------------------------
    
    n_models,
    n_models_non_luad,
    
    median_gene_effect,
    median_gene_effect_non_luad,
    
    LUAD_selectivity_difference,
    
    p_value,
    FDR,
    
    strong_LUAD_dependency =
      functional_dependency,
    
    LUAD_selective,
    LUAD_selective_FDR05,
    
    # --------------------------------------------------------
    # Evidence
    # --------------------------------------------------------
    
    transcriptomic_evidence,
    independent_validation,
    functional_dependency,
    
    core_evidence_count,
    priority_tier
  ) %>%
  
  dplyr::arrange(
    
    dplyr::desc(
      core_evidence_count
    ),
    
    TCGA_FDR
  )


# ============================================================
# 29. FINAL SANITY CHECKS
# ============================================================

cat(
  "\n============================================\n"
)

cat(
  "FINAL SANITY CHECKS\n"
)

cat(
  "============================================\n\n"
)


cat(
  "Number of final genes:",
  nrow(final_evidence),
  "\n"
)


cat(
  "Unique genes:",
  dplyr::n_distinct(
    final_evidence$SYMBOL
  ),
  "\n"
)


cat(
  "GSE32863 validated:",
  sum(
    final_evidence$GSE32863_validated,
    na.rm = TRUE
  ),
  "/",
  nrow(final_evidence),
  "\n"
)


cat(
  "Strong LUAD dependencies:",
  sum(
    final_evidence$functional_dependency,
    na.rm = TRUE
  ),
  "\n"
)


cat(
  "Directionally LUAD-selective candidates:",
  sum(
    final_evidence$LUAD_selective,
    na.rm = TRUE
  ),
  "\n"
)


cat(
  "Statistically supported LUAD-selective candidates (FDR < 0.05):",
  sum(
    final_evidence$LUAD_selective_FDR05,
    na.rm = TRUE
  ),
  "\n"
)


cat(
  "\nDepMap coverage:\n"
)


print(
  final_evidence %>%
    dplyr::select(
      SYMBOL,
      n_models,
      n_models_non_luad,
      median_gene_effect,
      median_gene_effect_non_luad,
      LUAD_selectivity_difference,
      FDR
    ),
  n = Inf
)


# ============================================================
# 30. FINAL TARGET PRIORITIZATION
# ============================================================

cat(
  "\n============================================\n"
)

cat(
  "FINAL TARGET PRIORITIZATION\n"
)

cat(
  "============================================\n\n"
)


print(
  final_priority_table,
  n = Inf
)


# ============================================================
# 31. SAVE RESULTS
# ============================================================

saveRDS(
  final_evidence,
  file.path(
    processed_dir,
    "final_target_evidence_integrated.rds"
  )
)


saveRDS(
  final_priority_table,
  file.path(
    processed_dir,
    "final_target_priority_table.rds"
  )
)


write.csv(
  final_priority_table,
  file.path(
    processed_dir,
    "final_target_priority_table.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 32. SAVE DEPMAP COMPARISON SEPARATELY
# ============================================================

saveRDS(
  dependency_comparison,
  file.path(
    processed_dir,
    "DepMap_dependency_comparison.rds"
  )
)


write.csv(
  dependency_comparison,
  file.path(
    processed_dir,
    "DepMap_dependency_comparison.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 33. FINAL MESSAGE
# ============================================================

cat(
  "\n============================================\n"
)

cat(
  "RESULTS SAVED SUCCESSFULLY\n"
)

cat(
  "============================================\n\n"
)

cat(
  "Files created:\n"
)

cat(
  "- final_target_evidence_integrated.rds\n"
)

cat(
  "- final_target_priority_table.rds\n"
)

cat(
  "- final_target_priority_table.csv\n"
)

cat(
  "- DepMap_dependency_comparison.rds\n"
)

cat(
  "- DepMap_dependency_comparison.csv\n"
)

cat(
  "\nScript 12 completed successfully.\n"
)