# SCRIPT 14 FUNCTIONAL & NETWORK CHARACTERIZATION

# Characterize the 12 final target genes identified using Script 12.

# Inputs:
#   data/processed/final_target_priority_table.csv

# Main analyses:
#   1. Final target validation
#   2. Network/module characterization
#   3. Hallmark pathway enrichment
#   4. Target-pathway overlap analysis
#   5. Visualization


# 1. INITIALIZATION

cat("\n")
cat("SCRIPT 14: FUNCTIONAL & NETWORK CHARACTERIZATION\n")

options(stringsAsFactors = FALSE)
set.seed(1234)


# 2. PACKAGE CHECK

required_packages <- c(
  "dplyr",
  "readr",
  "tibble",
  "tidyr",
  "ggplot2",
  "msigdbr"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {
  
  stop(
    paste0(
      "Missing required R packages: ",
      paste(missing_packages, collapse = ", "),
      "\nPlease install them before running Script 14."
    )
  )
}

cat("Required packages detected.\n")


# 3. DIRECTORY SETUP

processed_dir <- "data/processed"
figures_dir <- file.path(processed_dir, "figures")

if (!dir.exists(processed_dir)) {
  dir.create(
    processed_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

if (!dir.exists(figures_dir)) {
  dir.create(
    figures_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

cat("Directories ready.\n")


# 4. LOAD SCRIPT 12 RESULTS

cat("\n")
cat("Loading Script 12 final target table...\n")

input_file <- file.path(
  processed_dir,
  "final_target_priority_table.csv"
)

if (!file.exists(input_file)) {
  
  stop(
    paste0(
      "Required file not found:\n",
      input_file,
      "\n\nRun Script 12 first or place ",
      "final_target_priority_table.csv in ",
      processed_dir,
      "."
    )
  )
}

targets <- readr::read_csv(
  input_file,
  show_col_types = FALSE,
  progress = FALSE
)

cat("Rows loaded:", nrow(targets), "\n")
cat("Columns loaded:", ncol(targets), "\n")


# 5. VALIDATE INPUT

required_columns <- c(
  "SYMBOL",
  "TCGA_shrunken_log2FC",
  "TCGA_FDR",
  "GSE32863_log2FC",
  "GSE32863_FDR",
  "GSE32863_validated",
  "degree",
  "module",
  "n_pathways",
  "n_models",
  "n_models_non_luad",
  "median_gene_effect",
  "median_gene_effect_non_luad",
  "LUAD_selectivity_difference",
  "p_value",
  "FDR",
  "strong_LUAD_dependency",
  "LUAD_selective",
  "LUAD_selective_FDR05",
  "transcriptomic_evidence",
  "independent_validation",
  "core_evidence_count",
  "priority_tier"
)

missing_columns <- setdiff(
  required_columns,
  colnames(targets)
)

if (length(missing_columns) > 0) {
  
  stop(
    paste0(
      "Missing required columns:\n",
      paste(missing_columns, collapse = ", ")
    )
  )
}

cat("Required columns: PASSED\n")


# 6. CLEAN TARGET TABLE

targets <- targets |>
  dplyr::mutate(
    SYMBOL = trimws(as.character(.data$SYMBOL)),
    degree = suppressWarnings(
      as.numeric(as.character(.data$degree))
    ),
    module = as.character(.data$module),
    n_pathways = suppressWarnings(
      as.numeric(as.character(.data$n_pathways))
    ),
    TCGA_shrunken_log2FC = suppressWarnings(
      as.numeric(as.character(.data$TCGA_shrunken_log2FC))
    ),
    GSE32863_log2FC = suppressWarnings(
      as.numeric(as.character(.data$GSE32863_log2FC))
    ),
    core_evidence_count = suppressWarnings(
      as.numeric(as.character(.data$core_evidence_count))
    )
  ) |>
  dplyr::filter(
    !is.na(.data$SYMBOL),
    .data$SYMBOL != ""
  ) |>
  dplyr::distinct(
    .data$SYMBOL,
    .keep_all = TRUE
  )

cat("Unique final targets:", nrow(targets), "\n")

if (nrow(targets) != 12) {
  
  warning(
    paste0(
      "Expected 12 final targets, but found ",
      nrow(targets),
      "."
    )
  )
}


# 7. FINAL TARGET LIST

final_genes <- unique(targets$SYMBOL)

cat("\n")
cat("Final target genes:\n")
cat(paste(final_genes, collapse = ", "), "\n")


# 8. NETWORK CHARACTERIZATION

cat("\n")
cat("NETWORK CHARACTERIZATION\n")
network_summary <- targets |>
  dplyr::select(
    .data$SYMBOL,
    .data$degree,
    .data$module,
    .data$n_pathways,
    .data$core_evidence_count,
    .data$priority_tier
  ) |>
  dplyr::arrange(
    dplyr::desc(.data$degree)
  )

valid_degree <- network_summary$degree[
  !is.na(network_summary$degree)
]

if (length(valid_degree) > 0) {
  
  cat(
    "Network degree range:",
    min(valid_degree),
    "to",
    max(valid_degree),
    "\n"
  )
  
} else {
  
  cat("Network degree range: unavailable\n")
}

cat(
  "Network modules represented:",
  length(unique(network_summary$module)),
  "\n"
)

cat("\nNetwork summary:\n")

# Print a controlled number of columns/rows.
print(
  as.data.frame(network_summary),
  row.names = FALSE
)


# 9. MODULE COMPOSITION

module_summary <- targets |>
  dplyr::count(
    .data$module,
    name = "n_targets"
  ) |>
  dplyr::arrange(
    dplyr::desc(.data$n_targets),
    .data$module
  )

cat("\n")
cat("Module composition:\n")

print(
  as.data.frame(module_summary),
  row.names = FALSE
)


# 10. HALLMARK GENE SETS

cat("\n")
cat("LOADING MSIGDB HALLMARK GENE SETS\n")

hallmark_raw <- msigdbr::msigdbr(
  species = "Homo sapiens",
  collection = "H"
)

cat(
  "Raw Hallmark rows:",
  nrow(hallmark_raw),
  "\n"
)

# Select only columns needed for this analysis.
hallmark <- hallmark_raw |>
  dplyr::select(
    .data$gene_symbol,
    .data$gs_name
  ) |>
  dplyr::mutate(
    gene_symbol = trimws(
      as.character(.data$gene_symbol)
    ),
    gs_name = as.character(.data$gs_name)
  ) |>
  dplyr::filter(
    !is.na(.data$gene_symbol),
    .data$gene_symbol != "",
    !is.na(.data$gs_name),
    .data$gs_name != ""
  ) |>
  dplyr::distinct()

cat(
  "Hallmark gene-set rows:",
  nrow(hallmark),
  "\n"
)

cat(
  "Hallmark pathways:",
  length(unique(hallmark$gs_name)),
  "\n"
)


# 11. TARGET-HALLMARK OVERLAP

cat("\n")
cat("Calculating final-target / Hallmark overlaps...\n")

target_set <- unique(final_genes)

hallmark_universe <- unique(
  hallmark$gene_symbol
)

target_overlap <- intersect(
  target_set,
  hallmark_universe
)

cat(
  "Final targets represented in Hallmark:",
  length(target_overlap),
  "/",
  length(target_set),
  "\n"
)

if (length(target_overlap) == 0) {
  
  stop(
    "None of the final target genes could be matched to MSigDB Hallmark gene sets."
  )
}


# 12. HALLMARK ENRICHMENT

cat("\n")
cat("Running Hallmark enrichment...\n")

n_target <- length(target_overlap)
n_universe <- length(hallmark_universe)

# Pathway sizes

pathway_sizes <- hallmark |>
  dplyr::group_by(.data$gs_name) |>
  dplyr::summarise(
    pathway_size = dplyr::n_distinct(
      .data$gene_symbol
    ),
    .groups = "drop"
  )


# Target overlap per pathway

overlap_results <- hallmark |>
  dplyr::filter(
    .data$gene_symbol %in% target_overlap
  ) |>
  dplyr::group_by(.data$gs_name) |>
  dplyr::summarise(
    overlap = dplyr::n_distinct(
      .data$gene_symbol
    ),
    target_genes = paste(
      sort(unique(.data$gene_symbol)),
      collapse = ", "
    ),
    .groups = "drop"
  )

hallmark_enrichment <- pathway_sizes |>
  dplyr::left_join(
    overlap_results,
    by = "gs_name"
  ) |>
  dplyr::mutate(
    overlap = ifelse(
      is.na(.data$overlap),
      0,
      .data$overlap
    ),
    target_genes = ifelse(
      is.na(.data$target_genes),
      "",
      .data$target_genes
    )
  )


hallmark_enrichment$p_value <- stats::phyper(
  q = hallmark_enrichment$overlap - 1,
  m = hallmark_enrichment$pathway_size,
  n = n_universe - hallmark_enrichment$pathway_size,
  k = n_target,
  lower.tail = FALSE
)

hallmark_enrichment$FDR <- stats::p.adjust(
  hallmark_enrichment$p_value,
  method = "BH"
)


# Order results

hallmark_enrichment <- hallmark_enrichment |>
  dplyr::arrange(
    .data$FDR,
    dplyr::desc(.data$overlap)
  ) |>
  dplyr::filter(
    .data$overlap > 0
  )

cat(
  "Hallmark pathways containing at least one final target:",
  nrow(hallmark_enrichment),
  "\n"
)

cat(
  "Significant pathways (FDR < 0.05):",
  sum(
    hallmark_enrichment$FDR < 0.05,
    na.rm = TRUE
  ),
  "\n"
)


# 13. DISPLAY ENRICHMENT RESULTS

cat("\n")
cat("Top Hallmark pathways containing final targets:\n")

top_hallmark <- hallmark_enrichment |>
  dplyr::select(
    .data$gs_name,
    .data$overlap,
    .data$pathway_size,
    .data$p_value,
    .data$FDR,
    .data$target_genes
  ) |>
  dplyr::slice_head(n = 20)

print(
  as.data.frame(top_hallmark),
  row.names = FALSE
)


# 14. SIGNIFICANT HALLMARK PATHWAYS

significant_hallmark <- hallmark_enrichment |>
  dplyr::filter(
    !is.na(.data$FDR),
    .data$FDR < 0.05
  ) |>
  dplyr::arrange(
    .data$FDR,
    dplyr::desc(.data$overlap)
  )

cat("\n")
cat("Significant Hallmark pathways:\n")

if (nrow(significant_hallmark) > 0) {
  
  print(
    as.data.frame(
      significant_hallmark |>
        dplyr::select(
          .data$gs_name,
          .data$overlap,
          .data$p_value,
          .data$FDR,
          .data$target_genes
        )
    ),
    row.names = FALSE
  )
  
} else {
  
  cat(
    "No Hallmark pathway reached FDR < 0.05.\n"
  )
}


# 15. TARGET-PATHWAY MATRIX

cat("\n")
cat("Creating target × Hallmark pathway matrix...\n")

relevant_pathways <- unique(
  hallmark_enrichment$gs_name
)

target_pathway_long <- hallmark |>
  dplyr::filter(
    .data$gene_symbol %in% target_overlap,
    .data$gs_name %in% relevant_pathways
  ) |>
  dplyr::distinct(
    .data$gene_symbol,
    .data$gs_name
  )

if (nrow(target_pathway_long) > 0) {
  
  target_pathway_matrix <- target_pathway_long |>
    dplyr::mutate(
      present = 1L
    ) |>
    tidyr::pivot_wider(
      names_from = .data$gs_name,
      values_from = .data$present,
      values_fill = 0
    )
  
} else {
  
  target_pathway_matrix <- tibble::tibble(
    gene_symbol = character(0)
  )
}

matrix_pathway_count <- max(
  ncol(target_pathway_matrix) - 1,
  0
)

cat(
  "Target-pathway matrix dimensions:",
  nrow(target_pathway_matrix),
  "×",
  matrix_pathway_count,
  "\n"
)


# 16. PATHWAY COVERAGE PER TARGET

target_pathway_summary <- target_pathway_long |>
  dplyr::group_by(.data$gene_symbol) |>
  dplyr::summarise(
    n_hallmark_pathways = dplyr::n_distinct(
      .data$gs_name
    ),
    hallmark_pathways = paste(
      sort(unique(.data$gs_name)),
      collapse = "; "
    ),
    .groups = "drop"
  )

targets_functional <- targets |>
  dplyr::left_join(
    target_pathway_summary,
    by = c("SYMBOL" = "gene_symbol")
  ) |>
  dplyr::mutate(
    n_hallmark_pathways = ifelse(
      is.na(.data$n_hallmark_pathways),
      0,
      .data$n_hallmark_pathways
    ),
    hallmark_pathways = ifelse(
      is.na(.data$hallmark_pathways),
      "",
      .data$hallmark_pathways
    )
  )

cat("\n")
cat("Hallmark pathway coverage per target:\n")

coverage_display <- targets_functional |>
  dplyr::select(
    .data$SYMBOL,
    .data$degree,
    .data$module,
    .data$n_pathways,
    .data$n_hallmark_pathways,
    .data$hallmark_pathways
  ) |>
  dplyr::arrange(
    dplyr::desc(.data$n_hallmark_pathways),
    dplyr::desc(.data$degree)
  )

print(
  as.data.frame(coverage_display),
  row.names = FALSE
)


# 17. NETWORK CENTRALITY FIGURE

cat("\n")
cat("Creating network centrality figure...\n")

network_plot_data <- targets |>
  dplyr::arrange(.data$degree) |>
  dplyr::mutate(
    SYMBOL = factor(
      .data$SYMBOL,
      levels = .data$SYMBOL
    )
  )

p_network <- ggplot2::ggplot(
  network_plot_data,
  ggplot2::aes(
    x = .data$SYMBOL,
    y = .data$degree
  )
) +
  ggplot2::geom_col() +
  ggplot2::coord_flip() +
  ggplot2::labs(
    title = "Network Centrality of Final LUAD Targets",
    x = "Target gene",
    y = "Network degree"
  ) +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      face = "bold"
    )
  )

ggplot2::ggsave(
  filename = file.path(
    figures_dir,
    "Script14_network_centrality.png"
  ),
  plot = p_network,
  width = 9,
  height = 6,
  dpi = 300
)

ggplot2::ggsave(
  filename = file.path(
    figures_dir,
    "Script14_network_centrality.pdf"
  ),
  plot = p_network,
  width = 9,
  height = 6
)


# 18. HALLMARK ENRICHMENT FIGURE

cat("Creating Hallmark enrichment figure...\n")

plot_enrichment <- hallmark_enrichment |>
  dplyr::arrange(
    .data$FDR,
    dplyr::desc(.data$overlap)
  ) |>
  dplyr::slice_head(n = 15) |>
  dplyr::mutate(
    pathway_label = gsub(
      "^HALLMARK_",
      "",
      .data$gs_name
    ),
    pathway_label = gsub(
      "_",
      " ",
      .data$pathway_label
    ),
    pathway_label = factor(
      .data$pathway_label,
      levels = rev(
        unique(.data$pathway_label)
      )
    ),
    neg_log10_FDR = -log10(
      pmax(
        .data$FDR,
        .Machine$double.xmin
      )
    )
  )

if (nrow(plot_enrichment) > 0) {
  
  p_enrichment <- ggplot2::ggplot(
    plot_enrichment,
    ggplot2::aes(
      x = .data$overlap,
      y = .data$pathway_label,
      size = .data$neg_log10_FDR
    )
  ) +
    ggplot2::geom_point() +
    ggplot2::labs(
      title = "Hallmark Pathway Representation Among Final Targets",
      x = "Number of final targets",
      y = "Hallmark pathway",
      size = "-log10(FDR)"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold"
      )
    )
  
  ggplot2::ggsave(
    filename = file.path(
      figures_dir,
      "Script14_Hallmark_enrichment.png"
    ),
    plot = p_enrichment,
    width = 10,
    height = 7,
    dpi = 300
  )
  
  ggplot2::ggsave(
    filename = file.path(
      figures_dir,
      "Script14_Hallmark_enrichment.pdf"
    ),
    plot = p_enrichment,
    width = 10,
    height = 7
  )
  
} else {
  
  cat(
    "No pathway-overlap data available for enrichment figure.\n"
  )
}


# 19. TARGET-PATHWAY HEATMAP

cat("Creating target-pathway heatmap...\n")

heatmap_data <- target_pathway_long |>
  dplyr::mutate(
    gene_symbol = factor(
      .data$gene_symbol,
      levels = final_genes
    )
  )

if (nrow(heatmap_data) > 0) {
  
  p_heatmap <- ggplot2::ggplot(
    heatmap_data,
    ggplot2::aes(
      x = .data$gs_name,
      y = .data$gene_symbol
    )
  ) +
    ggplot2::geom_tile() +
    ggplot2::labs(
      title = "Final Targets Across Hallmark Pathways",
      x = "Hallmark pathway",
      y = "Target gene"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5,
        size = 7
      ),
      axis.text.y = ggplot2::element_text(
        size = 9
      ),
      plot.title = ggplot2::element_text(
        face = "bold"
      )
    )
  
  ggplot2::ggsave(
    filename = file.path(
      figures_dir,
      "Script14_target_pathway_heatmap.png"
    ),
    plot = p_heatmap,
    width = 14,
    height = 7,
    dpi = 300
  )
  
  ggplot2::ggsave(
    filename = file.path(
      figures_dir,
      "Script14_target_pathway_heatmap.pdf"
    ),
    plot = p_heatmap,
    width = 14,
    height = 7
  )
  
} else {
  
  cat(
    "No target-pathway relationships available for heatmap.\n"
  )
}


# 20. INTEGRATED FUNCTIONAL SUMMARY

cat("\n")
cat("Creating integrated functional summary...\n")

functional_columns <- c(
  "SYMBOL",
  "priority_tier",
  "TCGA_shrunken_log2FC",
  "TCGA_FDR",
  "GSE32863_log2FC",
  "GSE32863_FDR",
  "GSE32863_validated",
  "degree",
  "module",
  "n_pathways",
  "n_hallmark_pathways",
  "hallmark_pathways",
  "n_models",
  "n_models_non_luad",
  "median_gene_effect",
  "median_gene_effect_non_luad",
  "LUAD_selectivity_difference",
  "p_value",
  "FDR",
  "strong_LUAD_dependency",
  "LUAD_selective",
  "LUAD_selective_FDR05",
  "core_evidence_count"
)

functional_summary <- targets_functional |>
  dplyr::select(
    dplyr::all_of(functional_columns)
  ) |>
  dplyr::arrange(
    .data$priority_tier,
    dplyr::desc(.data$core_evidence_count),
    dplyr::desc(.data$degree)
  )


# 21. SAVE RESULTS

cat("\n")
cat("Saving Script 14 outputs...\n")

saveRDS(
  functional_summary,
  file.path(
    processed_dir,
    "final_target_functional_summary.rds"
  )
)

readr::write_csv(
  functional_summary,
  file.path(
    processed_dir,
    "final_target_functional_summary.csv"
  )
)

saveRDS(
  hallmark_enrichment,
  file.path(
    processed_dir,
    "final_target_Hallmark_enrichment.rds"
  )
)

readr::write_csv(
  hallmark_enrichment,
  file.path(
    processed_dir,
    "final_target_Hallmark_enrichment.csv"
  )
)

saveRDS(
  significant_hallmark,
  file.path(
    processed_dir,
    "final_target_Hallmark_significant.rds"
  )
)

readr::write_csv(
  significant_hallmark,
  file.path(
    processed_dir,
    "final_target_Hallmark_significant.csv"
  )
)

saveRDS(
  target_pathway_matrix,
  file.path(
    processed_dir,
    "final_target_Hallmark_matrix.rds"
  )
)

readr::write_csv(
  target_pathway_matrix,
  file.path(
    processed_dir,
    "final_target_Hallmark_matrix.csv"
  )
)

saveRDS(
  module_summary,
  file.path(
    processed_dir,
    "final_target_module_summary.rds"
  )
)

readr::write_csv(
  module_summary,
  file.path(
    processed_dir,
    "final_target_module_summary.csv"
  )
)


# 22. FINAL QUALITY CHECKS

cat("\n")
cat("FINAL QUALITY CHECKS\n")

cat(
  "Final target count:",
  nrow(targets),
  "\n"
)

cat(
  "Unique target count:",
  length(unique(targets$SYMBOL)),
  "\n"
)

cat(
  "Targets represented in Hallmark:",
  length(target_overlap),
  "/",
  length(target_set),
  "\n"
)

cat(
  "Hallmark pathways tested:",
  length(unique(hallmark$gs_name)),
  "\n"
)

cat(
  "Hallmark pathways with target overlap:",
  nrow(hallmark_enrichment),
  "\n"
)

cat(
  "Significant Hallmark pathways:",
  nrow(significant_hallmark),
  "\n"
)

cat(
  "Network modules represented:",
  length(unique(targets$module)),
  "\n"
)


# Target uniqueness

if (anyDuplicated(targets$SYMBOL) > 0) {
  
  stop(
    "Duplicate target symbols detected."
  )
  
} else {
  
  cat(
    "Target uniqueness: PASSED\n"
  )
}


# Network degree completeness

if (anyNA(targets$degree)) {
  
  warning(
    "Some targets have missing network degree values."
  )
  
} else {
  
  cat(
    "Network degree completeness: PASSED\n"
  )
}


# Hallmark enrichment completeness

if (
  anyNA(hallmark_enrichment$p_value) ||
  anyNA(hallmark_enrichment$FDR)
) {
  
  warning(
    "Some Hallmark enrichment values are missing."
  )
  
} else {
  
  cat(
    "Hallmark enrichment completeness: PASSED\n"
  )
}


# Output existence checks

expected_output_files <- c(
  
  file.path(
    processed_dir,
    "final_target_functional_summary.rds"
  ),
  
  file.path(
    processed_dir,
    "final_target_functional_summary.csv"
  ),
  
  file.path(
    processed_dir,
    "final_target_Hallmark_enrichment.rds"
  ),
  
  file.path(
    processed_dir,
    "final_target_Hallmark_enrichment.csv"
  ),
  
  file.path(
    processed_dir,
    "final_target_Hallmark_significant.rds"
  ),
  
  file.path(
    processed_dir,
    "final_target_Hallmark_significant.csv"
  ),
  
  file.path(
    processed_dir,
    "final_target_Hallmark_matrix.rds"
  ),
  
  file.path(
    processed_dir,
    "final_target_Hallmark_matrix.csv"
  ),
  
  file.path(
    processed_dir,
    "final_target_module_summary.rds"
  ),
  
  file.path(
    processed_dir,
    "final_target_module_summary.csv"
  )
)

output_status <- file.exists(
  expected_output_files
)

if (all(output_status)) {
  
  cat(
    "Core output files: PASSED\n"
  )
  
} else {
  
  missing_outputs <- expected_output_files[
    !output_status
  ]
  
  warning(
    paste0(
      "Some expected output files were not created:\n",
      paste(
        missing_outputs,
        collapse = "\n"
      )
    )
  )
}


# 23. COMPLETION

cat("\n")
cat("SCRIPT 14 COMPLETED SUCCESSFULLY\n")

cat("Outputs saved:\n")
cat("- final_target_functional_summary.rds\n")
cat("- final_target_functional_summary.csv\n")
cat("- final_target_Hallmark_enrichment.rds\n")
cat("- final_target_Hallmark_enrichment.csv\n")
cat("- final_target_Hallmark_significant.rds\n")
cat("- final_target_Hallmark_significant.csv\n")
cat("- final_target_Hallmark_matrix.rds\n")
cat("- final_target_Hallmark_matrix.csv\n")
cat("- final_target_module_summary.rds\n")
cat("- final_target_module_summary.csv\n")
cat("- figures/Script14_network_centrality.png\n")
cat("- figures/Script14_network_centrality.pdf\n")
cat("- figures/Script14_Hallmark_enrichment.png\n")
cat("- figures/Script14_Hallmark_enrichment.pdf\n")
cat("- figures/Script14_target_pathway_heatmap.png\n")
cat("- figures/Script14_target_pathway_heatmap.pdf\n")

cat("\n")
cat("Next stage: Final target biological interpretation and selection\n")