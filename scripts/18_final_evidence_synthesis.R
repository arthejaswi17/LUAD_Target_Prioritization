# ============================================================
# SCRIPT 18: FINAL EVIDENCE SYNTHESIS AND REPORTING
# ============================================================
# Purpose:
#   Convert the validated Script 17 biological decision layer
#   into a report-ready evidence package.
#
# IMPORTANT:
#   Script 18 does not re-rank targets or alter Script 17 decisions.
# ============================================================

cat("\n============================================================\n")
cat("SCRIPT 18: FINAL EVIDENCE SYNTHESIS AND REPORTING\n")
cat("============================================================\n\n")

project_dir <- "C:/Users/arthe/Downloads/LUAD_Target_Prioritization/LUAD_Target_Prioritization"
processed_dir <- file.path(project_dir, "data", "processed")
results_dir <- file.path(project_dir, "results")
figures_dir <- file.path(project_dir, "figures")

if (!dir.exists(processed_dir)) stop("Processed directory not found: ", processed_dir)
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)
if (!dir.exists(figures_dir)) dir.create(figures_dir, recursive = TRUE)

required_packages <- c("dplyr", "ggplot2", "readr", "tibble", "tidyr")
missing <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing) > 0) {
  stop("Missing packages: ", paste(missing, collapse = ", "))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tibble)
  library(tidyr)
})


# ------------------------------------------------------------
# 1. LOAD SCRIPT 17 OUTPUT
# ------------------------------------------------------------

candidate_rds <- c(
  file.path(processed_dir, "Script17_final_biological_decision_table.rds"),
  file.path(processed_dir, "Script17_final_target_decision_table.rds"),
  file.path(processed_dir, "Script17_final_biological_interpretation.rds")
)

candidate_csv <- c(
  file.path(results_dir, "Script17_final_biological_decision_table.csv"),
  file.path(results_dir, "Script17_final_target_decision_table.csv"),
  file.path(processed_dir, "Script17_final_biological_decision_table.csv"),
  file.path(processed_dir, "Script17_final_target_decision_table.csv")
)

rds_found <- candidate_rds[file.exists(candidate_rds)]
csv_found <- candidate_csv[file.exists(candidate_csv)]

if (length(rds_found) > 0) {
  input_file <- rds_found[1]
  final_table <- readRDS(input_file)
  
} else if (length(csv_found) > 0) {
  input_file <- csv_found[1]
  final_table <- readr::read_csv(input_file, show_col_types = FALSE)
  
} else {
  stop("No Script 17 final decision table was found.")
}

if (!is.data.frame(final_table)) {
  stop("Script 17 output is not a data frame.")
}

cat("Using Script 17 output: ", input_file, "\n", sep = "")


# ------------------------------------------------------------
# 2. STANDARDIZE AND VALIDATE
# ------------------------------------------------------------

symbol_candidates <- c(
  "FINAL_SYMBOL",
  "SYMBOL",
  "gene_symbol",
  "symbol"
)

gene_col <- symbol_candidates[
  symbol_candidates %in% names(final_table)
][1]

if (is.na(gene_col)) {
  stop("Could not identify gene-symbol column.")
}

if (gene_col != "FINAL_SYMBOL") {
  
  final_table$FINAL_SYMBOL <- as.character(
    final_table[[gene_col]]
  )
  
} else {
  
  final_table$FINAL_SYMBOL <- as.character(
    final_table$FINAL_SYMBOL
  )
  
}


expected_targets <- c(
  "AURKA",
  "TOP2A",
  "CDK1",
  "PLK1",
  "CHEK1",
  "SGO1",
  "CDCA5",
  "SKA1",
  "BRCA1",
  "GMNN",
  "PRKDC",
  "DNA2"
)

actual_targets <- unique(
  final_table$FINAL_SYMBOL
)

if (length(actual_targets) != 12) {
  stop(
    "Expected 12 unique targets; found ",
    length(actual_targets),
    "."
  )
}

if (length(setdiff(expected_targets, actual_targets)) > 0) {
  stop(
    "Expected target(s) missing: ",
    paste(
      setdiff(expected_targets, actual_targets),
      collapse = ", "
    )
  )
}


# Script 17 has used two compatible output schemas across revisions.
# Standardize the validated Script 17 output here rather than requiring
# one exact column naming convention.

if (
  !"confidence" %in% names(final_table) &&
  "biological_confidence" %in% names(final_table)
) {
  
  final_table$confidence <- as.character(
    final_table$biological_confidence
  )
  
}


if (
  !"biological_interpretation" %in% names(final_table) &&
  "interpretation" %in% names(final_table)
) {
  
  final_table$biological_interpretation <- as.character(
    final_table$interpretation
  )
  
}


if (!"network_supported" %in% names(final_table)) {
  
  if ("degree" %in% names(final_table)) {
    
    final_table$network_supported <-
      !is.na(final_table$degree) &
      as.numeric(final_table$degree) > 0
    
  } else {
    
    stop(
      "Script 17 output lacks both network_supported and degree."
    )
    
  }
  
}


if (!"hallmark_supported" %in% names(final_table)) {
  
  if ("n_hallmark_pathways" %in% names(final_table)) {
    
    final_table$hallmark_supported <-
      !is.na(final_table$n_hallmark_pathways) &
      as.numeric(final_table$n_hallmark_pathways) > 0
    
  } else {
    
    stop(
      "Script 17 output lacks both hallmark_supported and n_hallmark_pathways."
    )
    
  }
  
}


if (!"evidence_domain_count" %in% names(final_table)) {
  
  evidence_cols <- c(
    "strong_TCGA_DE",
    "independent_validation",
    "strong_LUAD_dependency",
    "LUAD_selective_FDR05",
    "network_supported",
    "hallmark_supported"
  )
  
  missing_evidence <- setdiff(
    evidence_cols,
    names(final_table)
  )
  
  if (length(missing_evidence) > 0) {
    
    stop(
      "Script 17 output lacks evidence columns needed to calculate ",
      "evidence_domain_count: ",
      paste(
        missing_evidence,
        collapse = ", "
      )
    )
    
  }
  
  final_table$evidence_domain_count <- rowSums(
    sapply(
      final_table[evidence_cols],
      function(x) {
        !is.na(x) & as.logical(x)
      }
    )
  )
  
}


required_cols <- c(
  "FINAL_SYMBOL",
  "selection_group",
  "confidence",
  "biological_role",
  "biological_process",
  "evidence_domain_count",
  "strong_TCGA_DE",
  "very_strong_TCGA_DE",
  "independent_validation",
  "strong_LUAD_dependency",
  "LUAD_selective_FDR05",
  "network_supported",
  "hallmark_supported",
  "biological_interpretation"
)

missing_cols <- setdiff(
  required_cols,
  names(final_table)
)

if (length(missing_cols) > 0) {
  
  stop(
    "Script 17 output is missing required columns after schema ",
    "standardization: ",
    paste(
      missing_cols,
      collapse = ", "
    )
  )
  
}


primary_targets <- c(
  "AURKA",
  "TOP2A",
  "CDK1",
  "PLK1",
  "CHEK1"
)

secondary_targets <- c(
  "SGO1",
  "CDCA5",
  "SKA1"
)

supporting_targets <- c(
  "BRCA1",
  "GMNN",
  "PRKDC",
  "DNA2"
)


if ("final_biological_rank" %in% names(final_table)) {
  
  final_table$final_biological_rank <-
    as.numeric(
      final_table$final_biological_rank
    )
  
} else {
  
  final_table <- final_table %>%
    dplyr::arrange(
      match(
        FINAL_SYMBOL,
        expected_targets
      )
    ) %>%
    dplyr::mutate(
      final_biological_rank = dplyr::row_number()
    )
  
}


# ------------------------------------------------------------
# 3. REPORT-READY EVIDENCE TABLES
# ------------------------------------------------------------

evidence_summary <- final_table %>%
  dplyr::transmute(
    final_biological_rank,
    target = FINAL_SYMBOL,
    selection_group,
    confidence,
    biological_role,
    biological_process,
    evidence_domains = evidence_domain_count,
    TCGA_differential_expression = strong_TCGA_DE,
    TCGA_very_strong_expression = very_strong_TCGA_DE,
    independent_validation,
    strong_LUAD_dependency,
    LUAD_selective_dependency = LUAD_selective_FDR05,
    network_support = network_supported,
    Hallmark_support = hallmark_supported,
    biological_interpretation
  ) %>%
  dplyr::arrange(
    final_biological_rank
  )


primary_table <- evidence_summary %>%
  dplyr::filter(
    target %in% primary_targets
  ) %>%
  dplyr::arrange(
    final_biological_rank
  )


# Long-format matrix for plotting and downstream use.
evidence_matrix <- final_table %>%
  dplyr::select(
    FINAL_SYMBOL,
    selection_group,
    strong_TCGA_DE,
    very_strong_TCGA_DE,
    independent_validation,
    strong_LUAD_dependency,
    LUAD_selective_FDR05,
    network_supported,
    hallmark_supported
  ) %>%
  dplyr::rename(
    target = FINAL_SYMBOL
  )


project_summary <- tibble::tibble(
  
  metric = c(
    
    "Final target count",
    "Primary targets",
    "Strong secondary targets",
    "Supporting targets",
    "Targets with strong TCGA differential expression",
    "Targets with very strong TCGA differential expression",
    "Targets independently validated in GSE32863",
    "Targets with strong LUAD dependency",
    "Targets with FDR-supported LUAD selectivity",
    "Targets with network support",
    "Targets with Hallmark support"
    
  ),
  
  value = c(
    
    nrow(evidence_summary),
    
    sum(
      evidence_summary$selection_group ==
        "Primary target",
      na.rm = TRUE
    ),
    
    sum(
      evidence_summary$selection_group ==
        "Strong secondary target",
      na.rm = TRUE
    ),
    
    sum(
      evidence_summary$selection_group ==
        "Supporting target",
      na.rm = TRUE
    ),
    
    sum(
      evidence_summary$TCGA_differential_expression,
      na.rm = TRUE
    ),
    
    sum(
      evidence_summary$TCGA_very_strong_expression,
      na.rm = TRUE
    ),
    
    sum(
      evidence_summary$independent_validation,
      na.rm = TRUE
    ),
    
    sum(
      evidence_summary$strong_LUAD_dependency,
      na.rm = TRUE
    ),
    
    sum(
      evidence_summary$LUAD_selective_dependency,
      na.rm = TRUE
    ),
    
    sum(
      evidence_summary$network_support,
      na.rm = TRUE
    ),
    
    sum(
      evidence_summary$Hallmark_support,
      na.rm = TRUE
    )
    
  )
  
)


limitations_table <- tibble::tibble(
  
  domain = c(
    "Transcriptomic analysis",
    "Independent validation",
    "DepMap dependency",
    "Dependency selectivity",
    "Network analysis",
    "Pathway analysis",
    "Biological prioritization",
    "Experimental validation"
  ),
  
  limitation = c(
    
    "TCGA-LUAD matched tumor-normal analysis establishes disease-associated transcriptional differences but does not by itself establish causality.",
    
    "GSE32863 provides independent expression validation, but external validation remains expression-based rather than functional.",
    
    "DepMap dependency data are derived from cancer cell-line models and may not fully reproduce tumor microenvironment or patient biology.",
    
    "Selectivity estimates depend on available model coverage; limited LUAD representation can reduce statistical power for individual genes.",
    
    "Network centrality reflects the selected interaction resource and network construction and should not be interpreted as proof of biological causality.",
    
    "Pathway enrichment reflects coordinated gene-set behavior and does not demonstrate that every prioritized target is independently responsible for pathway activation.",
    
    "The final target groups integrate multiple evidence domains but remain a computational prioritization framework rather than proof of therapeutic efficacy.",
    
    "Experimental assays are required to establish target perturbation effects, therapeutic selectivity and mechanistic causality."
    
  )
  
)


primary_rationale <- evidence_summary %>%
  
  dplyr::filter(
    target %in% primary_targets
  ) %>%
  
  dplyr::select(
    final_biological_rank,
    target,
    biological_role,
    biological_process,
    evidence_domains,
    TCGA_differential_expression,
    TCGA_very_strong_expression,
    independent_validation,
    strong_LUAD_dependency,
    LUAD_selective_dependency,
    network_support,
    Hallmark_support,
    confidence,
    biological_interpretation
  )


# ------------------------------------------------------------
# 4. SAVE TABLES
# ------------------------------------------------------------

readr::write_csv(
  evidence_summary,
  file.path(
    results_dir,
    "Script18_final_target_evidence_summary.csv"
  )
)

readr::write_csv(
  primary_table,
  file.path(
    results_dir,
    "Script18_primary_target_evidence_summary.csv"
  )
)

readr::write_csv(
  evidence_matrix,
  file.path(
    results_dir,
    "Script18_evidence_domain_matrix.csv"
  )
)

readr::write_csv(
  project_summary,
  file.path(
    results_dir,
    "Script18_project_summary.csv"
  )
)

readr::write_csv(
  limitations_table,
  file.path(
    results_dir,
    "Script18_project_limitations.csv"
  )
)

readr::write_csv(
  primary_rationale,
  file.path(
    results_dir,
    "Script18_primary_target_rationale.csv"
  )
)


saveRDS(
  evidence_summary,
  file.path(
    processed_dir,
    "Script18_final_target_evidence_summary.rds"
  )
)

saveRDS(
  primary_table,
  file.path(
    processed_dir,
    "Script18_primary_target_evidence_summary.rds"
  )
)

saveRDS(
  evidence_matrix,
  file.path(
    processed_dir,
    "Script18_evidence_domain_matrix.rds"
  )
)

saveRDS(
  project_summary,
  file.path(
    processed_dir,
    "Script18_project_summary.rds"
  )
)

saveRDS(
  limitations_table,
  file.path(
    processed_dir,
    "Script18_project_limitations.rds"
  )
)


# ------------------------------------------------------------
# 5. FIGURE 1 — PRIMARY TARGET EVIDENCE CONVERGENCE
# ------------------------------------------------------------

primary_plot_data <- primary_table %>%
  
  dplyr::select(
    target,
    TCGA_differential_expression,
    TCGA_very_strong_expression,
    independent_validation,
    strong_LUAD_dependency,
    LUAD_selective_dependency,
    network_support,
    Hallmark_support
  ) %>%
  
  tidyr::pivot_longer(
    cols = -target,
    names_to = "evidence_domain",
    values_to = "supported"
  ) %>%
  
  dplyr::mutate(
    
    supported = as.numeric(supported),
    
    evidence_domain = dplyr::recode(
      
      evidence_domain,
      
      TCGA_differential_expression =
        "TCGA DE",
      
      TCGA_very_strong_expression =
        "Very strong TCGA DE",
      
      independent_validation =
        "GSE32863 validation",
      
      strong_LUAD_dependency =
        "LUAD dependency",
      
      LUAD_selective_dependency =
        "LUAD selectivity",
      
      network_support =
        "Network support",
      
      Hallmark_support =
        "Hallmark support"
      
    )
    
  )


p1 <- ggplot2::ggplot(
  primary_plot_data,
  ggplot2::aes(
    x = evidence_domain,
    y = target,
    fill = supported
  )
) +
  
  ggplot2::geom_tile(
    color = "white",
    linewidth = 0.5
  ) +
  
  ggplot2::scale_fill_gradient(
    limits = c(0, 1),
    low = "grey92",
    high = "grey25",
    breaks = c(0, 1),
    labels = c(
      "Not supported",
      "Supported"
    )
  ) +
  
  ggplot2::labs(
    title = "Evidence convergence across primary LUAD targets",
    subtitle = "Evidence domains retained from the validated Script 17 decision layer",
    x = NULL,
    y = NULL,
    fill = NULL
  ) +
  
  ggplot2::theme_minimal(
    base_size = 11
  ) +
  
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 35,
      hjust = 1
    ),
    panel.grid = ggplot2::element_blank()
  )


ggplot2::ggsave(
  file.path(
    figures_dir,
    "Script18_primary_target_evidence_convergence.png"
  ),
  p1,
  width = 10,
  height = 5.5,
  dpi = 300
)

ggplot2::ggsave(
  file.path(
    figures_dir,
    "Script18_primary_target_evidence_convergence.pdf"
  ),
  p1,
  width = 10,
  height = 5.5
)


# ------------------------------------------------------------
# 6. FIGURE 2 — EVIDENCE DOMAIN COVERAGE BY TARGET
# ------------------------------------------------------------
# Uses ASCII "YES" instead of the Unicode checkmark.
# This avoids Windows graphics-device encoding warnings.

evidence_heatmap <- evidence_summary %>%
  
  dplyr::select(
    final_biological_rank,
    target,
    TCGA_differential_expression,
    TCGA_very_strong_expression,
    independent_validation,
    strong_LUAD_dependency,
    LUAD_selective_dependency,
    network_support,
    Hallmark_support
  ) %>%
  
  dplyr::arrange(
    final_biological_rank
  ) %>%
  
  tidyr::pivot_longer(
    cols = -c(
      final_biological_rank,
      target
    ),
    names_to = "evidence_domain",
    values_to = "supported"
  ) %>%
  
  dplyr::mutate(
    
    supported = as.numeric(
      !is.na(supported) &
        as.logical(supported)
    ),
    
    evidence_domain = dplyr::recode(
      
      evidence_domain,
      
      TCGA_differential_expression =
        "TCGA DE",
      
      TCGA_very_strong_expression =
        "Very strong TCGA DE",
      
      independent_validation =
        "GSE32863 validation",
      
      strong_LUAD_dependency =
        "LUAD dependency",
      
      LUAD_selective_dependency =
        "LUAD selectivity",
      
      network_support =
        "Network support",
      
      Hallmark_support =
        "Hallmark support"
      
    ),
    
    target = factor(
      target,
      levels = rev(
        evidence_summary$target[
          order(
            evidence_summary$final_biological_rank
          )
        ]
      )
    ),
    
    evidence_domain = factor(
      
      evidence_domain,
      
      levels = c(
        "TCGA DE",
        "Very strong TCGA DE",
        "GSE32863 validation",
        "LUAD dependency",
        "LUAD selectivity",
        "Network support",
        "Hallmark support"
      )
      
    )
    
  )


p2 <- ggplot2::ggplot(
  evidence_heatmap,
  ggplot2::aes(
    x = evidence_domain,
    y = target,
    fill = supported
  )
) +
  
  ggplot2::geom_tile(
    color = "white",
    linewidth = 0.6
  ) +
  
  ggplot2::scale_fill_gradient(
    limits = c(0, 1),
    low = "grey92",
    high = "grey30",
    breaks = c(0, 1),
    labels = c(
      "Not supported",
      "Supported"
    )
  ) +
  
  ggplot2::labs(
    title = "Evidence-domain coverage across final LUAD targets",
    subtitle = "Target-level support across the seven evidence domains retained from Script 17",
    x = NULL,
    y = NULL,
    fill = NULL
  ) +
  
  ggplot2::theme_minimal(
    base_size = 11
  ) +
  
  ggplot2::theme(
    
    axis.text.x = ggplot2::element_text(
      angle = 35,
      hjust = 1
    ),
    
    axis.text.y = ggplot2::element_text(
      face = "bold"
    ),
    
    panel.grid = ggplot2::element_blank(),
    
    legend.position = "bottom"
    
  )


ggplot2::ggsave(
  file.path(
    figures_dir,
    "Script18_evidence_domain_coverage.png"
  ),
  p2,
  width = 11,
  height = 7,
  dpi = 300
)

ggplot2::ggsave(
  file.path(
    figures_dir,
    "Script18_evidence_domain_coverage.pdf"
  ),
  p2,
  width = 11,
  height = 7
)


# ===========================================================
# FIGURE 3: FINAL TARGET-GROUP COMPOSITION
# ===========================================================

group_composition <- tibble::tibble(
  selection_group = factor(
    c(
      "Primary targets",
      "Strong secondary targets",
      "Supporting targets"
    ),
    levels = c(
      "Primary targets",
      "Strong secondary targets",
      "Supporting targets"
    )
  ),
  target_count = c(
    5,
    3,
    4
  )
)

p3 <- ggplot2::ggplot(
  group_composition,
  ggplot2::aes(
    x = selection_group,
    y = target_count
  )
) +
  
  ggplot2::geom_col(
    width = 0.65
  ) +
  
  ggplot2::geom_text(
    ggplot2::aes(
      label = target_count
    ),
    vjust = -0.5,
    size = 5,
    fontface = "bold"
  ) +
  
  ggplot2::scale_y_continuous(
    breaks = 0:5,
    limits = c(0, 5.8),
    expand = c(0, 0)
  ) +
  
  ggplot2::labs(
    title = "Final LUAD target-group composition",
    subtitle = "Number of targets retained in each evidence-based selection group",
    x = NULL,
    y = "Number of targets"
  ) +
  
  ggplot2::theme_minimal(
    base_size = 11
  ) +
  
  ggplot2::theme(
    panel.grid.major.x =
      ggplot2::element_blank(),
    
    panel.grid.minor =
      ggplot2::element_blank(),
    
    axis.text.x =
      ggplot2::element_text(
        size = 10,
        face = "bold"
      ),
    
    axis.title.y =
      ggplot2::element_text(
        size = 10
      ),
    
    plot.title =
      ggplot2::element_text(
        face = "bold",
        size = 15
      ),
    
    plot.subtitle =
      ggplot2::element_text(
        size = 9.5
      ),
    
    plot.margin =
      ggplot2::margin(
        10,
        15,
        10,
        15
      )
  )

# -----------------------------------------------------------
# Save Figure 3
# -----------------------------------------------------------

ggplot2::ggsave(
  filename = file.path(
    figures_dir,
    "Script18_final_target_group_composition.png"
  ),
  plot = p3,
  width = 8,
  height = 6,
  dpi = 300,
  bg = "white"
)

ggplot2::ggsave(
  filename = file.path(
    figures_dir,
    "Script18_final_target_group_composition.pdf"
  ),
  plot = p3,
  width = 8,
  height = 6,
  device = "pdf",
  bg = "white"
)

cat(
  "\nFigure 3 saved successfully:\n",
  "  - Script18_final_target_group_composition.png\n",
  "  - Script18_final_target_group_composition.pdf\n"
)

# ------------------------------------------------------------
# 8. FINAL QUALITY CONTROL
# ------------------------------------------------------------

if (nrow(evidence_summary) != 12) {
  stop(
    "QC failed: expected 12 final targets."
  )
}

if (length(unique(evidence_summary$target)) != 12) {
  stop(
    "QC failed: target symbols are not unique."
  )
}

if (
  any(
    is.na(evidence_summary$target) |
    evidence_summary$target == ""
  )
) {
  
  stop(
    "QC failed: missing target symbol."
  )
  
}

if (
  any(
    is.na(evidence_summary$biological_role) |
    evidence_summary$biological_role == ""
  )
) {
  
  stop(
    "QC failed: missing biological role."
  )
  
}

if (
  any(
    is.na(evidence_summary$biological_interpretation) |
    evidence_summary$biological_interpretation == ""
  )
) {
  
  stop(
    "QC failed: missing biological interpretation."
  )
  
}

if (
  !all(
    primary_targets %in%
    evidence_summary$target
  )
) {
  
  stop(
    "QC failed: primary targets incomplete."
  )
  
}

if (
  !all(
    secondary_targets %in%
    evidence_summary$target
  )
) {
  
  stop(
    "QC failed: secondary targets incomplete."
  )
  
}

if (
  !all(
    supporting_targets %in%
    evidence_summary$target
  )
) {
  
  stop(
    "QC failed: supporting targets incomplete."
  )
  
}


cat("\n============================================================\n")
cat("SCRIPT 18 FINAL QUALITY CONTROL: PASSED\n")
cat("============================================================\n\n")

cat("Final targets: 12\n")

cat(
  "Primary: AURKA, TOP2A, CDK1, PLK1, CHEK1\n"
)

cat(
  "Strong secondary: SGO1, CDCA5, SKA1\n"
)

cat(
  "Supporting: BRCA1, GMNN, PRKDC, DNA2\n\n"
)

cat(
  "Output tables written to: ",
  results_dir,
  "\n",
  sep = ""
)

cat(
  "Processed RDS files written to: ",
  processed_dir,
  "\n",
  sep = ""
)

cat(
  "Figures written to: ",
  figures_dir,
  "\n",
  sep = ""
)

cat(
  "\n============================================================\n"
)

cat(
  "SCRIPT 18 COMPLETED SUCCESSFULLY\n"
)

cat(
  "============================================================\n"
)