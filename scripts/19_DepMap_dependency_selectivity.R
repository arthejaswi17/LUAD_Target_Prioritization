# ============================================================
# SCRIPT 19: DEPMAP DEPENDENCY AND LUAD SELECTIVITY FIGURE
# ============================================================
#
# Purpose:
#   Create a report-ready figure showing:
#
#   Panel A:
#     Median LUAD Gene Effect for the 12 final targets
#
#   Panel B:
#     LUAD selectivity difference:
#     median non-LUAD Gene Effect -
#     median LUAD Gene Effect
#
# Data source:
#   data/processed/DepMap_dependency_comparison.rds
#
# Outputs:
#   figures/Script19_DepMap_dependency_selectivity.pdf
#   figures/Script19_DepMap_dependency_selectivity.png
#   results/Script19_DepMap_dependency_figure_data.csv
#
# ============================================================


# ============================================================
# 1. REQUIRED PACKAGES
# ============================================================

required_packages <- c(
  "dplyr",
  "ggplot2",
  "readr",
  "grid"
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
    "The following packages are missing:\n",
    paste(
      missing_packages,
      collapse = ", "
    ),
    "\n\nInstall them using:\n",
    "install.packages(c(",
    paste0(
      '"',
      missing_packages,
      '"',
      collapse = ", "
    ),
    "))"
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(grid)
})


# ============================================================
# 2. PROJECT DIRECTORIES
# ============================================================

project_root <- normalizePath(
  getwd(),
  winslash = "/",
  mustWork = TRUE
)

processed_dir <- file.path(
  project_root,
  "data",
  "processed"
)

figures_dir <- file.path(
  project_root,
  "figures"
)

results_dir <- file.path(
  project_root,
  "results"
)

if (!dir.exists(processed_dir)) {
  stop(
    "The processed-data directory was not found:\n",
    processed_dir,
    "\n\nPlease run this script from:\n",
    "LUAD_Target_Prioritization/LUAD_Target_Prioritization"
  )
}

dir.create(
  figures_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  results_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 3. LOAD EXISTING DEPMAP RESULTS
# ============================================================

depmap_file <- file.path(
  processed_dir,
  "DepMap_dependency_comparison.rds"
)

if (!file.exists(depmap_file)) {
  stop(
    "The required DepMap file was not found:\n",
    depmap_file,
    "\n\nExpected filename:\n",
    "DepMap_dependency_comparison.rds"
  )
}

depmap_raw <- readRDS(
  depmap_file
)

depmap_raw <- as.data.frame(
  depmap_raw,
  stringsAsFactors = FALSE
)

cat("\n============================================================\n")
cat("SCRIPT 19: DEPMAP DEPENDENCY AND SELECTIVITY FIGURE\n")
cat("============================================================\n\n")

cat(
  "Project directory:\n",
  project_root,
  "\n\n"
)

cat(
  "DepMap input file:\n",
  depmap_file,
  "\n\n"
)

cat(
  "Available columns in the DepMap file:\n"
)

print(
  names(depmap_raw)
)

cat("\n")


# ============================================================
# 4. FINAL TARGET ORDER
# ============================================================

# Established final biological ranking from Scripts 16 and 17.

final_rank <- c(
  AURKA = 1,
  TOP2A = 2,
  CDK1  = 3,
  PLK1  = 4,
  CHEK1 = 5,
  SGO1  = 6,
  CDCA5 = 7,
  SKA1  = 8,
  BRCA1 = 9,
  GMNN  = 10,
  PRKDC = 11,
  DNA2  = 12
)

final_targets <- names(
  final_rank
)


# ============================================================
# 5. HELPER FUNCTION FOR COLUMN DETECTION
# ============================================================

find_column <- function(
    data,
    possible_names,
    required = TRUE
) {
  
  exact_match <- possible_names[
    possible_names %in% names(data)
  ]
  
  if (length(exact_match) > 0) {
    return(
      exact_match[1]
    )
  }
  
  lower_names <- tolower(
    names(data)
  )
  
  lower_possible <- tolower(
    possible_names
  )
  
  fuzzy_match <- names(data)[
    lower_names %in% lower_possible
  ]
  
  if (length(fuzzy_match) > 0) {
    return(
      fuzzy_match[1]
    )
  }
  
  if (required) {
    stop(
      "Could not identify a required column.\n\n",
      "Possible column names checked:\n",
      paste(
        possible_names,
        collapse = ", "
      ),
      "\n\nAvailable columns:\n",
      paste(
        names(data),
        collapse = ", "
      )
    )
  }
  
  return(NULL)
}


# ============================================================
# 6. IDENTIFY REQUIRED COLUMNS
# ============================================================

symbol_column <- find_column(
  depmap_raw,
  c(
    "SYMBOL",
    "symbol",
    "gene_symbol",
    "gene",
    "Gene",
    "Hugo_Symbol"
  )
)

luad_n_column <- find_column(
  depmap_raw,
  c(
    "n_models",
    "n_luad_models",
    "LUAD_n_models",
    "n_LUAD",
    "luad_models"
  )
)

non_luad_n_column <- find_column(
  depmap_raw,
  c(
    "n_models_non_luad",
    "n_non_luad_models",
    "non_luad_n_models",
    "n_non_LUAD",
    "non_LUAD_models"
  )
)

luad_effect_column <- find_column(
  depmap_raw,
  c(
    "median_gene_effect",
    "median_LUAD_gene_effect",
    "LUAD_median_gene_effect",
    "median_luad_effect",
    "luad_median_effect"
  )
)

non_luad_effect_column <- find_column(
  depmap_raw,
  c(
    "median_gene_effect_non_luad",
    "median_non_luad_gene_effect",
    "non_LUAD_median_gene_effect",
    "median_non_luad_effect",
    "non_luad_median_effect"
  )
)

selectivity_column <- find_column(
  depmap_raw,
  c(
    "LUAD_selectivity_difference",
    "luad_selectivity_difference",
    "selectivity_difference",
    "LUAD_selectivity",
    "luad_selectivity"
  ),
  required = FALSE
)

fdr_column <- find_column(
  depmap_raw,
  c(
    "FDR",
    "fdr",
    "adjusted_p_value",
    "padj",
    "p_adjusted",
    "wilcoxon_fdr"
  ),
  required = FALSE
)

strong_dependency_column <- find_column(
  depmap_raw,
  c(
    "strong_LUAD_dependency",
    "strong_luad_dependency",
    "strong_dependency",
    "LUAD_dependency",
    "is_strong_dependency"
  ),
  required = FALSE
)


# ============================================================
# 7. STANDARDIZE COLUMN NAMES
# ============================================================

depmap <- depmap_raw %>%
  mutate(
    SYMBOL = as.character(
      .data[[symbol_column]]
    ),
    
    n_models = suppressWarnings(
      as.numeric(
        .data[[luad_n_column]]
      )
    ),
    
    n_models_non_luad = suppressWarnings(
      as.numeric(
        .data[[non_luad_n_column]]
      )
    ),
    
    median_gene_effect = suppressWarnings(
      as.numeric(
        .data[[luad_effect_column]]
      )
    ),
    
    median_gene_effect_non_luad = suppressWarnings(
      as.numeric(
        .data[[non_luad_effect_column]]
      )
    )
  )


# ============================================================
# 7A. SELECTIVITY DIFFERENCE
# ============================================================

if (!is.null(selectivity_column)) {
  
  depmap$LUAD_selectivity_difference <-
    suppressWarnings(
      as.numeric(
        depmap_raw[[selectivity_column]]
      )
    )
  
} else {
  
  depmap$LUAD_selectivity_difference <-
    depmap$median_gene_effect_non_luad -
    depmap$median_gene_effect
  
}


# ============================================================
# 7B. FDR
# ============================================================

if (!is.null(fdr_column)) {
  
  depmap$FDR <-
    suppressWarnings(
      as.numeric(
        depmap_raw[[fdr_column]]
      )
    )
  
} else {
  
  depmap$FDR <- NA_real_
  
}


# ============================================================
# 7C. STRONG DEPENDENCY FLAG
# ============================================================

if (!is.null(strong_dependency_column)) {
  
  raw_flag <-
    depmap_raw[[strong_dependency_column]]
  
  depmap$strong_dependency <-
    tolower(
      as.character(
        raw_flag
      )
    ) %in% c(
      "true",
      "t",
      "yes",
      "y",
      "1"
    )
  
} else {
  
  depmap$strong_dependency <-
    !is.na(
      depmap$n_models
    ) &
    depmap$n_models >= 10 &
    !is.na(
      depmap$median_gene_effect
    ) &
    depmap$median_gene_effect <= -1
  
}


# ============================================================
# 8. FILTER THE 12 FINAL TARGETS
# ============================================================

depmap_final <- depmap %>%
  filter(
    SYMBOL %in% final_targets
  ) %>%
  distinct(
    SYMBOL,
    .keep_all = TRUE
  ) %>%
  mutate(
    final_rank = unname(
      final_rank[SYMBOL]
    ),
    
    FDR_supported_selectivity =
      !is.na(FDR) &
      !is.na(n_models) &
      !is.na(n_models_non_luad) &
      n_models >= 10 &
      n_models_non_luad >= 10 &
      !is.na(
        LUAD_selectivity_difference
      ) &
      LUAD_selectivity_difference > 0 &
      FDR < 0.05
  ) %>%
  arrange(
    final_rank
  )


missing_targets <- setdiff(
  final_targets,
  depmap_final$SYMBOL
)

if (length(missing_targets) > 0) {
  stop(
    "The following final targets were not found in the DepMap file:\n",
    paste(
      missing_targets,
      collapse = ", "
    )
  )
}

if (nrow(depmap_final) != 12) {
  stop(
    "Expected 12 final targets, but found ",
    nrow(depmap_final),
    "."
  )
}


# ============================================================
# 9. QUALITY CONTROL
# ============================================================

cat(
  "Final targets plotted:",
  nrow(depmap_final),
  "\n"
)

cat(
  "Strong LUAD dependencies:",
  sum(
    depmap_final$strong_dependency,
    na.rm = TRUE
  ),
  "/12\n"
)

cat(
  "FDR-supported LUAD-selective targets:",
  sum(
    depmap_final$FDR_supported_selectivity,
    na.rm = TRUE
  ),
  "/12\n\n"
)

cat(
  "Final DepMap table:\n\n"
)

print(
  depmap_final %>%
    select(
      final_rank,
      SYMBOL,
      n_models,
      n_models_non_luad,
      median_gene_effect,
      median_gene_effect_non_luad,
      LUAD_selectivity_difference,
      FDR,
      strong_dependency,
      FDR_supported_selectivity
    ),
  row.names = FALSE
)


# ============================================================
# 10. SAVE FIGURE DATA
# ============================================================

figure_data <- depmap_final %>%
  select(
    final_rank,
    SYMBOL,
    n_models,
    n_models_non_luad,
    median_gene_effect,
    median_gene_effect_non_luad,
    LUAD_selectivity_difference,
    FDR,
    strong_dependency,
    FDR_supported_selectivity
  )

figure_data_file <- file.path(
  results_dir,
  "Script19_DepMap_dependency_figure_data.csv"
)

write_csv(
  figure_data,
  figure_data_file
)


# ============================================================
# 11. PREPARE PLOT DATA
# ============================================================

plot_df <- figure_data %>%
  mutate(
    target = factor(
      SYMBOL,
      levels = rev(
        final_targets
      )
    )
  )


# ============================================================
# 12. PANEL A - LUAD FUNCTIONAL DEPENDENCY
# ============================================================

dependency_plot <- ggplot(
  plot_df,
  aes(
    x = median_gene_effect,
    y = target
  )
) +
  
  # Strong-dependency threshold
  geom_vline(
    xintercept = -1,
    linetype = "dashed",
    linewidth = 0.7
  ) +
  
  # Horizontal line from zero to observed dependency
  geom_segment(
    aes(
      x = 0,
      xend = median_gene_effect,
      y = target,
      yend = target
    ),
    linewidth = 0.8
  ) +
  
  # Open = not strong dependency
  # Filled = strong dependency
  geom_point(
    aes(
      shape = strong_dependency
    ),
    size = 3.5
  ) +
  
  # Only show model counts where coverage is limited.
  # This highlights PRKDC and DNA2 without cluttering
  # the ten genes with 53 LUAD models.
  geom_text(
    data = subset(
      plot_df,
      n_models < 10
    ),
    aes(
      x = median_gene_effect - 0.10,
      label = paste0(
        "n = ",
        n_models
      )
    ),
    hjust = 1,
    size = 3.0
  ) +
  
  scale_shape_manual(
    values = c(
      `FALSE` = 1,
      `TRUE` = 16
    ),
    labels = c(
      `FALSE` = "Not strong dependency",
      `TRUE` = "Strong dependency"
    ),
    name = NULL
  ) +
  
  scale_x_continuous(
    name = "Median LUAD Gene Effect",
    expand = expansion(
      mult = c(
        0.08,
        0.08
      )
    )
  ) +
  
  labs(
    title =
      "A. LUAD functional dependency",
    
    subtitle =
      "More negative Gene Effect values indicate stronger dependency",
    
    y = NULL
  ) +
  
  theme_classic(
    base_size = 11
  ) +
  
  theme(
    
    plot.title =
      element_text(
        face = "bold",
        size = 13,
        margin =
          ggplot2::margin(
            b = 5
          )
      ),
    
    plot.subtitle =
      element_text(
        size = 9.5,
        margin =
          ggplot2::margin(
            b = 10
          )
      ),
    
    axis.text.y =
      element_text(
        face = "bold",
        size = 10
      ),
    
    axis.title.x =
      element_text(
        size = 10,
        margin =
          ggplot2::margin(
            t = 8
          )
      ),
    
    legend.position = "bottom",
    
    legend.margin =
      ggplot2::margin(
        t = 5
      ),
    
    legend.text =
      element_text(
        size = 8.5
      ),
    
    panel.grid.major.y =
      element_blank(),
    
    panel.grid.minor =
      element_blank(),
    
    plot.margin =
      ggplot2::margin(
        t = 8,
        r = 12,
        b = 8,
        l = 8
      )
  )


# ============================================================
# 13. PANEL B - LUAD SELECTIVITY
# ============================================================

selectivity_plot <- ggplot(
  plot_df,
  aes(
    x = LUAD_selectivity_difference,
    y = target
  )
) +
  
  # Zero = no difference between groups
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.7
  ) +
  
  # Horizontal line from zero to observed selectivity
  geom_segment(
    aes(
      x = 0,
      xend = LUAD_selectivity_difference,
      y = target,
      yend = target
    ),
    linewidth = 0.8
  ) +
  
  # Open = not FDR-supported
  # Filled = FDR-supported
  geom_point(
    aes(
      shape = FDR_supported_selectivity
    ),
    size = 3.5
  ) +
  
  # FDR labels only for supported targets
  geom_text(
    data = subset(
      plot_df,
      FDR_supported_selectivity
    ),
    aes(
      x =
        LUAD_selectivity_difference + 0.025,
      
      label =
        "FDR < 0.05"
    ),
    hjust = 0,
    size = 3.0
  ) +
  
  scale_shape_manual(
    values = c(
      `FALSE` = 1,
      `TRUE` = 16
    ),
    labels = c(
      `FALSE` = "Not FDR-supported",
      `TRUE` = "FDR-supported selectivity"
    ),
    name = NULL
  ) +
  
  scale_x_continuous(
    name =
      "LUAD selectivity difference",
    
    expand =
      expansion(
        mult = c(
          0.08,
          0.18
        )
      )
  ) +
  
  labs(
    title =
      "B. LUAD-selective dependency",
    
    subtitle =
      "Positive values indicate stronger dependency in LUAD than non-LUAD models",
    
    y = NULL
  ) +
  
  theme_classic(
    base_size = 11
  ) +
  
  theme(
    
    plot.title =
      element_text(
        face = "bold",
        size = 13,
        margin =
          ggplot2::margin(
            b = 5
          )
      ),
    
    plot.subtitle =
      element_text(
        size = 9.5,
        margin =
          ggplot2::margin(
            b = 10
          )
      ),
    
    axis.text.y =
      element_text(
        face = "bold",
        size = 10
      ),
    
    axis.title.x =
      element_text(
        size = 10,
        margin =
          ggplot2::margin(
            t = 8
          )
      ),
    
    legend.position = "bottom",
    
    legend.margin =
      ggplot2::margin(
        t = 5
      ),
    
    legend.text =
      element_text(
        size = 8.5
      ),
    
    panel.grid.major.y =
      element_blank(),
    
    panel.grid.minor =
      element_blank(),
    
    plot.margin =
      ggplot2::margin(
        t = 8,
        r = 12,
        b = 8,
        l = 8
      )
  )


# ============================================================
# 14. OUTPUT FILES
# ============================================================

pdf_file <- file.path(
  figures_dir,
  "Script19_DepMap_dependency_selectivity.pdf"
)

png_file <- file.path(
  figures_dir,
  "Script19_DepMap_dependency_selectivity.png"
)


# ============================================================
# 15. COMBINED FIGURE LAYOUT
# ============================================================

# Only the title and the two plots are included.
# No bottom footnotes are used, preventing clipping.

combined_layout <- grid.layout(
  nrow = 2,
  ncol = 2,
  
  heights = unit(
    c(
      0.10,
      0.90
    ),
    "null"
  ),
  
  widths = unit(
    c(
      1,
      1
    ),
    "null"
  )
)


# ============================================================
# 16. SAVE COMBINED PDF
# ============================================================

pdf(
  pdf_file,
  width = 14,
  height = 9,
  useDingbats = FALSE
)

grid.newpage()

pushViewport(
  viewport(
    layout = combined_layout
  )
)


# ------------------------------------------------------------
# Overall title
# ------------------------------------------------------------

grid.text(
  "Functional Dependency and LUAD Selectivity of Final Targets",
  
  x = 0.5,
  y = 0.5,
  
  vp = viewport(
    layout.pos.row = 1,
    layout.pos.col = 1:2
  ),
  
  gp = gpar(
    fontface = "bold",
    fontsize = 17
  )
)


# ------------------------------------------------------------
# Panel A
# ------------------------------------------------------------

print(
  dependency_plot,
  
  vp = viewport(
    layout.pos.row = 2,
    layout.pos.col = 1
  )
)


# ------------------------------------------------------------
# Panel B
# ------------------------------------------------------------

print(
  selectivity_plot,
  
  vp = viewport(
    layout.pos.row = 2,
    layout.pos.col = 2
  )
)

popViewport()

dev.off()


# ============================================================
# 17. SAVE COMBINED PNG
# ============================================================

png(
  png_file,
  width = 4200,
  height = 2700,
  res = 300
)

grid.newpage()

pushViewport(
  viewport(
    layout = combined_layout
  )
)


# ------------------------------------------------------------
# Overall title
# ------------------------------------------------------------

grid.text(
  "Functional Dependency and LUAD Selectivity of Final Targets",
  
  x = 0.5,
  y = 0.5,
  
  vp = viewport(
    layout.pos.row = 1,
    layout.pos.col = 1:2
  ),
  
  gp = gpar(
    fontface = "bold",
    fontsize = 48
  )
)


# ------------------------------------------------------------
# Panel A
# ------------------------------------------------------------

print(
  dependency_plot,
  
  vp = viewport(
    layout.pos.row = 2,
    layout.pos.col = 1
  )
)


# ------------------------------------------------------------
# Panel B
# ------------------------------------------------------------

print(
  selectivity_plot,
  
  vp = viewport(
    layout.pos.row = 2,
    layout.pos.col = 2
  )
)

popViewport()

dev.off()


# ============================================================
# 18. FINAL OUTPUT CHECKS
# ============================================================

required_outputs <- c(
  pdf_file,
  png_file,
  figure_data_file
)

missing_outputs <- required_outputs[
  !file.exists(
    required_outputs
  )
]

if (length(missing_outputs) > 0) {
  
  stop(
    "The following expected output files were not created:\n",
    paste(
      missing_outputs,
      collapse = "\n"
    )
  )
  
}


# ============================================================
# 19. FINAL MESSAGE
# ============================================================

cat(
  "\n============================================================\n"
)

cat(
  "SCRIPT 19 COMPLETED SUCCESSFULLY\n"
)

cat(
  "============================================================\n\n"
)

cat(
  "Created files:\n\n"
)

cat(
  pdf_file,
  "\n"
)

cat(
  png_file,
  "\n"
)

cat(
  figure_data_file,
  "\n\n"
)

cat(
  "Strong LUAD dependencies: ",
  sum(
    depmap_final$strong_dependency,
    na.rm = TRUE
  ),
  "/12\n",
  sep = ""
)

cat(
  "FDR-supported LUAD-selective targets: ",
  sum(
    depmap_final$FDR_supported_selectivity,
    na.rm = TRUE
  ),
  "/12\n",
  sep = ""
)

cat(
  "\n============================================================\n"
)