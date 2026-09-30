# ============================================================
# SCRIPT 20: FINAL PROJECT AUDIT
# ============================================================
#
# Purpose:
#   Final non-analytical audit of the completed LUAD project.
#
# IMPORTANT:
#   This script does NOT perform new analysis.
#   It does NOT change rankings.
#   It does NOT regenerate figures.
#   It only checks the existing project outputs.
#
# ============================================================


# ============================================================
# 1. PROJECT DIRECTORIES
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

scripts_dir <- file.path(
  project_root,
  "scripts"
)


cat("\n============================================================\n")
cat("SCRIPT 20: FINAL PROJECT AUDIT\n")
cat("============================================================\n\n")

cat(
  "Project root:\n",
  project_root,
  "\n\n"
)


# ============================================================
# 2. AUDIT STORAGE
# ============================================================

audit_results <- data.frame(
  category = character(),
  check = character(),
  status = character(),
  details = character(),
  stringsAsFactors = FALSE
)


add_check <- function(
    category,
    check,
    status,
    details
) {
  
  audit_results <<- rbind(
    audit_results,
    data.frame(
      category = category,
      check = check,
      status = status,
      details = details,
      stringsAsFactors = FALSE
    )
  )
  
}


# ============================================================
# 3. FINAL TARGET DEFINITION
# ============================================================

final_targets <- c(
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


expected_rank <- setNames(
  seq_along(final_targets),
  final_targets
)


# ============================================================
# 4. DIRECTORY CHECK
# ============================================================

required_directories <- c(
  processed = processed_dir,
  figures = figures_dir,
  results = results_dir,
  scripts = scripts_dir
)


for (directory_name in names(required_directories)) {
  
  directory_path <- required_directories[[directory_name]]
  
  if (dir.exists(directory_path)) {
    
    add_check(
      "Project structure",
      paste(
        "Directory exists:",
        directory_name
      ),
      "PASS",
      directory_path
    )
    
  } else {
    
    add_check(
      "Project structure",
      paste(
        "Directory exists:",
        directory_name
      ),
      "FAIL",
      directory_path
    )
    
  }
  
}


# ============================================================
# 5. SCRIPT INVENTORY
# ============================================================
#
# We only check that the major analysis stages are represented.
# Exact filenames are NOT treated as critical because several
# scripts were renamed during development.
# ============================================================

script_files <- if (dir.exists(scripts_dir)) {
  list.files(
    scripts_dir,
    pattern = "\\.R$",
    full.names = FALSE
  )
} else {
  character()
}


if (length(script_files) > 0) {
  
  add_check(
    "Scripts",
    "R analysis scripts present",
    "PASS",
    paste(
      length(script_files),
      "R scripts found."
    )
  )
  
} else {
  
  add_check(
    "Scripts",
    "R analysis scripts present",
    "FAIL",
    "No R scripts were found."
  )
  
}


# Check the final scripts specifically.

final_script_patterns <- c(
  "16.*ranking",
  "17.*[Bb]iological",
  "18.*[Ee]vidence",
  "19.*DepMap",
  "20.*audit"
)


for (pattern in final_script_patterns) {
  
  matches <- grep(
    pattern,
    script_files,
    value = TRUE
  )
  
  if (length(matches) > 0) {
    
    add_check(
      "Scripts",
      paste(
        "Final-stage script:",
        pattern
      ),
      "PASS",
      paste(
        matches,
        collapse = ", "
      )
    )
    
  } else {
    
    add_check(
      "Scripts",
      paste(
        "Final-stage script:",
        pattern
      ),
      "WARN",
      "No matching filename found."
    )
    
  }
  
}


# ============================================================
# 6. CORE INPUT / OUTPUT FILE DISCOVERY
# ============================================================

find_files <- function(
    directory,
    patterns
) {
  
  if (!dir.exists(directory)) {
    return(character())
  }
  
  all_files <- list.files(
    directory,
    full.names = TRUE,
    recursive = FALSE
  )
  
  matches <- character()
  
  for (pattern in patterns) {
    
    matches <- c(
      matches,
      all_files[
        grepl(
          pattern,
          basename(all_files),
          ignore.case = TRUE
        )
      ]
    )
    
  }
  
  unique(matches)
  
}


# ============================================================
# 7. DEPMAP OUTPUT CHECK
# ============================================================

depmap_matches <- find_files(
  processed_dir,
  c(
    "^DepMap_dependency_comparison\\.rds$"
  )
)


if (length(depmap_matches) > 0) {
  
  add_check(
    "Core data",
    "DepMap dependency dataset",
    "PASS",
    basename(depmap_matches[1])
  )
  
} else {
  
  add_check(
    "Core data",
    "DepMap dependency dataset",
    "FAIL",
    "DepMap_dependency_comparison.rds not found."
  )
  
}


# ============================================================
# 8. SCRIPT 17 DECISION TABLE DISCOVERY
# ============================================================

script17_matches <- find_files(
  processed_dir,
  c(
    "^Script17_final_biological_decision_table\\.rds$",
    "Script17.*decision.*\\.rds"
  )
)


script17_table <- NULL


if (length(script17_matches) > 0) {
  
  script17_file <- script17_matches[1]
  
  script17_table <- tryCatch(
    readRDS(script17_file),
    error = function(e) NULL
  )
  
  
  if (!is.null(script17_table)) {
    
    script17_table <- as.data.frame(
      script17_table,
      stringsAsFactors = FALSE
    )
    
    add_check(
      "Script 17",
      "Final decision table readable",
      "PASS",
      paste(
        basename(script17_file),
        "-",
        nrow(script17_table),
        "rows"
      )
    )
    
  } else {
    
    add_check(
      "Script 17",
      "Final decision table readable",
      "FAIL",
      "RDS could not be read."
    )
    
  }
  
} else {
  
  add_check(
    "Script 17",
    "Final decision table readable",
    "FAIL",
    "No Script 17 decision-table RDS found."
  )
  
}


# ============================================================
# 9. IDENTIFY SCRIPT 17 SYMBOL COLUMN
# ============================================================

script17_symbol_column <- NULL


if (!is.null(script17_table)) {
  
  candidate_columns <- c(
    "FINAL_SYMBOL",
    "SYMBOL",
    "gene_symbol",
    "symbol",
    "target"
  )
  
  available_columns <- candidate_columns[
    candidate_columns %in%
      names(script17_table)
  ]
  
  if (length(available_columns) > 0) {
    
    script17_symbol_column <- available_columns[1]
    
    add_check(
      "Script 17",
      "Target symbol column",
      "PASS",
      script17_symbol_column
    )
    
  } else {
    
    add_check(
      "Script 17",
      "Target symbol column",
      "FAIL",
      paste(
        "Available columns:",
        paste(
          names(script17_table),
          collapse = ", "
        )
      )
    )
    
  }
  
}


# ============================================================
# 10. SCRIPT 17 FINAL TARGET MEMBERSHIP
# ============================================================

if (
  !is.null(script17_table) &&
  !is.null(script17_symbol_column)
) {
  
  observed_targets <- unique(
    as.character(
      script17_table[[script17_symbol_column]]
    )
  )
  
  observed_targets <- observed_targets[
    !is.na(observed_targets) &
      nzchar(observed_targets)
  ]
  
  
  missing_targets <- setdiff(
    final_targets,
    observed_targets
  )
  
  extra_targets <- setdiff(
    observed_targets,
    final_targets
  )
  
  
  if (
    length(missing_targets) == 0 &&
    length(extra_targets) == 0
  ) {
    
    add_check(
      "Script 17",
      "Final target membership",
      "PASS",
      "Exactly the established 12 final targets."
    )
    
  } else {
    
    add_check(
      "Script 17",
      "Final target membership",
      "FAIL",
      paste(
        "Missing:",
        paste(
          missing_targets,
          collapse = ", "
        ),
        "| Extra:",
        paste(
          extra_targets,
          collapse = ", "
        )
      )
    )
    
  }
  
}


# ============================================================
# 11. TARGET GROUP CHECK
# ============================================================

if (
  length(final_targets) == 12 &&
  length(primary_targets) == 5 &&
  length(secondary_targets) == 3 &&
  length(supporting_targets) == 4 &&
  length(
    unique(
      c(
        primary_targets,
        secondary_targets,
        supporting_targets
      )
    )
  ) == 12
) {
  
  add_check(
    "Final targets",
    "Final target group composition",
    "PASS",
    "Primary = 5; secondary = 3; supporting = 4."
  )
  
} else {
  
  add_check(
    "Final targets",
    "Final target group composition",
    "FAIL",
    "Target groups are inconsistent."
  )
  
}


# ============================================================
# 12. SCRIPT 17 COLUMN AUDIT
# ============================================================
#
# Only require fields that are actually necessary for the final
# decision layer. Do NOT require columns merely because an older
# version of Script 18 used them.
# ============================================================

required_decision_columns <- c(
  "selection_group",
  "biological_role",
  "biological_process"
)


if (!is.null(script17_table)) {
  
  missing_decision_columns <- setdiff(
    required_decision_columns,
    names(script17_table)
  )
  
  
  if (
    length(missing_decision_columns) == 0
  ) {
    
    add_check(
      "Script 17",
      "Core biological decision columns",
      "PASS",
      "Selection group, biological role and biological process are present."
    )
    
  } else {
    
    add_check(
      "Script 17",
      "Core biological decision columns",
      "FAIL",
      paste(
        "Missing:",
        paste(
          missing_decision_columns,
          collapse = ", "
        )
      )
    )
    
  }
  
}


# ============================================================
# 13. SCRIPT 17 GROUP VALUES
# ============================================================

if (
  !is.null(script17_table) &&
  "selection_group" %in%
  names(script17_table)
) {
  
  group_values <- unique(
    as.character(
      script17_table$selection_group
    )
  )
  
  group_values <- group_values[
    !is.na(group_values)
  ]
  
  
  expected_groups <- c(
    "Primary",
    "Strong secondary",
    "Supporting"
  )
  
  
  invalid_groups <- setdiff(
    group_values,
    expected_groups
  )
  
  
  if (length(invalid_groups) == 0) {
    
    add_check(
      "Script 17",
      "Selection-group labels",
      "PASS",
      "All selection groups use the established categories."
    )
    
  } else {
    
    add_check(
      "Script 17",
      "Selection-group labels",
      "FAIL",
      paste(
        "Unexpected groups:",
        paste(
          invalid_groups,
          collapse = ", "
        )
      )
    )
    
  }
  
}


# ============================================================
# 14. SCRIPT 17 BIOLOGICAL INFORMATION COMPLETENESS
# ============================================================

if (
  !is.null(script17_table) &&
  !is.null(script17_symbol_column)
) {
  
  final_rows <- script17_table[
    as.character(
      script17_table[[script17_symbol_column]]
    ) %in% final_targets,
    ,
    drop = FALSE
  ]
  
  
  biological_columns <- intersect(
    c(
      "selection_group",
      "biological_role",
      "biological_process"
    ),
    names(final_rows)
  )
  
  
  missing_biological_values <- character()
  
  
  for (column_name in biological_columns) {
    
    if (
      any(
        is.na(
          final_rows[[column_name]]
        )
      )
    ) {
      
      missing_biological_values <- c(
        missing_biological_values,
        column_name
      )
      
    }
    
  }
  
  
  if (
    length(missing_biological_values) == 0
  ) {
    
    add_check(
      "Script 17",
      "Final biological information completeness",
      "PASS",
      "No missing core biological annotations."
    )
    
  } else {
    
    add_check(
      "Script 17",
      "Final biological information completeness",
      "WARN",
      paste(
        "Missing values in:",
        paste(
          missing_biological_values,
          collapse = ", "
        )
      )
    )
    
  }
  
}


# ============================================================
# 15. DEPMAP FIGURE DATA DISCOVERY
# ============================================================

depmap_figure_matches <- find_files(
  results_dir,
  c(
    "Script19.*DepMap.*figure.*data.*\\.csv",
    "DepMap.*figure.*data.*\\.csv"
  )
)


script19_table <- NULL


if (length(depmap_figure_matches) > 0) {
  
  depmap_figure_file <- depmap_figure_matches[1]
  
  script19_table <- tryCatch(
    read.csv(
      depmap_figure_file,
      stringsAsFactors = FALSE
    ),
    error = function(e) NULL
  )
  
  
  if (!is.null(script19_table)) {
    
    add_check(
      "Script 19",
      "DepMap figure data readable",
      "PASS",
      paste(
        basename(depmap_figure_file),
        "-",
        nrow(script19_table),
        "rows"
      )
    )
    
  } else {
    
    add_check(
      "Script 19",
      "DepMap figure data readable",
      "FAIL",
      "DepMap figure CSV could not be read."
    )
    
  }
  
} else {
  
  add_check(
    "Script 19",
    "DepMap figure data export",
    "WARN",
    "CSV export not found; figure may still be valid."
  )
  
}


# ============================================================
# 16. SCRIPT 19 DEPENDENCY COUNTS
# ============================================================

if (!is.null(script19_table)) {
  
  required_depmap_columns <- c(
    "SYMBOL",
    "n_models",
    "n_models_non_luad",
    "median_gene_effect",
    "median_gene_effect_non_luad",
    "LUAD_selectivity_difference",
    "FDR"
  )
  
  
  missing_depmap_columns <- setdiff(
    required_depmap_columns,
    names(script19_table)
  )
  
  
  if (
    length(missing_depmap_columns) == 0
  ) {
    
    add_check(
      "Script 19",
      "Core DepMap columns",
      "PASS",
      "All core numerical fields are present."
    )
    
  } else {
    
    add_check(
      "Script 19",
      "Core DepMap columns",
      "FAIL",
      paste(
        "Missing:",
        paste(
          missing_depmap_columns,
          collapse = ", "
        )
      )
    )
    
  }
  
  
  # ----------------------------------------------------------
  # Target membership
  # ----------------------------------------------------------
  
  depmap_targets <- unique(
    as.character(
      script19_table$SYMBOL
    )
  )
  
  
  missing_depmap_targets <- setdiff(
    final_targets,
    depmap_targets
  )
  
  
  extra_depmap_targets <- setdiff(
    depmap_targets,
    final_targets
  )
  
  
  if (
    length(missing_depmap_targets) == 0 &&
    length(extra_depmap_targets) == 0 &&
    length(depmap_targets) == 12
  ) {
    
    add_check(
      "Script 19",
      "DepMap target membership",
      "PASS",
      "Exactly the established 12 final targets."
    )
    
  } else {
    
    add_check(
      "Script 19",
      "DepMap target membership",
      "FAIL",
      paste(
        "Missing:",
        paste(
          missing_depmap_targets,
          collapse = ", "
        ),
        "| Extra:",
        paste(
          extra_depmap_targets,
          collapse = ", "
        )
      )
    )
    
  }
  
  
  # ----------------------------------------------------------
  # Strong dependency count
  # ----------------------------------------------------------
  
  if (
    "strong_dependency" %in%
    names(script19_table)
  ) {
    
    strong_count <- sum(
      script19_table$strong_dependency,
      na.rm = TRUE
    )
    
    if (strong_count == 7) {
      
      add_check(
        "Script 19",
        "Strong LUAD dependency count",
        "PASS",
        "7/12 targets."
      )
      
    } else {
      
      add_check(
        "Script 19",
        "Strong LUAD dependency count",
        "FAIL",
        paste(
          "Expected 7; found",
          strong_count
        )
      )
      
    }
    
  }
  
  
  # ----------------------------------------------------------
  # Selectivity count
  # ----------------------------------------------------------
  
  if (
    "FDR_supported_selectivity" %in%
    names(script19_table)
  ) {
    
    selective_count <- sum(
      script19_table$FDR_supported_selectivity,
      na.rm = TRUE
    )
    
    if (selective_count == 5) {
      
      add_check(
        "Script 19",
        "FDR-supported LUAD selectivity count",
        "PASS",
        "5/12 targets."
      )
      
    } else {
      
      add_check(
        "Script 19",
        "FDR-supported LUAD selectivity count",
        "FAIL",
        paste(
          "Expected 5; found",
          selective_count
        )
      )
      
    }
    
  }
  
  
  # ----------------------------------------------------------
  # Limited model coverage
  # ----------------------------------------------------------
  
  limited_models <- script19_table$SYMBOL[
    script19_table$n_models < 10
  ]
  
  
  if (
    setequal(
      limited_models,
      c(
        "PRKDC",
        "DNA2"
      )
    )
  ) {
    
    add_check(
      "Script 19",
      "Limited LUAD model coverage",
      "PASS",
      "PRKDC and DNA2 have only 2 LUAD models."
    )
    
  } else {
    
    add_check(
      "Script 19",
      "Limited LUAD model coverage",
      "WARN",
      paste(
        "Targets with <10 LUAD models:",
        paste(
          limited_models,
          collapse = ", "
        )
      )
    )
    
  }
  
}


# ============================================================
# 17. EXPECTED FINAL FIGURES
# ============================================================
#
# These are checked using flexible filename matching.
# ============================================================

figure_patterns <- c(
  
  "Script17.*final.*biological.*ranking.*\\.pdf",
  
  "Script17.*evidence.*domain.*heatmap.*\\.pdf",
  
  "Script17.*primary.*target.*evidence.*\\.pdf",
  
  "Script18.*evidence.*domain.*coverage.*\\.pdf",
  
  "Script18.*final.*target.*group.*composition.*\\.pdf",
  
  "Script19.*DepMap.*dependency.*selectivity.*\\.pdf",
  
  "Script19.*DepMap.*dependency.*selectivity.*\\.png"
  
)


for (pattern in figure_patterns) {
  
  matches <- find_files(
    figures_dir,
    pattern
  )
  
  
  if (length(matches) > 0) {
    
    add_check(
      "Figures",
      pattern,
      "PASS",
      basename(matches[1])
    )
    
  } else {
    
    add_check(
      "Figures",
      pattern,
      "WARN",
      "No matching figure found."
    )
    
  }
  
}


# ============================================================
# 18. OUTPUT INVENTORY
# ============================================================

processed_files <- if (
  dir.exists(processed_dir)
) {
  list.files(
    processed_dir,
    full.names = FALSE
  )
} else {
  character()
}


result_files <- if (
  dir.exists(results_dir)
) {
  list.files(
    results_dir,
    full.names = FALSE
  )
} else {
  character()
}


figure_files <- if (
  dir.exists(figures_dir)
) {
  list.files(
    figures_dir,
    full.names = FALSE
  )
} else {
  character()
}


cat("\n============================================================\n")
cat("OUTPUT INVENTORY\n")
cat("============================================================\n\n")

cat(
  "Processed files:",
  length(processed_files),
  "\n"
)

cat(
  "Result files:",
  length(result_files),
  "\n"
)

cat(
  "Figure files:",
  length(figure_files),
  "\n\n"
)


# ============================================================
# 19. FINAL CROSS-CHECK
# ============================================================

if (
  !is.null(script17_table) &&
  !is.null(script19_table) &&
  !is.null(script17_symbol_column)
) {
  
  targets17 <- unique(
    as.character(
      script17_table[[script17_symbol_column]]
    )
  )
  
  targets17 <- targets17[
    targets17 %in% final_targets
  ]
  
  
  targets19 <- unique(
    as.character(
      script19_table$SYMBOL
    )
  )
  
  
  if (
    setequal(
      targets17,
      targets19
    )
  ) {
    
    add_check(
      "Cross-analysis consistency",
      "Script 17 vs Script 19 target set",
      "PASS",
      "Target sets are identical."
    )
    
  } else {
    
    add_check(
      "Cross-analysis consistency",
      "Script 17 vs Script 19 target set",
      "FAIL",
      "Target sets differ."
    )
    
  }
  
}


# ============================================================
# 20. FINAL RANK CROSS-CHECK
# ============================================================

if (
  !is.null(script19_table) &&
  all(
    c(
      "SYMBOL",
      "final_rank"
    ) %in%
    names(script19_table)
  )
) {
  
  observed_rank <- setNames(
    script19_table$final_rank,
    script19_table$SYMBOL
  )
  
  
  common_targets <- intersect(
    names(expected_rank),
    names(observed_rank)
  )
  
  
  if (
    length(common_targets) == 12 &&
    all(
      expected_rank[
        common_targets
      ] ==
      observed_rank[
        common_targets
      ]
    )
  ) {
    
    add_check(
      "Cross-analysis consistency",
      "Final ranking consistency",
      "PASS",
      "Script 19 follows the established Script 17 ranking."
    )
    
  } else {
    
    add_check(
      "Cross-analysis consistency",
      "Final ranking consistency",
      "FAIL",
      "Ranking mismatch detected."
    )
    
  }
  
}


# ============================================================
# 21. FINAL SUMMARY
# ============================================================

pass_count <- sum(
  audit_results$status == "PASS"
)


warn_count <- sum(
  audit_results$status == "WARN"
)


fail_count <- sum(
  audit_results$status == "FAIL"
)


info_count <- sum(
  audit_results$status == "INFO"
)


cat("\n============================================================\n")
cat("FINAL AUDIT SUMMARY\n")
cat("============================================================\n\n")


cat(
  "PASS:",
  pass_count,
  "\n"
)


cat(
  "WARN:",
  warn_count,
  "\n"
)


cat(
  "FAIL:",
  fail_count,
  "\n"
)


cat(
  "INFO:",
  info_count,
  "\n\n"
)


# ============================================================
# 22. FAILURES
# ============================================================

if (fail_count > 0) {
  
  cat(
    "------------------------------------------------------------\n"
  )
  
  cat(
    "FAILURES REQUIRING ATTENTION\n"
  )
  
  cat(
    "------------------------------------------------------------\n\n"
  )
  
  
  print(
    audit_results[
      audit_results$status == "FAIL",
      ,
      drop = FALSE
    ],
    row.names = FALSE
  )
  
  
  cat("\n")
  
}


# ============================================================
# 23. WARNINGS
# ============================================================

if (warn_count > 0) {
  
  cat(
    "------------------------------------------------------------\n"
  )
  
  cat(
    "WARNINGS TO REVIEW\n"
  )
  
  cat(
    "------------------------------------------------------------\n\n"
  )
  
  
  print(
    audit_results[
      audit_results$status == "WARN",
      ,
      drop = FALSE
    ],
    row.names = FALSE
  )
  
  
  cat("\n")
  
}


# ============================================================
# 24. SAVE AUDIT CSV
# ============================================================

audit_csv <- file.path(
  results_dir,
  "Script20_final_project_audit.csv"
)


write.csv(
  audit_results,
  audit_csv,
  row.names = FALSE
)


# ============================================================
# 25. SAVE AUDIT TXT
# ============================================================

audit_txt <- file.path(
  results_dir,
  "Script20_final_project_audit.txt"
)


sink(
  audit_txt
)


cat(
  "============================================================\n"
)

cat(
  "SCRIPT 20: FINAL PROJECT AUDIT\n"
)

cat(
  "============================================================\n\n"
)


cat(
  "Project root:\n",
  project_root,
  "\n\n"
)


cat(
  "PASS:",
  pass_count,
  "\n"
)


cat(
  "WARN:",
  warn_count,
  "\n"
)


cat(
  "FAIL:",
  fail_count,
  "\n"
)


cat(
  "INFO:",
  info_count,
  "\n\n"
)


print(
  audit_results,
  row.names = FALSE
)


sink()


# ============================================================
# 26. FINAL STATUS
# ============================================================

cat(
  "\n============================================================\n"
)


if (fail_count == 0) {
  
  cat(
    "FINAL AUDIT STATUS: PASSED\n"
  )
  
  cat(
    "No critical consistency failures were detected.\n"
  )
  
} else {
  
  cat(
    "FINAL AUDIT STATUS: REVIEW REQUIRED\n"
  )
  
  cat(
    "Critical failures were detected.\n"
  )
  
}


cat(
  "============================================================\n\n"
)


cat(
  "Audit CSV:\n",
  audit_csv,
  "\n\n"
)


cat(
  "Audit TXT:\n",
  audit_txt,
  "\n\n"
)


if (fail_count == 0) {
  
  cat(
    "The R analysis can now be frozen.\n"
  )
  
} else {
  
  cat(
    "Resolve FAIL items before freezing the R analysis.\n"
  )
  
}


cat(
  "\n============================================================\n"
)