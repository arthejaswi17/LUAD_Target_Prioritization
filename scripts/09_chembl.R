# ============================================================
# 09_chembl.R
# ChEMBL pharmacological evidence for prioritized LUAD targets
# ============================================================

library(data.table)

# ------------------------------------------------------------
# 1. Load previously retrieved ChEMBL data
# ------------------------------------------------------------

chembl_quantitative <- readRDS(
  "data/processed/ChEMBL_quantitative_activity.rds"
)

gmnn_activity <- readRDS(
  "data/processed/ChEMBL_GMNN_activity.rds"
)

# ------------------------------------------------------------
# 2. Convert pChEMBL values to numeric
# ------------------------------------------------------------

chembl_quantitative[
  ,
  pchembl_numeric := as.numeric(pchembl_value)
]

gmnn_activity[
  ,
  pchembl_numeric := as.numeric(pchembl_value)
]

# ------------------------------------------------------------
# 3. Build summary for canonical quantitative targets
# ------------------------------------------------------------

chembl_evidence_summary <- chembl_quantitative[
  !is.na(pchembl_numeric),
  .(
    n_records = .N,
    n_pchembl = sum(!is.na(pchembl_numeric)),
    n_compounds = uniqueN(molecule_chembl_id),
    n_assays = uniqueN(assay_chembl_id),
    n_documents = uniqueN(document_chembl_id),
    median_pchembl = median(
      pchembl_numeric,
      na.rm = TRUE
    ),
    p90_pchembl = as.numeric(
      quantile(
        pchembl_numeric,
        probs = 0.90,
        na.rm = TRUE
      )
    ),
    max_pchembl = max(
      pchembl_numeric,
      na.rm = TRUE
    )
  ),
  by = SYMBOL
]

# ------------------------------------------------------------
# 4. Build GMNN summary
# ------------------------------------------------------------

gmnn_summary <- gmnn_activity[
  !is.na(pchembl_numeric),
  .(
    n_records = .N,
    n_pchembl = sum(!is.na(pchembl_numeric)),
    n_compounds = uniqueN(molecule_chembl_id),
    n_assays = uniqueN(assay_chembl_id),
    n_documents = uniqueN(document_chembl_id),
    median_pchembl = median(
      pchembl_numeric,
      na.rm = TRUE
    ),
    p90_pchembl = as.numeric(
      quantile(
        pchembl_numeric,
        probs = 0.90,
        na.rm = TRUE
      )
    ),
    max_pchembl = max(
      pchembl_numeric,
      na.rm = TRUE
    )
  ),
  by = SYMBOL
]

# ------------------------------------------------------------
# 5. Combine GMNN with the other targets
# ------------------------------------------------------------

chembl_evidence_summary <- rbind(
  chembl_evidence_summary,
  gmnn_summary,
  fill = TRUE
)

# ------------------------------------------------------------
# 6. Flag API-window truncation
# ------------------------------------------------------------

chembl_evidence_summary[
  ,
  retrieval_limit := n_records >= 1000
]

# ------------------------------------------------------------
# 7. Sort by therapeutic shortlist
# ------------------------------------------------------------

therapeutic_candidates <- c(
  "AURKA",
  "CDCA5",
  "GMNN",
  "SKA1",
  "SGO1",
  "PRKDC",
  "BRCA1",
  "PLK1",
  "CDK1",
  "TOP2A",
  "DNA2",
  "CHEK1"
)

chembl_evidence_summary[
  ,
  sort_order := match(
    SYMBOL,
    therapeutic_candidates
  )
]

setorder(
  chembl_evidence_summary,
  sort_order
)

chembl_evidence_summary[
  ,
  sort_order := NULL
]

# ------------------------------------------------------------
# 8. Save final outputs
# ------------------------------------------------------------

saveRDS(
  chembl_quantitative,
  "data/processed/ChEMBL_quantitative_activity.rds"
)

saveRDS(
  gmnn_activity,
  "data/processed/ChEMBL_GMNN_activity.rds"
)

saveRDS(
  chembl_evidence_summary,
  "data/processed/ChEMBL_evidence_summary.rds"
)

write.csv(
  chembl_quantitative,
  "data/processed/ChEMBL_quantitative_activity.csv",
  row.names = FALSE
)

write.csv(
  chembl_evidence_summary,
  "data/processed/ChEMBL_evidence_summary.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 9. Final checks
# ------------------------------------------------------------

cat("\n==============================\n")
cat("ChEMBL stage completed\n")
cat("==============================\n\n")

cat(
  "Targets represented:\n",
  paste(
    chembl_evidence_summary$SYMBOL,
    collapse = ", "
  ),
  "\n\n"
)

print(chembl_evidence_summary)

cat(
  "\nGMNN saved records:",
  nrow(gmnn_activity),
  "\n"
)

cat("\nOutputs saved successfully.\n")