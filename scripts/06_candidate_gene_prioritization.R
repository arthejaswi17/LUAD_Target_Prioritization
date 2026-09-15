# ============================================================
# SCRIPT 06: CANDIDATE GENE PRIORITIZATION
# ============================================================
#
# Purpose:
# Integrate:
#   1. DESeq2 differential expression
#   2. Hallmark GSEA
#   3. Hallmark leading-edge genes
#
# Output:
#   - candidate genes derived from DE + GSEA
#   - candidate prioritization table
#   - leading-edge gene table
#
# ============================================================

cat("\n")
cat("============================================================\n")
cat("SCRIPT 06: CANDIDATE GENE PRIORITIZATION\n")
cat("============================================================\n\n")


# ============================================================
# 1. LOAD PACKAGES
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
})


# ============================================================
# 2. DIRECTORIES
# ============================================================

processed_dir <- "data/processed"
figures_dir   <- "figures"

if (!dir.exists(processed_dir)) {
  dir.create(processed_dir, recursive = TRUE)
}

if (!dir.exists(figures_dir)) {
  dir.create(figures_dir, recursive = TRUE)
}


# ============================================================
# 3. LOAD DESEQ2 RESULTS
# ============================================================

cat("Loading DESeq2 results...\n")

de_results <- readRDS(
  file.path(
    processed_dir,
    "TCGA_LUAD_DESeq2_shrunken_results.rds"
  )
)

cat("Genes loaded:", nrow(de_results), "\n")
cat("Columns:", paste(colnames(de_results), collapse = ", "), "\n\n")


# ============================================================
# 4. CHECK DESEQ2 INPUT
# ============================================================

required_de_columns <- c(
  "gene",
  "baseMean",
  "log2FoldChange",
  "lfcSE",
  "pvalue",
  "padj"
)

missing_de_columns <- setdiff(
  required_de_columns,
  colnames(de_results)
)

if (length(missing_de_columns) > 0) {
  
  stop(
    paste(
      "Missing DESeq2 columns:",
      paste(missing_de_columns, collapse = ", ")
    )
  )
}


# ============================================================
# 5. REMOVE INVALID VALUES
# ============================================================

de_results <- de_results %>%
  dplyr::filter(
    !is.na(log2FoldChange),
    !is.na(padj),
    !is.na(gene)
  )


# ============================================================
# 6. DEFINE PRACTICAL DIFFERENTIALLY EXPRESSED GENES
# ============================================================

de_practical <- de_results %>%
  dplyr::filter(
    padj < 0.05,
    abs(log2FoldChange) >= 1
  )

de_up <- de_practical %>%
  dplyr::filter(
    log2FoldChange >= 1
  )

de_down <- de_practical %>%
  dplyr::filter(
    log2FoldChange <= -1
  )


cat("Practical DE genes:\n")
cat(
  "FDR < 0.05 and |log2FC| >= 1:",
  nrow(de_practical),
  "\n"
)

cat(
  "Upregulated:",
  nrow(de_up),
  "\n"
)

cat(
  "Downregulated:",
  nrow(de_down),
  "\n\n"
)


# ============================================================
# 7. LOAD HALLMARK GSEA RESULTS
# ============================================================

cat("Loading Hallmark GSEA results...\n")

gsea_file <- file.path(
  processed_dir,
  "GSEA_Hallmark_results.rds"
)

if (!file.exists(gsea_file)) {
  
  stop(
    paste(
      "GSEA results not found:",
      gsea_file
    )
  )
}

gsea_results <- readRDS(gsea_file)

cat(
  "GSEA pathways loaded:",
  nrow(gsea_results),
  "\n\n"
)

cat("GSEA columns:\n")
print(colnames(gsea_results))
cat("\n")


# ============================================================
# 8. CHECK GSEA STRUCTURE
# ============================================================

required_gsea_columns <- c(
  "pathway",
  "NES",
  "pval",
  "padj",
  "leadingEdge"
)

missing_gsea_columns <- setdiff(
  required_gsea_columns,
  colnames(gsea_results)
)

if (length(missing_gsea_columns) > 0) {
  
  stop(
    paste(
      "Missing GSEA columns:",
      paste(
        missing_gsea_columns,
        collapse = ", "
      )
    )
  )
}


# ============================================================
# 9. SIGNIFICANT GSEA PATHWAYS
# ============================================================

significant_gsea <- gsea_results %>%
  dplyr::filter(
    padj < 0.05
  )

cat(
  "Significant Hallmark pathways:",
  nrow(significant_gsea),
  "\n\n"
)


# ============================================================
# 10. EXTRACT LEADING-EDGE GENES
# ============================================================

cat(
  "Extracting leading-edge genes from significant pathways...\n"
)

# Keep only the required GSEA columns
leading_edge_base <- significant_gsea %>%
  dplyr::select(
    pathway,
    NES,
    padj,
    direction,
    leadingEdge
  )

# Convert the fgsea leadingEdge list-column into
# one row per pathway-gene relationship.
leading_edge_table <- do.call(
  rbind,
  lapply(
    seq_len(nrow(leading_edge_base)),
    function(i) {
      
      genes <- leading_edge_base$leadingEdge[[i]]
      
      if (length(genes) == 0) {
        return(NULL)
      }
      
      data.frame(
        pathway = leading_edge_base$pathway[i],
        NES = leading_edge_base$NES[i],
        padj = leading_edge_base$padj[i],
        direction = leading_edge_base$direction[i],
        gene = as.character(genes),
        stringsAsFactors = FALSE
      )
    }
  )
)

# Convert to tibble
leading_edge_table <- tibble::as_tibble(
  leading_edge_table
)


# ============================================================
# 11. CLEAN LEADING-EDGE GENES
# ============================================================

leading_edge_table <- leading_edge_table %>%
  dplyr::filter(
    !is.na(gene),
    gene != ""
  ) %>%
  dplyr::distinct(
    pathway,
    gene,
    .keep_all = TRUE
  )


cat(
  "Leading-edge gene-pathway relationships:",
  nrow(leading_edge_table),
  "\n"
)

cat(
  "Unique leading-edge genes:",
  dplyr::n_distinct(
    leading_edge_table$gene
  ),
  "\n\n"
)


# ============================================================
# 12. COUNT PATHWAY SUPPORT PER GENE
# ============================================================

gene_pathway_support <- leading_edge_table %>%
  dplyr::group_by(
    gene
  ) %>%
  dplyr::summarise(
    n_GSEA_pathways = dplyr::n_distinct(pathway),
    GSEA_pathways = paste(
      unique(pathway),
      collapse = "; "
    ),
    .groups = "drop"
  )


# ============================================================
# 13. ADD GSEA INFORMATION TO DE RESULTS
# ============================================================

candidate_table <- de_practical %>%
  dplyr::left_join(
    gene_pathway_support,
    by = c("gene" = "gene")
  )


# ============================================================
# 14. FILL GENES WITHOUT GSEA SUPPORT
# ============================================================

candidate_table <- candidate_table %>%
  dplyr::mutate(
    n_GSEA_pathways = ifelse(
      is.na(n_GSEA_pathways),
      0,
      n_GSEA_pathways
    ),
    GSEA_pathways = ifelse(
      is.na(GSEA_pathways),
      "",
      GSEA_pathways
    )
  )


# ============================================================
# 15. ADD DE DIRECTION
# ============================================================

candidate_table <- candidate_table %>%
  dplyr::mutate(
    DE_direction = dplyr::case_when(
      log2FoldChange >= 1 ~ "Upregulated",
      log2FoldChange <= -1 ~ "Downregulated",
      TRUE ~ "Not_DE"
    )
  )


# ============================================================
# 16. CALCULATE PRIORITIZATION SCORE
# ============================================================
#
# The score combines:
#
#   - magnitude of expression change
#   - statistical significance
#   - number of GSEA pathways containing the gene
#
# This is an exploratory prioritization score.
# It is NOT a statistical test.
#
# ============================================================

candidate_table <- candidate_table %>%
  dplyr::mutate(
    
    abs_log2FC = abs(log2FoldChange),
    
    neg_log10_FDR = -log10(
      pmax(
        padj,
        .Machine$double.xmin
      )
    ),
    
    prioritization_score =
      abs_log2FC +
      log1p(n_GSEA_pathways) +
      (neg_log10_FDR / 10)
    
  )


# ============================================================
# 17. SORT CANDIDATES
# ============================================================

candidate_table <- candidate_table %>%
  dplyr::arrange(
    dplyr::desc(prioritization_score)
  )


# ============================================================
# 18. SUMMARY
# ============================================================

cat("\n")
cat("============================================================\n")
cat("CANDIDATE PRIORITIZATION SUMMARY\n")
cat("============================================================\n\n")

cat(
  "Practical DE genes:",
  nrow(de_practical),
  "\n"
)

cat(
  "Genes with GSEA leading-edge support:",
  sum(candidate_table$n_GSEA_pathways > 0),
  "\n"
)

cat(
  "Genes without GSEA leading-edge support:",
  sum(candidate_table$n_GSEA_pathways == 0),
  "\n\n"
)


# ============================================================
# 19. TOP CANDIDATES
# ============================================================

cat("\nTop candidate genes:\n\n")

# Select the top 30 candidates based on the
# integrated prioritization score.
top_candidates <- candidate_table %>%
  dplyr::slice_head(n = 30)

# Create a clean display table for the console.
top_display <- top_candidates %>%
  dplyr::select(
    gene,
    log2FoldChange,
    padj,
    DE_direction,
    n_GSEA_pathways,
    GSEA_pathways,
    prioritization_score
  ) %>%
  dplyr::mutate(
    log2FoldChange = round(log2FoldChange, 3),
    padj = signif(padj, 3),
    prioritization_score = round(
      prioritization_score,
      3
    )
  ) %>%
  as.data.frame()

# Print only a manageable number of rows.
print(
  top_display,
  row.names = FALSE
)


# ============================================================
# 20. SAVE FULL CANDIDATE TABLE
# ============================================================

saveRDS(
  candidate_table,
  file.path(
    processed_dir,
    "candidate_gene_prioritization.rds"
  )
)

readr::write_csv(
  candidate_table,
  file.path(
    processed_dir,
    "candidate_gene_prioritization.csv"
  )
)


# ============================================================
# 21. SAVE LEADING-EDGE TABLE
# ============================================================

saveRDS(
  leading_edge_table,
  file.path(
    processed_dir,
    "GSEA_leading_edge_genes.rds"
  )
)

readr::write_csv(
  leading_edge_table,
  file.path(
    processed_dir,
    "GSEA_leading_edge_genes.csv"
  )
)


# ============================================================
# 22. SAVE TOP CANDIDATES
# ============================================================

saveRDS(
  top_candidates,
  file.path(
    processed_dir,
    "top_30_candidate_genes.rds"
  )
)

readr::write_csv(
  top_candidates,
  file.path(
    processed_dir,
    "top_30_candidate_genes.csv"
  )
)


# ============================================================
# 23. FINAL MESSAGE
# ============================================================

cat("\n")
cat("============================================================\n")
cat("SCRIPT 06 COMPLETED SUCCESSFULLY\n")
cat("============================================================\n\n")

cat("Outputs saved:\n")

cat("- candidate_gene_prioritization.rds\n")
cat("- candidate_gene_prioritization.csv\n")
cat("- GSEA_leading_edge_genes.rds\n")
cat("- GSEA_leading_edge_genes.csv\n")
cat("- top_30_candidate_genes.rds\n")
cat("- top_30_candidate_genes.csv\n")

cat("\n")
cat("============================================================\n")