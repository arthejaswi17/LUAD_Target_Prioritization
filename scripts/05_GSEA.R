# SCRIPT 05
# HALLMARK GENE SET ENRICHMENT ANALYSIS (GSEA)
# Perform ranked GSEA on TCGA-LUAD differential expression results using MSigDB Hallmark gene sets.
# Input:
#   TCGA_LUAD_DESeq2_shrunken_results.rds
# Packages:
#   dplyr
#   readr
#   fgsea
#   msigdbr



# 1. PACKAGES

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(fgsea)
  library(msigdbr)
})


# 2. DIRECTORIES

project_root <- getwd()
processed_dir <- file.path(
  project_root,
  "data",
  "processed"
)
figures_dir <- file.path(
  project_root,
  "figures"
)
if (!dir.exists(figures_dir)) {
  dir.create(
    figures_dir,
    recursive = TRUE
  )
}


# 3. LOAD DESEQ2 RESULTS

de_results <- readRDS(
  file.path(
    processed_dir,
    "TCGA_LUAD_DESeq2_shrunken_results.rds"
  )
)
cat(
  "\nDESeq2 results loaded.\n"
)
cat(
  "Genes:",
  nrow(de_results),
  "\n"
)
cat(
  "Columns:",
  paste(
    colnames(de_results),
    collapse = ", "
  ),
  "\n"
)


# 4. CHECK REQUIRED COLUMNS

required_columns <- c(
  "gene",
  "log2FoldChange",
  "padj"
)
missing_columns <- setdiff(
  required_columns,
  colnames(de_results)
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


# 5. PREPARE RANKED GENE LIST

gsea_input <- de_results %>%
  dplyr::select(
    gene,
    log2FoldChange,
    padj
  ) %>%
  dplyr::filter(
    !is.na(gene),
    !is.na(log2FoldChange)
  ) %>%
  dplyr::mutate(
    gene = as.character(gene)
  ) %>%
  dplyr::arrange(
    dplyr::desc(
      abs(log2FoldChange)
    )
  ) %>%
  dplyr::distinct(
    gene,
    .keep_all = TRUE
  )


# 6. CREATE RANKED STATISTIC
# For GSEA we use the shrunken log2 fold change directly.
# Positive values:
#   genes upregulated in LUAD
# Negative values:
#   genes downregulated in LUAD

gene_ranks <- gsea_input$log2FoldChange
names(gene_ranks) <- gsea_input$gene
# Remove duplicated names if any
gene_ranks <- gene_ranks[
  !duplicated(names(gene_ranks))
]
# Sort decreasing
gene_ranks <- sort(
  gene_ranks,
  decreasing = TRUE
)
cat(
  "\nRanked genes:",
  length(gene_ranks),
  "\n"
)
cat(
  "Maximum log2FC:",
  max(
    gene_ranks,
    na.rm = TRUE
  ),
  "\n"
)
cat(
  "Minimum log2FC:",
  min(
    gene_ranks,
    na.rm = TRUE
  ),
  "\n"
)


# 7. LOAD MSIGDB HALLMARK GENE SETS

cat(
  "\nLoading MSigDB Hallmark gene sets...\n"
)
hallmark <- msigdbr::msigdbr(
  species = "Homo sapiens",
  collection = "H"
)
cat(
  "Hallmark rows:",
  nrow(hallmark),
  "\n"
)
cat(
  "Hallmark pathways:",
  dplyr::n_distinct(
    hallmark$gs_name
  ),
  "\n"
)


# 8. PREPARE PATHWAYS FOR FGSEA

hallmark_pathways <- split(
  hallmark$ensembl_gene,
  hallmark$gs_name
)
# Remove missing Ensembl IDs
hallmark_pathways <- lapply(
  hallmark_pathways,
  function(x) {
    x <- x[
      !is.na(x) &
        x != ""
    ]
    unique(x)
  }
)
# Remove empty pathways
hallmark_pathways <- hallmark_pathways[
  lengths(hallmark_pathways) > 0
]
cat(
  "Usable Hallmark pathways:",
  length(hallmark_pathways),
  "\n"
)


# 9. CHECK OVERLAP

overlap_genes <- intersect(
  names(gene_ranks),
  unique(
    hallmark$ensembl_gene
  )
)
cat(
  "\nGene overlap between DESeq2 and Hallmark:",
  length(overlap_genes),
  "\n"
)
if (length(overlap_genes) < 1000) {
  warning(
    "Low overlap between DESeq2 genes and Hallmark gene sets."
  )
}


# 10. RUN FGSEA

cat(
  "\nRunning fgsea...\n"
)
set.seed(123)
gsea_results <- fgsea::fgsea(
  pathways = hallmark_pathways,
  stats = gene_ranks,
  minSize = 15,
  maxSize = 500,
  eps = 0
)


# 11. FORMAT RESULTS

gsea_results <- as.data.frame(
  gsea_results
) %>%
  dplyr::arrange(
    padj,
    dplyr::desc(
      abs(NES)
    )
  )


# 12. ADD INTERPRETATION

gsea_results <- gsea_results %>%
  dplyr::mutate(
    direction = dplyr::case_when(
      NES > 0 ~
        "Enriched in LUAD / upregulated genes",
      NES < 0 ~
        "Enriched in LUAD / downregulated genes",
      TRUE ~
        "No directional enrichment"
    ),
    significant_FDR05 =
      !is.na(padj) &
      padj < 0.05
  )


# 13. PRINT RESULTS

cat(
  "HALLMARK GSEA RESULTS\n"
)
cat(
  "Total pathways tested:",
  nrow(gsea_results),
  "\n"
)
cat(
  "Significant pathways (FDR < 0.05):",
  sum(
    gsea_results$significant_FDR05,
    na.rm = TRUE
  ),
  "\n"
)
cat(
  "Positive NES pathways:",
  sum(
    gsea_results$NES > 0,
    na.rm = TRUE
  ),
  "\n"
)
cat(
  "Negative NES pathways:",
  sum(
    gsea_results$NES < 0,
    na.rm = TRUE
  ),
  "\n"
)


# 14. PRINT TOP PATHWAYS

cat(
  "\nTop enriched Hallmark pathways:\n\n"
)
print(
  gsea_results %>%
    dplyr::select(
      pathway,
      NES,
      pval,
      padj,
      size,
      direction
    ) %>%
    dplyr::slice_head(
      n = 20
    ),
  row.names = FALSE
)


# 15. SIGNIFICANT PATHWAYS

gsea_significant <- gsea_results %>%
  
  dplyr::filter(
    significant_FDR05
  ) %>%
  
  dplyr::arrange(
    padj
  )


# 16. SAVE RESULTS

# Convert list-columns such as leadingEdge into text
gsea_results_export <- gsea_results
for (i in seq_along(gsea_results_export)) {
  if (is.list(gsea_results_export[[i]])) {
    gsea_results_export[[i]] <- vapply(
      gsea_results_export[[i]],
      function(x) {
        paste(
          x,
          collapse = ";"
        )
      },
      character(1)
    )
  }
}
# Save complete GSEA results
saveRDS(
  gsea_results,
  file.path(
    processed_dir,
    "GSEA_Hallmark_results.rds"
  )
)
write.csv(
  gsea_results_export,
  file.path(
    processed_dir,
    "GSEA_Hallmark_results.csv"
  ),
  row.names = FALSE
)
# Significant pathways
gsea_significant <- gsea_results %>%
  dplyr::filter(
    significant_FDR05
  ) %>%
  dplyr::arrange(
    padj
  )
gsea_significant_export <- gsea_significant
for (i in seq_along(gsea_significant_export)) {
  if (is.list(gsea_significant_export[[i]])) {
    gsea_significant_export[[i]] <- vapply(
      gsea_significant_export[[i]],
      function(x) {
        paste(
          x,
          collapse = ";"
        )
      },
      character(1)
    )
  }
}
saveRDS(
  gsea_significant,
  file.path(
    processed_dir,
    "GSEA_Hallmark_significant.rds"
  )
)
write.csv(
  gsea_significant_export,
  file.path(
    processed_dir,
    "GSEA_Hallmark_significant.csv"
  ),
  row.names = FALSE
)
cat(
  "\nGSEA results exported successfully.\n"
)


# 17. TOP PATHWAY PLOT

if (nrow(gsea_significant) > 0) {
  top_plot <- gsea_significant %>%
    dplyr::arrange(
      NES
    ) %>%
    dplyr::slice_head(
      n = min(
        15,
        nrow(gsea_significant)
      )
    )
  # Convert pathway names to readable format
  top_plot$pathway <- gsub(
    "^HALLMARK_",
    "",
    top_plot$pathway
  )
  top_plot$pathway <- gsub(
    "_",
    " ",
    top_plot$pathway
  )
  top_plot$pathway <- factor(
    top_plot$pathway,
    levels = top_plot$pathway
  )
  p <- ggplot2::ggplot(
    top_plot,
    ggplot2::aes(
      x = NES,
      y = pathway
    )
  ) +
    ggplot2::geom_point(
      ggplot2::aes(
        size = -log10(padj)
      )
    ) +
    ggplot2::geom_vline(
      xintercept = 0,
      linetype = "dashed"
    ) +
    ggplot2::labs(
      title = "Hallmark GSEA in TCGA-LUAD",
      x = "Normalized Enrichment Score (NES)",
      y = NULL,
      size = "-log10(FDR)"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title =
        ggplot2::element_text(
          face = "bold"
        ),
      axis.text.y =
        ggplot2::element_text(
          size = 10
        )
    )
  ggplot2::ggsave(
    filename = file.path(
      figures_dir,
      "GSEA_Hallmark_top_pathways.png"
    ),
    plot = p,
    width = 9,
    height = 7,
    dpi = 300
  )
  ggplot2::ggsave(
    filename = file.path(
      figures_dir,
      "GSEA_Hallmark_top_pathways.pdf"
    ),
    plot = p,
    width = 9,
    height = 7
  )
}


# 18. SESSION INFORMATION

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    processed_dir,
    "GSEA_Hallmark_sessionInfo.txt"
  )
)


# 19. FINAL CHECK

cat(
  "GSEA ANALYSIS COMPLETED SUCCESSFULLY\n"
)
cat(
  "Hallmark pathways tested:",
  nrow(gsea_results),
  "\n"
)
cat(
  "Significant pathways:",
  nrow(gsea_significant),
  "\n"
)
cat(
  "\nOutputs saved:\n"
)
cat(
  "GSEA_Hallmark_results.rds\n"
)
cat(
  "GSEA_Hallmark_results.csv\n"
)
cat(
  "GSEA_Hallmark_significant.rds\n"
)
cat(
  "GSEA_Hallmark_significant.csv\n"
)
cat(
  "figures/GSEA_Hallmark_top_pathways.png\n"
)
cat(
  "figures/GSEA_Hallmark_top_pathways.pdf\n"
)
cat(
  "GSEA_Hallmark_sessionInfo.txt\n"
)