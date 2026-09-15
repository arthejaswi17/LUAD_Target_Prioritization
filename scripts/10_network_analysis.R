# ============================================================
# 10_network_analysis.R
# Network and module analysis of prioritized LUAD candidates
# ============================================================

library(data.table)
library(httr2)
library(gprofiler2)

# ------------------------------------------------------------
# 1. Load candidate and DESeq2 data
# ------------------------------------------------------------

candidate_priority <- readRDS(
  "data/processed/candidate_priority.rds"
)

res_shrunk <- readRDS(
  "data/processed/TCGA_LUAD_DESeq2_shrunken_results.rds"
)

# ------------------------------------------------------------
# 2. Basic checks
# ------------------------------------------------------------

cat("Candidate genes:", nrow(candidate_priority), "\n")
cat("DESeq2 genes:", nrow(res_shrunk), "\n")

cat(
  "Candidate symbols:",
  length(unique(candidate_priority$SYMBOL)),
  "\n"
)

# ------------------------------------------------------------
# 3. Prepare candidate table
# ------------------------------------------------------------

network_candidates <- copy(candidate_priority)

network_candidates <- network_candidates[
  !is.na(SYMBOL) &
    SYMBOL != ""
]

network_candidates <- unique(
  network_candidates,
  by = "SYMBOL"
)

# ------------------------------------------------------------
# 4. Prepare DESeq2 statistics
# ------------------------------------------------------------

de_table <- as.data.table(
  as.data.frame(res_shrunk)
)

de_table[
  ,
  gene := rownames(res_shrunk)
]

de_table <- de_table[
  ,
  .(
    gene,
    baseMean,
    lfcSE,
    pvalue,
    padj
  )
]

# ------------------------------------------------------------
# 5. Merge only the additional DESeq2 statistics
# ------------------------------------------------------------

network_candidates <- merge(
  network_candidates,
  de_table[
    ,
    .(
      gene,
      baseMean,
      lfcSE,
      pvalue
    )
  ],
  by = "gene",
  all.x = TRUE
)

# ------------------------------------------------------------
# 6. Check merge
# ------------------------------------------------------------

cat(
  "Candidates after DESeq2 merge:",
  nrow(network_candidates),
  "\n"
)

cat(
  "Candidates with DESeq2 statistics:",
  sum(!is.na(network_candidates$pvalue)),
  "\n"
)

# ------------------------------------------------------------
# 7. Save network candidate table
# ------------------------------------------------------------

saveRDS(
  network_candidates,
  "data/processed/network_candidates.rds"
)

write.csv(
  network_candidates,
  "data/processed/network_candidates.csv",
  row.names = FALSE
)

cat("\nStage 10 input preparation complete.\n")

# ------------------------------------------------------------
# 8. STRING network preparation
# ------------------------------------------------------------

# STRING API
string_url <- "https://string-db.org/api/tsv/network"

# Candidate gene symbols
string_genes <- unique(
  network_candidates$SYMBOL
)

cat(
  "\nNumber of genes for STRING:",
  length(string_genes),
  "\n"
)

# ------------------------------------------------------------
# 9. Query STRING
# ------------------------------------------------------------

string_response <- httr2::request(string_url) |>
  httr2::req_url_query(
    identifiers = paste(
      string_genes,
      collapse = "%0d"
    ),
    species = 9606,
    required_score = 700,
    caller_identity = "LUAD_Target_Prioritization"
  ) |>
  httr2::req_perform()

string_text <- httr2::resp_body_string(
  string_response
)

# Save raw STRING response
writeLines(
  string_text,
  "data/processed/STRING_network_raw.tsv"
)

# ------------------------------------------------------------
# 10. Parse STRING network
# ------------------------------------------------------------

string_network <- fread(
  text = string_text
)

cat(
  "STRING edges retrieved:",
  nrow(string_network),
  "\n"
)

# ------------------------------------------------------------
# 11. Inspect network structure
# ------------------------------------------------------------

cat(
  "STRING columns:\n"
)

print(
  names(string_network)
)

# ------------------------------------------------------------
# 12. Save network
# ------------------------------------------------------------

saveRDS(
  string_network,
  "data/processed/STRING_network.rds"
)

write.csv(
  string_network,
  "data/processed/STRING_network.csv",
  row.names = FALSE
)

cat(
  "\nSTRING network retrieval complete.\n"
)

# ------------------------------------------------------------
# 13. Build network edge list
# ------------------------------------------------------------

network_edges <- string_network[
  ,
  .(
    gene_A = preferredName_A,
    gene_B = preferredName_B,
    score = score
  )
]

# Keep only candidate-to-candidate interactions
network_edges <- network_edges[
  gene_A %in% string_genes &
    gene_B %in% string_genes
]

# Remove self-interactions
network_edges <- network_edges[
  gene_A != gene_B
]

cat(
  "\nCandidate-to-candidate edges:",
  nrow(network_edges),
  "\n"
)

# ------------------------------------------------------------
# 14. Build igraph object
# ------------------------------------------------------------

library(igraph)

network_graph <- graph_from_data_frame(
  network_edges,
  directed = FALSE
)

# Remove duplicate edges
network_graph <- simplify(
  network_graph,
  remove.multiple = TRUE,
  remove.loops = TRUE
)

cat(
  "Network nodes:",
  vcount(network_graph),
  "\n"
)

cat(
  "Network edges:",
  ecount(network_graph),
  "\n"
)

cat(
  "Network density:",
  round(
    edge_density(network_graph),
    4
  ),
  "\n"
)

cat(
  "Connected components:",
  components(network_graph)$no,
  "\n"
)

# ------------------------------------------------------------
# 15. Network degree
# ------------------------------------------------------------

degree_table <- data.table(
  SYMBOL = V(network_graph)$name,
  degree = degree(network_graph)
)

degree_table <- merge(
  degree_table,
  network_candidates[
    ,
    .(
      SYMBOL,
      gene,
      n_pathways,
      log2FoldChange,
      padj
    )
  ],
  by = "SYMBOL",
  all.x = TRUE
)

degree_table <- degree_table[
  order(-degree)
]

cat(
  "\nTop 20 genes by network degree:\n"
)

print(
  degree_table[1:min(20, .N)]
)

# ------------------------------------------------------------
# 16. Save topology table
# ------------------------------------------------------------

saveRDS(
  degree_table,
  "data/processed/STRING_network_degree.rds"
)

write.csv(
  degree_table,
  "data/processed/STRING_network_degree.csv",
  row.names = FALSE
)

cat(
  "\nNetwork topology analysis complete.\n"
)

# ------------------------------------------------------------
# 17. Community detection
# ------------------------------------------------------------

set.seed(123)

leiden_modules <- cluster_leiden(
  network_graph,
  objective_function = "modularity",
  resolution = 1
)

module_membership <- data.table(
  SYMBOL = names(membership(leiden_modules)),
  module = as.integer(
    membership(leiden_modules)
  )
)

module_sizes <- module_membership[
  ,
  .(n_genes = .N),
  by = module
][
  order(-n_genes)
]

cat(
  "\nNumber of network modules:",
  length(unique(module_membership$module)),
  "\n"
)

cat(
  "\nModule sizes:\n"
)

print(module_sizes)

# ------------------------------------------------------------
# 18. Add module assignments to candidate table
# ------------------------------------------------------------

network_candidates_modules <- merge(
  network_candidates,
  module_membership,
  by = "SYMBOL",
  all.x = TRUE
)

network_candidates_modules <- network_candidates_modules[
  order(module, -log2FoldChange)
]

cat(
  "\nCandidates assigned to a module:",
  sum(!is.na(network_candidates_modules$module)),
  "\n"
)

cat(
  "Candidates without a module:",
  sum(is.na(network_candidates_modules$module)),
  "\n"
)

# ------------------------------------------------------------
# 19. Save module assignments
# ------------------------------------------------------------

saveRDS(
  module_membership,
  "data/processed/STRING_module_membership.rds"
)

saveRDS(
  network_candidates_modules,
  "data/processed/network_candidates_modules.rds"
)

write.csv(
  module_membership,
  "data/processed/STRING_module_membership.csv",
  row.names = FALSE
)

write.csv(
  network_candidates_modules,
  "data/processed/network_candidates_modules.csv",
  row.names = FALSE
)

cat(
  "\nCommunity detection complete.\n"
)

# ------------------------------------------------------------
# 20. Functional enrichment of major network modules
# ------------------------------------------------------------

major_modules <- module_sizes[
  n_genes >= 10,
  module
]

cat(
  "\nModules selected for enrichment:",
  paste(major_modules, collapse = ", "),
  "\n"
)

module_enrichment <- list()

for (m in major_modules) {
  
  genes_m <- module_membership[
    module == m,
    SYMBOL
  ]
  
  cat(
    "\nAnalyzing module",
    m,
    "with",
    length(genes_m),
    "genes...\n"
  )
  
  gost_result <- gost(
    query = genes_m,
    organism = "hsapiens",
    sources = c(
      "GO:BP",
      "REAC"
    ),
    correction_method = "fdr",
    significant = TRUE
  )
  
  module_enrichment[[paste0("module_", m)]] <- gost_result
}

# ------------------------------------------------------------
# 21. Save enrichment results
# ------------------------------------------------------------

saveRDS(
  module_enrichment,
  "data/processed/STRING_module_enrichment.rds"
)

cat(
  "\nModule enrichment complete.\n"
)

# ------------------------------------------------------------
# 22. Summarize top enriched terms per module
# ------------------------------------------------------------

module_enrichment_summary <- rbindlist(
  lapply(
    names(module_enrichment),
    function(x) {
      
      result <- module_enrichment[[x]]$result
      
      if (is.null(result) || nrow(result) == 0) {
        return(NULL)
      }
      
      result <- as.data.table(result)
      
      result <- result[
        order(p_value)
      ]
      
      result <- result[
        1:min(10, .N),
        .(
          term_id,
          term_name,
          source,
          p_value,
          intersection_size
        )
      ]
      
      result[
        ,
        module := sub("module_", "", x)
      ]
      
      result
    }
  ),
  fill = TRUE
)

# ------------------------------------------------------------
# 23. Order results
# ------------------------------------------------------------

module_enrichment_summary[
  ,
  module_num := as.integer(module)
]

setorder(
  module_enrichment_summary,
  module_num,
  p_value
)

cat(
  "\nTop enriched terms by module:\n\n"
)

print(
  module_enrichment_summary[
    ,
    .(
      module,
      term_name,
      source,
      p_value,
      intersection_size
    )
  ]
)

# ------------------------------------------------------------
# 24. Save enrichment summary
# ------------------------------------------------------------

module_enrichment_summary[
  ,
  module_num := NULL
]

saveRDS(
  module_enrichment_summary,
  "data/processed/STRING_module_enrichment_summary.rds"
)

write.csv(
  module_enrichment_summary,
  "data/processed/STRING_module_enrichment_summary.csv",
  row.names = FALSE
)

cat(
  "\nModule enrichment summary saved.\n"
)