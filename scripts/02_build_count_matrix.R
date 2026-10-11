# 02_Build TCGA-LUAD Raw Count Matrix


library(data.table)


# 1. Load validated metadata and file list


paired_metadata <- readRDS(
  "data/processed/paired_metadata.rds"
)

downloaded_files <- readRDS(
  "data/processed/downloaded_files.rds"
)

cat("Samples:", length(downloaded_files), "\n")
cat("Metadata rows:", nrow(paired_metadata), "\n")

stopifnot(length(downloaded_files) == 116)
stopifnot(nrow(paired_metadata) == 116)


# 2. Match downloaded files to metadata

file_names <- basename(downloaded_files)

stopifnot("file_name" %in% names(paired_metadata))

file_index <- match(
  paired_metadata$file_name,
  file_names
)

stopifnot(!anyNA(file_index))

downloaded_files <- downloaded_files[file_index]

stopifnot(length(downloaded_files) == 116)
stopifnot(!anyDuplicated(paired_metadata$sample))


# 3. Read first file

cat("\nReading first file...\n")

first_file <- fread(
  downloaded_files[1],
  data.table = FALSE
)


# 4. Keep Ensembl gene rows
# Exclude PAR_Y duplicate representations

first_genes <- first_file[
  grepl("^ENSG", first_file$gene_id) &
    !grepl("_PAR_Y$", first_file$gene_id),
  c("gene_id", "gene_name", "gene_type", "unstranded")
]

# Remove Ensembl version suffix
first_genes$gene_id <- sub(
  "\\..*$",
  "",
  first_genes$gene_id
)

# Confirm uniqueness
stopifnot(!anyDuplicated(first_genes$gene_id))

gene_ids <- first_genes$gene_id

cat("Genes retained:", length(gene_ids), "\n")


# 5. Initialize count matrix

count_matrix <- matrix(
  nrow = length(gene_ids),
  ncol = length(downloaded_files)
)

rownames(count_matrix) <- gene_ids


# 6. Read all 116 files

for (i in seq_along(downloaded_files)) {
  
  cat(
    sprintf(
      "[%d/%d] %s\n",
      i,
      length(downloaded_files),
      basename(downloaded_files[i])
    )
  )
  
  x <- fread(
    downloaded_files[i],
    data.table = FALSE
  )
  
  # Keep Ensembl genes, excluding PAR_Y representations
  x <- x[
    grepl("^ENSG", x$gene_id) &
      !grepl("_PAR_Y$", x$gene_id),
    c("gene_id", "unstranded")
  ]
  
  # Remove Ensembl version suffix
  x$gene_id <- sub(
    "\\..*$",
    "",
    x$gene_id
  )
  
  # Confirm identical gene structure
  stopifnot(identical(x$gene_id, gene_ids))
  
  # Store raw unstranded counts
  count_matrix[, i] <- x$unstranded
}


# 7. Convert to integer matrix

count_matrix <- round(count_matrix)

storage.mode(count_matrix) <- "integer"

colnames(count_matrix) <- paired_metadata$sample


# 8. Validate count matrix

cat("COUNT MATRIX VALIDATION\n")

cat("\nDimensions:\n")
print(dim(count_matrix))

cat("\nMissing values:\n")
print(sum(is.na(count_matrix)))

cat("\nDuplicated genes:\n")
print(anyDuplicated(rownames(count_matrix)))

cat("\nDuplicated samples:\n")
print(anyDuplicated(colnames(count_matrix)))

cat("\nMinimum count:\n")
print(min(count_matrix))

cat("\nMaximum count:\n")
print(max(count_matrix))

cat("\nNon-integer counts:\n")
print(sum(count_matrix %% 1 != 0))

cat("\nNegative counts:\n")
print(sum(count_matrix < 0))

# Hard validation
stopifnot(
  ncol(count_matrix) == 116,
  sum(is.na(count_matrix)) == 0,
  anyDuplicated(rownames(count_matrix)) == 0,
  anyDuplicated(colnames(count_matrix)) == 0,
  min(count_matrix) >= 0
)


# 9. Save count matrix

saveRDS(
  count_matrix,
  "data/processed/TCGA_LUAD_raw_counts.rds"
)

write.csv(
  count_matrix,
  "data/processed/TCGA_LUAD_raw_counts.csv"
)

cat("Count matrix saved successfully.\n")