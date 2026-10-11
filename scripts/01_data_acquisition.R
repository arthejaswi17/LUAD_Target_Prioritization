# 01_TCGA-LUAD Data Acquisition


library(TCGAbiolinks)
library(httr2)


# 1. Define directories

raw_dir <- "data/raw/TCGA_LUAD"
table_dir <- "tables"

dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)


# 2. QUERY TCGA-LUAD STAR Counts

query_luad <- GDCquery(
  project = "TCGA-LUAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

query_results <- getResults(query_luad)


# 3. Identify patients with both Tumor and Normal samples

tumor_patients <- unique(
  query_results$cases.submitter_id[
    query_results$sample_type == "Primary Tumor"
  ]
)

normal_patients <- unique(
  query_results$cases.submitter_id[
    query_results$sample_type == "Solid Tissue Normal"
  ]
)

matched_patients <- intersect(
  tumor_patients,
  normal_patients
)

cat("Matched patients:", length(matched_patients), "\n")


# 4. Restrict query to matched patients

query_matched <- GDCquery(
  project = "TCGA-LUAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  barcode = matched_patients
)

matched_results <- getResults(query_matched)


# 5. Select one Tumor and one Normal sample per patient

tumor_samples <- matched_results[
  matched_results$sample_type == "Primary Tumor",
]

normal_samples <- matched_results[
  matched_results$sample_type == "Solid Tissue Normal",
]

# Remove duplicate patients.
# This reproduces the validated cohort used during acquisition.
tumor_samples <- tumor_samples[
  !duplicated(tumor_samples$cases.submitter_id),
]

normal_samples <- normal_samples[
  !duplicated(normal_samples$cases.submitter_id),
]

common_patients <- intersect(
  tumor_samples$cases.submitter_id,
  normal_samples$cases.submitter_id
)

tumor_samples <- tumor_samples[
  tumor_samples$cases.submitter_id %in% common_patients,
]

normal_samples <- normal_samples[
  normal_samples$cases.submitter_id %in% common_patients,
]


# 6. Create paired sample metadata

paired_metadata <- rbind(
  data.frame(
    patient = tumor_samples$cases.submitter_id,
    sample = tumor_samples$sample.submitter_id,
    condition = "Tumor",
    file_id = tumor_samples$file_id,
    file_name = tumor_samples$file_name,
    md5sum = tumor_samples$md5sum
  ),
  data.frame(
    patient = normal_samples$cases.submitter_id,
    sample = normal_samples$sample.submitter_id,
    condition = "Normal",
    file_id = normal_samples$file_id,
    file_name = normal_samples$file_name,
    md5sum = normal_samples$md5sum
  )
)

write.csv(
  paired_metadata,
  file.path(table_dir, "paired_sample_metadata.csv"),
  row.names = FALSE
)

cat("Tumor samples:", sum(paired_metadata$condition == "Tumor"), "\n")
cat("Normal samples:", sum(paired_metadata$condition == "Normal"), "\n")
cat("Unique patients:", length(unique(paired_metadata$patient)), "\n")


# 7. Download files from the GDC API

failed_files <- character(0)

for (i in seq_len(nrow(paired_metadata))) {
  
  file_id <- paired_metadata$file_id[i]
  file_name <- paired_metadata$file_name[i]
  
  output_file <- file.path(raw_dir, file_name)
  
  cat(sprintf(
    "[%d/%d] %s\n",
    i,
    nrow(paired_metadata),
    file_name
  ))
  
  # Do not download again if the file already exists.
  if (file.exists(output_file) &&
      file.info(output_file)$size > 0) {
    
    cat("   Already exists - skipping\n")
    next
  }
  
  tryCatch({
    
    response <- request(
      paste0(
        "https://api.gdc.cancer.gov/data/",
        file_id
      )
    ) |>
      req_perform()
    
    if (resp_status(response) != 200) {
      stop(
        paste(
          "HTTP status:",
          resp_status(response)
        )
      )
    }
    
    writeBin(
      resp_body_raw(response),
      output_file
    )
    
    if (!file.exists(output_file) ||
        file.info(output_file)$size == 0) {
      
      stop("File was not written correctly")
    }
    
    cat("Downloaded successfully\n")
    
  }, error = function(e) {
    
    cat(
      "   FAILED:",
      conditionMessage(e),
      "\n"
    )
    
    failed_files <<- c(
      failed_files,
      file_id
    )
  })
}


# 8. Basic download validation

downloaded_files <- list.files(
  raw_dir,
  pattern = "\\.tsv$",
  full.names = TRUE
)

cat(
  "Downloaded TSV files:",
  length(downloaded_files),
  "\n"
)

file_info <- file.info(downloaded_files)

cat(
  "Empty files:",
  sum(file_info$size == 0),
  "\n"
)


# 9. MD5 integrity check against GDC

md5_actual <- unname(
  tools::md5sum(downloaded_files)
)

md5_expected <- paired_metadata$md5sum[
  match(
    basename(downloaded_files),
    paired_metadata$file_name
  )
]

md5_match <- md5_actual == md5_expected

cat(
  "Files with matching MD5:",
  sum(md5_match, na.rm = TRUE),
  "/",
  length(downloaded_files),
  "\n"
)


# 10. Save acquisition metadata

saveRDS(
  paired_metadata,
  "data/processed/paired_metadata.rds"
)

saveRDS(
  downloaded_files,
  "data/processed/downloaded_files.rds"
)

cat("\nData acquisition complete.\n")