# ============================================================
# LUAD Target Prioritization Project
# From Transcriptional Dysregulation to Therapeutic
# Target Prioritization in Lung Adenocarcinoma
# ============================================================

# Project setup
# R version: 4.5.2
# Bioconductor version: 3.22
# Date: 2026-09-08


# ------------------------------------------------------------
# 1. Load required packages
# ------------------------------------------------------------

library(DESeq2)
library(limma)
library(GEOquery)
library(apeglm)
library(TCGAbiolinks)

# ------------------------------------------------------------
# 2. Create project directories if they don't exist
# ------------------------------------------------------------

dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

dir.create("results/QC", recursive = TRUE, showWarnings = FALSE)
dir.create("results/differential_expression",
           recursive = TRUE, showWarnings = FALSE)
dir.create("results/pathways",
           recursive = TRUE, showWarnings = FALSE)
dir.create("results/networks",
           recursive = TRUE, showWarnings = FALSE)
dir.create("results/dependency",
           recursive = TRUE, showWarnings = FALSE)
dir.create("results/target_prioritization",
           recursive = TRUE, showWarnings = FALSE)
dir.create("results/validation",
           recursive = TRUE, showWarnings = FALSE)

dir.create("figures", showWarnings = FALSE)
dir.create("tables", showWarnings = FALSE)
dir.create("docs", showWarnings = FALSE)


# ------------------------------------------------------------
# 3. Record software versions
# ------------------------------------------------------------

cat("R version:", R.version.string, "\n")
cat("DESeq2 version:", as.character(packageVersion("DESeq2")), "\n")
cat("limma version:", as.character(packageVersion("limma")), "\n")
cat("GEOquery version:", as.character(packageVersion("GEOquery")), "\n")
cat("apeglm version:", as.character(packageVersion("apeglm")), "\n")
cat("TCGAbiolinks version:",
    as.character(packageVersion("TCGAbiolinks")), "\n")