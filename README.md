# Lung Adenocarcinoma Target Prioritization

## About the Project

This project uses computational approaches to study gene expression in lung adenocarcinoma (LUAD), a type of lung cancer. The aim is to identify and prioritize candidate genes for further investigation using gene expression data and additional biological evidence.

## Project Workflow

- Identify differentially expressed genes using TCGA-LUAD data.
- Investigate biological pathways associated with the analysis.
- Prioritize candidate genes for further investigation.
- Examine protein-protein interaction networks.
- Compare relevant findings with the GSE32863 dataset.
- Explore DepMap data for additional evidence related to cancer cell dependency.
- Integrate the available evidence to rank candidate genes.

## Tools and Data

- **Programming language:** R
- **Main dataset:** TCGA-LUAD
- **Comparison dataset:** GSE32863
- **Additional data source:** DepMap
- **Methods:** Differential expression analysis, functional enrichment, network analysis, and evidence integration.

## Repository Structure

- `scripts/` – R scripts for the analysis workflow.
- `data/processed/` – Selected processed datasets.
- `results/` – Analysis outputs and evidence summaries.
- `figures/` – Figures produced during the project.

## Limitations

This project uses computational analyses to prioritize candidate genes. The findings do not establish that these genes are effective therapeutic targets. Further biological and experimental studies are needed to investigate their functional and therapeutic relevance.

Raw downloaded datasets are not included in this repository. Reproducing the analysis may require obtaining the original datasets and adjusting file paths for the local environment.
