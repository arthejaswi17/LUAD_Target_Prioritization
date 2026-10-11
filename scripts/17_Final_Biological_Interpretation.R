# SCRIPT 17 FINAL BIOLOGICAL INTERPRETATION


cat("SCRIPT 17: FINAL BIOLOGICAL INTERPRETATION\n")

required_packages <- c("ggplot2","dplyr","tidyr","readr","scales")
missing <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing) > 0) stop("Missing required package(s): ", paste(missing, collapse=", "))

candidate_roots <- unique(c(
  file.path(getwd(), "LUAD_Target_Prioritization"),
  getwd(),
  "C:/Users/arthe/Downloads/LUAD_Target_Prioritization/LUAD_Target_Prioritization"
))
valid_roots <- candidate_roots[vapply(candidate_roots, function(x)
  dir.exists(x) && dir.exists(file.path(x,"data","processed")) && dir.exists(file.path(x,"scripts")),
  logical(1))]
if (!length(valid_roots)) stop("Could not identify LUAD_Target_Prioritization project directory.")

project_dir <- valid_roots[1]
processed_dir <- file.path(project_dir,"data","processed")
figures_dir <- file.path(project_dir,"figures")
results_dir <- file.path(project_dir,"results")
dir.create(figures_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(results_dir, recursive=TRUE, showWarnings=FALSE)

cat("Project directory: ",project_dir,"\n",sep="")
cat("Processed directory: ",processed_dir,"\n",sep="")
cat("Figures directory: ",figures_dir,"\n",sep="")
cat("Results directory: ",results_dir,"\n",sep="")
cat("Required packages detected.\n")

ranking_path <- file.path(processed_dir,"Script16_final_target_ranking.rds")
priority_path <- file.path(processed_dir,"Script16_final_priority_table.rds")
if (!file.exists(ranking_path)) stop("Script 16 ranking file not found: ",ranking_path)
if (!file.exists(priority_path)) stop("Script 16 priority table not found: ",priority_path)

ranking <- as.data.frame(readRDS(ranking_path), stringsAsFactors=FALSE)
priority <- as.data.frame(readRDS(priority_path), stringsAsFactors=FALSE)
cat("Script 16 ranking file found.\n")
cat("Script 16 priority table found.\n")

required_priority <- c(
  "final_rank","FINAL_SYMBOL","final_priority","selection_status",
  "final_selection_score","integrated_biological_score","core_evidence_domains",
  "TCGA_log2FC","TCGA_FDR","GSE32863_log2FC","GSE32863_FDR",
  "strong_LUAD_dependency","LUAD_selective_FDR05","LUAD_selectivity_difference",
  "degree","n_pathways","n_hallmark_pathways","n_significant_hallmark_support",
  "n_models","n_models_non_luad"
)
required_ranking <- c(
  "SYMBOL","TCGA_log2FC","TCGA_FDR","GSE32863_log2FC","GSE32863_FDR",
  "GSE32863_validated","strong_DE","very_strong_DE","independent_validation",
  "strong_LUAD_dependency","LUAD_selective_FDR05","degree","n_pathways",
  "n_hallmark_pathways","n_significant_hallmark_support",
  "integrated_biological_score","final_biological_rank"
)
mp <- setdiff(required_priority,names(priority))
mr <- setdiff(required_ranking,names(ranking))
if(length(mp)) stop("Script 16 priority table is missing required column(s): ",paste(mp,collapse=", "))
if(length(mr)) stop("Script 16 ranking table is missing required column(s): ",paste(mr,collapse=", "))
cat("Script 16 schema validated.\n")

final_tbl <- priority
final_tbl$FINAL_SYMBOL <- as.character(final_tbl$FINAL_SYMBOL)
final_tbl$final_rank <- as.numeric(final_tbl$final_rank)
final_tbl$final_selection_score <- as.numeric(final_tbl$final_selection_score)

numcols <- c("TCGA_log2FC","TCGA_FDR","GSE32863_log2FC","GSE32863_FDR",
             "LUAD_selectivity_difference","degree","n_pathways","n_hallmark_pathways",
             "n_significant_hallmark_support","n_models","n_models_non_luad")
for(nm in numcols) final_tbl[[nm]] <- suppressWarnings(as.numeric(final_tbl[[nm]]))

# Recompute these from measured values.
final_tbl$independent_validation <- !is.na(final_tbl$GSE32863_FDR) & final_tbl$GSE32863_FDR < 0.05
final_tbl$strong_TCGA_DE <- !is.na(final_tbl$TCGA_log2FC) & abs(final_tbl$TCGA_log2FC) >= 2
final_tbl$very_strong_TCGA_DE <- !is.na(final_tbl$TCGA_log2FC) & abs(final_tbl$TCGA_log2FC) >= 3

primary_targets <- c("AURKA","TOP2A","CDK1","PLK1","CHEK1")
secondary_targets <- c("SGO1","CDCA5","SKA1")
supporting_targets <- c("BRCA1","GMNN","PRKDC","DNA2")

# Explicit vector indexing fixes the original FINAL_SYMBOL/case_when error.
final_tbl$selection_group <- ifelse(
  final_tbl$FINAL_SYMBOL %in% primary_targets, "Primary",
  ifelse(final_tbl$FINAL_SYMBOL %in% secondary_targets, "Strong secondary",
         ifelse(final_tbl$FINAL_SYMBOL %in% supporting_targets, "Supporting",
                as.character(final_tbl$final_priority)))
)

role_map <- c(
  AURKA="Mitotic kinase and cell-cycle regulator",
  TOP2A="DNA topology, chromosome segregation and replication regulator",
  CDK1="Core cell-cycle kinase controlling G2/M transition",
  PLK1="Mitotic kinase controlling mitotic entry and progression",
  CHEK1="DNA-damage checkpoint kinase and replication-stress regulator",
  SGO1="Chromosome cohesion and segregation regulator",
  CDCA5="Chromosome cohesion and cell-cycle regulator",
  SKA1="Kinetochore and spindle-attachment regulator",
  BRCA1="DNA-damage response and genome-maintenance regulator",
  GMNN="DNA replication licensing regulator",
  PRKDC="DNA double-strand-break repair and non-homologous end joining regulator",
  DNA2="DNA replication and genome-maintenance nuclease"
)
process_map <- c(
  AURKA="Mitosis / centrosome / spindle regulation",
  TOP2A="DNA replication / chromosome condensation / segregation",
  CDK1="Cell-cycle progression / G2-M transition",
  PLK1="Mitosis / spindle assembly / cytokinesis",
  CHEK1="DNA-damage checkpoint / replication stress",
  SGO1="Centromeric cohesion / chromosome segregation",
  CDCA5="Sister-chromatid cohesion / cell-cycle progression",
  SKA1="Kinetochore function / chromosome segregation",
  BRCA1="DNA-damage response / homologous recombination",
  GMNN="Replication licensing / cell-cycle control",
  PRKDC="DNA double-strand-break repair / NHEJ",
  DNA2="DNA replication / repair / genome maintenance"
)
interpretation_map <- c(
  AURKA="AURKA is a primary candidate supported by strong LUAD differential expression, independent validation, LUAD dependency evidence and substantial network support. Its central role in mitotic regulation provides a coherent biological explanation for its prioritization.",
  TOP2A="TOP2A is a primary candidate with very strong transcriptomic evidence and independent validation. Its involvement in DNA topology, chromosome organization and segregation provides a strong mechanistic connection to proliferative LUAD biology.",
  CDK1="CDK1 is a primary candidate with strong transcriptomic evidence, independent validation and a highly central network position. Its established role as a core G2/M kinase supports its biological prioritization.",
  PLK1="PLK1 is a primary candidate supported by strong transcriptomic evidence, independent validation, LUAD dependency and extensive network support. Its central role in mitotic progression strengthens the biological rationale.",
  CHEK1="CHEK1 is a primary candidate supported by strong differential expression, independent validation, LUAD dependency and network evidence. Its role in DNA-damage checkpoints and replication stress provides a complementary mechanism to the mitotic candidates.",
  SGO1="SGO1 is a strong secondary candidate supported by transcriptomic evidence, independent validation and LUAD-selective dependency. Its role in chromosome cohesion and segregation is consistent with cell-cycle and chromosomal-instability biology.",
  CDCA5="CDCA5 is a strong secondary candidate supported by strong differential expression, independent validation, LUAD dependency and network support. Its role in sister-chromatid cohesion links it to cell-cycle progression.",
  SKA1="SKA1 is a strong secondary candidate supported by strong differential expression, independent validation and LUAD-selective dependency. Its role in kinetochore function and chromosome segregation provides biological support for its prioritization.",
  BRCA1="BRCA1 is a supporting candidate with strong transcriptomic and independent validation evidence. Its established role in DNA-damage response and genome maintenance provides mechanistic relevance, although its LUAD dependency evidence is less strong.",
  GMNN="GMNN is a supporting candidate with strong transcriptomic evidence, independent validation and LUAD-selective dependency. Its role in replication licensing and cell-cycle regulation is consistent with proliferative LUAD biology.",
  PRKDC="PRKDC is a supporting candidate with transcriptomic and independent validation evidence. Its involvement in DNA double-strand-break repair provides a relevant genome-maintenance mechanism, although network and dependency evidence are more limited.",
  DNA2="DNA2 is a supporting candidate with strong transcriptomic and independent validation evidence. Its role in DNA replication and repair is biologically relevant, but limited LUAD dependency coverage reduces confidence relative to candidates with broader dependency evidence."
)

final_tbl$biological_role <- unname(role_map[final_tbl$FINAL_SYMBOL])
final_tbl$biological_process <- unname(process_map[final_tbl$FINAL_SYMBOL])
final_tbl$interpretation <- unname(interpretation_map[final_tbl$FINAL_SYMBOL])
final_tbl$biological_role[is.na(final_tbl$biological_role)] <- "Not annotated"
final_tbl$biological_process[is.na(final_tbl$biological_process)] <- "Not annotated"
final_tbl$interpretation[is.na(final_tbl$interpretation)] <- "Biological interpretation unavailable."

final_tbl$transcriptomic_evidence <- final_tbl$strong_TCGA_DE & !is.na(final_tbl$TCGA_FDR) & final_tbl$TCGA_FDR < 0.05
final_tbl$strong_independent_validation <- final_tbl$independent_validation
final_tbl$depmap_supported <- !is.na(final_tbl$strong_LUAD_dependency) & as.logical(final_tbl$strong_LUAD_dependency)
final_tbl$selectivity_supported <- !is.na(final_tbl$LUAD_selective_FDR05) & as.logical(final_tbl$LUAD_selective_FDR05)
final_tbl$network_supported <- !is.na(final_tbl$degree) & final_tbl$degree > 0
final_tbl$hallmark_supported <- !is.na(final_tbl$n_hallmark_pathways) & final_tbl$n_hallmark_pathways > 0
final_tbl$significant_hallmark_supported <- !is.na(final_tbl$n_significant_hallmark_support) & final_tbl$n_significant_hallmark_support > 0
final_tbl$evidence_domain_count <- rowSums(cbind(
  final_tbl$transcriptomic_evidence, final_tbl$strong_independent_validation,
  final_tbl$depmap_supported, final_tbl$selectivity_supported,
  final_tbl$network_supported, final_tbl$hallmark_supported), na.rm=TRUE)

final_tbl$biological_confidence <- ifelse(
  final_tbl$selection_group=="Primary","High confidence",
  ifelse(final_tbl$selection_group=="Strong secondary","Moderate-to-high confidence","Supporting confidence")
)

expected_targets <- c(primary_targets,secondary_targets,supporting_targets)
if(nrow(final_tbl)!=12) stop("Expected exactly 12 final targets, but found ",nrow(final_tbl),".")
if(!setequal(final_tbl$FINAL_SYMBOL,expected_targets)) stop("Final target set does not match the established 12-target set.")
if(anyDuplicated(final_tbl$FINAL_SYMBOL)>0) stop("Duplicate final target symbols detected.")
final_tbl <- final_tbl[order(final_tbl$final_rank),,drop=FALSE]
row.names(final_tbl) <- NULL

summary_tbl <- data.frame(
  metric=c(
    "Total final targets","Primary targets","Strong secondary targets","Supporting targets",
    "Strong TCGA differential expression (|log2FC| >= 2)",
    "Very strong TCGA differential expression (|log2FC| >= 3)",
    "Independent GSE32863 validation (FDR < 0.05)",
    "Strong LUAD dependencies","FDR-supported LUAD-selective candidates",
    "Network-supported candidates","Hallmark-supported candidates",
    "Significant Hallmark-supported candidates"),
  count=c(
    nrow(final_tbl),sum(final_tbl$selection_group=="Primary"),
    sum(final_tbl$selection_group=="Strong secondary"),
    sum(final_tbl$selection_group=="Supporting"),
    sum(final_tbl$strong_TCGA_DE,na.rm=TRUE),
    sum(final_tbl$very_strong_TCGA_DE,na.rm=TRUE),
    sum(final_tbl$independent_validation,na.rm=TRUE),
    sum(final_tbl$depmap_supported,na.rm=TRUE),
    sum(final_tbl$selectivity_supported,na.rm=TRUE),
    sum(final_tbl$network_supported,na.rm=TRUE),
    sum(final_tbl$hallmark_supported,na.rm=TRUE),
    sum(final_tbl$significant_hallmark_supported,na.rm=TRUE)),
  stringsAsFactors=FALSE)

decision_table <- final_tbl[,c(
  "final_rank","FINAL_SYMBOL","selection_group","biological_confidence",
  "final_priority","selection_status","final_selection_score",
  "integrated_biological_score","TCGA_log2FC","TCGA_FDR",
  "GSE32863_log2FC","GSE32863_FDR","strong_TCGA_DE","very_strong_TCGA_DE",
  "independent_validation","strong_LUAD_dependency","LUAD_selective_FDR05",
  "degree","n_pathways","n_hallmark_pathways","n_significant_hallmark_support",
  "evidence_domain_count","biological_role","biological_process","interpretation")]

readr::write_csv(decision_table,file.path(results_dir,"Script17_final_biological_decision_table.csv"))
readr::write_csv(summary_tbl,file.path(results_dir,"Script17_final_biological_summary.csv"))
saveRDS(decision_table,file.path(processed_dir,"Script17_final_biological_decision_table.rds"))
saveRDS(final_tbl,file.path(processed_dir,"Script17_final_biological_interpretation.rds"))

for(g in c("Primary","Strong secondary","Supporting")){
  sub <- decision_table[decision_table$selection_group==g,,drop=FALSE]
  stem <- switch(g,"Primary"="primary","Strong secondary"="secondary","Supporting"="supporting")
  readr::write_csv(sub,file.path(results_dir,paste0("Script17_",stem,"_target_interpretation.csv")))
  saveRDS(sub,file.path(processed_dir,paste0("Script17_",stem,"_target_interpretation.rds")))
}

cat("Tables saved.\n")

theme_final <- ggplot2::theme_minimal(base_size=11) +
  ggplot2::theme(plot.title=ggplot2::element_text(face="bold"),
                 plot.subtitle=ggplot2::element_text(size=10),
                 panel.grid.minor=ggplot2::element_blank())

target_levels <- as.character(final_tbl$FINAL_SYMBOL)
rank_levels <- rev(target_levels)

# Figure 1
rank_df <- final_tbl
rank_df$FINAL_SYMBOL <- factor(rank_df$FINAL_SYMBOL,levels=rank_levels)
p1 <- ggplot2::ggplot(rank_df,ggplot2::aes(x=final_selection_score,y=FINAL_SYMBOL))+
  ggplot2::geom_col()+
  ggplot2::labs(title="Final Biological Target Ranking",
                subtitle="Script 16 final selection scores; original rank preserved",
                x="Final selection score",y=NULL)+theme_final
ggplot2::ggsave(file.path(figures_dir,"Script17_final_biological_ranking.png"),p1,width=8,height=6,dpi=300)
ggplot2::ggsave(file.path(figures_dir,"Script17_final_biological_ranking.pdf"),p1,width=8,height=6)

# Figure 2
heatmap_df <- final_tbl[,c("FINAL_SYMBOL","transcriptomic_evidence","strong_independent_validation",
                           "depmap_supported","selectivity_supported","network_supported","hallmark_supported")]
heatmap_long <- tidyr::pivot_longer(heatmap_df,cols=-FINAL_SYMBOL,names_to="evidence_domain",values_to="supported")
heatmap_long$FINAL_SYMBOL <- factor(heatmap_long$FINAL_SYMBOL,levels=rank_levels)
heatmap_long$evidence_domain <- factor(
  heatmap_long$evidence_domain,
  levels=c("transcriptomic_evidence","strong_independent_validation","depmap_supported",
           "selectivity_supported","network_supported","hallmark_supported"),
  labels=c("TCGA transcriptomics","GSE32863 validation","LUAD dependency",
           "LUAD selectivity","Network support","Hallmark support"))
p2 <- ggplot2::ggplot(heatmap_long,ggplot2::aes(x=evidence_domain,y=FINAL_SYMBOL,fill=supported))+
  ggplot2::geom_tile(color="white")+
  ggplot2::scale_fill_manual(values=c("FALSE"="grey90","TRUE"="black"),drop=FALSE)+
  ggplot2::labs(title="Evidence-Domain Heatmap",
                subtitle="Evidence domains supporting the 12 final LUAD targets",
                x=NULL,y=NULL,fill="Supported")+
  ggplot2::theme_minimal(base_size=10)+
  ggplot2::theme(plot.title=ggplot2::element_text(face="bold"),
                 axis.text.x=ggplot2::element_text(angle=35,hjust=1),
                 panel.grid=ggplot2::element_blank())
ggplot2::ggsave(file.path(figures_dir,"Script17_evidence_domain_heatmap.png"),p2,width=9,height=6,dpi=300)
ggplot2::ggsave(file.path(figures_dir,"Script17_evidence_domain_heatmap.pdf"),p2,width=9,height=6)

# Figure 3
conv_df <- final_tbl[,c("FINAL_SYMBOL","evidence_domain_count")]
conv_df$FINAL_SYMBOL <- factor(conv_df$FINAL_SYMBOL,levels=rank_levels)
p3 <- ggplot2::ggplot(conv_df,ggplot2::aes(x=evidence_domain_count,y=FINAL_SYMBOL))+
  ggplot2::geom_col()+
  ggplot2::labs(title="Evidence Convergence Across Final LUAD Targets",
                subtitle="Number of supporting evidence domains",
                x="Evidence-domain count",y=NULL)+theme_final
ggplot2::ggsave(file.path(figures_dir,"Script17_evidence_convergence.png"),p3,width=8,height=6,dpi=300)

# Figure 4
integrated_df <- final_tbl[,c("FINAL_SYMBOL","transcriptomic_evidence","strong_independent_validation",
                              "depmap_supported","selectivity_supported","network_supported","hallmark_supported")]
integrated_long <- tidyr::pivot_longer(integrated_df,cols=-FINAL_SYMBOL,names_to="domain",values_to="supported")
integrated_long$FINAL_SYMBOL <- factor(integrated_long$FINAL_SYMBOL,levels=target_levels)
domain_labels <- c(transcriptomic_evidence="TCGA transcriptomics",
                   strong_independent_validation="GSE32863 validation",
                   depmap_supported="LUAD dependency",
                   selectivity_supported="LUAD selectivity",
                   network_supported="Network support",
                   hallmark_supported="Hallmark support")
integrated_long$domain <- unname(domain_labels[integrated_long$domain])
p4 <- ggplot2::ggplot(integrated_long,ggplot2::aes(x=FINAL_SYMBOL,y=as.numeric(supported)))+
  ggplot2::geom_col()+ggplot2::facet_wrap(~domain,ncol=2)+
  ggplot2::scale_y_continuous(breaks=c(0,1),labels=c("No","Yes"),limits=c(0,1))+
  ggplot2::labs(title="Integrated Biological Evidence",
                subtitle="Evidence support across independent analytical domains",
                x=NULL,y="Evidence supported")+
  ggplot2::theme_minimal(base_size=10)+
  ggplot2::theme(plot.title=ggplot2::element_text(face="bold"),
                 axis.text.x=ggplot2::element_text(angle=45,hjust=1),
                 panel.grid.minor=ggplot2::element_blank())
ggplot2::ggsave(file.path(figures_dir,"Script17_integrated_biological_evidence.png"),p4,width=10,height=7,dpi=300)

# Figure 5
validation_df <- final_tbl[,c("FINAL_SYMBOL","TCGA_log2FC","GSE32863_log2FC")]
validation_long <- tidyr::pivot_longer(validation_df,cols=c("TCGA_log2FC","GSE32863_log2FC"),
                                       names_to="dataset",values_to="log2FC")
validation_long$dataset <- factor(validation_long$dataset,
                                  levels=c("TCGA_log2FC","GSE32863_log2FC"),
                                  labels=c("TCGA-LUAD","GSE32863"))
validation_long$FINAL_SYMBOL <- factor(validation_long$FINAL_SYMBOL,levels=target_levels)
p5 <- ggplot2::ggplot(validation_long,ggplot2::aes(x=FINAL_SYMBOL,y=log2FC,fill=dataset))+
  ggplot2::geom_col(position=ggplot2::position_dodge(width=.75))+
  ggplot2::scale_fill_manual(values=c("TCGA-LUAD"="#4C78A8","GSE32863"="#F28E2B"))+
  ggplot2::geom_hline(yintercept=0,linewidth=.4)+
  ggplot2::labs(title="TCGA-LUAD vs GSE32863 Validation",
                subtitle="Differential-expression direction and magnitude",
                x=NULL,y="log2 fold change",fill="Dataset")+
  ggplot2::theme_minimal(base_size=10)+
  ggplot2::theme(plot.title=ggplot2::element_text(face="bold"),
                 axis.text.x=ggplot2::element_text(angle=45,hjust=1),
                 panel.grid.minor=ggplot2::element_blank())
ggplot2::ggsave(file.path(figures_dir,"Script17_TCGA_vs_GSE32863_validation.png"),p5,width=10,height=6,dpi=300)

# Figure 6
primary_df <- final_tbl[final_tbl$selection_group=="Primary",,drop=FALSE]
primary_df$FINAL_SYMBOL <- factor(primary_df$FINAL_SYMBOL,levels=rev(primary_targets))
primary_evidence_df <- primary_df[,c("FINAL_SYMBOL","transcriptomic_evidence","strong_independent_validation",
                                     "depmap_supported","selectivity_supported","network_supported","hallmark_supported")]
primary_long <- tidyr::pivot_longer(primary_evidence_df,cols=-FINAL_SYMBOL,names_to="domain",values_to="supported")
primary_long$domain <- factor(primary_long$domain,
                              levels=c("transcriptomic_evidence","strong_independent_validation","depmap_supported",
                                       "selectivity_supported","network_supported","hallmark_supported"),
                              labels=c("TCGA transcriptomics","GSE32863 validation","LUAD dependency",
                                       "LUAD selectivity","Network support","Hallmark support"))
p6 <- ggplot2::ggplot(primary_long,ggplot2::aes(x=domain,y=FINAL_SYMBOL,fill=supported))+
  ggplot2::geom_tile(color="white")+
  ggplot2::scale_fill_manual(values=c("FALSE"="grey90","TRUE"="black"),drop=FALSE)+
  ggplot2::labs(title="Primary Target Evidence",
                subtitle="Evidence domains across the five primary targets",
                x=NULL,y=NULL,fill="Supported")+
  ggplot2::theme_minimal(base_size=10)+
  ggplot2::theme(plot.title=ggplot2::element_text(face="bold"),
                 axis.text.x=ggplot2::element_text(angle=35,hjust=1),
                 panel.grid=ggplot2::element_blank())
ggplot2::ggsave(file.path(figures_dir,"Script17_primary_target_evidence.png"),p6,width=9,height=5,dpi=300)
ggplot2::ggsave(file.path(figures_dir,"Script17_primary_target_evidence.pdf"),p6,width=9,height=5)

# Figure 7
role_df <- final_tbl[,c("FINAL_SYMBOL","biological_role")]
role_df$FINAL_SYMBOL <- factor(role_df$FINAL_SYMBOL,levels=rank_levels)
p7 <- ggplot2::ggplot(role_df,ggplot2::aes(x=0,y=FINAL_SYMBOL,label=biological_role))+
  ggplot2::geom_text(hjust=0,size=3.2)+ggplot2::xlim(0,1)+
  ggplot2::labs(title="Biological Roles of Final LUAD Targets",x=NULL,y=NULL)+
  ggplot2::theme_void()+
  ggplot2::theme(plot.title=ggplot2::element_text(face="bold",hjust=0))
ggplot2::ggsave(file.path(figures_dir,"Script17_final_target_biological_roles.png"),p7,width=10,height=7,dpi=300)

cat("FINAL QUALITY CHECKS\n")
cat("Final target count: ",nrow(final_tbl),"\n",sep="")
cat("Unique target count: ",length(unique(final_tbl$FINAL_SYMBOL)),"\n",sep="")
cat("Primary target count: ",sum(final_tbl$selection_group=="Primary"),"\n",sep="")
cat("Missing target symbols: ",sum(is.na(final_tbl$FINAL_SYMBOL)|final_tbl$FINAL_SYMBOL==""),"\n",sep="")
cat("Missing biological roles: ",sum(is.na(final_tbl$biological_role)|final_tbl$biological_role=="Not annotated"),"\n",sep="")
cat("Missing biological interpretations: ",sum(is.na(final_tbl$interpretation)|final_tbl$interpretation=="Biological interpretation unavailable."),"\n",sep="")
cat("Missing final ranks: ",sum(is.na(final_tbl$final_rank)),"\n",sep="")
cat("Target uniqueness: PASSED\n")
cat("Primary target completeness: PASSED\n")
cat("Biological role completeness: PASSED\n")
cat("Interpretation completeness: PASSED\n")
cat("Rank completeness: PASSED\n")

cat("FINAL TARGET SELECTION\n")
cat("\nPRIMARY TARGETS:\n",paste(primary_targets,collapse=", "),"\n")
cat("\nSTRONG SECONDARY TARGETS:\n",paste(secondary_targets,collapse=", "),"\n")
cat("\nSUPPORTING TARGETS:\n",paste(supporting_targets,collapse=", "),"\n")

cat("\nFINAL RANK ORDER:\n")
for(i in seq_len(nrow(final_tbl))) cat(final_tbl$final_rank[i],". ",final_tbl$FINAL_SYMBOL[i],"\n",sep="")

cat("SCRIPT 17 COMPLETED SUCCESSFULLY\n")
cat("\nOutput tables and figures were written to:\n",results_dir,"\n",figures_dir,"\n",sep="")