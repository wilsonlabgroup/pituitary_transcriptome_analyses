#' ---
#' title: "miRNA-gene correlation analysis"
#' author: "Cadia Chan"
#' affiliation: "SickKids Research Institute & University of Toronto"
#' date: "January 18, 2022"
#' output:
#'  html_document:
#'    code_folding: hide
#'    toc: true
#'    toc_float: true
#' ---
#' This script correlates miRNA and mRNA expression (log2(normalized CPM + 1)) from matched pituitary samples.
#' miRNA-gene pairs are identified from miRTarBase (experimentally validated pairs) and TargetScan (computationally predicted pairs).
#' Spearman's correlation coefficient (rho) is calculated for each pair.
#' A miRNA-gene pair network is generated for sex comparisons.
#' Genes implicated in pituitary disease, pubertal disorders are curated from GWAS and genetic studies.
#' Connection to Cytoscape is required to visualize the networks.
#'
#' To run: Set working directory "To Source File Location".
#'

#' ## Updates
#' * 2020/06/12 -> added in DE miRNA aspect
#' * 2020/09/29 -> removed heatmap creation of non DE miRNA and non DE genes, added in Fisher's exact test to calculated enrichment for each DE pair
#' * 2021/01/24 -> Updated to rmd format
#' * 2021/03/09 -> Added novel miRNA target correlation separately
#' * 2021/07/07 -> Removed age comparison correlations, only have sex-biased comparison correlations now
#' * 2022/01/18 -> Removed novel miRNA target correlation, removed example gene/miRNA plots, commented out line calling for Cytoscape dependency (makes interaction network)
#' 
#' To run: Set Working Directory to Source File Location

#' ## Libraries and source scripts
#+ load_library, warning = F, message = F
library(dplyr)
library(stringr)
library(edgeR)
library(gProfileR)
library(openxlsx)
library(pheatmap)
library(reshape2)
library(RColorBrewer)
library(ggplot2)
library(LaCroixColoR)
library(ggpubr)
library(knitr)
library(scales)
library(RCy3)
library(colorspace)
library(EDASeq)
library(biomaRt)

source("pit_mirna_mrna_correlation_functions_2021-07-07_rmd.R")

#' ## Load in miRNA and UTR counts
#' Similar to miRNA DE analysis, low count miRNAs are filtered out (CPM > 2 in 10+ samples are retained)
#' miRNA and mRNA normalized counts (not CPM) are log2-transformed and used for downstream analysis.
#' Only matching samples (n=41) in both miRNA and mRNA profiling are used (ie. PD37 is removed from mRNA as this age is not sequenced for miRNA expression).
#+ load_counts, warning = F, message = F
#### Load in miRNA and mRNA expression data ####
mirna_obj <- readRDS("input_files/RUV_corrected_pit_srna-seq_counts_combined.rds")
utr_obj <- readRDS("input_files/pit_utr_2019__RUV_k2_set1_2019-07-03.rds")

cpmcounts <- cpm(normCounts(mirna_obj))
keepcounts <- which(apply(cpmcounts, 1, function(x) length(which(x > 2)) >= 10))
mirna_counts <- as.data.frame(log2(normCounts(mirna_obj)[keepcounts, ] + 1))

utr_counts <- as.data.frame(log2(normCounts(utr_obj) + 1))

mirna_meta <- pData(mirna_obj) %>% mutate(sample = paste0("P", group, rep))
utr_meta <- pData(utr_obj) %>% mutate(sample = paste0("P", group, rep))

colnames(mirna_counts) <- mirna_meta$sample
colnames(utr_counts) <- utr_meta$sample

use_samples <- intersect(mirna_meta$sample, utr_meta$sample)

mirna_counts <- dplyr::select(mirna_counts, one_of(use_samples))
utr_counts <- dplyr::select(utr_counts, one_of(use_samples))

#' ## Curate miRNA-gene pairs
#' TargetScanMouse (v7.2) is used to curate computationally predicted miRNA-gene pairs (cumulative weighted context score < -0.1). Quantile of score shows 50% is -0.143.
#' miRTarBase (Release 8.0) is used to curate experimentally validated miRNA-gene pairs. This step may take a while.
#+ get_gene_pairs, warning = F, message = F
## Curate gene pairs (miRTarBase + TargetScan) ####
targetscan <- read.table("input_files/mmu_TargetScan_Mouse_7.2_Summary_Counts.default_predictions.txt", sep = "\t", header = T)
quantile(targetscan$Cumulative.weighted.context...score)  
targetscan <- filter(targetscan, Cumulative.weighted.context...score < -0.1) %>%
  dplyr::select(c(14, 2)) %>% mutate(database = "targetscan") # quantile of score shows 50% is -0.143
targetscan$Gene.Symbol <- str_to_sentence(targetscan$Gene.Symbol) # Make all genenames in sentence case
colnames(targetscan) <- c("mirna", "gene", "database")
targetscan <- unique(targetscan)

mirtarbase <- read.csv("input_files/mmu_MTI.txt", sep = "\t", header = T)[, c(2,4)] %>%
  mutate(database = "mirtarbase")
mirtarbase$Target.Gene <- str_to_sentence(mirtarbase$Target.Gene)
colnames(mirtarbase) <- c("mirna", "gene", "database")

mirnagene_pairs <- rbind(unique(targetscan), unique(mirtarbase)) %>%
  mutate(pair = paste0(mirna, "_", gene))
mirnagene_pairs <- aggregate(mirnagene_pairs$database, by = list(mirnagene_pairs$pair, mirnagene_pairs$mirna, mirnagene_pairs$gene), toString)
colnames(mirnagene_pairs) <- c("pair", "mirna", "gene", "database")
kable(mirnagene_pairs[1:5,], caption = "5 example miRNA-gene pairs curated from TargetScan and miRTarBase.")

#' ## Calculate Spearman's correlation for miRNA-gene pairs
#' Curated miRNA-gene pairs are first filtered for miRNAs and genes expressed in our samples.
#' Spearman's correlation coefficient (rho) and pvalue (pval) are then calculated for each pair.
#' FDR-adjusted p-value is then calculated to account for multiple testing.
#' "Positive" (rho > 0) and "negative" pairs (rho < 0) are filtered for (FDR < 0.1).
#+ calc_corr, message = F, warning = F
## Intersect with expressed genes and mirnas ####
use_pairs <- filter(mirnagene_pairs, gene %in% rownames(utr_counts)) %>%
  filter(mirna %in% rownames(mirna_counts))
use_pairs$database <- factor(ifelse(grepl(",", use_pairs$database), "both", use_pairs$database))

get_spearman_rho <- function(pair_idx, mirnadf, utrdf, pairdf) {
  mirna_expr <- as.matrix(mirnadf[as.character(pairdf[pair_idx, "mirna"]), ])
  utr_expr <- as.matrix(utrdf[as.character(pairdf[pair_idx, "gene"]), ])
  spearman_rho <- cor.test(mirna_expr, utr_expr, method = "spearman", exact = T)$estimate["rho"]
  spearman_pval <- cor.test(mirna_expr, utr_expr, method = "spearman", exact = T)$p.val
  # print(paste0(pair_idx, "_", spearman_rho))
  
  colnames(mirna_expr) <- paste0(colnames(mirna_expr), "_mirna")
  return(data.frame(pair_comb = paste0(as.character(pairdf[pair_idx, "pair"]), "_", spearman_rho, "_", spearman_pval)))
}

cor_pairs <- lapply(1:nrow(use_pairs), function(x) get_spearman_rho(x, mirna_counts, utr_counts, use_pairs))

mirna_counts_use <- mirna_counts
colnames(mirna_counts_use) <- paste0(colnames(mirna_counts_use), "_mirna")
mirna_counts_use <- mutate(mirna_counts_use, mirna = rownames(mirna_counts))

utr_counts_use <- mutate(utr_counts, gene = rownames(utr_counts))
# cor_pairs <- bind_rows(cor_pairs) # memory error
cor_pairs <- tibble(as.data.frame(unlist(cor_pairs)))
colnames(cor_pairs) <- c("pair_comb")
cor_pairs <- mutate(cor_pairs, mirna = sapply(pair_comb, function(x) strsplit(as.character(x), "_")[[1]][1]),
                    gene = sapply(pair_comb, function(x) strsplit(as.character(x), "_")[[1]][2]),
                    rho = sapply(pair_comb, function(x) strsplit(as.character(x), "_")[[1]][3]),
                    pval = sapply(pair_comb, function(x) strsplit(as.character(x), "_")[[1]][4])) %>%
  mutate(pair = paste0(mirna, "_", gene)) %>%
  merge(., mirna_counts_use, by = "mirna") %>%
  merge(., utr_counts_use, by = "gene") %>%
  merge(., dplyr::select(mirnagene_pairs, pair, database), by = "pair") %>%
  mutate(fdr = p.adjust(pval, "fdr"))
cor_pairs <- dplyr::select(cor_pairs, -pair_comb)

cor_pairs_df <- cor_pairs[, c(1, 4:5, 89, 3, 2, 88, 47:87, 6:46)]
cor_pairs_df$rho <- as.numeric(cor_pairs_df$rho)
cor_pairs_df$pval <- as.numeric(cor_pairs_df$pval)
print(paste0("Number of miRNAs included in analysis: ", length(unique(cor_pairs_df$mirna))))
print(paste0("Number of genes included in analysis: ", length(unique(cor_pairs_df$gene))))

dir.create("output_files")
saveRDS(cor_pairs_df, paste0("output_files/pit_mirna_mrna_pairs_w_correlation.rds"))

pos_pairs <- filter(cor_pairs_df, fdr < 0.1 & rho > 0)
neg_pairs <- filter(cor_pairs_df, fdr < 0.1 & rho < 0)

print(paste0("Number of total pairs: ", sum(nrow(pos_pairs), nrow(neg_pairs))))
print(paste0("Number of positive pairs: ", nrow(pos_pairs)))
print(paste0("Number of negative pairs: ", nrow(neg_pairs)))
print(paste0("Positive correlated pairs rho range from: ", min(pos_pairs$rho), " to ", max(pos_pairs$rho)))
print(paste0("Negative correlated pairs rho range from: -", min(abs(neg_pairs$rho)), " to -", max(abs(neg_pairs$rho))))

#' ## Make miRNA-gene interaction networks
#' The direction of age-biased and sex-biased genes and miRNAs are taken into account to further narrow down miRNA-gene pairs of interest.
#' Only negatively correlated miRNA-gene pairs are assessed.
#' Interaction networks generated are divided into (A) up-regulated genes with down-regulated miRNAs and (B) down-regulated genes with up-regulated miRNAs.
#' Puberty-related and pituitary disease genes are also highlighted in miRNA-gene pairs of interest.
#'
#' ### DE genes and miRNAs are loaded in
#+ load_de_genes_mirnas, message = F, warning = F
#### Load in DE genes ####
DE_utr <- readRDS("input_files/pit_utr_2019_de_result_list_2019-07-03.rds") 
DE_utr <- DE_utr[-(grep("vs", names(DE_utr)))] # DE genes from sex and consecutive age comparisons

# Filter for all up-regulated DE genes (consecutive age and sex comparisons)
DE_utr_up <- lapply(DE_utr, function(x) x[x$FDR < 0.05 & x$logFC > log2(1.5), ])
DE_utr_up <- DE_utr_up[lapply(DE_utr_up, nrow) > 0] # Filter out comparisons with 0 DE genes
names(DE_utr_up) <- paste0(names(DE_utr_up), "_up") 
DE_utr_up <- lapply(names(DE_utr_up), function(x) mutate(DE_utr_up[[x]], comparison = x))

# Filter for all down-regulated DE genes (consecutive age and sex comparisons)
DE_utr_down <- lapply(DE_utr, function(x) x[x$FDR < 0.05 & x$logFC < -log2(1.5), ])
DE_utr_down <- DE_utr_down[lapply(DE_utr_down, nrow) > 0] # Filter out comparisons with 0 DE genes
names(DE_utr_down) <- paste0(names(DE_utr_down), "_down")
DE_utr_down <- lapply(names(DE_utr_down), function(x) mutate(DE_utr_down[[x]], comparison = x))

#### Load in DE miRNAs ####
DE_mirna <- readRDS("input_files/mirna_pit_all_edgeR.RDS")
DE_mirna <- DE_mirna[c(1:7, 11:13)] # Select consecutive age comparisons and sex comparisons
names(DE_mirna) <- c("d12_sex", "d22_sex", "d27_sex", "d32_sex", "d22_M", "d27_M", "d32_M", "d22_F", "d27_F", "d32_F")

# Filter for all up-regulated DE miRNAs (consecutive age and sex comparisons)
DE_mirna_up <-  lapply(DE_mirna, function(x) x[x$FDR < 0.05 & x$logFC > log2(1.5), ])
DE_mirna_up <- DE_mirna_up[lapply(DE_mirna_up, nrow) > 0] # Filter out comparisons with 0 DE miRNAs
names(DE_mirna_up) <- paste0(names(DE_mirna_up), "_up")
DE_mirna_up <- lapply(names(DE_mirna_up), function(x) mutate(DE_mirna_up[[x]], comparison = x))

# Filter for all down-regulated DE miRNAs (consecutive age and sex comparisons)
DE_mirna_down <-  lapply(DE_mirna, function(x) x[x$FDR < 0.05 & x$logFC < -log2(1.5), ])
DE_mirna_down <- DE_mirna_down[lapply(DE_mirna_down, nrow) > 0]
names(DE_mirna_down) <- paste0(names(DE_mirna_down), "_down")
DE_mirna_down <- lapply(names(DE_mirna_down), function(x) mutate(DE_mirna_down[[x]], comparison = x))

#' ### Curated genes are loaded in
#' Puberty-related genes are curated from age at menarche (AAM) GWAS and voice breaking (VB) GWAS (Perry 2014, Day 2017, Hollis 2020) and gene implicated in idiopathic hypogonadotrophic hypogonadism (IHH) and Kallmann syndrome.
#' Pituitary disease genes are curated from combined pituitary hormone deficiency (CPHD) (Fang 2016), pituitary adenoma (PA) (Ye 2015, Hauser 2019) and hypopituitarism (Kurtoglu 2019) genetic studies.
#+ load_puberty_pituitary_genes, message = F, warning = F
#### Load in puberty and pituitary genes ####
# pub_genes <- read.table("input_files/GWAS_genes_nohiC_bloodeQTL_w_disease_2020-09-28.txt", sep = "\t", header = T)
# pub_genes <- filter(pub_genes, source != "Day2017_hiC", source != "Day2017_blood_eQTL", source != "Day2017_hiC;Day2017_blood_eQTL" )
# pub_genes$source <- gsub(";Day2017_hiC", "", pub_genes$source)
# pub_genes$source <- gsub("Day2017_hiC;", "", pub_genes$source)
# pub_genes$source <- gsub(";Day2017_blood_eQTL", "", pub_genes$source)
# pub_genes$source <- gsub("Day2017_blood_eQTL;", "", pub_genes$source)
# pit_genes <- read.table("input_files/pituitary_disease_GWAS_VB_gene_list_2020-09-22.txt", sep = "\t", header = T)
# 
# pub_pit_genes <- bind_rows(pub_genes[,-2], pit_genes)
# pub_pit_genes <- aggregate(pub_pit_genes$source, by = list(pub_pit_genes$MGI.symbol), toString)
# pub_pit_genes$x <- gsub(", ", ";", pub_pit_genes$x)
# pub_pit_genes$x <- gsub(",", ";", pub_pit_genes$x)
# colnames(pub_pit_genes) <- c("MGI.symbol", "source")
# 
# pub_pit_genes <- mutate(pub_pit_genes, source_group = gsub("IHH/Kallmann", "puberty_disease", source),
#                         source_group = gsub("Day2017_nearest", "puberty_GWAS", source_group),
#                         source_group = gsub("Perry2014", "puberty_GWAS", source_group),
#                         source_group = gsub("Hollis2020_GWAS_VB_FH", "puberty_GWAS", source_group),
#                         source_group = gsub("Day2015_GWAS_VB", "puberty_GWAS", source_group),
#                         source_group = gsub("Kurtoglu2019_hypopituitarism", "pituitary_disease", source_group),
#                         source_group = gsub("Fang2016_CPHD", "pituitary_disease", source_group),
#                         source_group = gsub("Hauser2019_PA", "pituitary_disease", source_group),
#                         source_group = gsub("Ye2015_PA_GWAS", "pituitary_disease_GWAS", source_group))
# 
# get_uniq_source <- function(use_string) {
#   trans_string <- unique(strsplit(use_string, ";")[[1]])
#   return(gsub(", ", ";", toString(trans_string)))
# }
# 
# pub_pit_genes$source_group <- sapply(pub_pit_genes$source_group,function(x) get_uniq_source(x))
# 
# write.table(pub_pit_genes, "pituitary_puberty_genes_combined_2020-09-28.txt", sep = "\t", quote = F, col.names = T, row.names = F)

pub_pit_genes <- read.table("input_files/pituitary_puberty_genes_combined_2020-09-28.txt", header = T, sep = "\t")
colnames(pub_pit_genes)[1] <- "gene"
kable(pub_pit_genes[1:5,], caption = "5 example genes from a curated list of puberty-related and pituitary disease genes.")

#' ### Make miRNA-gene interaction networks
#' #### Up-regulated genes and  down-regulated miRNAs from sex-biased comparisons
#+ up_genes_down_mirnas, warning = F, message = F, fig.width = 10, fig.height = 10, fig.cap = "Interaction heatmap of up-regulated genes and  down-regulated miRNAs from consecutive age-biased and sex-biased comparisons."
# gene up/ miRNA down ####
DE_up_df <- bind_rows(DE_utr_up)
DE_up_df_cut <- aggregate(DE_up_df$comparison, by = list(DE_up_df$genename), toString)
colnames(DE_up_df_cut) <- c("gene", "comparison")

DE_mirna_down_df <- bind_rows(DE_mirna_down)
DE_mirna_down_df <- aggregate(DE_mirna_down_df$comparison, by = list(DE_mirna_down_df$genes), toString)
colnames(DE_mirna_down_df) <- c("mirna", "comparison_mirna")

# Filter for negatively correlated pairs that contain an upregulated gene and downregulated miRNA
neg_DE_up_pairs <- filter(neg_pairs, gene %in% DE_up_df_cut$gene) %>%
  merge(., DE_up_df_cut, by.all = "gene") %>%
  full_join(., DE_mirna_down_df, by = "mirna") %>%
  mutate(mirna = gsub("mmu-", "", mirna), pair = gsub("mmu-", "", pair))

# Filter for negatively correlated pairs that only contain the downregulated miRNA
neg_DE_up_pairs_nosel <- filter(neg_pairs, mirna %in% DE_mirna_down_df$mirna) %>%
  merge(., DE_mirna_down_df, by.all = "mirna") %>%
  mutate(mirna = gsub("mmu-", "", mirna), pair = gsub("mmu-", "", pair))



#' #### Down-regulated genes and up-regulated miRNAs from sex-biased comparisons
#+ down_genes_up_mirnas, warning = F, message = F, fig.width = 10, fig.height = 25, fig.cap = "Interaction heatmap of down-regulated genes and up-regulated miRNAs from consecutive age-biased and sex-biased comparisons."
# gene down/ miRNA up ####
DE_down_df <- bind_rows(DE_utr_down)
DE_down_df_cut <- aggregate(DE_down_df$comparison, by = list(DE_down_df$genename), toString)
colnames(DE_down_df_cut) <- c("gene", "comparison")

DE_mirna_up_df <- bind_rows(DE_mirna_up)
DE_mirna_up_df <- aggregate(DE_mirna_up_df$comparison, by = list(DE_mirna_up_df$genes), toString)
colnames(DE_mirna_up_df) <- c("mirna", "comparison_mirna")

# Filter for negatively correlated pairs that contain an upregulated gene
neg_DE_down_pairs <- filter(neg_pairs, gene %in% DE_down_df_cut$gene) %>%
  merge(., DE_down_df_cut, by.all = "gene") %>%
  full_join(., DE_mirna_up_df, by = "mirna") %>%
  mutate(mirna = gsub("mmu-", "", mirna), pair = gsub("mmu-", "", pair))

#' ### Make sex-bias interaction network
#+ sex_bias_network, message = F, warning = F, fig.width = 8, fig.height = 7, fig.cap = "Interaction heatmap of sex-biased genes and sex-biased miRNAs."
# Make sex-biased heatmap ####
neg_DE_F_htmap <- dplyr::select(neg_DE_up_pairs, pair, gene, mirna, rho, comparison, comparison_mirna) %>%
  filter(grepl("sex", comparison)) %>%
  filter(grepl("sex", comparison_mirna))
neg_DE_M_htmap <- dplyr::select(neg_DE_down_pairs, pair, gene, mirna, rho, comparison, comparison_mirna) %>%
  filter(grepl("sex", comparison)) %>%
  filter(grepl("sex", comparison_mirna))

neg_DE_sex_htmap <- rbind(neg_DE_F_htmap, neg_DE_M_htmap)
neg_DE_up_sex <- neg_DE_up_pairs[grep("sex", neg_DE_up_pairs$comparison), ]
neg_DE_down_sex <- neg_DE_down_pairs[grep("sex", neg_DE_down_pairs$comparison), ]
htmap_data <- make_sex_corr_heatmap(df = neg_DE_sex_htmap,
                      plot_title = "Sex-biased genes and miRNAs",
                      out_file = "neg_sex_gene_mirna_pairs",
                      useQuantileFilter = T,
                      genelist = pub_pit_genes,
                      plot_width = 8,
                      plot_height = 7)

# net_data <- make_sex_interaction_network(df = neg_DE_sex_htmap,
#                                     plot_title = "Sex-biased genes and miRNAs",
#                                     out_file = "neg_sex_gene_mirna_pairs",
#                                     useQuantileFilter = T,
#                                     genelist = pub_pit_genes)
tmp <- dev.off()
print(htmap_data[[1]])
tmp <- dev.off()
# kable(net_data[[1]][1:5,], caption = "5 example nodes from interaction network generated between sex-biased genes and sex-biased miRNAs.")
# kable(net_data[[2]][1:5,], caption = "5 example edges from interaction network generated between sex-biased genes and sex-biased miRNAs.")

write.xlsx(list(F_biased_neg_pairs = neg_DE_up_pairs, M_biased_neg_pairs = neg_DE_down_pairs), 
           paste0("output_files/negative_sex_DE_correlation_pairs_2022-01-18.xlsx"))

#' ## Session Info
#+ session_info, warning = F, message = F
sessionInfo()