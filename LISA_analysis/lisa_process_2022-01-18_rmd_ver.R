#' ---
#' title: "Processing Lisa output"
#' author: "Cadia Chan"
#' affiliation: "SickKids Research Institute & University of Toronto"
#' date: "January 18, 2022"
#' output:
#'  html_document:
#'    code_folding: hide
#'    toc: true
#'    toc_float: true
#' ---
#' This script processes and plots the output from Lisa (http://lisa.cistrome.org/).
#' Command line version of the full Lisa model was run:
#' TF ChIP-seq Peak-RP (regulatory potential) and ISD-RP (in silico deletion-regulatory potential)
#' for both motif and ChIP-seq” methods) using the DNase-seq and H3K27Ac ChIP-seq data
#' and 3000 genes which were randomly selected as the background gene set.
#' Genes which were female- or male-biased in 2 or more ages between PD27, PD32, PD37 were used as input.
#' Results combined from H3K27ac-ChIP-seq and DNase-seq ISD models, and TF ChIP-seq
#' peak-only models using the Cauchy combination test for the ChIP-seq model are plotted.
#' 
#' ## Update
#' * 2022/01/18: Created .Rmd version
#' 
#' To run: Set working directory "To Source File Location".
#'
#+ setup
knitr::opts_chunk$set(message = F, warning = F)

#' ## Libraries and source scripts
#+ load_library
library(ggpubr)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(knitr)

#' ## Load LISA output
#+ load_data
tfrank_F <- read.table("input/pit_sex_utr_lisa_F_bias_chipseq_cauchy_combine_dedup.csv",
                       header = T, sep = ",")
tfrank_F <- dplyr::select(tfrank_F, pval, TF) %>%
  dplyr::rename("pval_Fbias" = pval)
  
tfrank_M <- read.table("input/pit_sex_utr_lisa_M_bias_chipseq_cauchy_combine_dedup.csv",
                       header = T, sep = ",")
tfrank_M <- dplyr::select(tfrank_M, pval, TF) %>%
  dplyr::rename("pval_Mbias" = pval)

#' ## Select top TFs
#' Change in P-value is determined by (F-bias pvalue - M-bias pvalue).
#' TFs with change in -log10(P-value) +/- 2 standard deviations are indicated & labelled in scatterplot.
#+ select_tfs
tfrank_plot <- full_join(tfrank_F, tfrank_M, by = "TF") %>%
  mutate(pval_Fbias = -log10(pval_Fbias), pval_Mbias = -log10(pval_Mbias)) %>%
  mutate(pval_diff = abs(pval_Fbias - pval_Mbias)) %>%
  arrange(-(pval_diff))
kable(head(tfrank_plot), caption = "Example of LISA output")
sd2_cut <- mean(tfrank_plot$pval_diff) + 2*sd(tfrank_plot$pval_diff)
tfrank_plot <- mutate(tfrank_plot, top_diff = ifelse(pval_diff > sd2_cut, "Y", "N")) %>%
  mutate(top_lab = ifelse(top_diff == "Y", as.character(TF), NA))

tfrank_plot <- mutate(tfrank_plot, sex_dir = ifelse((pval_Fbias - pval_Mbias) > 0 & top_diff == "Y", "F_bias", "N"))
tfrank_plot <- mutate(tfrank_plot, sex_dir = ifelse((pval_Mbias - pval_Fbias) > 0 & top_diff == "Y" & sex_dir == "N", "M_bias", sex_dir))

plot_col <- c(F_bias = "tomato", M_bias = "steelblue", N = "gray50")

#+ plot_lisa, fig.width = 8, fig.height = 8, fig.cap = "Scatterplot comparing Lisa TF rankings with combined P-values"
dir.create("output/")
p <- ggplot(tfrank_plot, aes(x = pval_Fbias, y = pval_Mbias, colour = sex_dir)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(values = plot_col, name = "", label = c(F_bias = "> 2sd F-bias", M_bias = "> 2sd M-bias", N = "< 2sd")) +
  geom_text_repel(label = tfrank_plot$top_lab, min.segment.length = 0.2) +
  # geom_text_repel(label = tfrank_plot$M_lab, min.segment.length = 0.2) +
  labs(x = "-log10(pval_Fbias)", y  = "-log10(pval_Mbias)", title = "TF ChIP-seq rank no bg dedup") +
  geom_abline(color = "gray50", linetype = "dashed") +
  ylim(0, 17) +
  xlim(0, 17)+
  theme_bw() +
  theme(text = element_text(size = 14))
print(p)
ggsave("output/lisa_TF_chipseq_rank_nobg_dedup_sex_compare.pdf", p, device = "pdf", width = 9, height = 8, useDingbats = F)
write.table(tfrank_plot, "output/lisa_TF_chipseq_rank_nobg_dedup_sex_compare_table.txt", sep = "\t",
            quote = F, row.names = F, col.names = T)

#' ## R Session
#+ r_session
sessionInfo()
