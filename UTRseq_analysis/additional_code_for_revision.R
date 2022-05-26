

# try normalization with ERCCs

y <- DGEList(counts=count_table, genes=row.names(count_table), group=sampleCondition$group)
y <- calcNormFactors(y)
count_cpm <- cpm(y)

# get spike-in gene names
spikes <- rownames(count_table)[grep("^ERCC", row.names(count_table))]

# filter for only genes expressed at cpm > 2 in at least 10 samples. Also removed mitochondrial genes and spike-ins
count_table_fil <- count_table[which(rowSums(count_cpm > 2) >= 10), ]
count_table_fil <- count_table_fil[!row.names(count_table_fil) %in% chrMgenes,]
#count_table_fil <- count_table_fil[!row.names(count_table_fil) %in% spikes,]

# construct edger expression set again with filtered data
y_fil <- DGEList(counts=count_table_fil, genes=row.names(count_table_fil), group=sampleCondition$age)
y_fil <- calcNormFactors(y_fil)

# get log CPM counts
logCPM <- cpm(y_fil, log=TRUE, prior.count=1)


# Correct batch effect/unwanted variables using RUV seq 

row.names(sampleCondition) <- colnames(count_table_fil)
set <- newSeqExpressionSet(as.matrix(count_table_fil),
                           phenoData = sampleCondition)
# upper quartile normalzation
set <- betweenLaneNormalization(set, which="upper")

# normalization with spike-ins
set1 <- RUVg(set, spikes[spikes %in% row.names(count_table_fil)], k=2)
logCPMc <- log2(normCounts(set1)+1)


# run qPRC comparison
# get qPRC results
all_data_normed_info_noC <- readRDS("datasets/fluidigm_all_data_normed_info_noC.rds")
# fix gene names
all_data_normed_info_noC$Primer[all_data_normed_info_noC$Primer == "Rab7l1"] <- "Rab29"
fluidigm_genes <- unique(all_data_normed_info_noC$Primer)
all_genes <- fluidigm_genes

# reformat pcr results table
PCR_puberty <- subset(all_data_normed_info_noC, tissue == "P") %>%
  mutate(Sample = gsub("_", "", Sample)) %>%
  group_by(Sample, Primer, age, gender) %>%
  summarise(Ct= mean(Ct)) %>%
  as.data.frame()

# compare pcr results and utrseq results (before and after normalization, only post-normalization data shown)
exp_compare_cor_values <- compare_qPCR(logCPM, logCPMc, PCR_puberty, outname, method = "spearman", ifgenplot = F)

compare_sum <- compare_qPCR_summaryplot(exp_compare_cor_values$cor_values)

exp_compare_cor_values$plots_after_correction$Pd37M4

compare_sum[[1]]
ggsave("~/Dropbox/wilson_lab/active_manuscripts/pituitary_paper/RC_revision/revision_figures/normalization_with_spikeins_qPCR_compare.pdf", width = 7)


# PCA ---------------------------------------------------------------------


# PCA
logCPMc_PCAplot <- gen_pca_plot(logCPMc, npc=20, sample_label=F, label =F)
logCPMc_PCAplot[[1]]
ggsave("~/Dropbox/wilson_lab/active_manuscripts/pituitary_paper/RC_revision/revision_figures/normalization_with_spikeins_PCA.pdf", width = 7)

# spike-in perc
spike_perc <- round(100*(colSums(count_table[spikes,])/colSums(count_table)),2)
spike_perc_df <- sampleCondition
spike_perc_df$spike_perc <- spike_perc[row.names(spike_perc_df)]

ggplot(spike_perc_df) +
  geom_col(aes(x = rep, y = spike_perc, fill=sex))+
  facet_grid(age ~ sex) +
  scale_fill_manual(values = c("M" = "steelblue", "F" = "tomato")) +
  theme_bw() +
  ylab("Percent of reads mapping to spike-ins")

ggsave("~/Dropbox/wilson_lab/active_manuscripts/pituitary_paper/RC_revision/revision_figures/spike_in_reads_perc.pdf", width = 7, height = 5)


# ERCC concentration correlation ------------------------------------------


# calculate ERCC correlations
ERCC_con_v <- readRDS("/mnt/mdwilson/huayun/ganno/ERCC_concentrations_2ul.rds")

cpm_spikes <- count_cpm[spikes,]

# compare ERCC amount and CPM
cpm_spikes_df <- melt(log2(cpm_spikes+1))
cpm_spikes_df$ERCC_amount <- ERCC_con_v[as.character(cpm_spikes_df[,1])]

get_cor <- function(sample){
  tempdata <- subset(cpm_spikes_df, Var2 == sample & value > 0)
  model <- lm(value~ERCC_amount, tempdata)
  slp <- coef(model)[2]
  cor <- with(tempdata, cor(value, ERCC_amount))
  return(list(cor, slp, nrow(tempdata)))
}

ERCCamount_cors <- unlist(as.matrix(sapply(unique(cpm_spikes_df$Var2), get_cor))[1,])
ERCCamount_slopes <- unlist(as.matrix(sapply(unique(cpm_spikes_df$Var2), get_cor))[2,])
ERCCamount_ndetected <- unlist(as.matrix(sapply(unique(cpm_spikes_df$Var2), get_cor))[3,])
names(ERCCamount_ndetected) <- names(ERCCamount_cors) <- names(ERCCamount_slopes) <- unique(cpm_spikes_df$Var2)
sort(ERCCamount_cors)
sort(ERCCamount_slopes)

spike_perc_df$ERCC_coef <- ERCCamount_slopes[row.names(spike_perc_df)]
ggplot(spike_perc_df) +
  geom_col(aes(x = rep, y = ERCC_coef, fill=sex))+
  facet_grid(age ~ sex) +
  scale_fill_manual(values = c("M" = "steelblue", "F" = "tomato")) +
  theme_bw() +
  ylab("linear regression coef (spike in read counts vs. concentration)")


# testing between pd27 samples --------------------------------------------

# method1 using all samples
sample_info_d27 <- pData(set1)
sample_info_d27 <- sample_info_d27 %>% 
  unite(sample, group, rep, remove = F) %>% 
  rowwise() %>% 
  mutate(puberty = case_when(sample %in% c("d27F_2","d27F_4","d27F_6") ~"yes",
                             sample %in% c("d27F_1","d27F_3","d27F_5") ~ "no",
                             TRUE ~ "NA"))
design_d27 <- model.matrix(~0+puberty+W_1+W_2, sample_info_d27)

y_fil_d27 <- estimateDisp(y_fil, robust = TRUE, design = design_d27)

fit_d27 <- glmQLFit(y_fil_d27, design=design_d27, robust = TRUE)

qlf_d27 <- glmQLFTest(fit_d27, contrast = makeContrasts(d27_vo = pubertyyes - pubertyno, levels = design_d27))
topTg_d27 <- topTags(qlf_d27, n=nrow(y_fil$counts))

de_d27 <- as.data.frame(topTg_d27[[1]])
de_d27$genename <- genename[de_d27$genes]

# method2 using only d27 samples

sample_info_d27_only <- subset(sample_info_d27, age == "d27" & sex == "F")
design_d27_only <- model.matrix(~W_1+W_2+puberty, sample_info_d27_only)

y_fil_d27_only <- DGEList(counts=count_table_fil[, grepl("Pd27F", colnames(count_table_fil))], genes=row.names(count_table_fil))
y_fil_d27_only <- calcNormFactors(y_fil_d27_only)

y_fil_d27_only <- estimateDisp(y_fil_d27_only, robust = TRUE, design = design_d27_only)
fit_d27_only <- glmQLFit(y_fil_d27_only, design=design_d27_only, robust = TRUE)

qlf_d27_only <- glmQLFTest(fit_d27_only, coef = 4)
topTg_d27_only <- topTags(qlf_d27_only, n=nrow(y_fil$counts))
de_d27_only <- as.data.frame(topTg_d27_only[[1]])
de_d27_only$genename <- genename[de_d27_only$genes]

up <- subset(de_d27_only, logFC > log2(1.5) & FDR < 0.05)$genename
write_lines(up, "~/Desktop/d27_up_genes.txt")

# heatmap showing PD27 DE genes when using only PD27 samples. PD27F1, PD37F4 and PD32M3 cluster together. Some of the genes correspond to M7 in coexpression analysis, mapping to posterior, connective tissue, vascular etc. If showing only PD27 samples, its clear that DE mostly driven by PD27F1

pd27 <- plot_gene_hm(subset(de_d27_only, abs(logFC) > log2(1.5) & FDR < 0.05)$genename, 
                  #logCPMc[, grepl("Pd27F", colnames(logCPMc))], 
                  logCPMc,
                  #clust_cols = F, 
                   # row_cut = 1, 
                  anno_colors = c(anno_cols, anno_cols_sex), 
                 # anno_rows = sex_gene_df[, c("d12","d22","d27","d32","d37")],
                  colnames = T,
                  # fontsize_row = 8,
                  colors = hmcol_yp, 
                  sampleCon = sampleCondition)
dev.off()
print(pd27)


# explore neuron genes ----------------------------------------------------

neuron_part_genes <- unique(c(to_mouse(subset(enrich_list_direction$d27_sex_up, `term.name` == "neuron part")$intersection),to_mouse(subset(enrich_list_direction$d32_sex_up, `term.name` == "neuron part")$intersection)))
                   
neuron_part_genes <- intersect(to_mouse(subset(enrich_list_direction$d27_sex_up, `term.name` == "neuron part")$intersection),to_mouse(subset(enrich_list_direction$d32_sex_up, `term.name` == "neuron part")$intersection))

write_lines(neuron_part_genes, "~/Desktop/neuron_part_genes.txt")    
write_lines(names(genename[genename %in% neuron_part_genes]), "~/Desktop/neuron_part_genes_ids.txt")  



# compare DEG with UTRseq and qPCR ----------------------------------------

qpcr <- read.xlsx("/mnt/mdwilson/huayun/puberty/pit_rna/scripts/pituitary_transcriptome_analyses/UTRseq_analysis/datasets/fluidigm_gender_diff_t_test_2016-10-03.xlsx")

de_sig_list_wPCR <- lapply(de_sig_list, function(x) subset(x, genename %in% unique(PCR_puberty$Primer)))

de_sig_list_wPCR <- bind_rows(de_sig_list_wPCR, .id = "age")

de_sig_list_wPCR <- de_sig_list_wPCR %>% 
  mutate(age = gsub("d|_sex","", age)) %>% 
  left_join(subset(qpcr, tissue == "P"), by = c("genename"="Primer", "age"))

ggplot(de_sig_list_wPCR) +
  geom_bar(aes(x=age, fill=p_adj < 0.05))

pcr_puberty_df <- PCR_puberty %>% 
  mutate(sample= paste0("pcr_", Sample)) %>% 
  pivot_wider(names_from = "sample", values_from = "Ct", id_cols = "Primer") %>% 
  column_to_rownames("Primer")

pcr_sig_df <- subset(qpcr, tissue == "P" & p_adj < 0.05)[,c("Primer","age")] %>% 
  mutate(age = paste0("d",age)) %>% 
  setNames(c("value","age"))

factor_plot(pcr_puberty_df, gene_name = "Fshb", name = T, convert_genename = F, ylab = "delta Ct", sig_df = pcr_sig_df)

factor_plot(logCPMc, gene_name = c("Fshb"), sig_df=sig_df)

dis_genes <- unique(subset(de_sig_list_wPCR, p_adj > 0.05)$genename)

p1 <- factor_plot(pcr_puberty_df, gene_name = dis_genes, name = T, convert_genename = F, ylab = "delta Ct", sig_df = pcr_sig_df, ncols = 2) + ggtitle("qPCR")
p2 <- factor_plot(logCPMc, gene_name = dis_genes, sig_df=sig_df, ncols = 2) + ggtitle("UTRseq")

p1 + p2 +  plot_annotation(tag_levels = 'A')
ggsave("~/Dropbox/wilson_lab/active_manuscripts/pituitary_paper/RC_revision/revision_figures/utrseq_qpcr_disagree_genes.pdf", width = 12, height = 8)
ggsave("~/Dropbox/wilson_lab/active_manuscripts/pituitary_paper/RC_revision/revision_figures/utrseq_qpcr_disagree_genes.png", width = 12, height = 8, dpi = 150, units = "in")

