#' ---
#' title: "Sex-biased miRNA analysis"
#' author: "Cadia Chan"
#' affiliation: "SickKids Research Institute & University of Toronto"
#' date: "January 17, 2022"
#' output:
#'  html_document:
#'    code_folding: hide
#'    toc: true
#'    toc_float: true
#' ---
#' This script creates a count matrix for miRNAs identified from sRNA-seq by miRDeep2.
#' RUV-seq is then run to remove batch effects from the data.
#' QC is performed with batch-corrected count matrix by visualizing the sample correlation heatmap and PCA plot.
#' edgeR is used to perform differential expression (DE) analysis on the batch-corrected count matrix.
#' DE analysis is visualized with volcano plots, barplots, and Upset plots.
#' Sex-biased miRNAs are identified at each profiled age (PD12, 22, 27, 32, 37).  
#' 
#' To run: Set working directory "To Source File Location".
#'
#' ## Updates
#' * 2021/03/17: Added in correlation plot comparing PD22-PD12 fold change in females vs males
#' * 2021/06/04: Added significance indicator for miRNA expression plots
#' * 2022/01/17: Removed age comparisons from analysis
#' 
#' ## Load in libraries and source scripts
#+ libraries, warning = F, message = F
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(gProfileR)
library(dplyr)
library(reshape2)
library(EDASeq)
library(RUVSeq)
library(pheatmap)
library(RColorBrewer)
library(edgeR)
library(ggfortify)
library(dendextend)
library(gridExtra)
library(colorspace)
library(LaCroixColoR)
library(UpSetR)
library(scales)
library(openxlsx)
library(scales)
library(knitr)

source("pit_srna_puberty_functions_2020-09-22.R")

#' ## Make count matrix
#' For each sample, counts are loaded in as 1) known mature miRNA counts, 2) novel miRNA counts
#' A novel miRNA key table is made: novel miRNAs are renamed here as novel1, novel2,...noveln
#' The mature sequence of each novel miRNA is kept track in the table
#' A raw count matrix is then created containing both known and novel miRNAs (rows) with each sample (columns)
#+ make_raw_countmat, warning = F, message = F
#### Make count matrix ####
dir.create("output_files")
files <- list.files("input_files/")
mature_files <- files[which(sapply(files, function(x) grepl("mature", x)))]
novel_files <- files[which(sapply(files, function(x) grepl("novel", x)))]
files <- files[c(which(files %in% mature_files), which(files %in% novel_files))]
sample_id <- unique(substr(files, 1, 10))

full_novel_mirdeep <- bind_rows(lapply(novel_files, function(x) read.csv(paste0("input_files/", x), sep = "\t", header = T) %>%
                                           filter(miRDeep2.score >= 2, rfam.alert == "-") %>%
                                         mutate(sample = x)))
novel_mirdeep <- unique( dplyr::select(full_novel_mirdeep, consensus.mature.sequence, precursor.coordinate) %>%
                           mutate(precursor.coordinate = gsub("..", "-", precursor.coordinate, fixed = T)))
                                          
novel_mirdeep <- unique(bind_rows(lapply(novel_files, function(x) read.csv(paste0("input_files/", x), sep = "\t", header = T) %>%
                                           filter(miRDeep2.score >= 2, rfam.alert == "-") %>%
                                           dplyr::select(consensus.mature.sequence, precursor.coordinate) %>%
                                           mutate(precursor.coordinate = gsub("..", "-", precursor.coordinate, fixed = T)))))
novel_mirdeep <- aggregate(novel_mirdeep$precursor.coordinate, by = list(novel_mirdeep$consensus.mature.sequence), toString) %>%
  mutate(id = paste0("novel", 1:nrow(.)))
colnames(novel_mirdeep) <- c("consensus.mature.sequence", "precursor.coordinate", "id")

novel_mirdeep_test <- unique(inner_join(novel_mirdeep,
                                dplyr::select(full_novel_mirdeep, consensus.mature.sequence,
                                              miRDeep2.score, sample), by = "consensus.mature.sequence"))

kable(head(novel_mirdeep), caption = "First 5 novel miRNAs with mature sequence and precursor coordinates")
write.table(novel_mirdeep, "output_files/novel_mirna_list.txt", sep = "\t", quote = F, col.names = T, row.names = F)

novel_counts <- lapply(novel_files, function(x) read.csv(paste0("input_files/", x), sep = "\t", header = T) %>%
                         filter(miRDeep2.score >= 2) %>%
                         select(consensus.mature.sequence, total.read.count) %>%
                         inner_join(., novel_mirdeep, by = "consensus.mature.sequence")) 
novel_counts <- lapply(novel_counts, function(x) aggregate(x$total.read.count, by = list(x$consensus.mature.sequence, x$id), sum))
names(novel_counts) <- novel_files
novel_counts <- lapply(names(novel_counts), function(i) mutate(novel_counts[[i]], filename = i, type = "novel"))
novel_df <- bind_rows(novel_counts)

mature_counts <- lapply(mature_files, function(x) read.csv(paste0("input_files/", x), sep = "\t", header = T) %>% filter(miRDeep2.score >= 2) %>%
                          select(mature.miRBase.miRNA, total.read.count) %>% mutate(id = sapply(as.character(mature.miRBase.miRNA), function(x) strsplit(x, "_")[[1]][1])))
mature_counts <- lapply(mature_counts, function(x) aggregate(x$total.read.count, by = list(x$mature.miRBase.miRNA, x$id), sum))
names(mature_counts) <- mature_files
mature_counts <- lapply(names(mature_counts), function(i) mutate(mature_counts[[i]], filename = i, type = "mature"))
mature_df <- bind_rows(mature_counts)

colnames(novel_df)[1:3] <- colnames(mature_df)[1:3] <- c("fullid", "id", "value")
count_df <- rbind(novel_df, mature_df) %>% mutate(sample = substr(filename, 1, 7), batch = sapply(filename, function(x) strsplit(x, "_")[[1]][3]))

countmat <- dcast(count_df, id~sample)
countmat[is.na(countmat)] <- 0
kable(countmat[1:5, 1:6], caption = "Example of 5 miRNAs (alphabetical order) for 5 samples")
write.table(countmat, "output_files/pit_srna-seq_unfiltered_countmat.csv", sep = "\t", quote = F, row.names = F, col.names = T)

#' ## Make sample condition table (metadata) and miRNA information table (featuredata)
#' For each sample, the sex, age, replicate number, and batch is stored parsed from filename and stored in metadata
#' For each miRNA, store the long id and short id in featuredata
#+ make_meta_fData, message = F, warning = F
#### Make metadata ####
metadata <- data.frame(fullid = sapply(novel_files, function(x) strsplit(x, "_novel")[[1]][1])) %>%
  mutate(sample = substr(fullid, 1, 7),  sex = substr(fullid, 5, 5), age = substr(fullid, 3, 4), rep = substr(fullid, 7, 7), batch = substr(fullid, 8, 10))
metadata$batch <- ifelse(metadata$batch == "_RO", "2018-12-11", "2018-01-26")
metadata <- metadata[, -1]
kable(metadata[1:5,], caption = "Metadata of 5 samples showing sex, age, replicate and batch")
write.table(metadata, "output_files/pit_srna-seq_metadata.csv", sep = "\t", quote = F, row.names = F, col.names = T)
rownames(metadata) <- metadata$sample
metadata[, c(2,3,4,5)] <- lapply(metadata[,c(2:5)], factor)

featdata <- unique(count_df[, c(1:2)]) %>% arrange(id)
rownames(featdata) <- featdata$id

#' ## Normalization using RUVseq (sRNA-seq data was sequenced in 2 batches)
#' RUVs was used to remove unwanted batch effects using replicate samples
#+ run_ruvseqs, message = F, warning = F
#### Run RUVSeq ####
rownames(countmat) <- countmat$id
countobj <- newSeqExpressionSet(as.matrix(countmat[, -1]), phenoData = AnnotatedDataFrame(metadata),
                    featureData = AnnotatedDataFrame(featdata))
print(countobj)
saveRDS(countobj, "output_files/pit_srna-seq_unfiltered_count_expressionset.rds")

metadata <- pData(countobj) %>% mutate(group = paste0("d", age, sex))
rownames(metadata) <- metadata$sample

countmat <- counts(countobj)

filter <- apply(countmat, 1, function(x) length(x[x>2])>=5)
filt_countmat <- countmat[filter, ]
mirna <- rownames(filt_countmat)
newfData <- filter(fData(countobj), id %in% rownames(filt_countmat))
rownames(newfData) <- newfData$id
set <- newSeqExpressionSet(as.matrix(filt_countmat),
                           phenoData = metadata,
                           featureData = newfData)

colors <- brewer.pal(3, "Set2")
#RLE = log-ratio of read count to median read count across sample
# plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[metadata$batch], )
# plotPCA(set, col=colors[metadata$batch], cex=1.2)

set <- betweenLaneNormalization(set, which="upper")
# plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[metadata$batch])
# plotPCA(set, col=colors[metadata$batch], cex=1.2)

design <- model.matrix(~metadata$group, data=pData(set))
y <- DGEList(counts=counts(set), group=metadata$group)
y <- calcNormFactors(y, method="upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)

## Run RUVSeqS
fit <- glmFit(y, design)
differences <- makeGroups(paste0(metadata$age, metadata$sex))
ruv_set <- RUVs(set, mirna, k=3, differences)

pdf("output_files/pit_srna-seq_ruv-seq_pca.pdf", width = 6, height = 6)
plotRLE(ruv_set, outline=FALSE, ylim=c(-2, 2), col=colors[factor(metadata$batch)])
plotPCA(ruv_set, col=colors[factor(metadata$batch)], cex=1.2)
plotPCA(ruv_set, col=colors[factor(metadata$age)], cex=1.2)
plotPCA(ruv_set, col=colors[factor(metadata$sex)], cex=1.2)
discard <- dev.off()
saveRDS(ruv_set, "output_files/RUV_corrected_pit_srna-seq_counts_combined.rds")

ruv_counts <- log(normCounts(ruv_set) + 1)
ruv_cor <- cor(ruv_counts, method = "pearson")
ruv_melt <- melt(ruv_cor, na.rm = T)
pheatmap(ruv_cor, annotation = metadata[,c(2,3, 5)],
         filename="output_files/RUV_correlation_heatmap_pit_srna-seq_combined.pdf",
         width = 9, height = 8)

ruv_dists <- dist(t(ruv_counts))
ruv_distsMatrix <- as.matrix(ruv_dists)
rownames(ruv_distsMatrix) <- colnames(ruv_counts)
colnames(ruv_distsMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Greens")) )(255)
pheatmap(t(ruv_distsMatrix),
         clustering_distance_rows=ruv_dists,
         clustering_distance_cols=ruv_dists,
         # col=colors,
         # main = plot_title,
         annotation_col = metadata[,c(2,3, 5)],
         treeheight_row = 0,
         filename="output_files/RUV_distance_heatmap_pit_srna-seq_combined.pdf",
         width = 9, height = 8)


#' ## Data QC
#' Low counts are first filtered from batch-normalized count matrix (CPM > 2 in 10 or more samples).
#' edgeR is used to calculate normalization factors to scale raw library sizes.
#' Then QC is performed to look for sample outliers by PCA and heatmap.
#+ data_qc, message = F, warning = F

#### Data QC ####
ruv_obj <- ruv_set
ruv_countmat <- normCounts(ruv_obj)
ruv_metadata <- pData(ruv_obj)

sample_sort <- order(ruv_metadata$sex)
ruv_metadata <- ruv_metadata[sample_sort, ]
ruv_metadata[, c(2:5)] <- lapply(ruv_metadata[, c(2:5)], factor)
ruv_countmat <- ruv_countmat[, sample_sort]

y <- DGEList(counts = ruv_countmat, 
             samples = ruv_metadata,
             genes = rownames(ruv_obj))
cpmcounts <- cpm(y)
keepcounts <- which(apply(cpmcounts, 1, function(x) length(which(x > 2)) >= 10))

y <- y[keepcounts, ]
y <- calcNormFactors(y)

ruv_norm_filt <- log2(y$counts + 1)
ruv_logcpm <- log2(cpm(y) + 1)


# Make filtered novel miRNA table
filt_fdata <- filter(novel_mirdeep, id %in% rownames(ruv_logcpm))
filt_fdata <- inner_join(filt_fdata, full_novel_mirdeep, by = "consensus.mature.sequence")
filt_fdata_mirdeep_score <- aggregate(filt_fdata$miRDeep2.score, by = list(filt_fdata$id), mean)
colnames(filt_fdata_mirdeep_score) <- c("id", "mean_mirdeep_score")
filt_fdata <- unique(inner_join(filt_fdata, filt_fdata_mirdeep_score, by = "id") %>%
  select(id, consensus.mature.sequence, precursor.coordinate.x, mean_mirdeep_score, rfam.alert))
kable(filt_fdata, caption = "Novel miRNAs with CPM > 2 in 10+ samples.")
write.table(filt_fdata, "output_files/novel_mirna_list_ruv_logcpm_filtered.txt", col.names = T, row.names = F, sep = "\t", quote = F)


#' ### Sample global overview
#' PCA plot and correlation heatmap are generated using log2-transformed RUV-seq-corrected counts.
#+ make_pca, fig.cap = "PCA plot of all miRNA samples", message = F, warning = F, fig.width = 5, fig.height = 7
#### Make PCA ####
prcmp <- prcomp(t(ruv_norm_filt))

sex_lab_pca <- plot_pca(ruv_norm_filt, ruv_metadata, "Pituitary miRNA - Coloured by sex", "sex", "age",
                        c("tomato", "steelblue"), label_anno = ruv_metadata$rep)

sex_pca <- plot_pca(ruv_norm_filt, ruv_metadata, "Pituitary miRNA - Coloured by sex", "sex", "age",
                    c("tomato", "steelblue"))
ggsave("output_files/mirna_pit_combined_pca_sex_coloured.pdf", sex_pca, device = "pdf", width = 5, height = 7, useDingbats = F)

pal <- brewer.pal(n = 9, "BuGn")
age_pca <- plot_pca(ruv_norm_filt, ruv_metadata, "Pituitary miRNA - Coloured by age", "age", "sex",
                    pal[c(3,5,7,9)])
batch_pca <- plot_pca(ruv_norm_filt, ruv_metadata, "Pituitary miRNA - Coloured by batch", "batch", "age",
                      c("darkmagenta", "springgreen2"))
p_all <- ggarrange(plotlist = list(sex_lab_pca, sex_pca, age_pca, batch_pca))
ggsave("output_files/mirna_pit_combined_all_pcas.pdf", p_all, device = "pdf", width = 10.5, height = 15, useDingbats = F)

print(sex_pca)

#+ make_corhtmap, fig.cap = "Heatmap showing correlation between all miRNA samples", message = F, warning = F, fig.width = 9, fig.height = 8
#### Make sample correlation heatmap ####

Var1 <- c("tomato", "steelblue")
names(Var1) <- c("F", "M")
Var2 <- pal[c(1,3,5,7)]
names(Var2) <- levels(ruv_metadata[,3])
use_anno_colors <- list(sex = Var1, age = Var2)

p <- plot_sample_dist_mat(ruv_norm_filt, ruv_metadata[, c(3, 2)], "Pituitary gland miRNA RUV log2(normCounts)", 
                          use_anno_colors)
print(p)
save_phtmap_pdf(p, "output_files/mirna_pit_ruv_log_normcounts_corr_heatmap.pdf", height = 8, width = 9)

cor_mat <- cor(ruv_norm_filt,method="pearson")
all_cor_means <- c()
ruv_metadata$group <- factor(ruv_metadata$group)
for(i in levels(ruv_metadata$group)) {
  sub_df <- rownames(ruv_metadata[ruv_metadata$group == i, ])
  sub_df <- cor_mat[sub_df, sub_df]
  sub_melt <- melt(sub_df)
  sub_melt <- sub_melt[-which(sub_melt$Var1 == sub_melt$Var2), ]
  sub_mean <- mean(sub_melt$value)
  all_cor_means <- c(all_cor_means, sub_mean)
}
print(paste0("Range of PCC: ", min(all_cor_means), " to ", max(all_cor_means)))

#' ### Expression of signature pituitary miRNAs
#' miRNAs selected are previously shown to be expressed in the pituitary gland.
#+ make_example_mirna, fig.cap = "Expression of signature pituitary miRNAs", message = F, warning = F, width = 8, height = 7
pit_expr_mirna <- data.frame("mirna" = c("mmu-miR-7a-5p",  "mmu-miR-26b-5p", "mmu-miR-212-5p", "mmu-miR-132-3p"))
rownames(pit_expr_mirna) <- pit_expr_mirna$mirna
ruv_norm_pitmirna <- ruv_norm_filt[as.character(pit_expr_mirna$mirna),]
use_anno_colors <- list(sex = Var1, age = Var2)

p <- lapply(pit_expr_mirna$mirna, function(x) mirnaplot_hhtheme(x, mirna_data = ruv_norm_filt, ylabel = "log2(normCounts)"))
p_all <- ggarrange(plotlist = p, common.legend = T, legend = "bottom")
ggsave("output_files/literature_pituitary_mirnas_exprplot.pdf", plot = p_all, width = 8, height = 7, useDingbats = F)
print(p_all)

#' ## Differential expression analysis
#' edgeR is used to perform differential expression analysis (pval < 0.05, abs(log2FC) > log2(1.5)).
#+ run_edgeR, message = F, warning = F
#### Run edgeR ####
dir.create("output_files/edgeR")
setwd("output_files/edgeR")
knitr::opts_knit$set(root.dir = 'output_files/edgeR')
countmat_order <- match(ruv_metadata$sample, colnames(filt_countmat))
y_fil <- filt_countmat[keepcounts, countmat_order]
y_fil <- DGEList(counts = y_fil, 
             samples = ruv_metadata,
             genes = rownames(y_fil))
y_fil <- calcNormFactors(y_fil)
design <- model.matrix(~0+group+W_1+W_2+W_3, y_fil$samples )

y_fil <- estimateDisp(y_fil, robust = TRUE, design = design)
fit <- glmQLFit(y_fil, design=design, robust = TRUE)

pair_contrasts <- makeContrasts(d12_sex = groupd12F - groupd12M,
                                d22_sex = groupd22F - groupd22M,
                                d27_sex = groupd27F - groupd27M,
                                d32_sex = groupd32F - groupd32M,
                                levels = design)

de_result_list <- lapply(colnames(pair_contrasts), function(x) get_edger_result(pair_contrasts, x))
names(de_result_list) <- colnames(pair_contrasts)

pval_val <- 0.05
log2FC_val <- log2(1.5)

join_lognorm <- mutate(as.data.frame(ruv_norm_filt), genes = rownames(ruv_norm_filt))
out_de_result_list <- lapply(de_result_list, function(x) inner_join(x, join_lognorm, by = "genes"))
write.xlsx(out_de_result_list, "mirna_pit_sex_all_edgeR.xlsx")
saveRDS(out_de_result_list, "mirna_pit_sex_all_edgeR.RDS")

results_df <- lapply(out_de_result_list, function(x) subset(x, FDR < pval_val & abs(logFC) > log2FC_val))
names(results_df) <- colnames(pair_contrasts)
names(results_df) <- lapply(names(results_df), function(x) gsub("_sex", " F-M", x))
write.xlsx(results_df, "mirna_pit_de_sex_sig_edgeR.xlsx")

#' ### Sex-biased miRNAs
#' Number of sex-biased miRNAs at each age is summarized as a Upset plot (Fig 4B).
#' Expression of all sex-biased miRNAs (log2-transformed counts) are plotted (Fig 4D, Supp Fig S6C).
#+ make_sex_upset, fig.cap = "Upset plot summarizing number of sex-biased miRNAs at each profiled age", message = F, warning = F, fig.height = 6, fig.width = 8
#### Sex upset ####
mirna_up <- lapply(results_df[1:4], function(x) x[x$logFC > 0,])
names(mirna_up) <- paste0(names(mirna_up), "_up")
F_bias_mirnas <- bind_rows(lapply(names(mirna_up), function(x) dplyr::select(mirna_up[[x]],genes) %>% mutate(comparison = x) %>%
                                   filter(grepl("novel", genes))))
mirna_down <- lapply(results_df[1:4], function(x) x[x$logFC < 0,])
names(mirna_down) <- paste0(names(mirna_down), "_down")
M_bias_mirnas <- bind_rows(lapply(names(mirna_down), function(x) dplyr::select(mirna_down[[x]],genes) %>% mutate(comparison = x) %>%
                                    filter(grepl("novel", genes))))
mirnas <- c(lapply(mirna_up, function(x) x[,"genes"]), lapply(mirna_down, function(x) x[,"genes"]))
mirnas <- fromListx(mirnas)

sex_upset <- upset(mirnas, order.by = "freq", main.bar.color = "black",
           sets= rev(colnames(mirnas)), keep.order = T,
           text.scale = 2)

pdf("mirna_pit_sex_biased_upset.pdf", 
    width = 6.5, height = 5, useDingbats = F)
print(sex_upset)
discard <- dev.off()
print(sex_upset)

#+ make_sex_expr, fig.cap = "Expression of all sex-biased miRNAs", message = F, warning = F, fig.height = 9, fig.width = 14
#### Make sex expression plot ###

sex_df <- lapply(de_result_list[grep("sex", names(de_result_list))], function(x) subset(x, FDR < pval_val & abs(logFC) > log2FC_val))
sex_df <- bind_rows(lapply(names(sex_df), function(x) mutate(sex_df[[x]], comparison = substr(x, 1, 3))))
colnames(sex_df)[1] <- "mirna"
sex_df <- mutate(sex_df, factor = paste0(mirna, "_", comparison))

sex_mirnaplots <- mirnafactor_plot(unique(sex_df$mirna),
                 mirna_data = ruv_norm_filt,
                ylabel = "log2(normCounts)",
                input_meta = metadata,
                de_input = sex_df)
print(sex_mirnaplots)
ggsave("mirna_pit_sex_biased_exprplot.pdf", sex_mirnaplots, device = "pdf", width = 12, height = 7, useDingbats = F)


#+ make_sex_expr_sel, fig.cap = "Expression of sex-biased miRNA chosen per age comparison", message = F, warning = F, fig.height = 4, fig.width = 6
mirna_plotlist <- mirnafactor_plot(c("novel46", "mmu-miR-342-5p", "mmu-miR-383-5p", "mmu-miR-224-5p"),
                                   mirna_data = ruv_norm_filt,
                                   ylabel = "log2(normCounts)",
                                   input_meta = metadata,
                                   de_input = filter(sex_df, mirna %in% c("novel46", "mmu-miR-342-5p", "mmu-miR-383-5p", "mmu-miR-224-5p")))
print(mirna_plotlist)
ggsave("mirna_pit_sex_biased_examples_exprplot.pdf", mirna_plotlist, device = "pdf", width = 5, height = 3.5, useDingbats = F)

#' ## Novel DE miRNAs
#' Novel miRNAs which were identified as DE in either sex or age comparisons are saved.
#' Mature sequences from this table are used with miRBase blast function.
#+ novel_DE_mirna, message = F, warning = F
novel_DE_mirnas <- rbind(F_bias_mirnas, M_bias_mirnas)
novel_DE_mirnas <- aggregate(novel_DE_mirnas$comparison, by = list(novel_DE_mirnas$genes), toString) 
colnames(novel_DE_mirnas) <- c("id", "comparison")
novel_DE_mirnas <- inner_join(novel_DE_mirnas, novel_mirdeep, by = "id")
kable(novel_DE_mirnas, caption = "novel DE miRNAs with mature sequence, precursor coordinates, and DE comparison")
write.table(novel_DE_mirnas, "novel_sex_DE_mirnas.txt", sep = "\t", quote = F, col.names = T, row.names = F)


#' ## Session Info
#+ session_info, warning = F, message = F
sessionInfo()