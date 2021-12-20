setwd("~/Dropbox/Mike_Anna_Mark/active_manuscript/Pituitary paper/scripts/mRNA_analysis_rmd/")
source("~/Dropbox/Mike_Anna_Mark/active_manuscript/Pituitary paper/scripts/mRNA_analysis_rmd/pit_rna_analysis_functions_clean.r")
counts_file <- "datasets/pit_utrseq_gencode_vM21_refSeq201906.txt"

# Obtain gene count table and sample conditions
temp <- get_count_table(counts_file, fixname = T)
count_table <- temp[[1]]
sampleCondition <- temp[[2]]
gene_info <- temp[[3]]
rm(temp)

colnames(count_table) <- gsub("p", "P", colnames(count_table))

count_file_gencode <- "/mnt/mdwilson/huayun/puberty/pit_rna/utrseq_2019/count_GENCODE/pit_utrseq_gencodeVM21_ERCC.txt"

temp <- get_count_table(count_file_gencode, fixname = T)
count_table2 <- temp[[1]]
sampleCondition <- temp[[2]]
gene_info2 <- temp[[3]]
rm(temp)
colnames(count_table2) <- gsub("p", "P", colnames(count_table2))

ave_counts <- rowMeans(count_table)
ave_counts2 <- rowMeans(count_table2)

diff <- which(ave_counts != ave_counts2)

more <- names(which(ave_counts[diff] > ave_counts2[diff] ))
table(ave_counts2[more] == 0)

ave_counts2[more][which(ave_counts2[more]==0)]

diff_level <- sort(ave_counts[more] - ave_counts2[more], decreasing = T)
names(diff_level) <- genename[names(diff_level)]

diff_level_fold <- sort((ave_counts[more] - ave_counts2[more])/ave_counts2[more], decreasing = T)
names(diff_level_fold) <- genename[names(diff_level_fold)]

novel_end <- unique(subset(novel_utrs_all, type == "intergenic")$gene_id)
novel_intron <- unique(subset(novel_utrs_all, type == "intronic")$gene_id)

novel_genes_all <- as.character(unique(novel_utrs_all$gene_id))

subset(novel_utrs_all, gene_id %in% novel_genes_all[!novel_genes_all %in% names(diff)])

goi = "Pou1f1"
goi_id <- names(genename[genename == goi])
ave_counts[goi_id]
ave_counts2[goi_id]
