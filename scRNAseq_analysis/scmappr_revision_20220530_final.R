#' ## Libraries and source scripts
#+ load_library, warning = F, message = F
#### Libraries ####
library(dplyr)
library(Seurat)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(scMappR)
library(scales)
library(colorspace)
library(ggrepel)
library(openxlsx)
library(biomaRt)
library(pheatmap)
library(knitr)
library(EDASeq)
library(WGCNA)
library(RColorBrewer)

source("scmappr_revision_functions.R")

#' ## Load bulk data
#+ load_bulk, warning = F, message = F
# Load bulk ####
bulk_obj <- readRDS("input/pit_utr_2019__RUV_k2_set1_2019-07-03.rds")
bulk_counts <- normCounts(bulk_obj)
bulk_de <- readRDS("input/pit_utr_2019_de_result_list_2019-07-03.rds")
bulk_de <- lapply(bulk_de, function(x) x[abs(x$logFC) > log2(1.5) & x$FDR < 0.05, ])
meta <- pData(bulk_obj)

#### Ruf-Zamojski data ####
#' ## Run scMappR with Ruf-Zamojski 2021 adult male and female data
#' Gene signature matrix is prepared for scMappR using both positive and negative gene markers identified by Seurat FindAllMarkers function.
#' Signature matrix contains odds ratio calculated by scMappR generes_to_heatmap function.
#+ make_scmappr_sci, warning = F, message = F
# Make scmappr input (OR matrix) ####
dir.create("output/ruf_scmappr")

ruf <- readRDS("input/20211214_ruf_2021_pit_GSE151961_integrated_seurat.rds")

int_simple_celltypes <- 
  c(Lactotropes = "Lactotropes",
    Somatotropes_1 = "Somatotropes",
    Somatotropes_2 = "Somatotropes",
    Melanotropes = "Melanotropes",
    `Somato/Lacto` = "Somato_Lacto",
    Gonadotropes = "Gonadotropes",
    `Stem-cell_1` = "Stem_cell",
    `Stem-cell_2` = "Stem_cell",
    Corticotropes = "Corticotropes",
    Proliferating = "Proliferating",
    Thryotropes = "Thyrotropes", 
    Pericytes = "Pericytes",
    Endothelial = "Endothelial", #Meis2/Plvap
    `Pituicytes (posterior)` = "Pituicytes",
    Debris = "Debris", #
    Macrophages = "Macrophages"
  )
ruf$man_clusters <- Idents(ruf)
ruf <- RenameIdents(ruf, int_simple_celltypes)

ruf$int_clusters <- factor(ruf@active.ident,
                           levels = list(Somatotropes = "Somatotropes", 
                                         Lactotropes = "Lactotropes", 
                                         Corticotropes = "Corticotropes", 
                                         Gonadotropes = "Gonadotropes", 
                                         Stem_cell = "Stem_cell",
                                         Endothelial = "Endothelial", 
                                         Proliferating = "Proliferating", 
                                         Melanotropes = "Melanotropes", 
                                         Pericytes = "Pericytes", 
                                         Macrophages = "Macrophages", 
                                         Pituicytes = "Pituicytes",
                                         Thyrotropes = "Thyrotropes",
                                         `Somato/Lacto` = "Somato_Lacto",
                                         Debris = "Debris"))
Idents(ruf) <- ruf$int_clusters
ruf <- subset(ruf, idents = "Debris", invert = T) # remove debris cluster
# Note the return threshold for the following command is 0.1 whereas the FindMarkers command used in scMappR seurat_to_generes
# command has no threshold set, so that all genes are returned regardless of pvalue. This means that
# the command used here is more stringent.
ruf_markers_sct <- FindAllMarkers(ruf, assay = "SCT", only.pos = F, recorrect_umi = F)
generes <- lapply(levels(ruf_markers_sct$cluster), function(x) ruf_markers_sct[ruf_markers_sct$cluster == x, ])
names(generes) <- levels(ruf_markers_sct$cluster)
generes <- lapply(generes, function(x) {rownames(x) <- as.character(x$gene); x})
generes <- lapply(generes, function(x) x[, c(5, 2)])

wilcoxon_scmappr <- generes_to_heatmap(generes, species = "mouse", make_names = F)
kable(head(wilcoxon_scmappr$OR), caption = "Example gene OR in signature matrix.")

HM <- wilcoxon_scmappr$OR
signatureVar <- rowVars(HM)
signature_o <- HM[order(signatureVar, decreasing = TRUE),]
if(nrow(signature_o) > 3000) {
  signature_o <- signature_o[1:3000,]
}

DeconMethod <- compare_deconvolution_methods(as.data.frame(bulk_counts), as.data.frame(signature_o))
DeconRNAseq <- DeconMethod$cellType_proportions$DeconRNAseq
WGCNA <- DeconMethod$cellType_proportions$WGCNA
DCQ <- DeconMethod$cellType_proportions$DCQ
rownames(DeconRNAseq) <- colnames(bulk_counts)

#' ## Run scMappR with custom sig mat from sci-seq
#+ run_sci_scmappr
rufdata <- ruf
rm(ruf)
scrna <- "ruf"
ruf_genes <- rownames(rufdata)
# setwd("output_files/")
# knitr::opts_knit$set(root.dir = paste0('output/", scrna, "_scmappr/'))

contrast_df <- matrix(c("Pd12F", "Pd12M", "Pd22F", "Pd22M", "Pd27F", "Pd27M", "Pd32F", "Pd32M", "Pd37F", "Pd37M",
                        "Pd22M", "Pd12M", "Pd27M", "Pd22M", "Pd32M", "Pd27M", "Pd37M", "Pd32M",
                        "Pd22F", "Pd12F", "Pd27F", "Pd22F", "Pd32F", "Pd27F", "Pd37F", "Pd32F",
                        "Pd37M", "Pd12M", "Pd37F", "Pd12F", "Pd37M", "Pd22M", "Pd37F", "Pd22F"),
                      nrow = 17, ncol = 2, byrow = T, dimnames = list(names(bulk_de), c("case", "ctrl")))
contrast_df <- contrast_df[-grep("vs", rownames(contrast_df)), ]

run_scmappr("d37_sex", contrast_df["d37_sex", 1], contrast_df["d37_sex", 2], file_out = paste0("./output/", scrna, "_scmappr/"))
run_scmappr("d32_sex", contrast_df["d32_sex", 1], contrast_df["d32_sex", 2], file_out = paste0("./output/", scrna, "_scmappr/"))
run_scmappr("d27_sex", contrast_df["d27_sex", 1], contrast_df["d27_sex", 2], file_out = paste0("./output/", scrna, "_scmappr/"))
run_scmappr("d22_sex", contrast_df["d22_sex", 1], contrast_df["d22_sex", 2], file_out = paste0("./output/", scrna, "_scmappr/"))
run_scmappr("d12_sex", contrast_df["d12_sex", 1], contrast_df["d12_sex", 2], file_out = paste0("./output/", scrna, "_scmappr/"))


#' ## Ruf-Zamojski cwFC evaluation
#+ cwfc_eval_ruf
ruf_markers_sct_pos <- filter(ruf_markers_sct, avg_log2FC > 0.25, p_val < 0.01)
lapply(c(12, 22, 27, 32, 37), function(x) plot_cwfc(x, scrna, ruf_genes, rufdata, ruf_markers_sct_pos))

