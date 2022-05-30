#' ---
#' title: "Reprocessing of Ruf-Zamojski pituitary gland snRNA-seq data with data integration"
#' author: "Cadia Chan"
#' affiliation: "SickKids Research Institute & University of Toronto"
#' date: "May 17, 2021"
#' output:
#'  html_document:
#'    code_folding: hide
#'    toc: true
#'    toc_float: true
#' ---
#' This script reprocesses 10X single-nuclei RNA-seq data from Ruf-Zamojski et al. 2021 (GSE151961) using Seurat v4.
#' Only male and female snap-frozen tissues were reprocessed (this data is the main dataset used in the paper).
#' Seurat PBMC analysis workflow was followed except SCTransform was used instead to regress out
#' mitochondrial (mt-), and ribosomal (Rpl, Rps) genes.
#' Data integration between male and female samples was performed after SCTransform.
#' Workflow followed is found here: https://satijalab.org/seurat/articles/integration_introduction.html#performing-integration-on-datasets-normalized-with-sctransform-1.
#' Cell-type markers are identified and clusters were manually annotated
#' using pituitary cell-type markers curated from the original publication (Ruf-Zamojski et al. 2021).
#' scMappR is also run with the integrated dataset count matrix for bulk sex-biased genes from all ages.
#' 
#' ## Updates
#' * 2021-12-14: Added in male-female integration for CIBERSORT and scMappR.
#' * 2022-05-17: Removed single-cell deconvolution + scMappR scripts.
#' 
#' To run: Set working directory "To Source File Location".
#'
#+ setup
knitr::opts_chunk$set(message = FALSE, warning = F)

#' ## Libraries and source scripts
#+ load_library, warning = F, message = F
library(Seurat)
library(ggplot2)
library(sctransform)
library(dplyr)
library(ggpubr)
library(scales)
library(colorspace)
library(ggrepel)
library(biomaRt)
library(pheatmap)
library(knitr)

#' ## Load in data
#' Load in h5 data files downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE151961.
#' Only the first 6 samples are used (snap-frozen tissues); the other 3 samples are dissociated then cyropreserved nuclei.
#+ load_data
filelist <- list.files("input")

run_load_data <- function(filename) {
  sample <- substr(filename, 20, 22)
  data <- Read10X_h5(paste0("input/", filename), use.names = TRUE, unique.features = TRUE)
  data <- CreateSeuratObject(counts = data, project = sample, min.cells = 3, min.features = 200)
  return(data)
}

datalist <- lapply(filelist, function(x) run_load_data(x))

#' ## Merge datasets
#+ merge_data
names(datalist) <- substr(filelist, 20, 22)
data_f_comb <- merge(datalist[[4]], y=c(datalist[5:6]),
                     add.cell.ids = names(datalist)[4:6],
                     project = "ruf_f")
data_m_comb <- merge(datalist[[1]], y=c(datalist[2:3]),
                     add.cell.ids = names(datalist)[1:3],
                     project = "ruf_m")
unique(sapply(X = strsplit(colnames(data_f_comb), split = "_"), FUN = "[", 1))
table(data_f_comb$orig.ident)
unique(sapply(X = strsplit(colnames(data_m_comb), split = "_"), FUN = "[", 1))
table(data_m_comb$orig.ident)


#' ## Data filtering
#+ data_filtering
data_f_comb[["percent.mt"]] <- PercentageFeatureSet(data_f_comb, pattern = "^mt-")
data_f_comb[["percent.rpl"]] <- PercentageFeatureSet(data_f_comb, pattern = "^Rpl")
data_f_comb[["percent.rps"]] <- PercentageFeatureSet(data_f_comb, pattern = "^Rps")
data_m_comb[["percent.mt"]] <- PercentageFeatureSet(data_m_comb, pattern = "^mt-")
data_m_comb[["percent.rpl"]] <- PercentageFeatureSet(data_m_comb, pattern = "^Rpl")
data_m_comb[["percent.rps"]] <- PercentageFeatureSet(data_m_comb, pattern = "^Rps")

VlnPlot(data_f_comb, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rpl", "percent.rps"), ncol = 3)
VlnPlot(data_m_comb, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rpl", "percent.rps"), ncol = 3)

plot1 <- FeatureScatter(data_f_comb, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(data_f_comb, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(data_f_comb, feature1 = "nCount_RNA", feature2 = "percent.rpl")
plot4 <- FeatureScatter(data_f_comb, feature1 = "nCount_RNA", feature2 = "percent.rps")
plot1 + plot2 + plot3 + plot4

plot1 <- FeatureScatter(data_m_comb, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(data_m_comb, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(data_m_comb, feature1 = "nCount_RNA", feature2 = "percent.rpl")
plot4 <- FeatureScatter(data_m_comb, feature1 = "nCount_RNA", feature2 = "percent.rps")
plot1 + plot2 + plot3 + plot4

data_f_comb <- subset(data_f_comb, subset = nFeature_RNA >= 500 & percent.mt <= 15 & percent.rps <= 3 & percent.rpl <= 3)
VlnPlot(data_f_comb, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rpl", "percent.rps"), ncol = 3)

data_m_comb <- subset(data_m_comb, subset = nFeature_RNA > 500 & percent.mt <= 15 & percent.rps <= 3 & percent.rpl <= 3)
VlnPlot(data_m_comb, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rpl", "percent.rps"), ncol = 3)

#' ## Apply scTransform (replaces data normalization, variable feature id'n, and data scaling step)
#+ sctransform
data_f_comb <- SCTransform(data_f_comb, vars.to.regress = c("percent.mt", "percent.rpl", "percent.rps"), verbose = FALSE)
data_m_comb <- SCTransform(data_m_comb, vars.to.regress = c("percent.mt", "percent.rpl", "percent.rps"), verbose = FALSE)


#' ## Integrate data
#' Integrate female and male data.
#' IntegrateData step was run on Zhanglab gpu (see below for details).
#' Note that FindIntegrationAnchors step takes a long time.
#+ data_integration
ruf_list <- list(data_f_comb, data_m_comb)
names(ruf_list) <- c("f", "m")

features <- SelectIntegrationFeatures(object.list = ruf_list, nfeatures = 3000)
ruf_list <- PrepSCTIntegration(object.list = ruf_list, anchor.features = features)

ruf_anchors <- FindIntegrationAnchors(object.list = ruf_list, normalization.method = "SCT",
                                         anchor.features = features)

# Ran on zhanglab gpu
# Exported ruf_anchors as .rds and then ran the line below.
# Saved resulting integrated data file as .rds
# Inputed as ruf_list_sct below.

# ruf_list_sct <- IntegrateData(anchorset = ruf_anchors, normalization.method = "SCT")
ruf_list_sct <- readRDS("output/20211214_ruf_anchors_output.rds")

ruf_list_sct <- RunPCA(ruf_list_sct, verbose = FALSE)
ElbowPlot(ruf_list_sct)
ruf_list_sct <- RunUMAP(ruf_list_sct, reduction = "pca", dims = 1:16)

ruf_list_sct$sex <- substr(ruf_list_sct$orig.ident, 1, 1)

#' ## Find clusters
#+ find_clusters
ruf_list_sct <- FindNeighbors(ruf_list_sct, reduction = "pca", dims = 1:16)
ruf_list_sct <- FindClusters(ruf_list_sct, resolution = 0.4)

#' ## Make UMAP
#+ run_umap
p1 <- DimPlot(ruf_list_sct, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(ruf_list_sct, reduction = "umap", group.by = "sex")
p3 <- DimPlot(ruf_list_sct, reduction = "umap", group.by = "seurat_clusters", label = T)
p1 + p2
ggsave("output/20211214_integrated_umap_samples.pdf", width = 12, height = 8)
p3
ggsave("output/20211214_integrated_umap_clusters.pdf", width = 9, height = 8)


#' ## Find markers
#+ find_markers
int_markers <- FindAllMarkers(ruf_list_sct, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

#' ## Reassign cluster names
#+ cluster_name
int_celltypes <- 
  c(`0` = "Lactotropes",
    `1` = "Somatotropes_1",
    `2` = "Somatotropes_2",
    `3` = "Melanotropes",
    `4` = "Somato/Lacto",
    `5` = "Gonadotropes",
    `6` = "Stem-cell_1",
    `7` = "Stem-cell_2",
    `8` = "Corticotropes",
    `9` = "Proliferating",
    `10` = "Thyrotropes", 
    `11` = "Pericytes",
    `12` = "Endothelial", #Meis2/Plvap
    `13` = "Pituicytes (posterior)",
    `14` = "Debris", #
    `15` = "Macrophages"
  )
Idents(ruf_list_sct) <- ruf_list_sct$seurat_clusters
ruf_list_sct <- RenameIdents(ruf_list_sct, int_celltypes)
DimPlot(ruf_list_sct, reduction = "umap", label = T)
ggsave("output/20211214_integrated_umap_clusters_labelled.pdf", width = 12, height = 10)


#' ## Save objects and markers
#+ save_output
saveRDS(ruf_list_sct, "output/20211214_ruf_2021_pit_GSE151961_integrated_seurat.rds")

write.table(int_markers, "output/20211214_ruf_2021_pit_GSE151961_integrated_markers.txt",
            sep = "\t", col.names = T, row.names = F, quote = F)


#' ## Make signature matrix for scMappR
#' Gene signature matrix is prepared for scMappR using both positive and negative gene markers identified by Seurat FindAllMarkers function.
#' Signature matrix contains odds ratio calculated by scMappR generes_to_heatmap function.
#+ make_scmappr_OR, warning = F, message = F
# Make scmappr input (OR matrix) ####
dir.create("output/scmappr")
ruf_list_sct <- readRDS("output/20211214_ruf_2021_pit_GSE151961_integrated_seurat.rds")

int_simple_celltypes <- 
  c(`0` = "Lactotropes",
    `1` = "Somatotropes",
    `2` = "Somatotropes",
    `3` = "Melanotropes",
    `4` = "Somato_Lacto",
    `5` = "Gonadotropes",
    `6` = "Stem_cell",
    `7` = "Stem_cell",
    `8` = "Corticotropes",
    `9` = "Proliferating",
    `10` = "Thyrotropes", 
    `11` = "Pericytes",
    `12` = "Endothelial", #Meis2/Plvap
    `13` = "Pituicytes_posterior",
    `14` = "Debris", #
    `15` = "Macrophages"
  )
ruf_list_sct$man_clusters <- Idents(ruf_list_sct)
Idents(ruf_list_sct) <- ruf_list_sct$seurat_clusters
ruf_list_sct <- RenameIdents(ruf_list_sct, int_simple_celltypes)
DimPlot(ruf_list_sct, reduction = "umap", label = T)
ggsave("output/20211214_integrated_umap_clusters_simple_labelled.pdf", width = 12, height = 10)

#' ## R Session
#+ r_session
sessionInfo()
