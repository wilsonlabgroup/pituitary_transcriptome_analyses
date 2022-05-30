#' ---
#' title: "Reprocessing of pituitary gland scRNA-seq data"
#' author: "Cadia Chan"
#' affiliation: "SickKids Research Institute & University of Toronto"
#' date: "March 2, 2022"
#' output:
#'  html_document:
#'    code_folding: hide
#'    toc: true
#'    toc_float: true
#' ---
#' This script reprocesses 10X single-cell RNA-seq data from Cheung et al. 2018 (GSE120410) using Seurat v4.
#' Standard filtering and normalization of data is performed (following Seurat PBMC analysis workflow).
#' Cell-type markers are identified and clusters are manually annotated using pituitary cell-type markers curated from the original publication (Cheung et al. 2018)
#' as well as rat anterior pituitary single-cell study (Fletcher et al. 2019), and mouse posterior pituitary single-cell study (Chen et al. 2019).
#' Using the reprocessed scRNA-seq data, enrichment for genes in co-expression modules within different single-cell clusters is assessed.
#' Enrichment is calculated based on expression of a given gene within a cell cluster compared to all other cell clusters using a one-tailed Kolmogorov-Smirnov (KS) test.
#' A gene is considered enriched if FDR < 0.05 and FC > 0 in a cell cluster compared to other cell clusters.
#' 
#' ## Updates
#' * 2022/03/02: Changed FindClusters resolution to 0.2. This decreased subclustering in corticotropes, and lactotropes. However, somatotropes still subcluster into 4.;
#' Updated scaling factor for filling in mean counts for signature matrix: changed from times 1000 to times numcells in the single cell dataset
#' * 2022/03/03: Changed to using Seurat v 4.0.5 instead of 3.2.3 (So that it becomes compatible with sci-seq data which relies on v4)
#' To run: Set working directory "To Source File Location".
#'

#' ## Libraries and source scripts
#+ load_library, warning = F, message = F
#### Libraries ####
library(dplyr)
library(Seurat)
library(knitr)
library(biomaRt)
library(scales)
library(EDASeq)
library(pheatmap)
library(LaCroixColoR)
library(grid)
library(plyr)
library(ggpubr)
library(reshape2)
library(ggplot2)
library(colorspace)
library(ggrepel)
library(openxlsx)

source("cheung_sc_reprocess_2021-03-22_functions_rmd_ver.R")

#' ## Load in GSE120410 Cheung et al 2018 single-cell data
#+ load_data, warning = F, message = F
#### Load in scRNA-seq data ####
pit_data <- Read10X(data.dir = "input_files/GSE120410_Cheung_2018/")
# pit_data <- Read10X(data.dir = "~/Dropbox (Wilson Lab)/Mike_Anna_Mark/active_manuscript/Pituitary paper/scripts/single_cell_rmd/input_files/GSE120410_Cheung_2018/")
pit_seurat <- CreateSeuratObject(counts = pit_data,
                                 project = "cheung_pit",
                                 min.cells = 3,
                                 min.features = 200)
pit_seurat

#' ## Unwanted cell filtering
#' In each cell, number of counts, number of RNA features (unique genes), and mitochondrial content is assessed.
#' Based on the scatter plots, cells are filtered for: nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt < 20.
#' This results in 13565/13620 cells retained.
#+ data_prefilt, warning = F, message = F, fig.cap = "Visualization of QC metrics pre-filtering", fig.width = 10, fig.height = 8
#### Cell filtering ####
# Visualize QC metrics as a violin plot and use to filter cells by unique feature counts
pit_seurat[["percent.mt"]] <- PercentageFeatureSet(pit_seurat, pattern = "^mt-")
print(VlnPlot(pit_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(pit_seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pit_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
print(plot1+plot2)

#+ data_postfilt, warning = F, message = F, fig.cap = "Visualization of QC metrics after filtering", fig.width = 10, fig.height = 8
# Apply filter to data based on QC metrics
pit_seurat <- subset(pit_seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt < 20)
pit_seurat # After subsetting with filter
print(VlnPlot(pit_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
plot1 <- FeatureScatter(pit_seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pit_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
print(plot1+plot2)

#' ## Global data normalization
#' Apply global-scaling normalization of the feature expression measurements for each cell by the total expression.
#' Counts are then mulitplied by a scaling factor (10,000) and then log-transformed.
#+ data_norm, warning = F, message = F
#### Data global normalization ####
pit_seurat <- NormalizeData(pit_seurat, normalization.method = "LogNormalize", scale.factor = 10000)

#' ## Gene feature identification
#' Genes (referred to as features) exhibiting high cell-to-cell variation are identified.
#' These cells are likely highly expressed in some cells but lowly expressed in others.
#' Focusing on these genes downstream helps highlight important biological signal in the dataset.
#' Selected features will be used in downstream analyses like PCA.
#' Here, we used 3000 variable features based on above scatter plots.
#+ feature_sel, warning = F, message = F, fig.width = 10, fig.height = 8
#### Feature (gene) selection ####
pit_seurat <- FindVariableFeatures(pit_seurat, selection.method = "vst", nfeatures = 3000) # based on violin plot

# Identify the 20 most highly variable genes
top20 <- head(VariableFeatures(pit_seurat), 20)
top20

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pit_seurat)
plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
print(plot1 + plot2)

#' ## Data scaling
#' Linear transformation 'scaling' prior to dimensional reduction like PCA shifts expression of each gene so mean expression across cells is 0.
#' Scales gene expression for each gene so that variance across cells is 1 (so that highly expressed genes do not dominate).
#+ data_scaling, warning = F, message = F
#### Data scaling ####
all.genes <- rownames(pit_seurat)
pit_seurat <- ScaleData(pit_seurat, features = all.genes) # pit_seurat[["RNA"]]@scale.data

#' ## Linear dimensionality reduction
#' 3000 variable features chosen previously are used for linear dimensionality reduction with PCA.
#' Genes associated with the top 5 PCs are visualized.
#' Expression of genes associated with the top 15 PCs in 500 cells displaying "extreme" expression from both ends of the spectrum are shown in the heatmap.
#' Heatmap can be used to visualize sources heterogeneity in the dataset and inform the which PCs to include in downstream analyses.
#+ dim_red_pca, warning = F, message = F, fig.cap = "Genes associated with top 5 principal components.", fig.width = 10, fig.height = 10
#### Linear dimension reduction (PCA) ####
pit_seurat <- RunPCA(pit_seurat, features = VariableFeatures(object = pit_seurat))
print(pit_seurat[["pca"]], dims = 1:5, nfeatures = 5)
print(VizDimLoadings(pit_seurat, dims = 1:5, reduction = "pca"))

#+ pca_dim_plot, warning = F, message = F, fig.cap = "Top 2 principal components of data.", fig.width = 10, fig.height = 10
DimPlot(pit_seurat, reduction = "pca")

#+ pca_dim_htmap, warning = F, message = F, fig.cap = "Expression of genes associated with the first 15 principal components.", fig.width = 10, fig.height = 10
DimHeatmap(pit_seurat, dims = 1:15, cells = 500, balanced = TRUE)

#' ## Determining "dimensionality" of dataset
#' Cells are clustered based on their PCA scores where each PC is a "metafeature" that combines info across correlated feature set.
#' To determine the number of PCs to use, an Elbow Plot is generated.
#' The elbow plot ranks the principle components based on % variance explained by each one.
#' The point at which we observe an 'elbow' or decrease in steepness of the line suggests that the majority of true signal is captured before this point.
#' Here we observe the elbow around PC15, hence we will use 15 PCs for our downstream analysis.
#+ elbow_plot, warning = F, message = F, fig.cap = "Elbow plot ranking principal components."
#### Elbow plot ####
print(ElbowPlot(pit_seurat))

#' ## Cell clustering
#' Graph-based clustering approach where genes which similar feature expression patterns have edges drawn between them.
#' Partition of graph into "communities" where nodes/genes are highly interconnected is then applied after.
#' First, K-nearest neighbour (KNN) graph based on euclidean distance in PCA space is constructed.
#' Second, edge weights between cells are refined based on shared overlap in 'local neighbourhood' (Jaccard similarity).
#' This second step takes in as input the previously defined dimensionality of the dataset (first 15 PCs).
#+ cell_clustering, warning = F, message = F
#### Cell clustering ####
pit_seurat <- FindNeighbors(pit_seurat, dims = 1:15)

#' To cluster the cells, Louvain algorithm is applied.
#' Cells are grouped together iteratively, with goal of optimizing standard modularity function.
#' Resolution increases for larger datasets and increases the number of clusters resulting.
#' Here, 14 cell clusters result.
#' From the previous version, the resolution was changed to 0.2. 
#+ find_clusters, warning = F, message = F
#### Louvain algorithm ####
pit_seurat <- FindClusters(pit_seurat, resolution = 0.2)
head(Idents(pit_seurat), 5)

#' ## Non-linear dimensionality reduction
#' UMAP reduce dimensions of the data and allow for visualization of the dataset in 2D space.
#' These algorithms try to learn the underlying manifold of the data to try to place similar cells together in low-dim space.
#+ umap_dim_red, warning = F, message = F, fig.show="hold", out.width="50%", fig.cap = "UMAP visualization of cell clusters.", fig.width = 28, fig.height = 14
#### Non-linear dimensional reduction ####
pit_seurat <- RunUMAP(pit_seurat, dims = 1:15)
dir.create("output_files")
p1 <- DimPlot(pit_seurat, reduction = "umap", label = T)
p2 <- DimPlot(pit_seurat, reduction = "umap", label = F)
pdf("output_files/seurat_v4_cheungpit_umap.pdf", width = 7, height = 6)
print(p1)
print(p2)
tmp <- dev.off()
print(p1+p2)

#' ## Finding cluster markers
#' Differential expression of genes within different clusters is tested to identify gene markers for each cluster.
#' Positive gene markers are identified for each cluster. 
#' Only genes which are detected in >25% in either one of the two clusters of a given comparison are used.
#' The default Wilcoxon test is used with logFC threshold of 0.25.
#' Violin plots and UMAP are used to visualize known anterior hormone-producing cell-type markers in each cluster.
#' Other pituitary cell-type markers are also visualized (pituicytes, FSCs, endothelial cells).
#' Finally, the top 10 expressing genes in each cluster is plotted as a heatmap.
#+ find_markers, warning = F, message = F
#### Find cluster biomarkers ####
# pit_seurat <- readRDS("output_files/seurat_v3_pitCheung.rds")
Idents(pit_seurat) <- pit_seurat$seurat_clusters

pit_markers <- FindAllMarkers(pit_seurat, only.pos = T, min.pct = 0.25)
kable(head(pit_markers), caption = "Example cluster gene markers identified.")
write.table(pit_markers, "output_files/seurat_v4_cheungpit_allposmarker_genes.txt",
            sep = "\t", col.names = T, row.names = F, quote = F)

dir.create("output_files/seurat_marker_plots")

#+ plot_anterior_markers_violin, message = F, warning = F, fig.cap = "Violin plot of anterior pituitary gene marker expression.", fig.width = 13, fig.height = 8
# Gh = somatotrope (0,2)
# Pou1f1 = Pou1f1 expressing (1)
# Prl = lactotrope (1)
# Lhb, Fshb = gonadotrope (3)
# Pomc, Crhr1 = corticotrope (4)
# Tshb = thyrotrope (5)
# Sox2, Rbpms = stem-cells (6)
# Top2a, Mki67 = Proliferating (7)
# Emcn, Meis2, Plvap = Endothelial (8)
# Pomc, Pax7, Oacyl = melanotrope (9) intermediate lobe
# Dcn, Lama2 = Pericytes (10)
# C1qa, Ctss = Macrophages (11)
# Scn7a, Col25a1 = Pituicytes (12)

p <- VlnPlot(pit_seurat, features = c("Gh", "Prl", "Lhb", "Fshb", "Tshb", "Pomc", "Crhr1", "Pax7", "Oacyl", "Sox2", "Rbpms", "Pou1f1",
                                      "Top2a", "Mki67", "Plvap", "Meis2", "Emcn", "Dcn", "Lama2", "C1qa", "Ctss",
                                      "Scn7a", "Col25a1"), pt.size = 0.5)
print(p)
pdf("output_files/seurat_marker_plots/cheungpit_markers_violin.pdf", width = 20, height = 20)
print(p)
tmp <- dev.off()
ggsave("output_files/seurat_marker_plots/cheungpit_markers_violin.png", p, width = 20, height = 25)

#+ plot_markers_umap, message = F, warning = F, fig.cap = "UMAP of pituitary gene markers.", fig.width = 13, fig.height = 8
p <- FeaturePlot(pit_seurat, features = c("Gh", "Prl", "Lhb", "Fshb", "Tshb", "Pomc", "Crhr1", "Pax7", "Oacyl", "Sox2", "Rbpms", "Pou1f1",
                                          "Top2a", "Mki67", "Plvap", "Meis2", "Emcn", "Dcn", "Lama2", "C1qa", "Ctss",
                                          "Scn7a", "Col25a1"), reduction = "umap")
print(p)
pdf("output_files/seurat_marker_plots/cheungpit_markers_umap.pdf", width = 20, height = 25)
print(p)
tmp <- dev.off()
ggsave("output_files/seurat_marker_plots/cheungpit_markers_umap.png", p, width = 20, height = 25)

#+ plot_top10_heatmap, message = F, warning = F, fig.cap = "Heatmap of top 10 gene markers in each cell cluster.", fig.width = 12, fig.height = 12
top10 <- pit_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
p <- DoHeatmap(pit_seurat, features = top10$gene) + NoLegend()
print(p)
pdf("output_files/seurat_marker_plots/cheungpit_TOP10markers_heatmap.pdf", width = 12, height = 12)
print(p)
tmp <- dev.off()
ggsave("output_files/seurat_marker_plots/cheungpit_TOP10markers_heatmap.png", p, width = 20, height = 25)

#+ plot_thyro_heatmap, message = F, warning = F, fig.cap = "UMAP of thyrotropes gene markers based on list in Cheung et al 2018.", fig.width = 10, fig.height = 8
p <- FeaturePlot(pit_seurat, features = c("Tshb", "Trhr", "Dio2"), reduction = "umap")
print(p)
pdf("output_files/seurat_marker_plots/cheungpit_thryotropes_markers_umap.pdf", width = 10, height = 8)
print(p)
tmp <- dev.off()
ggsave("output_files/seurat_marker_plots/cheungpit_thryotropes_markers_umap.pdf", p, width = 10, height = 8)


#' ## Manual cluster labelling
#' Clusters on are manually labelled and grouped based on curated known marker genes. 
#+ manual_label, message = F, warning = F
#### Manual cluster labelling ####

new_cluster_ids <- c(`0`="Somatotropes_0", #0
                     `1`="Lactotropes",#1
                     `2`="Somatotropes_1",#2
                     `3`="Gonadotropes",#3
                     `4`="Corticotropes",#4
                     `5`= "Thyrotropes",#5
                     `6`= "Stem_cells_Sox2", #6
                     `7`="Proliferating", #7
                     `8`="Endothelial", #8
                     `9`="Melanotropes", #9
                     `10`="Pericytes", #10
                     `11`="Macrophages", #11
                     `12`="Pituicytes") #12
pit_seurat <- RenameIdents(pit_seurat, new_cluster_ids)

#' ## Visualization of manually labelled cell clusters
#' UMAPs are revisualized with manually labelled clusters.
#+ umap_manual, message = F, warning = F, fig.show="hold", out.width="50%", fig.cap = "UMAP visualization of manually labelled cell clusters.", fig.width = 28, fig.height = 14
#### Manual cell cluster UMAP ####
p1 <- DimPlot(pit_seurat, reduction = "umap", label = T)
p2 <- DimPlot(pit_seurat, reduction = "umap", label = F)
print(p1+p2)
pdf("output_files/cheungpit_umap_manual_celltypes.pdf", width = 10, height = 10, useDingbats = F)
print(p1)
print(p2)
tmp <- dev.off()


#' ## Merge manual cluster labelling
#' Manually labelled clusters are merged based on cell type. 
#+ merge_manual_label, message = F, warning = F
#### Merge manual cluster ####
pit_seurat$manual_anno <- Idents(pit_seurat)
merge_cluster_ids <- c(Somatotropes_0 = "Somatotropes", #0
                     Lactotropes = "Lactotropes", #1
                     Somatotropes_1 = "Somatotropes",#2
                     Gonadotropes = "Gonadotropes",#3
                     Corticotropes = "Corticotropes",#4
                     Thyrotropes = "Thyrotropes", #5
                     Stem_cells_Sox2 = "Stem_cells_Sox2",#6
                     Proliferating = "Proliferating", #7
                     Endothelial = "Endothelial", #8
                     Melanotropes = "Melanotropes", #9
                     Pericytes = "Pericytes", #10
                     Macrophages = "Macrophages", #11
                     Pituicytes = "Pituicytes") #12
pit_seurat <- RenameIdents(pit_seurat, merge_cluster_ids)

#' ## Visualization of merged manually labelled cell clusters
#' UMAPs are revisualized with merged manually labelled clusters.
#+ merge_umap_manual, message = F, warning = F, fig.show="hold", out.width="50%", fig.cap = "UMAP visualization of merged manually labelled cell clusters.", fig.width = 28, fig.height = 14
#### Merged manual cell cluster UMAP ####
p1 <- DimPlot(pit_seurat, reduction = "umap", label = T)
p2 <- DimPlot(pit_seurat, reduction = "umap", label = F)
print(p1+p2)
pdf("output_files/cheungpit_umap_merged_manual_celltypes.pdf", width = 10, height = 10, useDingbats = F)
print(p1)
print(p2)
tmp <- dev.off()

#' ## Manually labelled cluster gene marker identification
#' Gene markers are reidentified using manually grouped cell clusters.
#' The same parameters are run with FindAllMarkers.
#' Additionally, both positive and negative gene markers are identified for downstream analyses purposes (scMappR).
#+ manual_find_marker, message = F, warning = F
#### Manual marker identification ####
pos_markers <- FindAllMarkers(pit_seurat, only.pos = T, min.pct = 0.25)
kable(head(pos_markers), caption = "Example manual cluster gene markers identified.")
write.table(pos_markers, "output_files/seurat_v4_cheung_etal_manual_anno_cell_cluster_allposmarkers.txt", sep = "\t", quote = F,col.names = T, row.names = F)

man_markers <- FindAllMarkers(pit_seurat, only.pos = F, min.pct = 0.25)
write.table(man_markers, "output_files/seurat_v4_cheung_etal_manual_anno_cell_cluster_allmarkers.txt", sep = "\t", quote = F,col.names = T, row.names = F)
# saveRDS(man_markers, "seurat_v4_cheung_all_markers.rds")

#' ## Save manually labelled Seurat object
#' Overwrites previous object
#+ save_manual_seurat, message = F, warning = F
saveRDS(pit_seurat, "output_files/seurat_v4_pitCheung.rds")

#' # Enrichment of module genes in single cell populations
#' ## Setup data for enrichment test
#' Co-expression module gene data is loaded in.
#' A colour palette is made based on the cell-types.
#+ setup_modules_data, warning = F, message = F
#### Enrichment of module genes in single cell populations ####
seurat <- pit_seurat
# seuratman <- pit_seurat_manual
clusters <- seurat$RNA_snn_res.0.2
manclusters <- Idents(seurat)
mod_utr <- readRDS("input_files/pit_utrseq_coexpression_modules_annotation.rds")
mod_utr$modules <- factor(mod_utr$modules)

celltypes <- merge_cluster_ids

# Create colour palette for heatmap annotation of cell types
ggcolours <- hue_pal()(12)
names(ggcolours) <- levels(manclusters)
celltypes_col <- celltypes
for(i in 1:length(celltypes_col)) {
  celltypes_col[i] <- ggcolours[celltypes_col[i]]
}
names(celltypes_col) <- as.character(celltypes)

#' ## Run KS test for all cell clusters
#' KS test is run and FDR is calculated for each module gene in each cell cluster (not manually labelled).
#' Genes with enriched with FC > 0, FDR adj < 0.05 in at least one cell cluster compared to all other clusters are retained.
#' Heatmap is plotted with FDR value for each module gene retained in each cell cluster.
#+ run_ks_test, warning = F, message = F
#Run KS test and calculate FDR for each module gene
ks_res <- lapply(levels(mod_utr$modules), function(x) calc_stat_test(mod_utr[mod_utr$modules == x, "genename"], "ks", seurat, clusters))
ks_respval <- lapply(ks_res, function(x) na.exclude(x[["pval"]]))
num_comparisons <- sum(unlist(lapply(ks_respval, function(x) print(as.numeric(dim(x)[1]) * as.numeric(dim(x)[2])))))

#Calculate FDR adjusted pvalue and filter for minimum FDR adj < 0.05 in each row
ks_resfdr <- lapply(ks_respval, function(x) matrix(p.adjust(as.matrix(x), method = "fdr", n = sum(num_comparisons)), nrow = nrow(x), ncol = ncol(x), dimnames = dimnames(x)))
ks_fdrcut <- lapply(ks_resfdr, function(x) x[rowMin(x) <= 0.05, ])


#Filter FC calculated so that only genes with FC > 0 in at least one comparison is included
ks_fc <- lapply(ks_res, function(x) na.exclude(x[["FC"]]))
names(ks_fdrcut) <- names(ks_fc) <- paste0("M", 1:9)

ks_fccut <- lapply(names(ks_fc), function(x) ks_fc[[x]][rownames(ks_fc[[x]]) %in% rownames(ks_fdrcut[[x]]), ])
ks_fccut <- lapply(ks_fccut, function(x) x[abs(rowSums(x)) > 0, ])

#Combine FDR df from all modules into one df
all_ks <- bind_rows(lapply(ks_fdrcut, function(x) as.data.frame(cbind("gene" = rownames(x), x))))
rownames(all_ks) <- all_ks$gene
all_ks <- all_ks[, -1]
gene_names <- rownames(all_ks)
all_ks <- sapply(all_ks, function(x) as.numeric(as.character(x)))
rownames(all_ks) <- gene_names

#Annotate each gene to corresponding module
anno_row <- bind_rows(lapply(names(ks_fdrcut), function(x) as.data.frame(cbind("gene" = rownames(ks_fdrcut[[x]]), "module" = x) )))
rownames(anno_row) <- anno_row$gene
anno_col <- list(celltype = ggcolours)

#Calculate breaks in heatmap based on modules
indices <- sapply(ks_fdrcut, nrow)
for(i in 2:length(indices)) {
  x <- as.numeric(indices[i])
  indices[i] <- as.numeric(indices[i - 1]) + as.numeric(indices[i])
}

cell_anno <- as.data.frame(cbind("cluster" = 0:12, "celltype" = celltypes))
rownames(cell_anno) <- cell_anno$cluster

use_breaks <- get_htmap_breaks(all_ks)

#+ ks_heatmap_all_clusters, warning = F, message = F, width = 8, height = 10, fig.cap = "Cell clusters are clustered based on similarity."
p <- pheatmap(all_ks,
              # border_color = "black",
              color = colorRampPalette(c("darkmagenta", "snow2"))(250),
              annotation_col = cell_anno[, 2, drop = F],
              cluster_rows = F,
              cluster_cols = T,
              # scale = "row",
              annotation_colors = anno_col,
              gaps_row = indices,
              show_rownames = F,
              labels_col = paste(cell_anno$cluster, cell_anno$celltype) 
              # breaks = use_breaks2,
)
save_phtmap_pdf(p, "output_files/pit_utr_module_cheung_scrna_ks_FDR_cellenrichment.pdf", 8, 10)
print(p)
tmp <- dev.off()

#+ ks_heatmap_all_clusters_ordered, warning = F, message = F, width = 8, height = 10, fig.cap = "Cell clusters are clustered based on manual cell-type annotation."
all_ks_fdr <- bind_rows(lapply(ks_fdrcut, function(x) as.data.frame(cbind("gene" = rownames(x), x))))
rownames(all_ks_fdr) <- all_ks_fdr$gene
all_ks_fdr <- all_ks_fdr[, -1]
gene_names <- rownames(all_ks_fdr)
all_ks_fdr <- sapply(all_ks_fdr, function(x) as.numeric(as.character(x)))
rownames(all_ks_fdr) <- gene_names

indices <- sapply(ks_fdrcut, nrow)
for(i in 2:length(indices)) {
  x <- as.numeric(indices[i])
  indices[i] <- as.numeric(indices[i - 1]) + as.numeric(indices[i])
}

cell_anno <- as.data.frame(cbind("cluster" = 0:12, "celltype" = celltypes))
rownames(cell_anno) <- cell_anno$cluster
all_ks_fdr_order <- all_ks_fdr[, order(cell_anno$celltype)]
p <- pheatmap(all_ks_fdr_order,
              # border_color = "black",
              color = colorRampPalette(c("darkmagenta", "snow3"))(250),
              annotation_col = cell_anno[, 2, drop = F],
              cluster_rows = F,
              cluster_cols = F,
              # scale = "row",
              annotation_colors = anno_col,
              gaps_row = indices,
              show_rownames = F
              # breaks = use_breaks,
)
save_phtmap_pdf(p, "output_files/pit_utr_module_cheung_scrna_ks_FDR_cellenrichment_cellordered.pdf", 8, 10)
print(p)
tmp <- dev.off()

#' ## Run KS test for manually annotated cell-types
#' KS test is run and FDR is calculated for each module gene in each manually annotated cell-type.
#' Genes with enriched with FC > 0, FDR adj < 0.05 in at least one cell-type compared to all other cell-types are retained.
#' Heatmap is plotted with FDR value for each module gene retained in each cell-type.
#+ run_ks_test_manual, warning = F, message = F
#### KS test manual annotation ####
#Run KS test and calculate FDR for each module gene
ks_res <- lapply(levels(mod_utr$modules), function(x) calc_stat_test(mod_utr[mod_utr$modules == x, "genename"], "ks", seurat, manclusters))
ks_respval <- lapply(ks_res, function(x) na.exclude(x[["pval"]]))
num_comparisons <- sum(unlist(lapply(ks_respval, function(x) print(as.numeric(dim(x)[1]) * as.numeric(dim(x)[2])))))

#Calculate FDR adjusted pvalue and filter for minimum FDR adj < 0.05 in each row
ks_resfdr <- lapply(ks_respval, function(x) matrix(p.adjust(as.matrix(x), method = "fdr", n = sum(num_comparisons)), nrow = nrow(x), ncol = ncol(x), dimnames = dimnames(x)))
ks_fdrcut <- lapply(ks_resfdr, function(x) x[rowMin(x) <= 0.05, ])


#Filter FC calculated so that only genes with FC > 0 in at least one comparison is included
ks_fc <- lapply(ks_res, function(x) na.exclude(x[["FC"]]))
names(ks_fdrcut) <- names(ks_fc) <- paste0("M", 1:9)

ks_fccut <- lapply(names(ks_fc), function(x) ks_fc[[x]][rownames(ks_fc[[x]]) %in% rownames(ks_fdrcut[[x]]), ])
ks_fccut <- lapply(ks_fccut, function(x) x[abs(rowSums(x)) > 0, ])


#Combine FDR df from all modules into one df
all_ks <- bind_rows(lapply(ks_fdrcut, function(x) as.data.frame(cbind("gene" = rownames(x), x))))
rownames(all_ks) <- all_ks$gene
all_ks <- all_ks[, -1]
gene_names <- rownames(all_ks)
all_ks <- sapply(all_ks, function(x) as.numeric(as.character(x)))
rownames(all_ks) <- gene_names

#Annotate each gene to corresponding module
anno_row <- bind_rows(lapply(names(ks_fdrcut), function(x) as.data.frame(cbind("gene" = rownames(ks_fdrcut[[x]]), "module" = x) )))
rownames(anno_row) <- anno_row$gene
anno_col <- list(celltype = ggcolours)

#Calculate breaks in heatmap based on modules
indices <- sapply(ks_fdrcut, nrow)
for(i in 2:length(indices)) {
  x <- as.numeric(indices[i])
  indices[i] <- as.numeric(indices[i - 1]) + as.numeric(indices[i])
}

cell_anno <- as.data.frame(colnames(all_ks))
colnames(cell_anno) <- "celltype"
rownames(cell_anno) <- cell_anno$celltype

use_breaks <- get_htmap_breaks(all_ks)

#+ ks_heatmap_manual, warning = F, message = F, width = 8, height = 10, fig.cap = "Cell-types are manually annotated."
p <- pheatmap(all_ks,
              # border_color = "black",
              color = colorRampPalette(c("darkmagenta", "snow2"))(250),
              annotation_col = cell_anno[, 1, drop = F],
              cluster_rows = F,
              cluster_cols = T,
              # scale = "row",
              annotation_colors = anno_col,
              gaps_row = indices,
              show_rownames = F
              # labels_col = paste(cell_anno$cluster, cell_anno$celltype), 
              # breaks = use_breaks2,
)
save_phtmap_pdf(p, "output_files/pit_utr_module_cheung_scrna_ks_FDR_cellenrichment_manual_anno.pdf", 8, 10)
print(p)
tmp <- dev.off()

#' ## Make signature matrix for CIBERSORT
#' Signature matrix contains gene counts averaged across all the cells for a given manually annotated cell cluster.
#' A factor of num_cells is used to ensure all the counts are integers (as required for CIBERSORT input).
#' Signature matrix is then used with CIBERSORT script (Newman et al 2015) run on the CIBERSORT browser application (https://cibersort.stanford.edu/).
#' CIBERSORT script is run with UTR-seq raw counts (after RUV-seq correction) and signature matrix with the following parameters:
#' * permutations n = 100
#' * quantile normalization = F
#' * absolute mode = F
#' * absolute method = sig.score.
#' For simplicity of running this script in one run, the CIBERSORT results are included in the input folder.
#+ make_sig_mat, warning = F, message = F
#### Make signature matrix ####
# with all gene expression and manually annotated cell clusters

pitman_countmat <- pit_seurat[["RNA"]]@counts
colnames(pitman_countmat) <- pit_seurat@active.ident

celltypes <- levels(pit_seurat)
man_mean_mat <- sapply(celltypes, function(x) fill_seurat_mat_mean(x, pit_seurat))
# man_mean_mat <- man_mean_mat*1000
man_mean_mat <- man_mean_mat*length(pit_seurat$nCount_RNA)
man_mean_mat <- apply(man_mean_mat, c(1,2), as.integer)
colnames(man_mean_mat) <- gsub(" ", "_", colnames(man_mean_mat))
colnames(man_mean_mat) <- gsub("/", "_", colnames(man_mean_mat))
man_mean_mat <- cbind("NAME" = rownames(man_mean_mat), man_mean_mat)
kable(man_mean_mat[1:5,], caption = "First 5 rows of signature matrix")

write.table(man_mean_mat, "output_files/seurat_v4_Cheungpit_manual_cluster_MEAN_gene_signature_matrix_2022-03-03.txt",
            sep = "\t", quote = F, row.names = F, col.names = T)

#' ## Plot CIBERSORT cell proportions
#' Estimated cell proportions from CIBERSORT (https://cibersort.stanford.edu/) are plotted for each cell-type.
#' For each cell-type, proportions are separated by age and sex.
#' Wilcoxon test is performed between males and females at each age.
#+ plot_cibersort, message = F, warning = F, fig.width = 14, fig.height = 7, fig.cap = "CIBERSORT estimated cell propotions."
#### Plot CIBERSORT cell proportions ####
# CIBERSORT was run on the browser application

cib_res <- read.table("input_files/cibersort_browser_cheung_mean_rawcounts_2022-03-03/CIBERSORT.cibersort_browser_cheung_v4_mean_rawcounts_2022-03-03.txt",
                      sep = "\t", header = T)

cib_plot <- tbl_df(melt(dplyr::select(cib_res, -Pearson.Correlation, -RMSE, -`P.value`)))
colnames(cib_plot) <- c("sample", "cell_type", "prop")
cib_plot <- mutate(cib_plot, age = substr(sample, 10, 11),
                   sex = substr(sample, 12, 12), rep = substr(sample, 13, 13))
cib_plot$age <- as.numeric(cib_plot$age)

gg_def_pal <- hue_pal()(length(levels(cib_plot$cell_type)))
lighter_cluster_colours <- sapply(gg_def_pal, function(x) lighten(x, amount = 0.5))
cluster_colours <- sapply(gg_def_pal, function(x) darken(x, amount = 0.2))
use_colors <- lapply(1:length(cluster_colours), function(x) c("F" = lighter_cluster_colours[[x]], "M" = cluster_colours[[x]]))
names(use_colors) <- levels(cib_plot$cell_type)

cell_prop_plots <- lapply(levels(cib_plot$cell_type), function(x) get_cell_prop_plot(x, cib_plot, use_colors))
p <- ggarrange(plotlist = cell_prop_plots)
pdf("output_files/cibersort_browser_cell_proportions_raw_counts_cheung_v4_2022-03-03.pdf", width = 14, height = 7, useDingbats = F)
print(p)
dev.off()

print(p)

#' ## R Session Info
#+ r_session, message = F, warning = F
sessionInfo()
