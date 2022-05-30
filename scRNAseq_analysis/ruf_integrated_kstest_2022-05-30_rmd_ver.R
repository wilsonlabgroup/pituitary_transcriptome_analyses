#' ---
#' title: "Enrichment of coexpression module genes in Ruf-Zamojski 2021 sex-integrated single-nuclei RNA-seq dataset"
#' author: "Cadia Chan"
#' affiliation: "SickKids Research Institute & University of Toronto"
#' date: "May 30, 2022"
#' output:
#'  html_document:
#'    code_folding: hide
#'    toc: true
#'    toc_float: true
#' ---
#' 
#' Using the snRNA-seq data (GSE151961) from Ruf-Zamojski 2021 (https://www.nature.com/articles/s41467-021-22859-w) which we integrated between sexes,
#' enrichment for genes in co-expression modules within different cell-types identified by snRNA-seq data is assessed.
#' Enrichment is calculated based on expression of a given gene within a cell type compared to all other cell types using a one-tailed Kolmogorov-Smirnov (KS) test.
#' A gene is considered enriched if FDR < 0.05 and FC > 0 in a cell type compared to other cell types.
#' 
#' To run: Set working directory "To Source File Location".
#'

#' ## Libraries and source scripts
#+ load_library, warning = F, message = F
#### Libraries ####
library(dplyr)
library(Seurat)
library(knitr)
# library(biomaRt)
# library(igraph)
# library(ggraph)
# library(intergraph)
# library(rowr)
library(scales)
library(EDASeq)
library(pheatmap)
library(LaCroixColoR)
library(grid)
library(plyr)
# library(ggpubr)
library(reshape2)
library(ggplot2)
library(colorspace)
library(ggrepel)
# library(openxlsx)

source("ruf_integrated_kstest_2022-05-30_functions_rmd_ver.R")

#' ## Load in integrated dataset
#+ load_data, warning = F, message = F
#### Load in scRNA-seq data ####
ruf <- readRDS("input/20211214_ruf_2021_pit_GSE151961_integrated_seurat.rds")

# Group cell-types
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
ruf <- subset(ruf, idents = "Debris", invert = T)

#' # Enrichment of module genes in single cell populations
#' ## Setup data for enrichment test
#' Co-expression module gene data is loaded in.
#' A colour palette is made based on the cell-types.
#+ setup_modules_data, warning = F, message = F
#### Enrichment of module genes in single cell populations ####
clusters <- factor(ruf$int_clusters, levels = levels(ruf$int_clusters)[1:length(levels(ruf$int_clusters))-1])
mod_utr <- readRDS("input/coexpression_modules_202205.rds")
mod_utr$modules <- factor(mod_utr$modules)

celltypes <- levels(clusters)

# Create colour palette for heatmap annotation of cell types
ggcolours <- hue_pal()(length(celltypes))
names(ggcolours) <- celltypes

#' ## Run KS test for combined cell clusters
#' KS test is run and FDR is calculated for each module gene in for combined cell clusters.
#' Genes with enriched with FC > 0, FDR adj < 0.05 in at least one cell type compared to all other cell types are retained.
#' Heatmap is plotted with FDR value for each module gene retained in each cell type.
#+ run_ks_test, warning = F, message = F
#Run KS test and calculate FDR for each module gene
ks_res <- lapply(levels(mod_utr$modules), function(x) calc_stat_test(mod_utr[mod_utr$modules == x, "genename"], "ks", ruf, clusters))
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
anno_col <- list(celltypes = ggcolours)

#Calculate breaks in heatmap based on modules
indices <- sapply(ks_fdrcut, nrow)
for(i in 2:length(indices)) {
  x <- as.numeric(indices[i])
  indices[i] <- as.numeric(indices[i - 1]) + as.numeric(indices[i])
}

cell_anno <- as.data.frame(celltypes)
rownames(cell_anno) <- cell_anno$celltypes

use_breaks <- get_htmap_breaks(all_ks)

ks_org <- lapply(1:length(indices), function(x) reorg_ks(x, all_ks))
ks_org_df <- bind_rows(ks_org)

#+ num_enrich, warning = F, message = F
sapply(1:length(ks_org), function(x) num_enrich(x))


#+ ks_heatmap_all_clusters, warning = F, message = F, width = 8, height = 10, fig.cap = "Heatmap showing enrichment of module genes in cell-types based on FDR-adjusted P-values from a KS test."
p <- pheatmap(ks_org_df,
              # border_color = "black",
              color = colorRampPalette(c("darkmagenta", "snow2"))(250),
              annotation_col = cell_anno[, 1, drop = F],
              cluster_rows = F,
              cluster_cols = F,
              # scale = "row",
              annotation_colors = anno_col,
              gaps_row = indices,
              show_rownames = F,
              labels_col = paste(cell_anno$cluster, cell_anno$celltype) 
              # breaks = use_breaks2,
)
save_phtmap_pdf(p, "output/pit_utr_module_cheung_scrna_ks_FDR_cellenrichment_no_debris.pdf", 8, 10)
print(p)
tmp <- dev.off()

#+ add_genelabels, warning = F, message = F, width = 8, height = 10, fig.cap = "Heatmap with top 5 hub genes labelled showing enrichment of module genes in cell-types based on FDR-adjusted P-values from a KS test."
hub_genes <- readRDS("input/top10_hub_genes_202205.rds")
hub_genes <- melt(apply(hub_genes, 1, function(x) unlist(strsplit(x, ",")))[1:5,])
anno_row <- mutate(anno_row, gene = ifelse(gene %in% hub_genes$value, paste0("----", gene), ""))

p <- pheatmap(ks_org_df,
              # border_color = "black",
              color = colorRampPalette(c("darkmagenta", "snow2"))(250),
              annotation_col = cell_anno[, 1, drop = F],
              cluster_rows = F,
              cluster_cols = F,
              # scale = "row",
              annotation_colors = anno_col,
              gaps_row = indices,
              # show_rownames = F,
              labels_col = paste(cell_anno$cluster, cell_anno$celltype),
              labels_row = anno_row$gene
              # breaks = use_breaks2,
)
save_phtmap_pdf(p, "output/pit_utr_module_cheung_scrna_ks_FDR_cellenrichment_no_debris_top5_label.pdf", 8, 10)
print(p)
tmp <- dev.off()

#' ## R Session Info
#+ r_session, message = F, warning = F
sessionInfo()
