#' ## Libraries
#+ libraries
library(scMappR)
library(Seurat)
library(ggplot2)
library(ggpubr)
library(reshape)
library(dplyr)
library(knitr)
library(WGCNA)
library(colorspace)
library(scales)
library(EDASeq)
source("decon_compare_functions.R")

#' ## Determine reference cell-type proportions
#' Load and save scRNA-seq cell-type proportions as a reference dataset.
#' This dataset was integrated using the script located at:
#' ~/⁨Dropbox (Wilson Lab)⁩/⁨Mike_Anna_Mark⁩/active_manuscript⁩/Pituitary paper⁩/scripts⁩/⁨single_cell_rmd⁩/20211214_integrated_process_ruf_snrnaseq.R
#+ load_scrna
dir.create("output/ruf_integrate", recursive = T)
ruf <- readRDS("input/20211214_ruf_2021_pit_GSE151961_integrated_seurat.rds")
# ruf <- readRDS("~/Dropbox (Wilson Lab)/Mike_Anna_Mark/active_manuscript/Pituitary paper/scripts/single_cell_rmd/GSE151961_RAW/output/20211214_ruf_2021_pit_GSE151961_integrated_seurat.rds")

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
p1 <- DimPlot(ruf, reduction = "umap", group.by = "sex")
p2 <- DimPlot(ruf, reduction = "umap", group.by = "int_clusters", label = T, repel = T)

p1 + p2
ggsave("output/ruf_integrate/20220315_integrated_umap_clusters_simple_labelled.pdf", width = 20, height = 10,
       useDingbats = F)
ggsave("output/ruf_integrate/20220315_integrated_umap_clusters_simple_labelled.png", width = 20, height = 10)

ggsave("output/ruf_integrate/20220315_integrated_umap_clusters_cluster_only.pdf", p2, width = 10, height = 10,
       useDingbats = F)
ggsave("output/ruf_integrate/20220315_integrated_umap_clusters_cluster_only.png", p2, width = 10, height = 10)

markers.to.plot <- c("Gh", "Ghrhr", "Dlk1", "Prl", "Lhb", "Fshb", "Tshb", "Pomc", "Crhr1", "Pax7",
                     "Oacyl", "Sox2", "Rbpms", "Pou1f1", "Top2a", "Mki67", "Plvap",
                     "Meis2", "Emcn", "Dcn", "Lama2", "C1qa", "Ctss", "Scn7a", "Col25a1", "Sgcz")
DotPlot(ruf,
        # assay = "SCT",
        features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8, split.by = "sex") +
  RotatedAxis()
ggsave("output/ruf_integrate/ruf_integrated_dotplot_markers.pdf", width = 10, height = 10)

neuron_genes <- read.table("input/neuron_part_genes.txt", header = F)
DotPlot(pit,
        assay = "SCT",
        features = neuron_genes$V1, cols = c("tomato", "steelblue"), dot.scale = 8, split.by = "sex") +
  # coord_flip() +
  RotatedAxis()
ggsave("output/ruf_integrate/ruf_integrated_dotplot_neuron_genes_sct.pdf", width = 16, height = 10)
ggsave("output/ruf_integrate/ruf_integrated_dotplot_neuron_genes_sct.png", width = 16, height = 10)


# Use integrated data for deconvolution
# Get average proportions from scRNA-seq data for each dataset
pit <- ruf
rm(ruf)
pit <- subset(pit, idents = "Debris", invert = T) # remove debris cluster
ruf_F <- pit[,(pit$sex == "F")]
ruf_M <- pit[,(pit$sex == "M")] 

total_cells_scRNA <- table(Idents(pit))
total_cells_F <- table(Idents(ruf_F))
total_cells_M <- table(Idents(ruf_M))

props_scRNA <- total_cells_scRNA/sum(total_cells_scRNA)
props_scRNA_F <- total_cells_F/sum(total_cells_F)
props_scRNA_M <- total_cells_M/sum(total_cells_M)
rm(list=c("ruf_F","ruf_M"))

#### Run DeconRNA-seq, WGCNA, DCQ ####
# Generate signature matrix for Adapts
pit_utr <- readRDS("input/pit_utr_2019__RUV_k2_set1_2019-07-03.rds")
pit_norm <- normCounts(pit_utr)

pit_markers_sct <- FindAllMarkers(pit, assay = "SCT", only.pos = F, recorrect_umi = F)

gene_res <- lapply(levels(pit_markers_sct$cluster), function(x) pit_markers_sct[pit_markers_sct$cluster == x, ])
names(gene_res) <- levels(pit_markers_sct$cluster)
gene_res <- lapply(gene_res, function(x) {rownames(x) <- as.character(x$gene); x})
gene_res <- lapply(gene_res, function(x) x[, c(5, 2)])

wilcoxon_scmappr <- generes_to_heatmap(gene_res, species = "mouse", make_names = F)
kable(head(wilcoxon_scmappr$OR), caption = "Example gene OR in signature matrix.")

HM <- wilcoxon_scmappr$OR
signatureVar <- rowVars(HM)
signature_o <- HM[order(signatureVar, decreasing = TRUE),]
signature_top3k <- signature_o[1:3000,]

# Normalized bulk counts
DeconMethod <- compare_deconvolution_methods(as.data.frame(pit_norm), as.data.frame(signature_top3k))
DeconRNAseq <- DeconMethod$cellType_proportions$DeconRNAseq
WGCNA <- DeconMethod$cellType_proportions$WGCNA
DCQ <- DeconMethod$cellType_proportions$DCQ
rownames(DeconRNAseq) <- colnames(pit_norm)

#' ## Compare methods
#' Plot cell proportion figures for each method.
#' Compare estimated cell proportions with "ground truth" proportion calculated from single-cell data.
#+ compare_methods
celltypes <- list(WGCNA = colnames(WGCNA),
                  DCQ = colnames(DCQ),
                  DeconRNAseq=colnames(DeconRNAseq))

samples <- list(WGCNA = rownames(WGCNA),
                DCQ = rownames(DCQ),
                DeconRNAseq=rownames(DeconRNAseq))

# Get 'overlapping' samples and cell-types to be able to compare them.
cellEstimated <- Reduce(intersect,celltypes)
# Reorder to match UMAP order
cellEstimated <- cellEstimated[match(c("Somatotropes", "Lactotropes", "Corticotropes", "Gonadotropes",
                                       "Stem_cell", "Endothelial", "Proliferating", "Melanotropes",
                                       "Pericytes", "Macrophages", "Pituicytes", "Thyrotropes", "Somato_Lacto"), cellEstimated)]
samplesEstimated <- Reduce(intersect,samples) 

WGCNA_flt <- WGCNA[samplesEstimated,cellEstimated]
DCQ_flt <- DCQ[samplesEstimated,cellEstimated]
DeconRNAseq_flt <- DeconRNAseq[samplesEstimated,cellEstimated]
estimated_proportions <- list(WGCNA=as.data.frame(WGCNA_flt),
                              DCQ=as.data.frame(DCQ_flt),
                              DeconRNAseq=as.data.frame(DeconRNAseq_flt))
props_scRNA <- props_scRNA[cellEstimated]
props_scRNA_F <- props_scRNA_F[cellEstimated]
props_scRNA_M <- props_scRNA_M[cellEstimated]
saveRDS(list(estimated=estimated_proportions,
             sc_props=props_scRNA,
             f_props = props_scRNA_F,
             m_props = props_scRNA_M),"output/ruf_integrate/ruf_integrate_different_method_estimated_proportions.rds")

#' ### Plot proportions
#+ plot_props
gg_def_pal <- hue_pal()(length(cellEstimated))
lighter_cluster_colours <- sapply(gg_def_pal, function(x) lighten(x, amount = 0.5))
cluster_colours <- sapply(gg_def_pal, function(x) darken(x, amount = 0.2))
use_colors <- lapply(1:length(cluster_colours), function(x) c("F" = lighter_cluster_colours[[x]], "M" = cluster_colours[[x]]))
names(use_colors) <- cellEstimated

meta <- pData(pit_utr) %>%
  mutate(sample=rownames(pData(pit_utr)),
         age = gsub("d", "", age))

melt_props <- lapply(estimated_proportions, function(x) props_melt(x))
prop_plots <- lapply(melt_props, function(decon_method) ggarrange(plotlist = lapply(cellEstimated, function(celltype) get_cell_prop_plot(celltype, decon_method, use_colors))))
names(prop_plots) <- names(estimated_proportions)

lapply(names(prop_plots), function(x) ggsave(paste0("output/ruf_integrate/ruf_integrate_v4_proportion_plots_", x, ".pdf"),
                                             prop_plots[[x]], width = 14 , height = 9, useDingbats = F))

