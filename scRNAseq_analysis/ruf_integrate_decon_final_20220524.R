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
source("decon_compare_functions.R")

#' ## Determine reference cell-type proportions
#' Load and save scRNA-seq cell-type proportions as a reference dataset.
#' This dataset was integrated using the script located at:
#' ~/⁨Dropbox (Wilson Lab)⁩/⁨Mike_Anna_Mark⁩/active_manuscript⁩/Pituitary paper⁩/scripts⁩/⁨single_cell_rmd⁩/20211214_integrated_process_ruf_snrnaseq.R
#+ load_scrna
dir.create("output/ruf_integrate")
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

#' ## CIBERSORT and CIBERSORTx
#' See cheung_sc_reprocess_2022-03-02_rmd_ver.R for generation of signature matrix (bottom of script).
#' Also try running it with CIBERSORT generated signature matrix. To do this, need to generate
#' reference count matrix (for both) and phenotype classes file (only for CIBERSORT).
#' Problems with uploading reference count matrix to CIBERSORT (>500 MB limit), therefore, use
#' signature matrix generated from CIBERSORTx with CIBERSORT after.
# Cibersort and Cibersort X data processing:

#' Processing in the scRNA-seq count matrix for CIBERSORTx
# norm_counts <- as.data.frame(pit_seurat@assays$RNA@counts)
# colnames(norm_counts) <- Idents(pit_seurat)
# write.table(norm_counts, file = "output/cheung_et_al_v4_raw_nolog_sc_ref_mat.txt", quote = FALSE, row.names = TRUE, col.names = TRUE, sep = "\t")

#' Generate phenotype classes file.
# Making the pheno file, rows are celltypes, columns are cells.
# Phenotype file is not used with CIBERSORTx.
# 1 is Yes, 2 is no
ID <- Idents(pit)
pheno <- matrix(2, length(levels(ID)),ncol(pit))
rownames(pheno) <- levels(ID)
for(i in levels(ID)){
  pheno[i,which(i==ID)] <- 1
}
write.table(pheno, file = "output/ruf_integrate/ruf_integrated_phenotype_files_2022-03-16.txt", quote = FALSE, row.names = TRUE, col.names = FALSE, sep = "\t")

celltypes <- levels(pit)
# Also try with all markers (neg and pos) identified by FindAllMarkers function (5435 total)
# Try with integrated@data
ciber_mat <- as.data.frame(pit@assays$integrated@data[rownames(pit@assays$integrated@data) %in% rownames(HM),])

# ciber_mat <- as.data.frame(pit@assays$integrated@data[rownames(pit@assays$integrated@data) %in% rownames(HM),])
ciber_mat <- 2^ciber_mat
ciber_mat <- round(ciber_mat, 2)
colnames(ciber_mat) <- Idents(pit)
write.table(cbind("NAME"=rownames(ciber_mat),ciber_mat), "output/ruf_integrate/ruf_integrated_no_log_data_all_markers_sc_ref_2022-03-16.txt",
            sep = "\t", quote = F, row.names = F, col.names = T)
# rm(ciber_mat)

###
# We used the web tool with default parameters. The analysis continues from the output.
# >[Options] perm: 100
# >[Options] QN: FALSE
###

Cibersort <- read.table("cibersort/ruf_integrate/markers/data/CIBERSORT.Output_Job29.txt",header = T,as.is=T,sep = "\t")
rownames(Cibersort) <- Cibersort$Input.Sample
Cibersort_props <- Cibersort[,intersect(colnames(Cibersort),names(total_cells_scRNA))]

Cibersortx <- read.table("cibersortX/ruf_integrated/markers/data/CIBERSORTx_Job24_Results.txt",header=T,as.is=T,sep="\t")
rownames(Cibersortx) <- Cibersortx$Mixture
Cibersortx_props <- Cibersortx[,intersect(colnames(Cibersortx),names(total_cells_scRNA))] # excludes pValue RMSE columns

#' ## CPM Deconvolution
#' Set up files in this script locally then export and import to script on cluster.
#' Run this script which is commented out at the bottom on the cluster.
#' This script takes a long time (we ran it on a cluster), and it is how we got the results for CPM.
#' 
# Note that on cluster I could only install scBio v0.1.4 and so the CPM command does not have typeTransformation to specify.
# This did not run on the DC cluster, trying locally (with v0.1.6). This took about 8 hours to run.
#+ cpm_decon
library(scBio)
Labels <- tochr(Idents(pit))
norm_counts <- as.data.frame(pit@assays$integrated@data)
norm_counts <- 2^norm_counts
# colnames(raw_counts) <- Idents(pit_seurat)
# write.table(raw_counts, file = "output/cheung_et_al_v4_raw_nolog_sc_ref_mat.txt", quote = FALSE, row.names = TRUE, col.names = TRUE, sep = "\t")

bulk_norm <- normCounts(pit_utr)
# write.table(cbind(genename=rownames(bulk_norm), bulk_norm), "input/pit_utr_2019__RUV_k2_set1_2019-07-03_normcounts_countmat.txt",
# row.names = F, col.names = T, sep = "\t", quote = F)
SCcellspace <- pit@reductions$umap@cell.embeddings
# saveRDS(list(Labels, raw_counts, bulk_raw, SCcellspace), file = "output/cheung_v4_raw_CPM_2022-03-04.rds")

res <- CPM(SCData = norm_counts, SCLabels = Labels,
           BulkData = bulk_norm, cellSpace = SCcellspace, quantifyTypes=TRUE, typeTransformation=TRUE)
# res <- readRDS("output/ruf_integrate/ruf_integrate_v4_norm_CPM_props_2022-03-16.rds")
saveRDS(res, file = "output/ruf_integrate/ruf_integrate_v4_norm_CPM_props_2022-03-16.rds")
dir.create("CPM")
propCPM <- res$cellTypePredictions
propCPM <- propCPM[, match(colnames(propCPM), levels(Idents(pit)))]
saveRDS(res, "output/ruf_integrate/ruf_integrate_CPM_cell_props.rds")
write.table(cbind(sampleid=rownames(propCPM), propCPM), "CPM/integrate_v4_sctransform_CPM.txt",
            row.names = F, col.names = T, sep = "\t", quote = F)

#' ## MuSiC Deconvolution
#' This script doesn't take quite as long but requires a lot of dependencies and takes some computational space, and it is how we got the results for Music
#' we tried scaled and unscaled as well as SCTransform and integrated counts. i.a. worked the best
#' Dustin's script did not work for my datasets because it relies on package that specifies Seurat v3 object (mine are v4 to work with SCI-seq data),
#' therefore, I directly followed the MuSiC tutorial (https://xuranw.github.io/MuSiC/articles/MuSiC.html#music-deconvolution-1).
#' To run MuSiC I needed to split the dataset randomly into 2 pseudobulk datasets (this is only an issue if your single-cell dataset is one replicate).
# See https://github.com/xuranw/MuSiC/issues/49.
#+ music_decon
#source("Manuscript_code/Preprocessing/deconvolution_music.R")
library(BisqueRNA)
library(MuSiC)
library(SummarizedExperiment)
library(xbioc) 

# Convert bulk data into expression set for both raw and normalized counts
bulk_norm_eset <- ExpressionSet(assayData = normCounts(pit_utr),
                                phenoData = AnnotatedDataFrame(pData(pit_utr)))

sc_pdata <- as.data.frame(cbind(celltype=as.character(Idents(pit)),
                                cellid=names(Idents(pit)),
                                set = pit$orig.ident))
sc_pdata$celltype <- factor(sc_pdata$celltype)
sc_pdata$set <- factor(sc_pdata$set)

rownames(sc_pdata) <- sc_pdata$cellid
scdata <- 2^(as.data.frame(pit@assays$integrated@data))
scdata <- round(scdata, 2)
sc_norm_eset <- ExpressionSet(assayData = as.matrix(scdata),
                              phenoData = AnnotatedDataFrame(sc_pdata[,c(1,3), drop = F]))
# Run Music with norm_bulk
norm_bulk_prop <- music_prop(bulk.eset = bulk_norm_eset,
                             sc.eset = sc_norm_eset,
                             clusters = "celltype",
                             samples = "set")
saveRDS(norm_bulk_prop, "output/ruf_integrate/ruf_integrate_norm_music_nnls_cell_props.rds")
MuSiC <- norm_bulk_prop$Est.prop.weighted
MuSiC <- MuSiC[, match(colnames(MuSiC), levels(Idents(pit)))]

NNLS <- norm_bulk_prop$Est.prop.allgene
NNLS <- NNLS[, match(colnames(NNLS), levels(Idents(pit)))]

dir.create("music_nnls")

write.table(cbind(sampleid=rownames(MuSiC), MuSiC), "music_nnls/ruf_integrate_v4_sctransform_only_music.txt",
            row.names = F, col.names = T, sep = "\t", quote = F)
write.table(cbind(sampleid=rownames(NNLS), NNLS), "music_nnls/ruf_integrate_v4_sctransform_only_nnls.txt",
            row.names = F, col.names = T, sep = "\t", quote = F)

# Benchmark using "artificial" bulk data constructed
pit.construct.full = bulk_construct(sc_norm_eset, clusters = 'celltype', samples = 'set')
names(pit.construct.full)
pit.construct.full$Bulk.counts
head(pit.construct.full$num.real)

pit.construct.full$prop.real = relative.ab(pit.construct.full$num.real, by.col = FALSE)
head(pit.construct.full$prop.real)

# Estimate cell type proportions of artificial bulk data
Est.prop.pit = music_prop(bulk.eset = pit.construct.full$Bulk.counts, sc.eset = sc_norm_eset,
                          clusters = 'celltype', samples = 'set')

Eval_multi(prop.real = data.matrix(pit.construct.full$prop.real),
           prop.est = list(data.matrix(Est.prop.pit$Est.prop.weighted),
                           data.matrix(Est.prop.pit$Est.prop.allgene)),
           method.name = c('MuSiC', 'NNLS'))

library(cowplot)
prop.comp.fig = Prop_comp_multi(prop.real = data.matrix(pit.construct.full$prop.real),
                                prop.est = list(data.matrix(Est.prop.pit$Est.prop.weighted),
                                                data.matrix(Est.prop.pit$Est.prop.allgene)),
                                method.name = c('MuSiC', 'NNLS'),
                                title = 'Heatmap of Real and Est. Prop' )

abs.diff.fig = Abs_diff_multi(prop.real = data.matrix(pit.construct.full$prop.real),
                              prop.est = list(data.matrix(Est.prop.pit$Est.prop.weighted),
                                              data.matrix(Est.prop.pit$Est.prop.allgene)),
                              method.name = c('MuSiC', 'NNLS'),
                              title = 'Abs.Diff between Real and Est. Prop' )
pdf("output/ruf_integrate/music_nnls_benchmark.pdf", width = 12, height = 9)
plot_grid(prop.comp.fig, abs.diff.fig, labels = "auto", rel_widths = c(4,3))
dev.off()

#' ## Compare methods
#' Plot cell proportion figures for each method.
#' Compare estimated cell proportions with "ground truth" proportion calculated from single-cell data.
#+ compare_methods
celltypes <- list(MuSiC = colnames(MuSiC),
                  NNLS = colnames(NNLS),
                  CPM=colnames(propCPM),
                  Cibersort=colnames(Cibersort_props),
                  Cibersortx=colnames(Cibersortx_props),
                  WGCNA = colnames(WGCNA),
                  DCQ = colnames(DCQ),
                  DeconRNAseq=colnames(DeconRNAseq))

samples <- list(MuSiC = rownames(MuSiC),
                NNLS = rownames(NNLS),
                CPM=rownames(propCPM),
                Cibersort=rownames(Cibersort_props),
                Cibersortx=rownames(Cibersortx_props),
                WGCNA = rownames(WGCNA),
                DCQ = rownames(DCQ),
                DeconRNAseq=rownames(DeconRNAseq))

# Get 'overlapping' samples and cell-types to be able to compare them.
cellEstimated <- Reduce(intersect,celltypes)
# Reorder to match UMAP order
cellEstimated <- cellEstimated[match(c("Somatotropes", "Lactotropes", "Corticotropes", "Gonadotropes",
                                       "Stem_cell", "Endothelial", "Proliferating", "Melanotropes",
                                       "Pericytes", "Macrophages", "Pituicytes", "Thyrotropes", "Somato_Lacto"), cellEstimated)]
samplesEstimated <- Reduce(intersect,samples) 

MuSiC_flt <- MuSiC[samplesEstimated,cellEstimated]
NNLS_flt <- NNLS[samplesEstimated,cellEstimated]
CPM_flt <- propCPM[samplesEstimated,cellEstimated]
Cibersort_flt <- Cibersort[samplesEstimated,cellEstimated]
Cibersortx_flt <- Cibersortx[samplesEstimated,cellEstimated]
WGCNA_flt <- WGCNA[samplesEstimated,cellEstimated]
DCQ_flt <- DCQ[samplesEstimated,cellEstimated]
DeconRNAseq_flt <- DeconRNAseq[samplesEstimated,cellEstimated]
estimated_proportions <- list(MuSiC=as.data.frame(MuSiC_flt),
                              NNLS=as.data.frame(NNLS_flt),
                              CPM=as.data.frame(CPM_flt),
                              Cibersort=as.data.frame(Cibersort_flt),
                              Cibersortx=as.data.frame(Cibersortx_flt),
                              WGCNA=as.data.frame(WGCNA_flt),
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

#' ### Number of detected cells
#+ num_detected_cells
unlist(lapply(estimated_proportions, function(x) length(which(colSums(x) > 0))))

#### From decon_comparison_correlation_final_20220518.R ####
#+ load_props
#### Load proportions ####
ruf <- readRDS("output/ruf_integrate/ruf_integrate_different_method_estimated_proportions.rds")

#' Use mean value for estimated cell proportion from PD37 split by male and female
#' and correlate with proportions from female and male adult pituitary single-cell data (Ruf-Zamojski 2021).
#+ correlation
#### Correlation ####
ruf_scrna_F <- ruf$f_props
ruf_scrna_M <- ruf$m_props

dir.create("output/corr_plots")

ruf_corr_37 <- plot_corr(ruf, "Ruf-Zamojski_2021", 37, ruf_scrna_F, ruf_scrna_M)

#### Print correlation values ####
print(lapply(ruf_corr_37$cor, function(sex) sapply(sex, function(decon_method) decon_method$estimate)))
