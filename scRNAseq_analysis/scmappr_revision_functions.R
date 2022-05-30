run_scmappr <- function(comparison, case, ctrl, file_out) {
  print(case)
  print(ctrl)
  use_bulk_DE <- bulk_de[[comparison]][, c(7, 6, 2)]
  
  if(nrow(use_bulk_DE) > 0) {
    scmappr_out <- scMappR_and_pathway_analysis(bulk_counts, wilcoxon_scmappr$OR,
                                                use_bulk_DE, case_grep = case,
                                                control_grep = ctrl,
                                                max_proportion_change = 10, print_plots = TRUE,
                                                plot_names = paste0("scMappR_revision_pituitary_", comparison), theSpecies = "mouse",
                                                output_directory = paste0("scMappR_revision_pituitary_", comparison),
                                                up_and_downregulated = TRUE, 
                                                internet = TRUE, toSave = TRUE,
                                                path = file_out,
                                                deconMethod = "WGCNA")
  }
  else {print(paste0(comparison, " has no bulk DE genes"))}
}

plot_cwfc <- function(compar, scrna, genes_expressed, scrna_data, scrna_markers) {
  cwfc <- get(load(paste0("output/", scrna, "_scmappr/scMappR_revision_pituitary_d", compar, "_sex/scMappR_revision_pituitary_d", compar, "_sex_cellWeighted_Foldchanges.Rdata")))
  props <- get(load(paste0("output/", scrna, "_scmappr/scMappR_revision_pituitary_d", compar, "_sex/scMappR_revision_pituitary_d", compar, "_sex_celltype_proportions.RData")))
  dds_significant <- bulk_de[[paste0("d", compar, "_sex")]]
  rownames(dds_significant) <- dds_significant$genename
  # format bulk DE gene list to: c("gene_name", "padj", "log2fc")
  dds_significant <- dplyr::select(dds_significant, genename, FDR, logFC)
  cwFC_eval <- cwFoldChange_evaluate(cwfc, props, dds_significant)
  
  cwFC_genes_assigned <- unique(unlist(lapply(cwFC_eval$cwFoldchange_gene_assigned, function(x) names(x))))
  cwFC_genes_FP <- unique(unlist(lapply(cwFC_eval$cwFoldchange_gene_flagged_FP, function(x) names(x))))
  cwFC_genes_expressed <- which(rownames(cwFC_eval$cwFoldChange_normalized) %in% genes_expressed) # Filter for genes which are expressed in the single cell dataset
  cwFC_df <- cwFC_eval$cwFoldChange_normalized[cwFC_genes_expressed,]
  if(length(cwFC_genes_FP) > 0) { # Remove genes flagged as false positive by scMappR
    cwFC_genes_FP <- which(rownames(cwFC_df) %in% cwFC_genes_FP)
    cwFC_df <- cwFC_df[-cwFC_genes_FP,]
  }
  cwFC_df <- cwFC_df[rownames(cwFC_df) %in% cwFC_genes_assigned,drop=F,]
  
  htmap_fxn(cwFC_df, paste0("output/", scrna, "_scmappr/scMappR_revision_pituitary_d", compar, "_sex/scMappR_revision_pituitary_d", compar, "_sex_scmappr_cwfc_heatmap_full.pdf"))
  
  cwFC_df_05 <- cwFC_df[which(apply(cwFC_df, 1, function(x) max(abs(x)) > 0.5)),drop = F,]
  htmap_fxn(cwFC_df_05,
            paste0("output/", scrna, "_scmappr/scMappR_revision_pituitary_d", compar, "_sex/scMappR_revision_pituitary_d", compar,
                   "_sex_scmappr_cwfc_heatmap_cwFC_05.pdf"))
  
  cwFC_df_hormone <- cwFC_df[, colnames(cwFC_df) %in% c("Somatotropes", "Lactotropes", "Corticotropes",
                                                        "Gonadotropes", "Melanotropes",
                                                        "Thyrotropes"), drop = F]
  cwFC_df_hormone <- cwFC_df_hormone[which(abs(rowSums(cwFC_df_hormone)) > 0.00001),drop = F,] # Filter out genes which had high cwFC in the non-hormone producing cell types which are filtered out
  
  htmap_fxn(cwFC_df_hormone,
            paste0("output/", scrna, "_scmappr/scMappR_revision_pituitary_d", compar, "_sex/scMappR_revision_pituitary_d", compar,
                   "_sex_scmappr_cwfc_heatmap_hormone.pdf"))
  
  cwFC_df_hormone_05 <- cwFC_df_hormone[which(apply(cwFC_df_hormone, 1, function(x) max(abs(x)) > 0.5)),drop = F,]
  htmap_fxn(cwFC_df_hormone_05,
            paste0("output/", scrna, "_scmappr/scMappR_revision_pituitary_d", compar, "_sex/scMappR_revision_pituitary_d", compar,
                   "_sex_scmappr_cwfc_heatmap_hormone_cwFC_05.pdf"))
 
  cwFC_expression_20 <- which(rowSums(scrna_data@assays$SCT[rownames(cwFC_df),]) >= 20)  # Filter for genes with sum(SCT_expr) >= 20
  
  cwFC_df <- cwFC_df[cwFC_expression_20,drop = F,]
  htmap_fxn(cwFC_df,
            paste0("output/", scrna, "_scmappr/scMappR_revision_pituitary_d", compar, "_sex/scMappR_revision_pituitary_d", compar,
                   "_sex_scmappr_cwfc_heatmap_cwFC_sct_expr_20.pdf"))
  
  cwFC_df_expr_05 <- cwFC_df[which(apply(cwFC_df, 1, function(x) max(abs(x)) > 0.5)),drop = F,]
  htmap_fxn(cwFC_df_expr_05,
            paste0("output/", scrna, "_scmappr/scMappR_revision_pituitary_d", compar, "_sex/scMappR_revision_pituitary_d", compar,
                   "_sex_scmappr_cwfc_heatmap_cwFC_sct_expr_20_cwFC_05.pdf"))
  
  out_cwFC <- lapply(cwFC_eval$cwFoldchange_gene_assigned, function(x) as.data.frame(x))
  out_cwFC <- lapply(out_cwFC, function(x) x[rownames(x) %in% rownames(cwFC_df),,drop=F])
  out_cwFC <- lapply(out_cwFC, function(i) dplyr::rename(i, "cwFC" = x) %>%
                       mutate(genename = rownames(i)))
  out_cwFC <- bind_rows(lapply(names(out_cwFC), function(x) mutate(out_cwFC[[x]], celltype = x)))
  write.table(out_cwFC, paste0("output/", scrna, "_scmappr/scMappR_revision_pituitary_d", compar, "_sex/scMappR_revision_pituitary_d", compar, "_sex_scmappr_cwfc_vals.txt"),
              quote = F, sep = "\t", col.names = T, row.names = F)
  
  # Filter for genes which are cell type markers as defined by FindAllMarkers function
  # in Seurat (with only positive markers)
  # Filter also for cwFC > 0.5
  # ct_specific <- lapply(levels(scrna_markers$cluster), function(x) prop_scmappr(x, filter(out_cwFC, abs(cwFC) > 0.5), dds_significant))
  # names(ct_specific) <- levels(scrna_markers$cluster)
  # 
  # ct_specific_use <- lapply(ct_specific, function(x) x[[1]])
  # ct_specific_use <- lapply(names(ct_specific_use), function(x) mutate(ct_specific_use[[x]], celltype = x))
  # ct_specific_use <- bind_rows(ct_specific_use)
  # 
  # htmap_fxn(cwFC_df_05[unique(ct_specific_use$genename),],
  #           paste0("output/", scrna, "_scmappr/scMappR_revision_pituitary_d", compar, "_sex/scMappR_revision_pituitary_d", compar,
  #                  "_sex_scmappr_cwfc_heatmap_cwFC_sct_cwFC_05_celltype_specific_genes.pdf"))
  # ct_specific_use <- mutate(ct_specific_use, celltype_gene = paste0(genename, "_", celltype))
  # out_cwFC_use <- mutate(out_cwFC, celltype_gene = paste0(genename, "_", celltype))
  # out_cwFC_use <- filter(out_cwFC_use, celltype_gene %in% ct_specific_use$celltype_gene)
  # write.table(out_cwFC_use, paste0("output/", scrna, "_scmappr/scMappR_revision_pituitary_d", compar, "_sex/scMappR_revision_pituitary_d", compar, "_sex_scmappr_cwfc_05_vals_celltype_marker.txt"),
  #             quote = F, sep = "\t", col.names = T, row.names = F)
  # cwFC_df_05_celltype_specific <- cwFC_df_05[unique(ct_specific_use$genename),c("Somatotropes", "Lactotropes", "Gonadotropes")]
  # cwFC_df_05_celltype_specific <- cwFC_df_05_celltype_specific[which(abs(rowSums(cwFC_df_05_celltype_specific)) > 0.15),drop = F,] # Filter out genes which had high cwFC in the non-hormone producing cell types which are filtered out
  # htmap_fxn(cwFC_df_05_celltype_specific,
  #           paste0("output/", scrna, "_scmappr/scMappR_revision_pituitary_d", compar, "_sex/scMappR_revision_pituitary_d", compar,
  #                  "_sex_scmappr_cwfc_heatmap_cwFC_sct_cwFC_05_celltype_specific_genes_ossd.pdf"))
  # 
}

prop_scmappr <- function(cell_type, cwfc_df, de_bulk) {
  cwfc_sub <- filter(cwfc_df, celltype == cell_type)
  marker_sub <- filter(ruf_markers_sct_pos, cluster == cell_type)
  bulk_sub <- filter(de_bulk, genename %in% marker_sub$gene)
  
  cw_in_cell <- cwfc_sub[which(cwfc_sub$genename %in% bulk_sub$genename),]
  cw_not_cell <- cwfc_sub[-which(cwfc_sub$genename %in% bulk_sub$genename),]
  
  return(list("CT_specific" = cw_in_cell, "not_CT_specific" = cw_not_cell))
}


get_htmap_breaks <- function(plotdata, use_pal) {
  paletteLength <- length(use_pal)
  myBreaks <- c(seq(min(plotdata), 0, length.out=ceiling(paletteLength/2) + 1), 
                seq(max(plotdata)/paletteLength, max(plotdata), length.out=floor(paletteLength/2)))
  custome_breaks <- unique(c(min(myBreaks), seq(quantile(myBreaks, 0.2), 0, length = paletteLength/2), seq(0, quantile(myBreaks, 0.8), length=paletteLength/2), max(myBreaks)))
  return(custome_breaks)
}

htmap_fxn <- function(df, outname) {
  if(nrow(df) == 0) {
    return(NA)
  }
  use_breaks <- get_htmap_breaks(df,
                                 colorRampPalette(c("orange", "white", "darkmagenta"))(255))
  if(nrow(df) == 1) {
    use_cluster <- F
  } else{use_cluster <- T}
  htmap <- pheatmap(t(df),
                    color = colorRampPalette(c("orange", "white", "darkmagenta"))(255),
                    # border = NA,
                    cluster_rows = F,
                    cluster_cols = use_cluster,
                    breaks = use_breaks,
                    filename = outname,
                    cellheight = 10, cellwidth = 8)
  return(htmap)
}
