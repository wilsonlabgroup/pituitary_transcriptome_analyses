get_cell_prop_plot <- function(cell_id, cell_prop_data, palette) {
  # print(cell_id)
  plot_prop <- filter(cell_prop_data, cell_type == cell_id)
  use_pal <- palette[[cell_id]]
  # plot_mean <- filter(cell_prop_mean_data, cell_type == cell_id)
  
  g <- ggline(plot_prop, x = "age", y = "prop", group = "sex",
              add = c("mean_se", "point"),
              col = "sex",
              palette = use_pal,
              size = 1,
              linetype = "sex",
              point.size = 2.5,
              shape = "sex",
              numeric.x.axis = T) +
    ggpubr::stat_compare_means(aes(group = sex), label = "p.signif", method = "wilcox.test", hide.ns = T, size = 6) +
    scale_x_continuous(breaks = c(12, 22, 27, 32, 37)) +
    facet_grid(.~cell_type) +
    labs(x = "Age", y = "Cell Proportion"
         # title = cluster
         # title = paste0("Cluster: ", substr(i, 2, 4)), subtitle = new_cluster_ids[index] # Use with unbiased clustering
    ) +
    theme_bw() +
    theme(
      strip.background.x = element_rect(fill = "grey85"),
      strip.text = element_text(size=14),
      axis.text = element_text(size=14 - 2, color="black"),
      axis.title = element_text(size=14),
      legend.position = "right", text = element_text(size = 14))
  g$layers[[1]]$aes_params$alpha <- 0.8
  g$layers[[2]]$aes_params$alpha <- 0.7
  g$layers[[3]]$aes_params$alpha  <- 0.7
  g$layers[[4]]$aes_params$alpha  <- 0.8
  # print(g)
  return(g)
}

make_decon_plots <- function(method, file_out) {
  decon_res <- decon_method$cellType_proportions[[method]]
  rownames(decon_res) <- rownames(decon_method$cellType_proportions$DCQ)
  decon_plot <- tbl_df(melt(decon_res))
  colnames(decon_plot) <- c("sample", "cell_type", "prop")
  decon_plot <- mutate(decon_plot, age = substr(sample, 10, 11),
                       sex = substr(sample, 12, 12), rep = substr(sample, 13, 13))
  decon_plot$age <- as.numeric(decon_plot$age)
  
  cell_prop_plots <- lapply(levels(decon_plot$cell_type), function(x) get_cell_prop_plot(x, decon_plot, use_colors))
  p <- ggarrange(plotlist = cell_prop_plots)
  pdf(paste0(file_out, method, "_cell_proportions.pdf"), width = 14, height = 7, useDingbats = F)
  print(p)
  dev.off()
  return(p)
}

run_scmappr <- function(comparison, case, ctrl, file_out) {
  print(case)
  print(ctrl)
  use_bulk_DE <- bulk_DE[[comparison]][, c(7, 6, 2)]
  
  if(nrow(use_bulk_DE) > 0) {
    scmappr_out <- scMappR_and_pathway_analysis(bulk_normalized, odds_ratio_in,
                                                use_bulk_DE, case_grep = case,
                                                control_grep = ctrl,
                                                max_proportion_change = 10, print_plots = TRUE,
                                                plot_names = paste0("scMappR_pituitary_", comparison), theSpecies = "mouse",
                                                output_directory = paste0("scMappR_pituitary_", comparison),
                                                up_and_downregulated = TRUE,
                                                internet = TRUE, toSave = TRUE,
                                                path = file_out)
  }
  else {print(paste0(comparison, " has no bulk DE genes"))}
}

cwFC_plot <- function(plotdata) {
  p <- ggplot(plotdata, aes(x = logFC, y = cwFC, color = genecolour, shape = geneshape)) +
    geom_point(size = 2, alpha = 0.8) +
    scale_color_manual(values = gg_def_pal_scmappr, guide = "none") +
    scale_shape_manual(values = c(2, 19)) +
    facet_wrap(.~celltype, scales = "free", ncol = 3) +
    labs(x = "Bulk logFC", y = "Scaled cwFC") +
    geom_text_repel(label = plotdata$genelab, force = 2, max.overlaps = 100) +
    theme_bw() +
    theme(
      text = element_text(size = 14))
  return(p)
}

get_top_cwFC <- function(use_data, cutoff) {
  sort_data <- as.data.frame(use_data)
  sort_data <- mutate(sort_data, genename = rownames(sort_data))
  colnames(sort_data) <- c("cwFC", "genename")
  sort_data <- filter(sort_data, abs(cwFC)>cutoff)
  return(sort_data)
}
make_cwfc_htmap <- function(plot_data){
  htmap_data <- bind_rows(lapply(names(htmap_data), function(x) mutate(htmap_data[[x]], celltype = x)))
  htmap_data$celltype <- factor(htmap_data$celltype)
  htmap_data <- dcast(htmap_data, genename~celltype, value.var = "cwFC")
  rownames(htmap_data) <- htmap_data$genename
  htmap_data <- htmap_data[,-1]
  htmap_data <- htmap_data[,which(colnames(htmap_data) %in% use_celltypes)]
  htmap_data[is.na(htmap_data)] <- 0
  return(htmap_data)
}

get_cwFC_genes <- function(use_data, use_celltype) {
  # print(use_celltype)
  sort_data <- as.data.frame(use_data)
  colnames(sort_data) <- "cwFC"
  sort_data <- mutate(sort_data, genename = rownames(sort_data), celltype=use_celltype)
  return(sort_data)
}

calc_cwfc_plots <- function(compar) {  
  cwFC <- get(load(paste0("scMappR_pituitary_d", compar, "_sex/scMappR_pituitary_d", compar, "_sex_cellWeighted_Foldchanges.Rdata")))
  props <- get(load(paste0("scMappR_pituitary_d", compar, "_sex/scMappR_pituitary_d", compar, "_sex_celltype_proportions.RData")))
  dds_significant <- bulk_DE[[paste0("d", compar, "_sex")]]
  
  # format bulk DE gene list to: c("gene_name", "padj", "log2fc")
  dds_significant <- dplyr::select(dds_significant, genename, FDR, logFC)
  # cwFC_eval <- cwFoldChange_evaluate_new(cwFC, props, dds_significant)
  cwFC_eval <- cwFoldChange_evaluate(cwFC, props, dds_significant)
  dir.create(paste0("scMappR_pituitary_d", compar, "_sex/cwFC_val_results"), showWarnings = F)
  write.table(cbind(genename = rownames(cwFC_eval$cwFoldchange_vs_bulk_rank_change),
                    cwFC_eval$cwFoldchange_vs_bulk_rank_change),
              paste0("scMappR_pituitary_d", compar, "_sex/cwFC_val_results/d",
                     compar, "bulk_vs_cwFC_rank_change.txt"),
              sep = "\t", quote = F, col.names = T, row.names = F)
  toOut <- lapply(cwFC_eval$cwFoldchange_gene_assigned, function(x) cbind(genename = names(x), as.data.frame(x)))
  write.xlsx(toOut,
             paste0("scMappR_pituitary_d", compar, "_sex/cwFC_val_results/d",
                    compar, "celltype_specific_DE_genes.xlsx"))
  
  cwFC_vs_bulk_rank <- apply(cwFC_eval$cwFoldchange_vs_bulk_rank_change, MARGIN = 2, function(x) list(rank_change=sort(abs(x-cwFC_eval$cwFoldchange_vs_bulk_rank_change$bulk),decreasing = T)))
  cwFC_vs_bulk_rank <- bind_rows(lapply(names(cwFC_vs_bulk_rank), function(x) cbind(celltype = x,genename = rownames(data.frame(cwFC_vs_bulk_rank[x][[1]])),
                                                                                    data.frame(cwFC_vs_bulk_rank[x][[1]])) %>%
                                          mutate(celltype_gene = paste0(celltype, "_", genename))))
  cwFC_vs_bulk_rank <- filter(cwFC_vs_bulk_rank, !celltype == "bulk")
  
  toOut <- toOut[which(lapply(toOut, nrow) > 0)]
  if(length(toOut) > 0) {
    cellspecific_DE <- bind_rows(lapply(names(toOut), function(x) arrange(cbind(celltype = x, toOut[[x]]), -abs(x)) %>%
                                          mutate(celltype_gene = paste0(celltype, "_", genename))))
    writeOut <- mutate(cellspecific_DE, direction = ifelse(x > 0, "F-biased", "M-biased")) %>%
      filter(abs(x) > 0.5) %>%
      dplyr::select(celltype, genename, direction)
    writeOut <- aggregate(writeOut$genename, by = list(writeOut$celltype, writeOut$direction), toString)
    colnames(writeOut) <- c("Cell-type", "DE_direction", "Genes")
    write.table(writeOut, paste0("scMappR_pituitary_d", compar, "_sex/cwFC_val_results/d",
                                 compar, "celltype_specific_DE_genenorm_cwFC_05_genes.txt"), quote = F,
                sep = "\t", row.names = F)
    combine_rank_cwDE <- full_join(dplyr::select(cellspecific_DE, x, celltype_gene, celltype),
                                   dplyr::select(cwFC_vs_bulk_rank, rank_change, celltype_gene),
                                   by = "celltype_gene") %>%
      dplyr::rename("cwFC" = x) %>%
      filter(!is.na(cwFC))
    
    # here the genes are first filtered for cwFC > 0.5
    # then the top 20 genes with the highest gene-normalized cwFC are labelled
    combine_rank_cwDE_label <- group_by(combine_rank_cwDE, celltype) %>%
      filter(abs(cwFC) > 0.5) %>%
      arrange(desc(rank_change), .by_group = T) %>%
      dplyr::slice(1:20, .by_group = T)
    
    cwFC_melt <- melt(cwFC_eval$cwFoldChange_normalized) %>%
      dplyr::rename("genename" = Var1, "celltype" = Var2, "cwFC" = value)
    cwFC_bulkFC <- full_join(cwFC_melt, dds_significant, by = "genename") %>% 
      filter(!is.na(celltype)) %>%
      mutate(celltype_gene = paste0(celltype, "_", genename)) %>%
      mutate(genelab = ifelse(celltype_gene %in% combine_rank_cwDE_label$celltype_gene, genename, "")) %>%
      mutate(genecolour = ifelse(celltype_gene %in% cellspecific_DE$celltype_gene, as.character(celltype), "not_cw_DE")) %>%
      mutate(geneshape = ifelse(genelab == "", F, T))
    
    # Plot all cell-types ##
    g <- cwFC_plot(cwFC_bulkFC)
    ggsave(paste0("scMappR_pituitary_d", compar, "_sex/cwFC_val_results/d",
                  compar, "_normcwFC_vs_bulkFC_plot.pdf"),
           g, width = 12, height = 12)
    
    # Only plot hormone-producing and precursor cell-types ##
    cwFC_bulkFC_sel <- filter(cwFC_bulkFC, celltype %in% c("Somatotropes", "Lactotropes", "Gonadotropes",
                                                           "Stem-cells_Sox2_FSC", "Proliferating_Pou1f1",
                                                           "Corticotropes", "Melanotropes_intermediate",
                                                           "Thyrotropes"))
    g <- cwFC_plot(cwFC_bulkFC_sel)
    ggsave(paste0("scMappR_pituitary_d", compar, "_sex/cwFC_val_results/d",
                  compar, "_normcwFC_vs_bulkFC_select_plot.pdf"),
           g, width = 12, height = 9)
    
    # Only plot sex-biased cell-types ##
    cwFC_bulkFC_sel2 <- filter(cwFC_bulkFC, celltype %in% c("Somatotropes", "Lactotropes", "Gonadotropes"))
    g <- cwFC_plot(cwFC_bulkFC_sel2)
    ggsave(paste0("scMappR_pituitary_d", compar, "_sex/cwFC_val_results/d",
                  compar, "_normcwFC_vs_bulkFC_select2_plot.pdf"),
           g, width = 10, height = 3.5)
    
    use_celltypes <- c("Somatotropes", "Lactotropes", "Gonadotropes",
                       "Corticotropes", "Melanotropes_intermediate",
                       "Thyrotropes")
    
    htmap_data <- cwFC_eval$cwFoldchange_gene_assigned[names(cwFC_eval$cwFoldchange_gene_assigned) %in% use_celltypes]
    htmap_data_sel <- lapply(htmap_data, function(x) get_top_cwFC(x, 0.5))
    htmap_data_sel <- bind_rows(lapply(names(htmap_data_sel), function(x) mutate(htmap_data_sel[[x]], celltype = x)))
    htmap_cwFC_plot <- cwFC_eval$cwFoldChange_normalized[which(rownames(cwFC_eval$cwFoldChange_normalized) %in% htmap_data_sel$genename),which(colnames(cwFC_eval$cwFoldChange_normalized) %in% use_celltypes), drop = F]
    
    # plot_cwfc <- apply(htmap_cwFC_plot, c(1,2), function(x) if(abs(x)>3){si=sign(x); x=si*3}else x = x)

    use_breaks <- get_htmap_breaks(htmap_cwFC_plot,
                                   colorRampPalette(c("orange", "white", "darkmagenta"))(255))
    
    if(nrow(htmap_cwFC_plot) > 1) {
      pcwfc <- pheatmap(data.frame(t(htmap_cwFC_plot)),
                        # scale = "column",
                        cluster_cols = T, cluster_rows = F, color = colorRampPalette(c("orange", "white", "darkmagenta") )(255),
                        breaks = use_breaks,
                        cellheight = 10, cellwidth = 8)
    }
    if(nrow(htmap_cwFC_plot) <= 1) {
      pcwfc <- pheatmap(data.frame(t(htmap_cwFC_plot)),
                        # scale = "column",
                        cluster_cols = F, cluster_rows = F, color = colorRampPalette(c("orange", "white", "darkmagenta") )(255),
                        breaks = use_breaks,
                        cellheight = 10, cellwidth = 8)
    }


    # pzscore <- pheatmap(t(htmap_cwFC_plot),
    #                  scale = "column",
    #                  cluster_cols = T, cluster_rows = F, color = colorRampPalette(c("orange", "white", "darkmagenta") )(255),
    #                  # breaks = use_breaks,
    #                  cellheight = 10, cellwidth = 8)
    write.table(mutate(as.data.frame(htmap_cwFC_plot), gene = rownames(htmap_cwFC_plot)), paste0("scMappR_pituitary_d", compar, "_sex/cwFC_val_results/d",
                                        compar, "_cwFC_table.txt"), quote = F, col.names = T, row.names = F, sep = "\t")
    
    sel_cwFC_genes <- bind_rows(lapply(names(cwFC_eval$cwFoldchange_gene_assigned[names(cwFC_eval$cwFoldchange_gene_assigned) %in% use_celltypes]),
                                       function(x) get_cwFC_genes(cwFC_eval$cwFoldchange_gene_assigned[[x]], x)))
    sel_cwFC_genes <- mutate(sel_cwFC_genes, cellgene = paste0(celltype, "_", genename))
    # View(sel_cwFC_genes)
    ortho <- getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol",
                    values = sel_cwFC_genes$genename, mart = mouse,
                    attributesL = c("rgd_symbol"), martL = rat,)
    colnames(ortho) <- c("genename", "rgd_symbol")
    ortho <- full_join(ortho, sel_cwFC_genes, by = "genename")
    colnames(ortho) <- c("mgi_symbol", "genename", "cwFC", "celltype", "cellgene")
    cwFC_fletcher <- full_join(ortho, fletcher_DE, by = "genename") %>%
      mutate(cwFC_sex_bias = ifelse(cwFC < 0, "M", "F")) %>%
      mutate(cellmatch = ifelse(celltype == cell_type & cwFC_sex_bias == sex_bias, "Y", "N")) %>%
      arrange(desc(cellmatch))
    
    cwFC_pit_pub <- full_join(sel_cwFC_genes, pit_pub_genes, by = "genename") %>%
      mutate(cwFC_sex_bias = ifelse(cwFC < 0, "M", "F")) %>%
      mutate(sort_col = ifelse((is.na(source_group) | is.na(cwFC_sex_bias)), NA, paste0(cwFC_sex_bias, "_", source_group))) %>%
      arrange(sort_col)
    
    gtex_cwfc <- filter(gtex_genes, variable == paste0("d", compar)) %>%
      filter(!is.na(value))
    
    gtex_cwfc <- full_join(gtex_cwfc, sel_cwFC_genes, by = "genename") %>%
      mutate(dir_match = ifelse(sign(MASH.beta) == sign(cwFC), "Y", NA)) %>%
      arrange(dir_match)
    
    return(list(heatmap_cwfc = pcwfc, fletcher = cwFC_fletcher, pit_pub = cwFC_pit_pub, gtex = gtex_cwfc,
                table_cwfc = htmap_cwFC_plot))
  }
  else{print(paste0("d", compar, ": No cell-type specific DE genes"))}
}

save_pheatmap_pdf <- function(x, filename, width, height) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width, height)
  for(i in 1:length(x)) {
    grid::grid.newpage()
    grid::grid.draw(x[[i]]$gtable)
  }
  tmp <- dev.off()
}

utrplot_hhtheme <- function(genelist, count_data, metadata,
                            plotrows = NULL,
                            type = "log2", pal_cols = c("F" = "tomato", "M" = "steelblue")) {
  use_genelist <- genelist[which(genelist %in% rownames(count_data))]
  if(length(genelist[!(genelist %in% rownames(count_data))]) > 0) {
    print(paste0("Gene not found in count data: ", genelist[!(genelist %in% rownames(count_data))]))
  }
  genecounts <- as.data.frame((count_data[which(rownames(count_data) %in% use_genelist),]))
  colnames(genecounts) <- use_genelist
  genecounts$sample <- rownames(genecounts)
  
  genecounts <- melt(genecounts)
  
  genecounts <- cbind(genecounts, "sex" = as.character(metadata$sex), "time" = gsub("d", "", metadata$age))
  
  genecounts$time <- as.numeric(as.character(genecounts$time))
  
  genecounts$reps <- paste0(metadata$sex, metadata$rep)
  
  gene_med <- aggregate(genecounts$value, by=list(genecounts$sex, genecounts$time, genecounts$variable), median)
  
  colnames(gene_med) <- c("sex", "time", "variable", "median")
  
  p <- ggplot() +
    geom_point(data = gene_med, aes(x = time, y = median, fill = sex, shape = sex, group = sex), color = "black", size = 4, alpha = 0.8) +
    geom_jitter(data = genecounts, aes(x = time, y = value, shape = sex, color = sex), width=0.2, height = 0) + 
    geom_line(data=gene_med, aes(x = time, y = median, color = sex, group = sex))  +
    xlab("Age (postnatal days)") + ylab("log2(normCounts)") +
    scale_fill_manual(values = pal_cols, name = "Sex") +
    scale_color_manual(values = pal_cols, name = "Sex") +
    scale_shape_manual(values = c(21, 24)) + 
    theme_bw() +
    facet_wrap( .~variable, scales = "free", nrow = plotrows) +
    theme(strip.text = element_text(size=18),
          axis.text = element_text(size=18 - 2, color="black"),
          axis.title = element_text(size=18),
          legend.position = "right", text = element_text(size = 18))+
    # strip.background = element_rect(colour = "red", fill = alpha("blue",0.2) )) +
    scale_x_continuous(breaks = c(12, 22, 27, 32, 37)) +
    guides(shape = F)
  return(p)
}

get_htmap_breaks <- function(plotdata, use_pal) {
  paletteLength <- length(use_pal)
  myBreaks <- c(seq(min(plotdata), 0, length.out=ceiling(paletteLength/2) + 1), 
                seq(max(plotdata)/paletteLength, max(plotdata), length.out=floor(paletteLength/2)))
  custome_breaks <- unique(c(min(myBreaks), seq(quantile(myBreaks, 0.2), 0, length = paletteLength/2), seq(0, quantile(myBreaks, 0.8), length=paletteLength/2), max(myBreaks)))
  return(custome_breaks)
}

# get_htmap_breaks <- function(plotdata){
#   m <- apply(plotdata, 1, mean, na.rm = T)
#   s <- apply(plotdata, 1, sd, na.rm = T)
#   use_data <- ((plotdata - m) / s)
#   all_values <- as.vector(as.matrix(use_data))
#   custome_breaks <- unique(c(min(all_values), seq(quantile(all_values, 0.01), 0, length = 125), seq(0, quantile(all_values, 0.99), length=125), max(all_values)))
#   # View(custome_breaks)
#   return(custome_breaks)
# }


abs_log <- function(x) {
  si <- sign(x)
  return(si*(log2(abs(x)+1)))
}
