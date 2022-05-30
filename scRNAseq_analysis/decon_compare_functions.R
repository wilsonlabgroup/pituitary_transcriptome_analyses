# Outputs mean counts for each gene in each cell type
fill_seurat_mat_mean <- function(cellid, seuratobj) {
  # print(cellid)
  cellbarcodes <- Idents(seuratobj) == cellid
  cellexpr <- seuratobj[["SCT"]]@data[, cellbarcodes]
  return(rowMeans(as.matrix(cellexpr)))
}

fill_seurat_mat_mean_test <- function(cellid, countdf, seuratobj) {
  # print(cellid)
  cellbarcodes <- Idents(seuratobj) == cellid
  cellexpr <- countdf[, cellbarcodes]
  return(rowMeans(as.matrix(cellexpr)))
}

# Reformats cell proportion dataframes for plotting
props_melt <- function(props_df) {
  props_df <- mutate(props_df, sample = rownames(props_df))
  props_df <- tbl_df(melt(props_df))
  colnames(props_df) <- c("sample", "cell_type", "prop")
  props_df <- inner_join(dplyr::select(meta, sample, age, sex, rep), props_df, by = "sample")
  props_df$age <- as.numeric(props_df$age)
  return(props_df)
}

# Plots cell proportions from deconvolution output
get_cell_prop_plot <- function(cell_id, cell_prop_data, palette) {
  print(cell_id)
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
    scale_x_continuous(breaks = c(12,22,27,32,37)) +
    facet_grid(.~cell_type) +
    labs(x = "Age", y = "Cell Proportion"
         # title = decon_method
         # title = paste0("Cluster: ", substr(i, 2, 4)), subtitle = new_cluster_ids[index] # Use with unbiased clustering
    ) +
    theme_bw() +
    theme(
      strip.background.x = element_rect(fill = "grey85"),
      strip.text = element_text(size=14, face="bold.italic"),
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

#### Functions for decon_comparison_correlation.R ####
plot_corr <- function(decon, scrna_input, age, scrna_props_F, scrna_props_M) {
  celltypes <- as.character(unique(unlist(lapply(decon$estimated, function(x) unique(colnames(x))))))
  
  # Check for missing cell types in certain deconvolution methods
  # If missing, give that cell type proportion values of 0 across all samples
  use_props <- lapply(decon$estimated, function(x) check_celltypes(x, celltypes))
  
  if(age == 37) {
    pd37_props_F <- lapply(use_props, function(x) melt(x[grep("Pd37F", rownames(x)),]) %>%
                             dplyr::rename("celltype" = variable))
    pd37_props_M <- lapply(use_props, function(x) melt(x[grep("Pd37M", rownames(x)),]) %>%
                             dplyr::rename("celltype" = variable))
  }
  
  if(age == 22) {
    pd37_props_F <- lapply(use_props, function(x) melt(x[grep("Pd22F", rownames(x)),]) %>%
                             dplyr::rename("celltype" = variable))
    pd37_props_M <- lapply(use_props, function(x) melt(x[grep("Pd22M", rownames(x)),]) %>%
                             dplyr::rename("celltype" = variable))
  }
  
  num_rep_F <- nrow(filter(pd37_props_F$MuSiC, celltype == "Somatotropes"))
  num_rep_M <- nrow(filter(pd37_props_M$MuSiC, celltype == "Somatotropes"))
  
  ruf_scrna_F <- scrna_props_F
  ruf_scrna_M <- scrna_props_M
  
  ruf_scrna_F_sub <- rep(ruf_scrna_F[names(ruf_scrna_F) %in% celltypes], each = num_rep_F)
  ruf_scrna_F_sub <- as.data.frame(cbind(names(ruf_scrna_F_sub), ruf_scrna_F_sub))
  ruf_scrna_M_sub <- rep(ruf_scrna_M[names(ruf_scrna_M) %in% celltypes], each = num_rep_M)
  ruf_scrna_M_sub <- as.data.frame(cbind(names(ruf_scrna_M_sub), ruf_scrna_M_sub))
  colnames(ruf_scrna_M_sub) <- colnames(ruf_scrna_F_sub) <- c("celltype", "value")
  
  F_cor <- lapply(pd37_props_F, function(decon_method) cor.test(decon_method$value, as.numeric(ruf_scrna_F_sub$value), method = "spearman"))
  M_cor <- lapply(pd37_props_M, function(decon_method) cor.test(decon_method$value, as.numeric(ruf_scrna_M_sub$value)))
  
  F_plot <- lapply(names(pd37_props_F), function(x) mutate(pd37_props_F[[x]], decon_method = x)) %>%
    bind_rows()
  F_plot <- rbind(unique(mutate(ruf_scrna_F_sub, decon_method = "scRNAseq_ref")),
                  F_plot)
  F_plot$value <- as.numeric(F_plot$value)
  celltype_named <- lapply(celltypes, function(x) x)
  names(celltype_named) <- celltypes
  F_plot$celltype <- factor(F_plot$celltype, levels = celltype_named)
  F_plot$decon_method <- factor(F_plot$decon_method, levels = list(Cibersort = "Cibersort",
                                                                   Cibersortx = "Cibersortx",
                                                                   CPM = "CPM",
                                                                   DCQ = "DCQ",
                                                                   DeconRNAseq = "DeconRNAseq",
                                                                   MuSiC = "MuSiC",
                                                                   NNLS = "NNLS",
                                                                   WGCNA = "WGCNA",
                                                                   scRNAseq_ref = "scRNAseq_ref"))
  
  gg_def_pal <- c(hue_pal()(length(levels(F_plot$decon_method)) - 1), "gray30")
  names(gg_def_pal) <- levels(F_plot$decon_method)
  
  M_plot <- lapply(names(pd37_props_M), function(x) mutate(pd37_props_M[[x]], decon_method = x)) %>%
    bind_rows()
  M_plot <- rbind(unique(mutate(ruf_scrna_M_sub, decon_method = "scRNAseq_ref")),
                  M_plot)
  
  M_plot$value <- as.numeric(M_plot$value)
  M_plot$celltype <- factor(M_plot$celltype, levels = celltype_named)
  M_plot$decon_method <- factor(M_plot$decon_method, levels = list(Cibersort = "Cibersort",
                                                                   Cibersortx = "Cibersortx",
                                                                   CPM = "CPM",
                                                                   DCQ = "DCQ",
                                                                   DeconRNAseq = "DeconRNAseq",
                                                                   MuSiC = "MuSiC",
                                                                   NNLS = "NNLS",
                                                                   WGCNA = "WGCNA",
                                                                   scRNAseq_ref = "scRNAseq_ref"))
  plot_list <- list(F_plot, M_plot)
  cor_list <- list(F_cor, M_cor)
  names(plot_list) <- names(cor_list) <- c(paste0("pd", age, "_F"), paste0("pd", age, "_M"))
  return(list(plots=lapply(names(plot_list), function(x) make_plots(x, plot_list[[x]],
                                                                    decon_method, gg_def_pal,
                                                                    scrna_input,
                                                                    cor_list[[x]])),
              cor=cor_list))
}

make_plots <- function(sex, F_plot, decon_method,
                       gg_def_pal, scrna_input,
                       F_cor) {
  g1 <- ggplot(F_plot, aes(x = decon_method, y = value,  color = decon_method)) +
    geom_boxplot() +
    geom_jitter() +
    labs(x = "Estimated cell proportions", y = "Deconvolution method",
         title = scrna_input,
         subtitle = sex) +
    scale_color_manual(values = gg_def_pal) +
    facet_wrap(.~celltype) +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          text = element_text(size = 14))
  g1a <- ggplot(F_plot, aes(x = decon_method, y = value,  color = decon_method)) +
    geom_boxplot() +
    geom_jitter() +
    labs(x = "Estimated cell proportions", y = "Deconvolution method",
         title = scrna_input,
         subtitle = sex) +
    scale_color_manual(values = gg_def_pal) +
    facet_wrap(.~celltype, scales = "free_y") +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          text = element_text(size = 14))
  F_cor_plot <- melt(t(as.data.frame(lapply(F_cor, function(x) x$estimate))))
  F_cor_plot$Var1 <- factor(F_cor_plot$Var1, levels = list(Cibersort = "Cibersort",
                                                           Cibersortx = "Cibersortx",
                                                           CPM = "CPM",
                                                           DCQ = "DCQ",
                                                           DeconRNAseq = "DeconRNAseq",
                                                           MuSiC = "MuSiC",
                                                           NNLS = "NNLS",
                                                           WGCNA = "WGCNA"))
  g2 <- ggplot(F_cor_plot, aes(x = Var1, y = value, color = Var1)) +
    geom_point(size = 3) +
    scale_color_manual(values = gg_def_pal[-length(gg_def_pal)]) +
    labs(x = "Deconvolution method", y = "Pearson corr") +
    theme_bw() +
    theme(legend.position = "none",
          text = element_text(size = 14))
  
  g1_comb <- ggarrange(plotlist=list(g1, g2),
                       nrow = 2,
                       heights = c(3,1))
  g1a_comb <- ggarrange(plotlist=list(g1a, g2),
                        nrow = 2,
                        heights = c(3,1))
  
  ggsave(paste0("output/corr_plots/", scrna_input, "_", sex, "_proportions_correlation.pdf"),
         g1_comb, width = 10, height = 10,
         useDingbats = F)
  ggsave(paste0("output/corr_plots/", scrna_input, "_", sex, "_scaled_proportions_correlation.pdf"),
         g1a_comb, width = 10, height = 10,
         useDingbats = F)
  return(list(g1, g1a_comb))
}

check_celltypes <- function(df, celltypes) {
  missing <- ncol(df) != length(celltypes)
  if(missing) {
    missing_celltype <- celltypes[-which(celltypes %in% colnames(df))]
    for(i in missing_celltype) {
      df$new <- 0
      colnames(df)[ncol(df)] <- i
    }
    return(df)
  } else { return(df) }
}

combine_prop_plots <- function(scidata, rufdata, decon_method) {
  scidata_sub <- scidata[, colnames(scidata) %in% colnames(rufdata)]
  rufdata_sub <- rufdata[, colnames(rufdata) %in% colnames(scidata)]
  
  scidata_sub <- melt(mutate(scidata_sub, sample=rownames(scidata_sub))) %>%
    mutate(scrna = "sci")
  rufdata_sub <- melt(mutate(rufdata_sub, sample=rownames(rufdata_sub))) %>%
    mutate(scrna = "ruf")
  
  combined <- rbind(scidata_sub, rufdata_sub) %>%
    full_join(., metadata, by = "sample") %>%
    mutate(age = as.numeric(gsub("d", "", age)))
  combined$scrna <- factor(combined$scrna, levels = list(sci = "sci",
                                                         ruf = "ruf"))
  combined_mean <- group_by(combined, group, variable) %>%
    summarise(mean=mean(value))
  combined_mean <- mutate(combined_mean, age = as.numeric(substr(group, 2, 3)), sex = substr(group,4, 4))
  combined_mean <- unique(full_join(combined_mean, dplyr::select(combined, group, scrna)))
  combined_mean$scrna <- factor(combined_mean$scrna, levels = list(sci = "sci",
                                                                   ruf = "ruf"))
  
  gg_def_pal <- c(hue_pal()(length(levels(combined$variable))))
  names(gg_def_pal) <- levels(combined$variable)
  
  # g <- ggplot(combined, aes(x = age, y = value, color = variable, fill = variable, group = group, shape = sex)) +
  #   # geom_boxplot() +
  #   geom_point(size = 2, alpha = 0.5) +
  #   scale_shape_manual(values = c(21,24)) +
  #   labs(x = "Postnatal age (days)", y = "Estimated cell proportions",
  #        title = decon_method) +
  #   # geom_jitter() +
  #   facet_wrap(variable~scrna, scales = "free", ncol = 6) +
  #   scale_x_continuous(breaks = c(12,22,27,32,37)) +
  #   theme_bw() +
  #   theme(text = element_text(size = 14))
  # g <- g + geom_line(data = combined_mean, aes(x = age, y = mean, group = sex, linetype = sex), size = 1, alpha = 0.8)
  # g <- g + geom_point(data = combined_mean, aes(x = age, y = mean, group = sex, shape = sex), size = 3, alpha = 0.8)
  g_all <- ggarrange(plotlist = lapply(levels(combined$variable), function(x) combined_plot_eachcell(x, filter(combined, variable == x),
                                                                                                     filter(combined_mean, variable == x),
                                                                                                     gg_def_pal[x])),
                     align = "hv") + ggtitle(decon_method)
  difference <- cbind(filter(combined_mean, sex == "F"), filter(combined_mean, sex == "M")) %>%
    mutate(diff = `mean...3`-`mean...9`) %>%
    dplyr::select(`variable...2`, `age...4`, diff)
  colnames(difference) <- c("celltype", "age", "diff")
  difference$age <- as.character(difference$age)
  pvals_list <- lapply(levels(combined$variable), function(x) calc_pvals(x, filter(combined, variable == x)))
  names(pvals_list) <- levels(combined$variable)
  pvals_df <- melt(bind_rows(unlist(pvals_list))) %>%
    mutate(celltype = unlist(lapply(as.character(variable), function(x) strsplit(x, ".", fixed = T)[[1]][1])),
           age = unlist(lapply(as.character(variable), function(x) strsplit(x, ".", fixed = T)[[1]][2])))
  difference <- full_join(difference, pvals_df, by = c("celltype", "age"))
  colnames(difference) <- c("celltype", "age", "diff", "variable", "p.val")
  celltype_factor <- levels(combined$variable)
  names(celltype_factor) <- levels(combined$variable)
  difference$celltype <- factor(difference$celltype, levels = celltype_factor)
  g1 <- ggplot(difference, aes(x = age, y = diff, fill = celltype)) +
    geom_bar(stat = "identity") +
    labs(x = "Postnatal age (days)", y = "Difference in cell proportions") +
    facet_wrap(.~celltype) +
    geom_text(label = difference$p.val, nudge_y = 0.01, check_overlap = T) +
    ylim(-max(abs(difference$diff))-0.01, max(abs(difference$diff))+0.01) +
    geom_hline(yintercept = c(-0.05, 0.05), linetype = "dashed", color = "gray50") +
    theme_bw() +
    theme(text = element_text(size = 14))
  ggsave(paste0("output/sci_ruf_combined_proportions/", decon_method, "_combined_proportions.pdf"),
         g_all, width = 14, height = 8, useDingbats = F)
  ggsave(paste0("output/sci_ruf_combined_proportions/", decon_method, "_proportion_difference.pdf"),
         g1, width = 14, height = 8)
  
  return(list(combined = combined, 
              combined_mean = combined_mean,
              difference = difference,
              plots = list(prop = g_all, diff = g1)))
}

calc_pvals <- function(celltype, combined_sub) {
  pvals <- sapply(unique(combined_sub$age), function(x) round(wilcox.test(filter(combined_sub, sex == "F", age == x)$value,
                                                                          filter(combined_sub, sex == "M", age == x)$value)$p.value, 3))
  names(pvals) <- unique(combined_sub$age)
  return(pvals)
}

combined_plot_eachcell <- function(celltype, combined_sub,
                                   combined_mean_sub, pal) {
  pvals <- sapply(unique(combined_sub$age), function(x) wilcox.test(filter(combined_sub, sex == "F", age == x)$value,
                                                                    filter(combined_sub, sex == "M", age == x)$value)$p.value)
  # fc <- sapply(unique(combined_sub$age), function(x) (log2(mean(filter(combined_sub, sex == "F", age == x)$value))-log2(mean(filter(combined_sub, sex == "M", age == x)$value))))
  # names(fc) <- names(pvals) <- unique(combined_sub$age)
  # 
  g <- ggplot(combined_sub, aes(x = age, y = value, group = group, color = variable, fill = variable, shape = sex)) +
    # geom_boxplot() +
    geom_point(size = 1, alpha = 0.5) +
    scale_shape_manual(values = c(21,24)) +
    scale_fill_manual(values = pal) +
    scale_color_manual(values = pal) +
    labs(x = "Postnatal age (days)", y = "Estimated cell proportions",
         # title = decon_method,
         subtitle = celltype) +
    # geom_jitter() +
    facet_wrap(.~scrna, scales = "free_x", ncol = 2) +
    scale_x_continuous(breaks = c(12,22,27,32,37)) +
    theme_bw()
  g <- g + ggpubr::stat_compare_means(data = combined_sub, aes(group = sex),
                                      label = "p.format", method = "wilcox.test",
                                      hide.ns = F, size = 3,
                                      vjust = 2)
  g <- g + geom_line(data = combined_mean_sub, aes(x = age, y = mean, group = sex, linetype = sex), size = 1, alpha = 0.8)
  g <- g + geom_point(data = combined_mean_sub, aes(x = age, y = mean, group = sex, shape = sex), size = 4, alpha = 0.8)
  g <- g + theme(text = element_text(size = 14),
                 legend.position = "none")
  return(g)
}
