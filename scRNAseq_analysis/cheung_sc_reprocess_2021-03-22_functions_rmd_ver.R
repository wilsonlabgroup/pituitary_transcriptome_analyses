#### Functions ####
calc_stat_test <- function(genelist, stat_test, seuratobj, use_clusters) {
  df <- df_FC <- as.data.frame(matrix(NA, nrow = length(genelist), ncol = length(levels(use_clusters)), dimnames = list(genelist, levels(use_clusters))))
  for(g in genelist) {
    tryCatch({
      scalex <- seuratobj@assays$RNA@scale.data[g, ]
      for(c in levels(use_clusters)) {
        cind <- which(names(use_clusters) %in% names(use_clusters[use_clusters == c]))
        c1 <- scalex[cind]
        c2 <- scalex[-cind]
        # c_plot <- cbind.fill(c1, c2, fill = NA)
        # colnames(c_plot) <- c("c1", "c2")
        # c_melt <- na.exclude(melt(c_plot))
        # p <- ggplot(c_melt) + geom_density(aes(x=value, fill = variable, color = variable), alpha = 0.7) +
        #   scale_fill_manual(values = c("darkmagenta", "gray70")) +
        #   scale_color_manual(values = c("darkmagenta", "gray70")) +
        #   labs(x = "Scaled expression value")
        #   theme_bw() +
        #   theme(text = element_text(size = 14))
        # print(p)
        if(stat_test == "ks") {
          pval <- ks.test(c1, c2, alternative = "less")$p.value
        }
        if(stat_test == "wilcox") {
          pval <- wilcox.test(c1, c2, alternative = "less")$p.value
        }
        FC <- median(c1) - median(c2)
        df[g, c] <- pval
        df_FC[g, c] <- FC
      }
    } , error = function(e) return(NA))
  }
  return(list("pval" = df, "FC" = df_FC))
}

get_htmap_breaks <- function(plotdata){
  m <- apply(plotdata, 1, mean, na.rm = T)
  s <- apply(plotdata, 1, sd, na.rm = T)
  use_data <- ((plotdata - m) / s)
  all_values <- as.vector(as.matrix(use_data))
  custome_breaks <- unique(c(min(all_values), seq(quantile(all_values, 0.001), 0, length = 125), seq(0, quantile(all_values, 0.999), length=125), max(all_values)))
  return(custome_breaks)
}

save_phtmap_pdf <- function(obj, filename, width, height) {
  stopifnot(!missing(obj))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(obj$gtable)
  tmp <- dev.off()
}

# Function takes cell ids/cluster ids and mean cell counts for each group into a count matrix accordingly 
fill_seurat_mat_mean <- function(cellid, seuratobj) {
  # print(cellid)
  cellbarcodes <- Idents(seuratobj) == cellid
  cellexpr <- seuratobj[["RNA"]]@data[, cellbarcodes]
  return(rowMeans(as.matrix(cellexpr)))
}

# Plots cell proportions from CIBERSORT output
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
    scale_x_continuous(breaks = c(12, 22, 27, 32, 37)) +
    facet_grid(.~cell_type) +
    labs(x = "Age", y = "Cell Proportion"
         # title = cluster
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


