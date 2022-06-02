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