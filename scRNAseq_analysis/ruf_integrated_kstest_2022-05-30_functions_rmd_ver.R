#### Functions ####
calc_stat_test <- function(genelist, stat_test, seuratobj, use_clusters) {
  df <- df_FC <- as.data.frame(matrix(NA, nrow = length(genelist), ncol = length(levels(use_clusters)), dimnames = list(genelist, levels(use_clusters))))
  for(g in genelist) {
    tryCatch({
      scalex <- ruf@assays$SCT@data[g, ]
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

reorg_ks <- function(index, df) {
  if(index == 1) {
    use_df <- df[1:indices[index],]
  } else{use_df <- df[(indices[index - 1]+1):indices[index],]} # +1 to skip that row (otherwise it will become duplicated)
  # row_rank <- order(rowMin(use_df))
  row_rank <- order(rowSums(use_df))
  use_df <- as.data.frame(use_df[row_rank,])
  # use_df$genes <- rownames(use_df)
  return(as.data.frame(use_df))
}

num_enrich <- function(module_num) {
  print(" ")
  for(i in 1:length(ks_org[[module_num]])) {
    celltype <- names(ks_org[[module_num]])[i]
    txt <- paste0("M", module_num, " in ", celltype, ": ",
                 length(which(ks_org[[module_num]][[i]] < 0.05)), "/",
                 length(ks_org[[module_num]][[i]]),
                 " of ", table(mod_utr$modules)[[module_num]], " genes")
    print(txt)
    return(txt)
  }
}

hypergeo_test <- function(module, celltype) {
  m <- nrow(ks_resfdr[[module]]) # Genes in given module
  n <- total_genes - m # Genes not in given module
  k <- sum(sapply(ks_resfdr, function(x) length(which(x[, celltype] < 0.05)))) # Genes with FDR < 0.05 in given cell type (across all modules)
  x <- c(0:m) # Vector for 0:m to run the test on all combinations
  
  # Use dhyper built-in function for hypergeometric density
  probabilities <- dhyper(x, m, n, k, log = F)
  # probabilities
  
  # Calculate the one-sided p-value for 64 or more genes both FDR < 0.05 and IN M7.
  num_test <- length(which(ks_resfdr[[module]][,celltype] < 0.05))# Genes both IN the given module and has FDR < 0.05 in the given cell type
  pvalue <- sum(probabilities[num_test:(m+1)])
  return(pvalue)
}

