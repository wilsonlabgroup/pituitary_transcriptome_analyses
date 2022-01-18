
count_obj <- readRDS("~/Dropbox (Wilson Lab)/Mike_Anna_Mark/active_manuscript/Pituitary paper/scripts/miRNA_analysis/output_files/RUV_corrected_pit_srna-seq_counts_combined.rds")

mirna_counts <- normCounts(count_obj)
mirna_log <- log2(mirna_counts + 1)

mirnaplot_hhtheme <- function(mirna, mirna_data = mirna_log, ylabel = "Expression level (log)") {
  
  if(mirna %in% rownames(mirna_data)) {
    mirnax <- mirna_data[which(rownames(mirna_data) == mirna),]
    # mirnax <- t(as.data.frame(mirnax))
    mirnax <- as.data.frame(mirnax)
    mirnax1 <- cbind(mirnax, "sex" = substr(rownames(mirnax), 5, 5), "time" = substr(rownames(mirnax), 3, 4))
    colnames(mirnax1)[1] <- "value"
    mirnax1$time <- as.numeric(as.character(mirnax1$time))
    mirnax1$mirna <- paste0(mirnax1$sex, substr(rownames(mirnax1), 7, 7))
    if(substr(mirna, 1, 3) == "mmu") {
      mirna_title <- substr(mirna, 5, 50)
    }
    else {
      mirna_title <- mirna
    }
    mirnax1_sum <- aggregate(mirnax1$value, by=list(mirnax1$sex, mirnax1$time), median)
    mirnax1$mmu <- mirna_title
    colnames(mirnax1_sum) <- c("sex", "time", "median")
    p <- ggplot() +
      geom_point(data = mirnax1_sum, aes(x = time, y = median, fill = sex, shape = sex, group = sex), color = "black",  size = 4, alpha = 0.8) +
      geom_jitter(data = mirnax1, aes(x = time, y = value, shape = sex, color = sex), width=0.2, height = 0) + 
      geom_line(data=mirnax1_sum, aes(x = time, y = median, color = sex, group = sex))  +
      xlab("Age (postnatal days)") + ylab(ylabel) +
      scale_fill_manual(values = c("F" = "tomato", "M" = "steelblue"), name = "Sex") +
      scale_color_manual(values = c("F"= "tomato", "M" ="steelblue"), name = "Sex") +
      scale_shape_manual(values = c(21, 24)) + 
      theme_bw() +
      facet_grid( .~mmu) +
      theme(strip.text = element_text(size=18),
            axis.text = element_text(size=18 - 2, color="black"),
            axis.title = element_text(size=18),
            legend.position = "right", text = element_text(size = 18))+
      # strip.background = element_rect(colour = "red", fill = alpha("blue",0.2) )) +
      scale_x_continuous(breaks = c(12, 22, 27, 32)) +
      scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
      guides(shape = F)
    return(p)
  }
  else{
    print("miRNA not found.")
  }
}


DE_number_counts <- function(comparison, DE_list) {
  df <- DE_list[[comparison]]
  nup <- nrow(df[df$logFC >= log2FC_val & df$FDR < pval_val, ])
  ndown <- nrow(df[df$logFC <= -log2FC_val & df$FDR < pval_val, ])
  return(c("up" = nup, "down" = ndown))
}

plot_pca <- function(counts, annotation_df, plot_title, colour_label = NULL, shape_label = NULL, anno_colors, label_anno = NULL) {
  p <- autoplot(prcomp(t(counts)), data = annotation_df,
                colour = colour_label,
                shape = shape_label,
                size = 7,
                alpha = 0.8) +
    scale_colour_manual(values = anno_colors) +
    ggtitle(plot_title) +
    theme_bw() +
    theme(text = element_text(size = 14), legend.position = "bottom",
          legend.box = "vertical", legend.title.align = 0)
  if(!is.null(label_anno)) {
    p <- p + geom_text_repel(aes(label = label_anno))
  }
  return(p)
}

get_htmap_breaks <- function(plotdata){
  m <- apply(plotdata, 1, mean, na.rm = T)
  s <- apply(plotdata, 1, sd, na.rm = T)
  use_data <- ((plotdata - m) / s)
  all_values <- as.vector(as.matrix(use_data))
  custome_breaks <- unique(c(min(all_values), seq(quantile(all_values, 0.005), 0, length = 125), seq(0, quantile(all_values, 0.995), length=125), max(all_values)))
  return(custome_breaks)
}


plot_sample_dist_mat <- function(counts, annotation_df, plot_title, anno_colors = NULL,
                                 colors = colorRampPalette(c("orange", "white", "darkmagenta") )(255)) {
  # sampleDists <- dist(t(counts))
  # sampleDistMatrix <- as.matrix(sampleDists)
  # rownames(sampleDistMatrix) <- colnames(counts)
  # colnames(sampleDistMatrix) <- NULL
  
  # colors <- colorRampPalette(c("orange", "white", "darkmagenta") )(255)
  p <- pheatmap(cor(counts,method="pearson"),
                clustering_distance_rows="correlation",
                clustering_distance_cols="correlation",
                col=colors,
                border_color = NA,
                main = plot_title,
                annotation_col = annotation_df,
                show_colnames = F,
                annotation_colors = anno_colors)
  return(p)
}

get_edger_result <- function(contrasts=pair_contrasts, name){
  qlf <- glmQLFTest(fit, contrast = contrasts[, name])
  topTg <- topTags(qlf, n=nrow(y_fil$counts))
  
  de <- as.data.frame(topTg[[1]])
  return(de)
}

plot_volc <- function(comparison, resultdf, num_label = 10) {
  #### Input parameters ####
  print(comparison)
  DE_res <- resultdf[[comparison]]
  log2_label <- paste0("log2FC > ", log2FC_val)
  pval_label <- paste0("FDR < ", pval_val)
  
  #### Add annotations for miRNAs ####
  input <- as.data.frame(DE_res) %>% mutate(label = "NS") %>%
    mutate(label = ifelse(logFC > log2FC_val & FDR < pval_val & label == "NS", "Up", label)) %>%
    mutate(label = ifelse(logFC < -log2FC_val & FDR < pval_val & label == "NS", "Down", label))
  
  #### Get list of DE mirnas (log2FC > 1 & padj < 0.05) ####
  DE_mirnas <- subset(input, input$label != "NS") %>% arrange(FDR, abs(logFC)) %>% mutate(genes = gsub("mmu-", "", genes)) %>%
    dplyr::slice(1:num_label)
  
  #### Volcano plot colours ####
  volc_colours <- c("darkmagenta", lacroix_palette("Orange", 6, type = "discrete")[1], "grey50")
  names(volc_colours) <- c("Up", "Down", "NS")
  
  #### Plot volcano plot ####
  input$label <- factor(input$label, levels = c("Up" = "Up", "Down" = "Down", "NS" = "NS"))
  volc <- ggplot(input, aes(x = logFC, y = -log10(FDR), color = label)) +
    ggtitle(paste0(comparison)) +
    labs(subtitle=paste0("Up = ", nrow(filter(input, label == "Up")), "; Down = ", nrow(filter(input, label == "Down")))) +
    ylab("-log10 FDR") + xlab("log2 FC") +
    geom_point(size = 2, alpha = 0.8) + 
    scale_color_manual(values = volc_colours, drop = F) +
    theme_bw() +
    theme(legend.title=element_blank(), legend.position = "bottom",
          text = element_text(size = 16), axis.text = element_text(size=16)) +
    geom_vline(xintercept = c(-log2FC_val, log2FC_val), alpha = 0.5, linetype = 2) +
    geom_hline(yintercept = -log10(pval_val), alpha = 0.5, linetype = 2) +
    xlim(-max(abs(input$logFC)), max(abs(input$logFC)))
  
  #### Label top DE miRNAs ####
  if (nrow(DE_mirnas) > 0) {
    volc <- volc + geom_text_repel(data=DE_mirnas, aes(label=DE_mirnas$genes, size = 5),
                                   box.padding = 0.4, min.segment.length = 0.4, show.legend = F,
                                   colour = "black")
  }
  return(volc)
}

annotate_direction <- function(comparison, df) {
  results <- df[[comparison]]
  type <- "age"
  if(grepl("F-M", comparison)) {
    type <- "sex"
  }
  if(type == "age") {
    results$direction <- factor(ifelse(results$logFC < 0, "younger-bias", "older-bias"))
  }
  if(type == "sex") {
    results$direction <- factor(ifelse(results$logFC < 0, "male-bias", "female-bias"))
  }
  return(results)
}


save_phtmap_pdf <- function(obj, filename, width, height) {
  stopifnot(!missing(obj))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(obj$gtable)
  dev.off()
}

fromListx <- function (input) 
{
  elements <- unique(unlist(input))
  data <- unlist(lapply(input, function(x) {
    x <- as.vector(match(elements, x))
  }))
  data[is.na(data)] <- as.integer(0)
  data[data != 0] <- as.integer(1)
  data <- data.frame(matrix(data, ncol = length(input), byrow = F))
  data <- data[which(rowSums(data) != 0), ]
  names(data) <- names(input)
  rownames(data) <- elements
  return(data)
}

myfun <- function(row, element){
  data <- sum(as.numeric(row[grepl(element, names(row))])) > 0
}

make_upset <- function(mirnalist, setorder = NULL) {
  mirna_up <- lapply(mirnalist, function(x) x[x$logFC > 0, "genes"])
  names(mirna_up) <- paste0(names(mirna_up), "_up")
  mirna_down <- lapply(mirnalist, function(x) x[x$logFC < 0, "genes"])
  names(mirna_down) <- paste0(names(mirna_down), "_down")
  mirnas <- c(mirna_up, mirna_down)
  mirnas <- fromListx(mirnas)
  x <- upset(mirnas, order.by = "freq", main.bar.color = "black",
             sets= setorder, keep.order = T,
             # queries = list(list(query = myfun, params = list("_*up"), color = "black", active = F),
             #                list(query = myfun, params = list("_*down"), color = "grey", active = T)),
             text.scale = 2)
  return(x)
}
save_pdf <- function(plot, path, width, height) {
  pdf(path, width, height)
  print(plot)
  dev.off()
}

DE_number_counts_hh <- function(sig_list_data, comparison){
  de_data <- sig_list_data[[comparison]]
  up <- as.numeric(table(de_data$logFC > 0)[2])
  down <- as.numeric(table(de_data$logFC > 0)[1])
  conditions <- unlist(strsplit(comparison, "_"))[1:2]
  out <- as.data.frame(rbind(c(conditions, up), c(rev(conditions), down)), stringsAsFactors=F)
  colnames(out) <- c("con1", "con2", "DE")
  return(out)
}

gen_de_number_hm <- function(sig_list_data, high_color = "steelblue", pu_genes = "", fisher=F){
  DE_count_table <- do.call(rbind, lapply(names(sig_list_data), function(x) DE_number_counts_hh(sig_list_data, x)))
  DE_count_table$DE[is.na(DE_count_table$DE)] <- 0
  DE_count_table$DE <- as.numeric(as.character(DE_count_table$DE))
  sig_list_data_pu <- lapply(sig_list_data, function(x) subset(x, genes %in% pu_genes))
  DE_count_table_pu <- do.call(rbind, lapply(names(sig_list_data_pu), function(x) DE_number_counts_hh(sig_list_data_pu, x)))
  DE_count_table$puberty <- as.numeric(as.character(DE_count_table_pu$DE))
  if(fisher){
    DE_count_table$fisherP <- apply(DE_count_table, 1, function(x) get_fisher(length(all_detected_genes), length(all_detected_genes[all_detected_genes %in% pu_genes]), as.numeric(x[3]), as.numeric(x[4])))
  }else{
    DE_count_table$fisherP <- 1
  }
  DE_count_table$Padjust <- p.adjust(DE_count_table$fisherP, method="BH")
  DE_count_table$sig <- ifelse(DE_count_table$Padjust < 0.1, "yes", "no")
  DE_count_table$sig[is.na(DE_count_table$sig)] <- "no"
  DE_p <- ggplot(DE_count_table) +
    geom_tile(aes(x=con2, y=con1, fill=DE), color="grey") +
    geom_text(aes(x=con2, y=con1, label=DE, color=sig)) +
    #geom_point(data=subset(DE_count_table, Padjust < 5), aes(x=con2, y=con1), shape=8) +
    scale_fill_gradient2(high=high_color, low="white", name="Number of DE miRNAs") +
    scale_color_manual(values=c("yes"="red", "no"="black")) +
    theme(axis.title = element_blank())
  return(list(DE_p, DE_count_table)) 
}


mirnafactor_plot <- function(mirna_list,
                             mirna_data,
                             ylabel,
                             input_meta,
                             de_input = NULL) {
  use_mirna_data <- melt(mirna_data[rownames(mirna_data) %in% mirna_list,])
  colnames(use_mirna_data) <- c("mirna", "sample", "value")
  use_mirna_data <- inner_join(use_mirna_data, input_meta, by = "sample") %>%
    mutate(factor = paste0(mirna, "_d", age))
  if(!is.null(de_input)){
    # find maximum y axis position
    label_loc <- use_mirna_data %>% 
      group_by(mirna) %>% 
      summarise(ymax = max(value), ymin = min(value)) %>% 
      mutate(yloc = ymin + 1.1*(ymax - ymin))
    # filter for age at which mirna is significant
    label_loc <- label_loc %>%
      right_join(., dplyr::select(de_input, comparison, factor, mirna), by = "mirna") %>%
      mutate(time = substr(comparison, 2,3))
    label_loc$time <- as.numeric(label_loc$time)
    
    # create fake group for plotting purposes
    # "group" needs to be present for plotting to work
    label_loc$group<- paste0(label_loc$comparison, "F")
  }
  mirna_med <- use_mirna_data %>%
    mutate(factor = paste0(factor, sex))
  
  mirna_med <- unique(mirna_med %>%
                        group_by(factor) %>%
                        summarise(med = median(value)) %>%
                        inner_join(., dplyr::select(mirna_med, mirna, sex, age, group, factor), by = "factor"))
  colnames(mirna_med) <- c("factor", "median", "mirna", "sex", "time", "group")
  mirna_med$time <- as.numeric(as.character(mirna_med$time))
  
  colnames(use_mirna_data) <- c("mirna", "sample", "value", "sex", "time", "rep", "batch", "group", "factor")
  use_mirna_data$time <- as.numeric(as.character(use_mirna_data$time))
  
  g <- ggplot() +
    geom_point(data = mirna_med, aes(x = time, y = median, fill = sex, shape = sex, group = sex), color = "black",  size = 4, alpha = 0.8) +
    geom_jitter(data = use_mirna_data, aes(x = time, y = value, shape = sex, color = sex), width=0.2, height = 0) + 
    geom_line(data=mirna_med, aes(x = time, y = median, color = sex, group = sex))  +
    xlab("Age (postnatal days)") + ylab(ylabel) +
    scale_fill_manual(values = c("F" = "tomato", "M" = "steelblue"), name = "Sex") +
    scale_color_manual(values = c("F"= "tomato", "M" ="steelblue"), name = "Sex") +
    scale_shape_manual(values = c(21, 24)) + 
    theme_bw() +
    facet_wrap( .~mirna, scales = "free_y") +
    theme(strip.text = element_text(size=18),
          axis.text = element_text(size=18 - 2, color="black"),
          axis.title = element_text(size=18),
          legend.position = "right", text = element_text(size = 18))+
    # strip.background = element_rect(colour = "red", fill = alpha("blue",0.2) )) +
    scale_x_continuous(breaks = c(12, 22, 27, 32)) +
    scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
    guides(shape = F)
  
  if(!is.null(de_input)) {
    g <- g + 
      geom_point(data = label_loc, aes(x = time, y = yloc), shape = 8) 
  }
  return(g)
}
