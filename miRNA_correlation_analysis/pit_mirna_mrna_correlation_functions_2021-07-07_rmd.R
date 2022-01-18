#### Palettes ####
age_comp_col <- c(`F` = "tomato", M = "steelblue")
age_compx_col <- c(`F` = "tomato", M = "steelblue", both = "orange")
age_mirna_col <- c(`F` = "tomato", M = "steelblue")
age_mirnax_col <- c(`F` = "tomato", M = "steelblue", both = "orange")
sex_comp_col <- c(F_bias = "tomato", M_bias = "steelblue")
sex_mirna_col <- c(F_bias = "tomato", M_bias = "steelblue")
pituitary_dis_col <- c(Y = lacroix_palette("PassionFruit")[1], N = "white")
pituitary_GWAS_col <- c(Y = lacroix_palette("PassionFruit")[3], N = "white")
puberty_GWAS_col <- c(Y = lacroix_palette("PassionFruit")[4], N = "white")
puberty_dis_col <- c(Y = lacroix_palette("PassionFruit")[6], N = "white")

genelist_col <- c(pituitary_disease = lacroix_palette("PassionFruit")[1],
                  pituitary_disease_GWAS = lacroix_palette("PassionFruit")[3],
                  puberty_GWAS = lacroix_palette("PassionFruit")[4],
                  puberty_disease = lacroix_palette("PassionFruit")[6],
                  N = "white")

htmap_anno_col <- list(age_comp = age_comp_col, age_compx = age_compx_col,
                       age_mirna = age_mirna_col, age_mirnax = age_mirnax_col,
                       pituitary_disease = pituitary_dis_col,
                       pituitary_disease_GWAS = pituitary_GWAS_col,
                       puberty_GWAS = puberty_GWAS_col,
                       puberty_disease = puberty_dis_col,
                       pub_pit_genes = genelist_col,
                       sex_comp = sex_comp_col,
                       sex_mirna = sex_mirna_col)
htmap_pal <- brewer.pal(n = 9, name = "BuPu")[3:9]


setup_interaction_df <- function(df) {
  df <- filter(df, grepl("_F", comparison) | grepl("_M", comparison)) %>%
    mutate(age_comp = ifelse(grepl("_F", comparison), "F", "M")) %>%
    mutate(age_compx = ifelse(grepl("_M", comparison) & age_comp == "F", "both", "NA")) %>%
    mutate(age_compx = ifelse(age_compx == "NA", age_comp, age_compx)) %>%
    filter(grepl("_F", comparison_mirna) | grepl("_M", comparison_mirna)) %>%
    mutate(age_mirna = ifelse(grepl("_F", comparison_mirna), "F", "M")) %>%
    mutate(age_mirnax = ifelse(grepl("_M", comparison_mirna) & age_mirna == "F", "both", "NA")) %>%
    mutate(age_mirnax = ifelse(age_mirnax == "NA", age_mirna, age_mirnax))
  return(df)
}

add_genelist_anno <- function(gene_anno, use_genes) {
  pub_pit_genes <- use_genes
  gene_anno <- left_join(gene_anno, pub_pit_genes, by = "gene")
  gene_anno <- mutate(gene_anno, pituitary_disease = ifelse(grepl("pituitary_disease", source_group), "Y", "N"),
                                   pituitary_disease_GWAS = ifelse(grepl("pituitary_disease_GWAS", source_group), "Y", "N"),
                                   puberty_GWAS = ifelse(grepl("puberty_GWAS", source_group), "Y", "N"),
                                   puberty_disease = ifelse(grepl("puberty_disease", source_group), "Y", "N"))
  gene_anno <- mutate(gene_anno, pub_pit_genes = ifelse(!is.na(source), as.character(source_group), "N"))
  rownames(gene_anno) <- gene_anno$gene
  return(gene_anno)
}
#### Interaction heatmap function ####
make_sex_corr_heatmap <- function(df, plot_title = "", out_file, useQuantileFilter = F, genelist, plot_width=8, plot_height=8) {
  df <- mutate(df, sex_comp = ifelse(grepl("_sex_up", comparison), "F_bias", "M_bias")) %>%
    mutate(sex_mirna = ifelse(grepl("_sex_up", comparison_mirna), "F_bias", "M_bias"))
  
  neg_DE_htmap_mirna_anno <- unique(dplyr::count(df, mirna) %>% 
                                      inner_join(., dplyr::select(df, mirna, comparison_mirna, sex_mirna), by = "mirna")) %>%
    arrange(desc(n))
  
  if(useQuantileFilter == T) {
    filter_mirna <- ceiling(quantile(neg_DE_htmap_mirna_anno$n)[2])
    neg_DE_htmap_mirna_anno <- filter(neg_DE_htmap_mirna_anno, n >= filter_mirna)
    df <- filter(df, mirna %in% neg_DE_htmap_mirna_anno$mirna)
    plot_title <- paste0(plot_title, "\nnTargets >= ", filter_mirna)
  }
  rownames(neg_DE_htmap_mirna_anno) <- neg_DE_htmap_mirna_anno$mirna
  neg_DE_htmap_gene_anno <- unique(dplyr::count(df, gene) %>%
                                     inner_join(., dplyr::select(df, gene, comparison, sex_comp), by = "gene"))
  rownames(neg_DE_htmap_gene_anno) <- neg_DE_htmap_gene_anno$gene
  plot_df <- dcast(df, gene~mirna, value.var = "rho")
  
  # Rearrange miRNA order in heatmap so that miRNA with top number of gene targets is first
  df <- merge(neg_DE_htmap_mirna_anno, df, by.all = "mirna", sort = F)
  neg_DE_htmap_gene_anno <- neg_DE_htmap_gene_anno[match(unique(df$gene), neg_DE_htmap_gene_anno$gene), ]
  rownames(neg_DE_htmap_gene_anno) <- neg_DE_htmap_gene_anno$gene
  plot_df <- plot_df[match(rownames(neg_DE_htmap_gene_anno), plot_df$gene), ]
  rownames(plot_df) <- plot_df$gene
  plot_df <- plot_df[, -1]
  
  # Add in genelist annotation
  neg_DE_htmap_gene_anno <- add_genelist_anno(neg_DE_htmap_gene_anno, genelist)
  
  # Rearrange gene order in heatmap so that genes are ordered by top miRNA targeting it
  plot_df <- t(plot_df[, match(rownames(neg_DE_htmap_mirna_anno), colnames(plot_df))])
  
  p <- pheatmap(plot_df,
           cluster_cols = F,
           cluster_rows = F,
           annotation_col = as.data.frame(neg_DE_htmap_gene_anno)[, c(4,7:10), drop = F],
           annotation_row = as.data.frame(neg_DE_htmap_mirna_anno)[, 4, drop = F], 
           color = rev(colorRampPalette(htmap_pal)(255)),
           annotation_colors = htmap_anno_col,
           na_col = "white", 
           cellwidth = 8,
           cellheight = 8,
           border_color = "gray90",
           main = plot_title,
           # width = 6, height = 11,
           # filename = paste0("output_files/", out_file, "_heatmap.pdf")
  )
  
  p_simple <- pheatmap(plot_df,
           cluster_cols = F,
           cluster_rows = F,
           annotation_col = as.data.frame(neg_DE_htmap_gene_anno)[, c(4,11), drop = F],
           annotation_row = as.data.frame(neg_DE_htmap_mirna_anno)[, 4, drop = F], 
           color = rev(colorRampPalette(htmap_pal)(255)),
           annotation_colors = htmap_anno_col,
           na_col = "white", 
           cellwidth = 8,
           cellheight = 8,
           border_color = "gray90",
           main = plot_title,
           # width = 6, height = 11,
           # filename = paste0("output_files/", out_file, "_simple_anno_heatmap.pdf")
  )
  save_pheatmap_pdf(list(p, p_simple), paste0("output_files/", out_file, "_heatmap.pdf"), plot_width, plot_height)
  return(list(p, p_simple))
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
#### Interaction network function ####
make_sex_interaction_network <- function(df, plot_title = "", out_file, useQuantileFilter = F, genelist) {
  df <- mutate(df, sex_comp = ifelse(grepl("_sex_up", comparison), "F_bias", "M_bias")) %>%
    mutate(sex_mirna = ifelse(grepl("_sex_up", comparison_mirna), "F_bias", "M_bias"))
  
  neg_DE_htmap_mirna_anno <- unique(dplyr::count(df, mirna) %>% 
                                      inner_join(., dplyr::select(df, mirna, comparison_mirna, sex_mirna), by = "mirna")) %>%
    arrange(desc(n))
  
  if(useQuantileFilter == T) {
    filter_mirna <- ceiling(quantile(neg_DE_htmap_mirna_anno$n)[2])
    neg_DE_htmap_mirna_anno <- filter(neg_DE_htmap_mirna_anno, n >= filter_mirna)
    df <- filter(df, mirna %in% neg_DE_htmap_mirna_anno$mirna)
    plot_title <- paste0(plot_title, "\nnTargets >= ", filter_mirna)
  }
  rownames(neg_DE_htmap_mirna_anno) <- neg_DE_htmap_mirna_anno$mirna
  neg_DE_htmap_gene_anno <- unique(dplyr::count(df, gene) %>%
                                     inner_join(., dplyr::select(df, gene, comparison, sex_comp), by = "gene"))
  rownames(neg_DE_htmap_gene_anno) <- neg_DE_htmap_gene_anno$gene
  plot_df <- dcast(df, gene~mirna, value.var = "rho")
  
  # Rearrange miRNA order in heatmap so that miRNA with top number of gene targets is first
  df <- merge(neg_DE_htmap_mirna_anno, df, by.all = "mirna", sort = F)
  neg_DE_htmap_gene_anno <- neg_DE_htmap_gene_anno[match(unique(df$gene), neg_DE_htmap_gene_anno$gene), ]
  rownames(neg_DE_htmap_gene_anno) <- neg_DE_htmap_gene_anno$gene
  
  # Add in genelist annotation
  neg_DE_htmap_gene_anno <- add_genelist_anno(neg_DE_htmap_gene_anno, genelist)
  
  library(ggraph)
  nodes <- unique(dplyr::select(df, gene, comparison, sex_comp) %>%
                    mutate(nodetype = "gene"))
  nodes <- inner_join(nodes, dplyr::select(neg_DE_htmap_gene_anno, gene, source_group), by = "gene") %>%
    mutate(source_group = ifelse(is.na(source_group), "N", as.character(source_group))) %>%
    dplyr::rename("id" = gene)
  
  nodes <- unique(rbind(nodes, dplyr::select(df, mirna, comparison_mirna, sex_mirna) %>%
                          dplyr::rename("comparison" = comparison_mirna, "sex_comp" = sex_mirna) %>%
                          mutate(nodetype = "mirna", source_group = "mirna") %>%
                          dplyr::rename("id" = mirna)))
  
  edges <- dplyr::select(df, mirna, gene, rho) %>%
    dplyr::rename("from" = mirna, "to" = gene)
  rho_lab <- unique(mutate(edges, label = "Y") %>%
                      dplyr::rename("id" = to) %>%
                      dplyr::select(id, label))
  
  rho_lab <- aggregate(rho_lab$label, by = list(rho_lab$id), toString)
  colnames(rho_lab) <- c("id", "label")
  rho_lab <- mutate(rho_lab, label = ifelse(label == "" | label == "Y", label, "Y"))
  # View(rho_lab)
  
  nodes <- unique(full_join(nodes, dplyr::select(rho_lab, id, label), by = "id")) %>%
    mutate(label = ifelse(is.na(label), id, label)) %>%
    mutate(label = ifelse(source_group != "N" & label == "", id, label)) %>%
    mutate(label = ifelse(label == "", label, id))
  library(igraph)
  net <- graph_from_data_frame(d = edges, vertices = nodes, directed = T)
  
  V(net)$connections <- degree(net)
  
  gene_sel <- degree(net, v = V(net), mode = c("all", "out", "in", "total"),
                     loops = TRUE, normalized = FALSE)
  gene_sel <- as.data.frame(gene_sel)
  # # gene_sel <- as.data.frame(gene_sel[order(gene_sel, decreasing = T)])
  gene_sel$gene <- rownames(gene_sel)
  colnames(gene_sel)[1] <- "degree"
  nodes <- merge(nodes, gene_sel, by.x = "id", by.y = "gene")
  
  net <- graph_from_data_frame(d = edges, vertices = nodes, directed = T)
  
  if(cytoscapePing() == "You are connected to Cytoscape!") {
    createNetworkFromIgraph(net, paste0(plot_title, "_igraph"))
  }
  return(list(nodes, edges))
}

#### Bargraph functions ####
calc_prop <- function(i) {
  # print(i)
  prop <- filter(plot_percent, dir_comp == i) %>%
    dplyr::select(num_genes)
  if(prop[1,1] != 0) {
    return(list(target = prop[1,1]/prop[2,1], all = 1-(prop[1,1]/prop[2,1])))
  } else{ return (list(target = 0, all = 1))}
}

#### gProfiler functions ####
count_intersection <- function(x){
  return(length(unlist(strsplit(x, ","))))
}
process_plotdata <- function(plotdata, exclude_domain="", testclusters=""){
  plotdata <- plotdata[order(plotdata$domain, plotdata$p.value),]
  plotdata$term.name <- factor(plotdata$term.name, levels=unique(as.character(plotdata$term.name)))
  #plotdata$domain <- factor(plotdata$domain, levels=c("BP", "MF", "CC", "keg", "hp"))
  plotdata <- subset(plotdata, !domain %in% exclude_domain)
  plotdata$ngenes <- paste0("n=",sapply(plotdata$intersection, count_intersection))
  plotdata$cluster <- factor(plotdata$cluster, levels=testclusters)
  
  textpos <- 0.90*(-log10(min(plotdata$p.value)))
  
  plotdata <- plotdata %>%
    dplyr::arrange(cluster, desc(domain), -p.value) %>%
    dplyr::mutate(order = row_number())
  return(plotdata)
}
plot_cluster_enrichment <- function(enrich_list, outname, width=8, text_angle=0, horizontal =F, topn=5, ol_size=1, exclude_domain="", cluster_order=NULL, dot=F, colors = "Set2"){
  testclusters <- names(enrich_list[sapply(enrich_list, nrow) >0])
  enrich_data_list <- lapply(testclusters, function(x){
    print(x)
    enrich_data <- enrich_list[[x]]
    enrich_data$domain <- factor(enrich_data$domain, levels = as.character(unique(enrich_data$domain)))
    enrich_data <- do.call("rbind", lapply(split(enrich_data, f=enrich_data$domain), function(x) {
      x <- subset(x, `overlap.size` >= ol_size)
      x <- na.omit(x[order(x$p.value),][1:min(nrow(x), topn),])
      return(x)
    }))
    enrich_data$term.name <- factor(enrich_data$term.name, levels=rev(unique(as.character(enrich_data$term.name))))
    enrich_data$cluster <- rep(x, nrow(enrich_data))
    return(enrich_data)
  })
  
  plotdata <- do.call("rbind", enrich_data_list)
  plotdata <- process_plotdata(plotdata, exclude_domain=exclude_domain, testclusters = testclusters)
  
  plotdata_all <- do.call("rbind", lapply(testclusters, function(x) {
    enrich_list[[x]]$cluster <- x
    return(enrich_list[[x]])}))
  plotdata_all <- process_plotdata(plotdata_all,exclude_domain=exclude_domain, testclusters = testclusters)
  
  
  if(!is.null(cluster_order)){
    plotdata$cluster <- factor(plotdata$cluster, levels=cluster_order)
    plotdata_all$cluster <- factor(plotdata_all$cluster, levels=cluster_order)
  }
  
  pdf(paste0(outname, "_gprofiler_enrichment_top", topn, "_", Sys.Date(), ".pdf"), width=width, height = nrow(plotdata)*0.15+1)
  textpos <- 0.90*(-log10(min(plotdata$p.value)))
  p <- ggplot(plotdata) +
    #geom_bar(aes(x=`term.name`, y=-log10(p.value), fill=domain), stat="identity") +
    geom_bar(aes(x=order, y=-log10(p.value), fill=domain), stat="identity") +
    geom_text(aes(x=order, y=textpos, label=ngenes)) +
    coord_flip()+
    theme_bw() +
    theme(axis.text = element_text(color="black", size=10)) +
    scale_x_continuous(
      breaks = plotdata$order,
      labels = plotdata$`term.name`,
      expand = c(0,0))
  p <- p +
    facet_grid(cluster~., scales = "free", space = "free") +
    theme(strip.text.y = element_text(angle=text_angle))
  if(horizontal){
    used_data <- subset(plotdata_all, `term.id` %in% plotdata$term.id)
    
    used_data <-  arrange(used_data, desc(domain),-p.value, cluster)
    idorder <- c(1:length(unique(used_data$term.name)))
    names(idorder) <- unique(used_data$term.name)
    used_data$order <- idorder[used_data$term.name]
    
    used_data$term.name <- factor(used_data$term.name, levels = unique(used_data$term.name))
    
    textpos <- 0.90*(-log10(min(used_data$p.value)))
    p <- ggplot(used_data) +
      #geom_bar(aes(x=`term.name`, y=-log10(p.value), fill=domain), stat="identity") +
      geom_bar(aes(x=`term.name`, y=-log10(p.value), fill=domain), stat="identity") +
      geom_text(aes(x=`term.name`, y=textpos, label=ngenes)) +
      coord_flip()+
      theme_bw() +
      theme(axis.text = element_text(color="black", size=10),
            axis.text.x = element_text(angle=text_angle)) +
      facet_grid(~cluster) +
      scale_fill_brewer(palette = colors)
    if(dot){
      p <- ggplot(used_data) +
        geom_point(aes(x=`term.name`, y=cluster, size=-log10(p.value), color=domain)) +
        coord_flip() +
        scale_color_brewer(palette = colors) +
        theme_bw() +
        theme(axis.text = element_text(color="black", size=10),
              axis.text.x = element_text(angle=text_angle)) 
    }
  }
  
  print(p)
  dev.off()
  return(p)
}

#### Pathway table functions ####
match_term <- function(term_list, use_df) {
  # print(unlist(term_list))
  match_df <- filter(use_df, term %in% unlist(term_list))
  topmirnas <- dplyr::count(match_df, mirna) %>% arrange(desc(n))
  # dplyr::count(match_df, gene) %>% arrange(desc(n))
  match_df <- unique(  inner_join(match_df, topmirnas, by = "mirna") %>%
                         arrange(desc(n)) %>%
                         dplyr::select(mirna, gene)) %>%
    mutate(gene = str_to_sentence(gene))
  
  match_df_aggr <- aggregate(match_df$gene, by = list(match_df$mirna), toString)
  colnames(match_df_aggr) <- c("mirna", "target_genes")
  match_df_aggr <- inner_join(match_df_aggr, topmirnas, by = "mirna") %>%
    arrange(desc(n)) %>%
    dplyr::select(mirna, target_genes)
  return(match_df_aggr)
}

genelist_match <- function(term_list, use_df) {
  match_df <- filter(use_df, term %in% unlist(term_list))
  match_df <- mutate(match_df, gene = str_to_sentence(gene))
  match_df <- unique(inner_join(match_df, pub_pit_genes, by = "gene") %>%
                       dplyr::select(mirna, gene, source, source_group))
  return(match_df)
}

#### miRNA and gene target expression plots ####
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


mirnaplot_hhtheme <- function(mirna, mirna_data, ylabel = "Expression level (log)") {
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
