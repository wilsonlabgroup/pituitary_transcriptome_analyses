# Load libraries
used_libraries <- c("edgeR","rtracklayer","GenomicRanges","RColorBrewer","gplots","pcaMethods","tidyverse","ggpubr", "ggrepel","RUVSeq","reshape2","gridExtra","cowplot","ggdendro", "dendextend", "fgsea", "gProfileR", "tidygraph", "ggraph", "scatterpie", "patchwork","UpSetR", "CEMiTool")
lapply(used_libraries, require, character.only=TRUE,  quietly=TRUE)

# set heatmap colors
hmcol<- colorRampPalette(brewer.pal(9, 'GnBu'))(100)
hmcol2 <- colorRampPalette(brewer.pal(9, 'RdBu'))(250)
hmcol_yp <- colorRampPalette(c("orange", "white", "darkmagenta") )(250)
# load mouse breeding info
#mouse_breeding <- readRDS("/mnt/mdwilson/huayun/puberty/phase1_mouse_breeding_info.rds")

# load named vector of gene id and name match
genename_ori <- genename <- readRDS("datasets/mm10_gencode_vM21_geneName.rds")
genetype <- readRDS("datasets/mm10_gencode_vM21_geneType.rds")

protein_coding <- names(genetype[genetype == "protein_coding"])

sexchrgenes <- c("Xist","Ddx3y","Uty","Gm29650","Kdm5d", "Eif2s3y", "Gm27733")

# load gmt files
#mouse_GO <- readRDS("/mnt/mdwilson/huayun/reference/gmt/Mm.c5.all.v7.1.symbol.rds")
#mouse_CP <- readRDS("/mnt/mdwilson/huayun/reference/gmt/Mm.c2.cp.v7.1.symbol.rds")

#################################################

# load information about GWAS
#alltfs <- readRDS("/mnt/mdwilson//huayun/puberty/all_tfs_msymbol.rds")

#HHgenes <- scan("~/Data/wilsonlab/puberty/fluidigm/HHgenes_mouse", what="character")
alltfs <- scan("datasets/mouse_TFs_animalTFDB3.txt", what = "character")

disease_genes <- read.table("datasets/pituitary_puberty_genes_combined_2020-09-28.txt", header = T, as.is = T, stringsAsFactors = F)

#GWAS_genes_human_mouse_anno <- readRDS("~/Dropbox/wilson_lab/Huayun/puberty/datasets/GWAS_genes_human_mouse_anno.rds")

# do not run
# gene_list <- lapply(names(GWAS_genes_human), function(x) {
#   dat = GWAS_genes_human[[x]]
#   return(data.frame(gene=dat, type=rep(x, length(dat))))})
# 
# GWAS_genes_human_anno <- do.call("rbind", gene_list) %>%
#   group_by(gene) %>%
#   summarize(source=paste(unique(sort(type)), collapse=";")) %>%
#   as.data.frame()

#GWAS_genes_human_mouse_anno <- merge(GWAS_genes_human_mouse_match, GWAS_genes_human_anno, by.x="HGNC.symbol", by.y="gene", all.x=T, all.y=F)

# all_puberty_genes <- subset(GWAS_genes_human_mouse_anno, source !="Day2017_blood_eQTL")$MGI.symbol
# all_puberty_genes_types <- subset(GWAS_genes_human_mouse_anno, source !="Day2017_blood_eQTL")$source
# all_puberty_genes_types <- gsub(";Day2017_blood_eQTL", "", all_puberty_genes_types)
# names(all_puberty_genes_types) <- all_puberty_genes
# 
# puberty_genes_nearest <- subset(GWAS_genes_human_mouse_anno, !source %in% c("Day2017_blood_eQTL","Day2017_hiC") )$MGI.symbol

#################################################
# get ERCC amount 
#ERCCMix1and2 <- read.table("/mnt/mdwilson/huayun/ganno/ERCC_concentration.txt", header = T, stringsAsFactors = F, as.is = T, sep="\t")

# ERCC_volume <- 2 # ul
# ERCC_dilution <- 1000 # 1/1000
# RNA_volume <- 0.1 # ng
# 
# ERCC_con <- ERCCMix1and2[, c(2,4)]
# names(ERCC_con)[2] <- "Mix1Conc.Attomoles_ul"
# 
# #ERCC_con$amount <- log2((ERCC_con$Mix1Conc.Attomoles_ul*ERCC_volume/ERCC_dilution)/RNA_volume)
# ERCC_con$amount <- log2(ERCC_con$Mix1Conc.Attomoles_ul*ERCC_volume/ERCC_dilution)
# ERCC_con <- ERCC_con[, -2]
# 
# ERCC_con_v <- ERCC_con$amount
# names(ERCC_con_v) <- ERCC_con$ERCC.ID

#saveRDS(ERCC_con_v, "/mnt/mdwilson/huayun/ganno/ERCC_concentrations_2ul.rds")

##################################################a

#library(devtools)
#source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

transform_df_samplename <- function(samplenames_all, field=7, list=TRUE, batch=FALSE, return_tissue=FALSE, repdouble=FALSE){
  transformed_list <- lapply(samplenames_all, function(samplenames){
    name_elements <- lapply(as.character(samplenames), function(x) unlist(strsplit(x, "_")))
    samplenames <- sapply(name_elements, "[[", field)
    
    tissue <- sapply(samplenames, function(x) substr(unlist(strsplit(x, "-"))[1],1,1))
    age <- sapply(samplenames, function(x) substr(unlist(strsplit(x, "-"))[1],2,4))
    sex <- sapply(samplenames, function(x) substr(unlist(strsplit(x, "-"))[1],5,5))
    #rep <- sapply(samplenames, function(x) unlist(strsplit(x, "-"))[2])
    if(repdouble){
      rep <- sapply(samplenames, function(x) substr(x, 6, nchar(x)))
    } else {
      rep <- sapply(samplenames, function(x) substr(x, nchar(x), nchar(x)))
    }
    
    
    ID <- sapply(name_elements, "[[", 1)
    # if(ID %in% c("WL551","WL552","WL553","WL554")){
    #   batch <- "third"
    # } else if(grepl("WL5", ID)){
    #   batch <- "first"
    # } else {
    #   batch <- "second"
    # }
    
    if(list){
      return(list(age, sex, rep, ID))
    } 
    if(return_tissue) {
      return(data.frame(age=age, sex=sex, rep=rep, id=ID, tissue=tissue))
    } else{
      return(data.frame(age=age, sex=sex, rep=rep, id=ID))
    }
  })
  outdf <- do.call("rbind", transformed_list)
  if(batch){
    outdf$batch <- batches[as.character(outdf$id)]
  }
  outdf$group <- with(outdf, paste0(age, sex))
  return(outdf)
}


orinames <- paste0("Hd12F", 2:7)
newnames <- paste0("Hd12F", 1:6)

fixname <- function(x){
  x_sample <- strsplit(x, "_")[[1]][2]
  if (x_sample %in% orinames){
    n <- which(orinames == x_sample)
    new <- gsub(orinames[n], newnames[n], x)
    return(new)
  } else{
    return(x)
  }
}

get_count_table <- function(counts_file, fixname=F){
  # load count table from featureCount (one run), also process colnames to get a sample condition table
  count_table <- read.table(counts_file, sep="\t", header=T, as.is = T, stringsAsFactors = F)
  row.names(count_table) <- count_table[,1]
  gene_info <- count_table[,1:6]
  count_table <- count_table[,7:ncol(count_table)]
  # obtain sample information
  sampleFiles <- names(count_table)
  samplenames <- unname(sapply(sampleFiles, function(x) {
    temp <- unlist(strsplit(x, "\\."))
    return(temp[grepl("^WL", temp)])}))
  samplenames <- unname(sapply(samplenames, function(x) paste(unlist(strsplit(basename(x), "_"))[c(1,7)], collapse="_")))
  if(fixname){
    samplenames <- unname(sapply(samplenames, fixname))
  }
  # rename columns of the count table
  names(count_table) <- samplenames
  
  # build a sample condition table from sample names 
  sampleCondition <- transform_df_samplename(samplenames, field = 2, list = FALSE)
  sampleCondition$rep <- sapply(row.names(sampleCondition), function(x) substring(x, nchar(x)))
  #row.names(sampleCondition) <- colnames(count_table)
  
  return(list(count_table, sampleCondition, gene_info))
}

get_count_stats <- function(count_table, cpm, minreads=5, select_genes=row.names(count_table)){
  count_stats <- data.frame(totalReads = apply(count_table, 2, sum),
                            genesDetected_counts = apply(count_table[select_genes,], 2, function(x) length(x[x>=minreads])),
                            genesDetected_cpm = apply(cpm[select_genes,], 2, function(x) length(x[x>=1]))) 
  count_stats$sample <- row.names(count_stats)
  count_stats.m <- melt(count_stats, id.vars = c("totalReads", "sample"))
  count_stats.m <- cbind(count_stats.m, transform_df_samplename(count_stats.m$sample, field = 2, list = F, return_tissue = T))
  
  return(count_stats.m)
}




#########################################################################

merge_RNA_PCR <- function(logCPM, PCR_puberty = PCR_puberty, use_genename=F){
  RNAseq <- as.data.frame(logCPM)
 # 
  if(use_genename){
    RNAseq$Primer <- row.names(RNAseq)
  } else{
    RNAseq$Primer <- unname(genename[row.names(RNAseq)])
  }
 
  RNAseq_puberty <- subset(RNAseq, Primer %in% all_genes)
  
  RNAseq_puberty_m <- melt(RNAseq_puberty, id.vars = "Primer")
  RNAseq_puberty_m$Sample <- sapply(RNAseq_puberty_m$variable, function(x) unlist(strsplit(as.character(x), "_"))[2])
  exp_compare <- merge(PCR_puberty, RNAseq_puberty_m, by=c("Sample", "Primer"), all=T)
  return(exp_compare)
}

get_cor_plot <- function(sample, exp_compare, type="plot", method="spearman"){
  plot_data <- subset(exp_compare, Sample == sample)
  min_value <-  min(na.omit(plot_data)$value)
  # num primers detected with qPCR
  npcr <- as.numeric(table(!is.na(plot_data$Ct))["TRUE"])
  nrna <- as.numeric(table(na.omit(plot_data)$value > min_value)["TRUE"])
  
  cor_value <- round(with(plot_data, cor(Ct, value, method=method)),2)
  cor_value2 <- round(with(na.omit(subset(plot_data, value != min_value)), cor(Ct, value, method=method)),2)
  variable <- plot_data[1,6]
  n <- nrow(na.omit(subset(plot_data, value != min_value)))

  plot_range <- plot_data[!is.na(plot_data$value),]
  ypos <- floor(max(na.omit(plot_range$value)))
  xpos <- round(mean(range(na.omit(plot_range$Ct))))
  xmin <- floor(range(na.omit(plot_range$Ct))[[1]])
  xmax <- ceiling(range(na.omit(plot_range$Ct))[[2]])
  
  p <- ggplot(plot_data, aes(x=Ct, y=value)) +
    #geom_point(alpha=0.7) +
    geom_text(aes(label = Primer), alpha=0.8) +
    annotate(geom="text", label = paste(variable, " : ", method ,"cor =", cor_value2, "\n   n=", n), x=xpos, y=ypos, color="blue") +
    ylab("RNA seq expression levels (logCPM)") +
    xlab("qPRC expression levels (delta Ct)") +
    xlim(c(xmin, xmax)) +
    #ylim(c(-1, 14)) +
    theme_bw()
  if(type != "plot"){
    return(list(cor_value2, npcr, nrna))
  } else {
    return(p)
  }
}



get_pca_var <- function(res){
  all_var <- R2cum(res)
  var1 <- round(all_var[1], 3)
  var2 <- round(all_var[2] - all_var[1], 3)
  var3 <- round(all_var[3] - all_var[2], 3)
  return(c(var1, var2))
}


get_weight_plot_2d <- function(temp_plot, method=method){
  # function to plot weight after unwanted factor estimation
  # color by condition/mouse breeding info
  p1 <- ggplot(temp_plot)+
    geom_point(aes(x=W_1, y=W_2, color=sex, shape=age), size=5) +
    geom_text_repel(aes(x=W_1, y=W_2, label=rep)) +
    #scale_shape_manual(values=21:25) +
  #  scale_color_manual(values=c("white", "black")) +
    scale_color_manual(values = c("red","steelblue")) +
    guides(color = guide_legend(override.aes = list(fill=c("F"="red","M"="steelblue")))) +
    theme_bw()
  
  num_colors <- length(unique(temp_plot$pregnancy.cage.ID))
  fill_cols <- brewer.pal(num_colors, "Paired")
  names(fill_cols) <- levels(factor(temp_plot$pregnancy.cage.ID))
  
  p2 <- ggplot(temp_plot)+
    geom_point(aes(x=W_1, y=W_2, color=as.factor(`pregnancy.cage.ID`), shape=age), size=5) +
    geom_text_repel(aes(x=W_1, y=W_2, label=rep)) +
  #  scale_shape_manual(values=21:25) +
    scale_color_manual(values = brewer.pal(num_colors, "Paired"), name="breeding_cage") +
    #guides(fill = guide_legend(override.aes = list(fill=fill_cols))) +
    #scale_color_manual(values = c("red","steelblue")) +
    theme_bw()
  
  # p3 <- ggplot(temp_plot) +
  #   geom_point(aes(x= as.factor(`pregnancy.cage.ID`), y=1, fill=as.factor(`pregnancy.cage.ID`)), shape=21, color="black", size=5) +
  #   scale_fill_manual(values = brewer.pal(12, "Paired"), name="breeding_cage") +
  #   theme_bw()
  
  pdf(paste0("weight_adjust_", method, "_", Sys.Date(),".pdf"))
  print(p1)
  print(p2)
  #print(p3)
  dev.off()
}

# get pca plot data using pcaMethod
get_pca_plotdata <- function(cpm_table, top=500, npc=15, field=2){
  sds <- apply(cpm_table, 1, sd)
  
  variable_genes <- names(sort(-abs(sds)))[1:top]
  #variable_genes <- names(sds)
  variable_counts <- cpm_table[variable_genes, ]
  
  center_data <- prep(t(variable_counts), scale="none", center=TRUE)
  resPCA <- pca(center_data, method="svd", center=FALSE, nPcs=npc)
  
  plot_data <- scores(resPCA)
  if(!is.null(field)){
    plot_data <- cbind(plot_data, transform_df_samplename(row.names(plot_data), field = field, list = FALSE))
    plot_data$age <- factor(plot_data$age, levels=c("d12","d22","d27","d32","d37"))
  }
  
  return(list(plot_data, resPCA))
}

# generate and save plots
get_pca_plots <- function(logCPM, logCPMc, method="ruv", npc=15){
  # removing chry_genes
  chry_genes <- gene_info$Geneid[grepl("chrY", gene_info$Chr)]
  
  # do pca without chrY genes
  PCA_removeBatch <- get_pca_plotdata(logCPMc[!(row.names(logCPMc) %in% chrXYgenes),], npc=npc)
  PCA_original <- get_pca_plotdata(logCPM[!(row.names(logCPM) %in% chrXYgenes),],npc=npc)
  
  removeBatch_var <- get_pca_var(PCA_removeBatch[[2]])
  original_var <- get_pca_var(PCA_original[[2]])
  
  p_removeBatch <- ggplot(PCA_removeBatch[[1]], aes(x=PC1, y=PC2, shape = age)) +
    geom_point(aes(color=sex), size=6, alpha=0.7) +
    theme_bw() +
    geom_text_repel(aes(label=id)) +
    scale_color_manual(values = c("red", "blue"))+
    xlab(paste("PC1 (", removeBatch_var[1], ")")) +
    ylab(paste("PC2 (", removeBatch_var[2], ")")) +
   # scale_shape_manual(values = c(1, 2, 0, 3, 7)) +
    theme(legend.position="none")
  
  p_removeBatch_color_by_batch <- ggplot(PCA_removeBatch[[1]], aes(x=PC1, y=PC2, shape = age)) +
    geom_point(aes(color=batch), size=6, alpha=0.7) +
    theme_bw() +
    geom_text_repel(aes(label=id)) +
    scale_color_manual(values = c("black", "chocolate", "springgreen4", "brown2"))+
    theme(legend.position="none") +
    xlab(paste("PC1 (", removeBatch_var[1], ")")) +
    ylab(paste("PC2 (", removeBatch_var[2], ")")) 
  
  
  p_ori <- ggplot(PCA_original[[1]], aes(x=PC1, y=PC2, shape = age)) +
    geom_point(aes(color=sex), size=6, alpha=0.7) +
    theme_bw() +
    geom_text_repel(aes(label=id)) +
    scale_color_manual(values = c("red", "blue")) +
    theme(legend.position="none") +
    xlab(paste("PC1 (", original_var[1], ")")) +
    ylab(paste("PC2 (", original_var[2], ")")) 
  
  p_ori_color_by_batch <- ggplot(PCA_original[[1]], aes(x=PC1, y=PC2, shape = age)) +
    geom_point(aes(color=batch), size=6, alpha=0.7) +
    theme_bw() +
    geom_text_repel(aes(label=id)) +
    scale_color_manual(values = c("black", "chocolate", "springgreen4", "brown2"))+
    theme(legend.position="none") +
    xlab(paste("PC1 (", original_var[1], ")")) +
    ylab(paste("PC2 (", original_var[2], ")")) 
  
  
  pdf(paste0("PCA_CPM_original_", method, "_", Sys.Date(), ".pdf"))
  print(p_ori)
 # print(p_ori_color_by_batch)
  dev.off()
  
  pdf(paste0("PCA_CPM_removeBatch_", method, "_", Sys.Date(),".pdf"))
  print(p_removeBatch)
  #print(p_removeBatch_color_by_batch)
  dev.off()
  
  PCA_all_removebatch <- get_pca_plotdata(logCPMc[!(row.names(logCPMc) %in% chrXYgenes),], top=nrow(logCPMc[!(row.names(logCPMc) %in% chrXYgenes),]), npc=npc)
  all_removeBatch_var <- get_pca_var(PCA_all_removebatch[[2]])
  
  pca_all <- ggplot(PCA_all_removebatch[[1]], aes(x=PC1, y=PC2, shape = age)) +
    geom_point(aes(color=sex), size=6, alpha=0.7) +
    theme_bw() +
    geom_text_repel(aes(label=rep)) +
    scale_color_manual(values = c("red", "blue")) +
    #theme(legend.position="none") +
    xlab(paste("PC1 (", all_removeBatch_var[1], ")")) +
    ylab(paste("PC2 (", all_removeBatch_var[2], ")")) 
  pdf(paste0("PCA_CPM_removeBatch_all_genes_", method, "_", Sys.Date(),".pdf"))
  print(pca_all)
  dev.off()
}

factor_plot <- function(data, gene_name, name=FALSE, legend=TRUE, ncols=0, fix_scale=FALSE, labelSize=18, field=2, multi_factor=F, sig_df = NULL, convert_genename = T,ylab = "log2(normCounts)"){
  # sig_df is a data frame with 2 columns: genename and the age at which it is significant
  # example: 
  # Fshb d12
  # Fshb d22
  # Lhb d22
  if(ncols==0){
    ncols=length(gene_name)
  }
    if(name){
    gene_id <- gene_name
  } else{
    gene_id <- names(genename[genename %in% gene_name])
  }
  data <- as.data.frame(data)
  gene_data <- data[row.names(data) %in% gene_id,]
  if(convert_genename){
    row.names(gene_data) <- genename[row.names(gene_data)]
  }
  
  #gene_data <- data.frame(sample=names(gene_data), logCPM=as.numeric(as.character(gene_data)))
  gene_data <- melt(as.matrix(gene_data))
  names(gene_data) <- c("factor", "sample", "logCPM")
  gene_data <- cbind(gene_data, do.call("rbind",lapply(gene_data$sample, function(x) transform_df_samplename(x, field = field, list = FALSE))))
  gene_data$age <- as.numeric(sapply(gene_data$age, function(x) substr(x, 2, 3)))
  
  gene_data_sum <- gene_data %>%
    group_by(factor, age, sex) %>%
    #summarise(mean = mean(logCPM), sd=sd(logCPM)) %>%
    summarise(median = median(logCPM, na.rm = T), sd=sd(logCPM)) %>%
    as.data.frame()
  
  gene_data_sum$factor <- factor(gene_data_sum$factor, levels = unique(as.character(gene_name)))
  gene_data$factor <- factor(gene_data$factor, levels = unique(as.character(gene_name)))
  
  # process data frame with significant genes
  if(!is.null(sig_df)){
    # find maximum y axis position
    label_loc <- gene_data %>% 
      group_by(factor) %>% 
      summarise(ymax = max(logCPM, na.rm = T), ymin = min(logCPM, na.rm = T)) %>% 
      mutate(yloc = ymin + 1.1*(ymax - ymin))
    
    # filter for significant gene
    names(sig_df) <- c("factor", "age")
    
    sig_df_used <- sig_df %>% 
      rowwise() %>% 
      mutate(age = as.numeric(gsub("d","", age))) %>% 
      filter(factor %in% gene_name) %>% 
      left_join(label_loc)
  }
  
  
  p_factor <- ggplot()+
    geom_point(data=gene_data_sum, aes(x=age, y=median, fill=sex, shape=sex), size=4, alpha=0.8) +
    geom_jitter(data=gene_data, aes(x=age, y=logCPM, color=sex, shape=sex), width=0.2, height = 0) +
    geom_line(data=gene_data_sum, aes(x=age, y=median, color=sex, group=sex)) +
    scale_color_manual( values=c("F" = "tomato", "M" = "steelblue")) +
    scale_fill_manual(values=c("F" ="tomato","M" = "steelblue")) +
    scale_shape_manual(values=c("F" = 21,"M" = 24)) +
   # scale_x_continuous(breaks=seq(10,40, by=5)) +
    scale_x_continuous(breaks=sort(as.numeric(unique(gene_data$age)))) +
    xlab("Age (postnatal days)") +
    ylab(ylab) +
    theme_bw() +
    theme(strip.text = element_text(size=labelSize, face="bold"),
          axis.text = element_text(size=labelSize - 2, color="black"),
          axis.title = element_text(size=labelSize))
  
  if(!is.null(sig_df)){
    p_factor <- p_factor + 
      geom_point(data = sig_df_used, aes(x = age, y = yloc), shape = 8) 
  }
  
  if(fix_scale){
    p_factor <- p_factor + 
      facet_wrap(~factor, ncol=ncols)
  } else{
    p_factor <- p_factor +
    facet_wrap(~factor, scales = "free_y", ncol=ncols)
  }
  
  if(!legend){
    p_factor <- p_factor +
      theme(legend.position = "none")
  }
  
  if(multi_factor){
    factor_cols <- colorRampPalette(brewer.pal(8, "Dark2"))(length(gene_name))
    p_factor <- ggplot(data=gene_data_sum, aes(x=age, y=median))+
      geom_point(aes(shape=sex, color=factor), size=4, alpha=0.8)+
      geom_line(aes(group=interaction(factor, sex), color=factor, linetype=sex)) +
      scale_x_continuous(breaks=sort(as.numeric(unique(gene_data$age)))) +
      xlab("Age (postnatal days)") +
      ylab(ylab) +
      theme_bw()+
      scale_shape_manual(values=c("F" = 16,"M" = 17)) +
      scale_color_manual(values = factor_cols)+
      scale_fill_brewer(palette = "Paired")+
      theme(strip.text = element_text(size=labelSize, face="bold"),
            axis.text = element_text(size=labelSize - 2, color="black"),
            axis.title = element_text(size=labelSize))
  }
  
  p_factor <- p_factor +  theme(strip.text = element_text(face="bold.italic"))
  return(p_factor)
}




library(RColorBrewer)
hmcol<- colorRampPalette(brewer.pal(9, 'Blues'))(100)

plotHM <- function(counts, method="pearson", cols=colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100), anno=sampleCondition_fil, show_col_names=T, show_row_names=T,new_col_names=colnames(counts), clust_method="complete"){
  mat<- as.matrix(cor(counts, method=method))
  
 # hc <- hclust(distsRL)
  anno_df <- anno[, c("age","sex")]
  rownames(anno_df) <- rownames(mat) <- colnames(mat) <- new_col_names
  pheatmap(mat, 
           clustering_distance_cols = "correlation",
           clustering_distance_rows = "correlation",
           clustering_method = clust_method,
           annotation_col = anno_df,
           annotation_colors = anno_cols,
           show_rownames = show_row_names,
           show_colnames = show_col_names,
           color = cols,
           border_color = NA)
}



library(RColorBrewer)
#greencols <- colorRampPalette(brewer.pal(11, 'PRGn'))(11)
greencols <- brewer.pal(n = 9, "BuGn")[c(1,3,5,7,9)]
agecolors <- greencols
names(agecolors) <- c("d12", "d22", "d27", "d32", "d37")
sexcolors <- c("tomato","steelblue")


# colors for pheatmap sidebars
age <- greencols
names(age) <- c("d12", "d22", "d27", "d32", "d37")
sex <- c("tomato","steelblue")
names(sex) <- c("F", "M")
anno_cols <- list(age = age, sex = sex)


plotHM2 <- function(counts, method="pearson", sampleCondition = sampleCondition_UTR){
  counts[counts==0] <- NA
  counts <- na.omit(counts)
  distsRL <- as.dist(1-cor(counts, method=method))
  #distsRL <- dist(t(counts))
  mat<- as.matrix(cor(counts, method=method))
 # disfun <- function(x) {as.dist(1-x)}
  #mat <- as.matrix(distsRL)
  rownames(mat) <- colnames(mat) <- colnames(counts)
  hc <- hclust(distsRL)

  hmap <- main_heatmap(mat, colors=hmcol) %>%
    add_row_labels(side="right") %>%
    add_col_labels() %>%
    add_col_annotation(sampleCondition[, c("age","sex")],
                       colors = list("age"=agecolors, "sex"=sexcolors)) %>%
    add_row_dendro(hc) %>%
    add_col_dendro(hc)
  return(hmap)  
}

hmcol2 <- colorRampPalette(brewer.pal(9, 'RdYlBu'))(100)
plotHM_expr <- function(plotdata, sampleCondition = sampleCondition_UTR, col_clust=F, k=4){
  #plotdata is the subset of logCPM
  #row.names(plotdata) <- genename[row.names(plotdata)]
  sampleCondition <- sampleCondition[colnames(plotdata),]
  hmap <- main_heatmap(plotdata, colors=rev(hmcol2)) %>%
    add_row_clustering(k=k) %>%
    add_row_labels(side="right", size=0.5) %>%
    add_col_labels() %>%
    add_col_annotation(sampleCondition[, c("age","sex")],
                       colors = list("age"=agecolors, "sex"=sexcolors)) %>%
    add_col_summary(groups=TRUE)
    
  if(col_clust){
    hmap <- hmap %>%
      add_col_clustering()
  }
  return(hmap)  
}



get_fisher <- function(all_detected_genes, all_detected_pu_genes, sig_genes, sig_pu_genes){
  #print(c(all_detected_genes, all_detected_pu_genes, sig_genes, sig_pu_genes))
  if(is.na(sig_genes)|is.na(sig_pu_genes)){
    return(NA)
  } else{
    no_sig_pu_genes <- all_detected_pu_genes - sig_pu_genes
    sig_no_pu_genes <- sig_genes - sig_pu_genes
    no_pu_no_sig_genes <- all_detected_genes - all_detected_pu_genes - sig_genes + sig_pu_genes
    test_mat <- matrix(c(sig_pu_genes,
                         no_sig_pu_genes, 
                         sig_no_pu_genes,
                         no_pu_no_sig_genes), nrow=2, byrow=T)
    res <- fisher.test(test_mat, alternative = "greater")
    return(res$p.value)
  }
}

gen_pca_plot <- function(count_table, npc=5, top=nrow(count_table), label=T, sample_label=T, colorF="tomato", colorM="steelblue", field=2){
  PCA <- get_pca_plotdata(count_table, npc=npc, top=top, field=field)
  var <- get_pca_var(PCA[[2]])
  
  p <- ggplot(as.data.frame(PCA[[1]]), aes(x=PC1, y=PC2, shape = age)) +
    geom_point(aes(color=sex), size=6, alpha=0.7) +
    theme_bw() +
    #geom_text_repel(aes(label=row.names(PCA[[1]]))) +
    scale_color_manual(values = c(colorF, colorM))+
    xlab(paste("PC1 (", 100*var[1], "%)")) +
    ylab(paste("PC2 (", 100*var[2], "%)")) 
    # scale_shape_manual(values = c(1, 2, 0, 3, 7)) +
    #theme(legend.position="none")
  if(label & sample_label){
    p <- p + geom_text_repel(aes(label=row.names(PCA[[1]])))
  } else if (label){
    p <- p + geom_text_repel(aes(label=rep))
  }
  
  # explore PCs
  # looking at PCA
  PCAdata <- melt(PCA[[1]], id.vars=c("id","group","age","sex", "rep"))
  
  idorder <- as.character(PCA[[1]][with(PCA[[1]], order(age, sex, rep)),]$id)
  PCAdata$id <- factor(PCAdata$id, levels=idorder)
  
  ppcs <- ggplot(PCAdata) +
    geom_bar(aes(x=id, y=value, fill=sex), stat="identity") +
    facet_grid(variable~age, scales="free") +
    theme(axis.text.x = element_text(angle=45, hjust=1))

  return(list(p, ppcs, PCAdata, PCA[[2]]))
}



get_sum_plots <- function(de_result_list, all_puberty_genes, sex_sample_range = 1:5, age_sample_range = 6:13, colorF="tomato", colorM="steelblue", fisher=T){
  gender_de_sig_list <- lapply(de_result_list[sex_sample_range], function(x) subset(x, FDR < 0.05 & abs(logFC) > log2(1.5)))
  
  gender_de_sig_df <- subset(melt(gender_de_sig_list), variable == "logFC")
  gender_de_sig_df$label <- ifelse(gender_de_sig_df$value > 0, "F_biased", "M_biased")
  
  gender_de_sig_df$puberty_genes <- ifelse(gender_de_sig_df$genename %in% all_puberty_genes, "yes", "no")
  
  sex_diff_count_table <- as.data.frame(with(gender_de_sig_df, table(L1, label)))
  names(sex_diff_count_table) <- c("age", "direction", "counts")
  write.table(sex_diff_count_table, paste0("sex_diff_geneCounts_table", Sys.Date(), ".txt"), sep="\t", col.names = T, row.names=F, quote=F)
  
  sex_sig_puberty_counts <- gender_de_sig_df %>%
    filter(puberty_genes == "yes") %>%
    group_by(L1, label) %>%
    summarise(counts=length(puberty_genes)) %>%
    as.data.frame()
  
  sex_counts_merged <- merge(sex_diff_count_table, sex_sig_puberty_counts, by.x=c("age", "direction"), by.y=c("L1", "label"), all.x=T, all.y=F)
  sex_counts_merged$counts.y[is.na(sex_counts_merged$counts.y)] <- 0
  
  if(fisher){
    sex_counts_merged$fishersP <- apply(sex_counts_merged, 1, function(x) get_fisher(length(all_detected_genes), length(all_detected_genes[all_detected_genes %in% all_puberty_genes]), as.numeric(x[3]), as.numeric(x[4])))
  }else{
    sex_counts_merged$fishersP <- 1
  }
  
  sex_counts_merged$fishersP_ori <- sex_counts_merged$fishersP
  sex_counts_merged$fishersP <- p.adjust(sex_counts_merged$fishersP, method="BH")
  
  ymax <- max(sex_counts_merged$counts.x)+10
  sex_counts_merged$age <- factor(gsub("_sex", "", sex_counts_merged$age), levels=c("d12","d22","d27","d32","d37"))
  p_sex_bias_count <- ggplot(sex_counts_merged) +
    geom_bar(aes(x=direction, y=counts.x, fill=direction), position="identity", stat="identity") +
    geom_point(aes(x=direction, y=counts.y), shape=17, size=1)+
    scale_fill_manual(values = c(colorF, colorM)) +
    scale_color_manual(values = c("white", "black")) +
    theme_bw() +
    facet_grid(.~age) +
    ylab("Number of sex-biased genes") +
    theme(axis.title.x=element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "top",
          legend.title = element_blank())
  
  pdata <- subset(na.omit(sex_counts_merged), fishersP < 0.1)
  if(nrow(pdata)>0){
    p_sex_bias_count <- p_sex_bias_count +
      geom_point(data=pdata, aes(x=direction, y=ymax), shape=8, size=2)
  }
  
  #ggsave(file="p_sex_biased_counts.pdf", p_sex_biase_count, width=5, height=3.5)
  ########################################################
  age_de_sig_list <- lapply(de_result_list[age_sample_range], function(x) subset(x, FDR < 0.05 & abs(logFC) > log2(1.5)))
  
  age_de_sig_df <- subset(melt(age_de_sig_list), variable == "logFC") %>%
    separate(L1, "_", into=c("age", "sex")) %>%
    as.data.frame()
  age_de_sig_df$label <- ifelse(age_de_sig_df$value > 0, "older_biased", "younger_biased")
  age_de_sig_df$puberty_genes <- ifelse(age_de_sig_df$genename %in% all_puberty_genes, "yes", "no")
  #gender_de_sig_df$puberty_genes <- ifelse(gender_de_sig_df$genename %in% all_puberty_genes, "yes", "no")
  
  age_diff_count_table <- age_de_sig_df %>%
    group_by(age, sex, label)%>%
    summarise(counts=length(value)) %>%
    as.data.frame()
  #colnames(sex_diff_count_table) <- c("age", "sex", "direction", "counts")
  write.table(age_diff_count_table, paste0("age_diff_geneCounts_table", Sys.Date(), ".txt"), sep="\t", col.names = T, row.names=F, quote=F)
  
  age_sig_puberty_counts <- subset(age_de_sig_df, puberty_genes == "yes")  %>%
    group_by(age, sex, label) %>%
    summarise(counts = length(puberty_genes)) %>%
    as.data.frame()
  
  age_counts_merged <- merge(age_diff_count_table, age_sig_puberty_counts, by=c("age", "sex", "label"), all.x=T, all.y=F)
  age_counts_merged$counts.y[is.na(age_counts_merged$counts.y)] <- 0
  if(fisher){
  age_counts_merged$fishersP <- apply(age_counts_merged, 1, function(x) get_fisher(length(all_detected_genes), length(all_detected_genes[all_detected_genes %in% all_puberty_genes]), as.numeric(x[4]), as.numeric(x[5])))
  } else{
    age_counts_merged$fishersP <- 1
  }
  
  age_counts_merged$fishersP_ori <- age_counts_merged$fishersP
  age_counts_merged$fishersP <- p.adjust(age_counts_merged$fishersP, method="BH")
  ymax <- max(sex_counts_merged$counts.x)+10
  
  p_age_bias_count <- ggplot(age_counts_merged) +
    geom_bar(aes(x=label, fill=label, y=counts.x), stat="identity", position="identity") +
    geom_point(aes(x=label, y=counts.y), shape=17, size=1)+
    scale_fill_manual(values = c("forestgreen", "orange")) +
    #scale_color_manual(values = c("white", "black")) +
    theme_bw() +
    facet_grid(sex~age, space="free") +
    ylab("Number of age-biased genes") +
    theme(axis.title.x=element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "top",
          legend.title = element_blank())
  if(nrow(subset(age_counts_merged, fishersP < 0.1)) >0){
    p_age_bias_count <- p_age_bias_count +
      geom_point(data=subset(age_counts_merged, fishersP < 0.1), aes(x=label, y=ymax), shape=8, size=2)
  }
  return(list(p_sex_bias_count, p_age_bias_count, sex_counts_merged, age_counts_merged))
}


# reshape normalized read counts table
ave_count_table <- function(logCPMc_used, scale=T){
  #logCPMc <- log(normCounts(set1) + 1)
  if(scale){
    logCPMc_used <- t(scale(t(logCPMc_used)))
  }
  logCPMc.m <- melt(logCPMc_used)
  
  logCPMc_ave_df <- logCPMc.m %>%
    separate(Var2, "_", into=c("ID", "sample")) %>%
    mutate(age = substr(sample, 2, 4), sex = substr(sample, 5,5), rep=substr(sample, 6, 6)) %>%
    group_by(Var1, age, sex) %>%
    summarize(median=median(value)) %>%
    spread(age, median) %>%
    as.data.frame()
  names(logCPMc_ave_df) <- c("gene","sex","d12","d22","d27","d32","d37")[1:ncol(logCPMc_ave_df)]
  logCPMc_ave_df$gene <- as.character(logCPMc_ave_df$gene)
  logCPMc_ave_df$sex <- as.character(logCPMc_ave_df$sex)
  return(logCPMc_ave_df)
}



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
    arrange(cluster, desc(domain), -p.value) %>%
    mutate(order = row_number())
  return(plotdata)
}



plot_cluster_enrichment <- function(enrich_list, outname, width=8, text_angle=0, horizontal =F, topn=5, ol_size=1, exclude_domain="", cluster_order=NULL, dot=F, colors = "Set2", save = TRUE){
  testclusters <- names(enrich_list[sapply(enrich_list, nrow) >0])
  enrich_data_list <- lapply(testclusters, function(x){
   # print(x)
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
      expand = c(0,0)) +
    scale_fill_brewer(palette = colors)
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
  if(save){
    pdf(paste0(outname, "_gprofiler_enrichment_top", topn, "_", Sys.Date(), ".pdf"), width=width, height = nrow(plotdata)*0.15+1)
    print(p)
    dev.off()
  } else{
    return(p)
  }
  
}

plot_cor <- function(de_result_list1, de_result_list2,condition1="d22_M", condition2="d22_F", merge_by = c("genes", "genename"), method="pearson", label1="", label2= ""){
  male22 <- de_result_list1[[condition1]]
  female22 <- de_result_list2[[condition2]] 
  merged <- merge(male22, female22, by=merge_by) %>% 
    mutate(sig1 = case_when(abs(logFC.x) >= log2FC_cutoff & FDR.x < FDR_cutoff ~ "sig",
                            TRUE ~ "noSig"),
           sig2 = case_when(abs(logFC.y) >= log2FC_cutoff & FDR.y < FDR_cutoff ~ "sig",
                            TRUE ~ "noSig")) %>% 
    mutate(sig = paste0(sig1, "_", sig2))
  
  cor <- round(with(merged, cor(logFC.x, logFC.y, method=method)), 2)
  corp <- round(with(merged, cor.test(logFC.x, logFC.y, method=method))$p.value, 2)
  
  
  pmerged <- ggplot() +
    geom_point(data=subset(merged, sig == "noSig_noSig", aes(x=logFC.x, y=logFC.y), alpha=0.1)) +
    geom_point(data=subset(merged, sig != "noSig_noSig"), aes(x=logFC.x, y=logFC.y, color = sig), alpha=0.6) +
    geom_point(data=subset(merged, abs(logFC.y) >= log2FC_cutoff & FDR.y < FDR_cutoff), aes(x=logFC.x, y=logFC.y), color="tomato", alpha=0.6) +
    geom_point(data=subset(merged, (abs(logFC.y) >= log2FC_cutoff & FDR.y < FDR_cutoff) & abs(logFC.x) >= log2FC_cutoff & FDR.x < FDR_cutoff), aes(x=logFC.x, y=logFC.y), color="purple", alpha=0.6) +
    annotate(geom = "text", x = 0.7*max(merged$logFC.x), y=0.9*max(merged$logFC.y), label=paste0(method, " cor=", cor, "\nP= ", corp), size=6) +
    theme_bw() +
    theme(axis.text = element_text(size=15),
          axis.title = element_text(size=15)) +
    xlab(paste0(condition1, label1, " log fold change")) +
    ylab(paste0(condition2, label2, "log fold change"))
  
  return(pmerged)
}

compare_qPCR <- function(logCPM, logCPMc, PCR_puberty, outname, method="spearman", ifgenplot=T){
 # PCR_puberty$Sample <- gsub("-","",PCR_puberty$Sample)
  exp_compare_ori <- merge_RNA_PCR(logCPM, PCR_puberty)
  exp_compare_corrected <- merge_RNA_PCR(logCPMc, PCR_puberty)
  
  all_variables <- sapply(strsplit(as.character(unique(colnames(logCPM))), "_"), "[[", 2)
  #all_samples <- as.character(unique(na.omit(PCR_puberty)$Sample))
  #sample_WL_match <- unlist(sapply(all_samples, function(x) all_variables[grepl(x, all_variables)]))
  
  cor_values_ori <- sapply(all_variables, get_cor_plot, exp_compare_ori, "number", method)
  cor_values_corrected <- sapply(all_variables, get_cor_plot, exp_compare_corrected, "number", method)
  
  exp_compare_cor_values <- as.data.frame(t(rbind(cor_values_ori, cor_values_corrected[1,])))
  colnames(exp_compare_cor_values) <- c("cor_original","n_PCR","n_RNA","cor_corrected")
  exp_compare_cor_values <- apply(exp_compare_cor_values, 2, as.numeric)
  row.names(exp_compare_cor_values) <- all_variables
  ### make a scatter plot for each samples, comparing UTR to qPCR. 
  
  
  plot_list_ori <- lapply(unique(na.omit(exp_compare_ori)$Sample), get_cor_plot, exp_compare_ori, method=method)
  plot_list_corrected <- lapply(unique(na.omit(exp_compare_corrected)$Sample), get_cor_plot, exp_compare_corrected, method=method)
  #save.image(paste0("pit_pairwise_analysis_",Sys.Date(),".RData"))
  names(plot_list_ori) <- names(plot_list_corrected) <- unique(na.omit(exp_compare_ori)$Sample)
  npage <- ceiling(length(all_variables)/9)
  
  if(ifgenplot){
    # output plots to file
    filename <- paste0(outname, "compare to PCR before batch correction ", Sys.Date(), ".pdf")
    pdf(filename, width=15, height=3*5, useDingbats = FALSE)
    for (i in 1:npage){
      if (i*9 >= length(plot_list_ori)){
        range <- c(((i-1)*9+1): length(plot_list_ori))
      } else{
        range <- c(((i-1)*9+1):(i*9))
      }
      #plot_list <- lapply(all_genes[range], function(x) factorplot(all_data_for_plot, x, ylab="delta_Ct", median=T))
      plot_list_temp <- plot_list_ori[range]
      args.list <- c(plot_list_temp,list(nrow=3,ncol=3))
      do.call(grid.arrange, args.list)
    }
    dev.off()
    
    filename <- paste0(outname, "compare to PCR after batch correction ", Sys.Date(), ".pdf")
    pdf(filename, width=15, height=3*5, useDingbats = FALSE)
    for (i in 1:npage){
      if (i*9 >= length(plot_list_corrected)){
        range <- c(((i-1)*9+1): length(plot_list_corrected))
      } else{
        range <- c(((i-1)*9+1):(i*9))
      }
      #plot_list <- lapply(all_genes[range], function(x) factorplot(all_data_for_plot, x, ylab="delta_Ct", median=T))
      plot_list_temp <- plot_list_corrected[range]
      args.list <- c(plot_list_temp,list(nrow=3,ncol=3))
      do.call(grid.arrange, args.list)
    }
    dev.off()
    return(exp_compare_cor_values)
  } else{
    return(list("cor_values" = exp_compare_cor_values,
                "plots_before_correction" = plot_list_ori,
                "plots_after_correction" = plot_list_corrected))
  }
 
}

compare_qPCR_summaryplot <- function(exp_compare_cor_values){
  exp_compare_cor_values.m <- melt(as.matrix(exp_compare_cor_values))
  sampleorder <- sort(row.names(exp_compare_cor_values))
  exp_compare_cor_values.m$Var1 <- factor(exp_compare_cor_values.m$Var1, levels = sampleorder)
  
  p0 <- ggplot(subset(exp_compare_cor_values.m, Var2 %in% c("cor_original","cor_corrected"))) +
    geom_density(aes(value, fill=Var2),  alpha=0.5) +
    scale_fill_manual(values=c("cor_original"="black", "cor_corrected"="red"))+
    theme_classic()+
   # ylab("cor coefficient") +
    theme(axis.text.x = element_text(size=12, color="black"),
          axis.text.y = element_text(size=12, color="black"),
          axis.title.x = element_blank(),
          legend.position = "none")
  
 # p0 <- ggdensity(subset(exp_compare_cor_values.m, Var2 %in% c("cor_original","cor_corrected")), x="value", fill="Var2", palette = c("black", "red"))
  
  p1 <- ggplot(subset(exp_compare_cor_values.m, Var2 %in% c("cor_original","cor_corrected"))) +
    geom_point(aes(x=Var1, y=value, color=Var2), alpha=0.6) +
    # geom_segment(aes(x=Var1, xend=Var1, y=0, yend=value, color=Var2))
    #ylim(c(0,1)) +
    scale_color_manual(values=c("cor_original"="black", "cor_corrected"="red"))+
    theme_bw()+
    ylab("cor coefficient") +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(size=12, color="black"),
          axis.title.x = element_blank(),
          legend.title = element_blank())
  p2 <- ggplot(subset(exp_compare_cor_values.m, Var2 %in% c("n_PCR","n_RNA"))) +
    geom_bar(aes(x=Var1, y=value, fill=Var2), color="gray30", stat="identity", position="identity", alpha=0.3) +
    scale_fill_manual(values=c("n_PCR"="white", "n_RNA"="red")) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, size=8, color="black"),
          axis.text.y = element_text(size=12, color="black"),
          legend.title = element_blank()) +
    ylab("detected genes") +
    xlab("samples")
  
 p0_v2 <- ggplot(subset(exp_compare_cor_values.m, Var2 %in% c("cor_corrected"))) +
    geom_density(aes(value, fill=Var2),  alpha=0.5) +
    scale_fill_manual(values=c("cor_original"="black", "cor_corrected"="red"))+
    theme_classic()+
    # ylab("cor coefficient") +
    theme(axis.text.x = element_text(size=12, color="black"),
          axis.text.y = element_text(size=12, color="black"),
          axis.title.x = element_blank(),
          legend.position = "none")
 p1_v2 <- ggplot(subset(exp_compare_cor_values.m, Var2 %in% c("cor_corrected"))) +
   geom_point(aes(x=Var1, y=value, color=Var2), alpha=0.6) +
   # geom_segment(aes(x=Var1, xend=Var1, y=0, yend=value, color=Var2))
   #ylim(c(0,1)) +
   scale_color_manual(values=c("cor_original"="black", "cor_corrected"="red"))+
   theme_bw()+
   ylab("cor coefficient") +
   theme(axis.text.x = element_blank(),
         axis.text.y = element_text(size=12, color="black"),
         axis.title.x = element_blank(),
         legend.title = element_blank())
  
  p <- plot_grid(p0,p1,p2, ncol=1, align="v", axis = "lr", rel_heights = c(1,1,2))
  p_correlatedonly <- plot_grid(p0_v2,p1_v2,p2, ncol=1, align="v", axis = "lr", rel_heights = c(1,1,2))
  return(list(p, p_correlatedonly))
}

get_ERCC_plots <- function(counts_cpm, outname){
  cpm_spikes <- count_cpm[spikes,]
  
  plotHM2(log(cpm_spikes+1), sampleCondition = sampleCondition) %>% save_iheatmap(paste0(outname, "ERCC_expr_cor_", Sys.Date(), ".png"),vwidth = 1100, vheight=992)  
  
  # compare ERCC amount and CPM
  cpm_spikes_df <- melt(log2(cpm_spikes+1))
  cpm_spikes_df$ERCC_amount <- ERCC_con_v[as.character(cpm_spikes_df[,1])]
  cpm_spikes_plot <- cbind(cpm_spikes_df, transform_df_samplename(samplenames_all = cpm_spikes_df$Var2, field = 2, list = F))
  
  
  get_cor <- function(sample){
    tempdata <- subset(cpm_spikes_df, Var2 == sample & value > 0)
    cor <- with(tempdata, cor(value, ERCC_amount))
    return(cor)
  }
  
  ERCCamount_cors <- sapply(unique(cpm_spikes_df$Var2), get_cor)
  names(ERCCamount_cors) <- unique(cpm_spikes_df$Var2)
  sort(ERCCamount_cors)
  
  cpm_spikes_plot$cor <- paste0("R=", round(ERCCamount_cors[as.character(cpm_spikes_plot$Var2)],2))
  cors_df <- unique(cpm_spikes_plot[, c("age", "sex", "rep", "cor")])
  
  ERCC_plot <- ggplot() +
    geom_point(data=subset(cpm_spikes_plot, value > 0), aes(x=ERCC_amount, y=value, color=group), alpha=0.6) +
    #geom_s(aes(group=id))
    geom_smooth(data=subset(cpm_spikes_plot, value > 0), aes(x=ERCC_amount, y=value, color=group, group=id), method="lm", alpha=0.6, se=F) +
    scale_x_continuous(breaks=scales::pretty_breaks(n=8))+
    scale_y_continuous(breaks=scales::pretty_breaks(n=8))+
    #facet_wrap(~group, ncol=5) +
    scale_color_brewer(palette = "Spectral")+
    ylab("ERCC_cpm")
  
  ggsave(paste0(outname, "ERCC_amount_reads_plot_", Sys.Date(), ".pdf"), ERCC_plot, width=5, height=4)
  
  ERCC_plot_all <- ERCC_plot + 
    geom_text(data=cors_df, aes(x=-6, y=10, label=cor), size=3) +
    facet_grid(age~sex+rep) +
    theme_bw()
  
  ggsave(paste0(outname, "ERCC_amount_reads_plot_all_", Sys.Date(), ".pdf"), ERCC_plot_all, width=15, height=8)
  return(cors_df)
}

ggHM <- function(plotdata, row_clust= T, col_clust=F, k = 0, midpoint=999, color_palette = "none"){
  library(dendextend)
  plotdata.m <- melt(as.matrix(plotdata))
  if(row_clust){
    rclust <- hclust(dist(plotdata))
    row.order <- as.character(row.names(plotdata)[rclust$order])
    plotdata.m$Var1 <- factor(plotdata.m$Var1, levels=row.order)
  }
  if(col_clust){
    cclust <- hclust(dist(t(plotdata)))
    col.order <- as.character(colnames(plotdata)[cclust$order])
    plotdata.m$Var2 <- factor(plotdata.m$Var2, levels=col.order)
  }
  
  #row_dendro <- ggdendrogram(as.dendrogram(rclust), rotate = T, labels = F, leaf_labels = F)
  
  
  rden <- rclust %>% as.dendrogram %>% 
    set("labels_cex", 0.5) 
  
  if(k!=0){
    rden <- rden %>% set("branches_k_color", k=k)
  }
  
  # trun into ggplot data
  rden_gg <- as.ggdend(rden) 
  
  # get dendro data for gene positioning 
  dend_data <- dendro_data(as.dendrogram(rclust))
  
  # Setup the data, so that the layout is inverted (this is more 
  # "clear" than simply using coord_flip())
  segment_data <- with(
    segment(dend_data), 
    data.frame(x = y, y = x, xend = yend, yend = xend))
  # Use the dendrogram label data to position the gene labels
  gene_pos_table <- with(
    dend_data$labels, 
    data.frame(y_center = x, gene = as.character(label), height = 1))
  
  # Table to position the samples
  sample_names <- colnames(plotdata)
  sample_pos_table <- data.frame(sample = sample_names) %>%
    mutate(x_center = (1:n()), 
           width = 1)
  
  # Neglecting the gap parameters
  heatmap_data <- as.matrix(plotdata) %>% 
    reshape2::melt(value.name = "expr", varnames = c("gene", "sample")) %>%
    left_join(gene_pos_table) %>%
    left_join(sample_pos_table)
  
  # Limits for the vertical axes
  gene_axis_limits <- with(
    gene_pos_table, 
    c(min(y_center - 0.5 * height), max(y_center + 0.5 * height))
  ) + 
    0.1 * c(-1, 1) # extra spacing: 0.1
  
  # Heatmap plot
  if(midpoint==999){
    midpoint=median(heatmap_data$expr)
  } 
  plt_hmap <- ggplot(heatmap_data, 
                     aes(x = x_center, y = y_center, fill = expr, 
                         height = height, width = width)) + 
    geom_tile() +
    scale_x_continuous(breaks = sample_pos_table$x_center, 
                       labels = sample_pos_table$sample, 
                       expand = c(0, 0)) + 
    # For the y axis, alternatively set the labels as: gene_position_table$gene
    scale_y_continuous(breaks = gene_pos_table$y_center, 
                       labels = rep("", nrow(gene_pos_table)),
                       limits = gene_axis_limits, 
                       expand = c(0, 0)) + 
    labs(x = "Sample", y = "") +
    theme_bw() +
    theme(axis.text.x = element_text(size = rel(1), hjust = 1, angle = 45), 
          # margin: top, right, bottom, and left
          plot.margin = unit(c(1, 0.2, 0.2, -0.7), "cm"), 
          panel.grid.minor = element_blank())
  
  if(color_palette != "none"){
    plt_hmap <- plt_hmap +
      scale_fill_distiller(palette = color_palette) 
  }else{
      scale_fill_gradient2("expr", high = "darkred", low = "darkblue", midpoint = midpoint) 
  }
  
  # Dendrogram plot using ggplot
  # plt_dendr <- ggplot(segment_data) + 
  #   geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
  #   scale_x_reverse(expand = c(0, 0.5)) + 
  #   scale_y_continuous(breaks = gene_pos_table$y_center, 
  #                      labels = gene_pos_table$gene, 
  #                      limits = gene_axis_limits, 
  #                      expand = c(0, 0)) + 
  #   labs(x = "Distance", y = "", colour = "", size = "") +
  #   theme_bw() + 
  #   theme(panel.grid = element_blank())
  
  # using dendextend
  plt_dendr <- ggplot(rden_gg, horiz = T) +
    scale_y_reverse(expand = c(0, 0.5)) + 
    scale_x_continuous(breaks = gene_pos_table$y_center, 
                       labels = gene_pos_table$gene, 
                       limits = gene_axis_limits, 
                       expand = c(0, 0)) + 
    labs(x = "Distance", y = "", colour = "", size = "") +
    theme_bw() + 
    theme(panel.grid = element_blank())
  
  library(cowplot)
  p <- plot_grid(plt_dendr, plt_hmap, align = 'h', rel_widths = c(1, 1))
  return(list(rclust, p))
}

# ave_count_table <- function(set1){
#   logCPMc <- log(normCounts(set1) + 1)
#   logCPMc.m <- melt(logCPMc)
#   
#   logCPMc_ave_df <- logCPMc.m %>%
#     separate(Var2, "_", into=c("ID", "sample")) %>%
#     mutate(age = substr(sample, 2, 4), sex = substr(sample, 5,5), rep=substr(sample, 6, 6)) %>%
#     group_by(Var1, age, sex) %>%
#     summarize(median=median(value)) %>%
#     spread(age, median) %>%
#     as.data.frame()
#   names(logCPMc_ave_df) <- c("gene","sex","d12","d22","d27","d32","d37")
#   logCPMc_ave_df$gene <- as.character(logCPMc_ave_df$gene)
#   logCPMc_ave_df$sex <- as.character(logCPMc_ave_df$sex)
#   return(logCPMc_ave_df)
# }

to_mouse <- function(gene){
  # This function converts gprofiler reported gene symbols (mouse symbols written in all caps) to mouse symbols used in the dataset 
  genes <- unlist(strsplit(gene, ","))
  out <- sapply(genes, function(x) paste(substr(x, 1, 1), tolower(substr(x, 2, nchar(x))), sep=""))
  return(unname(out))
}

library(pheatmap)


plot_gene_hm <- function(goi, used_data = logCPMc, name = F, ...){
  row.names(used_data) <- genename[row.names(used_data)]
  used_data <- used_data[, c(grep("F", colnames(used_data)), grep("M", colnames(used_data)))]
  if(name){
    goi <- genename[goi]
  }
  used_data_plot <- t(scale(t(used_data[goi[goi %in% row.names(used_data)],])))
  
  plot_phm(used_data_plot, ...)
}

anno_cols$type <- c(both="purple", "F_only"="tomato", "M_only"="steelblue")

# GWAS_cols <- colorRampPalette(brewer.pal(8,"Dark2"))(7)
# names(GWAS_cols) <- rev(sort(unique(all_puberty_genes_types)))
# anno_cols$GWAS <- GWAS_cols

darks <- colorRampPalette(brewer.pal(8,"Dark2"))(8)
names(darks) <- 1:8
anno_cols$clust <- darks

plot_phm <- function(plotdata, sampleCon = sampleCondition_fil, colors=rev(hmcol2), clust_cols=T,clust_rows=T, rownames=T, colnames=T, row_cut = 5, anno_rows=NULL, clust_distance_rows="euclidean", clust_distance_cols="euclidean", clust_method="complete", anno_colors = anno_cols, anno_df_cols =c("age", "sex"),...){

 
  anno_df <- sampleCon[colnames(plotdata), ]
  anno_df <- anno_df[anno_df_cols]
  
  all_values <- as.vector(as.matrix(plotdata))
  custome_breaks <- unique(c(min(all_values), seq(quantile(all_values, 0.005), 0, length = 125), seq(0, quantile(all_values, 0.995), length=125), max(all_values)))
  
  ptest <- pheatmap(plotdata,
                    color = colors,
                    breaks = custome_breaks,
                    cutree_rows = row_cut,
                    border_color = NA,
                    annotation_col = anno_df, 
                    annotation_colors = anno_cols,
                    cluster_cols = clust_cols,
                    cluster_rows = clust_rows,
                    show_rownames = rownames, 
                    show_colnames = colnames,
                    clustering_distance_cols = clust_distance_cols,
                    clustering_distance_rows = clust_distance_rows,
                    clustering_method = clust_method,
                    silent = TRUE)

  if(!is.null(anno_rows)){
    ptest <-  pheatmap(plotdata,
                       color = colors,
                       breaks = custome_breaks,
                       cutree_rows = row_cut,
                       border_color = NA,
                       annotation_col = anno_df, 
                       annotation_row = anno_rows,
                       annotation_colors = anno_colors,
                       cluster_cols = clust_cols,
                       cluster_rows = clust_rows,
                       show_rownames = rownames, 
                       show_colnames = colnames,
                       clustering_distance_cols = clust_distance_cols,
                       clustering_distance_rows = clust_distance_rows,
                       clustering_method = clust_method,
                       silent = TRUE,
                       ...)
  }
  
  return(ptest)
} 

# function to summarize DE genes numbers from different comparison
DE_number_counts <- function(sig_list_data, comparison){
  de_data <- sig_list_data[[comparison]]
  up <- as.numeric(table(de_data$logFC > 0)[2])
  down <- as.numeric(table(de_data$logFC > 0)[1])
  conditions <- unlist(strsplit(comparison, "_"))[1:2]
  out <- as.data.frame(rbind(c(conditions, up), c(rev(conditions), down)), stringsAsFactors=F)
  colnames(out) <- c("con1", "con2", "DE")
  return(out)
}

# function to generate a heatmap of all pair-wise comparisons
gen_de_number_hm <- function(sig_list_data, high_color = "steelblue", pu_genes, fisher=F){
  DE_count_table <- do.call(rbind, lapply(names(sig_list_data), function(x) DE_number_counts(sig_list_data, x)))
  DE_count_table$DE[is.na(DE_count_table$DE)] <- 0
  DE_count_table$DE <- as.numeric(as.character(DE_count_table$DE))
  
  sig_list_data_pu <- lapply(sig_list_data, function(x) subset(x, genes %in% pu_genes))
  DE_count_table_pu <- do.call(rbind, lapply(names(sig_list_data_pu), function(x) DE_number_counts(sig_list_data_pu, x)))
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
    scale_fill_gradient2(high=high_color, low="white", name="Number of DE genes") +
    scale_color_manual(values=c("yes"="red", "no"="black")) +
    theme(axis.title = element_blank())
  return(list(DE_p, DE_count_table)) 
}


get_sig_genes_df_from_list <- function(sig_list, pubery_genes = ""){
  de_merged <- melt(sig_list, d.vars=c("genes", "logFC"), measure.vars="logCPM") %>%
    select(c("genes", "logFC","L1")) %>%
    mutate(direction=ifelse(logFC>0, "up", "down")) %>%
    mutate(label = paste(L1, direction,sep="_")) %>%
    as.data.frame()
  temp_de_gene_list <- lapply(split(de_merged, de_merged$label), "[[", 1)
  
  sig_genes_all <- unique(unlist(temp_de_gene_list))
  sig_genes_df <- as.data.frame(do.call("cbind", lapply(temp_de_gene_list, function(x) ifelse(sig_genes_all %in% x, 1, 0))))
  row.names(sig_genes_df) <- sig_genes_all
  sig_genes_df$puberty <- ifelse(row.names(sig_genes_df) %in% pubery_genes, "yes", "no")
  return(sig_genes_df)
}

plot_clusters <- function(gene_clusters, ave, all_puberty_geneIDs=all_puberty_geneIDs, ncols=3, cluster_order=NULL){
  # pit_clusters is a vector where genenames are names and clusters are numbers
  cluster_names <- names(table(gene_clusters))
  clusters <- paste0("cluster", cluster_names)
  clustercounts <- paste0("n=", as.vector(table(gene_clusters)))
 # clustercounts <- clustercounts[cluster_names]
  names(clustercounts) <- clusters
  
  cluster_label <- function(variable, value){
    paste(value, ":",clustercounts[as.character(value)])
  }
  sig_genes <- names(gene_clusters)
  ave_plot <- subset(ave, gene %in% sig_genes)
  #row.names(pit_ave_sexdiff) <- pit_ave_sexdiff$gene
  #pit_ave_plot <- cbind(pit_ave_plot[,1:2], (pit_ave_plot[,3:7]-pit_ave_plot[,3]))
  ave_plot$cluster <- paste0("cluster", gene_clusters[ave_plot$gene])
  
  ave_plot <- melt(ave_plot, id.vars = c("gene","sex","cluster"))
  names(ave_plot)[4:5] <- c("age","logCPM")
  ave_plot$age <- as.numeric(gsub("d","",ave_plot$age))
  ave_plot_median <- ave_plot %>%
    group_by(sex, cluster, age) %>%
    summarise(median=mean(logCPM), sd = sd(logCPM), min=min(logCPM), max=max(logCPM)) %>%
    as.data.frame()
  ave_plot$pubertyGenes <- ifelse(ave_plot$gene %in% all_puberty_geneIDs, "yes", "no")
  
  if(!is.null(cluster_order)){
    cluster_order <- paste0("cluster", cluster_order)
    ave_plot$cluster <- factor(ave_plot$cluster, levels=cluster_order)
    ave_plot_median$cluster <- factor(ave_plot_median$cluster, levels = cluster_order)
  }
  
  ave_plot_median$sex <- factor(ave_plot_median$sex, levels = c("M","F"))
  
  pcluster <- ggplot() +
    #  geom_line(data=subset(pit_ave_plot, pubertyGenes == "no"), aes(x=age, y=logCPM, group=interaction(gene, sex), color=sex), alpha=0.5,size=0.5, linetype=2)+
    # geom_line(data=subset(pit_ave_plot, pubertyGenes == "yes"), aes(x=age, y=logCPM, group=interaction(gene, sex), color=sex), size=0.5, linetype=2)+
    #geom_line(data=pit_ave_plot, aes(x=age, y=logCPM, group=interaction(gene, sex), color=sex, alpha=pubertyGenes), size=0.5) +
    geom_line(data=ave_plot, aes(x=age, y=logCPM, group=interaction(gene, sex), color=sex), alpha=0.2,size=0.3, linetype=5) +
    #geom_ribbon(data=ave_plot_median, aes(age, ymin=median-sd, ymax=median+sd, fill=sex, group=sex), alpha=0.5) +
    geom_line(data=ave_plot_median, aes(x=age, y=median, color=sex, group=sex), size=1, alpha=0.7) +
    geom_point(data=ave_plot_median, aes(x=age, y=median, color=sex, shape=sex), size=2) +
    facet_wrap(~cluster,ncol = ncols, scales = "free_y", labeller = cluster_label) +
    theme_bw() +
    scale_x_continuous(breaks=c(12,22,27,32,37))  +
    scale_shape_manual(values = c("F"=16, "M"=17)) +
    #  scale_alpha_discrete(range=c(0.05, 0.5), guide=guide_legend(title="puberty gene")) +
    #  scale_linetype_manual(values = c("no"=2, "yes"=1), guide=guide_legend(title="puberty gene")) +
    scale_color_manual( values=c("F" = "tomato", "M" = "steelblue"))+
    scale_fill_manual( values=c("F" = "tomato", "M" = "steelblue"))
  return(pcluster)
}


calc_jaccard <- function(x, y){
  return(length(intersect(x, y))/length(union(x, y)))
}

enrichmap_like <- function(enrich_list_used, jaccard_cutoff=0.2, color_by="comparison", shape_by=NULL, size_by="pvalue", label=T, net=NULL, topn=NULL, ol_size=3, exclude_domain="", colors=NULL, low_limit = 0.2){
  # requires packages "tidygraph" "ggraph"
  if(is.null(net)){
    testclusters <- names(enrich_list_used[sapply(enrich_list_used, nrow) >0])
    if(!is.null(topn)){
      enrich_data_list <- lapply(testclusters, function(x){
        print(x)
        enrich_data <- enrich_list_used[[x]]
        #enrich_data <- subset(enrich_data, domain %in% used_domains)
        enrich_data$domain <- factor(enrich_data$domain, levels = as.character(unique(enrich_data$domain)))
        enrich_data <- do.call("rbind", lapply(split(enrich_data, f=enrich_data$domain), function(x) {
          x <- subset(x, `overlap.size` >= ol_size)
          x <- na.omit(x[order(x$p.value),][1:min(nrow(x), topn),])
          #x <- x[order(nrow(x):1),]
          return(x)
        }))
        # enrich_data$term.name <- factor(enrich_data$term.name, levels=rev(unique(as.character(enrich_data$term.name))))
        enrich_data$cluster <- rep(x, nrow(enrich_data))
        return(enrich_data)
      })
    } else{
      enrich_data_list <- lapply(testclusters, function(x){
        print(x)
        enrich_data <- enrich_list_used[[x]]
        enrich_data$cluster <- rep(x, nrow(enrich_data))
        return(enrich_data)
      })
    }
   
    plotdata <- do.call("rbind", enrich_data_list)
    plotdata <- process_plotdata(plotdata, exclude_domain=exclude_domain, testclusters = testclusters)
    plotdata$label <- with(plotdata, paste(`term.name`, cluster, domain, sep=";"))
    temp <- combn(unique(plotdata$label), m=2, simplify = F)
    
    term_jaccard <- function(labels){
      data1 <- to_mouse(subset(plotdata, label == labels[1])$intersection)
      data2 <- to_mouse(subset(plotdata, label == labels[2])$intersection)
      return(calc_jaccard(data1, data2))
    }
    
    nodes <- select(plotdata, id=label, termname=term.name, comparison=cluster, ngenes, pvalue=p.value, domain) %>%
      mutate(ngenes=as.numeric(gsub("n=","",ngenes)), pvalue=-log10(pvalue))
    
    jaccards <- sapply(temp, term_jaccard)
    edges <- as.data.frame(t(combn(unique(plotdata$label), m=2, simplify = T)))
    edges$jaccard <- jaccards
    names(edges) <- c("from", "to", "jaccard")
    edges$from <- as.character(edges$from)
    edges$to <- as.character(edges$to)
    
    edges_new <- subset(edges, jaccard > jaccard_cutoff) 
    nodes_new <- subset(nodes, id %in% c(edges_new[,1], edges_new[,2]))
    
    net <- tbl_graph(nodes = nodes, edges = edges_new, directed = F)
  }
  
  termnet <- ggraph(net) +
    geom_edge_link(aes(width=jaccard), alpha=0.3) +
    geom_node_point(aes_string(size=size_by, color=color_by, alpha=0.7)) +
    scale_size_continuous(range=c(1,8)) +
    #scale_fill_manual(values=c("miRNA"="tomato", "mRNA"="steelblue")) +
    #  scale_color_manual(values=c("No"= "black", "Yes"="red")) +
    #scale_shape_manual(values=c(16,17,18,15)) +
    scale_edge_width_continuous(limits = c(low_limit, 1), range=c(0.1,1.5)) +
    theme_bw()+
    theme(panel.background = element_blank(),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())
  if(label){
    termnet <- termnet +
      geom_node_text(aes_string(label="termname", color=color_by), repel = T) 
  }
  if(!is.null(colors)){
    termnet <- termnet +
      scale_color_manual(values = colors)
  }
  if(is.null(net)){
    return(list(net, termnet))
  } else{
    return(termnet)
  }
}

enrichmap_like_pie <- function(enrich_list_used, jaccard_cutoff=0.2, color_by="comparison", label=T, net_used=NULL, topn=NULL, ol_size=3, exclude_domain="", colors=NULL, low_limit = 0.2, radius_scale = 0.1, layout_used = "auto"){
  
  merge_genes <- function(x){
    x <- paste(x, collapse = ",")
    newstr <- paste(unique(unlist(strsplit(x, ","))), collapse=",")
    return(newstr)
  }
  testclusters <- names(enrich_list_used[sapply(enrich_list_used, nrow) >0])
# requires packages "tidygraph" "ggraph"
  if(is.null(net_used)){
    if(!is.null(topn)){
       enrich_data_list <- lapply(testclusters, function(x){
       print(x)
       enrich_data <- enrich_list_used[[x]]
      #enrich_data <- subset(enrich_data, domain %in% used_domains)
       enrich_data$domain <- factor(enrich_data$domain, levels = as.character(unique(enrich_data$domain)))
       enrich_data <- do.call("rbind", lapply(split(enrich_data, f=enrich_data$domain), function(x) {
        x <- subset(x, `overlap.size` >= ol_size)
        x <- na.omit(x[order(x$p.value),][1:min(nrow(x), topn),])
        #x <- x[order(nrow(x):1),]
        return(x)
      }))
      # enrich_data$term.name <- factor(enrich_data$term.name, levels=rev(unique(as.character(enrich_data$term.name))))
      enrich_data$cluster <- rep(x, nrow(enrich_data))
      return(enrich_data)
    })
  } else{
    enrich_data_list <- lapply(testclusters, function(x){
      print(x)
      enrich_data <- enrich_list_used[[x]]
      enrich_data$cluster <- rep(x, nrow(enrich_data))
      return(enrich_data)
    })
  }

  
  plotdata <- do.call("rbind", enrich_data_list)
  plotdata <- process_plotdata(plotdata, exclude_domain=exclude_domain, testclusters = testclusters)
  # if scatterpie
  term_with_genes <- group_by(plotdata, term.name) %>% 
    summarise(intersection = merge_genes(intersection)) %>% 
    rowwise() %>% 
    mutate(n_genes = length(unlist(strsplit(intersection, ","))))
  
  term_and_cluster <- plotdata %>% 
    dplyr::select(term.name, cluster) %>% 
    unique() %>% 
    group_by(term.name) %>% 
    mutate(size = 1/length(cluster)) %>% 
    spread(cluster, size) 
  term_and_cluster[is.na(term_and_cluster)] <- 0
  
  nodes <- left_join(term_and_cluster, term_with_genes, by = "term.name")
  
  temp <- combn(as.character(nodes$term.name), m = 2, simplify = FALSE)
  
  term_jaccard <- function(labels){
    data1 <- to_mouse(subset(nodes, term.name == labels[1])$intersection)
    data2 <- to_mouse(subset(nodes, term.name == labels[2])$intersection)
    return(calc_jaccard(data1, data2))
  }
  
  jaccards <- sapply(temp, term_jaccard)
  edges <- as.data.frame(t(combn(as.character(nodes$term.name), m=2, simplify = T)))
  edges$jaccard <- jaccards
  names(edges) <- c("from", "to", "jaccard")
  edges$from <- as.character(edges$from)
  edges$to <- as.character(edges$to)
  
  edges_new <- subset(edges, jaccard > jaccard_cutoff) 
  nodes_new <- subset(nodes, term.name %in% c(edges_new[,1], edges_new[,2]))
  
  net <- tbl_graph(nodes = nodes, edges = edges_new, directed = F)  
  }  
  
  if(!is.null(net_used)){
    net <- net_used
  }
  
  p <- ggraph(net, layout = layout_used) +
    geom_edge_link(aes(width=jaccard), alpha=0.3)  
  
  aa <- p$data
  aa$radius <- radius_scale*log2(aa$n_genes)
  
  relabel <- function(x){
    round(2^(x/radius_scale))
  }
  
  p <- p + geom_scatterpie(aes_(x=~x,y=~y, group = ~term.name, r = ~radius, alpha = .7), data=aa, cols=testclusters,color=NA,alpha=.7) +
    geom_scatterpie_legend(aa$radius, x = min(aa[, "x"]) -1, y = max(aa[, "y"])+1, labeller = relabel) +
    coord_equal() 
  p <- p + scale_edge_width_continuous(limits = c(low_limit, 1), range=c(0.1,1.5)) +
    theme_bw() +
    theme(panel.background = element_blank(),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank()) +
    scale_fill_brewer(palette = "Set2")
  
  

  if(!is.null(colors)){
    p <- p +
      scale_fill_manual(values = colors)
  }
  
  if(label){
    plabel <- p +
      geom_node_text(aes_string(label="term.name"), repel = T)
  }
  
  if(is.null(net_used)){
    return(list(net, p, plabel))
  } else{
    return(list(p, plabel))
  }
}

get_test_result <- function(contrasts=pair_contrasts, name){
  qlf <- glmQLFTest(fit, contrast = contrasts[, name])
  topTg <- topTags(qlf, n=nrow(y_fil$counts))
  
  de <- as.data.frame(topTg[[1]])
  de$genename <- genename[de$genes]
  return(de)
}

get_test_result_coef <- function(fit, y_fil, coef=NULL, contrast=NULL){
  # function to use multiple coefs or contrats
  if(!is.null(coef)){
    qlf <- glmQLFTest(fit, coef=coef)
  }else{
    qlf <- glmQLFTest(fit, contrast=contrast)
  }
  topTg <- topTags(qlf, n=nrow(y_fil$counts))
  
  de <- as.data.frame(topTg[[1]])
  de$genename <- genename[de$genes]
  return(de)
}


masigpro_runthrough <- function(data, set1,ndegree=3, counts=T){
  # data <- as.data.frame(normCounts(set1), stringsAsFactors = F)
  pdata <- pData(set1) 
  #pdata <- pdata[colnames(data),]
  design <- data.frame(sample = row.names(pdata)) %>%
    separate(sample, "_", into=c("ID", "condition")) %>%
    mutate(Time=substr(condition, 3,4)) %>%
    as.data.frame()
  
  sampleCondition <- data.frame(sample = colnames(data)) %>%
    separate(sample, "_", into=c("ID", "condition")) %>%
    mutate(age=substr(condition, 2, 4), sex=substr(condition, 5,5)) %>%
    as.data.frame()
  
  design$male <- ifelse(grepl("M", design$condition), 1, 0)
  design$female <- ifelse(grepl("F", design$condition), 1, 0)
  design$Replicate <- as.numeric(as.factor(substr(design$condition, 1, 5)))
  
  design <- design[, c("Time", "Replicate","male","female")]
  design <- apply(design, 2, as.numeric)
  row.names(design) <- row.names(pdata)
  
  designM <- make.design.matrix(design, degree = ndegree)
  
  # remove ERCC spike ins
  data <- data[!grepl("ERCC", row.names(data)), ]
  
  colnames(data) <- row.names(design)
  
  testp <- p.vector(data, designM, counts = counts)
  testT <- T.fit(testp)
  get <- get.siggenes(testT, vars="all",rsq=0.6)
  get$summary
  return(list(design, get, testp, testT))
}

get_time_results <- function(design, y_fil){
  y_fil <- estimateDisp(y_fil, robust = TRUE, design = design)
  #y_fil <- estimateGLMCommonDisp(y_fil, design)
  #y_fil <- estimateGLMTagwiseDisp(y_fil, design)
  fit <- glmQLFit(y_fil, robust=TRUE, design=design)
  
  # Make contrasts for all pair-wise age comparisions within female or male samples
  fconditions <- colnames(fit)[grepl("F", colnames(fit))]
  fconditions_pairs <- apply(combn(rev(fconditions), 2),2, paste, collapse="-")
  fcontrasts <- makeContrasts(contrasts=fconditions_pairs, levels = design)
  colnames(fcontrasts) <- gsub("group", "", colnames(fcontrasts))
  
  mconditions <- colnames(fit)[grepl("M", colnames(fit))]
  mconditions_pairs <- apply(combn(rev(mconditions), 2),2, paste, collapse="-")
  mcontrasts <- makeContrasts(contrasts=mconditions_pairs, levels = design)
  colnames(mcontrasts) <- gsub("group", "", colnames(mcontrasts))
  
  # Get result testing all differences between any pairs for male and female separately
  all_f_age_diff_result <- get_test_result_coef(fit=fit, y_fil=y_fil, contrast=fcontrasts)
  
  # same for male samples 
  all_m_age_diff_result <- get_test_result_coef(fit=fit, y_fil=y_fil, contrast=mcontrasts)
 
  return(list(all_f_age_diff_result, all_m_age_diff_result))
}



gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


get_gsea <- function(de_result,pathway, FC_col = "logFC") {
  colnames(de_result)[colnames(de_result) == FC_col] <- "logFC"
  rank_human_genes <- de_result %>%  
    arrange(-logFC) %>% 
    na.omit() %>% 
    dplyr::select("logFC", "genename") %>% 
    unique() 
  rank_human_genes <- rank_human_genes[!duplicated(rank_human_genes$genename),]
  
  rank_genes <- setNames(rank_human_genes$logFC, rank_human_genes$genename)
  
  gsea <- fgsea(pathway, rank_genes, minSize = 50, maxSize = 500, nperm = 1000)
  gsea <- gsea %>% 
    filter(padj < 0.1) %>% 
    arrange(-ES)
}

convertHumanGeneList <- function(x){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(genesV2)
}
