# load libraries
# suppressMessages(library(edgeR))
# library(derfinder)
# library(derfinderPlot)
# library(regionReport)
# library(GenomicRanges)
# library(GenomicFeatures)
# library(Biostrings)
# library(bumphunter)
# library(gplots)
# library(RColorBrewer)

used_libraries <- c("edgeR", "derfinder", "derfinderPlot", "regionReport", "GenomicRanges", "GenomicFeatures", "Biostrings", "bumphunter", "gplots","ggplot2","ggrepel", "gridExtra", "RColorBrewer", "dplyr", "tidyr", "UpSetR","reshape2")
library(BSgenome.Mmusculus.UCSC.mm10)

lapply(used_libraries, require, character.only = T, quietly = T)

refsequtr <- import("~/mdwilson/huayun/ganno/mm10/transcriptome/mm10_UCSC_RefSeq_3utr.bed", format="bed")

GWAS_genes_human_mouse_anno <- readRDS("~/Dropbox/wilson_lab/Huayun/puberty/datasets/GWAS_genes_human_mouse_anno.rds")

all_puberty_genes <- subset(GWAS_genes_human_mouse_anno, source !="Day2017_blood_eQTL")$MGI.symbol
all_puberty_genes_types <- subset(GWAS_genes_human_mouse_anno, source !="Day2017_blood_eQTL")$source
names(all_puberty_genes_types) <- all_puberty_genes
#library(TxDb.)

# gtf_file <- "~/mdwilson/genomes/mmus/gtf/gencode.vM11.annotation.gtf"
# txdb <- makeTxDbFromGFF(gtf_file, 
#                         format = "gtf",
#                         dataSource = "mouse vM11 transcript annotation",
#                         organism = "Mus musculus",
#                         chrominfo = chromInfo)
# alltx <- transcripts(txdb)
# 
# txgene <- select(txdb, keys=alltx$tx_name[1:10], columns=c("TXNAME","GENEID"), keytype="TXNAME")
# 
# genetype <- readRDS("~/mdwilson/huayun/ganno/mm10/transcriptome/mm10_gencode_vM11_geneType.rds")

hmcol<- colorRampPalette(brewer.pal(9, 'GnBu'))(100)

#batches <- readRDS("~/mdwilson/huayun/puberty/pit_rna/batches.rds")

# plotHM <- function(counts, method="pearson"){
#   distsRL <- as.dist(1-cor(counts, method=method))
#   #distsRL <- dist(t(counts))
#   mat<- as.matrix(cor(counts))
#   disfun <- function(x) {as.dist(1-x)}
#   #mat <- as.matrix(distsRL)
#   rownames(mat) <- colnames(mat) <- colnames(counts)
#   hc <- hclust(distsRL)
#   
#   #row_order <- hc$labels[hc$order]
#   # condition_df <- transform_df_samplename(colnames(mat), field = 2, list = F)
#   # sex_colors <- sapply(condition_df$sex, function(x) color_list[[as.character(x)]])
#   # age_colors <- sapply(condition_df$age, function(x) color_list[[as.character(x)]])
#   # batch_colors <- sapply(condition_df$batch, function(x) color_list[[as.character(x)]])
#   # colors <- cbind(sex_colors, age_colors, batch_colors)
#   
#   heatmap.2(mat, Rowv=as.dendrogram(hc),
#             Colv = as.dendrogram(hc),
#             distfun = disfun,
#             symm=TRUE, trace='none',
#             col = hmcol, margin=c(10, 13))#,
#   #           ColSideColors = colors,
#   #           ColSideColorsSize = 2)
#   # used_colors <- c(unique(as.character(condition_df$sex)),unique(as.character(condition_df$batch)),unique(as.character(condition_df$age))) 
#   # 
#   # used_color_list <- color_list[used_colors]
#   # 
#   # legend("topright",legend=names(used_color_list),
#   #        fill=unlist(used_color_list), border=TRUE, bty="n", y.intersp = 0.7, cex=0.7)
# }

single_hue_hm <- colorRampPalette(brewer.pal(n = 7, name ="Blues"))(100)

plotHM <- function(counts, method="pearson", cols=colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100), anno_df, show_col_names=T, show_row_names=T,new_col_names=colnames(counts), clust_method="complete"){
  mat<- as.matrix(cor(counts, method=method))
  
  # hc <- hclust(distsRL)
  #anno_df <- anno[, c("age","sex")]
  
  rownames(anno_df) <- rownames(mat) <- colnames(mat) <- new_col_names
  pheatmap(mat, 
           clustering_distance_cols = "correlation",
           clustering_distance_rows = "correlation",
           clustering_method = clust_method,
           annotation_col = anno_df,
           #annotation_colors = anno_cols,
           show_rownames = show_row_names,
           show_colnames = show_col_names,
           color = cols,
           border_color = NA)
}



get_samplesfiles <- function(sample_dir, pattern){
  samplefiles <- list.files(path=sample_dir,
                            pattern = glob2rx(pattern),
                            full.names = T)
  files <- rawFiles(datadir = NULL,
                    sampledirs = samplefiles,
                    fileterm = NULL)
  samplenames <- sapply(basename(files), function(x) paste(unlist(strsplit(x, "_"))[c(1,7)], collapse = "_"))
  names(files) <- samplenames
  return(files)
}



transform_df_samplename <- function(samplenames_all, field=7, list=TRUE){
  transformed_list <- lapply(samplenames_all, function(samplenames){
    name_elements <- lapply(as.character(samplenames), function(x) unlist(strsplit(x, "_")))
    samplenames <- sapply(name_elements, "[[", field)
    
    age <- sapply(samplenames, function(x) substr(unlist(strsplit(x, "-"))[1],2,4))
    sex <- sapply(samplenames, function(x) substr(unlist(strsplit(x, "-"))[1],5,5))
    rep <- sapply(samplenames, function(x) substr(x, nchar(x), nchar(x)))
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
    } else {
      return(data.frame(age=age, sex=sex, rep=rep, id=ID))
    }
  })
  outdf <- do.call("rbind", transformed_list)
  outdf$batch <- batches[as.character(outdf$id)]
  return(outdf)
}

get_txid <- function(anno, type="transcript_id"){
  anno_list <- unlist(strsplit(anno, "; "))
  output <- gsub(paste0(type, " "), "", anno_list[grepl(type, anno_list)])
  if(identical(output, character(0))){
    return(NA)
  } else{
    print(output)
    return(output)
  }
}

grange_add <- function(grange, add){
  meta1 <- as.data.frame(elementMetadata(grange))
  meta2 <- as.data.frame(add)
  meta <- cbind(meta1, meta2)
  elementMetadata(grange) <- meta
  return(grange)
}



grange_extend <- function(x, upstream=0, downstream=0){
  if (any(strand(x) == "*"))
    warning("'*' ranges were treated as '+'")
  on_plus <- strand(x) == "+" | strand(x) == "*"
  new_start <- start(x) - ifelse(on_plus, upstream, downstream)
  new_end <- end(x) + ifelse(on_plus, downstream, upstream)
  ranges(x) <- IRanges(new_start, new_end)
  trim(x)
}

extend_reduce <- function(anno_additional, left=40, right=20, ifgenename = T, ifgeneid=T){
  #function to reduce the extanded regions and match the gene id back
  anno_extended <- grange_extend(anno_additional, left, right)
  anno_ends_extended_reduced <- reduce(anno_extended)
  
  temp <- findOverlaps(anno_ends_extended_reduced, anno_extended, select="first")
  if(ifgeneid){
    mcols(anno_ends_extended_reduced)$gene_id <- mcols(anno_extended[temp])$gene_id 
  }
  if(ifgenename){
  mcols(anno_ends_extended_reduced)$gene_name <- mcols(anno_extended[temp])$gene_name
  }
 
  return(anno_ends_extended_reduced)
}


get_signal <- function(region,full_signal){
  chr <- as.character(seqnames(region))
  chr_signal <- full_signal[[chr]]
  signal <- colSums(as.data.frame(chr_signal[start(region):end(region),]))
  return (signal)
}

reformat_anno <- function(map_to_ends, txends, type="ends", nonUTR_regions=nonUTR_regions, additional=T, hitObject=T, distance){
  if(hitObject){
    anno_ends <- nonUTR_regions[from(map_to_ends)]
    ends_meta <- txends[to(map_to_ends)]
  } else{
    anno_ends <- nonUTR_regions
    ends_meta <- txends[map_to_ends]
  }
  
  #anno_end_txnames <- as.character(sapply(anno_ends$tx_name, function(x) unlist(strsplit(x, ";"))[1]))
  if(additional){
    anno_ends <- cbind(as.data.frame(anno_ends, row.names=NULL, stringsAsFactors=F), as.data.frame(ends_meta, row.names = NULL, stringsAsFactors=F)[, c("tx_id","tx_name")]) 
    anno_ends$gene_id <- txgenemap[as.character(anno_ends$tx_name)]
    anno_ends$gene_name <- genename[anno_ends$gene_id]
    anno_ends$type <- rep(type, nrow(anno_ends))
    anno_ends$distance <-distance
  } else{
    anno_ends <- cbind(as.data.frame(anno_ends, row.names=NULL, stringsAsFactors=F), gene_name = as.data.frame(ends_meta, row.names = NULL, stringsAsFactors=F)[, c("gene_name")]) 
  }
  
  # anno_ends2 <- anno_ends %>%
  #   group_by(seqnames, start, end, strand) %>%
  #   summarise(tx_name=paste(sort(tx_name), collapse=";"), gene_id = gene_id[1]) %>%
  #   mutate(gene_name = genename[gene_id], type=type) %>%
  #   as.data.frame()
  anno_ends <- makeGRangesFromDataFrame(anno_ends, keep.extra.columns = T)
  return(anno_ends)
}

grange_to_gtf <- function(refsequtr, type="RefSeq", outname){
  outgtf <- data.frame(chr=seqnames(refsequtr), 
                       type=type,
                       feature="exon",
                       start=start(refsequtr),
                       end=end(refsequtr),
                       noname1=".",
                       strand=strand(refsequtr),
                       noname2=".",
                       anno=paste('gene_id "', refsequtr$gene_id, '";',' gene_name "', refsequtr$gene_name, '";', sep=""))
  
  write.table(unique(outgtf), outname, quote=F, row.names=F, col.names = F, sep="\t")
}



# explore sequence composiiton of intronic signals
downstream_seq_sum_plot <- function(tx_regions, title=""){
  # tx_regions should have "region_name" column which is the index from all regions
  tx_downstream_seq <- downstream_seq[as.numeric(tx_regions$region_name)]
  letter_freq <- consensusMatrix(tx_downstream_seq)[1:4,]
  letter_freq <- as.data.frame(melt(letter_freq), stringAsFactors=F)
  names(letter_freq) <- c("Nucleotide", "pos","counts")
  letter_freq$pos <- as.numeric(letter_freq$pos) - 20
  
  # p <- ggplot(letter_freq) +
  #   geom_bar(aes(x=pos, y=counts, fill=letter), stat="identity", position="stack") +
  #   scale_x_continuous(breaks = seq(-20, 250, by = 20), expand = c(0,0)) +
  #   scale_y_continuous(expand = c(0,0)) +
  #   theme_linedraw() +
  #   xlab("relative position from ends of regions")
  
  p <- ggplot(letter_freq) +
    geom_line(aes(x=pos, y=counts, color=Nucleotide, group=Nucleotide)) +
    scale_x_continuous(breaks = seq(-20, 250, by = 20), expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    theme_bw() +
    xlab("relative position from ends of regions") +
    ggtitle(title) +
    theme(legend.position = c(0.9, 0.8),
          plot.title = element_text(hjust=0.5, margin=margin(t=10, b=-20)))
    #theme(legend.position = c(max(letter_freq$pos)-10, max(letter_freq$counts)-10))
  
  return(p)
}

