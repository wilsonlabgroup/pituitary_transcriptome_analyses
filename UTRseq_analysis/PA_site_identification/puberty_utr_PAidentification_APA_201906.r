library(derfinder)

#source("~/mdwilson/huayun/puberty/pit_rna/rh_pipeline/APA/derfinder_mouse_pit_functions.r")
source("~/Dropbox/wilson_lab/Huayun/puberty/scripts/puberty_utr_PAidentification_functions.r")

chrominfo_file <- "~/mdwilson/huayun/ganno/mm10/ChromInfo_noRandom.txt"
chromInfo <- read.table(chrominfo_file, sep="\t", as.is = T)
names(chromInfo) <- c("chrom", "length")

gtf_file <- "~/mdwilson/genomes/mmus/gtf/gencode.vM21.annotation.gtf"
#gtf_file <- "~/mdwilson/huayun/ganno/mm10/transcriptome/mm10_gencode_VM11_plus_UCSCrefseq_3utr_withTxNames.gtf"
txdb <- makeTxDbFromGFF(gtf_file, 
                        format = "gtf",
                        dataSource = "mouse vM21 transcript annotation",
                        organism = "Mus musculus",
                        chrominfo = chromInfo)
# generate chromatin state object using the txdb object. Function from derfinder package
#genomic_state_mm10 <- makeGenomicState(txdb, chrs = chromInfo$chrom)
#genomic_state <- readRDS("~/mdwilson/huayun/ganno/mm10/mm10_genomic_states_derfinder.rds")

#txgenemap <- readRDS("~/mdwilson/huayun/ganno/mm10/mm10_vM11_txgenemap.rds")

genename <- readRDS("~/mdwilson/huayun/ganno/mm10/transcriptome/mm10_gencode_vM21_geneName.rds")
genetype <- readRDS("~/mdwilson/huayun/ganno/mm10/transcriptome/mm10_gencode_vM21_geneType.rds")
rmsk <- readRDS("~/mdwilson/huayun/ganno/mm10/mm10_rmsk.rds")

blacklist <- read.table("~/mdwilson/genomes/mmus/mm10.blacklist.bed", sep="\t", as.is=T, stringsAsFactors = F)
names(blacklist) <- c("seqnames", "start", "end")
blacklist_gr <- makeGRangesFromDataFrame(blacklist)

genelength <- width(genes(txdb))
names(genelength) <- names(genes(txdb))


setwd("~/mdwilson/huayun/puberty/pit_rna/utrseq_2019/PAmethod/")
outname <- "pit_utr_regions_anno_"

#setwd("~/mdwilson/huayun/puberty/hypo_rna/utr_merge3prime/PAmethod/")
#outname <- "hypo_utr_regions_anno_to_intron_and_ends_counts_"
#sample_dir <- "~/mdwilson/huayun/puberty/hypo_rna/utr_merge3prime/raw_bw/"

# load sample files
sample_dir <- "~/mdwilson/huayun/puberty/pit_rna/utrseq_2019/raw_bw/"

files_plus <- get_samplesfiles(sample_dir, "*str1*")
files_neg <- get_samplesfiles(sample_dir, "*str2*")

#sample conditions
# parse sample names to get a sample information table
samplenames <- sapply(basename(files_plus), function(x) paste(unlist(strsplit(x, "_"))[c(1,7)], collapse = "_"))
#pheno <- data.frame(sample = sapply(samplenames, function(x) unlist(strsplit(x, "_"))[2]))
#row.names(pheno) <- samplenames

# get full coverage for all samples
all_chroms <- chromInfo$chrom
all_chroms <- all_chroms[all_chroms != "chrM"]


fullCov_pos <- fullCoverage(files = files_plus, chrs = all_chroms, verbose = FALSE)
fullCov_neg <- fullCoverage(files = files_neg, chrs = all_chroms, verbose = FALSE)

totalmapped_pos <- sapply(files_plus, getTotalMapped, chrs = all_chroms)
totalmapped_neg <- sapply(files_neg, getTotalMapped, chrs = all_chroms)

# # Get region matrix 
# DO NOT RUN
# test
#regionMat_pos <- regionMatrix(fullCov_pos, cutoff = 2, L = 1)
#regionMat_pos_lowcut <- regionMatrix(fullCov_pos, cutoff = 1, L = 1)
regionMat_pos <- regionMatrix(fullCov_pos, cutoff = 1, L = 1, totalMapped = totalmapped_pos, targetSize = 1000000)
regionMat_neg <- regionMatrix(fullCov_neg, cutoff = 1, L = 1, totalMapped = totalmapped_neg, targetSize = 1000000)

save(regionMat_pos, regionMat_neg, file = paste0(outname, Sys.Date(),".RData"))

region_list <- lapply(regionMat_pos, "[[", 1)
region_pos <- unlist(GRangesList(region_list), recursive = T, use.names = T)
strand(region_pos) <- "+"
#width_index_pos <- which(width(region_pos) > 10)
#region_pos <- region_pos[width_index_pos]

region_list_neg <- lapply(regionMat_neg, "[[", 1)
region_neg <- unlist(GRangesList(region_list_neg), recursive = T, use.names = T)
strand(region_neg) <- "-"

# get a list of all regions
regions_all_ori <- c(region_pos, region_neg)
regions_extend <- grange_extend(regions_all_ori, upstream = 50, downstream = 0)

# merge regions within 10bp of each other
regions_all <- reduce(regions_extend)

gaps_all <- width(gaps(regions_all))

as.data.frame(regions_all[which(gaps_all>10 & gaps_all < 20)[1:1000]])

# filter for blacklist regions
ol_bl <- findOverlaps(regions_all, blacklist_gr, select = "all")
regions_all <- regions_all[-unique(from(ol_bl))]

# rename regions
names(regions_all) <- seq(1,length(regions_all))
regions_all$region_name <- as.character(names(regions_all))
##############################################################
# get annotations
alltx <- transcripts(txdb)
txends <- resize(alltx, width=1, fix="end")

temp <- AnnotationDbi::select(txdb, keys=alltx$tx_name, keytype="TXNAME", columns=c("TXNAME", "GENEID"))
txgenemap <- temp$GENEID
names(txgenemap) <- temp$TXNAME

# extract genomic features
#txbygene <- transcriptsBy(txdb, "gene", use.names=T)


UTRs <- threeUTRsByTranscript(txdb, use.names=T)
UTRs <- unlist(UTRs, use.names = T)
mcols(UTRs)$gene_name <- genename[txgenename[names(UTRs)]]
#mcols(UTRs) <- mcols(UTRs)["gene_name"]

exons <- exonsBy(txdb, by="tx", use.names=T)
exons <- unlist(exons)

#last_exons <- endoapply(exons, function(x) x[length(x)])
#last_exons <- unlist(last_exons, use.names = T)

fiveUTR <- fiveUTRsByTranscript(txdb, use.names=T)
introns <- intronsByTranscript(txdb, use.names=T)



##############################################################
# refseq UTRs
refUTR <- read.table("~/mdwilson/huayun/ganno/mm10/transcriptome/mm10_UCSC_refSeq_3utr_1906.gtf", sep = "\t", as.is = T, stringsAsFactors = F)

refUTR_df <- data.frame(seqnames=refUTR$V1,
                        start=refUTR$V4,
                        end=refUTR$V5,
                        strand=refUTR$V7,
                        gene_name = gsub(";","",gsub("gene_name ", "", refUTR$V9)),
                        stringsAsFactors = F)

refUTR_gr <- makeGRangesFromDataFrame(refUTR_df, keep.extra.columns = T)

# save genome annotation results 
save(txgenemap,alltx, UTRs, exons, fiveUTR, introns, refUTR_gr, file = "~/mdwilson/huayun/ganno/mm10/transcriptome/mm10_vM21_genomic_annotation.RData" )

# merge 3utr annotation from gencode and UCSC(refseq)
# UTR_all <- c(UTRs, refUTR_gr)
# UTR_all <- UTR_all[!grepl("random", seqnames(UTR_all))]
# 
# # temp
# UTR_all <- UTRs
#################################################################################

# get downstream sequences
regions_all_rev <- invertStrand(regions_all)

# obtain regions 250bp downstream and 20bp upstream from the end of identified regions. 
downstream_regions <- promoters(regions_all_rev, upstream = 250, downstream =20)
#downstream_regions <- flank(regions_all, width=50, start=FALSE)
downstream_regions <- trim(downstream_regions) 
downstream_regions <- invertStrand(downstream_regions)
#strand(downstream_regions) <- unname(unlist(rev_strand[as.character(strand(downstream_regions))]))

# get sequences and filter for 18A matches
downstream_seq <- getSeq(BSgenome.Mmusculus.UCSC.mm10, downstream_regions)


#################################################################3
# annotate genes to 3utrs, exons, 5utrs, introns, and intergenic regions
map_to_gene_list <- lapply(list(UTRs, exons, fiveUTR, introns), function(x) unique(from(findOverlaps(regions_all, x))))
names(map_to_gene_list) <- c("threeUTR", "exon", "fiveUTR", "intron")
non_intergenic <- unique(unlist(map_to_gene_list))  
all_regions_index <- 1:length(regions_all)
map_to_gene_list$intergenic <- all_regions_index[!all_regions_index %in% non_intergenic]

# get a summary df of region annotations
regions_anno_counts_plot <- as.data.frame(do.call("cbind", lapply(map_to_gene_list, function(x) ifelse(all_regions_index %in% x, 1, 0))))

# make plots summarizing sequece composition downstream from each type of regions 
sequence_composiiton_plots <- lapply(names(map_to_gene_list), function(x) downstream_seq_sum_plot(regions_all[map_to_gene_list[[x]]], paste0("Regions mapped to ", x)))
library(gridExtra)
args.list <- c(sequence_composiiton_plots[c(1,4,5)],list(nrow=1,ncol=3))
# output plots 
pdf(paste0(outname, "_sequence_composition_plots_", Sys.Date(), ".pdf"), width=12, height=3)
do.call(grid.arrange, args.list)
dev.off()

##############################################################
# Labelling internal polyA events
# Count As up- and downstream from regions ends

# Search 150bp downstream 
downstream_seq_2 <- subseq(downstream_seq, start = 21, end = 170)
Aseq <- strrep("A", 18)
polyAcounts <- vcountPattern(Aseq, downstream_seq_2, max.mismatch = 6)

### Filter 2:
# Examine 20bp downtream 
# Any 6nt As within those regions?
Aseq <- strrep("A", 7)
sixAcounts <- vcountPattern(Aseq, subseq(downstream_seq, 21, 40), max.mismatch = 1)

A_freq <- letterFrequency(subseq(downstream_seq, 21, 40), "A", as.prob=TRUE)

### Filter 3: 
# look upstream for polyA signal, only considering AATAAA
upstream_regions <- promoters(regions_all_rev, upstream = 0, downstream = 70)
upstream_regions <- invertStrand(trim(upstream_regions))
upstream_seq <- getSeq(BSgenome.Mmusculus.UCSC.mm10, upstream_regions)
polyAsignal <- "AATAAA"
polyAsignal_counts <- vcountPattern(polyAsignal, upstream_seq, max.mismatch = 0) + vcountPattern("ATTAAA", upstream_seq, max.mismatch = 0)
regions_anno_counts_plot$polyAsignal <- polyAsignal_counts

# attach results from different filters to the summary table 
regions_anno_counts_plot$polyA_match <- polyAcounts
regions_anno_counts_plot$sixA_match <- sixAcounts
regions_anno_counts_plot$A_freq <- A_freq
#regions_anno_counts_plot$A_freq_pass <- A_freq_pass


myfun2 <- function(row, element, cutoff){
  data <- (row[element] > cutoff)
}

upset(regions_anno_counts_plot, sets=names(regions_anno_counts_plot)[1:5], order.by = "freq",
      queries = list(list(query = intersects, params = list("intron"), color="red")), expression = "sixA_match > 0")

# Plot the summary of annotation table, and all relavent filtering criteria. 
# plotting not working

# pdf(paste0(outname, "summary of region annotation ", Sys.Date(), ".pdf"), width=6, height=6, useDingbats = F)
# upset(regions_anno_counts_plot, sets = names(regions_anno_counts_plot)[1:5], 
#       query.legend = "top",
#       order.by = "freq",
#       queries = list(list(query = myfun2, params = list("threeUTR", 1), color="red", active=T, query.name = "3UTR"),
#                      list(query = myfun2, params = list("sixA_match", 1), color="purple", active=F, query.name = "six A_match > 0"),
#                      list(query = myfun2, params = list("polyA_match", 1), color="blue", active=F, query.name = "18-nt polyAs downstream"),
#                      list(query = myfun2, params = list("polyAsignal", 0), color="orange", active=T, query.name = "polyA signal upstream")),
#       #  boxplot.summary = c("AveSignal", "A_freq"))
#       boxplot.summary = "A_freq")
# dev.off()

# regions map to repeats 
map_to_rmsk <- unique(from(findOverlaps(regions_all, rmsk, select = "all")))
rmsk_ends <- resize(rmsk, width=1, fix="end")
distance_to_rmsk <- distanceToNearest(regions_all, rmsk_ends, select="arbitrary")
# either map to repeats or within 10bp downstream of a repeat 
map_to_rmsk <- unique(c(map_to_rmsk, from(subset(distance_to_rmsk, distance < 20))))

valid_index <- with(regions_anno_counts_plot, which(polyA_match == 0 & sixA_match==0 & A_freq < 0.5))


# region annotation
# all regions mapped to refseq UTRs but not GENCODE UTRs should be annotated
map_to_refseq <- findOverlaps(regions_all, refUTR_gr)
map_to_refseq_only <- unique(from(map_to_refseq)[!from(map_to_refseq) %in% map_to_gene_list$exon])
map_to_refseq_regions <- regions_all[map_to_refseq_only]

temp <- findOverlaps(map_to_refseq_regions, refUTR_gr)

# get the region:gene association based on refseq annotation
anno_to_refseq <- unique(reformat_anno(temp, refUTR_gr, type="refseq", map_to_refseq_regions, additional = F))

# Find the nearest gene ends
refseq_to_ends <- follow(map_to_refseq_regions, txends, select="last", ignore.strand=FALSE)
distance_to_ends_refseq <- distance(map_to_refseq_regions, txends[refseq_to_ends])
#distance_to_ends_refseq <- distance(map_to_refseq_regions[queryHits(refseq_to_ends)], txends[subjectHits(refseq_to_ends)])
#mcols(refseq_to_ends)$distance <- distance_to_ends_refseq
anno_refseq_to_ends <- reformat_anno(refseq_to_ends, txends, type="refseq_ends", map_to_refseq_regions, distance = distance_to_ends_refseq, hitObject = F)

###################################################################################################################
# map to transcripts
###################################################################################################################
# Annotate regions that map to refseq 3UTR annotation but do not map to GENCODE
refseq_to_introns <- findOverlaps(map_to_refseq_regions, alltx, select = "all")
anno_refseq_to_introns <- reformat_anno(refseq_to_introns, alltx, type="refseq_introns", map_to_refseq_regions, distance=NA)

temp <- sapply(1:length(anno_refseq_to_ends), function(x) anno_refseq_to_ends[x]$gene_name == anno_to_refseq[x]$gene_name)
refseq_to_ends_valid <- anno_refseq_to_ends[temp]

temp2 <- as.data.frame(anno_to_refseq) %>%
  select(c("region_name", "gene_name")) %>%
  left_join(unique(as.data.frame(anno_refseq_to_introns)[, c("region_name", "gene_name")]), by="region_name") %>%
  filter(gene_name.x == gene_name.y)

refseq_to_introns_valid <- subset(anno_refseq_to_introns, region_name %in% temp2$region_name)
refseq_valid <- rbind(as.data.frame(refseq_to_ends_valid[, c("region_name", "gene_id", "gene_name", "type")]),
                             as.data.frame(refseq_to_introns_valid[, c("region_name", "gene_id", "gene_name", "type")])) %>%
  group_by(seqnames, start, end, width, strand, region_name, gene_id, gene_name) %>%
  summarise(type=paste(type, collapse=";"))
refseq_valid_index <- refseq_valid$region_name
###################################################################################################################
# Map regions to UTRs
map_to_UTR_genes <- findOverlaps(regions_all[map_to_gene_list$exon], alltx)
anno_UTR <- reformat_anno(map_to_UTR_genes, alltx, type="tx", regions_all[map_to_gene_list$exon], distance = NA)
UTR_genes <- unique(c(refseq_valid$gene_id, as.character(anno_UTR$gene_id)))

# annotate intronic peaks
#map_to_intron <- map_to_gene_list$intron[!map_to_gene_list$intron %in% c(unlist(map_to_gene_list[c("threeUTR", "exon", "fiveUTR")]), map_to_rmsk)]
map_to_intron <- map_to_gene_list$intron[!map_to_gene_list$intron %in% c(unlist(map_to_gene_list[c("threeUTR", "exon", "fiveUTR")]), map_to_rmsk, refseq_valid_index) & map_to_gene_list$intron %in% valid_index]
map_to_intron_regions <- regions_all[map_to_intron]
map_to_tx <- findOverlaps(map_to_intron_regions, alltx, select="all")
distance_to_nearst_exon <- GenomicRanges::distanceToNearest(map_to_intron_regions[(from(map_to_tx))], exons, select=c("arbitrary"))
#mcols(map_to_tx)$distance <- mcols(distance_to_nearst_exon)$distance

anno_tx <- reformat_anno(map_to_tx, alltx, type="tx", map_to_intron_regions, distance = mcols(distance_to_nearst_exon)$distance)
# filter by distance 
#anno_tx_valid <- subset(anno_tx, distance > 40)
tx_genes <- unique(anno_tx$gene_id)

# annotate to intergenic regions
map_to_intergenic <- map_to_gene_list$intergenic[!map_to_gene_list$intergenic %in% c(map_to_rmsk, refseq_valid_index) & map_to_gene_list$intergenic %in% valid_index]
map_to_intergenic_regions <- regions_all[map_to_intergenic]
map_to_ends <- follow(map_to_intergenic_regions, txends, select="all", ignore.strand=FALSE)
distance_to_ends_downstream <- distance(map_to_intergenic_regions[queryHits(map_to_ends)], txends[subjectHits(map_to_ends)])
#mcols(map_to_ends)$distance <- distance_to_ends_downstream

# only considering regions within 1kb from its upstream transcript ends
#map_to_ends_valid <- subset(map_to_ends, distance > 40 & distance < 10000)

anno_ends <- reformat_anno(map_to_ends, txends, type="ends", map_to_intergenic_regions, distance= distance_to_ends_downstream)

mcols(anno_ends)$geneLength <- genelength[mcols(anno_ends)$gene_id]
mcols(anno_ends)$geneType <- genetype[mcols(anno_ends)$gene_id]
mcols(anno_ends)$perc <- with(mcols(anno_ends), distance/geneLength)

anno_ends_valid <- subset(anno_ends, geneType == "protein_coding" & distance < 5000)
#anno_ends_valid <- subset(anno_ends_valid, perc < 0.2 | distance < 5000)

ends_genes <- unique(anno_ends_valid$gene_id)

datalist <- list(UTR_genes,tx_genes, ends_genes)
names(datalist) <- c("Exonic signal", "Intronic signal", "Downstream signal")
pdf(paste0(outname, "_overlapping_genes_detected_diff_methods_", Sys.Date(),".pdf"), width=4, height = 4)
plot(venn(datalist))
dev.off()
# 
# intron_only_genes <- tx_genes[!tx_genes %in% c(UTR_genes, ends_genes)]
# end_only_genes <- ends_genes[!ends_genes %in% c(UTR_genes, tx_genes)]
# 
# 
# # plot the sequence composition after filter 
# p1 <- downstream_seq_sum_plot(unique(regions_all[anno_tx_valid$region_name]))
# p2 <- downstream_seq_sum_plot(unique(regions_all[anno_ends_valid$region_name]))
# 
# args.list <- c(list(p1, p2),list(nrow=1,ncol=2))
# # output plots 
# pdf(paste0(outname, "_after_filter_sequence_composition_plots_", Sys.Date(), ".pdf"), width=6, height=3)
# do.call(grid.arrange, args.list)
# dev.off()


# get the end of novel signal
# resize intronic regions
intronic_all <- c(makeGRangesFromDataFrame(subset(refseq_valid, type =="refseq_introns"), keep.extra.columns = T),
                  anno_tx[,c("region_name", "gene_id", "gene_name", "type")])


intron_novel <- unique(grange_extend(intronic_all, -45, 5)[, c("gene_id", "gene_name")])

ends_all <- c(makeGRangesFromDataFrame(subset(refseq_valid, type !="refseq_introns"), keep.extra.columns = T),
              anno_ends_valid[,c("region_name", "gene_id", "gene_name", "type")])

ends_of_novel <- resize(ends_all, width = 1, fix="end")
ends_of_novel_toUpstream_txend <- follow(ends_of_novel, txends, select="last", ignore.strand=FALSE)
ends_of_novel_anno <- cbind(as.data.frame(ends_of_novel), txend = start(txends[ends_of_novel_toUpstream_txend]))

new_df <- data.frame(seqnames=ends_of_novel_anno$seqnames, 
                     start=apply(ends_of_novel_anno, 1, function(x) as.integer(min(x["start"], x["txend"]))),
                     end=apply(ends_of_novel_anno, 1, function(x) as.integer(max(x["end"], x["txend"]))),
                     strand=ends_of_novel_anno$strand,
                     gene_id = ends_of_novel_anno$gene_id,
                     gene_name = ends_of_novel_anno$gene_name)

ends_novel <- makeGRangesFromDataFrame(new_df, keep.extra.columns = T)
ends_novel <- extend_reduce(ends_novel, 0,0)

intron_novel$type <- "intronic"
ends_novel$type <- "intergenic"

novel_utrs_all <- c(intron_novel, ends_novel)
saveRDS(novel_utrs_all, paste0(outname, "_novel_utrs_all_", Sys.Date(), ".rds"))
save.image(paste0(outname, "analysis_", Sys.Date(), ".RData"))


# output grange object to gtf 
grange_to_gtf(novel_utrs_all, "refseq_novel", paste0("mm10_GENCODE_Vm21_refSeq1906_pit_novel.gtf"))

#######################################################################################################
