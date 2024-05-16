#Author: Mats Pettersson
#mats.pettersson@imbim.uu.se

#The two haplotype-level assemblies
require(Biostrings)
require(GenomicRanges)
IPA_pri <- readDNAStringSet("/Users/mapet205/Projects/Sprat/data/assemblies/IPA/final.p_ctg.fasta.gz")
IPA_alt <- readDNAStringSet("/Users/mapet205/Projects/Sprat/data/assemblies/IPA/final.a_ctg.fasta.gz")

load("~/Projects/Sprat/data/genotypes/Sprat_pool_freq_post_lo_GR_v2.RData")


#Milestone1 vs IPA alternative
Sprat_pin_curated <- readDNAStringSet("~/Projects/Sprat/data/HiC_satsuma/pin_HiC_it6/curation/Sprat_pin_milestone1.fasta.gz")
Sprat_pin_curated <- Sprat_pin_curated[order(Sprat_pin_curated@ranges@width, decreasing = T)]
Sprat_pin_curated_main_scaffs <- names(Sprat_pin_curated)[1:50]

IPA_alt_v_M1_satsuma <- read.table("~/Projects/Sprat/data/HiC_satsuma/Sprat_pin_m1_v_IPA_alt_satsuma.out_sorted.gz", stringsAsFactors=F, sep = "\t", comment.char = "")

IPA_alt_v_M1_satsuma_GR <- GRanges(seqnames = IPA_alt_v_M1_satsuma$V4, ranges = IRanges(start = IPA_alt_v_M1_satsuma$V5+1, end = IPA_alt_v_M1_satsuma$V6), seqinfo = seqinfo(Sprat_pin_curated))
IPA_alt_v_M1_red_GR <- GenomicRanges::reduce(IPA_alt_v_M1_satsuma_GR, min.gapwidth = 1e5)

IPA_alt_v_M1_gaps_GR <- GenomicRanges::gaps(IPA_alt_v_M1_red_GR)
IPA_alt_v_M1_gaps_GR <- IPA_alt_v_M1_gaps_GR[IPA_alt_v_M1_gaps_GR@strand == "*"]
IPA_alt_v_M1_gaps_df <- as.data.frame(IPA_alt_v_M1_gaps_GR)
IPA_alt_v_M1_gaps_agg <- aggregate(IPA_alt_v_M1_gaps_df$width ~ IPA_alt_v_M1_gaps_df$seqnames, FUN = "sum")
names(IPA_alt_v_M1_gaps_agg) <- c("scaffold", "gap_width")
#not the same order/size
all(as.character(IPA_alt_v_M1_gaps_agg[,1]) == names(Sprat_pin_curated))

IPA_alt_v_M1_gaps_agg$scaff_size <- width(Sprat_pin_curated)[match(as.character(IPA_alt_v_M1_gaps_agg[,1]), names(Sprat_pin_curated))]

IPA_alt_v_M1_gaps_agg$scaff_gap_frac <- IPA_alt_v_M1_gaps_agg$gap_width/IPA_alt_v_M1_gaps_agg$scaff_size
IPA_alt_v_M1_gaps_agg$high_gap_scaff <- F
IPA_alt_v_M1_gaps_agg$high_gap_scaff[IPA_alt_v_M1_gaps_agg$scaff_gap_frac > 0.8] <- T

sum(IPA_alt_v_M1_gaps_agg$scaff_size[IPA_alt_v_M1_gaps_agg$high_gap_scaff])/1e6


avg_depth <- rowSums(as.matrix(sprat_post_lo_GR$lo_GRs@elementMetadata[,grep("_DP", names(sprat_post_lo_GR$lo_GRs@elementMetadata))]), na.rm = T)/length(grep("_DP", names(sprat_post_lo_GR$lo_GRs@elementMetadata)))
smooth_win <- 100
roll_avg_depth <- stats::filter(avg_depth, rep(1/smooth_win, smooth_win))


sprat_HiC_GR <- GRanges(seqnames = sprat_post_lo_GR$lo_GRs$HiC_scaffold, ranges = IRanges(start = sprat_post_lo_GR$lo_GRs$HiC_pos, end = sprat_post_lo_GR$lo_GRs$HiC_pos))
low_depth_GR <- GenomicRanges::reduce(sprat_HiC_GR[which(roll_avg_depth < median(roll_avg_depth, na.rm =T))], min.gapwidth = 1e4)
low_depth_GR <- low_depth_GR[width(low_depth_GR) > 5e4]


low_depth_v_gaps <- pintersect(findOverlapPairs(low_depth_GR, IPA_alt_v_M1_gaps_GR))

depth_v_gaps_df <- as.data.frame(low_depth_v_gaps)
depth_v_gaps_agg <- aggregate(depth_v_gaps_df$width ~ depth_v_gaps_df$seqnames, FUN = "sum")
names(depth_v_gaps_agg) <- c("scaffold", "gap_width")
#not same order
all(as.character(depth_v_gaps_agg[,1]) == names(Sprat_pin_curated))

depth_v_gaps_agg$scaff_size <- width(Sprat_pin_curated)[match(depth_v_gaps_agg$scaffold, names(Sprat_pin_curated))]
depth_v_gaps_agg$scaff_gap_frac <- depth_v_gaps_agg$gap_width/depth_v_gaps_agg$scaff_size
depth_v_gaps_agg$high_gap_scaff <- F
depth_v_gaps_agg$high_gap_scaff[depth_v_gaps_agg$scaff_gap_frac > 0.5] <- T

sum(depth_v_gaps_agg$scaff_size[depth_v_gaps_agg$high_gap_scaff])/1e6
hist(depth_v_gaps_agg$scaff_size[depth_v_gaps_agg$high_gap_scaff]/1e6, breaks = 20)




#Duplicated BUSCOS
full_table_sprat_df <- read.table("~/Projects/Sprat/data/HiC_satsuma/pin_HiC_it6/curation/DeDup/full_table_sprat_final.fasta.out.tsv", sep = "\t", stringsAsFactors = F, fill = T, header = T)
sprat_dup_BUSCO_df <- full_table_sprat_df[full_table_sprat_df$Status == "Duplicated",]
sprat_dup_contig_mat <- table(sprat_dup_BUSCO_df[,c(1,3)])

f_BUSCO_overlap <- Vectorize(FUN = function(X,Y, data_mat){sum(data_mat[,X] == 1 & data_mat[,Y] == 1)},vectorize.args = c("X", "Y"))

sprat_dup_contig_raw_pairs <- outer(X = (1:dim(sprat_dup_contig_mat)[2]), Y = (1:dim(sprat_dup_contig_mat)[2]), FUN = "f_BUSCO_overlap", data_mat = sprat_dup_contig_mat)
colnames(sprat_dup_contig_raw_pairs) <- colnames(sprat_dup_contig_mat)
rownames(sprat_dup_contig_raw_pairs) <- colnames(sprat_dup_contig_mat)

sprat_dup_contig_pairs <- sprat_dup_contig_raw_pairs

for(i in 1:dim(sprat_dup_contig_pairs)[2]){
  tmp_dup_vec <- sprat_dup_contig_pairs[,i]/sprat_dup_contig_pairs[i,i]
  sprat_dup_contig_pairs[,i] <- tmp_dup_vec
}
  
sprat_dup_contig_comp <- sprat_dup_contig_pairs - t(sprat_dup_contig_pairs)
sprat_dup_contig_pair_rank <- sprat_dup_contig_pairs
sprat_dup_contig_pair_rank[sprat_dup_contig_comp < 0] <- "Major"
sprat_dup_contig_pair_rank[sprat_dup_contig_comp > 0] <- "Minor"
sprat_dup_contig_minor_vec <- colSums(sprat_dup_contig_pair_rank == "Minor")
sprat_dup_contig_major_vec <- colSums(sprat_dup_contig_pair_rank == "Major")

sum(sprat_dup_contig_major_vec == 0 & sprat_dup_contig_minor_vec > 0)

ambigous_contigs <- colnames(sprat_dup_contig_pair_rank)[sprat_dup_contig_major_vec > 0 & sprat_dup_contig_minor_vec > 0]

major_contigs <- colnames(sprat_dup_contig_pair_rank)[sprat_dup_contig_major_vec > 0 & sprat_dup_contig_minor_vec == 0]
minor_contigs <- colnames(sprat_dup_contig_pair_rank)[sprat_dup_contig_major_vec == 0 & sprat_dup_contig_minor_vec > 0]


sprat_dup_contig_pair_rank[sprat_dup_contig_pair_rank[,ambigous_contigs[1]] %in% c("Major", "Minor"),sprat_dup_contig_pair_rank[,ambigous_contigs[1]] %in% c("Major", "Minor")]
IPA_pri["ctg.000001F"]

#It seems we need an iterative process, where the state is re-evaluated, maybe even after each contig removal, and also a cutoff for the number of pairs needed to excise

#Skipping the above, and going on high fraction of missing pieces in the Alt & low depth instead.

Sprat_DeDup_HiC <- Sprat_pin_curated[!(names(Sprat_pin_curated) %in% depth_v_gaps_agg$scaffold[depth_v_gaps_agg$high_gap_scaff]) & width(Sprat_pin_curated) > 1e5]
Sprat_excised_HiC <- Sprat_pin_curated[!(!(names(Sprat_pin_curated) %in% depth_v_gaps_agg$scaffold[depth_v_gaps_agg$high_gap_scaff]) & width(Sprat_pin_curated) > 1e5)]


plot(sprat_dup_contig_major_vec - sprat_dup_contig_minor_vec)

for(pin_scaff in Sprat_pin_curated_main_scaffs){
  pin_scaff_filter <- sprat_post_lo_GR$lo_GRs$HiC_scaffold == pin_scaff
  #png(filename = paste0("~/Projects/Sprat/doc/assembly_curation/", pin_scaff,"_SNP_depth.png"), height = 500, width = 1000)
  png(filename = paste0("~/Projects/Sprat/doc/assembly_curation/", pin_scaff,"_SNP_depth_replot.png"), height = 500, width = 1000)
  plot(x =  sprat_post_lo_GR$lo_GRs$HiC_pos[pin_scaff_filter][1], y = avg_depth[pin_scaff_filter][1], pch = 20, cex = 0.5, col = "grey60", ylim = c(0, max(avg_depth, na.rm = T)), xlim = c(0, max(sprat_post_lo_GR$lo_GRs$HiC_pos[pin_scaff_filter], na.rm = T)), type = "n")
  #abline(h = median(roll_avg_depth, na.rm = T), col = "red", lwd = 3)
  #gap_scaff_filter <- IPA_alt_v_M1_gaps_GR@seqnames == pin_scaff
  #segments(x0 = start(IPA_alt_v_M1_gaps_GR[gap_scaff_filter]), x1 = end(IPA_alt_v_M1_gaps_GR[gap_scaff_filter]), y0 = par("usr")[4]*0.95, col = "steelblue", lwd = 3)
  gap_scaff_filter <- low_depth_v_gaps@seqnames == pin_scaff
  rect(xleft = start(low_depth_v_gaps[gap_scaff_filter]), xright = end(low_depth_v_gaps[gap_scaff_filter]), ytop = par("usr")[4]*0.95, ybottom = 0, border = NA, col = "grey70")
  points(x =  sprat_post_lo_GR$lo_GRs$HiC_pos[pin_scaff_filter], y = avg_depth[pin_scaff_filter], pch = 20, cex = 0.5, col = "grey30")
  
  points(x =  sprat_post_lo_GR$lo_GRs$HiC_pos[pin_scaff_filter], y = roll_avg_depth[pin_scaff_filter], pch = 20, cex = 0.7, col ="darkorchid")
  abline(h = median(avg_depth), col = "black", lwd = 2, lty = "dashed")
  
  dev.off()
}

#New BUSCO lists
after_BUSCOs <- read.csv("~/Projects/Sprat/data/HiC_satsuma/pin_HiC_it6/curation/DeDup/BUSCO_v5.sprat_after_deduplication.full_table.tsv", sep = "\t", header = T, stringsAsFactors = F, skip = 2) 
names(after_BUSCOs)[1] <- "Busco_id"
before_BUSCOs <- read.csv("~/Projects/Sprat/data/HiC_satsuma/pin_HiC_it6/curation/DeDup/BUSCO_v5.sprat_before_deduplication.full_table.tsv", sep = "\t", header = T, stringsAsFactors = F, skip = 2) 
names(before_BUSCOs)[1] <- "Busco_id"

sum(after_BUSCOs$Busco_id %in% before_BUSCOs$Busco_id)
sum(!(after_BUSCOs$Busco_id %in% before_BUSCOs$Busco_id))
lost_BUSCOs <- after_BUSCOs$Busco_id[after_BUSCOs$Status %in% c("Missing", "Fragmented")]
recoverable_BUSCOs <- before_BUSCOs[(before_BUSCOs$Busco_id %in% lost_BUSCOs) & !(before_BUSCOs$Status %in% c("Missing", "Fragmented")),]
unique(recoverable_BUSCOs$Busco_id)
table(recoverable_BUSCOs$Sequence[recoverable_BUSCOs$Status == "Complete"])
table(recoverable_BUSCOs$Sequence[recoverable_BUSCOs$Status == "Duplicated"])

recoverable_BUSCO_GR <- GRanges(seqnames = recoverable_BUSCOs$Sequence, ranges = IRanges(start = recoverable_BUSCOs$Gene.Start, end = recoverable_BUSCOs$Gene.End))
recoverable_BUSCO_GR$Busco_id <- recoverable_BUSCOs$Busco_id

full_table_sprat_GR <- GRanges(seqnames = full_table_sprat_df$Contig[!is.na(full_table_sprat_df$Start)], ranges = IRanges(start = full_table_sprat_df$Start[!is.na(full_table_sprat_df$Start)], end = full_table_sprat_df$End[!is.na(full_table_sprat_df$Start)]))
full_table_sprat_GR$Busco_id <- full_table_sprat_df$Busco_id[!is.na(full_table_sprat_df$Start)]

full_v_rec <- findOverlaps(full_table_sprat_GR, recoverable_BUSCO_GR)

recoverable_BUSCO_GR$Alt_ID <- NA
recoverable_BUSCO_GR$Alt_ID[full_v_rec@to] <- full_table_sprat_GR$Busco_id[full_v_rec@from]

recoverable_BUSCO_GR[as.character(recoverable_BUSCO_GR@seqnames) %in% major_contigs]

recover_ctg_tab <- table(as.character(recoverable_BUSCO_GR@seqnames))

recover_ctg_tab <- recover_ctg_tab[order(recover_ctg_tab, decreasing = T)]
recover_contigs <- names(recover_ctg_tab)[1:50]
recover_contigs <- recover_contigs[!(recover_contigs %in% minor_contigs)]

agp_culled <- Sprat_m1_AGP[Sprat_m1_AGP$chr_name %in% names(Sprat_excised_HiC),]
agp_recovery <- agp_culled[agp_culled$contig %in% recover_contigs,]

Sprat_recovery_HiC <- Sprat_excised_HiC[names(Sprat_excised_HiC) %in% agp_recovery$chr_name]

#Outputting current snapshot
Sprat_m2_HiC <- c(Sprat_DeDup_HiC,Sprat_recovery_HiC)
Sprat_m2_HiC <- Sprat_m2_HiC[order(Sprat_m2_HiC@ranges@width, decreasing = T)]
writeXStringSet(Sprat_m2_HiC, filepath = "~/Projects/Sprat/data/assemblies/Sprat_m2_HiC.fasta.gz", compress = T)
##Note name change!
#rsync -axv --progress --no-group ./Sprat_m2_HiC.fasta.gz matsp@rackham.uppmax.uu.se:/proj/snic2020-2-19/webexport/Sprat_assembly/dedup/Sprat_DeDup_v2_HiC.fasta.gz

#Re-importing the current version
Sprat_DeDup_v2 <- readDNAStringSet("/Users/mapet205/Projects/Sprat/data/assemblies/Sprat_DeDup_v2_HiC.fasta.gz")
#DeDup_v2 assembly stats
DeDup_v2_w <- width(Sprat_DeDup_v2)
DeDup_v2_w <- DeDup_v2_w[order(DeDup_v2_w, decreasing = T)]
DeDup_v2_50 <- sum(DeDup_v2_w)/2
DeDup_v2_c50 <- min(which(cumsum(DeDup_v2_w) > DeDup_v2_50))
DeDup_v2_w[DeDup_v2_c50]


#Preparing liftover to the Herring assembly
HiC_DeDup_v2_v_Ch_v2.0.2_satsuma <- satsuma_processing_v2("~/Projects/Sprat/data/HiC_satsuma/Ch_v2.0.2_v_HiC_DeDup_v2_satsuma.out_sorted.gz")
HiC_DeDup_v2_v_Ch_v2.0.2_satsuma <- satsuma_direction_est(HiC_DeDup_v2_v_Ch_v2.0.2_satsuma)

names(HiC_DeDup_v2_v_Ch_v2.0.2_satsuma)[3:6] <- paste0("sprat_HiC_DeDup_v2_", names(HiC_DeDup_v2_v_Ch_v2.0.2_satsuma)[3:6])
names(HiC_DeDup_v2_v_Ch_v2.0.2_satsuma)[8:11] <- paste0("herring_v2.0.2_", names(HiC_DeDup_v2_v_Ch_v2.0.2_satsuma)[8:11])
names(HiC_DeDup_v2_v_Ch_v2.0.2_satsuma)[8:11] <- sub("T_", "", names(HiC_DeDup_v2_v_Ch_v2.0.2_satsuma)[8:11])

save(HiC_DeDup_v2_v_Ch_v2.0.2_satsuma, file = "~/Projects/Sprat/data/HiC_satsuma/Ch_v2.0.2_v_HiC_DeDup_v2_satsuma.RData")

load(file = "~/Projects/Sprat/data/genotypes/Sprat_pool_freq_DeDup_v2.RData")

Sprat_pool_freq_HiC_DeDup_v2_lo <- Sprat_HiC_liftover_v2(sprat_pos_data = Sprat_pool_freq, liftover_df = HiC_DeDup_v2_v_Ch_v2.0.2_satsuma, chr_size_df = Ch_v2.0.2_sizes, lo_cols = 3:dim(Sprat_pool_freq)[2])
save(Sprat_pool_freq_HiC_DeDup_v2_lo, file = "~/Projects/Sprat/data/HiC_satsuma/Sprat_pool_freq_HiC_DeDup_v2_lo.RData")

load(file = "~/Projects/Sprat/data/genotypes/Sprat_pool_DeDup_v2_GR.RData") ##Current, HiC_DeDup_v2 version

#Harmonizing the data
require(dplyr)
sprat_Ch_SNP_IDs <- paste(Sprat_pool_freq_HiC_DeDup_v2_lo[,"CHROM"], Sprat_pool_freq_HiC_DeDup_v2_lo[,"POS"], sep = "_")
sprat_ch_pos <- left_join(x = as.data.frame(sprat_lo_GR@elementMetadata[, c("CHROM", "POS")]), y = Sprat_pool_freq_HiC_DeDup_v2_lo[!duplicated(sprat_Ch_SNP_IDs),c("CHROM", "POS", "herring_v2.0.2_seqnames", "SNP_HiC_pos")], by = c("CHROM", "POS"), keep = T)
sprat_lo_GR@elementMetadata$Ch_v2.0.2_chr <- sprat_ch_pos$herring_v2.0.2_seqnames
sprat_lo_GR@elementMetadata$Ch_v2.0.2_pos <- sprat_ch_pos$SNP_HiC_pos

# sprat_site_filter - from Sprat_pool_contrast.R
save(sprat_lo_GR, sprat_site_filter, file = "~/Projects/Sprat/data/genotypes/Sprat_pool_DeDup_v2_GR_Ch_lo.RData") ##Current, HiC_DeDup_v2 version

png(file = "~/Projects/Sprat/doc/assembly_curation/DeDup_v2_lo_Chr18.png", width = 1000, height = 500)
plot(x = sprat_lo_GR$Ch_v2.0.2_pos[which(sprat_lo_GR$Ch_v2.0.2_chr == "chr18" &  sprat_site_filter)], y = sprat_lo_GR$oce_v_brack[which(sprat_lo_GR$Ch_v2.0.2_chr == "chr18" &  sprat_site_filter)], pch = 16, cex = 0.5, main = "Chr 18", xlab = "Position", ylab = "Oce v Brackish")
dev.off()

BS_LAND_freq_lo_cols <- grep("([.]BS|LAND15|LAND19)_Afreq", names(sprat_lo_GR@elementMetadata), value = T)
Baltic_freq_lo_cols <- grep("(AB|GOTB|BBS|GD)_Afreq", names(sprat_lo_GR@elementMetadata), value = T)
brackish_no_BlackSea_freq_lo_cols <- grep("(AB|GOTB|BBS|GD|LAND15|LAND19)_Afreq", names(sprat_lo_GR@elementMetadata), value = T)

sprat_lo_GR$Baltic_mean <- rowMeans(as.matrix(sprat_lo_GR@elementMetadata[,Baltic_freq_lo_cols]), na.rm = T)

NS_freq_lo_cols <- grep("(NS[0-9]+)_Afreq", names(sprat_lo_GR@elementMetadata), value = T)
S_Atl_freq_lo_cols <- grep("(BoB|CEL)_Afreq", names(sprat_lo_GR@elementMetadata), value = T)
SK_freq_lo_cols <- grep("(SK[0-9]+)_Afreq", names(sprat_lo_GR@elementMetadata), value = T)

sprat_lo_GR$brackish_v_BlackSea <- abs(rowMeans(as.matrix(sprat_lo_GR@elementMetadata[,brackish_no_BlackSea_freq_lo_cols]), na.rm = T) - sprat_lo_GR@elementMetadata[,"P30.BS_Afreq"])
sprat_lo_GR$BS_LAND_v_Baltic <- abs(rowMeans(as.matrix(sprat_lo_GR@elementMetadata[,BS_LAND_freq_lo_cols]), na.rm = T) - rowMeans(as.matrix(sprat_lo_GR@elementMetadata[,Baltic_freq_lo_cols]), na.rm = T))
sprat_lo_GR$NS_v_S_Atl <- abs(rowMeans(as.matrix(sprat_lo_GR@elementMetadata[,NS_freq_lo_cols]), na.rm = T) - rowMeans(as.matrix(sprat_lo_GR@elementMetadata[,S_Atl_freq_lo_cols]), na.rm = T))
sprat_lo_GR$SK_v_S_Atl <- abs(rowMeans(as.matrix(sprat_lo_GR@elementMetadata[,SK_freq_lo_cols ]), na.rm = T) - rowMeans(as.matrix(sprat_lo_GR@elementMetadata[,S_Atl_freq_lo_cols]), na.rm = T))
sprat_lo_GR$N_v_S_Atl <- abs(rowMeans(as.matrix(sprat_lo_GR@elementMetadata[,c(SK_freq_lo_cols,NS_freq_lo_cols)]), na.rm = T) - rowMeans(as.matrix(sprat_lo_GR@elementMetadata[,S_Atl_freq_lo_cols]), na.rm = T))
sprat_lo_GR$Oce_v_Baltic <- abs(sprat_lo_GR$mean_oceanic - sprat_lo_GR$Baltic_mean)

save(sprat_lo_GR, sprat_site_filter, file = "~/Projects/Sprat/data/genotypes/Sprat_pool_DeDup_v2_GR_Ch_lo_ext.RData") ##Current, HiC_DeDup_v2 version


#Comparing with the DToL assembly

HiC_DeDup_v2_v_fSprSpr1.1_satsuma <- satsuma_processing_v2("~/Projects/Sprat/data/HiC_satsuma/fSprSpr1.1_v_HiC_DeDup_v2_satsuma.out_sorted.gz")
HiC_DeDup_v2_v_fSprSpr1.1_satsuma <- satsuma_direction_est(HiC_DeDup_v2_v_fSprSpr1.1_satsuma)

fSprSpr1.1 <- readDNAStringSet("/Users/mapet205/Projects/Sprat/data/assemblies/DToL/ncbi_dataset/ncbi_dataset/data/GCA_963457725.1/GCA_963457725.1_fSprSpr1.1_genomic.fna")



HiC_DeDup_main_scaffs <- names(table(HiC_DeDup_v2_v_fSprSpr1.1_satsuma$seqnames)[order(table(HiC_DeDup_v2_v_fSprSpr1.1_satsuma$seqnames),decreasing = T)][1:25])
fSprSpr1.1_main_scaffs <- grep("chrom", as.character(unique(HiC_DeDup_v2_v_fSprSpr1.1_satsuma$T_seqnames)), value = T)
fSprSpr1.1_main_scaffs <- fSprSpr1.1_main_scaffs[order(fSprSpr1.1_main_scaffs)]
DeDup_v2_v_fSprSpr1.1_satsuma_clean <- HiC_DeDup_v2_v_fSprSpr1.1_satsuma[(HiC_DeDup_v2_v_fSprSpr1.1_satsuma$seqnames %in% HiC_DeDup_main_scaffs) & (HiC_DeDup_v2_v_fSprSpr1.1_satsuma$T_seqnames %in% fSprSpr1.1_main_scaffs),]
DeDup_v2_v_fSprSpr1.1_satsuma_clean$seqnames <- as.character(DeDup_v2_v_fSprSpr1.1_satsuma_clean$seqnames)
DeDup_v2_v_fSprSpr1.1_satsuma_clean$T_seqnames <- as.character(DeDup_v2_v_fSprSpr1.1_satsuma_clean$T_seqnames)

DeDup_v2_v_fSprSpr1.1_mapping_summary <- table(DeDup_v2_v_fSprSpr1.1_satsuma_clean[,c("seqnames", "T_seqnames")])
DeDup_v2_v_fSprSpr1.1_mapping_summary["PGA_scaffold_1118__66_contigs__length_56234773", ]
DeDup_v2_v_fSprSpr1.1_mapping_summary["PGA_scaffold_1111__78_contigs__length_49772038", ]
DeDup_v2_v_fSprSpr1.1_mapping_summary[,"OY735285.1_Sprattus_sprattus_genome_assembly,_chromosome:_1"]
DeDup_v2_v_fSprSpr1.1_mapping_summary[,"OY735286.1_Sprattus_sprattus_genome_assembly,_chromosome:_2"]
DeDup_v2_v_fSprSpr1.1_mapping_summary[,"OY735287.1_Sprattus_sprattus_genome_assembly,_chromosome:_3"]
pdf(file = "~/Projects/Sprat/doc/assembly_curation/DeDup_v2_vs_fSprSpr1.1_hm.pdf", width = 12, height = 12)
heatmap(DeDup_v2_v_fSprSpr1.1_mapping_summary, scale = "col", margins = c(14, 14), cexRow = 0.6, cexCol = 0.6)
heatmap(DeDup_v2_v_fSprSpr1.1_mapping_summary, scale = "row", margins = c(14, 14), cexRow = 0.6, cexCol = 0.6)
dev.off()

#Looking at split Chr2
fSprSpr1.1[4]# Length: 71 378 072
#PGA_scaffold_1109__40_contigs__length_34818411 PGA_scaffold_1110__33_contigs__length_20919432 # length: 34818411 + 20919432 = 55737843
DeDup_v2_v_fSprSpr1.1_Chr2 <- DeDup_v2_v_fSprSpr1.1_satsuma_clean[DeDup_v2_v_fSprSpr1.1_satsuma_clean$T_seqnames == "OY735286.1_Sprattus_sprattus_genome_assembly,_chromosome:_2",]
DeDup_v2_v_fSprSpr1.1_Chr2 <- DeDup_v2_v_fSprSpr1.1_Chr2[DeDup_v2_v_fSprSpr1.1_Chr2$seqnames %in% c("PGA_scaffold_1109__40_contigs__length_34818411", "PGA_scaffold_1110__33_contigs__length_20919432"),]
DeDup_v2_v_fSprSpr1.1_Chr2$plot_col <- NA
DeDup_v2_v_fSprSpr1.1_Chr2$plot_col[DeDup_v2_v_fSprSpr1.1_Chr2$seqnames == "PGA_scaffold_1109__40_contigs__length_34818411"] <- "olivedrab"
DeDup_v2_v_fSprSpr1.1_Chr2$plot_col[DeDup_v2_v_fSprSpr1.1_Chr2$seqnames == "PGA_scaffold_1110__33_contigs__length_20919432"] <- "darkorchid"


DeDup_v2_v_fSprSpr1.1_Chr1 <- DeDup_v2_v_fSprSpr1.1_satsuma_clean[DeDup_v2_v_fSprSpr1.1_satsuma_clean$T_seqnames == "OY735285.1_Sprattus_sprattus_genome_assembly,_chromosome:_1",]
DeDup_v2_v_fSprSpr1.1_Chr1$plot_col <- "black"
DeDup_v2_v_fSprSpr1.1_Chr1$plot_col[DeDup_v2_v_fSprSpr1.1_Chr1$seqnames == "PGA_scaffold_1118__66_contigs__length_56234773"] <- "darkorchid"

DeDup_v2_v_fSprSpr1.1_Chr3 <- DeDup_v2_v_fSprSpr1.1_satsuma_clean[DeDup_v2_v_fSprSpr1.1_satsuma_clean$T_seqnames == "OY735287.1_Sprattus_sprattus_genome_assembly,_chromosome:_3",]
DeDup_v2_v_fSprSpr1.1_Chr3$plot_col <- "black"
DeDup_v2_v_fSprSpr1.1_Chr3$plot_col[DeDup_v2_v_fSprSpr1.1_Chr3$seqnames == "PGA_scaffold_1111__78_contigs__length_49772038"] <- "darkorchid"


pdf(file = "~/Projects/Sprat/doc/assembly_curation/DeDup_v2_vs_fSprSpr1.1_dot.pdf", width = 12, height = 12)

plot(y = DeDup_v2_v_fSprSpr1.1_Chr1$start, 
     x = DeDup_v2_v_fSprSpr1.1_Chr1$T_start,
     col = DeDup_v2_v_fSprSpr1.1_Chr1$plot_col, pch = 16, main = "Chr 1", xlab = "fSprSpr1.1 (DToL)", ylab = "DeDup_v2 (in house)")

plot(y = DeDup_v2_v_fSprSpr1.1_Chr2$start, 
     x = DeDup_v2_v_fSprSpr1.1_Chr2$T_start, 
     col = DeDup_v2_v_fSprSpr1.1_Chr2$plot_col, pch = 16, main = "Chr 2", xlab = "fSprSpr1.1 (DToL)", ylab = "DeDup_v2 (in house)")

plot(y = DeDup_v2_v_fSprSpr1.1_Chr3$start, 
     x = DeDup_v2_v_fSprSpr1.1_Chr3$T_start,
     col = DeDup_v2_v_fSprSpr1.1_Chr3$plot_col, pch = 16, main = "Chr 3", xlab = "fSprSpr1.1 (DToL)", ylab = "DeDup_v2 (in house)")

dev.off()

Sprat_HiC_liftover_v2 <- function(sprat_pos_data, liftover_df, chr_size_df, lo_cols){
  require(GenomicRanges)
  
  scaffold_GR <- GRanges(seqnames = sprat_pos_data$CHROM,ranges = IRanges(start = sprat_pos_data$POS, end = sprat_pos_data$POS))
  #BGI_v_Ilu_df
  liftover_GR <- GRanges(seqnames= liftover_df$sprat_HiC_DeDup_v2_seqnames,ranges=IRanges(start = liftover_df$sprat_HiC_DeDup_v2_start, end = liftover_df$sprat_HiC_DeDup_v2_end))
  
  #Matching SNP positions with entries in the BGI to Ilu satsuma alignment
  scaffold_matches <- findOverlaps(query = scaffold_GR, subject = liftover_GR, type = "any" )
  if(class(lo_cols) == "integer") lo_cols <- names(sprat_pos_data)[lo_cols]
  
  lo_df <- cbind(sprat_pos_data[scaffold_matches@from,c("CHROM","POS",lo_cols)], liftover_df[scaffold_matches@to,])
  
  #Adjusting positions within each matched interval
  target_SNPs <-  sign(lo_df$direction_est) == 1
  lo_df[target_SNPs,"SNP_HiC_pos"] <- lo_df[target_SNPs,"herring_v2.0.2_start"] + (lo_df[target_SNPs,"POS"] - lo_df[target_SNPs,"sprat_HiC_DeDup_v2_start"])
  target_SNPs <-  sign(lo_df$direction_est) != 1
  lo_df[target_SNPs,"SNP_HiC_pos"] <- lo_df[target_SNPs,"herring_v2.0.2_start"] + (lo_df[target_SNPs,"sprat_HiC_DeDup_v2_end"] - lo_df[target_SNPs,"POS"])
  
  #Add a cumulative position
  #HiC_chr_vec <- chromosome_list[order(as.numeric(sub("hic_scaffold_","",chromosome_list)))]
  pos_adj_vec <- c(0,cumsum(chr_size_df[,2]))
  pos_adj_df <- cbind(chr_size_df,pos_adj_vec[-length(pos_adj_vec)])
  lo_df[,"SNP_cumulative_pos"] <- lo_df[,"SNP_HiC_pos"] + pos_adj_df[match(lo_df[,"herring_v2.0.2_seqnames"], pos_adj_df[,1]),3]
  chr_col_vec <- c("grey30", "grey70")[match(lo_df[,"herring_v2.0.2_seqnames"], pos_adj_df[,1]) %% 2 + 1]
  lo_df[,"col"] <- chr_col_vec
  
  lo_df <- lo_df[!is.na(lo_df[,"SNP_cumulative_pos"]),]
  
  lo_df <- lo_df[order(lo_df[,"SNP_cumulative_pos"]),]
  return(lo_df)
}



