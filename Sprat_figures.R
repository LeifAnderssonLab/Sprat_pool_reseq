#Author: Mats Pettersson
#mats.pettersson@imbim.uu.se


#Figures - Sprat adaption to brackish

#Necessary data
load("~/Projects/Sprat/data/sprat_sample_info.RData")

#Older, outdeted versions
#load("~/Projects/Sprat/data/PacBio_satsuma/Sprat_pool_freq_lo.RData")
#load(file = "~/Projects/Sprat/data/PacBio_satsuma/sprat_site_filter_lo.RData")
#load(file = "~/Projects/Sprat/data/genotypes/Sprat_pool_freq_post_lo_GR.RData")
#load(file = "~/Projects/Sprat/data/genotypes/Sprat_pool_freq_post_lo_GR_ext2.RData")
#load(file = "~/Projects/Sprat/data/genotypes/sprat_site_filter_m1.RData")
#load("~/Projects/Sprat/data/genotypes/Sprat_pool_freq_post_lo_GR_Ch_pos_v2.RData")
###

#load("~/Projects/Sprat/data/genotypes/Sprat_pool_DeDup_v2_GR_Ch_lo.RData") #Current, HiC_DeDup_v2 version
load("~/Projects/Sprat/data/genotypes/Sprat_pool_DeDup_v2_GR_Ch_lo_ext.RData") #Current, HiC_DeDup_v2 version, with more contrasts pre-made

load("~/Projects/Herring/data/polygenic/selection_regions.RData")
regions_raw <- readRDS("~/Projects/Herring/data/British_populations/Top50_SNPs_Cluhar_72pools_50contrasts_sigBonf.minGapWidth.50K.minWidth.2.minNdAFgt04.1.minNreg.2.topSNPs.mindAF.037.minNdAFgt04.1.RDS")
regions_raw$clean_reg <- sub("([A-Za-z_]+[0-9]+[:][0-9]+[-][0-9]+)[.].+","\\1", regions_raw$region)
mixed_reg_GR <- GenomicRanges::reduce(GRanges(unique(regions_raw$clean_reg)), min.gapwidth = 1e3)

#plot order for heatmaps
pools_to_plot <- rev(sprat_sample_info$clean_name[c(10:11,9,8, 12,13,5, 14,1:4,19,6:7,15:18)])

#Fig 1 - general tree
#load(file = "~/Projects/Sprat/data/genotypes/Sprat_pool_raw_dist.RData") #Original, outdated version
load(file ="~/Projects/Sprat/data/genotypes/Sprat_pool_raw_dist_DeDup_v2.RData") #Current, HiC_DeDup_v2 version

require(Biostrings)
require(ape)
#HapDNA <- DNAStringSet(hap_vec)
#HapDNA_bin <- as.DNAbin(HapDNA)
Sprat_pool_dist <- as.dist(Sprat_raw_dist)
Sprat_tree <- bionj(Sprat_pool_dist)

#Original color scheme
#tip_col_vec <- rep("darkorange1", length(tree$tip.label))
#tip_col_vec[grep("Baltic", tree$tip.label)] <- "darkorange4"
#tip_col_vec[grep("Landvik", tree$tip.label)] <- "darkorange1"

#tip_col_vec[grep("HWS", tree$tip.label)] <- "olivedrab4"
#tip_col_vec[grep("PB8", tree$tip.label)] <- "olivedrab2"
#tip_col_vec[grep("HWS6", tree$tip.label)] <- "darkolivegreen"

#plot_tree <- tree 
#plot_tree$tip.label <- rep("\U2022", length(tree$tip.label))
tree_sprat_IDs <- sub("_Afreq", "", grep("_Afreq", Sprat_tree$tip.label, value = T))
tree_sprat_pool_names <- sprat_sample_info$clean_name[match(tree_sprat_IDs, sprat_sample_info$POOL)]
Sprat_tree$tip.label <- tree_sprat_pool_names

png(file = "~/Projects/Sprat/doc/pool_tree_rename_HiC_DeDup_v2.png", height = 1000, width = 1000)
plot.phylo(Sprat_tree, type="unrooted", lab4ut = "axial", cex = 3, edge.width = 2.5)
add.scale.bar(length = NULL, ask = FALSE,
              lwd = 1, lcol = "black")
dev.off()

#Figure 2 - DAF histogram
hist(sprat_lo_GR$oce_v_brack[sprat_site_filter])
quantile(x = sprat_lo_GR$oce_v_brack[sprat_site_filter], probs = 1-c(0.1, 0.02, 0.003))
pdf(file = "~/Projects/Sprat/doc/OvB_DAF_hist_HiC_DeDup_v2.pdf", width = 10, height = 7)
hist(sprat_lo_GR$oce_v_brack[sprat_site_filter])
dev.off()

#Figure 2.5 - Herring lift-over (GW)
#Sprat_lo_manhattan_GW(sprat_data = sprat_post_lo_GR$lo_GRs, scaff_list = paste0("chr", c(1:26)), daf_col = "oce_v_brack", png_file = "~/Projects/Sprat/doc/figure_sketches/oce_v_brack_Ch_lo_GW_replot.png", site_filt_vec = sprat_site_filter_m1, backbone = "herring")
#Sprat_lo_manhattan_GW(sprat_data = sprat_post_lo_GR$lo_GRs, scaff_list = paste0("chr", c(1:26)), daf_col = "oce_v_fjord", png_file = "~/Projects/Sprat/doc/figure_sketches/oce_v_fjord_Ch_lo_GW.png", site_filt_vec = sprat_site_filter_m1, backbone = "herring")
#Sprat_lo_manhattan_GW(sprat_data = sprat_post_lo_GR$lo_GRs, scaff_list = paste0("chr", c(1:26)), daf_col = "fjord_v_brack", png_file = "~/Projects/Sprat/doc/figure_sketches/fjord_v_brack_Ch_lo_GW.png", site_filt_vec = sprat_site_filter_m1, backbone = "herring")
Sprat_lo_manhattan_GW(sprat_data = sprat_lo_GR, scaff_list = paste0("chr", c(1:26)), daf_col = "oce_v_brack", png_file = "~/Projects/Sprat/doc/figure_sketches/oce_v_brack_Ch_lo_GW_HiC_DeDup_v2.png", site_filt_vec = sprat_site_filter, backbone = "herring")
Sprat_lo_manhattan_GW(sprat_data = sprat_lo_GR, scaff_list = paste0("chr", c(1:26)), daf_col = "oce_v_fjord", png_file = "~/Projects/Sprat/doc/figure_sketches/oce_v_fjord_Ch_lo_GW_HiC_DeDup_v2.png", site_filt_vec = sprat_site_filter, backbone = "herring")
Sprat_lo_manhattan_GW(sprat_data = sprat_lo_GR, scaff_list = paste0("chr", c(1:26)), daf_col = "fjord_v_brack", png_file = "~/Projects/Sprat/doc/figure_sketches/fjord_v_brack_Ch_lo_GW_HiC_DeDup_v2.png", site_filt_vec = sprat_site_filter, backbone = "herring")




sprat_lo_GR$oce_DDAF <- sprat_lo_GR$oce_v_brack - sprat_lo_GR$oce_v_fjord
Sprat_lo_manhattan_GW(sprat_data = sprat_lo_GR, scaff_list = paste0("chr", c(1:26)), daf_col = "oce_DDAF", png_file = "~/Projects/Sprat/doc/figure_sketches/oce_DDAF_Ch_lo_GW_HiC_DeDup_v2.png", site_filt_vec = sprat_site_filter, backbone = "herring", y_lim = c(-1,1))


#Figure 3 - Herring lift-over (inversions)

#Sprat_lo_manhattan(sprat_data = Sprat_pool_freq_lo, daf_col = "oce_v_brack", sel_reg_GR = mixed_reg_GR, plot_prefix = "~/Projects/Sprat/doc/Oce_v_Brack/Oce_v_Brack_mixed_reg_", site_filt_vec = sprat_site_filter_lo )
#Sprat_lo_manhattan_Ch(sprat_data = Sprat_pool_freq_lo, scaff_list = paste0("chr", 1:26), daf_col = "oce_v_brack", sel_reg_GR = sal_p20_regions_GR, plot_prefix = "~/Projects/Sprat/doc/draft_manhattan/Oce_v_Brack_replot_", site_filt_vec = sprat_site_filter_lo )
#Sprat_lo_manhattan_v3(sprat_data = sprat_post_lo_GR$lo_GRs, scaff_list = names(Sprat_pin_curated)[1], daf_col = "Oce_v_Baltic", sel_reg_GR = sal_p20_regions_GR, plot_prefix = "~/Projects/Sprat/doc/figure_sketches/Oce_v_Baltic_m1_replot_", site_filt_vec = sprat_site_filter_m1 )
#Sprat_lo_manhattan_v3(sprat_data = sprat_post_lo_GR$lo_GRs, scaff_list = paste0("chr", c(4,5,7,13, 18, 21)), daf_col = "oce_v_brack", sel_reg_GR = sal_p20_regions_GR, plot_prefix = "~/Projects/Sprat/doc/figure_sketches/oce_v_brack_replot2_Ch_", site_filt_vec = sprat_site_filter_m1, backbone = "herring")
#Sprat_lo_manhattan_v3(sprat_data = sprat_lo_GR, scaff_list = paste0("chr", c(4,5,7,13, 18, 21)), daf_col = "oce_v_brack", sel_reg_GR = sal_p20_regions_GR[1], plot_prefix = "~/Projects/Sprat/doc/figure_sketches/oce_v_brack_HiC_DeDup_v2_Ch_", site_filt_vec = sprat_site_filter, backbone = "herring")
Sprat_lo_manhattan_v3(sprat_data = sprat_lo_GR, scaff_list = paste0("chr", c(4,5,7,13, 18, 21)), daf_col = "oce_v_brack", sel_reg_GR = sal_p20_regions_GR[1], plot_prefix = "~/Projects/Sprat/doc/figure_sketches/replot_OvB_HiC_DeDup_v2_Ch_", site_filt_vec = sprat_site_filter, backbone = "herring", th_lwd = 1.5, point_cex = 0.9)


#PGA_scaffold_1129__24_contigs__length_18973662


#Figure 4(i) - Comparing the contrasts
png(file = "~/Projects/Sprat/doc/figure_sketches/OvB_vs_OvF_scatter.png", width = 1000, height = 1000)
plot(x = sprat_lo_GR$oce_v_brack[sprat_site_filter], y = sprat_lo_GR$oce_v_fjord[sprat_site_filter], pch = 16, col = "black", cex = 0.5)
dev.off()


fjord_assc_filter <- sprat_site_filter & sprat_lo_GR$oce_v_brack < 0.1 & sprat_lo_GR$oce_v_fjord > 0.4 & sprat_lo_GR$CHROM == "PGA_scaffold_1113__35_contigs__length_27024755" & sprat_lo_GR$POS < 6.5e6

simple_sprat_pool_hm_v4(snp_data = sprat_lo_GR, reg_filter = fjord_assc_filter, sample_info_df = sprat_sample_info, pdf_file =  "~/Projects/Sprat/doc/figure_sketches/PGA_scaffold_1113_DeDup_v2_fjord_HM.pdf", freq_ref_pools = sprat_sample_info$clean_name[15:18], ordered_plot_pools = pools_to_plot)



#Figure 4 (ii) - putative inversion
#tmp_filter <- Sprat_pool_freq_lo$herring_v2.0.2_seqnames == "chr18" & sprat_site_filter_lo & Sprat_pool_freq_lo$oce_v_brack > 0.7
#simple_sprat_pool_hm(snp_data = Sprat_pool_freq_lo, reg_filter = tmp_filter, pdf_file =  "~/Projects/Sprat/doc/reg_heatmaps/Chr18_inv.pdf")

chr18_inv_filter <- sprat_lo_GR$CHROM == "PGA_scaffold_40__17_contigs__length_6491151" & sprat_site_filter & sprat_lo_GR$oce_v_brack > 0.7 & sprat_lo_GR$POS > 1.5e6
#simple_sprat_pool_hm_v3(snp_data = sprat_post_lo_GR$lo_GRs, reg_filter = chr18_inv_filter, sample_info_df = sprat_sample_info, pdf_file =  "~/Projects/Sprat/doc/reg_heatmaps_m1/PGA_scaffold_40_rename.pdf")
simple_sprat_pool_hm_v4(snp_data = sprat_lo_GR, reg_filter = chr18_inv_filter, sample_info_df = sprat_sample_info, pdf_file =  "~/Projects/Sprat/doc/figure_sketches/PGA_scaffold_40_DeDup_v2_order_HM.pdf", freq_ref_pools = sprat_sample_info$clean_name[15:18], ordered_plot_pools = pools_to_plot)





#Figure 5

#Sprat_lo_manhattan_v3(sprat_data = sprat_post_lo_GR$lo_GRs, scaff_list = "chr12", daf_col = "Oce_v_Baltic", plot_prefix = "~/Projects/Sprat/doc/figure_sketches/Oce_v_Baltic_replot_Ch_", site_filt_vec = sprat_site_filter_m1, backbone = "herring", sel_reg_GR = sal_p20_regions_GR) 
#Sprat_lo_manhattan_v3(sprat_data = sprat_post_lo_GR$lo_GRs, scaff_list = "chr19", daf_col = "Oce_v_Baltic", plot_prefix = "~/Projects/Sprat/doc/figure_sketches/Oce_v_Baltic_replot_Ch_", site_filt_vec = sprat_site_filter_m1, backbone = "herring", sel_reg_GR = sal_p20_regions_GR)
Sprat_lo_manhattan_v3(sprat_data = sprat_lo_GR, scaff_list = "chr12", daf_col = "Oce_v_Baltic", plot_prefix = "~/Projects/Sprat/doc/figure_sketches/Oce_v_Baltic_DeDup_v2_Ch_", site_filt_vec = sprat_site_filter, backbone = "herring", sel_reg_GR = sal_p20_regions_GR) 
Sprat_lo_manhattan_v3(sprat_data = sprat_lo_GR, scaff_list = "chr19", daf_col = "Oce_v_Baltic", plot_prefix = "~/Projects/Sprat/doc/figure_sketches/Oce_v_Baltic_DeDup_v2_Ch_", site_filt_vec = sprat_site_filter, backbone = "herring", sel_reg_GR = sal_p20_regions_GR)


#chr12_ch_filter <- Sprat_pool_freq_lo$herring_v2.0.2_seqnames == "chr12" & sprat_site_filter_lo & Sprat_pool_freq_lo$oce_v_brack > 0.4 & Sprat_pool_freq_lo$SNP_HiC_pos > 2e6 & Sprat_pool_freq_lo$SNP_HiC_pos < 3e6
#simple_sprat_pool_hm_Ch(snp_data = Sprat_pool_freq_lo, reg_filter = chr12_ch_filter, sample_info_df = sprat_sample_info, pdf_file =  "~/Projects/Sprat/doc/reg_heatmaps/Chr12_2.3MB_rename.pdf", freq_ref_pools = sprat_sample_info$clean_name[15:18], ordered_plot_pools = pools_to_plot)

#chr12_spsp_filter <- 1:length(sprat_post_lo_GR$lo_GRs$Ch_v2.0.2_chr) %in% which(sprat_post_lo_GR$lo_GRs$Ch_v2.0.2_chr == "chr12" & sprat_site_filter_m1 & sprat_post_lo_GR$lo_GRs$oce_v_brack > 0.4 & sprat_post_lo_GR$lo_GRs$Ch_v2.0.2_pos > 2e6 & sprat_post_lo_GR$lo_GRs$Ch_v2.0.2_pos < 3e6)
#simple_sprat_pool_hm_v4(snp_data = sprat_post_lo_GR$lo_GRs, reg_filter = chr12_spsp_filter, sample_info_df = sprat_sample_info, pdf_file =  "~/Projects/Sprat/doc/reg_heatmaps_m1/Chr12_2.3MB_lo_join.pdf", freq_ref_pools = sprat_sample_info$clean_name[15:18], ordered_plot_pools = pools_to_plot)

chr12_spsp_filter <- 1:length(sprat_lo_GR$Ch_v2.0.2_chr) %in% which(sprat_lo_GR$Ch_v2.0.2_chr == "chr12" & sprat_site_filter & sprat_lo_GR$Oce_v_Baltic > 0.7 & sprat_lo_GR$Ch_v2.0.2_pos > 2e6 & sprat_lo_GR$Ch_v2.0.2_pos < 3e6)
simple_sprat_pool_hm_v4(snp_data = sprat_lo_GR, reg_filter = chr12_spsp_filter, sample_info_df = sprat_sample_info, pdf_file =  "~/Projects/Sprat/doc/figure_sketches/Chr12_2.3MB_DeDup_v2.pdf", freq_ref_pools = sprat_sample_info$clean_name[15:18], ordered_plot_pools = pools_to_plot)



#chr12_snp_idx <- Sprat_pool_freq_lo$herring_v2.0.2_seqnames == "chr12" & Sprat_pool_freq_lo$SNP_HiC_pos == 2273814
#Sprat_pool_freq_lo[chr12_snp_idx,1:30]
#IPA_pri[["ctg.000058F"]][1125866 + -40:40]

#chr19_ch_filter <- Sprat_pool_freq_lo$herring_v2.0.2_seqnames == "chr19" & sprat_site_filter_lo & Sprat_pool_freq_lo$oce_v_brack > 0.4 & Sprat_pool_freq_lo$SNP_HiC_pos > 6.3e6 & Sprat_pool_freq_lo$SNP_HiC_pos < 6.5e6
#simple_sprat_pool_hm_Ch(snp_data = Sprat_pool_freq_lo, reg_filter = chr19_ch_filter, sample_info_df = sprat_sample_info, pdf_file =  "~/Projects/Sprat/doc/reg_heatmaps/Chr19_6.3MB_rename.pdf", freq_ref_pools = sprat_sample_info$clean_name[15:18], ordered_plot_pools = pools_to_plot)

#chr19_spsp_filter <- 1:length(sprat_post_lo_GR$lo_GRs$Ch_v2.0.2_chr) %in% which(sprat_post_lo_GR$lo_GRs$Ch_v2.0.2_chr == "chr19" & sprat_site_filter_m1 & sprat_post_lo_GR$lo_GRs$oce_v_brack > 0.4 & sprat_post_lo_GR$lo_GRs$Ch_v2.0.2_pos >  6.3e6 & sprat_post_lo_GR$lo_GRs$Ch_v2.0.2_pos < 6.5e6)
#simple_sprat_pool_hm_v4(snp_data = sprat_post_lo_GR$lo_GRs, reg_filter = chr19_spsp_filter, sample_info_df = sprat_sample_info, pdf_file =  "~/Projects/Sprat/doc/reg_heatmaps_m1/Chr19_6.3MB_lo_join.pdf", freq_ref_pools = sprat_sample_info$clean_name[15:18], ordered_plot_pools = pools_to_plot)

chr19_spsp_filter <- 1:length(sprat_lo_GR$Ch_v2.0.2_chr) %in% which(sprat_lo_GR$Ch_v2.0.2_chr == "chr19" & sprat_site_filter & sprat_lo_GR$oce_v_brack > 0.4 & sprat_lo_GR$Ch_v2.0.2_pos >  6.3e6 & sprat_lo_GR$Ch_v2.0.2_pos < 6.5e6)
simple_sprat_pool_hm_v4(snp_data = sprat_lo_GR, reg_filter = chr19_spsp_filter, sample_info_df = sprat_sample_info, pdf_file =  "~/Projects/Sprat/doc/figure_sketches/Chr19_6.3MB_DeDup_v2_HM.pdf", freq_ref_pools = sprat_sample_info$clean_name[15:18], ordered_plot_pools = pools_to_plot)


#chr19_snp_idx <- Sprat_pool_freq_lo$herring_v2.0.2_seqnames == "chr19" & Sprat_pool_freq_lo$SNP_HiC_pos == 6380297
#Sprat_pool_freq_lo[chr19_snp_idx,1:30]
#IPA_pri[["ctg.001310F"]][114012 + -20:20]



####



#Fig 2b - tree after excluding inversions
#Needs to be re-done
#all inversions
#paste0("chr", c(4,5,7,13, 18, 21))
#chr4: 27 MB to Chr end
#chr5: 11 to 13 MB
#chr7: 23 MB to Chr end
#ch13: 8 to 10 MB
#chr18: 13 to 20 MB
#chr21: 2 to 5 MB



c4_inv_filter <- 1:length(sprat_lo_GR$Ch_v2.0.2_chr) %in% which(sprat_lo_GR$Ch_v2.0.2_chr == "chr4" & sprat_lo_GR$Ch_v2.0.2_pos > 27e6 & sprat_site_filter)
c5_inv_filter <- 1:length(sprat_lo_GR$Ch_v2.0.2_chr) %in% which(sprat_lo_GR$Ch_v2.0.2_chr == "chr5" & sprat_lo_GR$Ch_v2.0.2_pos > 11e6 & sprat_lo_GR$Ch_v2.0.2_pos < 13e6 &sprat_site_filter)
c7_inv_filter <- 1:length(sprat_lo_GR$Ch_v2.0.2_chr) %in% which(sprat_lo_GR$Ch_v2.0.2_chr == "chr7" & sprat_lo_GR$Ch_v2.0.2_pos > 23e6 & sprat_site_filter)
c13_inv_filter <- 1:length(sprat_lo_GR$Ch_v2.0.2_chr) %in% which(sprat_lo_GR$Ch_v2.0.2_chr == "chr13" & sprat_lo_GR$Ch_v2.0.2_pos > 8e6 & sprat_lo_GR$Ch_v2.0.2_pos < 10e6 &sprat_site_filter)
c18_inv_filter <- 1:length(sprat_lo_GR$Ch_v2.0.2_chr) %in% which(sprat_lo_GR$Ch_v2.0.2_chr == "chr18" & sprat_lo_GR$Ch_v2.0.2_pos > 13e6 & sprat_lo_GR$Ch_v2.0.2_pos < 20e6 &sprat_site_filter) 
c21_inv_filter <- 1:length(sprat_lo_GR$Ch_v2.0.2_chr) %in% which(sprat_lo_GR$Ch_v2.0.2_chr == "chr21" & sprat_lo_GR$Ch_v2.0.2_pos > 2e6 & sprat_lo_GR$Ch_v2.0.2_pos < 5e6 &sprat_site_filter)

putative_inv_GR <- GRanges(seqnames = c("chr4","chr4", "chr5", "chr7", "chr13", "chr18", "chr21"), ranges = IRanges(start = c(27e6, 31e6, 11e6, 23e6, 8e6, 13e6, 2e6), end = c(29e6, Ch_v2.0.2_sizes$size[Ch_v2.0.2_sizes$name == "chr4"],13e6, 28e6 ,10e6, 20e6, 5e6)))

combined_inv_filter <- c4_inv_filter | c5_inv_filter | c7_inv_filter | c13_inv_filter | c18_inv_filter | c21_inv_filter
high_daf_filter <- 1:length(sprat_lo_GR$Ch_v2.0.2_chr) %in% which((sprat_lo_GR$oce_v_brack > 0.5 | sprat_lo_GR$oce_v_fjord > 0.5 | sprat_lo_GR$fjord_v_brack > 0.5) & sprat_site_filter)
no_inv_no_diff_filter <- !(high_daf_filter | combined_inv_filter) & sprat_site_filter

no_inv_mat  <- as.matrix(sprat_lo_GR@elementMetadata[no_inv_no_diff_filter, grep("_Afreq", names(sprat_lo_GR@elementMetadata))])
f3 <- Vectorize(FUN = function(X,Y, data_mat){sum(abs(data_mat[,X] - data_mat[,Y]), na.rm = T)/sum(!is.na(data_mat[,X]) & !is.na(data_mat[,Y]))}, vectorize.args = c("X", "Y"))
no_inv_dist <- outer(FUN ="f3", X=1:dim(no_inv_mat)[2], Y=1:dim(no_inv_mat)[2], data_mat = no_inv_mat)
rownames(no_inv_dist) <- grep("_Afreq", names(sprat_lo_GR@elementMetadata), value = T)
colnames(no_inv_dist) <- grep("_Afreq", names(sprat_lo_GR@elementMetadata), value = T)
diag(no_inv_dist) <- NA
no_inv_dist <- no_inv_dist/sum(no_inv_no_diff_filter)
rm(no_inv_mat)

save(no_inv_dist, file = "~/Projects/Sprat/data/tree/no_inv_dist_DeDup_v2.RData") #- corrected
#load(file = "~/Projects/Sprat/data/tree/no_inv_dist.RData") - still contians the inversion regions

no_inv_tree <- bionj(as.dist(no_inv_dist))
no_inv_tree_IDs <- sub("_Afreq", "", grep("_Afreq", no_inv_tree$tip.label, value = T))
no_inv_tree_names <- sprat_sample_info$clean_name[match(no_inv_tree_IDs, sprat_sample_info$POOL)]
no_inv_tree$tip.label <- no_inv_tree_names


pdf(file = "~/Projects/Sprat/doc/figure_sketches/no_inv_no_diff_DeDup_v2_tree.pdf")
plot.phylo(no_inv_tree, type = "unrooted", lab4ut = "axial", no.margin = T)
add.scale.bar(length = NULL, ask = FALSE, lwd = 1, lcol = "black")
dev.off()

#Fig 2c - tree using only divergent regions
#load(file = "~/Projects/Sprat/data/tree/OvB_Reg_dist.RData") - also affected by issues above

inv_and_diff_filter <- (high_daf_filter | combined_inv_filter) & sprat_site_filter

inv_and_diff_mat  <- as.matrix(sprat_lo_GR@elementMetadata[inv_and_diff_filter, grep("_Afreq", names(sprat_lo_GR@elementMetadata))])
inv_and_diff_dist <- outer(FUN ="f3", X=1:dim(inv_and_diff_mat)[2], Y=1:dim(inv_and_diff_mat)[2], data_mat = inv_and_diff_mat)
rownames(inv_and_diff_dist) <- grep("_Afreq", names(sprat_lo_GR@elementMetadata), value = T)
colnames(inv_and_diff_dist) <- grep("_Afreq", names(sprat_lo_GR@elementMetadata), value = T)
diag(inv_and_diff_dist) <- NA
#inv_and_diff_dist <- inv_and_diff_dist/sum(inv_and_diff_filter)

save(inv_and_diff_dist, file = "~/Projects/Sprat/data/tree/inv_and_diff_dist_DeDup_v2.RData") #- corrected

inv_tree <- bionj(as.dist(inv_and_diff_dist))
inv_tree_IDs <- sub("_Afreq", "", grep("_Afreq", inv_tree$tip.label, value = T))
inv_tree_names <- sprat_sample_info$clean_name[match(inv_tree_IDs, sprat_sample_info$POOL)]
inv_tree$tip.label <- inv_tree_names


pdf(file = "~/Projects/Sprat/doc/figure_sketches/inv_high_diff_DeDup_v2_tree.pdf")
plot.phylo(inv_tree, type = "unrooted", lab4ut = "axial", no.margin = T)
add.scale.bar(length = NULL, ask = FALSE, lwd = 1, lcol = "black")
dev.off()


#diff_mat  <- as.matrix(sprat_lo_GR@elementMetadata[high_daf_filter, grep("_Afreq", names(sprat_lo_GR@elementMetadata))])
#diff_dist <- outer(FUN ="f3", X=1:dim(diff_mat)[2], Y=1:dim(diff_mat)[2], data_mat = diff_mat)
#rownames(diff_dist) <- grep("_Afreq", names(sprat_lo_GR@elementMetadata), value = T)
#colnames(diff_dist) <- grep("_Afreq", names(sprat_lo_GR@elementMetadata), value = T)
#diag(diff_dist) <- NA
#diff_dist <- diff_dist/sum(high_daf_filter)

#Supp Fig 3 - inversion regions
#Chr 4
c4_inv_diff_filter <- c4_inv_filter & sprat_lo_GR$oce_v_brack > 0.7
simple_sprat_pool_hm_v4(snp_data = sprat_lo_GR, reg_filter = c4_inv_diff_filter, sample_info_df = sprat_sample_info, pdf_file =  "~/Projects/Sprat/doc/figure_sketches/chr4_inv_DeDup_v2_HM.pdf", freq_ref_pools = sprat_sample_info$clean_name[15:18], ordered_plot_pools = pools_to_plot)

#Chr 5
c5_inv_diff_filter <- c5_inv_filter & sprat_lo_GR$oce_v_brack > 0.6
simple_sprat_pool_hm_v4(snp_data = sprat_lo_GR, reg_filter = c5_inv_diff_filter, sample_info_df = sprat_sample_info, pdf_file =  "~/Projects/Sprat/doc/figure_sketches/chr5_inv_DeDup_v2_HM.pdf", freq_ref_pools = sprat_sample_info$clean_name[15:18], ordered_plot_pools = pools_to_plot)

#Chr 7
c7_inv_diff_filter <- c7_inv_filter & sprat_lo_GR$oce_v_brack > 0.7
simple_sprat_pool_hm_v4(snp_data = sprat_lo_GR, reg_filter = c7_inv_diff_filter, sample_info_df = sprat_sample_info, pdf_file =  "~/Projects/Sprat/doc/figure_sketches/chr7_inv_DeDup_v2_HM.pdf", freq_ref_pools = sprat_sample_info$clean_name[15:18], ordered_plot_pools = pools_to_plot)

#Chr 13
c13_inv_diff_filter <- c13_inv_filter & sprat_lo_GR$oce_v_brack > 0.7
simple_sprat_pool_hm_v4(snp_data = sprat_lo_GR, reg_filter = c13_inv_diff_filter, sample_info_df = sprat_sample_info, pdf_file =  "~/Projects/Sprat/doc/figure_sketches/chr13_inv_DeDup_v2_HM.pdf", freq_ref_pools = sprat_sample_info$clean_name[15:18], ordered_plot_pools = pools_to_plot)

#Chr 21
c21_inv_diff_filter <- c21_inv_filter & sprat_lo_GR$oce_v_brack > 0.7
simple_sprat_pool_hm_v4(snp_data = sprat_lo_GR, reg_filter = c21_inv_diff_filter, sample_info_df = sprat_sample_info, pdf_file =  "~/Projects/Sprat/doc/figure_sketches/chr21_inv_DeDup_v2_HM.pdf", freq_ref_pools = sprat_sample_info$clean_name[15:18], ordered_plot_pools = pools_to_plot)

#require(ape)
#OvB_reg_tree <- bionj(as.dist(inv_and_diff_dist))
#OvB_reg_tree_IDs <- sub("_Afreq", "", grep("_Afreq", OvB_reg_tree$tip.label, value = T))
#OvB_reg_tree_names <- sprat_sample_info$clean_name[match(OvB_reg_tree_IDs, sprat_sample_info$POOL)]
#OvB_reg_tree_names <- paste0("___", OvB_reg_tree_names, "___")
#OvB_reg_tree$tip.label <- OvB_reg_tree_names

#pdf(file = "~/Projects/Sprat/doc/inv_and_diff_tree.pdf")
#plot.phylo(OvB_reg_tree, type = "unrooted", lab4ut = "axial", no.margin = T, cex = 0.9)
#dev.off()

#OvB_diff_tree <- bionj(as.dist(diff_dist))
#OvB_diff_tree_IDs <- sub("_Afreq", "", grep("_Afreq", OvB_diff_tree$tip.label, value = T))
#OvB_diff_tree_names <- sprat_sample_info$clean_name[match(OvB_diff_tree_IDs, sprat_sample_info$POOL)]
#OvB_diff_tree_names <- paste0("___", OvB_diff_tree_names, "___")
#OvB_diff_tree$tip.label <- OvB_diff_tree_names

#pdf(file = "~/Projects/Sprat/doc/OvB_diff_tree.pdf")
#plot.phylo(OvB_diff_tree, type = "unrooted", lab4ut = "axial", no.margin = T, cex = 0.9)
#add.scale.bar(length = NULL, ask = FALSE, lwd = 1, lcol = "black")
#dev.off()




#Supp Fig 4 - Scaff 1114 heatmap
#simple_sprat_pool_hm_v3(snp_data = sprat_post_lo_GR$lo_GRs, reg_filter = tmp_filter_1114_1, sample_info_df = sprat_sample_info, pdf_file = "~/Projects/Sprat/doc/reg_heatmaps_m1/PGA_scaffold_1114_6MB_rename.pdf")
tmp_filter_1114_1 <- sprat_lo_GR$CHROM == "PGA_scaffold_1114__20_contigs__length_23438440" & sprat_site_filter & sprat_lo_GR$oce_v_brack > 0.35 & sprat_lo_GR$POS > 5.0e6 & sprat_lo_GR$POS < 7.0e6
simple_sprat_pool_hm_v4(snp_data = sprat_lo_GR, reg_filter = tmp_filter_1114_1, sample_info_df = sprat_sample_info, freq_ref_pools = sprat_sample_info$clean_name[15:18], ordered_plot_pools = pools_to_plot, pdf_file = "~/Projects/Sprat/doc/figure_sketches//PGA_scaffold_1114_6MB_DeDup_v2.pdf")

#Supp Fig 5 - Scaff 4 heatmap
#simple_sprat_pool_hm_v3(snp_data = sprat_post_lo_GR$lo_GRs, reg_filter = tmp_filter_4, sample_info_df = sprat_sample_info, pdf_file =  "~/Projects/Sprat/doc/reg_heatmaps_m1/PGA_scaffold_4_19MB_rename.pdf")
#simple_sprat_pool_hm_v4(snp_data = sprat_post_lo_GR$lo_GRs, reg_filter = tmp_filter_4, sample_info_df = sprat_sample_info, freq_ref_pools = sprat_sample_info$clean_name[15:18],  pdf_file =  "~/Projects/Sprat/doc/reg_heatmaps_m1/PGA_scaffold_4_19MB_flip.pdf")
tmp_filter_4 <- sprat_lo_GR$CHROM == "PGA_scaffold_4__24_contigs__length_24387585" & sprat_site_filter & sprat_lo_GR$oce_v_brack > 0.2 & sprat_lo_GR$POS < 20.0e6 & sprat_lo_GR$POS > 18.5e6
simple_sprat_pool_hm_v4(snp_data = sprat_lo_GR, reg_filter = tmp_filter_4, sample_info_df = sprat_sample_info, freq_ref_pools = sprat_sample_info$clean_name[15:18], ordered_plot_pools = pools_to_plot, pdf_file =  "~/Projects/Sprat/doc/figure_sketches/PGA_scaffold_4_19MB_DeDup_v2.pdf")



##Possible extra supporting figures
##Ove v Baltic contrast - all chromsomes
#Sprat_lo_manhattan_v3(sprat_data = sprat_post_lo_GR$lo_GRs, scaff_list = paste0("chr", 1:26), daf_col = "Oce_v_Baltic", sel_reg_GR = sal_p20_regions_GR, plot_prefix = "~/Projects/Sprat/doc/Oce_v_Baltic_Ch/Oce_v_Baltic_Ch_", site_filt_vec = sprat_site_filter_m1, backbone = "herring")
#filter_chr3_31MB <- 1:length(sprat_post_lo_GR$lo_GRs$Ch_v2.0.2_chr) %in% which(sprat_post_lo_GR$lo_GRs$Ch_v2.0.2_chr == "chr3" & sprat_post_lo_GR$lo_GRs$Oce_v_Baltic > 0.60 & sprat_post_lo_GR$lo_GRs$Ch_v2.0.2_pos > 3e7 & sprat_site_filter_m1)
#simple_sprat_pool_hm_v4(snp_data = sprat_post_lo_GR$lo_GRs, reg_filter = filter_chr3_31MB, sample_info_df = sprat_sample_info, freq_ref_pools = sprat_sample_info$clean_name[15:18], ordered_plot_pools = pools_to_plot, pdf_file =  "~/Projects/Sprat/doc/figure_sketches/Chr3_31.9MB_order.pdf")

Sprat_lo_manhattan_v3(sprat_data = sprat_lo_GR, scaff_list = paste0("chr", c(3, 12, 19)), daf_col = "Oce_v_Baltic", sel_reg_GR = sal_p20_regions_GR, plot_prefix = "~/Projects/Sprat/doc/figure_sketches/oce_v_balt_HiC_DeDup_v2_Ch_", site_filt_vec = sprat_site_filter, backbone = "herring")
filter_chr3_31MB <- 1:length(sprat_lo_GR$Ch_v2.0.2_chr) %in% which(sprat_lo_GR$Ch_v2.0.2_chr == "chr3" & sprat_lo_GR$Oce_v_Baltic > 0.60 & sprat_lo_GR$Ch_v2.0.2_pos > 3e7 & sprat_site_filter)

simple_sprat_pool_hm_v4(snp_data = sprat_lo_GR, reg_filter = filter_chr3_31MB, sample_info_df = sprat_sample_info, freq_ref_pools = sprat_sample_info$clean_name[15:18], ordered_plot_pools = pools_to_plot, pdf_file =  "~/Projects/Sprat/doc/figure_sketches/Chr3_31.9MB_DeDup_v2.pdf")




#"Fjord region" on Chr 12
#chr12_fjord_filter <- (1:length(sprat_post_lo_GR$lo_GRs$oce_v_fjord)) %in% which(sprat_post_lo_GR$lo_GRs$oce_v_fjord > 0.5 & sprat_post_lo_GR$lo_GRs$oce_v_brack < 0.2 & sprat_site_filter_m1 & sprat_post_lo_GR$lo_GRs$Ch_v2.0.2_chr == "chr12" & sprat_post_lo_GR$lo_GRs$Ch_v2.0.2_pos < 1e6)
#simple_sprat_pool_hm_v4(snp_data = sprat_post_lo_GR$lo_GRs, reg_filter = chr12_fjord_filter, sample_info_df = sprat_sample_info, freq_ref_pools = sprat_sample_info$clean_name[15:18], ordered_plot_pools = pools_to_plot, pdf_file =  "~/Projects/Sprat/doc/reg_heatmaps_m1/Chr12_1Mb_flip_order.pdf")
chr12_fjord_filter <- (1:length(sprat_lo_GR$oce_v_fjord)) %in% which(sprat_lo_GR$oce_v_fjord > 0.5 & sprat_lo_GR$oce_v_brack < 0.2 & sprat_site_filter & sprat_lo_GR$Ch_v2.0.2_chr == "chr12" & sprat_lo_GR$Ch_v2.0.2_pos < 1e6)
simple_sprat_pool_hm_v4(snp_data = sprat_lo_GR, reg_filter = chr12_fjord_filter, sample_info_df = sprat_sample_info, freq_ref_pools = sprat_sample_info$clean_name[15:18], ordered_plot_pools = pools_to_plot, pdf_file =  "~/Projects/Sprat/doc/figure_sketches/Chr12_1Mb_DeDup_v2.pdf")



#Independent regions of selection

#Cumulative positions
#scaff_list <- paste0("chr", c(1:26))
scaff_list <- names(Sprat_DeDup_v2)[1:23]
site_filt_vec <- sprat_site_filter
#tmp_chr_filter <- (1:length(sprat_lo_GR$Ch_v2.0.2_chr) %in% which(sprat_lo_GR$Ch_v2.0.2_chr %in% scaff_list)) & site_filt_vec
tmp_chr_filter <- (1:length(sprat_lo_GR$CHROM) %in% which(sprat_lo_GR$CHROM %in% scaff_list)) & site_filt_vec

#scaff_maxes <- aggregate(sprat_lo_GR$Ch_v2.0.2_pos[tmp_chr_filter] ~ sprat_lo_GR$Ch_v2.0.2_chr[tmp_chr_filter], FUN = "max")
scaff_maxes <- aggregate(sprat_lo_GR$POS[tmp_chr_filter] ~ sprat_lo_GR$CHROM[tmp_chr_filter], FUN = "max")
scaff_maxes <- scaff_maxes[match(scaff_list, scaff_maxes[,1]),]

scaff_maxes$offset <- c(0,cumsum(scaff_maxes[,2])[-dim(scaff_maxes)[1]])
scaff_maxes$col <- rep(c("grey30", "grey70"), times = ceiling(length(scaff_maxes$offset)/2))[1:length(scaff_maxes$offset)]

#OvB_thresh_val <- qnorm(p = 1-(0.05/sum(sprat_site_filter_m1)), mean = mean(sprat_post_lo_GR$lo_GRs@elementMetadata[sprat_site_filter_m1, "oce_v_brack"]), sd  = sd(sprat_post_lo_GR$lo_GRs@elementMetadata[sprat_site_filter_m1, "oce_v_brack"]))
#OvB_reg_filter <- (1:length(sprat_post_lo_GR$lo_GRs$oce_v_brack)) %in% which(sprat_post_lo_GR$lo_GRs$oce_v_brack > OvB_thresh_val & sprat_site_filter_m1 & !is.na(sprat_post_lo_GR$lo_GRs$Ch_v2.0.2_chr))


#OvB_reg_GR_raw <- GRanges(seqnames = sprat_post_lo_GR$lo_GRs$Ch_v2.0.2_chr[OvB_reg_filter], IRanges(start = sprat_post_lo_GR$lo_GRs$Ch_v2.0.2_pos[OvB_reg_filter], end = sprat_post_lo_GR$lo_GRs$Ch_v2.0.2_pos[OvB_reg_filter]))
#OvB_regions_GR <- reduce(OvB_reg_GR_raw, min.gapwidth = 5e5)
#OvB_regions_GR$cumulative_start <- OvB_regions_GR@ranges@start + scaff_maxes$offset[match(as.character(OvB_regions_GR@seqnames), scaff_maxes[,1])]

fuse_thresh <- 5e5
OvB_thresh_val <- qnorm(p = 1-(0.05/sum(sprat_site_filter)), mean = mean(sprat_lo_GR@elementMetadata[sprat_site_filter, "oce_v_brack"]), sd  = sd(sprat_lo_GR@elementMetadata[sprat_site_filter, "oce_v_brack"]))
OvB_reg_filter <- (1:length(sprat_lo_GR$oce_v_brack)) %in% which(sprat_lo_GR$oce_v_brack > OvB_thresh_val & sprat_site_filter)

OvB_reg_DeDup_v2_GR_raw <- GRanges(seqnames = sprat_lo_GR$CHROM[OvB_reg_filter], IRanges(start = sprat_lo_GR$POS[OvB_reg_filter], end = sprat_lo_GR$POS[OvB_reg_filter]))
OvB_reg_DeDup_v2_GR <- reduce(OvB_reg_DeDup_v2_GR_raw, min.gapwidth = fuse_thresh)
OvB_reg_DeDup_v2_GR$cumulative_start <- OvB_reg_DeDup_v2_GR@ranges@start + scaff_maxes$offset[match(as.character(OvB_reg_DeDup_v2_GR@seqnames), scaff_maxes[,1])]
OvB_reg_DeDup_v2_GR$contrast <- "Oceanic_vs_Brackish" 

OvF_thresh_val <- qnorm(p = 1-(0.05/sum(sprat_site_filter)), mean = mean(sprat_lo_GR@elementMetadata[sprat_site_filter, "oce_v_fjord"]), sd  = sd(sprat_lo_GR@elementMetadata[sprat_site_filter, "oce_v_fjord"]))
OvF_reg_filter <- (1:length(sprat_lo_GR$oce_v_fjord)) %in% which(sprat_lo_GR$oce_v_fjord > OvF_thresh_val & sprat_site_filter)

OvF_reg_DeDup_v2_GR_raw <- GRanges(seqnames = sprat_lo_GR$CHROM[OvF_reg_filter], IRanges(start = sprat_lo_GR$POS[OvF_reg_filter], end = sprat_lo_GR$POS[OvF_reg_filter]))
OvF_reg_DeDup_v2_GR <- reduce(OvF_reg_DeDup_v2_GR_raw, min.gapwidth = fuse_thresh)
OvF_reg_DeDup_v2_GR$cumulative_start <- OvF_reg_DeDup_v2_GR@ranges@start + scaff_maxes$offset[match(as.character(OvF_reg_DeDup_v2_GR@seqnames), scaff_maxes[,1])]

OvF_reg_DeDup_v2_GR$contrast <- "Oceanic_vs_Coastal" 

reg_out_df <- rbind(as.data.frame(OvB_reg_DeDup_v2_GR)[,c("seqnames","start", "end", "contrast")], as.data.frame(OvF_reg_DeDup_v2_GR)[,c("seqnames","start", "end", "contrast")])
write.table(reg_out_df, file = "~/Projects/Sprat/data/sprat_selection_regions.txt", row.names = F, quote = F, sep = "\t")



#sum(grepl("chr", OvB_regions_GR@seqnames))
#hist(OvB_regions_GR@ranges@width)

hist(OvB_reg_DeDup_v2_GR@ranges@width)
hist(OvF_reg_DeDup_v2_GR@ranges@width)
length(OvB_reg_DeDup_v2_GR)
sum(OvB_reg_DeDup_v2_GR@ranges@width > 1)
sum(width(OvB_reg_DeDup_v2_GR))/1e6

length(OvF_reg_DeDup_v2_GR)
sum(OvF_reg_DeDup_v2_GR@ranges@width > 1)
sum(width(OvF_reg_DeDup_v2_GR))/1e6




#Sprat_lo_manhattan_GW(sprat_data = sprat_post_lo_GR$lo_GRs, scaff_list = paste0("chr", c(1:26)), daf_col = "oce_v_brack", png_file = "~/Projects/Sprat/doc/figure_sketches/oce_v_brack_Ch_GW_regions.png", site_filt_vec = sprat_site_filter_m1, backbone = "herring", sel_reg = OvB_regions_GR)
Sprat_DeDup_manhattan_GW(sprat_data = sprat_lo_GR, scaff_list = names(Sprat_DeDup_v2)[1:23], daf_col = "oce_v_brack", png_file = "~/Projects/Sprat/doc/figure_sketches/oce_v_brack_DeDup_v2_GW_regions.png", site_filt_vec = sprat_site_filter, backbone = "sprat", sel_reg = OvB_reg_DeDup_v2_GR)
Sprat_DeDup_manhattan_GW(sprat_data = sprat_lo_GR, scaff_list = names(Sprat_DeDup_v2)[1:23], daf_col = "oce_v_fjord", png_file = "~/Projects/Sprat/doc/figure_sketches/oce_v_fjord_DeDup_v2_GW_regions.png", site_filt_vec = sprat_site_filter, backbone = "sprat", sel_reg = OvF_reg_DeDup_v2_GR)

OvF_v_OvB <- findOverlaps(OvF_reg_DeDup_v2_GR, OvB_reg_DeDup_v2_GR)
OvF_v_OvB_pairs <- findOverlapPairs(OvF_reg_DeDup_v2_GR, OvB_reg_DeDup_v2_GR)
OvF_v_OvB_intersect <- pintersect(OvF_v_OvB_pairs)

sum(width(OvB_reg_DeDup_v2_GR[unique(OvF_v_OvB@to)]))/1e6
sum(width(OvF_v_OvB_intersect))/1e6

OvB_reg_size_by_scaff <- aggregate(width(OvB_reg_DeDup_v2_GR)~as.character(seqnames(OvB_reg_DeDup_v2_GR)), FUN = "sum")
OvB_reg_size_by_scaff <- OvB_reg_size_by_scaff[order(OvB_reg_size_by_scaff[,2], decreasing = T),]



#OvF_thresh_val <- qnorm(p = 1-(0.05/sum(sprat_site_filter_m1)), mean = mean(sprat_post_lo_GR$lo_GRs@elementMetadata[sprat_site_filter_m1, "oce_v_fjord"]), sd  = sd(sprat_post_lo_GR$lo_GRs@elementMetadata[sprat_site_filter_m1, "oce_v_fjord"]))
#OvF_reg_filter <- (1:length(sprat_post_lo_GR$lo_GRs$oce_v_fjord)) %in% which(sprat_post_lo_GR$lo_GRs$oce_v_fjord > OvF_thresh_val & sprat_site_filter_m1 & !is.na(sprat_post_lo_GR$lo_GRs$Ch_v2.0.2_chr))

#OvF_reg_GR_raw <- GRanges(seqnames = sprat_post_lo_GR$lo_GRs$Ch_v2.0.2_chr[OvF_reg_filter], IRanges(start = sprat_post_lo_GR$lo_GRs$Ch_v2.0.2_pos[OvF_reg_filter], end = sprat_post_lo_GR$lo_GRs$Ch_v2.0.2_pos[OvF_reg_filter]))
#OvF_regions_GR <- reduce(OvF_reg_GR_raw, min.gapwidth = 5e5)
#OvF_regions_GR$cumulative_start <- OvF_regions_GR@ranges@start + scaff_maxes$offset[match(as.character(OvF_regions_GR@seqnames), scaff_maxes[,1])]

#OvF_reg_m1_GR_raw <- GRanges(seqnames = sprat_post_lo_GR$lo_GRs$HiC_scaffold[OvF_reg_filter], IRanges(start = sprat_post_lo_GR$lo_GRs$HiC_pos[OvF_reg_filter], end = sprat_post_lo_GR$lo_GRs$HiC_pos[OvF_reg_filter]))
#OvF_regions_m1_GR <- reduce(OvF_reg_m1_GR_raw, min.gapwidth = 5e5)
#Sprat_lo_manhattan_GW(sprat_data = sprat_post_lo_GR$lo_GRs, scaff_list = paste0("chr", c(1:26)), daf_col = "oce_v_fjord", png_file = "~/Projects/Sprat/doc/figure_sketches/oce_v_fjord_Ch_GW_regions.png", site_filt_vec = sprat_site_filter_m1, backbone = "herring", sel_reg = OvF_regions_GR)

OvF_v_OvB <- findOverlaps(OvF_regions_GR, OvB_regions_GR)
OvF_v_OvB_pairs <- findOverlapPairs(OvF_regions_GR, OvB_regions_GR)
OvF_v_OvB_intersect <- pintersect(OvF_v_OvB_pairs)

sum(width(OvB_regions_GR[unique(OvF_v_OvB@to)]))


OvB_regions_GR[grepl("chr",OvB_regions_GR@seqnames) & OvB_regions_GR@ranges@width > 1]
sum(width(OvB_regions_GR[grepl("chr",OvB_regions_GR@seqnames) & OvB_regions_GR@ranges@width > 1]))

OvF_regions_GR[grepl("chr",OvF_regions_GR@seqnames)]
OvF_regions_GR[grepl("chr",OvF_regions_GR@seqnames) & OvF_regions_GR@ranges@width > 1]
sum(width(OvF_regions_GR[grepl("chr",OvF_regions_GR@seqnames) & OvF_regions_GR@ranges@width > 1]))

#The inversions 
OvB_regions_GR[grepl("chr(4|5|7|13|18|21)",OvB_regions_GR@seqnames) & OvB_regions_GR@ranges@width > 1 & 1:length(OvB_regions_GR@seqnames) %in% OvF_v_OvB@to]
#non-inversion overlaps
non_inv_hits <- OvF_v_OvB_intersect %outside% OvB_regions_GR[grepl("chr(4|5|7|13|18|21)",OvB_regions_GR@seqnames) & OvB_regions_GR@ranges@width > 1 & 1:length(OvB_regions_GR@seqnames) %in% OvF_v_OvB@to]
sum(width(OvF_v_OvB_intersect[non_inv_hits & grepl("chr",OvF_v_OvB_intersect@seqnames)]))


#cumulative_pos_vec <- sprat_post_lo_GR$lo_GRs$Ch_v2.0.2_pos + scaff_maxes$offset[match(sprat_post_lo_GR$lo_GRs$Ch_v2.0.2_chr, scaff_maxes[,1])]
#cumulative_col_vec <- scaff_maxes$col[match(sprat_data$Ch_v2.0.2_chr[tmp_chr_filter], scaff_maxes[,1])]

#png(filename = "~/Projects/Sprat/doc/figure_sketches/OvB_OvF_scatterplot.png", height = 1000, width = 1000)
#plot(x = sprat_lo_GR$oce_v_brack[sprat_site_filter_m1], y = sprat_lo_GR$oce_v_fjord[sprat_site_filter])
#abline(v = OvB_thresh_val, col = "red")
#abline(h = OvF_thresh_val, col = "darkorchid")
#dev.off()

#Summary statistics
summary(sprat_lo_GR$oce_v_brack[sprat_site_filter])
summary(sprat_lo_GR$oce_v_fjord[sprat_site_filter])
Ovb_OvF_cor <- cor.test(x = sprat_lo_GR$oce_v_brack[sprat_site_filter], y = sprat_lo_GR$oce_v_fjord[sprat_site_filter])
Ovb_OvF_cor$estimate^2 #r^2

Ovb_OvF_diff_cor <- cor.test(x = sprat_lo_GR$oce_v_brack[inv_and_diff_filter], y = sprat_lo_GR$oce_v_fjord[inv_and_diff_filter])
Ovb_OvF_diff_cor$estimate^2 #r^2

#Genome wide liftover summary
#Sprat_scaff_vs_Ch_chr <- aggregate(sprat_lo_GR$Ch_v2.0.2_chr[!is.na(sprat_lo_GR$Ch_v2.0.2_chr)], by = sprat_lo_GR$CHROM[!is.na(sprat_lo_GR$Ch_v2.0.2_chr)], FUN = "sum")
Sprat_scaff_vs_Ch_chr_df  <- data.frame(sprat = as.character(sprat_lo_GR$CHROM[!is.na(sprat_lo_GR$Ch_v2.0.2_chr)]), sprat_pos = as.character(sprat_lo_GR$POS[!is.na(sprat_lo_GR$Ch_v2.0.2_chr)]), herring = as.character(sprat_lo_GR$Ch_v2.0.2_chr[!is.na(sprat_lo_GR$Ch_v2.0.2_chr)]), herring_pos = as.character(sprat_lo_GR$Ch_v2.0.2_pos[!is.na(sprat_lo_GR$Ch_v2.0.2_chr)]), stringsAsFactors = F)
Sprat_scaff_vs_Ch_chr_df$sprat_pos <- as.numeric(Sprat_scaff_vs_Ch_chr_df$sprat_pos)
Sprat_scaff_vs_Ch_chr_df$herring_pos <- as.numeric(Sprat_scaff_vs_Ch_chr_df$herring_pos)

Sprat_scaff_vs_Ch_chr_df <- Sprat_scaff_vs_Ch_chr_df[grep("chr", Sprat_scaff_vs_Ch_chr_df$herring),]
Sprat_scaff_vs_Ch_chr_tab <- table(Sprat_scaff_vs_Ch_chr_df[,c("sprat", "herring")])
Sprat_scaff_vs_Ch_chr_tab <- Sprat_scaff_vs_Ch_chr_tab[rowSums(Sprat_scaff_vs_Ch_chr_tab) > 5e3,]
Sprat_scaff_vs_Ch_chr_tab[Sprat_scaff_vs_Ch_chr_tab[,"chr8"] > 5000,]

heatmap(Sprat_scaff_vs_Ch_chr_tab, Rowv = NULL, Colv = NULL, scale = "none")

pdf(file = "~/Projects/Sprat/doc/assembly_curation/Sprat_vs_Ch_hm.pdf", height = 10, width = 15)
heatmap(Sprat_scaff_vs_Ch_chr_tab[,order(as.numeric(sub("chr", "", colnames(Sprat_scaff_vs_Ch_chr_tab))))], Colv = NA, scale = "col", margins = c(4,14))
dev.off()

Sprat_vs_Ch_PGA_scaffold_1118 <- Sprat_scaff_vs_Ch_chr_df[Sprat_scaff_vs_Ch_chr_df$sprat == "PGA_scaffold_1118__66_contigs__length_56234773",]
Sprat_vs_Ch_PGA_scaffold_1118[,"col"] <- NA
Sprat_vs_Ch_PGA_scaffold_1118[Sprat_vs_Ch_PGA_scaffold_1118$herring == "chr8","col"] <- "darkorchid"
Sprat_vs_Ch_PGA_scaffold_1118[Sprat_vs_Ch_PGA_scaffold_1118$herring == "chr24","col"] <- "olivedrab"
Sprat_vs_Ch_PGA_scaffold_1118[Sprat_vs_Ch_PGA_scaffold_1118$herring == "chr26","col"] <- "firebrick"

Sprat_vs_Ch_PGA_scaffold_1118[,"lvl"] <- 1
Sprat_vs_Ch_PGA_scaffold_1118[Sprat_vs_Ch_PGA_scaffold_1118$herring == "chr24","lvl"] <- 1.3
Sprat_vs_Ch_PGA_scaffold_1118[Sprat_vs_Ch_PGA_scaffold_1118$herring == "chr26","lvl"] <- 0.7

pdf(file = "~/Projects/Sprat/doc/assembly_curation/Sprat_vs_Ch_PGA_scaffold_1118.pdf", height = 6, width = 10)
plot(x = Sprat_vs_Ch_PGA_scaffold_1118$sprat_pos, y = Sprat_vs_Ch_PGA_scaffold_1118$lvl, col = Sprat_vs_Ch_PGA_scaffold_1118$col, pch = 20, xlab = "Position", main = "PGA_scaffold_1118", ylab = "", ylim = c(0.5, 1.5))
legend(x = "topleft", legend = c("chr8", "chr24", "chr26"), col = c("darkorchid", "olivedrab", "firebrick"), pch = 20, cex = 1.5)
dev.off()

#Sprat_vs_Ch_PGA_scaffold_48 <- Sprat_scaff_vs_Ch_chr_df[Sprat_scaff_vs_Ch_chr_df$sprat == "PGA_scaffold_48__33_contigs__length_23346322",]
#Sprat_vs_Ch_PGA_scaffold_48 <- Sprat_vs_Ch_PGA_scaffold_48[Sprat_vs_Ch_PGA_scaffold_48$herring == "chr9",]

#pdf(file = "~/Projects/Sprat/doc/assembly_curation/Sprat_scaffold_48_vs_Ch_chr9.pdf", height = 6, width = 10)
#plot(x = Sprat_vs_Ch_PGA_scaffold_48$sprat_pos[1], y = 1, pch = 20, xlab = "Position", main = "PGA_scaffold_48 v chr9", ylab = "", ylim = c(0.4, 1.6), xlim = c(0, max(c(Sprat_vs_Ch_PGA_scaffold_48$sprat_pos, Sprat_vs_Ch_PGA_scaffold_48$herring_pos))), type = "n")
#segments(x0 = 0, x1 = c(max(Sprat_vs_Ch_PGA_scaffold_48$herring_pos), max(Sprat_vs_Ch_PGA_scaffold_48$sprat_pos)), y0 = c(1.5, 0.5), col = c("firebrick", "darkorchid"), lwd = 2)
#segments(x0 = Sprat_vs_Ch_PGA_scaffold_48$herring_pos, x1 = Sprat_vs_Ch_PGA_scaffold_48$sprat_pos, y0 = 1.5, y1 =  0.5, col = "slateblue", lwd = 0.3)

#legend(x = "topleft", legend = c("chr8", "chr24", "chr26"), col = c("darkorchid", "olivedrab", "firebrick"), pch = 20, cex = 1.5)
#dev.off()

Sp_scaff_v_Ch_chr_plot(ch_chr = "chr9", png_file = "~/Projects/Sprat/doc/assembly_curation/Sprat_vs_Ch_chr9.png", in_data = Sprat_scaff_vs_Ch_chr_df, n_scaffs = 2)
for(ch_chr in paste0("chr", 1:26)) Sp_scaff_v_Ch_chr_plot(ch_chr = ch_chr, png_file = paste0("~/Projects/Sprat/doc/assembly_curation/Sprat_vs_Ch_", ch_chr, ".png"), in_data = Sprat_scaff_vs_Ch_chr_df, n_scaffs = 2)

Sp_scaff_v_Ch_chr_plot(sprat_scaff = "PGA_scaffold_1118__66_contigs__length_56234773", png_file = "~/Projects/Sprat/doc/assembly_curation/Sprat_PGA_scaffold_1118_vs_Ch.png", in_data = Sprat_scaff_vs_Ch_chr_df, n_scaffs = 2, ch_sizes = Ch_v2.0.2_sizes)
Sp_scaff_v_Ch_chr_plot(sprat_scaff = "PGA_scaffold_1081__1_contigs__length_1176958", png_file = "~/Projects/Sprat/doc/assembly_curation/Sprat_PGA_scaffold_1081_vs_Ch.png", in_data = Sprat_scaff_vs_Ch_chr_df, n_scaffs = 2, ch_sizes = Ch_v2.0.2_sizes) #Has hits to only
for(sprat_scaff in rownames(Sprat_scaff_vs_Ch_chr_tab)) Sp_scaff_v_Ch_chr_plot(sprat_scaff = sprat_scaff, png_file = paste0("~/Projects/Sprat/doc/assembly_curation/Sprat_", sprat_scaff, "_vs_Ch_.png"), in_data = Sprat_scaff_vs_Ch_chr_df, n_scaffs = 2, ch_sizes = Ch_v2.0.2_sizes)

#Manhattan plots using the spart backbone instead
#Sprat_DeDup_manhattan_GW(sprat_data = sprat_lo_GR, scaff_list = rownames(Sprat_scaff_vs_Ch_chr_tab), daf_col = "fjord_v_brack", png_file = "~/Projects/Sprat/doc/figure_sketches/fjord_v_brack_HiC_DeDup_v2_GW.png", site_filt_vec = sprat_site_filter, backbone = "sprat")
#Sprat_DeDup_manhattan_GW(sprat_data = sprat_lo_GR, scaff_list = rownames(Sprat_scaff_vs_Ch_chr_tab), daf_col = "oce_v_brack", png_file = "~/Projects/Sprat/doc/figure_sketches/oce_v_brack_HiC_DeDup_v2_GW.png", site_filt_vec = sprat_site_filter, backbone = "sprat")
#Sprat_DeDup_manhattan_GW(sprat_data = sprat_lo_GR, scaff_list = rownames(Sprat_scaff_vs_Ch_chr_tab), daf_col = "oce_v_fjord", png_file = "~/Projects/Sprat/doc/figure_sketches/oce_v_fjord_HiC_DeDup_v2_GW.png", site_filt_vec = sprat_site_filter, backbone = "sprat")
Sprat_DeDup_manhattan_GW(sprat_data = sprat_lo_GR, scaff_list = rownames(Sprat_scaff_vs_Ch_chr_tab), daf_col = "fjord_v_brack", png_file = "~/Projects/Sprat/doc/figure_sketches/fjord_v_brack_HiC_DeDup_v2_GW_replot.png", site_filt_vec = sprat_site_filter, backbone = "sprat", plot_w = 10000, plot_h = 2000, th_lwd = 16, point_cex = 3) #Improved resolution & threshhold line
Sprat_DeDup_manhattan_GW(sprat_data = sprat_lo_GR, scaff_list = rownames(Sprat_scaff_vs_Ch_chr_tab), daf_col = "oce_v_brack", png_file = "~/Projects/Sprat/doc/figure_sketches/oce_v_brack_HiC_DeDup_v2_GW_replot.png", site_filt_vec = sprat_site_filter, backbone = "sprat", plot_w = 10000, plot_h = 2000, th_lwd = 16, point_cex = 3) #Improved resolution & threshhold line
Sprat_DeDup_manhattan_GW(sprat_data = sprat_lo_GR, scaff_list = rownames(Sprat_scaff_vs_Ch_chr_tab), daf_col = "oce_v_fjord", png_file = "~/Projects/Sprat/doc/figure_sketches/oce_v_fjord_HiC_DeDup_v2_GW_replot.png", site_filt_vec = sprat_site_filter, backbone = "sprat", plot_w = 10000, plot_h = 2000, th_lwd = 16, point_cex = 3) #Improved resolution & threshhold line



#Inversions on the sprat assembly

#PGA_scaffold_1110__33_contigs__length_20919432: Chr4 inv
#PGA_scaffold_1114__20_contigs__length_23438440:Chr 5 inv 
#PGA_scaffold_13__40_contigs__length_28874545: Chr 7 inv
#PGA_scaffold_173__1_contigs__length_182285, PGA_scaffold_374__1_contigs__length_149044: Chr 13 inv
#PGA_scaffold_40__17_contigs__length_6491151: Chr 18 inv
#PGA_scaffold_1111__78_contigs__length_49772038: Chr 21 inv 

inv_scaffs <- c("PGA_scaffold_1110__33_contigs__length_20919432","PGA_scaffold_1114__20_contigs__length_23438440",
              "PGA_scaffold_13__40_contigs__length_28874545", "PGA_scaffold_173__1_contigs__length_182285", "PGA_scaffold_374__1_contigs__length_149044",
              "PGA_scaffold_40__17_contigs__length_6491151", "PGA_scaffold_1111__78_contigs__length_49772038")
#Sprat_lo_manhattan_v4(sprat_data = sprat_lo_GR, scaff_list = inv_scaffs, daf_col = "oce_v_brack", sel_reg_GR = sal_p20_regions_GR[1], plot_prefix = "~/Projects/Sprat/doc/figure_sketches/oce_v_brack_HiC_DeDup_inv_", site_filt_vec = sprat_site_filter, backbone = "sprat")
Sprat_lo_manhattan_v4(sprat_data = sprat_lo_GR, scaff_list = inv_scaffs, daf_col = "oce_v_brack", sel_reg_GR = sal_p20_regions_GR[1], plot_prefix = "~/Projects/Sprat/doc/figure_sketches/replot_OvB_HiC_DeDup_inv_", site_filt_vec = sprat_site_filter, backbone = "sprat", th_lwd = 8, point_cex = 2.5, point_col = "grey20")



#Examining some in-scaffold cross-mappings
tpm_sp_v_ch <- Sprat_scaff_vs_Ch_chr_df[Sprat_scaff_vs_Ch_chr_df$herring == "chr26",]
table(tpm_sp_v_ch$sprat)
hist(tpm_sp_v_ch$sprat_pos[tpm_sp_v_ch$sprat == "PGA_scaffold_1118__66_contigs__length_56234773"])
sum(tpm_sp_v_ch$sprat_pos[tpm_sp_v_ch$sprat == "PGA_scaffold_1118__66_contigs__length_56234773"] <= 25e6)
#[1] 5
sum(tpm_sp_v_ch$sprat_pos[tpm_sp_v_ch$sprat == "PGA_scaffold_1118__66_contigs__length_56234773"] >= 45e6)
#[1] 149

scaffold_1111_AGP <- as.data.frame(Sprat_m1_AGP_GR[Sprat_m1_AGP_GR$HiC_scaffold == "PGA_scaffold_1111__78_contigs__length_49772038"])
#tail(scaffold_1111_AGP)

tpm_sp_v_ch <- Sprat_scaff_vs_Ch_chr_df[Sprat_scaff_vs_Ch_chr_df$herring == "chr23",]
table(tpm_sp_v_ch$sprat)
hist(tpm_sp_v_ch$sprat_pos[tpm_sp_v_ch$sprat == "PGA_scaffold_1111__78_contigs__length_49772038"])
sum(tpm_sp_v_ch$sprat_pos[tpm_sp_v_ch$sprat == "PGA_scaffold_1111__78_contigs__length_49772038"] >= 11e6 & tpm_sp_v_ch$sprat_pos[tpm_sp_v_ch$sprat == "PGA_scaffold_1111__78_contigs__length_49772038"] <= 15e6)
#[1] 1179
snps_to_plot <- tpm_sp_v_ch$sprat == "PGA_scaffold_1111__78_contigs__length_49772038" & tpm_sp_v_ch$sprat_pos >= 8e6 & tpm_sp_v_ch$sprat_pos <= 15e6
pdf(file = "~/Projects/Sprat/doc/assembly_curation/PGA_scaffold_1111_v_Chr23_contig_breaks.pdf")
plot(x = tpm_sp_v_ch$sprat_pos[snps_to_plot], y = rep(1, sum(snps_to_plot)), xlim = c(8e6, 15e6), col = "darkorchid", pch = 20, ylab = "", xlab = "Sprat position")
segments(x0 = scaffold_1111_AGP$HiC_start, x1 = scaffold_1111_AGP$HiC_end, y0 = 0.8, col = c("olivedrab", "black"), lwd = c(1,6))
abline(v =scaffold_1111_AGP$HiC_start[scaffold_1111_AGP$seqnames == "HiC_gap"], lwd = 1, lty = "dashed", col = "grey40")
dev.off()


tpm_sp_v_ch <- Sprat_scaff_vs_Ch_chr_df[Sprat_scaff_vs_Ch_chr_df$herring == "chr21",]
table(tpm_sp_v_ch$sprat)
hist(tpm_sp_v_ch$sprat_pos[tpm_sp_v_ch$sprat == "PGA_scaffold_1111__78_contigs__length_49772038"])
sum(tpm_sp_v_ch$sprat_pos[tpm_sp_v_ch$sprat == "PGA_scaffold_1111__78_contigs__length_49772038"] >= 45e6)
#[1] 4028

pdf(file = "~/Projects/Sprat/doc/assembly_curation/PGA_scaffold_1111_v_Chr21_contig_breaks.pdf")
snps_to_plot <- tpm_sp_v_ch$sprat == "PGA_scaffold_1111__78_contigs__length_49772038" & tpm_sp_v_ch$sprat_pos >= 45e6
plot(x = tpm_sp_v_ch$sprat_pos[snps_to_plot], y = rep(1, sum(snps_to_plot)), xlim = c(46e6, 50e6), col = "darkorchid", pch = 20, ylab = "", xlab = "Sprat position")
segments(x0 = scaffold_1111_AGP$HiC_start, x1 = scaffold_1111_AGP$HiC_end, y0 = 0.8, col = c("olivedrab", "black"), lwd = c(1,6))
abline(v =scaffold_1111_AGP$HiC_start[scaffold_1111_AGP$seqnames == "HiC_gap"], lwd = 1, lty = "dashed", col = "grey40")
dev.off()


#Support functions
Sp_scaff_v_Ch_chr_plot <- function(sprat_scaff = NULL, ch_chr = NULL, ch_sizes = NULL, png_file, in_data, n_scaffs = 1){
  if(is.null(sprat_scaff) & is.null(ch_chr)){
    print("At least one target needed!")
    return(0)
  } else if(is.null(sprat_scaff)){
    tmp_sprat_tab <- table(in_data$sprat[in_data$herring == ch_chr])
    sprat_scaff <-  names(tmp_sprat_tab)[order(tmp_sprat_tab, decreasing = T)][1:n_scaffs]
  } else if(is.null(ch_chr)){
    tmp_herring_tab <- table(in_data$herring[in_data$sprat == sprat_scaff])
    ch_chr <-  names(tmp_herring_tab)[order(tmp_herring_tab, decreasing = T)][1:n_scaffs]
  } 
  print(ch_chr)
  
  Sp_scaff_v_Ch_chr_df <- in_data[which(in_data$sprat %in% sprat_scaff),]
  Sp_scaff_v_Ch_chr_df <- Sp_scaff_v_Ch_chr_df[which(Sp_scaff_v_Ch_chr_df$herring %in% ch_chr),]
  png(file = png_file, height = 600, width = 1000)
  
  if(n_scaffs == 1){
    print("Type 1")
    flush.console()
    
    plot(x = in_data$sprat_pos[1], y = 1, pch = 20, xlab = "Position", main = "", ylab = "", ylim = c(0.4, 1.6), xlim = c(0, max(c(Sp_scaff_v_Ch_chr_df$sprat_pos, Sp_scaff_v_Ch_chr_df$herring_pos))), type = "n")
    segments(x0 = 0, x1 = c(max(Sp_scaff_v_Ch_chr_df$herring_pos), max(Sp_scaff_v_Ch_chr_df$sprat_pos)), y0 = c(1.5, 0.5), col = c("firebrick", "darkorchid"), lwd = 2)
    text(labels = c(sprat_scaff, ch_chr), x = rep(mean(c(par("usr")[1], par("usr")[2])), 2), y = c(0.4, 1.6))
    segments(x0 = Sp_scaff_v_Ch_chr_df$herring_pos, x1 = Sp_scaff_v_Ch_chr_df$sprat_pos, y0 = 1.5, y1 =  0.5, col = "slateblue", lwd = 0.3)
  }
  
  if(n_scaffs == 2 & length(ch_chr) == 1){
    #One herring chr vs 2 sprat scaffolds
    print("Type 2")
    flush.console()
    
    sprat_scaff_1 <-  sprat_scaff[1]
    sprat_max_1 <- as.numeric(sub(".+_length_([0-9]+)", "\\1", sprat_scaff_1))
   
    sprat_scaff_2 <-  sprat_scaff[2]
    sprat_max_2 <- as.numeric(sub(".+_length_([0-9]+)", "\\1", sprat_scaff_2))
    
    plot(x = in_data$sprat_pos[1], y = 1, pch = 20, xlab = "Position", main = "", ylab = "", ylim = c(0.4, 1.6), xlim = c(0, max(c(sprat_max_1, sprat_max_2, in_data$herring_pos*1.1))), type = "n")
    
    segments(x0 = Sp_scaff_v_Ch_chr_df$herring_pos[Sp_scaff_v_Ch_chr_df$sprat == sprat_scaff_1], x1 = Sp_scaff_v_Ch_chr_df$sprat_pos[Sp_scaff_v_Ch_chr_df$sprat == sprat_scaff_1], y0 = 1, y1 =  0.5, col = "slateblue", lwd = 0.3)
    segments(x0 = 0, x1 = sprat_max_1, y0 = 0.5, col = "darkorchid", lwd = 5)
    
    
    segments(x0 = Sp_scaff_v_Ch_chr_df$herring_pos[Sp_scaff_v_Ch_chr_df$sprat == sprat_scaff_2], x1 = Sp_scaff_v_Ch_chr_df$sprat_pos[Sp_scaff_v_Ch_chr_df$sprat == sprat_scaff_2], y0 = 1, y1 =  1.5, col = "steelblue", lwd = 0.3)
    segments(x0 = 0, x1 = sprat_max_2, y0 = 1.5, col = "salmon", lwd = 5)
    
    segments(x0 = 0, x1 = c(max(Sp_scaff_v_Ch_chr_df$herring_pos)), y0 = 1, col = c("firebrick"), lwd = 8)
    
    text(labels = c(sprat_scaff_1,sprat_scaff_2, ch_chr), x = c(rep(mean(c(par("usr")[1], par("usr")[2])), 2), max(Sp_scaff_v_Ch_chr_df$herring_pos*1.05)), y = c(0.4, 1.6, 1), cex = 1.5)
  }
  
  if(n_scaffs == 2 & length(sprat_scaff) == 1){
    #One sprat scaffold vs 2 herring Chrs
    print("Type 3")
    flush.console()
    
    ch_chr_1 <-  ch_chr[1]
    ch_max_1 <- ch_sizes$size[ch_sizes$name == ch_chr_1]
    ch_max_2 <- 0 #in case there is no second Chr
    
    ch_chr_2 <-  ch_chr[2]
    if(!is.na(ch_chr_2)) ch_max_2 <- ch_sizes$size[ch_sizes$name == ch_chr_2]
    
    
    #print(paste0("first max: ", ch_max_1, "; second max: ", ch_max_2))
    #flush.console()
    
    plot(x = Sp_scaff_v_Ch_chr_df$herring_pos[1], y = 1, pch = 20, xlab = "Position", main = "", ylab = "", ylim = c(0.4, 1.6), xlim = c(0, max(c(ch_max_1, ch_max_2, Sp_scaff_v_Ch_chr_df$sprat_pos*1.25))), type = "n")
    segments(x0 = Sp_scaff_v_Ch_chr_df$sprat_pos[Sp_scaff_v_Ch_chr_df$herring == ch_chr_1], x1 = Sp_scaff_v_Ch_chr_df$herring_pos[Sp_scaff_v_Ch_chr_df$herring == ch_chr_1], y0 = 1, y1 =  0.5, col = "slateblue", lwd = 0.3)
    segments(x0 = 0, x1 = ch_max_1, y0 = 0.5, col = "darkorchid", lwd = 5)
    
    segments(x0 = Sp_scaff_v_Ch_chr_df$sprat_pos[Sp_scaff_v_Ch_chr_df$herring == ch_chr_2], x1 = Sp_scaff_v_Ch_chr_df$herring_pos[Sp_scaff_v_Ch_chr_df$herring == ch_chr_2], y0 = 1, y1 =  1.5, col = "steelblue", lwd = 0.3)
    segments(x0 = 0, x1 = ch_max_2, y0 = 1.5, col = "salmon", lwd = 5)
 
    segments(x0 = 0, x1 = c(max(Sp_scaff_v_Ch_chr_df$sprat_pos)), y0 = 1, col = c("firebrick"), lwd = 8)
    
    text(labels = c(ch_chr_1,ch_chr_2, sub("__.+", "", sprat_scaff)), x = c(rep(mean(c(par("usr")[1], par("usr")[2])), 2), max(Sp_scaff_v_Ch_chr_df$sprat_pos*1.05)),pos =c(NULL, NULL, 4), y = c(0.4, 1.6, 1), cex = 1.5)
  }
  dev.off()
}

Sprat_lo_manhattan_v4 <- function(sprat_data, scaff_list, daf_col, site_filt_vec,  sel_reg_GR, plot_prefix, backbone = "sprat", th_lwd = 1.5, point_cex = 0.9, point_col = "gray80"){
  daf_thresh_val <- qnorm(p = 1-(0.05/sum(site_filt_vec)), mean = mean(sprat_data@elementMetadata[site_filt_vec, daf_col]), sd  = sd(sprat_data@elementMetadata[site_filt_vec, daf_col]))
  for(chr in scaff_list){
    tmp_png_file <- paste0(plot_prefix,chr, ".png")
    tmp_chr_filter <- NULL
    png(filename = tmp_png_file, height = 1000, width = 2000)
    
    if(backbone == "sprat"){
      tmp_chr_filter <- sprat_data$CHROM == chr
      plot(x = 1, y = 0, col = "black", xlim = c(0, max(sprat_data$POS[site_filt_vec & tmp_chr_filter])), ylim = c(0,1), xlab  = "Position", ylab = "DAF", main = chr, type = "n")
      points(x = sprat_data$POS[site_filt_vec & tmp_chr_filter], y = sprat_data@elementMetadata[site_filt_vec & tmp_chr_filter, daf_col], col = point_col, pch = 16, cex = point_cex)
    }
    
    
    if(backbone == "herring"){
      tmp_chr_filter <- 1:length(sprat_data$Ch_v2.0.2_chr) %in% which(sprat_data$Ch_v2.0.2_chr == chr)
      tmp_reg_filter <- as.character(sel_reg_GR@seqnames) == chr
      plot(x = 1, y = 0, col = "black", xlim = c(0, max(sprat_data$Ch_v2.0.2_pos[site_filt_vec & tmp_chr_filter])), ylim = c(0,1), xlab  = "Position", ylab = "DAF", main = chr, type = "n")
      if(sum(tmp_reg_filter) > 0) rect(xleft = sel_reg_GR@ranges@start[tmp_reg_filter], xright = sel_reg_GR@ranges@start[tmp_reg_filter] + sel_reg_GR@ranges@width[tmp_reg_filter], ybottom = -2, ytop = 2, col = "grey90", border = NA)
      points(x = sprat_data$Ch_v2.0.2_pos[site_filt_vec & tmp_chr_filter], y = sprat_data@elementMetadata[site_filt_vec & tmp_chr_filter, daf_col], col = point_col, pch = 16, cex = point_cex)
    }
    
    abline(h = daf_thresh_val, col = "red", lwd = th_lwd)
    
    dev.off()
  }
}


simple_sprat_pool_hm_Ch <- function(snp_data, reg_filter, sample_info_df, pdf_file, freq_ref_pools, ordered_plot_pools = NULL){
  tmp_hm_freq <- as.matrix(snp_data[reg_filter, grep("_Afreq", names(snp_data))])
  rownames(tmp_hm_freq) <- paste(snp_data$herring_v2.0.2_seqnames[reg_filter], snp_data$SNP_HiC_pos[reg_filter], sep = "_")
  
  tmp_sprat_IDs <- sub("_Afreq", "", grep("_Afreq", colnames(tmp_hm_freq), value = T))
  tmp_sprat_pool_names <- sample_info_df$clean_name[match(tmp_sprat_IDs, sample_info_df$POOL)]
  colnames(tmp_hm_freq) <- tmp_sprat_pool_names
  
  freq_flip_vec <- rowMeans(tmp_hm_freq[,freq_ref_pools]) > 0.5
  tmp_hm_freq[freq_flip_vec,] <- (1 - tmp_hm_freq[freq_flip_vec,])
  if(is.null(ordered_plot_pools)){
    tmp_hm_freq <- tmp_hm_freq[order(snp_data$SNP_HiC_pos[reg_filter], decreasing = F),]
  } else{
    pool_order <- match(ordered_plot_pools, tmp_sprat_pool_names)
    tmp_hm_freq <- tmp_hm_freq[order(snp_data$SNP_HiC_pos[reg_filter], decreasing = F),pool_order]
  }
  
  pdf(file = pdf_file)
  if(is.null(ordered_plot_pools)){
    heatmap(t(tmp_hm_freq), scale = "none", Colv = NA, margins = c(5,10))
  } else{
    heatmap(t(tmp_hm_freq), scale = "none", Colv = NA, Rowv = NA, margins = c(5,10))
  }
  dev.off()
}

simple_sprat_pool_hm_v4 <- function(snp_data, reg_filter, sample_info_df, pdf_file, freq_ref_pools, ordered_plot_pools = NULL){
  tmp_hm_freq <- as.matrix(snp_data@elementMetadata[reg_filter, grep("_Afreq", names(snp_data@elementMetadata))])
  #rownames(tmp_hm_freq) <- paste(snp_data$HiC_scaffold[reg_filter], snp_data$HiC_pos[reg_filter], sep = "_")
  rownames(tmp_hm_freq) <- paste(snp_data$CHROM[reg_filter], snp_data$POS[reg_filter], sep = "_")
  tmp_sprat_IDs <- sub("_Afreq", "", grep("_Afreq", colnames(tmp_hm_freq), value = T))
  tmp_sprat_pool_names <- sample_info_df$clean_name[match(tmp_sprat_IDs, sample_info_df$POOL)]
  colnames(tmp_hm_freq) <- tmp_sprat_pool_names
  freq_flip_vec <- rowMeans(tmp_hm_freq[,freq_ref_pools]) > 0.5
  tmp_hm_freq[freq_flip_vec,] <- (1 - tmp_hm_freq[freq_flip_vec,])
  if(is.null(ordered_plot_pools)){
    #tmp_hm_freq <- tmp_hm_freq[order(snp_data$HiC_pos[reg_filter], decreasing = F),]
    tmp_hm_freq <- tmp_hm_freq[order(snp_data$POS[reg_filter], decreasing = F),]
  } else{
    pool_order <- match(ordered_plot_pools, tmp_sprat_pool_names)
    #tmp_hm_freq <- tmp_hm_freq[order(snp_data$HiC_pos[reg_filter], decreasing = F),pool_order]
    tmp_hm_freq <- tmp_hm_freq[order(snp_data$POS[reg_filter], decreasing = F),pool_order]
  }
  
  pdf(file = pdf_file)
  if(is.null(ordered_plot_pools)){
    heatmap(t(tmp_hm_freq), scale = "none", Colv = NA, margins = c(5,10))
  } else{
    heatmap(t(tmp_hm_freq), scale = "none", Colv = NA, Rowv = NA, margins = c(5,10))
  }
  dev.off()
}

Sprat_lo_manhattan_v3 <- function(sprat_data, scaff_list, daf_col, site_filt_vec,  sel_reg_GR, plot_prefix, backbone = "sprat", th_lwd = 1.5, point_cex = 0.9){
  daf_thresh_val <- qnorm(p = 1-(0.05/sum(site_filt_vec)), mean = mean(sprat_data@elementMetadata[site_filt_vec, daf_col]), sd  = sd(sprat_data@elementMetadata[site_filt_vec, daf_col]))
  for(chr in scaff_list){
    tmp_png_file <- paste0(plot_prefix,chr, ".png")
    tmp_chr_filter <- NULL
    png(filename = tmp_png_file, height = 1000, width = 2000)
    
    if(backbone == "sprat"){
      tmp_chr_filter <- sprat_data$HiC_scaffold == chr
      plot(x = 1, y = 0, col = "black", xlim = c(0, max(sprat_data$HiC_pos[site_filt_vec & tmp_chr_filter])), ylim = c(0,1), xlab  = "Position", ylab = "DAF", main = chr, type = "n")
      points(x = sprat_data$HiC_pos[site_filt_vec & tmp_chr_filter], y = sprat_data@elementMetadata[site_filt_vec & tmp_chr_filter, daf_col], col = "gray80", pch = 16, cex = point_cex)
    }
    
  
    if(backbone == "herring"){
      tmp_chr_filter <- 1:length(sprat_data$Ch_v2.0.2_chr) %in% which(sprat_data$Ch_v2.0.2_chr == chr)
      tmp_reg_filter <- as.character(sel_reg_GR@seqnames) == chr
      plot(x = 1, y = 0, col = "black", xlim = c(0, max(sprat_data$Ch_v2.0.2_pos[site_filt_vec & tmp_chr_filter])), ylim = c(0,1), xlab  = "Position", ylab = "DAF", main = chr, type = "n")
      if(sum(tmp_reg_filter) > 0) rect(xleft = sel_reg_GR@ranges@start[tmp_reg_filter], xright = sel_reg_GR@ranges@start[tmp_reg_filter] + sel_reg_GR@ranges@width[tmp_reg_filter], ybottom = -2, ytop = 2, col = "grey90", border = NA)
      points(x = sprat_data$Ch_v2.0.2_pos[site_filt_vec & tmp_chr_filter], y = sprat_data@elementMetadata[site_filt_vec & tmp_chr_filter, daf_col], col = "gray30", pch = 16, cex = point_cex)
    }
    
    abline(h = daf_thresh_val, col = "red", lwd = th_lwd)

    dev.off()
  }
}

Sprat_lo_manhattan_GW <- function(sprat_data, scaff_list, daf_col, site_filt_vec, png_file, backbone = "sprat", y_lim =  c(0,1), sel_reg = NULL){
  daf_thresh_val <- qnorm(p = 1-(0.05/sum(site_filt_vec)), mean = mean(sprat_data@elementMetadata[site_filt_vec, daf_col]), sd  = sd(sprat_data@elementMetadata[site_filt_vec, daf_col]))
  #for(chr in scaff_list){
  #tmp_png_file <- paste0(plot_prefix,chr, ".png")
  #tmp_chr_filter <- NULL
  png(filename = png_file, height = 1000, width = 10000)
  
  if(backbone == "sprat"){
    tmp_chr_filter <- (1:length(sprat_data$HiC_scaffold) %in% which(sprat_data$HiC_scaffold %in% scaff_list)) & site_filt_vec
    #TO ADD: cumulative position of included scaffolds.
    scaff_maxes <- aggregate(sprat_data$HiC_pos[tmp_chr_filter]~sprat_data$HiC_scaffold[tmp_chr_filter], FUN = "max")
    scaff_maxes$offset <- c(0,cumsum(scaff_maxes[,2])[-dim(scaff_maxes)[1]])
    scaff_maxes$col <- rep(c("grey30", "grey70"), times = ceiling(length(scaff_maxes$offset)/2))[1:length(scaff_maxes$offset)]
    cumulative_pos_vec <- sprat_data$HiC_pos[tmp_chr_filter] + scaff_maxes$offset[match(sprat_data$CHROM[tmp_chr_filter], scaff_maxes[,1])]
    cumulative_col_vec <- scaff_maxes$col[match(sprat_data$CHROM[tmp_chr_filter], scaff_maxes[,1])]
    
    plot(x = 1, y = 0, col = "black", xlim = c(0, max(cumulative_pos_vec)), ylim = y_lim, xlab  = "Position", ylab = "DAF", main = "", type = "n")
    points(x = cumulative_pos_vec, y = sprat_data@elementMetadata[tmp_chr_filter, daf_col], col = cumulative_col_vec, pch = 16, cex = 0.9)
  }
  
  
  if(backbone == "herring"){
    tmp_chr_filter <- (1:length(sprat_data$Ch_v2.0.2_chr) %in% which(sprat_data$Ch_v2.0.2_chr %in% scaff_list)) & site_filt_vec
    scaff_maxes <- aggregate(sprat_data$Ch_v2.0.2_pos[tmp_chr_filter]~sprat_data$Ch_v2.0.2_chr[tmp_chr_filter], FUN = "max")
    scaff_maxes <- scaff_maxes[match(scaff_list, scaff_maxes[,1]),]
    
    scaff_maxes$offset <- c(0,cumsum(scaff_maxes[,2])[-dim(scaff_maxes)[1]])
    scaff_maxes$col <- rep(c("grey30", "grey70"), times = ceiling(length(scaff_maxes$offset)/2))[1:length(scaff_maxes$offset)]
    cumulative_pos_vec <- sprat_data$Ch_v2.0.2_pos[tmp_chr_filter] + scaff_maxes$offset[match(sprat_data$Ch_v2.0.2_chr[tmp_chr_filter], scaff_maxes[,1])]
    cumulative_col_vec <- scaff_maxes$col[match(sprat_data$Ch_v2.0.2_chr[tmp_chr_filter], scaff_maxes[,1])]
    #tmp_reg_filter <- as.character(sel_reg_GR@seqnames) == chr
    plot(x = 1, y = 0, col = "black", xlim = c(0, max(cumulative_pos_vec)), ylim = y_lim, xlab  = "Position", ylab = "DAF", main = "", type = "n")
    #if(sum(tmp_reg_filter) > 0) rect(xleft = sel_reg_GR@ranges@start[tmp_reg_filter], xright = sel_reg_GR@ranges@start[tmp_reg_filter] + sel_reg_GR@ranges@width[tmp_reg_filter], ybottom = -2, ytop = 2, col = "grey90", border = NA)
    points(x = cumulative_pos_vec, y = sprat_data@elementMetadata[tmp_chr_filter, daf_col], col = cumulative_col_vec, pch = 16, cex = 0.9)
  }
  
  abline(h = daf_thresh_val, col = "red", lwd = 1.5)
  if(!is.null(sel_reg)){
    segments(x0 = sel_reg$cumulative_start, x1 = sel_reg$cumulative_start + as.numeric(sel_reg@ranges@width), y0 = 1.0, y1 = 1.0, col = "darkorchid", lwd = 8)
  }
  
  dev.off()
  #}
  return(scaff_maxes)
}

Sprat_DeDup_manhattan_GW <- function(sprat_data, scaff_list, daf_col, site_filt_vec, png_file, backbone = "sprat", y_lim =  c(0,1), sel_reg = NULL, plot_w = 10000, plot_h = 2000, th_lwd = 2.5, point_cex = 0.9){
  daf_thresh_val <- qnorm(p = 1-(0.05/sum(site_filt_vec)), mean = mean(sprat_data@elementMetadata[site_filt_vec, daf_col]), sd  = sd(sprat_data@elementMetadata[site_filt_vec, daf_col]))
  #for(chr in scaff_list){
  #tmp_png_file <- paste0(plot_prefix,chr, ".png")
  #tmp_chr_filter <- NULL
  png(filename = png_file, height = plot_h, width = plot_w)
  
  if(backbone == "sprat"){
    tmp_chr_filter <- (1:length(sprat_data$CHROM) %in% which(sprat_data$CHROM %in% scaff_list)) & site_filt_vec
    #TO ADD: cumulative position of included scaffolds.
    scaff_maxes <- aggregate(sprat_data$POS[tmp_chr_filter]~sprat_data$CHROM[tmp_chr_filter], FUN = "max")
    scaff_maxes <- scaff_maxes[match(scaff_list, scaff_maxes[,1]),]
    scaff_maxes$offset <- c(0,cumsum(scaff_maxes[,2])[-dim(scaff_maxes)[1]])
    scaff_maxes$col <- rep(c("grey30", "grey70"), times = ceiling(length(scaff_maxes$offset)/2))[1:length(scaff_maxes$offset)]
    cumulative_pos_vec <- sprat_data$POS[tmp_chr_filter] + scaff_maxes$offset[match(sprat_data$CHROM[tmp_chr_filter], scaff_maxes[,1])]
    cumulative_col_vec <- scaff_maxes$col[match(sprat_data$CHROM[tmp_chr_filter], scaff_maxes[,1])]
    
    plot(x = 1, y = 0, col = "black", xlim = c(0, max(cumulative_pos_vec)), ylim = y_lim, xlab  = "Position", ylab = "DAF", main = "", type = "n")
    points(x = cumulative_pos_vec, y = sprat_data@elementMetadata[tmp_chr_filter, daf_col], col = cumulative_col_vec, pch = 16, cex = point_cex)
  }
  
  
  if(backbone == "herring"){
    tmp_chr_filter <- (1:length(sprat_data$Ch_v2.0.2_chr) %in% which(sprat_data$Ch_v2.0.2_chr %in% scaff_list)) & site_filt_vec
    scaff_maxes <- aggregate(sprat_data$Ch_v2.0.2_pos[tmp_chr_filter]~sprat_data$Ch_v2.0.2_chr[tmp_chr_filter], FUN = "max")
    scaff_maxes <- scaff_maxes[match(scaff_list, scaff_maxes[,1]),]
    
    scaff_maxes$offset <- c(0,cumsum(scaff_maxes[,2])[-dim(scaff_maxes)[1]])
    scaff_maxes$col <- rep(c("grey30", "grey70"), times = ceiling(length(scaff_maxes$offset)/2))[1:length(scaff_maxes$offset)]
    cumulative_pos_vec <- sprat_data$Ch_v2.0.2_pos[tmp_chr_filter] + scaff_maxes$offset[match(sprat_data$Ch_v2.0.2_chr[tmp_chr_filter], scaff_maxes[,1])]
    cumulative_col_vec <- scaff_maxes$col[match(sprat_data$Ch_v2.0.2_chr[tmp_chr_filter], scaff_maxes[,1])]
    #tmp_reg_filter <- as.character(sel_reg_GR@seqnames) == chr
    plot(x = 1, y = 0, col = "black", xlim = c(0, max(cumulative_pos_vec)), ylim = y_lim, xlab  = "Position", ylab = "DAF", main = "", type = "n")
    #if(sum(tmp_reg_filter) > 0) rect(xleft = sel_reg_GR@ranges@start[tmp_reg_filter], xright = sel_reg_GR@ranges@start[tmp_reg_filter] + sel_reg_GR@ranges@width[tmp_reg_filter], ybottom = -2, ytop = 2, col = "grey90", border = NA)
    points(x = cumulative_pos_vec, y = sprat_data@elementMetadata[tmp_chr_filter, daf_col], col = cumulative_col_vec, pch = 16, cex = 0.9)
  }
  
  abline(h = daf_thresh_val, col = "firebrick", lwd = th_lwd)
  if(!is.null(sel_reg)){
    segments(x0 = sel_reg$cumulative_start, x1 = sel_reg$cumulative_start + as.numeric(sel_reg@ranges@width), y0 = 1.0, y1 = 1.0, col = "darkorchid", lwd = 8)
  }
  
  dev.off()
  #}
  return(scaff_maxes)
}

Sprat_lo_manhattan_Ch <- function(sprat_data, scaff_list, daf_col, site_filt_vec,  sel_reg_GR, plot_prefix){
  daf_thresh_val <- qnorm(p = 1-(1/sum(site_filt_vec)), mean = mean(sprat_data[site_filt_vec, daf_col]), sd  = sd(sprat_data[site_filt_vec, daf_col]))
  for(chr in scaff_list){
    tmp_png_file <- paste0(plot_prefix,chr, ".png")
    tmp_chr_filter <- sprat_data$herring_v2.0.2_seqnames == chr
    png(filename = tmp_png_file, height = 1000, width = 2000)
    plot(x = sprat_data$SNP_HiC_pos[site_filt_vec & tmp_chr_filter], y = sprat_data[site_filt_vec & tmp_chr_filter, daf_col], col = "black", ylim = c(0,1), pch = 16, cex = 0.8, xlab  = "Position", ylab = "DAF", main = chr)
    tmp_reg_filter <- as.character(sel_reg_GR@seqnames) == chr
    if(sum(tmp_reg_filter) > 0) rect(xleft = sel_reg_GR@ranges@start[tmp_reg_filter], xright = sel_reg_GR@ranges@start[tmp_reg_filter] + sel_reg_GR@ranges@width[tmp_reg_filter], ybottom = 0, ytop = 1, col = rgb(1,0,0,0.2), border = NULL)
    abline(h = daf_thresh_val, col = "red")
    
    dev.off()
  }
}


