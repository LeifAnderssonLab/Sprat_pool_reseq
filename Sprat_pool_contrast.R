#Author: Mats Pettersson
#mats.pettersson@imbim.uu.se

oceanic_freq <- names(Sprat_pool_freq)[grepl("_Afreq", names(Sprat_pool_freq)) & grepl("NS|BoB|SK|CEL", names(Sprat_pool_freq))]
brackish_freq <- names(Sprat_pool_freq)[grepl("_Afreq", names(Sprat_pool_freq)) & grepl("AB|BS|GOTB|BBS|LAND|GD", names(Sprat_pool_freq))]
fjord_freq <- names(Sprat_pool_freq)[grepl("_Afreq", names(Sprat_pool_freq)) & grepl("NOR|HAR|SOG", names(Sprat_pool_freq))]

require(matrixStats)
Sprat_pool_freq$mean_oceanic <- rowMeans(Sprat_pool_freq[,oceanic_freq], na.rm = T)
Sprat_pool_freq$mean_brackish <- rowMeans(Sprat_pool_freq[,brackish_freq], na.rm = T)
Sprat_pool_freq$mean_fjord <- rowMeans(Sprat_pool_freq[,fjord_freq], na.rm = T)

Sprat_pool_freq$oce_v_brack <- abs(Sprat_pool_freq$mean_brackish - Sprat_pool_freq$mean_oceanic)
Sprat_pool_freq$fjord_v_brack <- abs(Sprat_pool_freq$mean_brackish - Sprat_pool_freq$mean_fjord)
Sprat_pool_freq$oce_v_fjord <- abs(Sprat_pool_freq$mean_fjord - Sprat_pool_freq$mean_oceanic)

#save(Sprat_pool_freq, sprat_site_filter, file = "~/Projects/Sprat/data/genotypes/IPA_pri_genos/Sprat_pool_freq_ext.RData") ##Original, IPA-pri version
save(Sprat_pool_freq, file = "~/Projects/Sprat/data/genotypes/Sprat_pool_freq_DeDup_v2_ext.RData") ##Current, HiC_DeDup_v2 version



#Adding some site filtering a la Multispecies SNP chip
sprat_opp_fix_filter <- (rowSums(Sprat_pool_freq[,grep("_Afreq", names(Sprat_pool_freq))] == 0, na.rm = T) > 0) & (rowSums(Sprat_pool_freq[,grep("_Afreq", names(Sprat_pool_freq))] == 1, na.rm = T) > 0)
#Attmepting to eliminate SNPs where the large majority of counts are exactly 0 or 1
head(Sprat_pool_freq[sprat_opp_fix_filter,c(1,2,grep("_Afreq", names(Sprat_pool_freq)))])
sprat_pool_dp_filter <- !(rowSums(Sprat_pool_freq[,grep("_DP", names(Sprat_pool_freq))] <= 10, na.rm = T) > 0)
head(Sprat_pool_freq[sprat_pool_dp_filter,])
sprat_site_filter <- !sprat_opp_fix_filter & sprat_pool_dp_filter

sprat_lo_GR <- GRanges(seqnames = Sprat_pool_freq$CHROM, ranges = IRanges(start = Sprat_pool_freq$POS, end = Sprat_pool_freq$POS))
for(i in 1:dim(Sprat_pool_freq)[2]){
  sprat_lo_GR@elementMetadata[,i] <- Sprat_pool_freq[,i]
  names(sprat_lo_GR@elementMetadata)[i] <- names(Sprat_pool_freq)[i]
}

save(sprat_lo_GR, sprat_site_filter, file = "~/Projects/Sprat/data/genotypes/Sprat_pool_DeDup_v2_GR.RData") ##Current, HiC_DeDup_v2 version


##Checking an expected inversion
##Target scaffold: PGA_scaffold_40__17_contigs__length_6491151
png(filename = "~/Projects/Sprat/doc/oce_v_brack_DeDup_v2_scaff40.png", width = 1000)
plot_filt <- sprat_lo_GR$CHROM == "PGA_scaffold_40__17_contigs__length_6491151" & sprat_site_filter
plot(x = sprat_lo_GR$POS[plot_filt], y = sprat_lo_GR$oce_v_brack[plot_filt], pch = 16, cex = 0.3, ylab = "DAF (Oce v Brackish)", main = "PGA_scaffold_40__17_contigs", xlab = "Position")
dev.off()

png(filename = "~/Projects/Sprat/doc/oce_v_brack.png", width = 1000)
plot(x = 1:sum(sprat_site_filter), y = Sprat_pool_freq$oce_v_brack[sprat_site_filter], pch = 16, cex = 0.3, ylab = "DAF")
lines(x = 1:sum(sprat_site_filter), y = stats::filter(Sprat_pool_freq$oce_v_brack[sprat_site_filter], filter = rep(1/100, 100)), col = "red")
dev.off()

hist(stats::filter(Sprat_pool_freq$oce_v_brack[sprat_site_filter], filter = rep(1/100, 100)))

png(filename = "~/Projects/Sprat/doc/oce_v_fjord.png", width = 1000)
plot(x = 1:sum(sprat_site_filter), y = Sprat_pool_freq$oce_v_fjord[sprat_site_filter], pch = 16, cex = 0.3, ylab = "DAF")
lines(x = 1:sum(sprat_site_filter), y = stats::filter(Sprat_pool_freq$oce_v_fjord[sprat_site_filter], filter = rep(1/100, 100)), col = "red")
dev.off()

png(filename = "~/Projects/Sprat/doc/brack_v_fjord.png", width = 1000)
plot(x = 1:sum(sprat_site_filter), y = Sprat_pool_freq$fjord_v_brack[sprat_site_filter], pch = 16, cex = 0.3, ylab = "DAF")
lines(x = 1:sum(sprat_site_filter), y = stats::filter(Sprat_pool_freq$fjord_v_brack[sprat_site_filter], filter = rep(1/100, 100)), col = "red")
dev.off()

png(filename = "~/Projects/Sprat/doc/oce_contrast_cor.png", width = 1000)
plot(x =Sprat_pool_freq$oce_v_brack[sprat_site_filter], y = Sprat_pool_freq$oce_v_fjord[sprat_site_filter], pch = 16, cex = 0.3, ylab = "Oce_v_Fjord", xlab = "Oce_v_Brack")
dev.off()

oce_v_brack_high_DAF <- Sprat_pool_freq$oce_v_brack > 0.7 & sprat_site_filter 

tmp_filter <- Sprat_pool_freq$CHROM == "ctg.001058F" & sprat_site_filter  
png(filename = "~/Projects/Sprat/doc/oce_v_brack_001058F.png", width = 1000)
plot(x = Sprat_pool_freq$POS[tmp_filter], y = Sprat_pool_freq$oce_v_brack[tmp_filter], pch = 16, cex = 0.3, ylab = "DAF")
lines(x = Sprat_pool_freq$POS[tmp_filter], y = stats::filter(Sprat_pool_freq$oce_v_brack[tmp_filter], filter = rep(1/100, 100)), col = "red")
dev.off()

tmp_filter <- Sprat_pool_freq$CHROM == "ctg.000449F" & sprat_site_filter  
png(filename = "~/Projects/Sprat/doc/oce_v_brack_000449F.png", width = 1000)
plot(x = Sprat_pool_freq$POS[tmp_filter], y = Sprat_pool_freq$oce_v_brack[tmp_filter], pch = 16, cex = 0.3, ylab = "DAF")
lines(x = Sprat_pool_freq$POS[tmp_filter], y = stats::filter(Sprat_pool_freq$oce_v_brack[tmp_filter], filter = rep(1/100, 100)), col = "red")
dev.off()

tmp_hm_freq <- as.matrix(test_sprat_SNPs_v2$SNP_set$picked_SNPs[test_sprat_SNPs_v2$SNP_set$picked_SNPs$CHROM == "ctg.000449F", -c(1,2)])
rownames(tmp_hm_freq) <- paste(test_sprat_SNPs_v2$SNP_set$picked_SNPs[test_sprat_SNPs_v2$SNP_set$picked_SNPs$CHROM == "ctg.000449F", 1], test_sprat_SNPs_v2$SNP_set$picked_SNPs[test_sprat_SNPs_v2$SNP_set$picked_SNPs$CHROM == "ctg.000449F", 2], sep = "_")
pdf(file = "~/Projects/Sprat/doc/Ctg_000449F_hm.pdf")
  heatmap(t(tmp_hm_freq), scale = "none", Colv = NA, margins = c(5,10))
dev.off()



#Fjord-unique
table(Sprat_pool_freq$CHROM[Sprat_pool_freq$oce_v_fjord > 0.4 & Sprat_pool_freq$oce_v_brack < 0.1 & sprat_site_filter])
#ctg.000013F: 1153

tmp_filter <- Sprat_pool_freq$CHROM == "ctg.000013F" & sprat_site_filter  
png(filename = "~/Projects/Sprat/doc/ctg.000013F_oce_v_brack.png", width = 1000)
plot(x = Sprat_pool_freq$POS[tmp_filter], y = Sprat_pool_freq$oce_v_brack[tmp_filter], pch = 16, cex = 0.3, ylab = "DAF")
lines(x = Sprat_pool_freq$POS[tmp_filter], y = stats::filter(Sprat_pool_freq$oce_v_brack[tmp_filter], filter = rep(1/100, 100)), col = "red")
dev.off()

png(filename = "~/Projects/Sprat/doc/ctg.000013F_oce_v_fjord.png", width = 1000)
plot(x = Sprat_pool_freq$POS[tmp_filter], y = Sprat_pool_freq$oce_v_fjord[tmp_filter], pch = 16, cex = 0.3, ylab = "DAF")
lines(x = Sprat_pool_freq$POS[tmp_filter], y = stats::filter(Sprat_pool_freq$oce_v_fjord[tmp_filter], filter = rep(1/100, 100)), col = "red")
dev.off()

tmp_hm_freq <- as.matrix(test_sprat_SNPs_v2$SNP_set$picked_SNPs[test_sprat_SNPs_v2$SNP_set$picked_SNPs$CHROM == "ctg.000013F", -c(1,2)])
rownames(tmp_hm_freq) <- paste(test_sprat_SNPs_v2$SNP_set$picked_SNPs[test_sprat_SNPs_v2$SNP_set$picked_SNPs$CHROM == "ctg.000013F", 1], test_sprat_SNPs_v2$SNP_set$picked_SNPs[test_sprat_SNPs_v2$SNP_set$picked_SNPs$CHROM == "ctg.000013F", 2], sep = "_")
pdf(file = "~/Projects/Sprat/doc/ctg_000013F_hm.pdf")
heatmap(t(tmp_hm_freq), scale = "none", Colv = NA, margins = c(5,10))
dev.off()

tmp_filter <- Sprat_pool_freq$CHROM == "ctg.000013F" & sprat_site_filter & Sprat_pool_freq$oce_v_fjord > 0.4
#tmp_filter <- Sprat_pool_freq$CHROM == "ctg.000013F" & sprat_site_filter
tmp_hm_freq <- as.matrix(Sprat_pool_freq[tmp_filter, grep("_Afreq", names(Sprat_pool_freq))])
rownames(tmp_hm_freq) <- paste(Sprat_pool_freq[tmp_filter, 1], Sprat_pool_freq[tmp_filter, 2], sep = "_")
#pdf(file = "~/Projects/Sprat/doc/ctg_000013F_full_hm.pdf")
pdf(file = "~/Projects/Sprat/doc/ctg_000013F_high_daf_hm.pdf")
heatmap(t(tmp_hm_freq), scale = "none", Colv = NA, margins = c(5,10))
dev.off()



Sprat_pool_freq_lo <- Sprat_HiC_liftover(sprat_pos_data = Sprat_pool_freq, liftover_df = IPA_pri_v_Ch_v2.0.2_satsuma, chr_size_df = Ch_v2.0.2_sizes, lo_cols = 3:dim(Sprat_pool_freq)[2])
save(Sprat_pool_freq_lo, file = "~/Projects/Sprat/data/PacBio_satsuma/Sprat_pool_freq_lo.RData")

#Adding some site filtering a la Multispecies SNP chip
sprat_opp_fix_filter_lo <- (rowSums(Sprat_pool_freq_lo[,grep("_Afreq", names(Sprat_pool_freq_lo))] == 0, na.rm = T) > 0) & (rowSums(Sprat_pool_freq_lo[,grep("_Afreq", names(Sprat_pool_freq_lo))] == 1, na.rm = T) > 0)
#Attmepting to eliminate SNPs where the large majority of counts are exactly 0 or 1
head(Sprat_pool_freq_lo[sprat_opp_fix_filter_lo,c(1,2,grep("_Afreq", names(Sprat_pool_freq_lo)))])
sprat_pool_dp_filter_lo <- !(rowSums(Sprat_pool_freq_lo[,grep("_DP", names(Sprat_pool_freq_lo))] <= 20, na.rm = T) > 0)
head(Sprat_pool_freq_lo[sprat_pool_dp_filter_lo,])
sprat_site_filter_lo <- !sprat_opp_fix_filter_lo & sprat_pool_dp_filter_lo
save(sprat_site_filter_lo, file = "~/Projects/Sprat/data/PacBio_satsuma/sprat_site_filter_lo.RData" )


#Loading herring regions, for comparison
#sal_p20_regions_GR - older salinity regions
regions_raw <- readRDS("~/Projects/Herring/data/British_populations/Top50_SNPs_Cluhar_72pools_50contrasts_sigBonf.minGapWidth.50K.minWidth.2.minNdAFgt04.1.minNreg.2.topSNPs.mindAF.037.minNdAFgt04.1.RDS")
regions_raw$clean_reg <- sub("([A-Za-z_]+[0-9]+[:][0-9]+[-][0-9]+)[.].+","\\1", regions_raw$region)

mixed_reg_GR <- GenomicRanges::reduce(GRanges(unique(regions_raw$clean_reg)), min.gapwidth = 1e3)


png(filename = "~/Projects/Sprat/doc/Sprat_on_Ch_v2.0.2_Oce_v_Brack.png", height = 1000, width = 2000)
plot(x = Sprat_pool_freq_lo$SNP_cumulative_pos[sprat_site_filter_lo], y = Sprat_pool_freq_lo$oce_v_brack[sprat_site_filter_lo], col = Sprat_pool_freq_lo$col[sprat_site_filter_lo], ylim = c(0,1.1), pch = 16, cex = 0.5, xlab  = "Cumulative position", ylab = "DAF", main = "Oceanic vs Brackish")
sel_reg_x_vec <- sal_p20_regions_GR@ranges@start + Ch_v2.0.2_sizes$offset[match(as.character(sal_p20_regions_GR@seqnames), Ch_v2.0.2_sizes$name)]
segments(x0 = sel_reg_x_vec, x1 = sel_reg_x_vec + sal_p20_regions_GR@ranges@width, y0 = 1.05, col = "red", lwd = 5)
dev.off()

#plot_GR <- mixed_reg_GR 
#plot_GR <- sal_p20_regions_GR
Sprat_lo_manhattan(sprat_data = Sprat_pool_freq_lo, daf_col = "oce_v_brack", sel_reg_GR = mixed_reg_GR, plot_prefix = "~/Projects/Sprat/doc/Oce_v_Brack/Oce_v_Brack_mixed_reg_", site_filt_vec = sprat_site_filter_lo )
Sprat_lo_manhattan(sprat_data = Sprat_pool_freq_lo, daf_col = "oce_v_fjord", sel_reg_GR = mixed_reg_GR, plot_prefix = "~/Projects/Sprat/doc/Oce_v_Fjord/Oce_v_Fjord_mixed_reg_", site_filt_vec = sprat_site_filter_lo )
Sprat_lo_manhattan(sprat_data = Sprat_pool_freq_lo, daf_col = "fjord_v_brack", sel_reg_GR = mixed_reg_GR, plot_prefix = "~/Projects/Sprat/doc/Fjord_v_Brack/Fjord_v_Brack_mixed_reg_", site_filt_vec = sprat_site_filter_lo )



plot_GR <- sal_p20_regions_GR
thresh_val_ob <- qnorm(p = 1-(1/sum(sprat_site_filter_lo)), mean = mean(Sprat_pool_freq_lo$oce_v_brack[sprat_site_filter_lo]), sd  = sd(Sprat_pool_freq_lo$oce_v_brack[sprat_site_filter_lo]))
for(chr in paste0("chr", 1:26)){
  tmp_png_file <- paste0("~/Projects/Sprat/doc/Oce_v_Brack/Oce_v_Brack_",chr, ".png")
  tmp_chr_filter <- Sprat_pool_freq_lo$herring_v2.0.2_seqnames == chr
  png(filename = tmp_png_file, height = 1000, width = 2000)
  plot(x = Sprat_pool_freq_lo$SNP_HiC_pos[sprat_site_filter_lo & tmp_chr_filter], y = Sprat_pool_freq_lo$oce_v_brack[sprat_site_filter_lo & tmp_chr_filter], col = "black", ylim = c(0,1), pch = 16, cex = 0.8, xlab  = "Position", ylab = "DAF", main = chr)
  tmp_reg_filter <- as.character(plot_GR@seqnames) == chr
  if(sum(tmp_reg_filter) > 0) rect(xleft = plot_GR@ranges@start[tmp_reg_filter], xright = plot_GR@ranges@start[tmp_reg_filter] + plot_GR@ranges@width[tmp_reg_filter], ybottom = 0, ytop = 1, col = rgb(1,0,0,0.2), border = NULL)
  abline(h = thresh_val_ob, col = "red")
  if(chr == "chr4") abline(v = c(11217660,11217516), col = "darkorchid") 
  dev.off()
}

thresh_val_of <- qnorm(p = 1-(1/sum(sprat_site_filter_lo)), mean = mean(Sprat_pool_freq_lo$oce_v_fjord[sprat_site_filter_lo]), sd  = sd(Sprat_pool_freq_lo$oce_v_fjord[sprat_site_filter_lo]))
for(chr in paste0("chr", 1:26)){
  tmp_png_file <- paste0("~/Projects/Sprat/doc/Oce_v_Fjord/Oce_v_Fjord_",chr, ".png")
  tmp_chr_filter <- Sprat_pool_freq_lo$herring_v2.0.2_seqnames == chr
  png(filename = tmp_png_file, height = 1000, width = 2000)
  plot(x = Sprat_pool_freq_lo$SNP_HiC_pos[sprat_site_filter_lo & tmp_chr_filter], y = Sprat_pool_freq_lo$oce_v_fjord[sprat_site_filter_lo & tmp_chr_filter], col = "black", ylim = c(0,1), pch = 16, cex = 0.8, xlab  = "Position", ylab = "DAF", main = chr)
  tmp_reg_filter <- as.character(plot_GR@seqnames) == chr
  if(sum(tmp_reg_filter) > 0) rect(xleft = plot_GR@ranges@start[tmp_reg_filter], xright = plot_GR@ranges@start[tmp_reg_filter] + plot_GR@ranges@width[tmp_reg_filter], ybottom = 0, ytop = 1, col = rgb(1,0,0,0.2), border = NULL)
  abline(h = thresh_val_of, col = "red") 
  dev.off()
}

thresh_val_fb <- qnorm(p = 1-(1/sum(sprat_site_filter_lo)), mean = mean(Sprat_pool_freq_lo$fjord_v_brack[sprat_site_filter_lo]), sd  = sd(Sprat_pool_freq_lo$fjord_v_brack[sprat_site_filter_lo]))
for(chr in paste0("chr", 1:26)){
  tmp_png_file <- paste0("~/Projects/Sprat/doc/Fjord_v_brack/Fjord_v_brack_",chr, ".png")
  tmp_chr_filter <- Sprat_pool_freq_lo$herring_v2.0.2_seqnames == chr
  png(filename = tmp_png_file, height = 1000, width = 2000)
  plot(x = Sprat_pool_freq_lo$SNP_HiC_pos[sprat_site_filter_lo & tmp_chr_filter], y = Sprat_pool_freq_lo$fjord_v_brack[sprat_site_filter_lo & tmp_chr_filter], col = "black", ylim = c(0,1), pch = 16, cex = 0.8, xlab  = "Position", ylab = "DAF", main = chr)
  tmp_reg_filter <- as.character(plot_GR@seqnames) == chr
  if(sum(tmp_reg_filter) > 0) rect(xleft = plot_GR@ranges@start[tmp_reg_filter], xright = plot_GR@ranges@start[tmp_reg_filter] + plot_GR@ranges@width[tmp_reg_filter], ybottom = 0, ytop = 1, col = rgb(1,0,0,0.2), border = NULL)
  abline(h = thresh_val_fb, col = "red") 
  dev.off()
}

Sprat_lo_GR <- GRanges(seqnames = Sprat_pool_freq_lo$herring_v2.0.2_seqnames, ranges = IRanges(start = Sprat_pool_freq_lo$SNP_HiC_pos, end = Sprat_pool_freq_lo$SNP_HiC_pos))
Sprat_lo_GR$site_filter <- sprat_site_filter_lo
snps_v_sel_reg <- findOverlaps(Sprat_lo_GR, sal_p20_regions_GR)
Sprat_lo_GR$sel_reg <- F
Sprat_lo_GR$sel_reg[snps_v_sel_reg@from] <- T
Sprat_lo_GR$oce_v_brack <- Sprat_pool_freq_lo$oce_v_brack
boxplot(Sprat_lo_GR$oce_v_brack[Sprat_lo_GR$site_filter]~Sprat_lo_GR$sel_reg[Sprat_lo_GR$site_filter])
t.test(Sprat_lo_GR$oce_v_brack[Sprat_lo_GR$site_filter]~Sprat_lo_GR$sel_reg[Sprat_lo_GR$site_filter])

qnorm(p = 1-(1/sum(Sprat_lo_GR$site_filter)), mean = mean(Sprat_lo_GR$oce_v_brack[Sprat_lo_GR$site_filter]), sd  = sd(Sprat_lo_GR$oce_v_brack[Sprat_lo_GR$site_filter]))
#[1] 0.5826694
plot(y = dnorm(x=seq(from = 0, to  = 1, by = 0.001), mean = mean(Sprat_lo_GR$oce_v_brack[Sprat_lo_GR$site_filter]), sd  = sd(Sprat_lo_GR$oce_v_brack[Sprat_lo_GR$site_filter])), x = seq(from = 0, to  = 1, by = 0.001))
abline(v = 0.5826694, col = "red" ) 


#Additional SNP-chip candidates
chip_filter <- sprat_site_filter & (Sprat_pool_freq$oce_v_fjord > 0.65 | Sprat_pool_freq$oce_v_brack > 0.70 | Sprat_pool_freq$fjord_v_brack > 0.55)
Sprat_contrast_chip <- Sprat_pool_freq[chip_filter, 1:2]
Sprat_contrast_chip$SNP_ID <- paste(Sprat_pool_freq[chip_filter, 1], Sprat_pool_freq[chip_filter, 2], sep = "_")


#Fjord-associated SNPS
fjord_snps <- sprat_site_filter_lo & Sprat_pool_freq_lo$fjord_v_brack > 0.5 & Sprat_pool_freq_lo$oce_v_fjord > 0.5
fjord_snp_chr <- table(Sprat_pool_freq_lo$herring_v2.0.2_seqnames[fjord_snps])
head(fjord_snp_chr[order(fjord_snp_chr, decreasing = T)])
hist(Sprat_pool_freq_lo$SNP_HiC_pos[fjord_snps & Sprat_pool_freq_lo$herring_v2.0.2_seqnames == "chr12"])


#Regional heatmaps
#Putative inversions
tmp_filter <- Sprat_pool_freq_lo$herring_v2.0.2_seqnames == "chr4" & sprat_site_filter_lo & (Sprat_pool_freq_lo$oce_v_brack > 0.6) & Sprat_pool_freq_lo$SNP_HiC_pos > 25e6
simple_sprat_pool_hm(snp_data = Sprat_pool_freq_lo, reg_filter = tmp_filter, pdf_file = "~/Projects/Sprat/doc/reg_heatmaps/Chr4_end.pdf")

tmp_filter <- Sprat_pool_freq_lo$herring_v2.0.2_seqnames == "chr5" & sprat_site_filter_lo & (Sprat_pool_freq_lo$oce_v_brack > 0.6)
simple_sprat_pool_hm(snp_data = Sprat_pool_freq_lo, reg_filter = tmp_filter, pdf_file = "~/Projects/Sprat/doc/reg_heatmaps/Chr5_12_Mb.pdf")

tmp_filter <- Sprat_pool_freq_lo$herring_v2.0.2_seqnames == "chr7" & sprat_site_filter_lo & (Sprat_pool_freq_lo$oce_v_brack > 0.6) & Sprat_pool_freq_lo$SNP_HiC_pos > 22e6
simple_sprat_pool_hm(snp_data = Sprat_pool_freq_lo, reg_filter = tmp_filter, pdf_file = "~/Projects/Sprat/doc/reg_heatmaps/Chr7_end.pdf")

tmp_filter <- Sprat_pool_freq_lo$herring_v2.0.2_seqnames == "chr11" & sprat_site_filter_lo & Sprat_pool_freq_lo$oce_v_brack > 0.4 & Sprat_pool_freq_lo$SNP_HiC_pos < 7e6
simple_sprat_pool_hm(snp_data = Sprat_pool_freq_lo, reg_filter = tmp_filter, pdf_file =  "~/Projects/Sprat/doc/reg_heatmaps/Chr11_5_Mb.pdf")

tmp_filter <- Sprat_pool_freq_lo$herring_v2.0.2_seqnames == "chr12" & sprat_site_filter_lo & (Sprat_pool_freq_lo$oce_v_fjord > 0.35 | Sprat_pool_freq_lo$fjord_v_brack > 0.35) & Sprat_pool_freq_lo$SNP_HiC_pos < 2e6
simple_sprat_pool_hm(snp_data = Sprat_pool_freq_lo, reg_filter = tmp_filter, pdf_file =  "~/Projects/Sprat/doc/reg_heatmaps/Chr12_start.pdf")

tmp_filter <- Sprat_pool_freq_lo$herring_v2.0.2_seqnames == "chr13" & sprat_site_filter_lo & Sprat_pool_freq_lo$oce_v_brack > 0.6 & Sprat_pool_freq_lo$SNP_HiC_pos < 10e6
simple_sprat_pool_hm(snp_data = Sprat_pool_freq_lo, reg_filter = tmp_filter, pdf_file =  "~/Projects/Sprat/doc/reg_heatmaps/Chr13_7_Mb.pdf")

tmp_filter <- Sprat_pool_freq_lo$herring_v2.0.2_seqnames == "chr18" & sprat_site_filter_lo & Sprat_pool_freq_lo$oce_v_brack > 0.7
simple_sprat_pool_hm(snp_data = Sprat_pool_freq_lo, reg_filter = tmp_filter, pdf_file =  "~/Projects/Sprat/doc/reg_heatmaps/Chr18_inv.pdf")

tmp_filter <- Sprat_pool_freq_lo$herring_v2.0.2_seqnames == "chr21" & sprat_site_filter_lo & Sprat_pool_freq_lo$oce_v_brack > 0.6 & Sprat_pool_freq_lo$SNP_HiC_pos < 5e6
simple_sprat_pool_hm(snp_data = Sprat_pool_freq_lo, reg_filter = tmp_filter, pdf_file =  "~/Projects/Sprat/doc/reg_heatmaps/Chr21_3_Mb.pdf")

#Signals overlapping the herring
tmp_filter <- Sprat_pool_freq_lo$herring_v2.0.2_seqnames == "chr3" & sprat_site_filter_lo & Sprat_pool_freq_lo$oce_v_brack > 0.6 & Sprat_pool_freq_lo$SNP_HiC_pos > 30e6 
simple_sprat_pool_hm(snp_data = Sprat_pool_freq_lo, reg_filter = tmp_filter, pdf_file =  "~/Projects/Sprat/doc/reg_heatmaps/Chr6_33_MB.pdf")
Sprat_pool_freq_lo[tmp_filter, 1:30]
IPA_pri[["ctg.000027F"]][529537 + -40:40]

tmp_filter <- Sprat_pool_freq_lo$herring_v2.0.2_seqnames == "chr12" & sprat_site_filter_lo & Sprat_pool_freq_lo$oce_v_brack > 0.4 & Sprat_pool_freq_lo$SNP_HiC_pos > 2e6 & Sprat_pool_freq_lo$SNP_HiC_pos < 3e6
simple_sprat_pool_hm(snp_data = Sprat_pool_freq_lo, reg_filter = tmp_filter, pdf_file =  "~/Projects/Sprat/doc/reg_heatmaps/Chr12_2.3_MB.pdf")
chr12_snp_idx <- Sprat_pool_freq_lo$herring_v2.0.2_seqnames == "chr12" & Sprat_pool_freq_lo$SNP_HiC_pos == 2273814
Sprat_pool_freq_lo[chr12_snp_idx,1:30]
IPA_pri[["ctg.000058F"]][1125866 + -40:40]


tmp_filter <- Sprat_pool_freq_lo$herring_v2.0.2_seqnames == "chr19" & sprat_site_filter_lo & Sprat_pool_freq_lo$oce_v_brack > 0.4 & Sprat_pool_freq_lo$SNP_HiC_pos > 6.3e6 & Sprat_pool_freq_lo$SNP_HiC_pos < 6.5e6
simple_sprat_pool_hm(snp_data = Sprat_pool_freq_lo, reg_filter = tmp_filter, pdf_file =  "~/Projects/Sprat/doc/reg_heatmaps/Chr19_6.3_MB.pdf")
chr19_snp_idx <- Sprat_pool_freq_lo$herring_v2.0.2_seqnames == "chr19" & Sprat_pool_freq_lo$SNP_HiC_pos == 6380297
Sprat_pool_freq_lo[chr19_snp_idx,1:30]
IPA_pri[["ctg.001310F"]][114012 + -20:20]


#Support functions
simple_sprat_pool_hm <- function(snp_data, reg_filter, pdf_file){
  tmp_hm_freq <- as.matrix(snp_data[reg_filter, grep("_Afreq", names(snp_data))])
  rownames(tmp_hm_freq) <- paste(snp_data$herring_v2.0.2_seqnames[reg_filter], snp_data$SNP_HiC_pos[reg_filter], sep = "_")
  pdf(file = pdf_file)
  heatmap(t(tmp_hm_freq), scale = "none", Colv = NA, margins = c(5,10))
  dev.off()
}

Sprat_lo_manhattan <- function(sprat_data, daf_col, site_filt_vec, sel_reg_GR, plot_prefix){
  daf_thresh_val <- qnorm(p = 1-(1/sum(site_filt_vec)), mean = mean(sprat_data[site_filt_vec, daf_col]), sd  = sd(sprat_data[site_filt_vec, daf_col]))
  for(chr in paste0("chr", 1:26)){
    tmp_png_file <- paste0(plot_prefix,chr, ".png")
    tmp_chr_filter <- sprat_data$herring_v2.0.2_seqnames == chr
    png(filename = tmp_png_file, height = 1000, width = 2000)
    plot(x = sprat_data$SNP_HiC_pos[site_filt_vec & tmp_chr_filter], y = sprat_data[site_filt_vec & tmp_chr_filter, daf_col], col = "black", ylim = c(0,1), pch = 16, cex = 0.8, xlab  = "Position", ylab = "DAF", main = chr)
    tmp_reg_filter <- as.character(sel_reg_GR@seqnames) == chr
    if(sum(tmp_reg_filter) > 0) rect(xleft = sel_reg_GR@ranges@start[tmp_reg_filter], xright = sel_reg_GR@ranges@start[tmp_reg_filter] + sel_reg_GR@ranges@width[tmp_reg_filter], ybottom = 0, ytop = 1, col = rgb(1,0,0,0.2), border = NULL)
    abline(h = daf_thresh_val, col = "red")
    #if(chr == "chr4") abline(v = c(11217660,11217516), col = "darkorchid") 
    dev.off()
  }
}

##TBD!!!
Sprat_ctg_plot <- function(target_ctg, spart_freq_data = Sprat_pool_freq, pdf_file ){
  tmp_filter <- spart_freq_data$CHROM == target_ctg & sprat_site_filter  
  png(filename = paste0("~/Projects/Sprat/doc/oce_v_brack_", target_ctg,".png"), width = 1000)
  plot(x = spart_freq_data$POS[tmp_filter], y = spart_freq_data$oce_v_brack[tmp_filter], pch = 16, cex = 0.3, ylab = "DAF")
  lines(x = spart_freq_data$POS[tmp_filter], y = stats::filter(spart_freq_data$oce_v_brack[tmp_filter], filter = rep(1/100, 100)), col = "red")
  dev.off()
  
  tmp_filter <- spart_freq_data$CHROM == "ctg.000013F" & sprat_site_filter & Sprat_pool_freq$oce_v_fjord > 0.4
  #tmp_filter <- Sprat_pool_freq$CHROM == "ctg.000013F" & sprat_site_filter
  tmp_hm_freq <- as.matrix(spart_freq_data[tmp_filter, grep("_Afreq", names(spart_freq_data))])
  rownames(tmp_hm_freq) <- paste(spart_freq_data[tmp_filter, 1], spart_freq_data[tmp_filter, 2], sep = "_")
  #pdf(file = "~/Projects/Sprat/doc/ctg_000013F_full_hm.pdf")
  pdf(file = pdf_file)
  heatmap(t(tmp_hm_freq), scale = "none", Colv = NA, margins = c(5,10))
  dev.off()
}
##TBD!!!