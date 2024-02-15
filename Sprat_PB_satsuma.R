#Author: Mats Pettersson
#mats.pettersson@imbim.uu.se


IPA_pri_satsuma <- read.table("~/Projects/Sprat/data/PacBio_satsuma/Sprat_IPA_primary_v_Ch_v2.0.2.chained.out_sorted.gz", stringsAsFactors=F, sep = "\t", comment.char = "")
IPA_pri_v_Ch_v2.0.2_satsuma <- satsuma_processing_v2("~/Projects/Sprat/data/PacBio_satsuma/Sprat_IPA_primary_v_Ch_v2.0.2.chained.out_sorted.gz")
IPA_pri_v_Ch_v2.0.2_satsuma <- satsuma_direction_est(IPA_pri_v_Ch_v2.0.2_satsuma)

names(IPA_pri_v_Ch_v2.0.2_satsuma)[3:6] <- paste0("sprat_IPA_", names(IPA_pri_v_Ch_v2.0.2_satsuma)[3:6])
names(IPA_pri_v_Ch_v2.0.2_satsuma)[8:11] <- paste0("herring_v2.0.2_", names(IPA_pri_v_Ch_v2.0.2_satsuma)[8:11])
names(IPA_pri_v_Ch_v2.0.2_satsuma)[8:11] <- sub("T_", "", names(IPA_pri_v_Ch_v2.0.2_satsuma)[8:11])

save(IPA_pri_v_Ch_v2.0.2_satsuma, file = "~/Projects/Sprat/data/PacBio_satsuma/Sprat_IPA_primary_v_Ch_v2.0.2_processed.RData")

IPA_pri_chr8 <- IPA_pri_satsuma[IPA_pri_satsuma[,4] == "chr8",]



IPA_alt_satsuma <- read.table("~/Projects/Sprat/data/PacBio_satsuma/Sprat_IPA_alternative_v_Ch_v2.0.2.chained.out_sorted", stringsAsFactors=F, sep = "\t", comment.char = "")
IPA_alt_chr8 <- IPA_alt_satsuma[IPA_alt_satsuma[,4] == "chr8",]

IPA_pri_v_chrY_satsuma <- read.table("~/Projects/Sprat/data/PacBio_satsuma/Sprat_IPA_pri_v_Ch_chrY_satsuma.out_sorted", stringsAsFactors=F, sep = "\t", comment.char = "")
IPA_pri_v_chrY <- IPA_pri_v_chrY_satsuma[IPA_pri_v_chrY_satsuma[,4] == "chr8_Y",]

IPA_alt_v_chrY_satsuma <- read.table("~/Projects/Sprat/data/PacBio_satsuma/Sprat_IPA_alt_v_Ch_chrY_satsuma.out_sorted", stringsAsFactors=F, sep = "\t", comment.char = "")
IPA_alt_v_chrY <- IPA_alt_v_chrY_satsuma[IPA_alt_v_chrY_satsuma[,4] == "chr8_Y",]

#HiC assmeblies
#Sprat_3D_v_Ch_v2.0.2 <- read_and_plot_fish_satsuma2("~/Projects/Sprat/data/HiC_satsuma/Sprat_3d_DNA_v_Ch_v2.0.2_satsuma.out_sorted", fish_name = "Sprat_3D_Ch_v2.0.2", plot_dir = "~/Projects/Sprat/doc/")
#Sprat_3D_v_Ch_v2.0.2_syn <- satsuma_synteny_blocks2("~/Projects/Sprat/data/HiC_satsuma/Sprat_3d_DNA_v_Ch_v2.0.2_satsuma.out_sorted", "~/Projects/Sprat/doc/Sprat_3D_Ch_v2.0.2_synteny.pdf", hit_co = 20000)
Sprat_Salsa_v_Ch_v2.0.2 <- read_and_plot_fish_satsuma2("~/Projects/Sprat/data/HiC_satsuma/Sprat_Salsa_v_Ch_v2.0.2_satsuma.out_sorted", fish_name = "Sprat_Salsa_Ch_v2.0.2", plot_dir = "~/Projects/Sprat/doc/")
Sprat_Salsa_v_Ch_v2.0.2_syn <- satsuma_synteny_blocks2("~/Projects/Sprat/data/HiC_satsuma/Sprat_Salsa_v_Ch_v2.0.2_satsuma.out_sorted", "~/Projects/Sprat/doc/Sprat_Salsa_Ch_v2.0.2_synteny.pdf", hit_co = 20000)

Sprat_3Dv2_v_Ch_v2.0.2 <- read_and_plot_fish_satsuma2("~/Projects/Sprat/data/HiC_satsuma/Sprat_3d_DNAv2_v_Ch_v2.0.2_satsuma.out_sorted", fish_name = "Sprat_3Dv2_Ch_v2.0.2", plot_dir = "~/Projects/Sprat/doc/")
Sprat_3Dv3_v_Ch_v2.0.2_syn <- satsuma_synteny_blocks2("~/Projects/Sprat/data/HiC_satsuma/Sprat_3d_DNAv2_v_Ch_v2.0.2_satsuma.out_sorted", "~/Projects/Sprat/doc/Sprat_3Dv2_Ch_v2.0.2_synteny.pdf", hit_co = 20000)

#Chr12 inversion
Sprat_3Dv2_v_Ch_chr12 <- Sprat_3Dv2_v_Ch_v2.0.2[Sprat_3Dv2_v_Ch_v2.0.2$V4 == "chr12",] # Sprat_3Dv2_v_Ch_v2.0.2$V6 > 4.5e6 & Sprat_3Dv2_v_Ch_v2.0.2$V5 < 25.9e6
Sprat_3Dv2_v_Ch_chr12$col <- "grey50"
Sprat_3Dv2_v_Ch_chr12$col[Sprat_3Dv2_v_Ch_chr12$V1 == "HiC_scaffold_19"] <- "red"
Sprat_3Dv2_v_Ch_chr12$col[Sprat_3Dv2_v_Ch_chr12$V1 == "HiC_scaffold_9"] <- "blue"
Sprat_3Dv2_v_Ch_chr12$col[Sprat_3Dv2_v_Ch_chr12$V1 == "HiC_scaffold_10"] <- "darkorchid"
Sprat_3Dv2_v_Ch_chr12$cex <- 0.1
Sprat_3Dv2_v_Ch_chr12$cex[Sprat_3Dv2_v_Ch_chr12$V1 %in% c("HiC_scaffold_19","HiC_scaffold_9","HiC_scaffold_10")] <-1



png("~/Projects/Sprat/doc/Herring_chr12_v_Sprat_3Dv2.png", height = 1000, width = 1000)
plot(x = Sprat_3Dv2_v_Ch_chr12$V5, y = Sprat_3Dv2_v_Ch_chr12$V2, col = Sprat_3Dv2_v_Ch_chr12$col, cex = Sprat_3Dv2_v_Ch_chr12$cex, pch = 20, ylab = "Sprat HiC position", xlab = "Herring Chr 12 position") 
abline(v = 17.826e6)
abline(v = 25.60e6)
#abline(h = 42.0e6)
#abline(h = 48.5e6)
legend(x = "right", pch = 20, cex = 1.2, col = c("red", "blue", "darkorchid","grey50"), legend = c("HiC_scaffold_19", "HiC_scaffold_9", "HiC_scaffold_10", "other"))
dev.off()

#Gap on Sprat HiC_scaffold_9 - not so interesting mostly an extra box of the flanking repeats, it seems, see above.
#Sprat_HiC9_v_Ch_v2.0.2 <- Sprat_3Dv2_v_Ch_v2.0.2[Sprat_3Dv2_v_Ch_v2.0.2$V4 == "chr12" & Sprat_3Dv2_v_Ch_v2.0.2$V1 == "HiC_scaffold_9" & Sprat_3Dv2_v_Ch_v2.0.2$V3 > 42.0e6 & Sprat_3Dv2_v_Ch_v2.0.2$V2 < 48.0e6,]
#Sprat_HiC9_v_Ch_v2.0.2$col <- "black"
#Sprat_HiC9_v_Ch_v2.0.2$col[Sprat_HiC9_v_Ch_v2.0.2$V1 == "HiC_scaffold_19"] <- "red"
#Sprat_HiC9_v_Ch_v2.0.2$col[Sprat_HiC9_v_Ch_v2.0.2$V1 == "HiC_scaffold_9"] <- "blue"

#plot(x = Sprat_HiC9_v_Ch_v2.0.2$V5, y = Sprat_HiC9_v_Ch_v2.0.2$V2, col = Sprat_HiC9_v_Ch_v2.0.2$col)

#Chr6 inversion
Sprat_3Dv2_v_Ch_chr6 <- Sprat_3Dv2_v_Ch_v2.0.2[Sprat_3Dv2_v_Ch_v2.0.2$V4 == "chr6",] # & Sprat_3Dv2_v_Ch_v2.0.2$V6 > 20.0e6 & Sprat_3Dv2_v_Ch_v2.0.2$V5 < 26.0e6
Sprat_3Dv2_v_Ch_chr6$col <- "grey50"
Sprat_3Dv2_v_Ch_chr6$cex <- 0.3

png("~/Projects/Sprat/doc/Herring_chr6_v_Sprat_3Dv2.png", height = 1000, width = 1000)
plot(x = Sprat_3Dv2_v_Ch_chr6$V5, y = Sprat_3Dv2_v_Ch_chr6$V2, col = Sprat_3Dv2_v_Ch_chr6$col, cex = Sprat_3Dv2_v_Ch_chr6$cex, pch = 20, ylab = "Sprat HiC position", xlab = "Herring Chr 6 position") 
abline(v = 22.283e6)
abline(v = 24.869e6)
dev.off()

#Chr17 inversion
Sprat_3Dv2_v_Ch_chr17 <- Sprat_3Dv2_v_Ch_v2.0.2[Sprat_3Dv2_v_Ch_v2.0.2$V4 == "chr17",] # & Sprat_3Dv2_v_Ch_v2.0.2$V6 > 20.0e6 & Sprat_3Dv2_v_Ch_v2.0.2$V5 < 26.0e6
Sprat_3Dv2_v_Ch_chr17$col <- "grey50"
Sprat_3Dv2_v_Ch_chr17$cex <- 0.3

png("~/Projects/Sprat/doc/Herring_chr17_v_Sprat_3Dv2.png", height = 1000, width = 1000)
plot(x = Sprat_3Dv2_v_Ch_chr17$V5, y = Sprat_3Dv2_v_Ch_chr17$V2, col = Sprat_3Dv2_v_Ch_chr17$col, cex = Sprat_3Dv2_v_Ch_chr17$cex, pch = 20, ylab = "Sprat HiC position", xlab = "Herring Chr 17 position") 
abline(v = 25.805e6)
abline(v = 27.569e6)
dev.off()

#Chr23 inversion
Sprat_3Dv2_v_Ch_chr23 <- Sprat_3Dv2_v_Ch_v2.0.2[Sprat_3Dv2_v_Ch_v2.0.2$V4 == "chr23",] # & Sprat_3Dv2_v_Ch_v2.0.2$V6 > 20.0e6 & Sprat_3Dv2_v_Ch_v2.0.2$V5 < 26.0e6
Sprat_3Dv2_v_Ch_chr23$col <- "grey50"
Sprat_3Dv2_v_Ch_chr23$cex <- 0.3

png("~/Projects/Sprat/doc/Herring_chr23_v_Sprat_3Dv2.png", height = 1000, width = 1000)
plot(x = Sprat_3Dv2_v_Ch_chr23$V5, y = Sprat_3Dv2_v_Ch_chr23$V2, col = Sprat_3Dv2_v_Ch_chr23$col, cex = Sprat_3Dv2_v_Ch_chr23$cex, pch = 20, ylab = "Sprat HiC position", xlab = "Herring Chr 23 position") 
abline(v = 16.226e6)
abline(v = 17.604e6)
dev.off()

#satsuma_target_GR <- GRanges(seqnames=sub_satsuma[,14], ranges=IRanges(start=satsuma_in[,5], end = satsuma_in[,6]))
plot(y = IPA_pri_chr8[,2], x = IPA_pri_chr8[,5], xlab = "Herring position", ylab = "Query position", type = "n", main = "Sprat IPA pri vs Herring chr8", xlim = c(0, 3.5e7))
points(y = IPA_pri_chr8[,2], x = IPA_pri_chr8[,5], pch = 16, col = as.integer(as.factor(IPA_pri_chr8[,1])) %% 8 + 1)

plot(y = IPA_pri_chr8[,2], x = IPA_pri_chr8[,5], xlab = "Herring position", ylab = "Query position", type = "n", main = "Sprat IPA pri vs Herring chr8", xlim = c(1.95e7, 2.20e7))
points(y = IPA_pri_chr8[,2], x = IPA_pri_chr8[,5], pch = 16, col = as.integer(as.factor(IPA_pri_chr8[,1])) %% 8 + 1)

plot(y = IPA_alt_chr8[,2], x = IPA_alt_chr8[,5], xlab = "Herring position", ylab = "Query position", type = "n", main = "Sprat IPA alt vs Herring chr8", xlim = c(1.95e7, 2.20e7))
points(y = IPA_alt_chr8[,2], x = IPA_alt_chr8[,5], pch = 16, col = as.integer(as.factor(IPA_alt_chr8[,1])) %% 8 + 1)

plot(y = IPA_pri_v_chrY[,2], x = IPA_pri_v_chrY[,5], xlab = "Herring position", ylab = "Query position", type = "n", main = "Sprat IPA pri vs Herring chr8_Y", xlim = c(1.95e7, 2.20e7))
points(y = IPA_pri_v_chrY[,2], x = IPA_pri_v_chrY[,5], pch = 16, col = as.integer(as.factor(IPA_pri_v_chrY[,1])) %% 8 + 1)

plot(y = IPA_alt_v_chrY[,2], x = IPA_alt_v_chrY[,5], xlab = "Herring position", ylab = "Query position", type = "n", main = "Sprat IPA alt vs Herring chr8_Y", xlim = c(1.95e7, 2.20e7))
points(y = IPA_alt_v_chrY[,2], x = IPA_alt_v_chrY[,5], pch = 16, col = as.integer(as.factor(IPA_alt_v_chrY[,1])) %% 8 + 1)


table(IPA_pri_chr8[IPA_pri_chr8[,5] > 21.0e6 & IPA_pri_chr8[,5] < 21.4e6,1])

pdf(file = "~/Projects/Sprat/doc/Sprat_IPA_v_SDR.pdf", width = 10)
plot(y = IPA_pri_chr8[,2], x = IPA_pri_chr8[,5], xlab = "Herring position", ylab = "Query position", type = "n", main = "Sprat IPA pri vs Herring chr8_X", xlim = c(2.05e7, 2.16e7), ylim = c(0, 7.5e5))
points(y = IPA_pri_chr8[IPA_pri_chr8[,1] == "ctg.000592F",2], x = IPA_pri_chr8[IPA_pri_chr8[,1] == "ctg.000592F",5], pch = 16, col = "darkorchid")
points(y = IPA_pri_chr8[IPA_pri_chr8[,1] != "ctg.000592F",2], x = IPA_pri_chr8[IPA_pri_chr8[,1] != "ctg.000592F",5], pch = 16, col = "grey50", cex = 0.7)

plot(y = IPA_pri_v_chrY[,2], x = IPA_pri_v_chrY[,5], xlab = "Herring position", ylab = "Query position", type = "n", main = "Sprat IPA pri vs Herring chr8_Y", xlim = c(2.05e7, 2.17e7), ylim = c(0, 7.5e5))
points(y = IPA_pri_v_chrY[IPA_pri_v_chrY[,1] == "ctg.000592F",2], x = IPA_pri_v_chrY[IPA_pri_v_chrY[,1] == "ctg.000592F",5], pch = 16, col = "darkorchid")
points(y = IPA_pri_v_chrY[IPA_pri_v_chrY[,1] != "ctg.000592F",2], x = IPA_pri_v_chrY[IPA_pri_v_chrY[,1] != "ctg.000592F",5], pch = 16, col = "grey50", cex = 0.7)

plot(y = IPA_alt_chr8[,2], x = IPA_alt_chr8[,5], xlab = "Herring position", ylab = "Query position", type = "n", main = "Sprat IPA alt vs Herring chr8_X", xlim = c(2.05e7, 2.16e7), ylim = c(0, 7.5e5))
points(y = IPA_alt_chr8[IPA_alt_chr8[,1] == "hap_ctg.001096F_HAPLOTIG",2], x = IPA_alt_chr8[IPA_alt_chr8[,1] == "hap_ctg.001096F_HAPLOTIG",5], pch = 16, col = "darkorchid")
points(y = IPA_alt_chr8[IPA_alt_chr8[,1] == "hap_ctg.001499F_HAPLOTIG",2], x = IPA_alt_chr8[IPA_alt_chr8[,1] == "hap_ctg.001499F_HAPLOTIG",5], pch = 16, col = "steelblue")
points(y = IPA_alt_chr8[!(IPA_alt_chr8[,1] %in% c("hap_ctg.001096F_HAPLOTIG", "hap_ctg.001499F_HAPLOTIG")),2], x = IPA_alt_chr8[!(IPA_alt_chr8[,1] %in% c("hap_ctg.001096F_HAPLOTIG", "hap_ctg.001499F_HAPLOTIG")),5], pch = 16, col = "grey50", cex = 0.7)

plot(y = IPA_alt_v_chrY[,2], x = IPA_alt_v_chrY[,5], xlab = "Herring position", ylab = "Query position", type = "n", main = "Sprat IPA alt vs Herring chr8_Y", xlim = c(2.05e7, 2.17e7), ylim = c(0, 7.5e5))
points(y = IPA_alt_v_chrY[IPA_alt_v_chrY[,1] == "hap_ctg.001096F_HAPLOTIG",2], x = IPA_alt_v_chrY[IPA_alt_v_chrY[,1] == "hap_ctg.001096F_HAPLOTIG",5], pch = 16, col = "darkorchid")
points(y = IPA_alt_v_chrY[IPA_alt_v_chrY[,1] == "hap_ctg.001499F_HAPLOTIG",2], x = IPA_alt_v_chrY[IPA_alt_v_chrY[,1] == "hap_ctg.001499F_HAPLOTIG",5], pch = 16, col = "steelblue")
points(y = IPA_alt_v_chrY[!(IPA_alt_v_chrY[,1] %in% c("hap_ctg.001096F_HAPLOTIG", "hap_ctg.001499F_HAPLOTIG")),2], x = IPA_alt_v_chrY[!(IPA_alt_v_chrY[,1] %in% c("hap_ctg.001096F_HAPLOTIG", "hap_ctg.001499F_HAPLOTIG")),5], pch = 16, col = "grey50", cex = 0.7)
dev.off()

#Blast of Herring SDR vs Sprat assemblies
require(Biostrings)
chr8_Y <- readDNAStringSet("/Users/mapet205/Projects/Herring/data/SexDetermination/Ch2.0.2_CHR_with_chr8_Y.fa.gz")
writeXStringSet(x = subseq(chr8_Y["chr8_Y"], start = 20.8e6, end = 21.6e6), file = "~/Projects/Sprat/data/Herring_SDR.fa")
#blastn -task blastn -query ./Herring_SDR.fa -outfmt 6 -culling_limit 5 -db ../assemblies/IPA/Sprat_IPA -out Herring_SDR_v_Sprat_IPA.blastout
SDR_blastout <- read.table("~/Projects/Sprat/data/SDR/Herring_SDR_v_Sprat_IPA.blastout", sep = "\t", stringsAsFactors = F)

plot(x = range(SDR_blastout[,7:8]), y = range(SDR_blastout[,9:10]), type = "n", ylim = c(0, 1e6))
hit_filter <- SDR_blastout$V3 > 80 & SDR_blastout$V4 > 100
segments(x0 =SDR_blastout[hit_filter,7], x1 = SDR_blastout[hit_filter,8], y0 = SDR_blastout[hit_filter,10], y1 = SDR_blastout[hit_filter,9])

####
IPA_pri <- readDNAStringSet("/Users/mapet205/Projects/Sprat/data/assemblies/IPA/final.p_ctg.fasta.gz")

#IPA_pri assembly stats
IPA_pri_w <- width(IPA_pri)
IPA_pri_w <- IPA_pri_w[order(IPA_pri_w, decreasing = T)]
IPA_pri_50 <- sum(width(IPA_pri))/2
IPA_pri_c50 <- min(which(cumsum(IPA_pri_w) > IPA_pri_50))
IPA_pri_w[IPA_pri_c50]


IPA_alt <- readDNAStringSet("/Users/mapet205/Projects/Sprat/data/assemblies/IPA/final.a_ctg.fasta.gz")
HiC_3D_DNA <- readDNAStringSet("/Users/mapet205/Projects/Sprat/data/assemblies/3D_DNA/sprat_3DDNA_20210222.fasta.gz")

summary(IPA_pri@ranges@width)

writeXStringSet(x = subseq(IPA_pri["ctg.000003F"], start = 1, end = 1.5e5), file = "~/Projects/Sprat/data/ctg.000003F_0_150kb.fa")

#Testing the lift-over function
Sprat_pool_freq_lo <- Sprat_HiC_liftover(sprat_pos_data = Sprat_pool_freq, liftover_df = IPA_pri_v_Ch_v2.0.2_satsuma, chr_size_df = Ch_v2.0.2_sizes, lo_cols = 3:dim(Sprat_pool_freq)[2])
save(Sprat_pool_freq_lo, file = "~/Projects/Sprat/data/PacBio_satsuma/Sprat_pool_freq_lo.RData")
png(filename = "~/Projects/Sprat/doc/Sprat_IPA_on_Ch_v2.0.2_manhattan.png", height = 1000, width = 2000)
plot(x = Sprat_pool_freq_lo$SNP_cumulative_pos, y = Sprat_pool_freq_lo$oce_v_brack, col = Sprat_pool_freq_lo$col, pch = 16, cex = 0.5)
dev.off()


#Tree comparing Sprat to Herring, over one example block
#imbim36-193:tree mapet205$ samtools faidx ~/Projects/Herring/data/v2.0.2_reference/Ch_v2.0.2.fasta chr8:27700000-27800000 | bcftools consensus -H 1 --sample BM15_HastKar_Baltic_Spring  ~/Projects/Herring/data/v2.0.2_genotypes/phased_79_ind_v2.0.2.vcf.gz > BM15_HastKar_Baltic_Spring_chr8_27.7Mb.fa
#Applied 222 variants
#imbim36-193:tree mapet205$ samtools faidx ~/Projects/Herring/data/v2.0.2_reference/Ch_v2.0.2.fasta chr8:27700000-27800000 | bcftools consensus -H 1 --sample Gavle54_Gavle_Baltic_Autumn  ~/Projects/Herring/data/v2.0.2_genotypes/phased_79_ind_v2.0.2.vcf.gz > Gavle54_Gavle_Baltic_Autumn_chr8_27.7Mb.fa
#Applied 253 variants
#imbim36-193:tree mapet205$ samtools faidx ~/Projects/Herring/data/v2.0.2_reference/Ch_v2.0.2.fasta chr8:27700000-27800000 | bcftools consensus -H 1 --sample NSSH34_Norway_Atlantic_Spring  ~/Projects/Herring/data/v2.0.2_genotypes/phased_79_ind_v2.0.2.vcf.gz > NSSH34_Norway_Atlantic_Spring_chr8_27.7Mb.fa
#Applied 284 variants
#imbim36-193:tree mapet205$ samtools faidx ~/Projects/Herring/data/v2.0.2_reference/Ch_v2.0.2.fasta chr8:27700000-27800000 | bcftools consensus -H 1 --sample Z14_IsleofMan_Atlantic_Autumn  ~/Projects/Herring/data/v2.0.2_genotypes/phased_79_ind_v2.0.2.vcf.gz > Z14_IsleofMan_Atlantic_Autumn_chr8_27.7Mb.fa
#Applied 233 variants
#imbim36-193:tree mapet205$ samtools faidx ~/Projects/Herring/data/v2.0.2_reference/Ch_v2.0.2.fasta chr8:27700000-27800000 | bcftools consensus -H 1 --sample HWS22_PechoraSea_BarentsSea  ~/Projects/Herring/data/v2.0.2_genotypes/phased_79_ind_v2.0.2.vcf.gz > HWS22_PechoraSea_BarentsSea_chr8_27.7Mb.fa
#Applied 378 variants
#imbim36-193:tree mapet205$ samtools faidx ~/Projects/Herring/data/v2.0.2_reference/Ch_v2.0.2.fasta chr8:27700000-27800000 | bcftools consensus -H 1 --sample HWS12_Japan_SeaofJapan  ~/Projects/Herring/data/v2.0.2_genotypes/phased_79_ind_v2.0.2.vcf.gz > HWS12_Japan_SeaofJapan_chr8_27.7Mb.fa
#Applied 377 variants
#imbim36-193:tree mapet205$ samtools faidx ~/Projects/Herring/data/v2.0.2_reference/Ch_v2.0.2.fasta chr8:27700000-27800000 | bcftools consensus -H 1 --sample HWS32_WhiteSea_WhiteSea  ~/Projects/Herring/data/v2.0.2_genotypes/phased_79_ind_v2.0.2.vcf.gz > HWS32_WhiteSea_WhiteSea_chr8_27.7Mb.fa
#Applied 370 variants
#imbim36-193:tree mapet205$ samtools faidx ~/Projects/Herring/data/v2.0.2_reference/Ch_v2.0.2.fasta chr8:27700000-27800000 | bcftools consensus -H 1 --sample HWS42_KandalakshaBay_WhiteSea_Spring  ~/Projects/Herring/data/v2.0.2_genotypes/phased_79_ind_v2.0.2.vcf.gz > HWS42_KandalakshaBay_WhiteSea_Spring_chr8_27.7Mb.fa
#pplied 336 variants
#imbim36-193:tree mapet205$ samtools faidx ~/Projects/Herring/data/v2.0.2_reference/Ch_v2.0.2.fasta chr8:27700000-27800000 | bcftools consensus -H 1 --sample HWS53_KandalakshaBay_WhiteSea_Summer  ~/Projects/Herring/data/v2.0.2_genotypes/phased_79_ind_v2.0.2.vcf.gz > HWS53_KandalakshaBay_WhiteSea_Summer_chr8_27.7Mb.fa
#Applied 316 variants
#imbim36-193:tree mapet205$ samtools faidx ~/Projects/Herring/data/v2.0.2_reference/Ch_v2.0.2.fasta chr8:27700000-27800000 | bcftools consensus -H 1 --sample HWS63_Balsfjord_Atlantic  ~/Projects/Herring/data/v2.0.2_genotypes/phased_79_ind_v2.0.2.vcf.gz > HWS63_Balsfjord_Atlantic_chr8_27.7Mb.fa
#Applied 291 variants
#imbim36-193:tree mapet205$ samtools faidx ~/Projects/Herring/data/v2.0.2_reference/Ch_v2.0.2.fasta chr8:27700000-27800000 | bcftools consensus -H 1 --sample Pacific3_Vancouver_Pacific  ~/Projects/Herring/data/v2.0.2_genotypes/phased_79_ind_v2.0.2.vcf.gz > Pacific3_Vancouver_Pacific_chr8_27.7Mb.fa
#Applied 405 variants

BM15_chr8_27.7Mb <- read.FASTA("~/Projects/Sprat/data/tree/BM15_HastKar_Baltic_Spring_chr8_27.7Mb.fa")
HWS22_chr8_27.7Mb <- read.FASTA("~/Projects/Sprat/data/tree/HWS22_PechoraSea_BarentsSea_chr8_27.7Mb.fa")
Pacific3_chr8_27.7Mb <- read.FASTA("~/Projects/Sprat/data/tree/Pacific3_Vancouver_Pacific_chr8_27.7Mb.fa")
NSSH34_chr8_27.7Mb <- read.FASTA("~/Projects/Sprat/data/tree/NSSH34_Norway_Atlantic_Spring_chr8_27.7Mb.fa")

Sprat_ctg.000003F_0_150kb <- read.FASTA("~/Projects/Sprat/data/ctg.000003F_0_150kb.fa")


chr8_27.7Mb_DNA <- c(HWS22_chr8_27.7Mb,BM15_chr8_27.7Mb, Pacific3_chr8_27.7Mb, NSSH34_chr8_27.7Mb, Sprat_ctg.000003F_0_150kb)
chr8_27.7Mb_clustalo <- clustalomega(x = chr8_27.7Mb_DNA)
save(chr8_27.7Mb_clustalo, file = "~/Projects/Sprat/data/tree/chr8_27.7Mb_clustalo.RData")

#Smaller  subsections
Sprat_ctg.000003F_0_25kb <- Sprat_ctg.000003F_0_150kb
Sprat_ctg.000003F_0_25kb[[1]] <- Sprat_ctg.000003F_0_150kb[[1]][1:25000]

BM15_chr8_27.7Mb_0_25kb <- BM15_chr8_27.7Mb
BM15_chr8_27.7Mb_0_25kb[[1]] <- BM15_chr8_27.7Mb[[1]][1:25000]

HWS22_chr8_27.7Mb_0_25kb <- HWS22_chr8_27.7Mb
HWS22_chr8_27.7Mb_0_25kb[[1]] <- HWS22_chr8_27.7Mb[[1]][1:25000]

Pacific3_chr8_27.7Mb_0_25kb <- Pacific3_chr8_27.7Mb
Pacific3_chr8_27.7Mb_0_25kb[[1]] <- Pacific3_chr8_27.7Mb[[1]][1:25000]

NSSH34_chr8_27.7Mb_0_25kb <- NSSH34_chr8_27.7Mb
NSSH34_chr8_27.7Mb_0_25kb[[1]] <- NSSH34_chr8_27.7Mb[[1]][1:25000]

chr8_27.7Mb_DNA_0_25kb <- c(HWS22_chr8_27.7Mb_0_25kb,BM15_chr8_27.7Mb_0_25kb, Pacific3_chr8_27.7Mb_0_25kb, NSSH34_chr8_27.7Mb_0_25kb, Sprat_ctg.000003F_0_25kb)
chr8_27.7Mb_clustalo_0_25kb <- clustalomega(x = chr8_27.7Mb_DNA_0_25kb)
save(chr8_27.7Mb_clustalo_0_25kb, file = "~/Projects/Sprat/data/tree/chr8_27.7Mb_0_25kb_clustalo.RData")

Sprat_ctg.000003F_75_120kb <- Sprat_ctg.000003F_0_150kb
Sprat_ctg.000003F_75_120kb[[1]] <- Sprat_ctg.000003F_0_150kb[[1]][75001:120000]

BM15_chr8_27.7Mb_75_120kb <- BM15_chr8_27.7Mb
BM15_chr8_27.7Mb_75_120kb[[1]] <- BM15_chr8_27.7Mb[[1]][70001:100000]

HWS22_chr8_27.7Mb_75_120kb <- HWS22_chr8_27.7Mb
HWS22_chr8_27.7Mb_75_120kb[[1]] <- HWS22_chr8_27.7Mb[[1]][70001:100000]

Pacific3_chr8_27.7Mb_75_120kb <- Pacific3_chr8_27.7Mb
Pacific3_chr8_27.7Mb_75_120kb[[1]] <- Pacific3_chr8_27.7Mb[[1]][70001:100000]

NSSH34_chr8_27.7Mb_75_120kb <- NSSH34_chr8_27.7Mb
NSSH34_chr8_27.7Mb_75_120kb[[1]] <- NSSH34_chr8_27.7Mb[[1]][70001:100000]

chr8_27.7Mb_DNA_75_120kb <- c(HWS22_chr8_27.7Mb_75_120kb,BM15_chr8_27.7Mb_75_120kb, Pacific3_chr8_27.7Mb_75_120kb, NSSH34_chr8_27.7Mb_75_120kb, Sprat_ctg.000003F_75_120kb)
chr8_27.7Mb_clustalo_75_120kb <- clustalomega(x = chr8_27.7Mb_DNA_75_120kb)
save(chr8_27.7Mb_clustalo_75_120kb, file = "~/Projects/Sprat/data/tree/chr8_27.7Mb_75_120kb_clustalo.RData")

#Example regions:
table(IPA_pri_satsuma[IPA_pri_satsuma$V4 == "chr1" & IPA_pri_satsuma$V5 > 1.00e7 & IPA_pri_satsuma$V5 < 1.01e7, "V1"])
chr1_10.0_Mb <- Sprat_v_herring_aln(herring_GR = GRanges("chr1:10.02e6-10.07e6"), sprat_GR = GRanges("ctg.002579F:2.3e4-7e4"))
chr1_10.0_Mb_clustalo <- clustalomega(x = chr1_10.0_Mb)
save(chr1_10.0_Mb_clustalo, file = "~/Projects/Sprat/data/tree/chr1_10.0_Mb_clustalo.RData")
checkAlignment(chr1_10.0_Mb_clustalo)
SvH_tree <- bionj(dist.dna(chr1_10.0_Mb_clustalo[,2000:18000]))
SvH_tree_rooted <- root(SvH_tree, outgroup = "ctg.002579F", resolve.root = T)

reg_filter <- IPA_pri_satsuma$V4 == "chr13" & IPA_pri_satsuma$V5 > 22.05e6 & IPA_pri_satsuma$V5 < 22.15e6
table(IPA_pri_satsuma[reg_filter, "V1"])
plot(y = IPA_pri_satsuma[IPA_pri_satsuma$V4 == "chr13" & IPA_pri_satsuma$V5 > 22.05e6 & IPA_pri_satsuma$V5 < 22.15e6, "V2"], x = IPA_pri_satsuma[IPA_pri_satsuma$V4 == "chr13" & IPA_pri_satsuma$V5 > 22.05e6 & IPA_pri_satsuma$V5 < 22.15e6, "V5"])

reg_filter_alt <- IPA_alt_satsuma$V4 == "chr13" & IPA_alt_satsuma$V5 > 22.05e6 & IPA_alt_satsuma$V5 < 22.15e6
table(IPA_alt_satsuma[reg_filter_alt, "V1"])
plot(y = IPA_alt_satsuma[reg_filter_alt, "V2"], x = IPA_alt_satsuma[reg_filter_alt, "V5"])

chr13_22.1_Mb <- Sprat_v_herring_aln(herring_GR = GRanges("chr13:22.10e6-22.14e6"), sprat_pri_GR = GRanges("ctg.000944F:0.4e4-4.5e4"), sprat_pri_dir = "neg", sprat_alt_GR = GRanges("hap_ctg.003798F JUNK:2.0e4-5.1e4"), sprat_alt_dir = "pos")
chr13_22.1_Mb_clustalo <- clustalomega(x = chr13_22.1_Mb)
save(chr13_22.1_Mb_clustalo, file = "~/Projects/Sprat/data/tree/chr13_22.1_Mb_clustalo.RData")
aln_basic_plots(clustal_obj = chr13_22.1_Mb_clustalo, pdf_file = "~/Projects/Sprat/doc/chr13_22.10e6_22.14e6_trees.pdf")


reg_filter <- IPA_pri_satsuma$V4 == "chr15" & IPA_pri_satsuma$V5 > 9.8e6 & IPA_pri_satsuma$V5 < 10.0e6
table(IPA_pri_satsuma[reg_filter, "V1"])
plot(y = IPA_pri_satsuma[reg_filter, "V2"], x = IPA_pri_satsuma[reg_filter, "V5"])

reg_filter_alt <- IPA_alt_satsuma$V4 == "chr15" & IPA_alt_satsuma$V5 > 9.8e6 & IPA_alt_satsuma$V5 < 10.0e6
table(IPA_alt_satsuma[reg_filter_alt, "V1"])
plot(y = IPA_alt_satsuma[reg_filter_alt, "V2"], x = IPA_alt_satsuma[reg_filter_alt, "V5"])


chr15_9.9_Mb <- Sprat_v_herring_aln(herring_GR = GRanges("chr15:9.85e6-9.90e6"), sprat_pri_GR = GRanges("ctg.000078F:0.9e5-1.4e5"), sprat_pri_dir = "pos", sprat_alt_GR = GRanges("hap_ctg.000731F HAPLOTIG:0.90e5-1.4e5"), sprat_alt_dir = "pos")
chr15_9.9_Mb_clustalo <- clustalomega(x = chr15_9.9_Mb)
save(chr15_9.9_Mb_clustalo, file = "~/Projects/Sprat/data/tree/chr15_9.9_Mb_clustalo.RData")
aln_basic_plots(clustal_obj = chr15_9.9_Mb_clustalo, pdf_file = "~/Projects/Sprat/doc/chr15_9.85e6_9.90e6_trees.pdf")



reg_filter <- IPA_pri_satsuma$V4 == "chr15" & IPA_pri_satsuma$V5 > 8.8e6 & IPA_pri_satsuma$V5 < 9.0e6
table(IPA_pri_satsuma[reg_filter, "V1"])
plot(y = IPA_pri_satsuma[reg_filter, "V2"], x = IPA_pri_satsuma[reg_filter, "V5"])

reg_filter_alt <- IPA_alt_satsuma$V4 == "chr15" & IPA_alt_satsuma$V5 > 8.8e6 & IPA_alt_satsuma$V5 < 9.0e6
table(IPA_alt_satsuma[reg_filter_alt, "V1"])
plot(y = IPA_alt_satsuma[reg_filter_alt, "V2"], x = IPA_alt_satsuma[reg_filter_alt, "V5"])

chr15_8.88_Mb <- Sprat_v_herring_aln(herring_GR = GRanges("chr15:8.86e6-8.92e6"), sprat_pri_GR = GRanges("ctg.002308F:2.5e4-7.5e4"), sprat_pri_dir = "neg", sprat_alt_GR = GRanges("hap_ctg.002584F HAPLOTIG:4.0e4-8.7e4"), sprat_alt_dir = "neg")
chr15_8.88_Mb_clustalo <- clustalomega(x = chr15_8.88_Mb)
save(chr15_8.88_Mb_clustalo, file = "~/Projects/Sprat/data/tree/chr15_8.88_Mb_clustalo.RData")
aln_basic_plots(clustal_obj = chr15_8.88_Mb_clustalo, pdf_file = "~/Projects/Sprat/doc/chr15_8.86e6_8.90e6_trees.pdf")

chr15_8.93_Mb <- Sprat_v_herring_aln(herring_GR = GRanges("chr15:8.92e6-8.945e6"), sprat_pri_GR = GRanges("ctg.002308F:1-2.5e4"), sprat_pri_dir = "neg", sprat_alt_GR = GRanges("hap_ctg.002584F HAPLOTIG:1-3.5e4"), sprat_alt_dir = "neg")
chr15_8.93_Mb_clustalo <- clustalomega(x = chr15_8.93_Mb)
save(chr15_8.93_Mb_clustalo, file = "~/Projects/Sprat/data/tree/chr15_8.93_Mb_clustalo.RData")
aln_basic_plots(clustal_obj = chr15_8.93_Mb_clustalo, pdf_file = "~/Projects/Sprat/doc/chr15_8.92e6_8.94e6_trees.pdf")


reg_filter <- IPA_pri_satsuma$V4 == "chr8" & IPA_pri_satsuma$V5 > 20.8e6 & IPA_pri_satsuma$V5 < 21.6e6
table(IPA_pri_satsuma[reg_filter, "V1"])
plot(y = IPA_pri_satsuma[reg_filter, "V2"], x = IPA_pri_satsuma[reg_filter, "V5"])

reg_filter_alt <- IPA_alt_satsuma$V4 == "chr8" & IPA_alt_satsuma$V5 > 20.8e6 & IPA_alt_satsuma$V5 < 21.6e6
table(IPA_alt_satsuma[reg_filter_alt, "V1"])
plot(y = IPA_alt_satsuma[reg_filter_alt, "V2"], x = IPA_alt_satsuma[reg_filter_alt, "V5"])

chr8_SDR <- Sprat_v_herring_aln(herring_GR = GRanges("chr8:21.15e6-21.23e6"), sprat_pri_GR = GRanges("ctg.000592F:2.70e5-3.5e5"), sprat_pri_dir = "neg", sprat_alt_GR = GRanges("hap_ctg.001096F HAPLOTIG:2.70e5-3.50e5"), sprat_alt_dir = "pos")
chr8_SDR_clustalo <- clustalomega(x = chr8_SDR)
save(chr8_SDR_clustalo, file = "~/Projects/Sprat/data/tree/chr8_SDR_clustalo.RData")
aln_basic_plots(clustal_obj = chr8_SDR_clustalo, pdf_file = "~/Projects/Sprat/doc/chr8_SDR_trees.pdf")


reg_filter <- IPA_pri_satsuma$V4 == "chr12" & IPA_pri_satsuma$V5 > 18.1e6 & IPA_pri_satsuma$V5 < 18.2e6
table(IPA_pri_satsuma[reg_filter, "V1"])
plot(y = IPA_pri_satsuma[reg_filter, "V2"], x = IPA_pri_satsuma[reg_filter, "V5"])

reg_filter_alt <- IPA_alt_satsuma$V4 == "chr12" & IPA_alt_satsuma$V5 > 18.1e6 & IPA_alt_satsuma$V5 < 18.2e6
table(IPA_alt_satsuma[reg_filter_alt, "V1"])
plot(y = IPA_alt_satsuma[reg_filter_alt, "V2"], x = IPA_alt_satsuma[reg_filter_alt, "V5"])

chr12_18.15_Mb <- Sprat_v_herring_aln(herring_GR = GRanges("chr12:18.13e6-18.17e6"), sprat_pri_GR = GRanges("ctg.000421F:6.45e5-6.95e5"), sprat_pri_dir = "neg", sprat_alt_GR = GRanges("hap_ctg.001344F HAPLOTIG:1.00e5-1.50e5"), sprat_alt_dir = "pos")
chr12_18.15_Mb_clustalo <- clustalomega(x = chr12_18.15_Mb)
save(chr12_18.15_Mb_clustalo, file = "~/Projects/Sprat/data/tree/chr12_18.15_Mb_clustalo.RData")
aln_basic_plots(clustal_obj = chr12_18.15_Mb_clustalo, pdf_file = "~/Projects/Sprat/doc/chr12_18.15_Mb_trees.pdf")


reg_filter <- IPA_pri_satsuma$V4 == "chr6" & IPA_pri_satsuma$V5 > 10.1e6 & IPA_pri_satsuma$V5 < 10.2e6
chr6_10.15_Mb <- Sprat_v_herring_aln(herring_GR = GRanges("chr6:10.11e6-10.15e6"), sprat_pri_GR = GRanges("ctg.001110F:1.45e5-1.85e5"), sprat_pri_dir = "pos", sprat_alt_GR = sprt_alt_GR, sprat_alt_dir = "pos")
chr6_10.15_Mb_clustalo <- clustalomega(x = chr6_10.15_Mb)
table(IPA_pri_satsuma[reg_filter, "V1"])
plot(y = IPA_pri_satsuma[reg_filter, "V2"], x = IPA_pri_satsuma[reg_filter, "V5"])

reg_filter_alt <- IPA_alt_satsuma$V4 == "chr6" & IPA_alt_satsuma$V5 > 10.1e6 & IPA_alt_satsuma$V5 < 10.2e6
table(IPA_alt_satsuma[reg_filter_alt, "V1"])
plot(y = IPA_alt_satsuma[reg_filter_alt, "V2"], x = IPA_alt_satsuma[reg_filter_alt, "V5"])
sprt_alt_GR <- GRanges(seqnames = "hap_ctg.001144F REPEAT", ranges = IRanges(start=1.30e5, end = 1.70e5))
save(chr6_10.15_Mb_clustalo, file = "~/Projects/Sprat/data/tree/chr6_10.15_Mb_clustalo.RData")
aln_basic_plots(clustal_obj = chr6_10.15_Mb_clustalo, pdf_file = "~/Projects/Sprat/doc/chr6_10.15_Mb_trees.pdf")


reg_filter <- IPA_pri_satsuma$V4 == "chr25" & IPA_pri_satsuma$V5 > 8.1e6 & IPA_pri_satsuma$V5 < 8.2e6
table(IPA_pri_satsuma[reg_filter, "V1"])
plot(y = IPA_pri_satsuma[reg_filter, "V2"], x = IPA_pri_satsuma[reg_filter, "V5"])

reg_filter_alt <- IPA_alt_satsuma$V4 == "chr25" & IPA_alt_satsuma$V5 > 8.1e6 & IPA_alt_satsuma$V5 < 8.2e6
table(IPA_alt_satsuma[reg_filter_alt, "V1"])
plot(y = IPA_alt_satsuma[reg_filter_alt, "V2"], x = IPA_alt_satsuma[reg_filter_alt, "V5"])

chr25_8.15_Mb <- Sprat_v_herring_aln(herring_GR = GRanges("chr25:8.10e6-8.14e6"), sprat_pri_GR = GRanges("ctg.000437F:0.25e5-0.65e5"), sprat_pri_dir = "neg", sprat_alt_GR = GRanges("hap_ctg.000522F HAPLOTIG:6.65e5-7.10e5"), sprat_alt_dir = "pos")
chr25_8.15_Mb_clustalo <- clustalomega(x = chr25_8.15_Mb)
save(chr25_8.15_Mb_clustalo, file = "~/Projects/Sprat/data/tree/chr25_8.15_Mb_clustalo.RData")
aln_basic_plots(clustal_obj = chr25_8.15_Mb_clustalo, pdf_file = "~/Projects/Sprat/doc/chr25_8.15_Mb_trees.pdf")

reg_filter <- IPA_pri_satsuma$V4 == "chr4" & IPA_pri_satsuma$V5 > 30.8e6 & IPA_pri_satsuma$V5 < 30.815e6
table(IPA_pri_satsuma[reg_filter, "V1"])
plot(y = IPA_pri_satsuma[reg_filter, "V2"], x = IPA_pri_satsuma[reg_filter, "V5"])

reg_filter_alt <- IPA_alt_satsuma$V4 == "chr4" & IPA_alt_satsuma$V5 > 30.8e6 & IPA_alt_satsuma$V5 < 30.815e6
table(IPA_alt_satsuma[reg_filter_alt, "V1"])
plot(y = IPA_alt_satsuma[reg_filter_alt, "V2"], x = IPA_alt_satsuma[reg_filter_alt, "V5"])

chr4_3.85_Mb <- Sprat_v_herring_aln(herring_GR = GRanges("chr4:30.80e6-30.815e6"), sprat_pri_GR = GRanges("ctg.001953F:1.45e5-1.60e5"), sprat_pri_dir = "neg", sprat_alt_GR = GRanges("hap_ctg.002132F OVLP:0.14e5-0.29e5"), sprat_alt_dir = "pos")
chr4_3.85_Mb_clustalo <- clustalomega(x = chr4_3.85_Mb)
save(chr4_3.85_Mb_clustalo, file = "~/Projects/Sprat/data/tree/chr4_3.85_Mb_clustalo.RData")
aln_basic_plots(clustal_obj = chr4_3.85_Mb_clustalo, pdf_file = "~/Projects/Sprat/doc/chr4_3.85_Mb_trees.pdf")


#Region on Ch_v2.0.2 chr17-14.7 Mb (Introgressed into Arctic herring)
IPA_pri_chr17_14.7 <- IPA_pri_satsuma[IPA_pri_satsuma[,4] == "chr17" & IPA_pri_satsuma[,6] > 14.68e6 & IPA_pri_satsuma[,5] < 14.70e6,]
plot(y = IPA_pri_chr17_14.7[, "V2"], x = IPA_pri_chr17_14.7[, "V5"])
IPA_pri_chr17_14.7 <- IPA_pri_satsuma[IPA_pri_satsuma[,4] == "chr17" & IPA_pri_satsuma[,6] > 14.738e6 & IPA_pri_satsuma[,5] < 14.743e6,]
plot(y = IPA_pri_chr17_14.7[, "V2"], x = IPA_pri_chr17_14.7[, "V5"])

IPA_alt_chr17_14.7 <- IPA_alt_satsuma[IPA_alt_satsuma[,4] == "chr17" & IPA_alt_satsuma[,6] > 14.68e6 & IPA_alt_satsuma[,5] < 14.70e6,]
plot(y = IPA_alt_chr17_14.7[, "V2"], x = IPA_alt_chr17_14.7[, "V5"])


chr17_14.69_Mb <- Sprat_v_herring_aln(herring_GR = GRanges("chr17:14.685e6-14.700e6"), sprat_pri_GR = GRanges("ctg.003110F:0.22e5-0.37e5"), sprat_pri_dir = "neg", sprat_alt_GR = GRanges("hap_ctg.003053F HAPLOTIG:0.39e5-0.55e5"), sprat_alt_dir = "pos")
chr17_14.69_clustalo <- clustalomega(x = chr17_14.69_Mb)
save(chr17_14.69_clustalo, file = "~/Projects/Sprat/data/tree/chr17_14.69_Mb_clustalo.RData")
aln_basic_plots(clustal_obj = chr17_14.69_clustalo, pdf_file = "~/Projects/Sprat/doc/chr17_14.69_Mb_trees.pdf")


#Whole-genome self-blast

#main HiC_scaffolds
Sprat_HiC_1MB_scaff <-names(HiC_3D_DNA[order(HiC_3D_DNA@ranges@width, decreasing = T)][1:35])

#Sprat_HiC_blast <- read.table("~/Projects/Sprat/data/3D_DNA_blast/sprat_3D_DNA_selfblast.out.gz", stringsAsFactors=F)
Sprat_HiC_blast <- read.table("~/Projects/Sprat/data/3D_DNA_blast/sprat_3D_DNA_selfblast_cl10_e10.out", stringsAsFactors=F)
#Sprat_HiC_blast <- Sprat_HiC_blast[-1519715,] #Using megablast and higher culling limit, tmp fix for incomplete file
#Seems better, will replace
save(Sprat_HiC_blast, file = "~/Projects/Sprat/data/3D_DNA_blast/sprat_3D_DNA_selfblast.RData")

off_diag_co <- 500
length_co <- 1000
Sprat_blast_summary_df <- data.frame(row.names = Sprat_HiC_1MB_scaff, on_diag = rep(NA, length(Sprat_HiC_1MB_scaff)), off_diag = NA, stringsAsFactors = F)

for(q_seq in Sprat_HiC_1MB_scaff){
  png_file <- paste0("~/Projects/Sprat/doc/blast/", q_seq,"_selfblast.png")
  png(png_file, width = 2000, height = 2000)
  blast_subset <- Sprat_HiC_blast[Sprat_HiC_blast[,1] == q_seq & Sprat_HiC_blast[,4] > length_co, ]
  #main_t_seq <- names(table(blast_subset[,2]))[order(table(blast_subset[,2]), decreasing=T)][1]
  #if(main_t_seq == q_seq){
   # main_t_seq <- names(table(blast_subset[,2]))[order(table(blast_subset[,2]), decreasing=T)][2]
  #}
  #main_hit_filter <- blast_subset[,2] == main_t_seq
  #aln_fit <- lm(blast_subset[main_hit_filter,7]~blast_subset[main_hit_filter,9])
  #aln_direction <- sign(aln_fit$coefficients[2])
 
  self_off_diag <- blast_subset[,2] == q_seq & (pmin(abs(blast_subset[,7] - blast_subset[,9]), abs(blast_subset[,7] - blast_subset[,10])) > off_diag_co | pmin(abs(blast_subset[,8] - blast_subset[,9]), abs(blast_subset[,8] - blast_subset[,10])) > off_diag_co)
  self_on_diag <- blast_subset[,2] == q_seq & !self_off_diag
  off_filter <- (blast_subset[,2] != q_seq) & (blast_subset[,2] %in% Sprat_HiC_1MB_scaff)
  Sprat_blast_summary_df[q_seq, "on_diag"] <- sum(blast_subset[self_on_diag,4])
  Sprat_blast_summary_df[q_seq, "off_diag"] <- sum(blast_subset[self_off_diag,4])
  
  plot(x = 1, y = 1, type = "n", xlim = c(0, max(blast_subset[,7:8])), ylim = c(0, max(blast_subset[,7:8])), xlab = q_seq, ylab = q_seq)
  segments(x0 = blast_subset[off_filter,7], x1 = blast_subset[off_filter,8], y0 = blast_subset[off_filter,9], y1 = blast_subset[off_filter,10], col = "grey30", lty = "dotted", lwd = 0.5)
  segments(x0 = blast_subset[self_on_diag,7], x1 = blast_subset[self_on_diag,8], y0 = blast_subset[self_on_diag,9], y1 = blast_subset[self_on_diag,10], col = "blue", lwd = 2)
  segments(x0 = blast_subset[self_off_diag,7], x1 = blast_subset[self_off_diag,8], y0 = blast_subset[self_off_diag,9], y1 = blast_subset[self_off_diag,10], col = "red", lwd = 2)	
  dev.off()
}
Sprat_blast_summary_df[, "off_frac"] <- (Sprat_blast_summary_df[,2]/2)/HiC_3D_DNA[Sprat_HiC_1MB_scaff]@ranges@width
Sprat_blast_summary_df[, "tot_size"] <- HiC_3D_DNA[Sprat_HiC_1MB_scaff]@ranges@width


#Pin HiC re-assembly attempt
#Note the renaming of the final fasta
#rsync -axv --progress --no-group matsp@rackham.uppmax.uu.se:/proj/snic2020-2-19/private/sprat/users/mats/pin_HiC/pin_out/scaffolds_final.fa ./Sprat_pin_HiC_r1.fasta
require(Biostrings)
Sprat_pin_HiC <- readDNAStringSet("~/Projects/Sprat/data/assemblies/pin_HiC/Sprat_pin_HiC_r1.fasta.gz")
Sprat_pin_v_Ch_v2.0.2 <- read_and_plot_fish_satsuma2("~/Projects/Sprat/data/HiC_satsuma/Sprat_pin_HiC_v_Ch_v2.0.2_satsuma.out_sorted", fish_name = "Sprat_pin_v_Ch_v2.0.2", plot_dir = "~/Projects/Sprat/doc/")
Sprat_pin_v_Ch_v2.0.2_syn <- satsuma_synteny_blocks2("~/Projects/Sprat/data/HiC_satsuma/Sprat_pin_HiC_v_Ch_v2.0.2_satsuma.out_sorted", "~/Projects/Sprat/doc/Sprat_pin_v_Ch_v2.0.2_synteny.pdf", hit_co = 20000)
save(Sprat_pin_v_Ch_v2.0.2, Sprat_pin_v_Ch_v2.0.2_syn, file = "~/Projects/Sprat/data/HiC_satsuma/Sprat_pin_HiC_v_Ch_v2.0.2.RData")

Sprat_pin_HiC <- readDNAStringSet("~/Projects/Sprat/data/HiC_satsuma/pin_HiC_it6/Sprat_pin_HiC_it6.fasta.gz")
Sprat_pin_it6_v_Ch_v2.0.2_syn <- satsuma_synteny_blocks2("~/Projects/Sprat/data/HiC_satsuma/pin_HiC_it6/Sprat_pin_HiC_it6_v_Ch_v2.0.2_satsuma.out_sorted", "~/Projects/Sprat/doc/Sprat_pin_it6_v_Ch_v2.0.2_synteny.pdf", hit_co = 20000)
save(Sprat_pin_v_Ch_v2.0.2, Sprat_pin_v_Ch_v2.0.2_syn, Sprat_pin_it6_v_Ch_v2.0.2_syn, file = "~/Projects/Sprat/data/HiC_satsuma/Sprat_pin_HiC_v_Ch_v2.0.2.RData")



#Getting to juicebox
#Pruning non-paired reads
#module load bioinfo-tools samtools python3
#~/.local/bin/pairtools parse --assembly Sprat_IPA_pri -c /proj/snic2020-2-19/private/sprat/users/mats/pin_HiC/pruned_bams/Sprat_IPA_primary.sizes /proj/snic2020-2-19/private/sprat/users/mats/pin_HiC/pruned_bams/matlock_test.bam | gzip > /proj/snic2020-2-19/private/sprat/users/mats/pin_HiC/pruned_bams/prune_test_sbatch.pairsam.gz
#~/.local/bin/pairtools sort  --cmd-in gunzip ./prune_test_sbatch.pairsam.gz | gzip > prune_test_sorted.pairsam.gz
#~/.local/bin/pairtools dedup --mark-dups --cmd-in gunzip ./prune_test_sorted.pairsam.gz | gzip > prune_test_dedup.pairsam.gz
#~/.local/bin/pairtools select '(pair_type == "UU")' --cmd-in gunzip ./prune_test_dedup.pairsam.gz | gzip > prune_test_filter.pairsam.gz
#~/.local/bin/pairtools split --output-sam prune_test.filtered.bam --cmd-in gunzip ./prune_test_filter.pairsam.gz 

#/home/matsp/private/Software/PhaseGenomics/matlock/bin/matlock bam2 juicer ./prune_test_workflow_filtered.bam ./Sprat_test.links

#bash /proj/snic2020-2-19/private/herring/users/mats/three_d_dna/visualize/run-assembly-visualizer.sh -p false ../pin_out_it6/Sprat_pin_it6.assembly Sprat_test.links 

Sprat_pin_curated <- readDNAStringSet("~/Projects/Sprat/data/HiC_satsuma/pin_HiC_it6/curation/Sprat_pin_milestone1.fasta.gz")
Sprat_pin_curated <- Sprat_pin_curated[order(Sprat_pin_curated@ranges@width, decreasing = T)]
Sprat_pin_m1_v_Ch_v2.0.2_syn <- satsuma_synteny_blocks2("~/Projects/Sprat/data/HiC_satsuma/", "~/Projects/Sprat/doc/Sprat_pin_m1_v_Ch_v2.0.2_synteny.pdf", hit_co = 20000)
save(Sprat_pin_m1_v_Ch_v2.0.2_syn, file = "~/Projects/Sprat/data/HiC_satsuma/Sprat_pin_m1_v_Ch_v2.0.2.RData")

#Preparing ranges from AGP for liftover
Sprat_m1_AGP <- read.table("~/Projects/Sprat/data/HiC_satsuma/pin_HiC_it6/curation/Sprat_pin_milestone1.agp", stringsAsFactors = F, sep = "\t", header = F)
names(Sprat_m1_AGP) <- c("chr_name", "chr_start", "chr_end", "idx", "type", "raw_contig", "fragment_start", "fragment_end", "direction")

Sprat_m1_AGP[Sprat_m1_AGP$type == "U", "direction"] <- "+"
Sprat_m1_AGP[Sprat_m1_AGP$type == "U", "raw_contig"] <- "HiC_gap"

Sprat_m1_AGP[,"size"] <- Sprat_m1_AGP$chr_end - Sprat_m1_AGP$chr_start + 1
Sprat_m1_AGP[,"contig"] <- sub("([0-9A-Za-z_|.]+)[:]{3}.+", "\\1", Sprat_m1_AGP[,"raw_contig"])

for(contig in unique(Sprat_m1_AGP[,"contig"])){
  contig_sizes <- Sprat_m1_AGP[Sprat_m1_AGP[,"contig"] == contig,"size"]
  if(length(contig_sizes) > 1){
    frag_order <- as.integer(sub("([0-9A-Za-z_|.]+)[:]{3}fragment_([0-9+]).*", "\\2", Sprat_m1_AGP[Sprat_m1_AGP[,"contig"] == contig,"raw_contig"]))
  } else {
    frag_order <- 1
  }
  #Incorrect!
  #Sprat_m1_AGP[Sprat_m1_AGP[,"contig"] == contig,"fragment_start"] <-(c(0, cumsum(contig_sizes[frag_order])[-length(contig_sizes)]) + 1)[frag_order]
  #Sprat_m1_AGP[Sprat_m1_AGP[,"contig"] == contig,"fragment_end"] <- cumsum(contig_sizes[frag_order])[frag_order]
  Sprat_m1_AGP[Sprat_m1_AGP[,"contig"] == contig,"fragment_start"] <- (c(0, cumsum(contig_sizes[match(1:length(contig_sizes), frag_order)])[-length(contig_sizes)]) + 1)[frag_order]
  Sprat_m1_AGP[Sprat_m1_AGP[,"contig"] == contig,"fragment_end"] <- cumsum(contig_sizes[match(1:length(contig_sizes), frag_order)])[frag_order]
}

Sprat_m1_AGP[Sprat_m1_AGP$type == "U", "fragment_start"] <- 1
Sprat_m1_AGP[Sprat_m1_AGP$type == "U", "fragment_end"] <- 100
Sprat_m1_AGP[, "fragment_start"] <- as.integer(Sprat_m1_AGP[, "fragment_start"])
Sprat_m1_AGP[, "fragment_end"] <- as.integer(Sprat_m1_AGP[, "fragment_end"])


Sprat_m1_AGP_GR <- GRanges(seqnames = Sprat_m1_AGP$contig, ranges = IRanges(start = Sprat_m1_AGP$fragment_start, end = Sprat_m1_AGP$fragment_end))
Sprat_m1_AGP_GR$direction <- Sprat_m1_AGP$direction
Sprat_m1_AGP_GR$HiC_scaffold <- Sprat_m1_AGP$chr_name
Sprat_m1_AGP_GR$HiC_start <- Sprat_m1_AGP$chr_start
Sprat_m1_AGP_GR$HiC_end <- Sprat_m1_AGP$chr_end

#Positions & frequencies to transfer
#load("~/Projects/Sprat/data/genotypes/Sprat_pool_freq_ext.RData")
#sprat_lo_test <- Sprat_pool_freq[sample.int(dim(Sprat_pool_freq)[1], 200),]
#sprat_lo_test_GR <- GRanges(seqnames = sprat_lo_test$CHROM, ranges = IRanges(start = sprat_lo_test$POS, end = sprat_lo_test$POS))
#for(i in 1:dim(sprat_lo_test)[2]){
#  sprat_lo_test_GR@elementMetadata[,i] <- sprat_lo_test[,i]
#  names(sprat_lo_test_GR@elementMetadata)[i] <- names(sprat_lo_test)[i]
#}
#sprat_lo_test_results <- agp_SNP_liftover(lo_GRs = sprat_lo_test_GR, AGP_PB_GR = Sprat_m1_AGP_GR)
#sprat_lo_test_results$lo_GRs@elementMetadata[1:5,c("CHROM", "POS", "HiC_scaffold", "HiC_start", "HiC_end" ,"direction", "PB_start", "HiC_pos")]


sprat_lo_GR <- GRanges(seqnames = Sprat_pool_freq$CHROM, ranges = IRanges(start = Sprat_pool_freq$POS, end = Sprat_pool_freq$POS))
for(i in 1:dim(Sprat_pool_freq)[2]){
  sprat_lo_GR@elementMetadata[,i] <- Sprat_pool_freq[,i]
  names(sprat_lo_GR@elementMetadata)[i] <- names(Sprat_pool_freq)[i]
}


save(sprat_lo_GR, file = "~/Projects/Sprat/data/genotypes/Sprat_pool_freq_ext_GR.RData")

sprat_post_lo_GR <- agp_SNP_liftover(lo_GRs = sprat_lo_GR, AGP_PB_GR = Sprat_m1_AGP_GR)

#HiC positions in this one is affected by the fragmenting ordering bug!
save(sprat_post_lo_GR, file = "~/Projects/Sprat/data/genotypes/Sprat_pool_freq_post_lo_GR.RData")
#

save(sprat_post_lo_GR, file = "~/Projects/Sprat/data/genotypes/Sprat_pool_freq_post_lo_GR_v2.RData")

#Site filtering
#Attmepting to eliminate SNPs where counts are both exactly 0 and 1
sprat_m1_freq_mat <- as.matrix(sprat_post_lo_GR$lo_GRs@elementMetadata[,grep("_Afreq", names(sprat_post_lo_GR$lo_GRs@elementMetadata))])
sprat_opp_fix_filter_m1 <- (rowSums(sprat_m1_freq_mat == 0, na.rm = T) > 0) & (rowSums(sprat_m1_freq_mat == 1, na.rm = T) > 0)
head(sprat_m1_freq_mat[sprat_opp_fix_filter_m1,])
rm(sprat_m1_freq_mat)                                                                           

#head(Sprat_pool_freq_lo[sprat_opp_fix_filter_lo,c(1,2,grep("_Afreq", names(Sprat_pool_freq_lo)))])
sprat_m1_DP_mat <- as.matrix(sprat_post_lo_GR$lo_GRs@elementMetadata[,grep("_DP", names(sprat_post_lo_GR$lo_GRs@elementMetadata))])
sprat_m1_dp_filter <- !(rowSums(sprat_m1_DP_mat  <= 20, na.rm = T) > 0)
head(sprat_m1_DP_mat[sprat_m1_dp_filter,])
rm(sprat_m1_DP_mat) 
sprat_site_filter_m1 <- !sprat_opp_fix_filter_m1 & sprat_m1_dp_filter
save(sprat_site_filter_m1, file = "~/Projects/Sprat/data/genotypes/sprat_site_filter_m1.RData")


#Looking into putative inversions
#Chr4
#Using the Herring liftover to find critical contigs
chr4_inv_filter <- sprat_pool_dp_filter_lo & Sprat_pool_freq_lo$herring_v2.0.2_seqnames == "chr4" & Sprat_pool_freq_lo$SNP_HiC_pos > 25e6
chr4_inv_ctg <- names(table(Sprat_pool_freq_lo$CHROM[chr4_inv_filter]))[table(Sprat_pool_freq_lo$CHROM[chr4_inv_filter]) > 50]
m1_chr4_inv_filter <- as.character(sprat_post_lo_GR$lo_GRs@seqnames) %in% chr4_inv_ctg
table(sprat_post_lo_GR$lo_GRs$HiC_scaffold[m1_chr4_inv_filter])
m1_chr4_inv_main_scaff <- "PGA_scaffold_1110__33_contigs__length_20919432" 
plot(x= sprat_post_lo_GR$lo_GRs$HiC_pos[sprat_post_lo_GR$lo_GRs$HiC_scaffold == m1_chr4_inv_main_scaff & sprat_site_filter_m1], y = sprat_post_lo_GR$lo_GRs$oce_v_brack[sprat_post_lo_GR$lo_GRs$HiC_scaffold == m1_chr4_inv_main_scaff & sprat_site_filter_m1])


#Chr5
#Using the Herring liftover to find critical contigs
chr5_inv_filter <- sprat_pool_dp_filter_lo & Sprat_pool_freq_lo$herring_v2.0.2_seqnames == "chr5" & Sprat_pool_freq_lo$SNP_HiC_pos > 10e6 & Sprat_pool_freq_lo$SNP_HiC_pos < 15e6
chr5_inv_ctg <- names(table(Sprat_pool_freq_lo$CHROM[chr5_inv_filter]))[table(Sprat_pool_freq_lo$CHROM[chr5_inv_filter]) > 50]
m1_chr5_inv_filter <- as.character(sprat_post_lo_GR$lo_GRs@seqnames) %in% chr5_inv_ctg
table(sprat_post_lo_GR$lo_GRs$HiC_scaffold[m1_chr5_inv_filter])
m1_chr5_inv_main_scaff <- "PGA_scaffold_1114__20_contigs__length_23438440" 
plot(x= sprat_post_lo_GR$lo_GRs$HiC_pos[sprat_post_lo_GR$lo_GRs$HiC_scaffold == m1_chr5_inv_main_scaff & sprat_site_filter_m1], y = sprat_post_lo_GR$lo_GRs$oce_v_brack[sprat_post_lo_GR$lo_GRs$HiC_scaffold == m1_chr5_inv_main_scaff & sprat_site_filter_m1])


#Chr7
#Using the Herring liftover to find critical contigs
chr7_inv_filter <- sprat_pool_dp_filter_lo & Sprat_pool_freq_lo$herring_v2.0.2_seqnames == "chr7" & Sprat_pool_freq_lo$SNP_HiC_pos > 23e6
chr7_inv_ctg <- names(table(Sprat_pool_freq_lo$CHROM[chr7_inv_filter]))[table(Sprat_pool_freq_lo$CHROM[chr7_inv_filter]) > 50]
m1_chr7_inv_filter <- as.character(sprat_post_lo_GR$lo_GRs@seqnames) %in% chr7_inv_ctg
table(sprat_post_lo_GR$lo_GRs$HiC_scaffold[m1_chr7_inv_filter])
m1_chr7_inv_main_scaff <- "PGA_scaffold_13__40_contigs__length_28874545" 
plot(x= sprat_post_lo_GR$lo_GRs$HiC_pos[sprat_post_lo_GR$lo_GRs$HiC_scaffold == m1_chr7_inv_main_scaff & sprat_site_filter_m1], y = sprat_post_lo_GR$lo_GRs$oce_v_brack[sprat_post_lo_GR$lo_GRs$HiC_scaffold == m1_chr7_inv_main_scaff & sprat_site_filter_m1])
m1_chr7_lo_col_filter <- sprat_post_lo_GR$lo_GRs$HiC_scaffold == m1_chr13_inv_main_scaff & sprat_site_filter_m1 & m1_chr7_inv_filter 
points(x= sprat_post_lo_GR$lo_GRs$HiC_pos[m1_chr7_lo_col_filter], y = sprat_post_lo_GR$lo_GRs$oce_v_brack[m1_chr7_lo_col_filter], pch = 20, col = "blue")

tmp_filter <- sprat_post_lo_GR$lo_GRs$HiC_scaffold == m1_chr7_inv_main_scaff & sprat_site_filter_m1 & sprat_post_lo_GR$lo_GRs$oce_v_brack > 0.5 & sprat_post_lo_GR$lo_GRs$HiC_pos > 23e6
simple_sprat_pool_hm_v2(snp_data = sprat_post_lo_GR$lo_GRs, reg_filter = tmp_filter, pdf_file =  "~/Projects/Sprat/doc/reg_heatmaps_m1/PGA_scaffold_13_25MB.pdf")



#Chr13
#Using the Herring liftover to find critical contigs
chr13_inv_filter <- sprat_pool_dp_filter_lo & Sprat_pool_freq_lo$herring_v2.0.2_seqnames == "chr13" & Sprat_pool_freq_lo$SNP_HiC_pos > 7e6 & Sprat_pool_freq_lo$SNP_HiC_pos < 9e6
chr13_inv_ctg <- names(table(Sprat_pool_freq_lo$CHROM[chr13_inv_filter]))[table(Sprat_pool_freq_lo$CHROM[chr13_inv_filter]) > 50]
m1_chr13_inv_filter <- as.character(sprat_post_lo_GR$lo_GRs@seqnames) %in% chr13_inv_ctg
table(sprat_post_lo_GR$lo_GRs$HiC_scaffold[m1_chr13_inv_filter])
m1_chr13_inv_main_scaff <- "PGA_scaffold_173__1_contigs__length_182285" 
plot(x= sprat_post_lo_GR$lo_GRs$HiC_pos[sprat_post_lo_GR$lo_GRs$HiC_scaffold == m1_chr13_inv_main_scaff & sprat_site_filter_m1], y = sprat_post_lo_GR$lo_GRs$oce_v_brack[sprat_post_lo_GR$lo_GRs$HiC_scaffold == m1_chr13_inv_main_scaff & sprat_site_filter_m1])

m1_chr13_lo_col_filter <- sprat_post_lo_GR$lo_GRs$HiC_scaffold == m1_chr13_inv_main_scaff & sprat_site_filter_m1 & m1_chr13_inv_filter 
points(x= sprat_post_lo_GR$lo_GRs$HiC_pos[m1_chr13_lo_col_filter], y = sprat_post_lo_GR$lo_GRs$oce_v_brack[m1_chr13_lo_col_filter], pch = 20, col = "red")


tmp_filter <- sprat_post_lo_GR$lo_GRs$HiC_scaffold == m1_chr13_inv_main_scaff & sprat_site_filter_m1 & sprat_post_lo_GR$lo_GRs$oce_v_brack > 0.5
simple_sprat_pool_hm_v2(snp_data = sprat_post_lo_GR$lo_GRs, reg_filter = tmp_filter, pdf_file =  "~/Projects/Sprat/doc/reg_heatmaps_m1/PGA_scaffold_173.pdf")



#Chr18
#Based on PB v HiC_m1 satsuma
plot(x= sprat_post_lo_GR$lo_GRs$HiC_pos[sprat_post_lo_GR$lo_GRs$HiC_scaffold == "PGA_scaffold_40__17_contigs__length_6491151"& sprat_site_filter_m1], y = sprat_post_lo_GR$lo_GRs$oce_v_brack[sprat_post_lo_GR$lo_GRs$HiC_scaffold == "PGA_scaffold_40__17_contigs__length_6491151"& sprat_site_filter_m1])




#All major scaffolds
Sprat_lo_manhattan_v2(sprat_data = sprat_post_lo_GR$lo_GRs, scaff_list = names(Sprat_pin_curated)[1:50], daf_col = "oce_v_brack", plot_prefix = "~/Projects/Sprat/doc/Oce_v_Brack_m1/Oce_vBrack_m1_", site_filt_vec = sprat_site_filter_m1 )
Sprat_lo_manhattan_v2(sprat_data = sprat_post_lo_GR$lo_GRs, scaff_list = names(Sprat_pin_curated)[1:50], daf_col = "oce_v_fjord", plot_prefix = "~/Projects/Sprat/doc/Oce_v_Fjord_m1/Oce_v_Fjord_m1_", site_filt_vec = sprat_site_filter_m1 )


tmp_filter_1113 <- sprat_post_lo_GR$lo_GRs$HiC_scaffold == "PGA_scaffold_1113__35_contigs__length_27024755" & sprat_site_filter_m1 & sprat_post_lo_GR$lo_GRs$oce_v_fjord > 0.5 & sprat_post_lo_GR$lo_GRs$HiC_pos < 5e6
simple_sprat_pool_hm_v2(snp_data = sprat_post_lo_GR$lo_GRs, reg_filter = tmp_filter_1113, pdf_file =  "~/Projects/Sprat/doc/reg_heatmaps_m1/PGA_scaffold_1113_3MB.pdf")

tmp_filter_4 <- sprat_post_lo_GR$lo_GRs$HiC_scaffold == "PGA_scaffold_4__24_contigs__length_24387585" & sprat_site_filter_m1 & sprat_post_lo_GR$lo_GRs$oce_v_brack > 0.2 & sprat_post_lo_GR$lo_GRs$HiC_pos < 20.0e6 & sprat_post_lo_GR$lo_GRs$HiC_pos > 18.5e6
simple_sprat_pool_hm_v2(snp_data = sprat_post_lo_GR$lo_GRs, reg_filter = tmp_filter_4, pdf_file =  "~/Projects/Sprat/doc/reg_heatmaps_m1/PGA_scaffold_4_19MB.pdf")

tmp_filter_69 <- sprat_post_lo_GR$lo_GRs$HiC_scaffold == "PGA_scaffold_69__9_contigs__length_3302248" & sprat_site_filter_m1 & sprat_post_lo_GR$lo_GRs$oce_v_brack > 0.6 & sprat_post_lo_GR$lo_GRs$HiC_pos > 2.5e6
simple_sprat_pool_hm_v2(snp_data = sprat_post_lo_GR$lo_GRs, reg_filter = tmp_filter_69, pdf_file =  "~/Projects/Sprat/doc/reg_heatmaps_m1/PGA_scaffold_69_3MB.pdf")

tmp_filter_1109 <- sprat_post_lo_GR$lo_GRs$HiC_scaffold == "PGA_scaffold_1109__40_contigs__length_34818411" & sprat_site_filter_m1 & sprat_post_lo_GR$lo_GRs$oce_v_fjord > 0.5 & sprat_post_lo_GR$lo_GRs$HiC_pos > 30.0e6
simple_sprat_pool_hm_v2(snp_data = sprat_post_lo_GR$lo_GRs, reg_filter = tmp_filter_1109, pdf_file =  "~/Projects/Sprat/doc/reg_heatmaps_m1/PGA_scaffold_1109_32MB.pdf")

tmp_filter_1111 <- sprat_post_lo_GR$lo_GRs$HiC_scaffold == "PGA_scaffold_1111__78_contigs__length_49772038" & sprat_site_filter_m1 & sprat_post_lo_GR$lo_GRs$oce_v_brack > 0.75 & sprat_post_lo_GR$lo_GRs$HiC_pos < 7.0e6
simple_sprat_pool_hm_v2(snp_data = sprat_post_lo_GR$lo_GRs, reg_filter = tmp_filter_1111, pdf_file =  "~/Projects/Sprat/doc/reg_heatmaps_m1/PGA_scaffold_1111_5MB.pdf")

tmp_filter_1114_1 <- sprat_post_lo_GR$lo_GRs$HiC_scaffold == "PGA_scaffold_1114__20_contigs__length_23438440" & sprat_site_filter_m1 & sprat_post_lo_GR$lo_GRs$oce_v_brack > 0.35 & sprat_post_lo_GR$lo_GRs$HiC_pos > 5.0e6 & sprat_post_lo_GR$lo_GRs$HiC_pos < 7.0e6
simple_sprat_pool_hm_v2(snp_data = sprat_post_lo_GR$lo_GRs, reg_filter = tmp_filter_1114_1, pdf_file =  "~/Projects/Sprat/doc/reg_heatmaps_m1/PGA_scaffold_1114_6MB.pdf")

tmp_filter_1114_2 <- sprat_post_lo_GR$lo_GRs$HiC_scaffold == "PGA_scaffold_1114__20_contigs__length_23438440" & sprat_site_filter_m1 & sprat_post_lo_GR$lo_GRs$oce_v_brack > 0.5 & sprat_post_lo_GR$lo_GRs$HiC_pos > 20.0e6
simple_sprat_pool_hm_v2(snp_data = sprat_post_lo_GR$lo_GRs, reg_filter = tmp_filter_1114_2, pdf_file =  "~/Projects/Sprat/doc/reg_heatmaps_m1/PGA_scaffold_1114_23MB.pdf")

tmp_filter_1132 <- sprat_post_lo_GR$lo_GRs$HiC_scaffold == "PGA_scaffold_1132__20_contigs__length_21373292" & sprat_site_filter_m1 & sprat_post_lo_GR$lo_GRs$oce_v_brack > 0.6
simple_sprat_pool_hm_v2(snp_data = sprat_post_lo_GR$lo_GRs, reg_filter = tmp_filter_1132, pdf_file =  "~/Projects/Sprat/doc/reg_heatmaps_m1/PGA_scaffold_1132_9MB.pdf")

tmp_filter_13_13Mb <- sprat_post_lo_GR$lo_GRs$HiC_scaffold == "PGA_scaffold_13__40_contigs__length_28874545" & sprat_site_filter_m1 & sprat_post_lo_GR$lo_GRs$oce_v_brack > 0.3 & sprat_post_lo_GR$lo_GRs$HiC_pos > 11.5e6 & sprat_post_lo_GR$lo_GRs$HiC_pos < 14.0e6
simple_sprat_pool_hm_v2(snp_data = sprat_post_lo_GR$lo_GRs, reg_filter = tmp_filter_13_13Mb, pdf_file =  "~/Projects/Sprat/doc/reg_heatmaps_m1/PGA_scaffold_13_13MB.pdf")

tmp_filter_13_25Mb <- sprat_post_lo_GR$lo_GRs$HiC_scaffold == "PGA_scaffold_13__40_contigs__length_28874545" & sprat_site_filter_m1 & sprat_post_lo_GR$lo_GRs$oce_v_brack > 0.5 & sprat_post_lo_GR$lo_GRs$HiC_pos > 23e6


##
##Distances/trees of selected regions
f3 <- Vectorize(FUN = function(X,Y, data_mat){sum(abs(data_mat[,X] - data_mat[,Y]))}, vectorize.args = c("X", "Y"))
#tmp_filter_comb <- tmp_filter_1111 | tmp_filter_1113 | tmp_filter_40 | tmp_filter_69 # More restricitve version
tmp_filter_comb <- tmp_filter_1111 | tmp_filter_1113 | tmp_filter_1109 | tmp_filter_69 | tmp_filter_13_13Mb | tmp_filter_1132 | tmp_filter_1114_1 | tmp_filter_1114_2 | tmp_filter_13_25Mb # Complements the "outside" version below



OvB_Reg_mat  <- as.matrix(sprat_post_lo_GR$lo_GRs@elementMetadata[tmp_filter_comb, grep("_Afreq", names(sprat_post_lo_GR$lo_GRs@elementMetadata))])
OvB_Reg_dist <- outer(FUN ="f3", X=1:dim(OvB_Reg_mat)[2], Y=1:dim(OvB_Reg_mat)[2], data_mat = OvB_Reg_mat)
rownames(OvB_Reg_dist) <- grep("_Afreq", names(sprat_post_lo_GR$lo_GRs@elementMetadata), value = T)
colnames(OvB_Reg_dist) <- grep("_Afreq", names(sprat_post_lo_GR$lo_GRs@elementMetadata), value = T)
diag(OvB_Reg_dist) <- NA
OvB_Reg_dist <- OvB_Reg_dist/sum(tmp_filter_comb)
save(OvB_Reg_dist, file = "~/Projects/Sprat/data/tree/OvB_Reg_dist.RData")
plot(hclust(as.dist(OvB_Reg_dist)))
require(ape)
OvB_Reg_tree <- bionj(as.dist(OvB_Reg_dist))

pdf(file = "~/Projects/Sprat/doc/comb_reg_tree.pdf")
plot.phylo(OvB_Reg_tree, type = "unrooted", lab4ut = "axial", no.margin = T)
dev.off()




#Removing inversions (and maybe some other divergent regions)
reg_filter_1113 <- sprat_post_lo_GR$lo_GRs$HiC_scaffold == "PGA_scaffold_1113__35_contigs__length_27024755" & sprat_site_filter_m1  & sprat_post_lo_GR$lo_GRs$HiC_pos < 5e6
reg_filter_4 <- sprat_post_lo_GR$lo_GRs$HiC_scaffold == "PGA_scaffold_4__24_contigs__length_24387585" & sprat_site_filter_m1 & sprat_post_lo_GR$lo_GRs$HiC_pos < 20.0e6 & sprat_post_lo_GR$lo_GRs$HiC_pos > 18.5e6
reg_filter_69 <- sprat_post_lo_GR$lo_GRs$HiC_scaffold == "PGA_scaffold_69__9_contigs__length_3302248" & sprat_site_filter_m1 & sprat_post_lo_GR$lo_GRs$HiC_pos > 2.5e6
reg_filter_1109 <- sprat_post_lo_GR$lo_GRs$HiC_scaffold == "PGA_scaffold_1109__40_contigs__length_34818411" & sprat_site_filter_m1 & sprat_post_lo_GR$lo_GRs$HiC_pos > 30.0e6
reg_filter_1111 <- sprat_post_lo_GR$lo_GRs$HiC_scaffold == "PGA_scaffold_1111__78_contigs__length_49772038" & sprat_site_filter_m1 & sprat_post_lo_GR$lo_GRs$HiC_pos < 7.0e6
reg_filter_1114_1 <- sprat_post_lo_GR$lo_GRs$HiC_scaffold == "PGA_scaffold_1114__20_contigs__length_23438440" & sprat_site_filter_m1 & sprat_post_lo_GR$lo_GRs$HiC_pos > 5.0e6 & sprat_post_lo_GR$lo_GRs$HiC_pos < 7.0e6
reg_filter_1114_2 <- sprat_post_lo_GR$lo_GRs$HiC_scaffold == "PGA_scaffold_1114__20_contigs__length_23438440" & sprat_site_filter_m1 & sprat_post_lo_GR$lo_GRs$HiC_pos > 20.0e6
reg_filter_1132 <- sprat_post_lo_GR$lo_GRs$HiC_scaffold == "PGA_scaffold_1132__20_contigs__length_21373292" & sprat_site_filter_m1
reg_filter_13_13Mb <- sprat_post_lo_GR$lo_GRs$HiC_scaffold == "PGA_scaffold_13__40_contigs__length_28874545" & sprat_site_filter_m1 & sprat_post_lo_GR$lo_GRs$HiC_pos > 11.5e6 & sprat_post_lo_GR$lo_GRs$HiC_pos < 14.0e6
reg_filter_13_25Mb <- sprat_post_lo_GR$lo_GRs$HiC_scaffold == "PGA_scaffold_13__40_contigs__length_28874545" & sprat_site_filter_m1  & sprat_post_lo_GR$lo_GRs$HiC_pos > 23e6

##Is currently not removing all (or even any?) inversions!
no_inv_filter <- sprat_site_filter_m1 & !(reg_filter_13_13Mb | reg_filter_1111 | reg_filter_1113 | reg_filter_1109 | reg_filter_69 | reg_filter_4 | reg_filter_1132 | reg_filter_1114_1 | reg_filter_1114_2 | reg_filter_13_25Mb)
no_inv_mat  <- as.matrix(sprat_post_lo_GR$lo_GRs@elementMetadata[no_inv_filter, grep("_Afreq", names(sprat_post_lo_GR$lo_GRs@elementMetadata))])
no_inv_dist <- outer(FUN ="f3", X=1:dim(no_inv_mat)[2], Y=1:dim(no_inv_mat)[2], data_mat = no_inv_mat)
rownames(no_inv_dist) <- grep("_Afreq", names(sprat_post_lo_GR$lo_GRs@elementMetadata), value = T)
colnames(no_inv_dist) <- grep("_Afreq", names(sprat_post_lo_GR$lo_GRs@elementMetadata), value = T)
diag(no_inv_dist) <- NA
no_inv_dist <- no_inv_dist/sum(no_inv_filter)
save(no_inv_dist, file = "~/Projects/Sprat/data/tree/no_inv_dist.RData")
###


no_inv_tree <- bionj(as.dist(no_inv_dist))
pdf(file = "~/Projects/Sprat/doc/no_inv_tree.pdf")
plot.phylo(no_inv_tree, type = "unrooted", lab4ut = "axial", no.margin = T)
dev.off()


#Within-group comparisons
brackish_no_BlackSea_freq_lo_cols <- grep("(AB|GOTB|BBS|LAND15|LAND19|GD)_Afreq", names(sprat_post_lo_GR$lo_GRs@elementMetadata), value = T)
#head(rowMeans(as.matrix(sprat_post_lo_GR$lo_GRs@elementMetadata[,brackish_no_BlackSea_freq_lo_cols]), na.rm = T))


BS_LAND_freq_lo_cols <- grep("([.]BS|LAND15|LAND19)_Afreq", names(sprat_post_lo_GR$lo_GRs@elementMetadata), value = T)
Baltic_freq_lo_cols <- grep("(AB|GOTB|BBS|GD)_Afreq", names(sprat_post_lo_GR$lo_GRs@elementMetadata), value = T)
sprat_post_lo_GR$lo_GRs$Baltic_mean <- rowMeans(as.matrix(sprat_post_lo_GR$lo_GRs@elementMetadata[,Baltic_freq_lo_cols]), na.rm = T)

NS_freq_lo_cols <- grep("(NS[0-9]+)_Afreq", names(sprat_post_lo_GR$lo_GRs@elementMetadata), value = T)
S_Atl_freq_lo_cols <- grep("(BoB|CEL)_Afreq", names(sprat_post_lo_GR$lo_GRs@elementMetadata), value = T)
SK_freq_lo_cols <- grep("(SK[0-9]+)_Afreq", names(sprat_post_lo_GR$lo_GRs@elementMetadata), value = T)

sprat_post_lo_GR$lo_GRs$brackish_v_BlackSea <- abs(rowMeans(as.matrix(sprat_post_lo_GR$lo_GRs@elementMetadata[,brackish_no_BlackSea_freq_lo_cols]), na.rm = T) - sprat_post_lo_GR$lo_GRs@elementMetadata[,"P30.BS_Afreq"])
sprat_post_lo_GR$lo_GRs$BS_LAND_v_Baltic <- abs(rowMeans(as.matrix(sprat_post_lo_GR$lo_GRs@elementMetadata[,BS_LAND_freq_lo_cols]), na.rm = T) - rowMeans(as.matrix(sprat_post_lo_GR$lo_GRs@elementMetadata[,Baltic_freq_lo_cols]), na.rm = T))
sprat_post_lo_GR$lo_GRs$NS_v_S_Atl <- abs(rowMeans(as.matrix(sprat_post_lo_GR$lo_GRs@elementMetadata[,NS_freq_lo_cols]), na.rm = T) - rowMeans(as.matrix(sprat_post_lo_GR$lo_GRs@elementMetadata[,S_Atl_freq_lo_cols]), na.rm = T))
sprat_post_lo_GR$lo_GRs$SK_v_S_Atl <- abs(rowMeans(as.matrix(sprat_post_lo_GR$lo_GRs@elementMetadata[,SK_freq_lo_cols ]), na.rm = T) - rowMeans(as.matrix(sprat_post_lo_GR$lo_GRs@elementMetadata[,S_Atl_freq_lo_cols]), na.rm = T))
sprat_post_lo_GR$lo_GRs$N_v_S_Atl <- abs(rowMeans(as.matrix(sprat_post_lo_GR$lo_GRs@elementMetadata[,c(SK_freq_lo_cols,NS_freq_lo_cols)]), na.rm = T) - rowMeans(as.matrix(sprat_post_lo_GR$lo_GRs@elementMetadata[,S_Atl_freq_lo_cols]), na.rm = T))
sprat_post_lo_GR$lo_GRs$Oce_v_Baltic <- abs(sprat_post_lo_GR$lo_GRs$mean_oceanic - sprat_post_lo_GR$lo_GRs$Baltic_mean)

###Unreliable! (bugged HiC positions!)
#save(sprat_post_lo_GR, file = "~/Projects/Sprat/data/genotypes/Sprat_pool_freq_post_lo_GR_ext.RData")
#save(sprat_post_lo_GR, file = "~/Projects/Sprat/data/genotypes/Sprat_pool_freq_post_lo_GR_ext2.RData")
###

save(sprat_post_lo_GR, file = "~/Projects/Sprat/data/genotypes/Sprat_pool_freq_post_lo_GR_ext3.RData")

Sprat_lo_manhattan_v2(sprat_data = sprat_post_lo_GR$lo_GRs, scaff_list = names(Sprat_pin_curated)[1:50], daf_col = "brackish_v_BlackSea", plot_prefix = "~/Projects/Sprat/doc/Brackish_v_BlackSea_m1/Brackish_v_BlackSea_m1_replot_", site_filt_vec = sprat_site_filter_m1 )
Sprat_lo_manhattan_v2(sprat_data = sprat_post_lo_GR$lo_GRs, scaff_list = names(Sprat_pin_curated)[1:50], daf_col = "BS_LAND_v_Baltic", plot_prefix = "~/Projects/Sprat/doc/BS_LAND_v_Baltic_m1/BS_LAND_v_Baltic_m1_", site_filt_vec = sprat_site_filter_m1 )
Sprat_lo_manhattan_v2(sprat_data = sprat_post_lo_GR$lo_GRs, scaff_list = names(Sprat_pin_curated)[1:50], daf_col = "NS_v_S_Atl", plot_prefix = "~/Projects/Sprat/doc/NS_v_S_Atl_m1/NS_v_S_Atl_m1_", site_filt_vec = sprat_site_filter_m1 )
Sprat_lo_manhattan_v2(sprat_data = sprat_post_lo_GR$lo_GRs, scaff_list = names(Sprat_pin_curated)[1:50], daf_col = "SK_v_S_Atl", plot_prefix = "~/Projects/Sprat/doc/SK_v_S_Atl_m1/SK_v_S_Atl_m1_", site_filt_vec = sprat_site_filter_m1 )
Sprat_lo_manhattan_v2(sprat_data = sprat_post_lo_GR$lo_GRs, scaff_list = names(Sprat_pin_curated)[1:50], daf_col = "N_v_S_Atl", plot_prefix = "~/Projects/Sprat/doc/N_v_S_Atl_m1/N_v_S_Atl_m1_", site_filt_vec = sprat_site_filter_m1 )
Sprat_lo_manhattan_v2(sprat_data = sprat_post_lo_GR$lo_GRs, scaff_list = names(Sprat_pin_curated)[1:30], daf_col = "Oce_v_Baltic", plot_prefix = "~/Projects/Sprat/doc/Oce_v_Baltic_m1/Oce_v_Baltic_m1_replot_", site_filt_vec = sprat_site_filter_m1 )



tmp_filter_48_15Mb <- sprat_post_lo_GR$lo_GRs$HiC_scaffold == "PGA_scaffold_48__33_contigs__length_23346322" & sprat_site_filter_m1 & sprat_post_lo_GR$lo_GRs$BS_LAND_v_Baltic > 0.6 & sprat_post_lo_GR$lo_GRs$HiC_pos > 15.0e6 & sprat_post_lo_GR$lo_GRs$HiC_pos < 15.8e6
simple_sprat_pool_hm_v2(snp_data = sprat_post_lo_GR$lo_GRs, reg_filter = tmp_filter_48_15Mb, pdf_file =  "~/Projects/Sprat/doc/reg_heatmaps_m1/PGA_scaffold_48_15MB.pdf")

tmp_filter_48_17Mb <- sprat_post_lo_GR$lo_GRs$HiC_scaffold == "PGA_scaffold_48__33_contigs__length_23346322" & sprat_site_filter_m1 & sprat_post_lo_GR$lo_GRs$BS_LAND_v_Baltic > 0.7 & sprat_post_lo_GR$lo_GRs$HiC_pos > 16.0e6 & sprat_post_lo_GR$lo_GRs$HiC_pos < 18.0e6
simple_sprat_pool_hm_v2(snp_data = sprat_post_lo_GR$lo_GRs, reg_filter = tmp_filter_48_17Mb, pdf_file =  "~/Projects/Sprat/doc/reg_heatmaps_m1/PGA_scaffold_48_17MB.pdf")

tmp_filter_1106_8Mb <- sprat_post_lo_GR$lo_GRs$HiC_scaffold == "PGA_scaffold_1106__41_contigs__length_25158786" & sprat_site_filter_m1 & sprat_post_lo_GR$lo_GRs$BS_LAND_v_Baltic > 0.6  & sprat_post_lo_GR$lo_GRs$HiC_pos > 7.0e6 & sprat_post_lo_GR$lo_GRs$HiC_pos < 9.0e6
simple_sprat_pool_hm_v2(snp_data = sprat_post_lo_GR$lo_GRs, reg_filter = tmp_filter_1106_8Mb, pdf_file =  "~/Projects/Sprat/doc/reg_heatmaps_m1/PGA_scaffold_1106_8MB.pdf")

tmp_filter_1107_10_20Mb <- sprat_post_lo_GR$lo_GRs$HiC_scaffold == "PGA_scaffold_1107__54_contigs__length_31732007" & sprat_site_filter_m1 & sprat_post_lo_GR$lo_GRs$BS_LAND_v_Baltic > 0.6  & sprat_post_lo_GR$lo_GRs$HiC_pos > 10.0e6 & sprat_post_lo_GR$lo_GRs$HiC_pos < 20.0e6
simple_sprat_pool_hm_v2(snp_data = sprat_post_lo_GR$lo_GRs, reg_filter = tmp_filter_1107_10_20Mb, pdf_file =  "~/Projects/Sprat/doc/reg_heatmaps_m1/PGA_scaffold_1107_10_20MB.pdf")

reg_filter_1107_16Mb <- sprat_post_lo_GR$lo_GRs$HiC_scaffold == "PGA_scaffold_1107__54_contigs__length_31732007" & sprat_site_filter_m1  & sprat_post_lo_GR$lo_GRs$HiC_pos > 15.0e6 & sprat_post_lo_GR$lo_GRs$HiC_pos < 18.0e6
tmp_filter_1107_16Mb <- reg_filter_1107_16Mb & sprat_post_lo_GR$lo_GRs$Oce_v_Baltic > 0.7
simple_sprat_pool_hm_v2(snp_data = sprat_post_lo_GR$lo_GRs, reg_filter = tmp_filter_1107_16Mb, pdf_file =  "~/Projects/Sprat/doc/reg_heatmaps_m1/PGA_scaffold_1107_16MB.pdf")
plot(x = sprat_post_lo_GR$lo_GRs$HiC_pos[reg_filter_1107_16Mb], y = sprat_post_lo_GR$lo_GRs$Oce_v_Baltic[reg_filter_1107_16Mb],  pch = 16, cex = 0.5)


tmp_filter_1108_35Mb <- sprat_post_lo_GR$lo_GRs$HiC_scaffold == "PGA_scaffold_1108__43_contigs__length_36098062" & sprat_site_filter_m1 & sprat_post_lo_GR$lo_GRs$BS_LAND_v_Baltic > 0.7 & sprat_post_lo_GR$lo_GRs$HiC_pos > 30.0e6
simple_sprat_pool_hm_v2(snp_data = sprat_post_lo_GR$lo_GRs, reg_filter = tmp_filter_1108_35Mb, pdf_file =  "~/Projects/Sprat/doc/reg_heatmaps_m1/PGA_scaffold_1108_35MB.pdf")

tmp_filter_1108_23Mb <- sprat_post_lo_GR$lo_GRs$HiC_scaffold == "PGA_scaffold_1108__43_contigs__length_36098062" & sprat_site_filter_m1 & sprat_post_lo_GR$lo_GRs$N_v_S_Atl > 0.7 & sprat_post_lo_GR$lo_GRs$HiC_pos > 22.0e6
simple_sprat_pool_hm_v2(snp_data = sprat_post_lo_GR$lo_GRs, reg_filter = tmp_filter_1108_23Mb, pdf_file =  "~/Projects/Sprat/doc/reg_heatmaps_m1/PGA_scaffold_1108_23MB.pdf")




tmp_filter_1113_5Mb <- sprat_post_lo_GR$lo_GRs$HiC_scaffold == "PGA_scaffold_1113__35_contigs__length_27024755" & sprat_site_filter_m1 & sprat_post_lo_GR$lo_GRs$BS_LAND_v_Baltic > 0.7 & sprat_post_lo_GR$lo_GRs$HiC_pos > 5.0e6 & sprat_post_lo_GR$lo_GRs$HiC_pos < 6.0e6
simple_sprat_pool_hm_v2(snp_data = sprat_post_lo_GR$lo_GRs, reg_filter = tmp_filter_1113_5Mb, pdf_file =  "~/Projects/Sprat/doc/reg_heatmaps_m1/PGA_scaffold_1113_5MB.pdf")

tmp_filter_1118_44Mb <- sprat_post_lo_GR$lo_GRs$HiC_scaffold == "PGA_scaffold_1118__66_contigs__length_56234773" & sprat_site_filter_m1 & sprat_post_lo_GR$lo_GRs$BS_LAND_v_Baltic > 0.7 & sprat_post_lo_GR$lo_GRs$HiC_pos > 43.0e6
simple_sprat_pool_hm_v2(snp_data = sprat_post_lo_GR$lo_GRs, reg_filter = tmp_filter_1118_44Mb, pdf_file =  "~/Projects/Sprat/doc/reg_heatmaps_m1/PGA_scaffold_1118_44MB.pdf")

tmp_filter_1119_14Mb <- sprat_post_lo_GR$lo_GRs$HiC_scaffold == "PGA_scaffold_1119__31_contigs__length_21549996" & sprat_site_filter_m1 & sprat_post_lo_GR$lo_GRs$SK_v_S_Atl > 0.6 & sprat_post_lo_GR$lo_GRs$HiC_pos > 13.0e6 & sprat_post_lo_GR$lo_GRs$HiC_pos < 15.0e6
simple_sprat_pool_hm_v2(snp_data = sprat_post_lo_GR$lo_GRs, reg_filter = tmp_filter_1119_14Mb, pdf_file =  "~/Projects/Sprat/doc/reg_heatmaps_m1/PGA_scaffold_1119_14MB.pdf")


#Polishing names
raw_sprat_IDs <- sub("_Afreq", "", grep("_Afreq", names(sprat_post_lo_GR$lo_GRs@elementMetadata), value = T))

sprat_sample_info <- read.table("~/Projects/Sprat/doc/Sprat_pools_sequenced2021.txt", sep = "\t", header = T, stringsAsFactors = F)
#Naming
sprat_sample_info$clean_name <- sprat_sample_info$Site
sprat_sample_info$clean_name[c(1,2,6,7,8,9)] <- paste(sprat_sample_info$Site[c(1,2,6,7,8,9)], sprat_sample_info$Year[c(1,2,6,7,8,9)])
sprat_sample_info$clean_name[c(12,13)] <- paste(sprat_sample_info$Site[c(12,13)], c("North", "South"))
sprat_sample_info$clean_name[c(14)] <- sub(" fjord", "fjorden", sprat_sample_info$clean_name[14])
sprat_sample_info$clean_name <- sub(" $", "", sprat_sample_info$clean_name)
sprat_sample_info$clean_name[5] <- sub("\xbf", "\u00f8", sprat_sample_info$clean_name[5])
sprat_sample_info$clean_name[16] <- sub(" S", "", sprat_sample_info$clean_name[16])


#Combing coordinates
sprat_sample_info$lat_lon <- paste0(sprat_sample_info$Latitude, " N ", sprat_sample_info$Longitude, " E")
sprat_sample_info$lat_lon[sprat_sample_info$Longitude < 0] <- paste0(sprat_sample_info$Latitude, " N ", -1*sprat_sample_info$Longitude, " W")[sprat_sample_info$Longitude < 0]
write(sprat_sample_info$lat_lon, ncolumns = 1, file = "~/Projects/Sprat/doc/SRA_sub/lat_lons.txt")


#Plot positioning
sprat_sample_info$name_pos <- 4
sprat_sample_info$name_pos[c(12, 13, 16)] <- 1
sprat_sample_info$name_pos[c(1,6,8,19)] <- 2
sprat_sample_info$name_pos[c(15,7)] <- 3

save(sprat_sample_info, file = "~/Projects/Sprat/data/sprat_sample_info.RData")

#sprat_pool_names <- sprat_sample_info$clean_name[match(raw_sprat_IDs, sprat_sample_info$POOL)]

#simple_sprat_pool_hm_v3(snp_data = sprat_post_lo_GR$lo_GRs, reg_filter = tmp_filter_1119_14Mb, sample_info_df = sprat_sample_info, pdf_file =  "~/Projects/Sprat/doc/reg_heatmaps_m1/PGA_scaffold_1119_14MB_rename.pdf")

#Cleaned-up versions
simple_sprat_pool_hm_v3(snp_data = sprat_post_lo_GR$lo_GRs, reg_filter = tmp_filter_48_15Mb, sample_info_df = sprat_sample_info, pdf_file =  "~/Projects/Sprat/doc/reg_heatmaps_m1/PGA_scaffold_48_15MB_rename.pdf")
simple_sprat_pool_hm_v3(snp_data = sprat_post_lo_GR$lo_GRs, reg_filter = tmp_filter_1114_1, sample_info_df = sprat_sample_info, pdf_file = "~/Projects/Sprat/doc/reg_heatmaps_m1/PGA_scaffold_1114_6MB_rename.pdf")
#chr7_inv_filter<- sprat_post_lo_GR$lo_GRs$HiC_scaffold == m1_chr7_inv_main_scaff & sprat_site_filter_m1 & sprat_post_lo_GR$lo_GRs$oce_v_brack > 0.5 & sprat_post_lo_GR$lo_GRs$HiC_pos > 23e6
simple_sprat_pool_hm_v3(snp_data = sprat_post_lo_GR$lo_GRs, reg_filter = chr7_inv_filter, sample_info_df = sprat_sample_info, pdf_file =  "~/Projects/Sprat/doc/reg_heatmaps_m1/PGA_scaffold_13_25MB_rename.pdf")
#chr13_inv_filter <- sprat_post_lo_GR$lo_GRs$HiC_scaffold == m1_chr13_inv_main_scaff & sprat_site_filter_m1 & sprat_post_lo_GR$lo_GRs$oce_v_brack > 0.5
simple_sprat_pool_hm_v3(snp_data = sprat_post_lo_GR$lo_GRs, reg_filter = chr13_inv_filter, sample_info_df = sprat_sample_info, pdf_file =  "~/Projects/Sprat/doc/reg_heatmaps_m1/PGA_scaffold_173_rename.pdf")
#chr18_inv_filter <- sprat_post_lo_GR$lo_GRs$HiC_scaffold == "PGA_scaffold_40__17_contigs__length_6491151" & sprat_site_filter_m1 & sprat_post_lo_GR$lo_GRs$oce_v_brack > 0.7 & sprat_post_lo_GR$lo_GRs$HiC_pos > 1.5e6
simple_sprat_pool_hm_v3(snp_data = sprat_post_lo_GR$lo_GRs, reg_filter = chr18_inv_filter, sample_info_df = sprat_sample_info, pdf_file =  "~/Projects/Sprat/doc/reg_heatmaps_m1/PGA_scaffold_40_rename.pdf")

simple_sprat_pool_hm_v3(snp_data = sprat_post_lo_GR$lo_GRs, reg_filter = tmp_filter_1113, sample_info_df = sprat_sample_info, pdf_file =  "~/Projects/Sprat/doc/reg_heatmaps_m1/PGA_scaffold_1113_3MB_rename.pdf")
simple_sprat_pool_hm_v3(snp_data = sprat_post_lo_GR$lo_GRs, reg_filter = tmp_filter_4, sample_info_df = sprat_sample_info, pdf_file =  "~/Projects/Sprat/doc/reg_heatmaps_m1/PGA_scaffold_4_19MB_rename.pdf")
simple_sprat_pool_hm_v3(snp_data = sprat_post_lo_GR$lo_GRs, reg_filter = tmp_filter_69, sample_info_df = sprat_sample_info, pdf_file =  "~/Projects/Sprat/doc/reg_heatmaps_m1/PGA_scaffold_69_3MB_rename.pdf")
simple_sprat_pool_hm_v3(snp_data = sprat_post_lo_GR$lo_GRs, reg_filter = tmp_filter_1109, sample_info_df = sprat_sample_info, pdf_file =  "~/Projects/Sprat/doc/reg_heatmaps_m1/PGA_scaffold_1109_32MB_rename.pdf")
simple_sprat_pool_hm_v3(snp_data = sprat_post_lo_GR$lo_GRs, reg_filter = tmp_filter_1111, sample_info_df = sprat_sample_info, pdf_file =  "~/Projects/Sprat/doc/reg_heatmaps_m1/PGA_scaffold_1111_5MB_rename.pdf")
simple_sprat_pool_hm_v3(snp_data = sprat_post_lo_GR$lo_GRs, reg_filter = tmp_filter_1114_1, sample_info_df = sprat_sample_info, pdf_file =  "~/Projects/Sprat/doc/reg_heatmaps_m1/PGA_scaffold_1114_6MB_rename.pdf")
simple_sprat_pool_hm_v3(snp_data = sprat_post_lo_GR$lo_GRs, reg_filter = tmp_filter_1114_2, sample_info_df = sprat_sample_info, pdf_file =  "~/Projects/Sprat/doc/reg_heatmaps_m1/PGA_scaffold_1114_23MB_rename.pdf")
simple_sprat_pool_hm_v3(snp_data = sprat_post_lo_GR$lo_GRs, reg_filter = tmp_filter_1132, sample_info_df = sprat_sample_info, pdf_file =  "~/Projects/Sprat/doc/reg_heatmaps_m1/PGA_scaffold_1132_9MB_rename.pdf")
simple_sprat_pool_hm_v3(snp_data = sprat_post_lo_GR$lo_GRs, reg_filter = tmp_filter_13_13Mb, sample_info_df = sprat_sample_info, pdf_file =  "~/Projects/Sprat/doc/reg_heatmaps_m1/PGA_scaffold_13_13MB_rename.pdf")

simple_sprat_pool_hm_v3(snp_data = sprat_post_lo_GR$lo_GRs, reg_filter = tmp_filter_48_15Mb, sample_info_df = sprat_sample_info, pdf_file =  "~/Projects/Sprat/doc/reg_heatmaps_m1/PGA_scaffold_48_15MB_rename.pdf")
simple_sprat_pool_hm_v3(snp_data = sprat_post_lo_GR$lo_GRs, reg_filter = tmp_filter_48_17Mb, sample_info_df = sprat_sample_info, pdf_file =  "~/Projects/Sprat/doc/reg_heatmaps_m1/PGA_scaffold_48_17MB_rename.pdf")
simple_sprat_pool_hm_v3(snp_data = sprat_post_lo_GR$lo_GRs, reg_filter = tmp_filter_1106_8Mb, sample_info_df = sprat_sample_info, pdf_file =  "~/Projects/Sprat/doc/reg_heatmaps_m1/PGA_scaffold_1106_8MB_rename.pdf")
simple_sprat_pool_hm_v3(snp_data = sprat_post_lo_GR$lo_GRs, reg_filter = tmp_filter_1107_10_20Mb, sample_info_df = sprat_sample_info, pdf_file =  "~/Projects/Sprat/doc/reg_heatmaps_m1/PGA_scaffold_1107_10_20MB_rename.pdf")
simple_sprat_pool_hm_v3(snp_data = sprat_post_lo_GR$lo_GRs, reg_filter = tmp_filter_1107_16Mb, sample_info_df = sprat_sample_info, pdf_file =  "~/Projects/Sprat/doc/reg_heatmaps_m1/PGA_scaffold_1107_16MB_rename.pdf")
simple_sprat_pool_hm_v3(snp_data = sprat_post_lo_GR$lo_GRs, reg_filter = tmp_filter_1108_35Mb, sample_info_df = sprat_sample_info, pdf_file =  "~/Projects/Sprat/doc/reg_heatmaps_m1/PGA_scaffold_1108_35MB_rename.pdf")
simple_sprat_pool_hm_v3(snp_data = sprat_post_lo_GR$lo_GRs, reg_filter = tmp_filter_1108_23Mb, sample_info_df = sprat_sample_info, pdf_file =  "~/Projects/Sprat/doc/reg_heatmaps_m1/PGA_scaffold_1108_23MB_rename.pdf")


tmp_filter_1113_5Mb <- sprat_post_lo_GR$lo_GRs$HiC_scaffold == "PGA_scaffold_1113__35_contigs__length_27024755" & sprat_site_filter_m1 & sprat_post_lo_GR$lo_GRs$BS_LAND_v_Baltic > 0.7 & sprat_post_lo_GR$lo_GRs$HiC_pos > 5.0e6 & sprat_post_lo_GR$lo_GRs$HiC_pos < 6.0e6
simple_sprat_pool_hm_v2(snp_data = sprat_post_lo_GR$lo_GRs, reg_filter = tmp_filter_1113_5Mb, pdf_file =  "~/Projects/Sprat/doc/reg_heatmaps_m1/PGA_scaffold_1113_5MB.pdf")

tmp_filter_1118_44Mb <- sprat_post_lo_GR$lo_GRs$HiC_scaffold == "PGA_scaffold_1118__66_contigs__length_56234773" & sprat_site_filter_m1 & sprat_post_lo_GR$lo_GRs$BS_LAND_v_Baltic > 0.7 & sprat_post_lo_GR$lo_GRs$HiC_pos > 43.0e6
simple_sprat_pool_hm_v2(snp_data = sprat_post_lo_GR$lo_GRs, reg_filter = tmp_filter_1118_44Mb, pdf_file =  "~/Projects/Sprat/doc/reg_heatmaps_m1/PGA_scaffold_1118_44MB.pdf")

tmp_filter_1119_14Mb <- sprat_post_lo_GR$lo_GRs$HiC_scaffold == "PGA_scaffold_1119__31_contigs__length_21549996" & sprat_site_filter_m1 & sprat_post_lo_GR$lo_GRs$SK_v_S_Atl > 0.6 & sprat_post_lo_GR$lo_GRs$HiC_pos > 13.0e6 & sprat_post_lo_GR$lo_GRs$HiC_pos < 15.0e6
simple_sprat_pool_hm_v2(snp_data = sprat_post_lo_GR$lo_GRs, reg_filter = tmp_filter_1119_14Mb, pdf_file =  "~/Projects/Sprat/doc/reg_heatmaps_m1/PGA_scaffold_1119_14MB.pdf")

simple_sprat_pool_hm_v3(snp_data = sprat_post_lo_GR$lo_GRs, reg_filter = tmp_filter_1113_5Mb, sample_info_df = sprat_sample_info, pdf_file =  "~/Projects/Sprat/doc/reg_heatmaps_m1/PGA_scaffold_1113_5MB_rename.pdf")
simple_sprat_pool_hm_v3(snp_data = sprat_post_lo_GR$lo_GRs, reg_filter = tmp_filter_1118_44Mb, sample_info_df = sprat_sample_info, pdf_file =  "~/Projects/Sprat/doc/reg_heatmaps_m1/PGA_scaffold_1118_44MB_rename.pdf")
simple_sprat_pool_hm_v3(snp_data = sprat_post_lo_GR$lo_GRs, reg_filter = tmp_filter_1119_14Mb, sample_info_df = sprat_sample_info, pdf_file =  "~/Projects/Sprat/doc/reg_heatmaps_m1/PGA_scaffold_1119_14MB_rename.pdf")




require(maptools)
require(raster)
wrld_map <- readShapePoly("~/Projects/Herring/data/SNP_chip/SEC16B/10m-admin-0-countries/10m_admin_0_countries.shp")

pdf("~/Projects/Sprat/doc/Sprat_BS_loc.pdf", width = 5, height = 4, bg = "white")
plot(x = 15, y = 57, type = "n", xlim = c(28,41), ylim = c(40,46), xlab = "Longitude", ylab = "Latitude", main = "", cex.axis = 1.5)
plot(wrld_map, add = T, col = "grey90", border = "grey70")
points(x= sprat_sample_info[,"Longitude"], y = sprat_sample_info[,"Latitude"], col = "black", pch = 20, cex = 2)
text(x= sprat_sample_info[,"Longitude"], y = sprat_sample_info[,"Latitude"], labels = sprat_sample_info$clean_name, pos = sprat_sample_info$name_pos, cex = 1.1)
dev.off()

pdf("~/Projects/Sprat/doc/Sprat_Pool_loc.pdf", width = 12, height = 12)
plot(x = 15, y = 57, type = "n", xlim = c(-9,25), ylim = c(47,65), xlab = "Longitude", ylab = "Latitude", main = "", cex.axis = 2)
plot(wrld_map, add = T, col = "grey90", border = "grey70")
points(x= sprat_sample_info[,"Longitude"], y = sprat_sample_info[,"Latitude"], col = "black", pch = 20, cex = 2)
text(x= sprat_sample_info[,"Longitude"], y = sprat_sample_info[,"Latitude"], labels = sprat_sample_info$clean_name, pos = sprat_sample_info$name_pos, cex = 1.1)
dev.off()




#Support functions

simple_sprat_pool_hm_v3 <- function(snp_data, reg_filter, sample_info_df, pdf_file){
  tmp_hm_freq <- as.matrix(snp_data@elementMetadata[reg_filter, grep("_Afreq", names(snp_data@elementMetadata))])
  rownames(tmp_hm_freq) <- paste(snp_data$HiC_scaffold[reg_filter], snp_data$HiC_pos[reg_filter], sep = "_")
  tmp_sprat_IDs <- sub("_Afreq", "", grep("_Afreq", colnames(tmp_hm_freq), value = T))
  tmp_sprat_pool_names <- sample_info_df$clean_name[match(tmp_sprat_IDs, sample_info_df$POOL)]
  colnames(tmp_hm_freq) <- tmp_sprat_pool_names
  tmp_hm_freq <- tmp_hm_freq[order(snp_data$HiC_pos[reg_filter], decreasing = F),]
  pdf(file = pdf_file)
  heatmap(t(tmp_hm_freq), scale = "none", Colv = NA, margins = c(5,10))
  dev.off()
}

simple_sprat_pool_hm_v2 <- function(snp_data, reg_filter, pdf_file){
  tmp_hm_freq <- as.matrix(snp_data@elementMetadata[reg_filter, grep("_Afreq", names(snp_data@elementMetadata))])
  rownames(tmp_hm_freq) <- paste(snp_data$HiC_scaffold[reg_filter], snp_data$HiC_pos[reg_filter], sep = "_")
  tmp_hm_freq <- tmp_hm_freq[order(snp_data$HiC_pos[reg_filter], decreasing = F),]
  pdf(file = pdf_file)
  heatmap(t(tmp_hm_freq), scale = "none", Colv = NA, margins = c(5,10))
  dev.off()
}

Sprat_lo_manhattan_v2 <- function(sprat_data, scaff_list, daf_col, site_filt_vec, plot_prefix){
  daf_thresh_val <- qnorm(p = 1-(1/sum(site_filt_vec)), mean = mean(sprat_data@elementMetadata[site_filt_vec, daf_col]), sd  = sd(sprat_data@elementMetadata[site_filt_vec, daf_col]))
  for(chr in scaff_list){
    tmp_png_file <- paste0(plot_prefix,chr, ".png")
    tmp_chr_filter <- sprat_data$HiC_scaffold == chr
    png(filename = tmp_png_file, height = 1000, width = 2000)
    plot(x = sprat_data$HiC_pos[site_filt_vec & tmp_chr_filter], y = sprat_data@elementMetadata[site_filt_vec & tmp_chr_filter, daf_col], col = "black", ylim = c(0,1), pch = 16, cex = 0.8, xlab  = "Position", ylab = "DAF", main = chr)
    #tmp_reg_filter <- as.character(sel_reg_GR@seqnames) == chr
    #if(sum(tmp_reg_filter) > 0) rect(xleft = sel_reg_GR@ranges@start[tmp_reg_filter], xright = sel_reg_GR@ranges@start[tmp_reg_filter] + sel_reg_GR@ranges@width[tmp_reg_filter], ybottom = 0, ytop = 1, col = rgb(1,0,0,0.2), border = NULL)
    abline(h = daf_thresh_val, col = "red")
    #if(chr == "chr4") abline(v = c(11217660,11217516), col = "darkorchid") 
    dev.off()
  }
}

Sprat_v_herring_aln <- function(sprat_pri_GR, sprat_alt_GR, herring_GR, herring_samples = c("BM15_HastKar_Baltic_Spring", "Gavle54_Gavle_Baltic_Autumn", "NSSH34_Norway_Atlantic_Spring", "Z14_IsleofMan_Atlantic_Autumn", "HWS22_PechoraSea_BarentsSea", "HWS42_KandalakshaBay_WhiteSea_Spring", "HWS53_KandalakshaBay_WhiteSea_Summer", "HWS12_Japan_SeaofJapan", "HWS32_WhiteSea_WhiteSea", "HWS63_Balsfjord_Atlantic", "Pacific3_Vancouver_Pacific"), sprat_pri = IPA_pri, sprat_alt = IPA_alt, sprat_pri_dir = "pos", sprat_alt_dir = "pos"){
  require(Biostrings)
  require(GenomicRanges)
  require(ape)
  old_wd <-  getwd()
  setwd("~/Projects/Sprat/data/tree/")
  if(file.exists("./tmp_SvH_aln.fasta")){
    file.remove("./tmp_SvH_aln.fasta")
    print("Cleared old temp file!")
  } 
  for(herring_sample in herring_samples){
    system_cmd <- paste0("export PATH=$PATH:/Users/mapet205/Software/local/bin; samtools faidx ~/Projects/Herring/data/v2.0.2_reference/Ch_v2.0.2.fasta ", herring_GR, " | bcftools consensus -H 1 --sample ", herring_sample, " ~/Projects/Herring/data/v2.0.2_genotypes/phased_79_ind_v2.0.2.vcf.gz >> ./tmp_SvH_aln.fasta")
    system(system_cmd)
  }
  if(sprat_pri_dir == "pos") writeXStringSet(x = sprat_pri[sprat_pri_GR], file = "./tmp_SvH_aln.fasta", append = T)
  if(sprat_pri_dir == "neg") writeXStringSet(x = reverseComplement(sprat_pri[sprat_pri_GR]), file = "./tmp_SvH_aln.fasta", append = T)
  
  if(sprat_alt_dir == "pos") writeXStringSet(x = sprat_alt[sprat_alt_GR], file = "./tmp_SvH_aln.fasta", append = T)
  if(sprat_alt_dir == "neg") writeXStringSet(x = reverseComplement(sprat_alt[sprat_alt_GR]), file = "./tmp_SvH_aln.fasta", append = T)
  
  
  raw_seq_list <- read.FASTA("./tmp_SvH_aln.fasta")
  names(raw_seq_list)[1:length(herring_samples)] <- herring_samples
  if(sprat_pri_dir == "pos") names(raw_seq_list)[length(herring_samples)+1] <- paste(sprat_pri_GR, "plus", sep = "_")
  if(sprat_pri_dir == "neg") names(raw_seq_list)[length(herring_samples)+1] <- paste(sprat_pri_GR, "minus", sep = "_")
  if(sprat_alt_dir == "pos") names(raw_seq_list)[length(herring_samples)+2] <- paste(sprat_alt_GR, "plus", sep = "_")
  if(sprat_alt_dir == "neg") names(raw_seq_list)[length(herring_samples)+2] <- paste(sprat_alt_GR, "minus", sep = "_")
  
  
  setwd(old_wd)
  return(raw_seq_list)
}

herring_vcf_extract <- function(herring_GR, herring_samples = c("BM15_HastKar_Baltic_Spring", "Gavle54_Gavle_Baltic_Autumn", "NSSH34_Norway_Atlantic_Spring", "Z14_IsleofMan_Atlantic_Autumn", "HWS22_PechoraSea_BarentsSea", "HWS42_KandalakshaBay_WhiteSea_Spring", "HWS53_KandalakshaBay_WhiteSea_Summer", "HWS12_Japan_SeaofJapan", "HWS32_WhiteSea_WhiteSea", "HWS63_Balsfjord_Atlantic", "Pacific3_Vancouver_Pacific")){
  require(Biostrings)
  require(GenomicRanges)
  require(ape)
  old_wd <-  getwd()
  setwd("~/Projects/Sprat/data/tree/")
  if(file.exists("./tmp_H_extract.fasta")){
    file.remove("./tmp_H_extract.fasta")
    print("Cleared old temp file!")
  } 
  for(herring_sample in herring_samples){
    system_cmd <- paste0("export PATH=$PATH:/Users/mapet205/Software/local/bin; samtools faidx ~/Projects/Herring/data/v2.0.2_reference/Ch_v2.0.2.fasta ", herring_GR, " | bcftools consensus -H 1 --sample ", herring_sample, " ~/Projects/Herring/data/v2.0.2_genotypes/phased_79_ind_v2.0.2.vcf.gz >> ./tmp_H_extract.fasta")
    system(system_cmd)
  }

  raw_seq_list <- read.FASTA("./tmp_SvH_aln.fasta")
  names(raw_seq_list)[1:length(herring_samples)] <- herring_samples
  setwd(old_wd)
  return(raw_seq_list)
}

aln_basic_plots <- function(clustal_obj, pdf_file){
  #old.par <- par(no.readonly = T)
  png(filename = sub("pdf", "png", pdf_file), width = 1000, height = 1000)
  checkAlignment(clustal_obj)
  dev.off()
  #par(old.par)
  pdf(pdf_file)
  SvH_tree <- bionj(dist.dna(clustal_obj))
  plot.phylo(SvH_tree, type = "unrooted", lab4ut = "axial")
  SvH_tree_rooted <- root(SvH_tree, outgroup = grep("ctg", SvH_tree$tip.label, value = T), resolve.root = T)
  plot.phylo(SvH_tree_rooted)
  dev.off()
}

#Primary only version
Sprat_v_herring_aln_v1 <- function(sprat_GR, herring_GR, herring_samples = c("BM15_HastKar_Baltic_Spring", "Gavle54_Gavle_Baltic_Autumn", "NSSH34_Norway_Atlantic_Spring", "Z14_IsleofMan_Atlantic_Autumn", "HWS22_PechoraSea_BarentsSea", "HWS42_KandalakshaBay_WhiteSea_Spring", "HWS53_KandalakshaBay_WhiteSea_Summer", "HWS12_Japan_SeaofJapan", "HWS32_WhiteSea_WhiteSea", "HWS63_Balsfjord_Atlantic", "Pacific3_Vancouver_Pacific"), sprat_genome = IPA_pri, sprat_dir = "pos"){
  require(Biostrings)
  require(GenomicRanges)
  require(ape)
  old_wd <-  getwd()
  setwd("~/Projects/Sprat/data/tree/")
  if(file.exists("./tmp_SvH_aln.fasta")){
    file.remove("./tmp_SvH_aln.fasta")
    print("Cleared old temp file!")
  } 
  for(herring_sample in herring_samples){
    system_cmd <- paste0("export PATH=$PATH:/Users/mapet205/Software/local/bin; samtools faidx ~/Projects/Herring/data/v2.0.2_reference/Ch_v2.0.2.fasta ", herring_GR, " | bcftools consensus -H 1 --sample ", herring_sample, " ~/Projects/Herring/data/v2.0.2_genotypes/phased_79_ind_v2.0.2.vcf.gz >> ./tmp_SvH_aln.fasta")
    system(system_cmd)
  }
  if(sprat_dir == "pos") writeXStringSet(x = sprat_genome[sprat_GR], file = "./tmp_SvH_aln.fasta", append = T)
  if(sprat_dir == "neg") writeXStringSet(x = reverseComplement(sprat_genome[sprat_GR]), file = "./tmp_SvH_aln.fasta", append = T)
  
  raw_seq_list <- read.FASTA("./tmp_SvH_aln.fasta")
  names(raw_seq_list)[1:length(herring_samples)] <- herring_samples
  if(sprat_dir == "pos") names(raw_seq_list)[length(herring_samples)+1] <- paste(sprat_GR, "plus", sep = "_")
  if(sprat_dir == "neg") names(raw_seq_list)[length(herring_samples)+1] <- paste(sprat_GR, "minus", sep = "_")
  setwd(old_wd)
  return(raw_seq_list)
}


read_and_plot_fish_satsuma2 <- function(satsuma_file, fish_name, plot_dir = "~/Projects/Herring/doc/HiC_satsuma/other_fish/"){
  #fish_name <- sub(,,satsuma_file)
  fish_v_BGI <- read.table(satsuma_file, sep = "\t", stringsAsFactors=F, comment.char="")
  fish_v_BGI[,4] <- paste(sub(".+ovlk_", "", fish_v_BGI[,4]), sep = "")
  
  fish_col_vec <- as.integer(as.factor(fish_v_BGI[,1])) %% 8 + 1
  
  pdf(paste(plot_dir, fish_name, "_v_Herring.pdf", sep = ""))
  for(HiC_chr in unique(fish_v_BGI[,4])){
    fish_plot_filter <- fish_v_BGI[,4] == HiC_chr
    plot(y = fish_v_BGI[fish_plot_filter,2], x = fish_v_BGI[fish_plot_filter,5], pch = 16, col = fish_col_vec[fish_plot_filter], type = "p", main = HiC_chr, xlab = "Herring position", ylab = paste(fish_name, "position"))
  }
  dev.off()
  return(invisible(fish_v_BGI))
}

satsuma_synteny_blocks2 <- function(satsuma_file, pdf_file, gap_allowance = 5e5, name_clean_re = NULL, hit_co = 100){
  satsuma_in <- read.table(satsuma_file, stringsAsFactors=F, sep = "\t", comment.char = "")
  satsuma_in[,4] <- sub(".+ovlk_", "", satsuma_in[,4])
  if(!is.null(name_clean_re)){
    satsuma_in[,1] <- sub(name_clean_re, "\\1", satsuma_in[,1])
  }
  
  output_list <- list()
  ##split by chromsome
  pdf(pdf_file)
  q_chr_list <- names(table(satsuma_in[,1])[table(satsuma_in[,1]) >= hit_co])
  for(q_chr in q_chr_list){
    sub_satsuma <- satsuma_in[satsuma_in[,1] == q_chr,]
    #satsuma_target_GR <- GRanges(seqnames=sub_satsuma[,14], ranges=IRanges(start=satsuma_in[,5], end = satsuma_in[,6]))
    
    ##Detect major chromosme hit
    main_target <- names(table(sub_satsuma[,4]))[which.max(table(sub_satsuma[,4]))]
    sub_satsuma_main <- sub_satsuma[sub_satsuma[,4] == main_target,]
    ##Find blocks of matching hits
    block_counter <- 1
    sub_satsuma_main[1,"block_idx"] <- block_counter
    for(i in 2:dim(sub_satsuma_main)[1]){
      if(min(c(abs(sub_satsuma_main[i-1, 5] - sub_satsuma_main[i, 6]), abs(sub_satsuma_main[i-1, 6] - sub_satsuma_main[i, 5]))) > gap_allowance){
        block_counter <- block_counter + 1
      }
      sub_satsuma_main[i,"block_idx"] <- block_counter
      
    }
    ##Define borders
    block_extent_df <- data.frame(block_idx = unique(sub_satsuma_main[,"block_idx"]), start = NA, end = NA)
    for(idx in block_extent_df[,"block_idx"]){
      block_extent_df[block_extent_df[,"block_idx"] == idx, c("start", "end")] <- range(sub_satsuma_main[sub_satsuma_main[,"block_idx"] == idx,c(5,6)])
    }
    sub_satsuma[sub_satsuma[,4] == main_target,"block_idx"] <- sub_satsuma_main[,"block_idx"] 
    output_list[[q_chr]] <- list(t_chr = main_target, q_chr = q_chr, block_locations = block_extent_df, sub_satsuma = sub_satsuma)
    plot(y = sub_satsuma[,2], x = sub_satsuma[,5], pch = 16, col = as.integer(as.factor(sub_satsuma[,4])) %% 8 + 1, xlab = "Herring position", ylab = "Query position", type = "n", main = paste("Target: ", main_target, "; Query: ", q_chr, sep  =""), xlim = c(0, 4.0e7))
    rect(xleft = block_extent_df[,"start"], ybottom = 0, xright = block_extent_df[,"end"], ytop = 4e7, col = c("darkgrey", "lightgrey"))
    points(y = sub_satsuma[,2], x = sub_satsuma[,5], pch = 16, col = as.integer(as.factor(sub_satsuma[,4])) %% 8 + 1)
  }
  ##Return findings
  dev.off()
  
  return(output_list)
}


##To be adjusted
Sprat_HiC_liftover <- function(sprat_pos_data, liftover_df, chr_size_df, lo_cols){
  require(GenomicRanges)
  
  scaffold_GR <- GRanges(seqnames = sprat_pos_data$CHROM,ranges = IRanges(start = sprat_pos_data$POS, end = sprat_pos_data$POS))
  #BGI_v_Ilu_df
  liftover_GR <- GRanges(seqnames= liftover_df$sprat_IPA_seqnames,ranges=IRanges(start = liftover_df$sprat_IPA_start, end = liftover_df$sprat_IPA_end))
  
  #Matching SNP positions with entries in the BGI to Ilu satsuma alignment
  scaffold_matches <- findOverlaps(query = scaffold_GR, subject = liftover_GR, type = "any" )
  if(class(lo_cols) == "integer") lo_cols <- names(sprat_pos_data)[lo_cols]
  
  lo_df <- cbind(sprat_pos_data[scaffold_matches@from,c("CHROM","POS",lo_cols)], liftover_df[scaffold_matches@to,])
  
  #Adjusting positions within each matched interval
  target_SNPs <-  sign(lo_df$direction_est) == 1
  lo_df[target_SNPs,"SNP_HiC_pos"] <- lo_df[target_SNPs,"herring_v2.0.2_start"] + (lo_df[target_SNPs,"POS"] - lo_df[target_SNPs,"sprat_IPA_start"])
  target_SNPs <-  sign(lo_df$direction_est) != 1
  lo_df[target_SNPs,"SNP_HiC_pos"] <- lo_df[target_SNPs,"herring_v2.0.2_start"] + (lo_df[target_SNPs,"sprat_IPA_end"] - lo_df[target_SNPs,"POS"])
  
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
###

satsuma_processing_v2 <- function(satsuma_file, reduce_gap = 1e3){
  require(GenomicRanges)
  satsuma_in <- read.table(satsuma_file, stringsAsFactors=F, sep = "\t")
  satsuma_target_GR <- GRanges(seqnames=satsuma_in[,4], ranges=IRanges(start=satsuma_in[,5], end = satsuma_in[,6]))
  satsuma_query_GR <- GRanges(seqnames=satsuma_in[,1], ranges=IRanges(start=satsuma_in[,2], end = satsuma_in[,3]))
  satsuma_red_target_GR <- GenomicRanges::reduce(satsuma_target_GR, min.gapwidth = reduce_gap, with.revmap = T)
  revmap <- mcols(satsuma_red_target_GR)$revmap
  satsuma_relist_query_GR <- relist(satsuma_query_GR[unlist(revmap), ], revmap)
  satsuma_red_query_GR <- GenomicRanges::reduce(satsuma_relist_query_GR, min.gapwidth = reduce_gap)
  satsuma_query_df <- as.data.frame(satsuma_red_query_GR)
  #names(satsuma_target_df ) <- paste("q_", names(satsuma_target_df), sep = "")
  satsuma_target_df <- as.data.frame(satsuma_red_target_GR)
  names(satsuma_target_df ) <- paste("T_", names(satsuma_target_df), sep = "")
  
  satsuma_query_df <- satsuma_query_df[order(satsuma_query_df[,"group"], -satsuma_query_df[,"width"]),]
  satsuma_query_df <- satsuma_query_df[!duplicated(satsuma_query_df[,"group"]),]
  #satsuma_query_df[,"seqnames"] <- as.character(satsuma_query_df[,"seqnames"])
  satsuma_comb_df <- cbind(satsuma_query_df, satsuma_target_df[match(satsuma_query_df[,1], rownames(satsuma_target_df)),])
  satsuma_comb_df <- satsuma_comb_df[order(satsuma_comb_df[, "seqnames"], satsuma_comb_df[, "start"]),]
}

satsuma_direction_est <- function(satsuma_df){
  block_ends <- c(0,which(!(satsuma_df[(2:dim(satsuma_df)[1])-1, "seqnames"] == satsuma_df[2:dim(satsuma_df)[1], "seqnames"] & satsuma_df[(2:dim(satsuma_df)[1])-1, "T_seqnames"] == satsuma_df[2:dim(satsuma_df)[1], "T_seqnames"])), dim(satsuma_df)[1])
  for(i in 1:(length(block_ends)-1)){
    satsuma_df[(block_ends[i]+1):block_ends[i+1],"aln_block"] <- i
  }
  
  for(aln_block in unique(satsuma_df[,"aln_block"])){	
    aln_block_df <- satsuma_df[satsuma_df[,"aln_block"] == aln_block,]
    satsuma_df[satsuma_df[,"aln_block"] == aln_block,"direction_est"] <- sum(sign(diff(aln_block_df[,"start"])) * sign(diff(aln_block_df[,"T_start"])))
  }
  return(satsuma_df)
}

#AGP-based liftover function
agp_SNP_liftover <- function(lo_GRs, AGP_PB_GR){
  lo_hits <- findOverlaps(lo_GRs, AGP_PB_GR)
  lo_GRs$HiC_scaffold[lo_hits@from] <- AGP_PB_GR$HiC_scaffold[lo_hits@to]
  lo_GRs$HiC_start[lo_hits@from] <- AGP_PB_GR$HiC_start[lo_hits@to]
  lo_GRs$HiC_end[lo_hits@from] <- AGP_PB_GR$HiC_end[lo_hits@to]
  lo_GRs$direction[lo_hits@from] <- AGP_PB_GR$direction[lo_hits@to]
  lo_GRs$PB_start[lo_hits@from] <- AGP_PB_GR@ranges@start[lo_hits@to]
  
  lo_GRs$HiC_pos <- NA
  lo_GRs$HiC_pos[which(lo_GRs$direction == "+")] <- (lo_GRs$HiC_start + lo_GRs@ranges@start - lo_GRs$PB_start)[which(lo_GRs$direction == "+")]
  lo_GRs$HiC_pos[which(lo_GRs$direction == "-")] <- (lo_GRs$HiC_end - lo_GRs@ranges@start + lo_GRs$PB_start)[which(lo_GRs$direction == "-")]
  
  #lo_GRs$HiC_end_pos <- NA
  #lo_GRs$HiC_end_pos[which(lo_GRs$direction == "+")] <- (lo_gff$HiC_start + lo_gff@ranges@start + lo_gff@ranges@width - lo_gff$PB_start)[which(lo_gff$direction == "+")] - 1
  #lo_GRs$HiC_end_pos[which(lo_GRs$direction == "-")] <- (lo_gff$HiC_end - lo_gff@ranges@start -lo_gff@ranges@width + lo_gff$PB_start)[which(lo_gff$direction == "-")] + 1
  
  #hic_start_vec <- pmin(lo_gff$HiC_start_pos, lo_gff$HiC_end_pos)
  #hic_end_vec <- pmax(lo_gff$HiC_start_pos, lo_gff$HiC_end_pos)
  #hic_strand_vec <- rep("-", times = length(lo_gff))
  #hic_strand_vec[lo_gff$direction == "+" & as.character(lo_gff@strand) == "+"] <- "+"
  #hic_strand_vec[lo_gff$direction == "-" & as.character(lo_gff@strand) == "-"] <- "+"
  #hic_strand_vec <- factor(hic_strand_vec, levels = c("+", "-"))
  #hic_gff <- GRanges(seqnames = lo_gff$HiC_scaffold, ranges = IRanges(start = hic_start_vec , end = hic_end_vec), strand = hic_strand_vec)

  #hic_gff@elementMetadata <- lo_gff@elementMetadata
  
  return(list(lo_GRs = lo_GRs))
}

