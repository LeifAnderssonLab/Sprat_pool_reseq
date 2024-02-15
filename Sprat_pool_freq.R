#Author: Mats Pettersson
#mats.pettersson@imbim.uu.se

## For pooled data, AD counts
#zgrep "CHROM" Sprat.filtered.vcf.gz | awk -F'\t' '{$1="CHROM"; $3=""; $6=""; $8=""; $9=""; print}' > Sprat_pool_header.txt
#bcftools query -i 'TYPE=="SNP" & N_ALT==1' -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t[%AD\t]\n' ./Sprat.filtered.vcf.gz > Sprat_pool_AD.txt
#cat Sprat_pool_header.txt Sprat_pool_AD.txt > Sprat_pool_AD_header.txt
#zgrep -E "CHROM|PASS" Sprat_pool_AD_header.txt.gz > Sprat_pool_AD_PASS.txt

###New run - HiC_DeDup_v2
#zgrep -m1 "CHROM" ../../../variants/Sprat.filtered.vcf.gz | awk -F'\t' '{$1="CHROM"; $3=""; $6=""; $8=""; $9=""; print}' > Sprat_pool_header.txt
#bcftools query -i 'TYPE=="SNP" & N_ALT==1' -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t[%AD\t]\n' ../../../variants/Sprat.filtered.vcf.gz > Sprat_pool_AD.txt
#cat Sprat_pool_header.txt Sprat_pool_AD.txt > Sprat_pool_AD_header.txt
#mv Sprat_pool_AD_header.txt Sprat_pool_AD_DeDup_v2_raw.txt 
#zgrep -E "CHROM|PASS" Sprat_pool_AD_DeDup_v2_raw.txt.gz > Sprat_pool_AD_DeDup_v2_PASS.txt


#Sprat_pool_freq <- read.table("~/Projects/Sprat/data/genotypes/IPA_pri_genos/Sprat_pool_AD_PASS.txt.gz", header = T, sep = "", stringsAsFactors = F, fill = T) ## Original, IPA version
Sprat_pool_freq <- read.table("~/Projects/Sprat/data/genotypes/Sprat_pool_AD_DeDup_v2_PASS.txt.gz", header = T, sep = "", stringsAsFactors = F, fill = T) ## Current, HiC_DeDup_v2 version

sprat_pools <- names(Sprat_pool_freq)[6:24]

for(sprat_p in sprat_pools){
  Sprat_pool_freq[, paste0(sprat_p, "_A")] <- as.numeric(sub("([0-9]+)[,]([0-9]+)", "\\1", Sprat_pool_freq[,sprat_p]))
  Sprat_pool_freq[, paste0(sprat_p, "_D")] <- as.numeric(sub("([0-9]+)[,]([0-9]+)", "\\2", Sprat_pool_freq[,sprat_p]))
  Sprat_pool_freq[, paste0(sprat_p, "_Afreq")] <- Sprat_pool_freq[, paste0(sprat_p, "_A")]/(Sprat_pool_freq[, paste0(sprat_p, "_A")] + Sprat_pool_freq[, paste0(sprat_p, "_D")])
}

#save(Sprat_pool_freq, file = "~/Projects/Sprat/data/genotypes/Sprat_pool_freq.RData")
save(Sprat_pool_freq, file = "~/Projects/Sprat/data/genotypes/Sprat_pool_freq_DeDup_v2.RData")

for(sprat_p in sprat_pools){
  Sprat_pool_freq[, paste0(sprat_p, "_DP")] <- Sprat_pool_freq[, paste0(sprat_p, "_D")] + Sprat_pool_freq[, paste0(sprat_p, "_A")]
}


f3 <- Vectorize(FUN = function(X,Y, data_mat){sum(abs(data_mat[,X] - data_mat[,Y]), na.rm = T)/sum(!is.na(data_mat[,X]) & !is.na(data_mat[,Y]))}, vectorize.args = c("X", "Y"))
#Sprat_tree_SNPs <- sample.int(n = dim(Sprat_pool_freq)[1], size = 1e6, replace = F)
Sprat_tree_SNPs <- which(sprat_site_filter)
Sprat_Afreq_col <- grep("_Afreq", names(Sprat_pool_freq))
Sprat_raw_dist <- outer(X = Sprat_Afreq_col, Y = Sprat_Afreq_col, FUN = "f3", data_mat = Sprat_pool_freq[Sprat_tree_SNPs,])
rownames(Sprat_raw_dist) <- names(Sprat_pool_freq)[Sprat_Afreq_col]
colnames(Sprat_raw_dist) <- names(Sprat_pool_freq)[Sprat_Afreq_col]

#save(Sprat_raw_dist, file = "~/Projects/Sprat/data/genotypes/IPA_pri_genos/Sprat_pool_raw_dist.RData") ## Original, IPA_pri version
save(Sprat_raw_dist, file = "~/Projects/Sprat/data/genotypes/Sprat_pool_raw_dist_DeDup_v2.RData") ## Current, HiC_DeDup_v2 version


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


png(file = "~/Projects/Sprat/doc/pool_tree_DeDup_v2.png", height = 1000, width = 1000)
plot.phylo(Sprat_tree, type="unrooted", lab4ut = "axial", cex = 3, edge.width = 2.5)
add.scale.bar(length = NULL, ask = FALSE,
              lwd = 1, lcol = "black")
dev.off()






