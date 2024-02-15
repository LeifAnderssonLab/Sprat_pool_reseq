#Author: Mats Pettersson
#mats.pettersson@imbim.uu.se

#The two haplotype-level assemblies
IPA_pri <- readDNAStringSet("/Users/mapet205/Projects/Sprat/data/assemblies/IPA/final.p_ctg.fasta.gz")
IPA_alt <- readDNAStringSet("/Users/mapet205/Projects/Sprat/data/assemblies/IPA/final.a_ctg.fasta.gz")

#Milestone1 vs IPA alternative
Sprat_pin_curated <- readDNAStringSet("~/Projects/Sprat/data/HiC_satsuma/pin_HiC_it6/curation/Sprat_pin_milestone1.fasta.gz")
Sprat_pin_curated <- Sprat_pin_curated[order(Sprat_pin_curated@ranges@width, decreasing = T)]
IPA_alt_v_M1_satsuma <- read.table("~/Projects/Sprat/data/HiC_satsuma/Sprat_pin_m1_v_IPA_alt_satsuma.out_sorted.gz", stringsAsFactors=F, sep = "\t", comment.char = "")





#Evaluation region
ctg.000003F_blastout <- read.table("~/Projects/Sprat/data/IPA_blast/ctg.000003F_0_150kb.blastout", sep = "\t", stringsAsFactors = F)
ctg.000003F_hap_blast <- ctg.000003F_blastout[grepl("hap_", ctg.000003F_blastout[,2]),]
head(ctg.000003F_hap_blast[order(ctg.000003F_hap_blast[,4], decreasing = T),])

pre_aln_ctgs <- c(IPA_pri["ctg.000003F"], IPA_alt["hap_ctg.002952F HAPLOTIG"]) #Selected from the dominant blast hits


Sp1_seqnames <- c("ctg.000003F", "hap_ctg.002952F HAPLOTIG")
Sp1_starts <- c(1, 1)
Sp1_ends <- c(88000, 88623)

Sp1_clustalo <- align_target_region_Sp(seqname_vec = Sp1_seqnames, start_vec = Sp1_starts, stop_vec = Sp1_ends, genome_set = pre_aln_ctgs, dir_vec = c(1,-1))

pdf(file = "~/Projects/Sprat/doc/diversity/ctg.000003F_0_88kb.pdf", width = 14, height = 7)
checkAlignment(Sp1_clustalo)
dev.off()
dist.dna(Sp1_clustalo)
#aln_diff_plot(reg1_clustalo)
Sp1_diff <- aln_diff_plot_Sp(Sp1_clustalo, win_size = 500, pdf_file = "~/Projects/Sprat/doc/diversity/ctg.000003F_0_88kb_diff.pdf", y_lim = c(-0.5,0.1))

#Using herring "translation" to find corresponding contigs
IPA_pri_satsuma <- read.table("~/Projects/Sprat/data/PacBio_satsuma/Sprat_IPA_primary_v_Ch_v2.0.2.chained.out_sorted.gz", stringsAsFactors=F, sep = "\t", comment.char = "")
IPA_alt_satsuma <- read.table("~/Projects/Sprat/data/PacBio_satsuma/Sprat_IPA_alternative_v_Ch_v2.0.2.chained.out_sorted.gz", stringsAsFactors=F, sep = "\t", comment.char = "")


reg_filter <- IPA_pri_satsuma$V4 == "chr11" & IPA_pri_satsuma$V5 > 8.8e6 & IPA_pri_satsuma$V5 < 9.0e6
table(IPA_pri_satsuma[reg_filter, "V1"])
plot(y = IPA_pri_satsuma[reg_filter, "V2"], x = IPA_pri_satsuma[reg_filter, "V5"])

reg_filter_alt <- IPA_alt_satsuma$V4 == "chr11" & IPA_alt_satsuma$V5 > 8.8e6 & IPA_alt_satsuma$V5 < 9.0e6
table(IPA_alt_satsuma[reg_filter_alt, "V1"])
plot(y = IPA_alt_satsuma[reg_filter_alt, "V2"], x = IPA_alt_satsuma[reg_filter_alt, "V5"])


Sp2_seqnames <- c("ctg.000470F", "hap_ctg.000688F HAPLOTIG")
Sp2_starts <- c(1.5e5, 1.5e5)
Sp2_ends <- c(2e5, 2e5)
pre_aln_ctgs_2 <- c(IPA_pri["ctg.000470F"], IPA_alt["hap_ctg.000688F HAPLOTIG"]) #Selected from the dominant satsuma hits
Sp2_clustalo <- align_target_region_Sp(seqname_vec = Sp2_seqnames, start_vec = Sp2_starts, stop_vec = Sp2_ends, genome_set = pre_aln_ctgs_2, dir_vec = c(1,1))

checkAlignment(Sp2_clustalo)
Sp2_diff <- aln_diff_plot_Sp(Sp2_clustalo, win_size = 500, pdf_file = "~/Projects/Sprat/doc/diversity/ctg.000470F_150_200kb_diff.pdf", y_lim = c(-0.5,0.5))

reg_filter <- IPA_pri_satsuma$V4 == "chr3" & IPA_pri_satsuma$V5 > 26.9e6 & IPA_pri_satsuma$V5 < 27.0e6
table(IPA_pri_satsuma[reg_filter, "V1"])
plot(y = IPA_pri_satsuma[reg_filter, "V2"], x = IPA_pri_satsuma[reg_filter, "V5"])

#reg_filter_alt <- IPA_alt_satsuma$V4 == "chr3" & IPA_alt_satsuma$V5 > 26.9e6 & IPA_alt_satsuma$V5 < 27.0e6
#table(IPA_alt_satsuma[reg_filter_alt, "V1"])
#plot(y = IPA_alt_satsuma[reg_filter_alt, "V2"], x = IPA_alt_satsuma[reg_filter_alt, "V5"])
#Both copies in primary

aln_ctgs_4 <- c(IPA_pri[GRanges("ctg.000018F:5.50e5-6.00e5")], IPA_pri[GRanges("ctg.000030F:23.95e5-24.45e5")]) #Selected from the dominant satsuma hits
aln_ctgs_4
Sp4_clustalo <- align_target_region_Sp(aln_seq = aln_ctgs_4, dir_vec = c(-1,1))
checkAlignment(Sp4_clustalo)
Sp4_diff <- aln_diff_plot_Sp(Sp4_clustalo, win_size = 500, pdf_file = "~/Projects/Sprat/doc/diversity/ctg.000018F_550_600kb_diff.pdf", y_lim = c(-1,0.1))



##Trying regions from inside inversions

##Chr18
reg_filter <- IPA_pri_satsuma$V4 == "chr18" & IPA_pri_satsuma$V5 > 14.8e6 & IPA_pri_satsuma$V5 < 15.0e6
table(IPA_pri_satsuma[reg_filter, "V1"])
plot(y = IPA_pri_satsuma[reg_filter, "V2"], x = IPA_pri_satsuma[reg_filter, "V5"])

reg_filter_alt <- IPA_alt_satsuma$V4 == "chr18" & IPA_alt_satsuma$V5 > 14.8e6 & IPA_alt_satsuma$V5 < 15.0e6
table(IPA_alt_satsuma[reg_filter_alt, "V1"])
plot(y = IPA_alt_satsuma[reg_filter_alt, "V2"], x = IPA_alt_satsuma[reg_filter_alt, "V5"])

aln_ctgs_3 <- c(IPA_pri[GRanges("ctg.000969F:2.27e5-2.73e5")], IPA_alt[GRanges("hap_ctg.002173F HAPLOTIG:1.0e5-1.455e5")]) #Selected from the dominant satsuma hits
aln_ctgs_3
Sp3_clustalo <- align_target_region_Sp(aln_seq = aln_ctgs_3, dir_vec = c(-1,1))
checkAlignment(Sp3_clustalo)
Sp3_diff <- aln_diff_plot_Sp(Sp3_clustalo, win_size = 500, pdf_file = "~/Projects/Sprat/doc/diversity/ctg.000969F_227_273kb_diff.pdf", y_lim = c(-1,0.1))

### Chr5
reg_filter <- IPA_pri_satsuma$V4 == "chr5" & IPA_pri_satsuma$V5 > 11.9e6 & IPA_pri_satsuma$V5 < 12.1e6
table(IPA_pri_satsuma[reg_filter, "V1"])
plot(y = IPA_pri_satsuma[reg_filter, "V2"], x = IPA_pri_satsuma[reg_filter, "V5"])

reg_filter_alt <- IPA_alt_satsuma$V4 == "chr5" & IPA_alt_satsuma$V5 > 11.9e6 & IPA_alt_satsuma$V5 < 12.1e6
table(IPA_alt_satsuma[reg_filter_alt, "V1"])
plot(y = IPA_alt_satsuma[reg_filter_alt, "V2"], x = IPA_alt_satsuma[reg_filter_alt, "V5"])

aln_ctgs_5 <- c(IPA_pri[GRanges("ctg.001169F:1.00e5-1.50e5")], IPA_alt[GRanges("hap_ctg.001595F HAPLOTIG:1.0e5-1.50e5")]) #Selected from the dominant satsuma hits
aln_ctgs_5
Sp5_clustalo <- align_target_region_Sp(aln_seq = aln_ctgs_5, dir_vec = c(1,1))
checkAlignment(Sp5_clustalo)
Sp5_diff <- aln_diff_plot_Sp(Sp5_clustalo, win_size = 500, pdf_file = "~/Projects/Sprat/doc/diversity/ctg.001169F_100_150kb_diff.pdf", y_lim = c(-1,0.1))

### Chr4
reg_filter <- IPA_pri_satsuma$V4 == "chr4" & IPA_pri_satsuma$V5 > 28.6e6 & IPA_pri_satsuma$V5 < 28.8e6
table(IPA_pri_satsuma[reg_filter, "V1"])
plot(y = IPA_pri_satsuma[reg_filter, "V2"], x = IPA_pri_satsuma[reg_filter, "V5"])

reg_filter_alt <- IPA_alt_satsuma$V4 == "chr4" & IPA_alt_satsuma$V5 > 28.6e6 & IPA_alt_satsuma$V5 < 28.8e6
table(IPA_alt_satsuma[reg_filter_alt, "V1"])
plot(y = IPA_alt_satsuma[reg_filter_alt, "V2"], x = IPA_alt_satsuma[reg_filter_alt, "V5"])

aln_ctgs_6 <- c(IPA_pri[GRanges("ctg.001598F:0.01e5-0.50e5")], IPA_alt[GRanges("hap_ctg.001802F HAPLOTIG:0.01e5-0.50e5")]) #Selected from the dominant satsuma hits
aln_ctgs_6
Sp6_clustalo <- align_target_region_Sp(aln_seq = aln_ctgs_6, dir_vec = c(1,1))
checkAlignment(Sp6_clustalo)
Sp6_diff <- aln_diff_plot_Sp(Sp6_clustalo, win_size = 500, pdf_file = "~/Projects/Sprat/doc/diversity/ctg.001598F_1_50kb_diff.pdf", y_lim = c(-1,0.1))


reg_filter <- IPA_pri_satsuma$V4 == "chr4" & IPA_pri_satsuma$V5 > 31.8e6 & IPA_pri_satsuma$V5 < 32.0e6
table(IPA_pri_satsuma[reg_filter, "V1"])
plot(y = IPA_pri_satsuma[reg_filter, "V2"], x = IPA_pri_satsuma[reg_filter, "V5"])

reg_filter_alt <- IPA_alt_satsuma$V4 == "chr4" & IPA_alt_satsuma$V5 > 31.8e6 & IPA_alt_satsuma$V5 < 32.0e6
table(IPA_alt_satsuma[reg_filter_alt, "V1"])
plot(y = IPA_alt_satsuma[reg_filter_alt, "V2"], x = IPA_alt_satsuma[reg_filter_alt, "V5"])

aln_ctgs_7 <- c(IPA_pri[GRanges("ctg.002332F:0.01e5-0.50e5")], IPA_alt[GRanges("hap_ctg.002422F HAPLOTIG:0.01e5-0.50e5")]) #Selected from the dominant satsuma hits
aln_ctgs_7
Sp7_clustalo <- align_target_region_Sp(aln_seq = aln_ctgs_7, dir_vec = c(1,1))
checkAlignment(Sp7_clustalo)
Sp7_diff <- aln_diff_plot_Sp(Sp7_clustalo, win_size = 500, pdf_file = "~/Projects/Sprat/doc/diversity/ctg.002332F_1_50kb_diff.pdf", y_lim = c(-1,0.1))
####

dim(IPA_alt_v_M1_satsuma[IPA_alt_v_M1_satsuma$V4 == "PGA_scaffold_40__17_contigs__length_6491151" & IPA_alt_v_M1_satsuma$V5 >= 1763627 & IPA_alt_v_M1_satsuma$V6 <= 3165171, ]) #Corresponds to (Herring) Chr 18 inversion
#[1] 8812    8

sp_reg_satsuma <- IPA_alt_v_M1_satsuma[IPA_alt_v_M1_satsuma$V4 == "PGA_scaffold_40__17_contigs__length_6491151" & IPA_alt_v_M1_satsuma$V5 >= 1785000 & IPA_alt_v_M1_satsuma$V6 <= 1763627 + 1e5, ] #Alignable block
IPA_alt_reg_GR <- GRanges(seqnames = sub("(.+)_([A-Z]+$)", "\\1 \\2", unique(sp_reg_satsuma$V1)), ranges = IRanges(start = min(sp_reg_satsuma$V2), end = max(sp_reg_satsuma$V3)))
M1_reg_GR <- GRanges(seqnames = unique(sp_reg_satsuma$V4), ranges = IRanges(start = min(sp_reg_satsuma$V5), end = max(sp_reg_satsuma$V6)))
aln_ctgs_satsuma <- c(Sprat_pin_curated[M1_reg_GR], IPA_alt[IPA_alt_reg_GR]) #Selected from the dominant satsuma hits
sp_reg_satsuma_clustalo <- align_target_region_Sp(aln_seq = aln_ctgs_satsuma, dir_vec = c(1,1))
checkAlignment(sp_reg_satsuma_clustalo)
sp_reg_satsuma_diff <- aln_diff_plot_Sp(sp_reg_satsuma_clustalo, win_size = 100, pdf_file = "~/Projects/Sprat/doc/diversity/sp_reg_satsuma_diff.pdf", y_lim = c(-1,0.1))

#More general approach
sp_reg_EnA <- extract_and_align_Sp(satsuma_df = IPA_alt_v_M1_satsuma, 
                                   target_scaff = "PGA_scaffold_40__17_contigs__length_6491151", 
                                   reg_start_pos = 1905000, reg_size = 0.4e5,
                                   HiC_genome = Sprat_pin_curated, query_genome = IPA_alt)

rnd_win_size <- 0.25e5
main_scaff_vec <- names(table(IPA_alt_v_M1_satsuma$V4))[order(table(IPA_alt_v_M1_satsuma$V4), decreasing = T)][1:75]

rnd_scaff_vec <- sample(main_scaff_vec, size = 100, replace = T)
aln_list <- list()
for(rnd_scaff in rnd_scaff_vec){
  if(width(Sprat_pin_curated[rnd_scaff]) - (rnd_win_size + 1) > 0){
    rnd_start <- sample.int(n = width(Sprat_pin_curated[rnd_scaff]) - (rnd_win_size + 1), size = 1)
    
    aln_list[[paste(rnd_scaff, rnd_start, sep = "_")]] <- extract_and_align_Sp(satsuma_df = IPA_alt_v_M1_satsuma, 
                                                                               target_scaff = rnd_scaff, 
                                                                               reg_start_pos = rnd_start, reg_size = rnd_win_size,
                                                                               HiC_genome = Sprat_pin_curated, query_genome = IPA_alt) 
  }else {
    print("Scaff too small!")
  }
}

sp_avg_div <- data.frame(block = names(aln_list), obs_div = NA, aln_length = NA)
for(aln_idx in 1:length(aln_list)){
  sp_avg_div[aln_idx, "aln_length"] <- sum(as.character(aln_list[[aln_idx]][2,]) != "-" & as.character(aln_list[[aln_idx]][1,]) != "-")
  sp_avg_div[aln_idx, "obs_div"] <- dist.dna(aln_list[[aln_idx]])
}


save(sp_avg_div, aln_list, file = "~/Projects/Sprat/data/diversity/sp_avg_div_run6.RData") #Filled slots: 6, 5, 4, 3, 2, 1

#sp_mean_obs_div <- sum((sp_avg_div$aln_length*sp_avg_div$obs_div))/sum(sp_avg_div$aln_length)
run_filter <- sp_avg_div$obs_div < 0.05 #Suppressing failed alignments
sp_mean_obs_div <- sum((sp_avg_div$aln_length[run_filter]*sp_avg_div$obs_div[run_filter]))/sum(sp_avg_div$aln_length[run_filter])
sp_mean_obs_div*100
sp_median_obs_div <- median(sp_avg_div$obs_div[run_filter])
sp_median_obs_div*100

#Compiling results
for(i in 1:6){
  div_df_file <- paste0("~/Projects/Sprat/data/diversity/sp_avg_div_run",i,".RData")
  load(div_df_file)
  if(i == 1) sp_div_compilation <- sp_avg_div
  if(i != 1){
    sp_div_compilation <- rbind(sp_div_compilation,sp_avg_div)
  }
}
compilation_filter <- sp_div_compilation$obs_div < 0.05 #Suppressing failed alignments
sp_mean_obs_div <- sum((sp_div_compilation$aln_length[compilation_filter]*sp_div_compilation$obs_div[compilation_filter]))/sum(sp_div_compilation$aln_length[compilation_filter])
sp_mean_obs_div*100

pdf(file = "~/Projects/Sprat/doc/diversity/compiled_region_hist.pdf")
hist(sp_div_compilation$obs_div[compilation_filter], breaks = 50, main = "Divergence between Pimary and Alternative contigs; 25 kb blocks", xlab = "Divergence", col = "salmon") 
abline(v = 0.003, col = "darkorchid", lty = "dashed", lwd = 2)
abline(v = median(sp_div_compilation$obs_div[compilation_filter]), lwd = 2, col = "grey70")
#abline(v = median(sp_div_compilation$obs_div), lwd = 2, col = "grey40")
abline(v = sp_mean_obs_div, lwd = 2, col = "black")
legend(x = "topright", lwd = 2, col = c("darkorchid","grey70","black"), lty = c(2,1,1), legend = c("Herring", "Median", "Weighted mean"))
dev.off()

extract_and_align_Sp <- function(satsuma_df, target_scaff, reg_start_pos, reg_size = 1e5, HiC_genome, query_genome){
  reg_satsuma <- satsuma_df[satsuma_df$V4 == target_scaff & satsuma_df$V5 >= reg_start_pos & IPA_alt_v_M1_satsuma$V6 <= reg_start_pos + reg_size, ] #Alignable block
  main_query <- names(table(reg_satsuma$V1))[order(table(reg_satsuma$V1), decreasing = T)][1]
  reg_satsuma <- reg_satsuma[reg_satsuma$V1 == main_query, ]
    
  #query_reg_GR <- GRanges(seqnames = sub("(.+)_([A-Z]+$)", "\\1 \\2", unique(reg_satsuma$V1)), ranges = IRanges(start = min(reg_satsuma$V2), end = max(reg_satsuma$V3)))
  aln_present <- F
  aln_ctgs <- NULL
  if(dim(reg_satsuma)[1] > 10){ #Checking for the amount of satsuma hits 
    query_reg_GR <- GRanges(seqnames = gsub("([^p])_", "\\1 ", unique(reg_satsuma$V1)), ranges = IRanges(start = min(reg_satsuma$V2), end = max(reg_satsuma$V3)))
    HiC_reg_GR <- GRanges(seqnames = unique(reg_satsuma$V4), ranges = IRanges(start = min(reg_satsuma$V5), end = max(reg_satsuma$V6)))
    aln_ctgs <- try(c(HiC_genome[HiC_reg_GR], query_genome[query_reg_GR])) #Selected from the dominant satsuma hits
    print(class(aln_ctgs))
  #plot(x = reg_satsuma$V5, y = reg_satsuma$V2)
    dir_vec = c(1,1)
    if(sign(lm(formula =  reg_satsuma$V2 ~ reg_satsuma$V5)$coef["reg_satsuma$V5"]) == -1) dir_vec = c(-1,1)
    print(lm(formula =  reg_satsuma$V2 ~ reg_satsuma$V5))
    print(paste("Aln directions. HiC:", dir_vec[1], "; Query:", dir_vec[2]))
    aln_present <- T
  }
  
  if(reg_size <= 1e5 & aln_present & class(aln_ctgs) == "DNAStringSet"){ 
    reg_clustalo <- align_target_region_Sp(aln_seq = aln_ctgs, dir_vec = dir_vec)
    chkAln <- checkAlignment(reg_clustalo)
    print(paste("Observed divergence:", round(dist.dna(reg_clustalo)*100, digits = 1), "%"))
    #sp_reg_satsuma_diff <- aln_diff_plot_Sp(sp_reg_satsuma_clustalo, win_size = 100, pdf_file = "~/Projects/Sprat/doc/diversity/sp_reg_satsuma_diff.pdf", y_lim = c(-1,0.1))
  } else{
    print(paste0("Region too large and/or inconclusive! Size: ", reg_size, "; Region in satsuma file; ", aln_present, "; Target contigs: ", gsub("([^p])_", "\\1 ", unique(reg_satsuma$V1)), " & ", unique(reg_satsuma$V4)))
    reg_clustalo <- NULL
  }
  return(reg_clustalo)
}


#Support functions
align_target_region_Sp <- function(aln_seq, dir_vec){
  require(ape)
  require(GenomicRanges)
  #aln_target <- GRanges(seqnames = seqname_vec, ranges = IRanges(start = start_vec, end = stop_vec))
  #aln_seq <- genome_set[aln_target]
  #names(aln_seq) <- paste(aln_target)
  for(i in 1:length(dir_vec)){
    if(dir_vec[i] == -1) {aln_seq[i] <- reverseComplement(aln_seq[i])}
  }
  aln_seq_bin <- as.DNAbin(aln_seq)
  reg_clustalo <- clustalomega(aln_seq_bin, exec="clustalo", quiet = F)
  #checkAlignment(reg_clustalo)
  return(reg_clustalo)
}

aln_diff_plot_Sp <- function(aln_obj, win_size = 1000, pdf_file = NULL, y_lim = c(-1,1), ...){
  aln_char_us <- as.character(aln_obj[1,])
  aln_char_chr <- as.character(aln_obj[2,])
  
  aln_pos <- which(aln_char_chr != "-" & aln_char_us != "-")
  gap_pos <- aln_char_chr == "-" | aln_char_us == "-"
  aln_pos_clustalo <- aln_obj[,aln_pos]
  #aln_pos_chr_gaps <- as.numeric(as.character(aln_pos_chr_clustalo[1,]) == "-")
  avg_gaps <- stats::filter(as.numeric(gap_pos), rep(1/win_size, win_size))
  
  #double_aln_pos <- which(as.character(aln_pos_chr_clustalo[1,]) != "-")
  diff_vec <- as.character(aln_obj[1,aln_pos]) != as.character(aln_obj[2,aln_pos])
  aln_avg_diff <- stats::filter(as.numeric(diff_vec), rep(1/win_size, win_size))
  
  #avg_diff_vec <- filter(as.numeric(diff_vec), rep(1/win_size, win_size))
  
  if(!is.null(pdf_file)) pdf(file = pdf_file, width = 10, height = 7)
  #if(is.null(y_lim)) y_lim <- c(-1,1)
  plot(x= aln_pos, y = aln_avg_diff, ylim = y_lim, type = "n")
  #gtf_tmp <- bmpr1ba_gtf[1:length(bmpr1ba_gtf) %in% bmpr1ba_gtf_subset & bmpr1ba_gtf$type == "CDS"]
  #rect(xleft = gtf_tmp@ranges@start, xright = gtf_tmp@ranges@start + gtf_tmp@ranges@width, ybottom = -1.5, ytop = 1.5, border = NA, lwd  =.5, col = "grey80")
  abline(h = 0, col = "grey70")
  abline(h = seq(from = 0.1, to = 0.5, by = 0.1), col = "grey80")
  points(x= 1:length(aln_obj[1,]) , y = -avg_gaps, pch = 16, cex = 0.5, col = "olivedrab")
  points(x= aln_pos, y = aln_avg_diff, pch = 16, cex = 0.5, col = "darkorchid")
  #plot(y = avg_diff_vec, x = 1:length(avg_diff_vec), xlab = "", ylab = "Difference (fraction)", ... )
  
  if(!is.null(pdf_file)) dev.off()
  return(invisible(list(avg_diff = aln_avg_diff, avg_gaps = avg_gaps, gap_pos = gap_pos, aln_pos = aln_pos)))
}
