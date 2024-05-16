#Author: Mats Pettersson
#mats.pettersson@imbim.uu.se

#IPA_pri <- readDNAStringSet("/Users/mapet205/Projects/Sprat/data/assemblies/IPA/final.p_ctg.fasta.gz")

mdm2_blastout_raw <- read.table("/Users/mapet205/Projects/Sprat/data/IPA_blast/Ch_MDM2vSprat_IPA_nt.blastout", sep = "\t", stringsAsFactors = F)
mdm2_blastout <- mdm2_blastout_raw[mdm2_blastout_raw[,3] > 90 & mdm2_blastout_raw[,4] > 50, ]
CluHar_mdm <- readDNAStringSet("/Users/mapet205/Projects/Sprat/data/IPA_blast/Ch_MDM2_nt.fasta")

#First block, seems like a complete one-to-one copy
mdm2_block1 <- mdm2_blastout[mdm2_blastout[,2] == "ctg.000255F",]
mdm2_block1 <- mdm2_block1[order(mdm2_block1[,9]),]
mdm2_block1[,"adj_q_start"] <- mdm2_block1[,7] 
#Removing overlaps to perserve frame, based on comparing breakpoints with the herring cDNA
mdm2_block1[2,"adj_q_start"] <- 60
mdm2_block1[3,"adj_q_start"] <- 137
mdm2_block1[4,"adj_q_start"] <- 267
mdm2_block1[5,"adj_q_start"] <- 326
mdm2_block1[6,"adj_q_start"] <- 386
mdm2_block1[7,"adj_q_start"] <- 443
mdm2_block1[9,"adj_q_start"] <- 932

mdm2_block1[,"adj_q_end"] <- mdm2_block1[,8]
mdm2_block1[7,"adj_q_end"] <- 600

mdm2_block1[,"start_shift"] <- mdm2_block1[,"adj_q_start"] - mdm2_block1[,7]
mdm2_block1[,"adj_start"] <- mdm2_block1[,"start_shift"] + mdm2_block1[,9]
mdm2_block1[,"end_shift"] <-  mdm2_block1[,8] - mdm2_block1[,"adj_q_end"]
mdm2_block1[,"adj_end"] <- mdm2_block1[,10] - mdm2_block1[,"end_shift"]
mdm2_block1[7,"adj_end"] <-  481435 #Adjusting 



mdm2_block1_GR <- GRanges(seqnames = mdm2_block1[,2], ranges = IRanges(start = mdm2_block1[,"adj_start"], end = mdm2_block1[,"adj_end"]))
#tp53_block1_GR <- tp53_block1_GR[order(tp53_block1_GR@ranges@start)]
Sprat_mdm2_CDS <- as.DNAbin(unlist(IPA_pri[mdm2_block1_GR]))
Sprat_est_mdm2 <- paste0(as.character(trans(Sprat_mdm2_CDS[[1]])), collapse = "")
Ch_mdm2_CDS <- as.DNAbin(CluHar_mdm[[1]])
Ch_mdm2 <- paste0(as.character(trans(Ch_mdm2_CDS[[1]])), collapse = "")

mdm2_cDNAs <- c(Sprat_mdm2_CDS, Ch_mdm2_CDS)
names(mdm2_cDNAs) <- c("Sprat", "Herring")
mdm2_DNA_clustalo <- clustalomega(mdm2_cDNAs)
image(mdm2_DNA_clustalo)

mdm2_AA <- c(trans(Sprat_mdm2_CDS[1]), trans(Ch_mdm2_CDS[1]))
names(mdm2_AA) <- c("Sprat", "Herring")
mdm2_AA_clustalo <- clustalomega(mdm2_AA)
image(mdm2_AA_clustalo)

write(">Sprattus_sprattus_putative_MDM2", file = "~/Projects/Sprat/doc/protein_comp/Sprat_MDM2.fa")
write(Sprat_est_mdm2, file = "~/Projects/Sprat/doc/protein_comp/Sprat_MDM2.fa", append = T)
alview(mdm2_AA_clustalo, file = "/Users/mapet205/Projects/Sprat/doc/protein_comp/MDM2.aln")
