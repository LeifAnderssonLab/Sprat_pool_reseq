Sprat_v_herring_aln <- function(sprat_GR, herring_GR, herring_samples = c("BM15_HastKar_Baltic_Spring", "Gavle54_Gavle_Baltic_Autumn", "NSSH34_Norway_Atlantic_Spring", "Z14_IsleofMan_Atlantic_Autumn", "HWS22_PechoraSea_BarentsSea", "HWS42_KandalakshaBay_WhiteSea_Spring", "HWS53_KandalakshaBay_WhiteSea_Summer", "HWS12_Japan_SeaofJapan", "HWS32_WhiteSea_WhiteSea", "HWS63_Balsfjord_Atlantic", "Pacific3_Vancouver_Pacific"), sprat_genome = IPA_pri, sprat_dir = "pos"){
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

require(Biostrings)
IPA_pri <- readDNAStringSet("/Users/mapet205/Projects/Sprat/data/assemblies/IPA/final.p_ctg.fasta")
chr15_8.93_Mb <- Sprat_v_herring_aln(herring_GR = GRanges("chr15:8.92e6-8.945e6"), sprat_GR = GRanges("ctg.002308F:1-2.5e4"), sprat_dir = "neg")
chr15_8.93_Mb_clustalo <- clustalomega(x = chr15_8.93_Mb)
save(chr15_8.93_Mb_clustalo, file = "~/Projects/Sprat/data/tree/chr15_8.93_Mb_clustalo.RData")



