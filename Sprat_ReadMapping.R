#Author: Mats Pettersson
#mats.pettersson@imbim.uu.se


#Sprat data processing
#find /proj/snic2020-2-19/private/sprat/reads/data_from_delivery04934/ -maxdepth 5 -type f -name '*.fastq.gz' | grep "R1_00" > /proj/snic2020-2-19/private/sprat/users/mats/mapping/Sprat_R1_files.txt
#find /proj/snic2020-2-19/private/sprat/reads/data_from_delivery04934/ -maxdepth 5 -type f -name '*.fastq.gz' | grep "R2_00" > /proj/snic2020-2-19/private/sprat/users/mats/mapping/Sprat_R2_files.txt

R1_files <- read.table("~/Projects/Sprat/data/mapping/Sprat_R1_files.txt", stringsAsFactors = F)
R1_files[,"ID"] <- sub(".+/([A-Za-z0-9_-]+)_R[12]_001.fastq.gz", "\\1", R1_files[,1])
R1_files[,"Short_ID"] <- sub("_S[0-9]+_L[0-9]+", "", R1_files[,"ID"])
R1_files[,"Short_ID"] <- sub("UA[-]2823[-]", "", R1_files[,"Short_ID"])
R2_files <- read.table("~/Projects/Sprat/data/mapping/Sprat_R2_files.txt", stringsAsFactors = F)
R2_files[,"ID"] <- sub(".+/([A-Za-z0-9_-]+)_R[12]_001.fastq.gz", "\\1", R2_files[,1])
all(R2_files[,"ID"] == R1_files[,"ID"])
#in order
write(R1_files[,"Short_ID"], file = "~/Projects/Sprat/data/mapping/Sprat_sample_IDs.txt")


#GATK
#find /proj/snic2020-2-19/private/sprat/alignments/ -maxdepth 5 -type f -name '*.MarkDup.bam' > /proj/snic2020-2-19/private/sprat/users/mats/mapping/Sprat_bams.txt
sprat_bams <- read.table("~/Projects/Sprat/data/mapping/Sprat_bams.txt", stringsAsFactors = F)
sprat_bams[,"ID"] <- sub(".+/([A-Za-z0-9_-]+).sort.MarkDup.bam", "\\1", sprat_bams[,1])
write(sprat_bams[,"ID"], file = "~/Projects/Sprat/data/mapping/Sprat_bam_IDs.txt")

## 2023-05-02 Update: recreating BAMs using this refrence:/proj/snic2020-2-19/private/sprat/assemblies/current_snapshot/Sprat_DeDup_v2_HiC.fasta
#Job 16 failed (P11-UV) - disk quota exceeded
# Re-launched after removing ~1 TB

#Preparing for SRA 
#Constructing symlinks
# while read line; do ln -s "$line" "${line##*/}" ; done <../../mapping/Sprat_R1_files.txt 

#File names,no path
R1_files$file_name <- sub(".+[/](UA-2823.+)", "\\1", R1_files$V1)

