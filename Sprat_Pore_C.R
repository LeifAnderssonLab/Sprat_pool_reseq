#Author: Mats Pettersson
#mats.pettersson@imbim.uu.se

eval_data_set <- read.table("~/Projects/Sprat/data/Pore_C/reads_with_sec_map_positions.txt", stringsAsFactors = F, header = F, sep = "\t")
length(unique(eval_data_set[,1]))

eval_read_set <- read.table("~/Projects/Sprat/data/Pore_C/reads_with_sec_map.txt", stringsAsFactors = F, header = F, sep = "\t")
length(unique(eval_read_set[,1]))

IPA_pri[order(IPA_pri@ranges@width, decreasing = T)]

eval_data_set_ctg.000000F <- eval_data_set[eval_data_set[,3] == "ctg.000000F", ]
length(unique(eval_data_set_ctg.000000F[,1]))
ctg.000000F_dup_reads <- names(table(eval_data_set_ctg.000000F[,1])[table(eval_data_set_ctg.000000F[,1]) > 1])
plot(x = 1, y = 1, xlim = c(0, max(eval_data_set_ctg.000000F[,4])), ylim = c(0,3), type = "n")
for(read_id in ctg.000000F_dup_reads){
  lines(x = eval_data_set_ctg.000000F[eval_data_set_ctg.000000F[,1] == read_id, 4], y = rep(x = 3*which(ctg.000000F_dup_reads == read_id)/length(ctg.000000F_dup_reads), times = sum(eval_data_set_ctg.000000F[,1] == read_id)))
}

plot(x = 1, y = 1, type = "n", xlim =c(3.5e6,4e6), ylim =c(3.5e6,4e6) ) # xlim = c(0, max(eval_data_set_ctg.000000F[,4])), ylim = c(0, max(eval_data_set_ctg.000000F[,4]))
abline(a = 0, b = 1)
for(read_id in ctg.000000F_dup_reads){
  tmp_points <- eval_data_set_ctg.000000F[eval_data_set_ctg.000000F[,1] == read_id, 4]
  for(i in 1:length(tmp_points)){
    if(i < length(tmp_points)) points(x = tmp_points[i], y = tmp_points[i+1], pch = 16, cex = 0.4, col = "red")
    if(i == length(tmp_points)) points(x = tmp_points[i], y = tmp_points[1], pch = 16, cex = 0.4, col = "red" )
  }
}

IPA_major_ctgs <- names(IPA_pri)[IPA_pri@ranges@width > 2e6]
poreC_diffs <- numeric()
for(m_ctg in IPA_major_ctgs){
  eval_data_set_ctg <- eval_data_set[eval_data_set[,3] == m_ctg, ]
  ctg_dup_reads <- names(table(eval_data_set_ctg[,1])[table(eval_data_set_ctg[,1]) > 1])
  for(read_id in ctg_dup_reads){
    tmp_map_filter <- eval_data_set_ctg[,1] == read_id & eval_data_set_ctg[,4] > 2e5 & eval_data_set_ctg[,4] < IPA_pri[m_ctg]@ranges@width - 2e5
    tmp_points <- eval_data_set_ctg[tmp_map_filter, 4]
    #for(i in 1:length(tmp_points)){
    #  if(i < length(tmp_points)) points(x = tmp_points[i], y = tmp_points[i+1], pch = 16, cex = 0.4, col = "red")
    #  if(i == length(tmp_points)) points(x = tmp_points[i], y = tmp_points[1], pch = 16, cex = 0.4, col = "red" )
    #}
    read_diffs <- numeric()
    if(length(tmp_points) > 1){ 
      read_diffs <- abs(diff(tmp_points))
      poreC_diffs <- c(poreC_diffs, read_diffs)
    }
  }
}

hist(poreC_diffs)
hist(log10(poreC_diffs), breaks = 50)





