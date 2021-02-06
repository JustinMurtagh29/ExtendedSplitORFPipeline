library("dplyr")
library("UpSetR")
bed <- read.table("E:/Justin_Backup_29.04.2020/Justin/Uni-Frankfurt/FP+Masterarbeit/PipeTest/Pipeline/Output/run_06.01.2021-11.01.29_NMD_with_longest_isoforms/BOWTIE/Heart/ERR3367797-iPScm/hs_iPScm_01_Ri_vs_NMD_intersect_counts_relative_sorted.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
list<- rep(NA, length(bed[,1]))
for(i in 1:length(bed[,1])){
  list[i]=paste0(bed[i,1],":",bed[i,2],":",bed[i,3])
}
rownames(bed)=list
bed = filter(bed, V4 > 5)
bed = filter(bed, V5 > 0.1)

bed2 <- read.table("E:/Justin_Backup_29.04.2020/Justin/Uni-Frankfurt/FP+Masterarbeit/PipeTest/Pipeline/Output/run_06.01.2021-11.01.29_NMD_with_longest_isoforms/BOWTIE/Heart/SRR6900519-RiboLace/RiboLace_vs_NMD_intersect_counts_relative_sorted.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
list2<- rep(NA, length(bed2[,1]))
for(i in 1:length(bed2[,1])){
  list2[i]=paste0(bed2[i,1],":",bed2[i,2],":",bed2[i,3])
}
rownames(bed2)=list2
bed2 = filter(bed2, V4 > 5)
bed2 = filter(bed2, V5 > 0.1)

bed3 <- read.table("E:/Justin_Backup_29.04.2020/Justin/Uni-Frankfurt/FP+Masterarbeit/PipeTest/Pipeline/Output/run_06.01.2021-11.01.29_NMD_with_longest_isoforms/BOWTIE/Heart/SRR9332878-NSC/NSC_vs_NMD_intersect_counts_relative_sorted.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
list3<- rep(NA, length(bed3[,1]))
for(i in 1:length(bed3[,1])){
  list3[i]=paste0(bed3[i,1],":",bed3[i,2],":",bed3[i,3])
}
rownames(bed3)=list3
bed3 = filter(bed3, V4 > 5)
bed3 = filter(bed3, V5 > 0.1)

bed4 <- read.table("E:/Justin_Backup_29.04.2020/Justin/Uni-Frankfurt/FP+Masterarbeit/PipeTest/Pipeline/Output/run_06.01.2021-11.01.29_NMD_with_longest_isoforms/BOWTIE/Adenocarcinom/ctrl1/Adeno_ctrl1_vs_NMD_intersect_counts_relative_sorted.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
list4<- rep(NA, length(bed4[,1]))
for(i in 1:length(bed4[,1])){
  list4[i]=paste0(bed4[i,1],":",bed4[i,2],":",bed4[i,3])
}
rownames(bed4)=list4
bed4 = filter(bed4, V4 > 5)
bed4 = filter(bed4, V5 > 0.1)

bed5 <- read.table("E:/Justin_Backup_29.04.2020/Justin/Uni-Frankfurt/FP+Masterarbeit/PipeTest/Pipeline/Output/run_06.01.2021-11.01.29_NMD_with_longest_isoforms/BOWTIE/Adenocarcinom/ctrl2/Adeno_ctrl2_vs_NMD_intersect_counts_relative_sorted.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
list5<- rep(NA, length(bed5[,1]))
for(i in 1:length(bed5[,1])){
  list5[i]=paste0(bed5[i,1],":",bed5[i,2],":",bed5[i,3])
}
rownames(bed5)=list5
bed5 = filter(bed5, V4 > 5)
bed5 = filter(bed5, V5 > 0.1)

bed6 <- read.table("E:/Justin_Backup_29.04.2020/Justin/Uni-Frankfurt/FP+Masterarbeit/PipeTest/Pipeline/Output/run_06.01.2021-11.01.29_NMD_with_longest_isoforms/BOWTIE/Adenocarcinom/sh1/Adeno_sh1_vs_NMD_intersect_counts_relative_sorted.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
list6<- rep(NA, length(bed6[,1]))
for(i in 1:length(bed6[,1])){
  list6[i]=paste0(bed6[i,1],":",bed6[i,2],":",bed6[i,3])
}
rownames(bed6)=list6
bed6 = filter(bed6, V4 > 5)
bed6 = filter(bed6, V5 > 0.1)

bed7 <- read.table("E:/Justin_Backup_29.04.2020/Justin/Uni-Frankfurt/FP+Masterarbeit/PipeTest/Pipeline/Output/run_06.01.2021-11.01.29_NMD_with_longest_isoforms/BOWTIE/Adenocarcinom/sh2/Adeno_sh2_vs_NMD_intersect_counts_relative_sorted.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
list7<- rep(NA, length(bed7[,1]))
for(i in 1:length(bed7[,1])){
  list7[i]=paste0(bed7[i,1],":",bed7[i,2],":",bed7[i,3])
}
rownames(bed7)=list7
bed7 = filter(bed7, V4 > 5)
bed7 = filter(bed7, V5 > 0.1)

bed8 <- read.table("E:/Justin_Backup_29.04.2020/Justin/Uni-Frankfurt/FP+Masterarbeit/PipeTest/Pipeline/Output/run_06.01.2021-11.01.29_NMD_with_longest_isoforms/BOWTIE/LTM/rep1/LTM_rep1_vs_NMD_intersect_counts_relative_sorted.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
list8<- rep(NA, length(bed8[,1]))
for(i in 1:length(bed8[,1])){
  list8[i]=paste0(bed8[i,1],":",bed8[i,2],":",bed8[i,3])
}
rownames(bed8)=list8
bed8 = filter(bed8, V4 > 5)
bed8 = filter(bed8, V5 > 0.1)

bed9 <- read.table("E:/Justin_Backup_29.04.2020/Justin/Uni-Frankfurt/FP+Masterarbeit/PipeTest/Pipeline/Output/run_06.01.2021-11.01.29_NMD_with_longest_isoforms/BOWTIE/LTM/rep2/LTM_rep2_vs_NMD_intersect_counts_relative_sorted.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
list9<- rep(NA, length(bed9[,1]))
for(i in 1:length(bed9[,1])){
  list9[i]=paste0(bed9[i,1],":",bed9[i,2],":",bed9[i,3])
}
rownames(bed9)=list9
bed9 = filter(bed9, V4 > 5)
bed9 = filter(bed9, V5 > 0.1)

bed10 <- read.table("E:/Justin_Backup_29.04.2020/Justin/Uni-Frankfurt/FP+Masterarbeit/PipeTest/Pipeline/Output/run_06.01.2021-11.01.29_NMD_with_longest_isoforms/BOWTIE/Harr/Harr_vs_NMD_intersect_counts_relative_sorted.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
list10<- rep(NA, length(bed10[,1]))
for(i in 1:length(bed10[,1])){
  list10[i]=paste0(bed10[i,1],":",bed10[i,2],":",bed10[i,3])
}
rownames(bed10)=list10
bed10 = filter(bed10, V4 > 5)
bed10 = filter(bed10, V5 > 0.1)

listInput=c('iPScm'=list(row.names(bed)),'RiboLace'=list(row.names(bed2)),'NSC'=list(row.names(bed3)),'Adeno_ctrl1'=list(row.names(bed4)),'Adeno_ctrl2'=list(row.names(bed5)),'Adeno_sh1'=list(row.names(bed6)),'Adeno_sh2'=list(row.names(bed7)),'LTM_rep1'=list(row.names(bed8)),'LTM_rep2'=list(row.names(bed9)),'Harr'=list(row.names(bed10)))

upset(fromList(listInput), order.by = "freq", nsets = 10, point.size = 2.5, line.size = 1.5, 
      mainbar.y.label = "Unique region Intersections", sets.x.label = "Matching unique regions per dataset",
      text.scale = c(1.5, 1.25, 1, 1.25, 1.25, 1.5))
