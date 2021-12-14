NMD_ORFs=read.table(
  file = 'G:/Justin_Backup_29.04.2020/Justin/Uni-Frankfurt/FP+Masterarbeit/PipeTest/Pipeline/Output/run_11.02.2021-13.21.38_new_NMD/ValidORF_DNA_Sequences_length.bed',
  sep = '\t', header = FALSE)
RI_ORFs=read.table(
  file = 'G:/Justin_Backup_29.04.2020/Justin/Uni-Frankfurt/FP+Masterarbeit/PipeTest/Pipeline/Output/run_11.02.2021-15.10.13_new_RI/ValidORF_DNA_Sequences_length.bed',
  sep = '\t', header = FALSE)
NMD_ORFs_Mouse=read.table(
  file = 'G:/Justin_Backup_29.04.2020/Justin/Uni-Frankfurt/FP+Masterarbeit/PipeTest/Pipeline/Output/run_15.07.2021-12.20.56_NMD_mouse/ValidORF_DNA_Sequences_length.bed',
  sep = '\t', header = FALSE)
RI_ORFs_Mouse=read.table(
  file = 'G:/Justin_Backup_29.04.2020/Justin/Uni-Frankfurt/FP+Masterarbeit/PipeTest/Pipeline/Output/run_15.07.2021-13.02.13_RI_mouse/ValidORF_DNA_Sequences_length.bed',
  sep = '\t', header = FALSE)


NMD_Unique_Protein=read.table(
  file = 'G:/Justin_Backup_29.04.2020/Justin/Uni-Frankfurt/FP+Masterarbeit/PipeTest/Pipeline/Output/run_11.02.2021-13.21.38_new_NMD/Unique_Protein_regions_length.bed',
  sep = '\t', header = FALSE)
RI_Unique_Protein=read.table(
  file = 'G:/Justin_Backup_29.04.2020/Justin/Uni-Frankfurt/FP+Masterarbeit/PipeTest/Pipeline/Output/run_11.02.2021-15.10.13_new_RI/Unique_Protein_regions_length.bed',
  sep = '\t', header = FALSE)
NMD_Unique_Protein_Mouse=read.table(
  file = 'G:/Justin_Backup_29.04.2020/Justin/Uni-Frankfurt/FP+Masterarbeit/PipeTest/Pipeline/Output/run_15.07.2021-12.20.56_NMD_mouse/Unique_Protein_regions_length.bed',
  sep = '\t', header = FALSE)
RI_Unique_Protein_Mouse=read.table(
  file = 'G:/Justin_Backup_29.04.2020/Justin/Uni-Frankfurt/FP+Masterarbeit/PipeTest/Pipeline/Output/run_15.07.2021-13.02.13_RI_mouse/Unique_Protein_regions_length.bed',
  sep = '\t', header = FALSE)


NMD_Unique_DNA=read.table(
  file = 'G:/Justin_Backup_29.04.2020/Justin/Uni-Frankfurt/FP+Masterarbeit/PipeTest/Pipeline/Output/run_11.02.2021-13.21.38_new_NMD/Unique_DNA_regions_length.bed',
  sep = '\t', header = FALSE)
RI_Unique_DNA=read.table(
  file = 'G:/Justin_Backup_29.04.2020/Justin/Uni-Frankfurt/FP+Masterarbeit/PipeTest/Pipeline/Output/run_11.02.2021-15.10.13_new_RI/Unique_DNA_regions_length.bed',
  sep = '\t', header = FALSE)
NMD_Unique_DNA_Mouse=read.table(
  file = 'G:/Justin_Backup_29.04.2020/Justin/Uni-Frankfurt/FP+Masterarbeit/PipeTest/Pipeline/Output/run_15.07.2021-12.20.56_NMD_mouse/Unique_DNA_regions_length.bed',
  sep = '\t', header = FALSE)
RI_Unique_DNA_Mouse=read.table(
  file = 'G:/Justin_Backup_29.04.2020/Justin/Uni-Frankfurt/FP+Masterarbeit/PipeTest/Pipeline/Output/run_15.07.2021-13.02.13_RI_mouse/Unique_DNA_regions_length.bed',
  sep = '\t', header = FALSE)

ORFlengths=c("NMD_Human"=list(NMD_ORFs[,2]),"RI_Human" = list(RI_ORFs[,2]), "NMD_Mouse" = list(NMD_ORFs_Mouse[,2]), "RI_Mouse" = list(RI_ORFs_Mouse[,2]))
boxplot(ORFlengths,ylab="length (bp)",outline = FALSE, main="Length of valid ORFs")

DNAlengths=c("NMD_Human"=list(NMD_Unique_DNA[,2]),"RI_Human" = list(RI_Unique_DNA[,2]), "NMD_Mouse" = list(NMD_Unique_DNA_Mouse[,2]), "RI_Mouse" = list(RI_Unique_DNA_Mouse[,2]))
boxplot(DNAlengths,ylab="length (bp)",outline = FALSE, main="Length of unique DNA Regions")

Protlengths=c("NMD_Human"=list(NMD_Unique_Protein[,2]),"RI_Human" = list(RI_Unique_Protein[,2]), "NMD_Mouse" = list(NMD_Unique_Protein_Mouse[,2]), "RI_Mouse" = list(RI_Unique_Protein_Mouse[,2]))
boxplot(Protlengths,ylab="length (aa)",outline = FALSE, main="Length of unique Protein Regions")









