Oldfile=read.csv("G:/Justin_Backup_29.04.2020/Justin/Uni-Frankfurt/FP+Masterarbeit/TestOldPipe/Output/NMD/Valid_ORF_Proteins.bed", header=FALSE, sep="\t")
Newfile=read.csv("G:/Justin_Backup_29.04.2020/Justin/Uni-Frankfurt/FP+Masterarbeit/PipeTest/Pipeline/Output/run_11.02.2021-13.21.38_new_NMD/Valid_ORF_Proteins.bed", header=FALSE, sep="\t")

write(setdiff(substr(Oldfile[,1],1,31),substr(Newfile[,1],1,31)), file="C:/Users/Justi/OneDrive/Desktop/Diff_Old_New_Pipe.txt")
