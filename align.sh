#This script calls the aligning scripts to align the Ribo-Seq reads to NMD and RI transcripts as well as the Genome and to determine their 
#overlap with the previously determined unique regions (main pipeline)
#The following are the main directories for the pipeline outputs and for the alignment outputs
nmd="Path/to/nmd/alignment/main/directory"
nmdmain="Path/to/nmd/main/directory"
ri="Path/to/ri/alignment/main/directory"
rimain="Path/to/ri/main/directory"
output="./Output"
outputBowtie="./Output/BOWTIE" #For genomic alignment
#The following are the paths to the riboseq reads for each sample
ERR3367797="./Raw-Data/Heart/ERR3367797/ERR3367797.fastq.gz"
SRR6900519="./Raw-Data/Heart/SRR6900519/SRR6900519.fastq.gz"
SRR9332878="./Raw-Data/Heart/SRR9332878/SRR9332878.fastq.gz"
SRR964946="./Raw-Data/Harringtonin_and_LTM/SRR964946_HEK293_Harringtonine/SRR964946.fastq.gz"
SRR618772="./Raw-Data/Harringtonin_and_LTM/SRR618772_HEK293_lactimidomycin_rep1/SRR618772.fastq.gz"
SRR618773="./Raw-Data/Harringtonin_and_LTM/SRR618773_HEK293_lactimidomycin_rep2/SRR618773.fastq.gz"
SRR8883211="./Raw-Data/Ribo-seq/SRR8883211_GSM3718424-RiboSeq_Ctrl1/SRR8883211.fastq.gz"
SRR8883212="./Raw-Data/Ribo-seq/SRR8883212_GSM3718425-RiboSeq_Ctrl2/SRR8883212.fastq.gz"
SRR8883213="./Raw-Data/Ribo-seq/SRR8883213_GSM3718426-RiboSeq_SH1/SRR8883213.fastq.gz"
SRR8883214="./Raw-Data/Ribo-seq/SRR8883214_GSM3718427-RiboSeq_SH2/SRR8883214.fastq.gz"
#Create a Logfile for the alignments in the output directory
exec > >(tee -i $output/AlignmentLogfile.txt)
exec 2>&1

#The following block calls the Bowtie_Align script, which creates a BOWTIE index for the NMD and RI transcripts and aligns the 
#Ribo-seq data against them before checking the overlap with the determined unique regions. For further analysis a file with random regions
#is also created and used to determine background overlap
echo "Starting alignment against transcripts"
source ./PipeTest/Pipeline/Bowtie_Align.sh -i ./Path/to/NMD_transcripts.fa 10 $nmd/NMD_Bowtie_index $ERR3367797 $nmd/Heart/ERR3367797-iPScm/hs_iPScm_01_Ri_vs_NMD $nmdmain/Unique_DNA_Regions_merged.bed ./Path/to/NMD_transcripts.fa
source ./PipeTest/Pipeline/Bowtie_Align.sh 10 $nmd/NMD_Bowtie_index $SRR6900519 $nmd/Heart/SRR6900519-RiboLace/RiboLace_vs_NMD $nmdmain/Unique_DNA_Regions_merged.bed ./Path/to/NMD_transcripts.fa
source ./PipeTest/Pipeline/Bowtie_Align.sh 10 $nmd/NMD_Bowtie_index $SRR9332878 $nmd/Heart/SRR9332878-NSC/NSC_vs_NMD $nmdmain/Unique_DNA_Regions_merged.bed ./Path/to/NMD_transcripts.fa
source ./PipeTest/Pipeline/Bowtie_Align.sh 10 $nmd/NMD_Bowtie_index $SRR964946 $nmd/Harr/Harr_vs_NMD $nmdmain/Unique_DNA_Regions_merged.bed ./Path/to/NMD_transcripts.fa
source ./PipeTest/Pipeline/Bowtie_Align.sh 10 $nmd/NMD_Bowtie_index $SRR618772 $nmd/LTM/rep1/LTM_rep1_vs_NMD $nmdmain/Unique_DNA_Regions_merged.bed ./Path/to/NMD_transcripts.fa
echo "=====...............	25%"
source ./PipeTest/Pipeline/Bowtie_Align.sh 10 $nmd/NMD_Bowtie_index $SRR618773 $nmd/LTM/rep2/LTM_rep2_vs_NMD $nmdmain/Unique_DNA_Regions_merged.bed ./Path/to/NMD_transcripts.fa
source ./PipeTest/Pipeline/Bowtie_Align.sh 10 $nmd/NMD_Bowtie_index $SRR8883211 $nmd/Adenocarcinom/ctrl1/Adeno_ctrl1_vs_NMD $nmdmain/Unique_DNA_Regions_merged.bed ./Path/to/NMD_transcripts.fa
source ./PipeTest/Pipeline/Bowtie_Align.sh 10 $nmd/NMD_Bowtie_index $SRR8883212 $nmd/Adenocarcinom/ctrl2/Adeno_ctrl2_vs_NMD $nmdmain/Unique_DNA_Regions_merged.bed ./Path/to/NMD_transcripts.fa
source ./PipeTest/Pipeline/Bowtie_Align.sh 10 $nmd/NMD_Bowtie_index $SRR8883213 $nmd/Adenocarcinom/sh1/Adeno_sh1_vs_NMD $nmdmain/Unique_DNA_Regions_merged.bed ./Path/to/NMD_transcripts.fa
source ./PipeTest/Pipeline/Bowtie_Align.sh 10 $nmd/NMD_Bowtie_index $SRR8883214 $nmd/Adenocarcinom/sh2/Adeno_sh2_vs_NMD $nmdmain/Unique_DNA_Regions_merged.bed ./Path/to/NMD_transcripts.fa
echo "==========..........	50%"
source ./PipeTest/Pipeline/Bowtie_Align.sh -i ./Path/to/RI_transcripts.fa 10 $ri/retained_intron_Bowtie_index $ERR3367797 $ri/Heart/ERR3367797-iPScm/hs_iPScm_01_Ri_vs_RI $rimain/Unique_DNA_Regions_merged.bed ./Path/to/RI_transcripts.fa
source ./PipeTest/Pipeline/Bowtie_Align.sh 10 $ri/retained_intron_Bowtie_index $SRR6900519 $ri/Heart/SRR6900519-RiboLace/RiboLace_vs_RI $rimain/Unique_DNA_Regions_merged.bed ./Path/to/RI_transcripts.fa
source ./PipeTest/Pipeline/Bowtie_Align.sh 10 $ri/retained_intron_Bowtie_index $SRR9332878 $ri/Heart/SRR9332878-NSC/NSC_vs_RI $rimain/Unique_DNA_Regions_merged.bed ./Path/to/RI_transcripts.fa
source ./PipeTest/Pipeline/Bowtie_Align.sh 10 $ri/retained_intron_Bowtie_index $SRR964946 $ri/Harr/Harr_vs_RI $rimain/Unique_DNA_Regions_merged.bed ./Path/to/RI_transcripts.fa
source ./PipeTest/Pipeline/Bowtie_Align.sh 10 $ri/retained_intron_Bowtie_index $SRR618772 $ri/LTM/rep1/LTM_rep1_vs_RI $rimain/Unique_DNA_Regions_merged.bed ./Path/to/RI_transcripts.fa
echo "===============.....	75%"
source ./PipeTest/Pipeline/Bowtie_Align.sh 10 $ri/retained_intron_Bowtie_index $SRR618773 $ri/LTM/rep2/LTM_rep2_vs_RI $rimain/Unique_DNA_Regions_merged.bed ./Path/to/RI_transcripts.fa
source ./PipeTest/Pipeline/Bowtie_Align.sh 10 $ri/retained_intron_Bowtie_index $SRR8883211 $ri/Adenocarcinom/ctrl1/Adeno_ctrl1_vs_RI $rimain/Unique_DNA_Regions_merged.bed ./Path/to/RI_transcripts.fa
source ./PipeTest/Pipeline/Bowtie_Align.sh 10 $ri/retained_intron_Bowtie_index $SRR8883212 $ri/Adenocarcinom/ctrl2/Adeno_ctrl2_vs_RI $rimain/Unique_DNA_Regions_merged.bed ./Path/to/RI_transcripts.fa
source ./PipeTest/Pipeline/Bowtie_Align.sh 10 $ri/retained_intron_Bowtie_index $SRR8883213 $ri/Adenocarcinom/sh1/Adeno_sh1_vs_RI $rimain/Unique_DNA_Regions_merged.bed ./Path/to/RI_transcripts.fa
source ./PipeTest/Pipeline/Bowtie_Align.sh 10 $ri/retained_intron_Bowtie_index $SRR8883214 $ri/Adenocarcinom/sh2/Adeno_sh2_vs_RI $rimain/Unique_DNA_Regions_merged.bed ./Path/to/RI_transcripts.fa
echo "====================	100%"

#The following block calls the Bowtie_Align_genomic and Bowtie_Align_genomic2 script, which creates a BOWTIE index for the genome and aligns 
#the Ribo-seq data against it before checking the overlap with the determined unique regions. For further analysis a file with random regions
#is also created and used to determine background overlap. As the alignment does not differ for nmd and ri, the Bowtie_Align_genomic2 script
#does not perform the alignment but uses the previously created alignment file for the overlap determination
echo "Starting Alignment againt genome"
echo "....................	0%"
source ./PipeTest/Pipeline/Bowtie_Align_genomic.sh 10 $output/GRCh38.p13_Bowtie_index $ERR3367797 $outputBowtie/Heart/ERR3367797-iPScm/hs_iPScm_01_Ri_vs_GRCh38_p13 $nmdmain/Unique_DNA_Regions_merged_genomic_with_chr.bed ./Raw-Data/ncbi-genomes-2021-01-27/GCF_000001405.39_GRCh38.p13_genomic_new_headers.fa
source ./PipeTest/Pipeline/Bowtie_Align_genomic.sh 10 $output/GRCh38.p13_Bowtie_index $SRR6900519 $outputBowtie/Heart/SRR6900519-RiboLace/RiboLace_vs_GRCh38_p13 $nmdmain/Unique_DNA_Regions_merged_genomic_with_chr.bed ./Raw-Data/ncbi-genomes-2021-01-27/GCF_000001405.39_GRCh38.p13_genomic_new_headers.fa
source ./PipeTest/Pipeline/Bowtie_Align_genomic.sh 10 $output/GRCh38.p13_Bowtie_index $SRR9332878 $outputBowtie/Heart/SRR9332878-NSC/NSC_vs_GRCh38_p13 $nmdmain/Unique_DNA_Regions_merged_genomic_with_chr.bed ./Raw-Data/ncbi-genomes-2021-01-27/GCF_000001405.39_GRCh38.p13_genomic_new_headers.fa
source ./PipeTest/Pipeline/Bowtie_Align_genomic.sh 10 $output/GRCh38.p13_Bowtie_index $SRR964946 $outputBowtie/Harr/Harr_vs_GRCh38_p13 $nmdmain/Unique_DNA_Regions_merged_genomic_with_chr.bed ./Raw-Data/ncbi-genomes-2021-01-27/GCF_000001405.39_GRCh38.p13_genomic_new_headers.fa
source ./PipeTest/Pipeline/Bowtie_Align_genomic.sh 10 $output/GRCh38.p13_Bowtie_index $SRR618772 $outputBowtie/LTM/rep1/LTM_rep1_vs_GRCh38_p13 $nmdmain/Unique_DNA_Regions_merged_genomic_with_chr.bed ./Raw-Data/ncbi-genomes-2021-01-27/GCF_000001405.39_GRCh38.p13_genomic_new_headers.fa
echo "=====...............	25%"
source ./PipeTest/Pipeline/Bowtie_Align_genomic.sh 10 $output/GRCh38.p13_Bowtie_index $SRR618773 $outputBowtie/LTM/rep2/LTM_rep2_vs_GRCh38_p13 $nmdmain/Unique_DNA_Regions_merged_genomic_with_chr.bed ./Raw-Data/ncbi-genomes-2021-01-27/GCF_000001405.39_GRCh38.p13_genomic_new_headers.fa
source ./PipeTest/Pipeline/Bowtie_Align_genomic.sh 10 $output/GRCh38.p13_Bowtie_index $SRR8883211 $outputBowtie/Adenocarcinom/ctrl1/Adeno_ctrl1_vs_GRCh38_p13 $nmdmain/Unique_DNA_Regions_merged_genomic_with_chr.bed ./Raw-Data/ncbi-genomes-2021-01-27/GCF_000001405.39_GRCh38.p13_genomic_new_headers.fa
source ./PipeTest/Pipeline/Bowtie_Align_genomic.sh 10 $output/GRCh38.p13_Bowtie_index $SRR8883212 $outputBowtie/Adenocarcinom/ctrl2/Adeno_ctrl2_vs_GRCh38_p13 $nmdmain/Unique_DNA_Regions_merged_genomic_with_chr.bed ./Raw-Data/ncbi-genomes-2021-01-27/GCF_000001405.39_GRCh38.p13_genomic_new_headers.fa
source ./PipeTest/Pipeline/Bowtie_Align_genomic.sh 10 $output/GRCh38.p13_Bowtie_index $SRR8883213 $outputBowtie/Adenocarcinom/sh1/Adeno_sh1_vs_GRCh38_p13 $nmdmain/Unique_DNA_Regions_merged_genomic_with_chr.bed ./Raw-Data/ncbi-genomes-2021-01-27/GCF_000001405.39_GRCh38.p13_genomic_new_headers.fa
source ./PipeTest/Pipeline/Bowtie_Align_genomic.sh 10 $output/GRCh38.p13_Bowtie_index $SRR8883214 $outputBowtie/Adenocarcinom/sh2/Adeno_sh2_vs_GRCh38_p13 $nmdmain/Unique_DNA_Regions_merged_genomic_with_chr.bed ./Raw-Data/ncbi-genomes-2021-01-27/GCF_000001405.39_GRCh38.p13_genomic_new_headers.fa
echo "==========..........	50%"
source ./PipeTest/Pipeline/Bowtie_Align_genomic2.sh $outputBowtie/Heart/ERR3367797-iPScm/hs_iPScm_01_Ri_vs_GRCh38_p13.bed $outputBowtie/Heart/ERR3367797-iPScm/hs_iPScm_01_Ri_vs_GRCh38_p13 $rimain/Unique_DNA_Regions_merged_genomic_with_chr.bed ./Raw-Data/ncbi-genomes-2021-01-27/GCF_000001405.39_GRCh38.p13_genomic_new_headers.fa
source ./PipeTest/Pipeline/Bowtie_Align_genomic2.sh $outputBowtie/Heart/SRR6900519-RiboLace/RiboLace_vs_GRCh38_p13.bed $outputBowtie/Heart/SRR6900519-RiboLace/RiboLace_vs_GRCh38_p13 $rimain/Unique_DNA_Regions_merged_genomic_with_chr.bed ./Raw-Data/ncbi-genomes-2021-01-27/GCF_000001405.39_GRCh38.p13_genomic_new_headers.fa
source ./PipeTest/Pipeline/Bowtie_Align_genomic2.sh $outputBowtie/Heart/SRR9332878-NSC/NSC_vs_GRCh38_p13.bed $outputBowtie/Heart/SRR9332878-NSC/NSC_vs_GRCh38_p13 $rimain/Unique_DNA_Regions_merged_genomic_with_chr.bed ./Raw-Data/ncbi-genomes-2021-01-27/GCF_000001405.39_GRCh38.p13_genomic_new_headers.fa
source ./PipeTest/Pipeline/Bowtie_Align_genomic2.sh $outputBowtie/Harr/Harr_vs_GRCh38_p13.bed $outputBowtie/Harr/Harr_vs_GRCh38_p13 $rimain/Unique_DNA_Regions_merged_genomic_with_chr.bed ./Raw-Data/ncbi-genomes-2021-01-27/GCF_000001405.39_GRCh38.p13_genomic_new_headers.fa
source ./PipeTest/Pipeline/Bowtie_Align_genomic2.sh $outputBowtie/LTM/rep1/LTM_rep1_vs_GRCh38_p13.bed $outputBowtie/LTM/rep1/LTM_rep1_vs_GRCh38_p13 $rimain/Unique_DNA_Regions_merged_genomic_with_chr.bed ./Raw-Data/ncbi-genomes-2021-01-27/GCF_000001405.39_GRCh38.p13_genomic_new_headers.fa
echo "===============.....	75%"
source ./PipeTest/Pipeline/Bowtie_Align_genomic2.sh $outputBowtie/LTM/rep2/LTM_rep2_vs_GRCh38_p13.bed $outputBowtie/LTM/rep2/LTM_rep2_vs_GRCh38_p13 $rimain/Unique_DNA_Regions_merged_genomic_with_chr.bed ./Raw-Data/ncbi-genomes-2021-01-27/GCF_000001405.39_GRCh38.p13_genomic_new_headers.fa
source ./PipeTest/Pipeline/Bowtie_Align_genomic2.sh $outputBowtie/Adenocarcinom/ctrl1/Adeno_ctrl1_vs_GRCh38_p13.bed $outputBowtie/Adenocarcinom/ctrl1/Adeno_ctrl1_vs_GRCh38_p13 $rimain/Unique_DNA_Regions_merged_genomic_with_chr.bed ./Raw-Data/ncbi-genomes-2021-01-27/GCF_000001405.39_GRCh38.p13_genomic_new_headers.fa
source ./PipeTest/Pipeline/Bowtie_Align_genomic2.sh $outputBowtie/Adenocarcinom/ctrl2/Adeno_ctrl2_vs_GRCh38_p13.bed $outputBowtie/Adenocarcinom/ctrl2/Adeno_ctrl2_vs_GRCh38_p13 $rimain/Unique_DNA_Regions_merged_genomic_with_chr.bed ./Raw-Data/ncbi-genomes-2021-01-27/GCF_000001405.39_GRCh38.p13_genomic_new_headers.fa
source ./PipeTest/Pipeline/Bowtie_Align_genomic2.sh $outputBowtie/Adenocarcinom/sh1/Adeno_sh1_vs_GRCh38_p13.bed $outputBowtie/Adenocarcinom/sh1/Adeno_sh1_vs_GRCh38_p13 $rimain/Unique_DNA_Regions_merged_genomic_with_chr.bed ./Raw-Data/ncbi-genomes-2021-01-27/GCF_000001405.39_GRCh38.p13_genomic_new_headers.fa
source ./PipeTest/Pipeline/Bowtie_Align_genomic2.sh $outputBowtie/Adenocarcinom/sh2/Adeno_sh2_vs_GRCh38_p13.bed $outputBowtie/Adenocarcinom/sh2/Adeno_sh2_vs_GRCh38_p13 $rimain/Unique_DNA_Regions_merged_genomic_with_chr.bed ./Raw-Data/ncbi-genomes-2021-01-27/GCF_000001405.39_GRCh38.p13_genomic_new_headers.fa
echo "====================	100%"

#This block corrects the variables for the rmd input and calls two RMD scripts, who give insight into the alignment of the Ribo-seq data
nmd=$nmd/
ri=$ri/
outputBowtie="./Output/BOWTIE/"
R -e "rmarkdown::render('RiboSeqReport.Rmd',output_file='./Output/RiboSeq_Report.html',params=list(args = c('$nmd','$ri')))"
R -e "rmarkdown::render('RiboSeqReportGenomic.Rmd',output_file='./Output/RiboSeq_Report_genomic.html',params=list(args = c('$outputBowtie','$outputBowtie')))"