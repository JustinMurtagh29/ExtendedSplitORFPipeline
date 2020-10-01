#!/bin/bash

#Help message:
usage="
Usage: ./run_pipeline.sh [-h] proteins.fa transcripts.fa annotation.bed proteincodingtranscripts.fa

proteins.fa 			should be a multi fasta file containing the amino acid sequences of the proteins that are used as reference (whole transcriptome).
transcripts.fa 			should be a multi fasta file containing the DNA sequences of the reads/transcripts that shall be analyzed.
annotation.bed 			should be a bedfile containing the annotations for the used genome build.
				The standard annotation files for human and mouse (ENSEMBL 95) can be found in the annotations directory in the SplitORF directory.
proteinCodingTranscripts.fa 	should be a multi fasta file containing the DNA sequences of the protein coding transcripts that are used as reference.

where:
-h	show this help"

#available options for the programm
while getopts ':h' option; do
  case "$option" in
    h) echo "$usage"
       exit
       ;;
   \?) printf "illegal option: -%s\n" "$OPTARG" >&2
       echo "$usage" >&2
       exit 1
       ;;
  esac
done

#Check that all arguments are provided and are in the correct file format. Give error messages according to errors in the call and exit the programm if something goes wrong
RED='\033[0;31m' #Red colour for error messages
NC='\033[0m'	 #No colour for normal messages

#check for right number of arguments
if [ "$#" -ne 4 ]; then #check for right number of arguments
  echo -e "${RED}
ERROR while executing the Pipeline!
Wrong number of arguments.${NC}"
  echo "$usage" >&2
  exit 1
#check if any argument is a directory
elif  [ -d "$1" ] || [ -d "$2" ] || [ -d "$3" ] || [ -d "$4" ]; then
  echo -e "${RED}
ERROR while executing the Pipeline!
One or more of the arguments are directories.${NC}"
  echo "$usage" >&2
  exit 1
#check if every argument has the correct file extension
elif ! [[ $1 =~ \.fa$ ]] || ! [[ $2 =~ \.fa$ ]] || ! [[ $3 =~ \.bed$ ]] || ! [[ $4 =~ \.fa$ ]]; then
  echo -e "${RED}
ERROR while executing the Pipeline!
One or more of the arguments are not in the specified file format.${NC}"
  echo "$usage" >&2
  exit 1
fi

#Assign the given arguments
proteins=$1
transcripts=$2
annotation=$3
proteinCodingTranscripts=$4

[ -d "./Output" ] && echo "Directory ./Output exists." || echo "Directory ./Output does not exist, creating ..." && mkdir ./Output && echo "Directory ./Output created."

# create run specific output folder
timestamp=$(date "+%d.%m.%Y-%H.%M.%S")
if [ -d "./Output/run_$timestamp" ]
then 
	echo "Directory ./Output/run_$timestamp already exists." >&2
	exit 1
else
	mkdir ./Output/run_$timestamp
fi

#activate conda
source $(conda info --base)/etc/profile.d/conda.sh

#activate py2 environment in order to execute runSplitOrfs.sh (see runSplitOrfs.sh for more information)
conda activate py2
echo 'Run SplitORFs'
source ./SplitOrfs-master/runSplitOrfs.sh ./Output/run_$timestamp $proteins $transcripts $annotation

#Create a report with basic statistics of the the pipeline results
R -e "rmarkdown::render('Split-ORF_Report.Rmd',output_file='./Output/run_$timestamp/Split-ORF_Report.html',params=list(args = c('/Output/run_$timestamp/ValidProteinORFPairs.txt','/Output/run_$timestamp/UniqueProteinORFPairs_annotated.txt')))"

#activate py3_7 environment in order to execute Select_SplitORF_Sequences.py and Find_Unique_Regions.py --> Dependency BioSeqIO which needs python3
conda activate py3_7

#Select_SplitORF_Sequences.py takes the output of runSplitOrfs.sh and extracts the annotated sequences from the given transcripts.fa
echo 'Select SplitORF Sequences'
python ./Select_SplitORF_Sequences.py $transcripts ./Output/run_$timestamp/UniqueProteinORFPairs_annotated.txt ./Output/run_$timestamp/UniqueProteinORFPairsSequences.fa 

#Determine Unique DNA regions by calling mummer maxmatch with a minimum length of 20, annotating the matches (non unique regions) in a bedfile and using bedtools subtract to get the non matching regions
#which are then annotated as the unique regions in another bedfile
echo "Align ORF-transcripts(DNA) to protein coding transcripts with mummer -maxmatch"
mummer -maxmatch $proteinCodingTranscripts ./Output/run_$timestamp/UniqueProteinORFPairsSequences.fa > ./Output/run_$timestamp/DNA_maxmatch.mums
echo "Select the non matching regions as unique regions and save as bedfile"
python ./Uniqueness_scripts/Find_Unique_Regions.py ./Output/run_$timestamp/DNA_maxmatch.mums ./Output/run_$timestamp/DNA_non_unique.bed ./Output/run_$timestamp/UniqueProteinORFPairsSequences.fa ./Output/run_$timestamp/UniqueProteinORFPairsSequences.bed ./Output/run_$timestamp/Unique_DNA_Regions.bed

#Determine Unique Protein regions by calling mummer maxmatch with a minimum length of 10, annotating the matches (non unique regions) in a bedfile and using bedtools subtract to get the non matching regions
#which are then annotated as the unique regions in another bedfile
echo "Align ORF-transcripts(Protein) to protein coding transcripts mummer -maxmatch -l 10"
mummer -maxmatch -l 10 $proteins ./Output/run_$timestamp/ORFProteins.fa > ./Output/run_$timestamp/Proteins_maxmatch_l10.mums
echo "Select the non matching regions as unique regions and save as bedfile"
python ./Uniqueness_scripts/Find_Unique_Regions.py ./Output/run_$timestamp/Proteins_maxmatch_l10.mums ./Output/run_$timestamp/Protein_non_unique.bed ./Output/run_$timestamp/ORFProteins.fa ./Output/run_$timestamp/OrfProteins.bed ./Output/run_$timestamp/Unique_Protein_Regions.bed

#Use bedtools getfasta to extract the fasta sequences of the unique regions annotated in the produced bedfiles
bedtools getfasta -fi ./Output/run_$timestamp/UniqueProteinORFPairsSequences.fa -fo ./Output/run_$timestamp/Unique_DNA_regions.fa -bed ./Output/run_$timestamp/Unique_DNA_Regions.bed
bedtools getfasta -fi ./Output/run_$timestamp/ORFProteins.fa -fo ./Output/run_$timestamp/Unique_Protein_regions.fa -bed ./Output/run_$timestamp/Unique_Protein_Regions.bed

#Create a report with basics statistics of the uniqueness scripts
R -e "rmarkdown::render('Extended_Pipeline.Rmd',output_file='./Output/run_$timestamp/Uniqueness_Report.html',params=list(args = c('/Output/run_$timestamp/Unique_DNA_Regions.fa','/Output/run_$timestamp/Unique_Protein_Regions.fa')))"