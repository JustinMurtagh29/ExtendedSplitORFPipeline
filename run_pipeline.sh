#!/bin/bash

#This is the main script to call the pipeline. It determines split-ORFs for the given transcript, determines their unique regions and creates files
#showing their position within a transcript and within the genome. It does so for DNA and protein sequences and extracts the relevant regions into 
#faste files. Moreover it provides general statistik reports on the findings. The output files can be further analysed by using the aligning scripts
#provided in the main directory

#Help message:
usage="
Usage: ./run_pipeline.sh [-h] proteins.fa transcripts.fa annotation.bed longestproteincodingtranscripts_protein.fa longestproteincodingtranscripts_DNA.fa GenomicPositions.bed

proteins.fa 			should be a multi fasta file containing the amino acid sequences of the proteins that are used as reference (whole transcriptome).
transcripts.fa 			should be a multi fasta file containing the DNA sequences of the reads/transcripts that shall be analyzed.
annotation.bed 			should be a bedfile containing the annotations for the used genome build.
				The standard annotation files for human and mouse (ENSEMBL 95) can be found in the annotations directory in the SplitORF directory.
longestproteincodingtranscripts_protein.fa 	should be a multi fasta file containing the protein sequences of the longest isoforms of the protein coding transcripts that are used as reference.
longestproteincodingtranscripts_DNA.fa 	should be a multi fasta file containing the DNA sequences of the longest isoforms of the protein coding transcripts that are used as reference.
Tipp: The longest Isoforms can be extracted using the getLongestIsoform.py script within the Uniqueness_scripts folder
GenomicPositions.bed		A bed file containing the chromosome specific positions of the transcripts
should be in the format: Gene stable ID	Transcript stable ID	Chromosome/scaffold name	Transcript start (bp)	Transcript end (bp)

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
if [ "$#" -ne 6 ]; then #check for right number of arguments
  echo -e "${RED}
ERROR while executing the Pipeline!
Wrong number of arguments.${NC}"
  echo "$usage" >&2
  exit 1
#check if any argument is a directory
elif  [ -d "$1" ] || [ -d "$2" ] || [ -d "$3" ] || [ -d "$4" ] || [ -d "$5" ] || [ -d "$6" ]; then
  echo -e "${RED}
ERROR while executing the Pipeline!
One or more of the arguments are directories.${NC}"
  echo "$usage" >&2
  exit 1
#check if every argument has the correct file extension
elif ! [[ $1 =~ \.fa$ ]] || ! [[ $2 =~ \.fa$ ]] || ! [[ $3 =~ \.bed$ ]] || ! [[ $4 =~ \.fa$ ]] || ! [[ $5 =~ \.fa$ ]]; then
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
proteins2=$5
genomicpositions=$6

if [ -d "./Output" ]
then
	echo "Directory ./Output already exists."
else
	echo "Directory ./Output does not exist, creating ..."
	mkdir ./Output
	echo "Directory ./Output created."
fi

# create run specific output folder
timestamp=$(date "+%d.%m.%Y-%H.%M.%S")
if [ -d "./Output/run_$timestamp" ]
then 
	echo "Directory ./Output/run_$timestamp already exists." >&2
	exit 1
else
	mkdir ./Output/run_$timestamp
fi

exec > >(tee -i ./Output/run_$timestamp/Logfile.txt)
exec 2>&1

echo "Log Location should be: [ ./Output/run_$timestamp ]"

#activate conda
source $(conda info --base)/etc/profile.d/conda.sh

#activate py2 environment in order to execute runSplitOrfs.sh (see runSplitOrfs.sh for more information)
conda activate py2
echo 'Run SplitORFs'
output="./Output/run_$timestamp"
echo "*********$output**********"


echo "run the SplitOrfs script on: " $output $proteins $transcripts $annotation $proteinCodingTranscripts $proteins2 $genomicpositions

#create Orf sequences
python ./SplitOrfs-master/OrfFinder.py $transcripts > $output/OrfProteins.fa

#activate py3_7 environment
conda activate py3_7

#Determine Unique Protein regions by calling mummer maxmatch with a minimum length of 10, annotating the matches (non unique regions) in a bedfile and using bedtools subtract to get the non matching regions
#which are then annotated as the unique regions in another bedfile
echo "Align ORF-transcripts(Protein) to protein coding transcripts mummer -maxmatch -l 10"
mummer -maxmatch -l 10 $proteins2 ./Output/run_$timestamp/ORFProteins.fa > ./Output/run_$timestamp/Proteins_maxmatch_l10.mums
echo "Select the non matching regions as unique regions and save as bedfile"
python ./Uniqueness_scripts/Find_Unique_Regions.py ./Output/run_$timestamp/Proteins_maxmatch_l10.mums ./Output/run_$timestamp/Protein_non_unique.bed ./Output/run_$timestamp/ORFProteins.fa ./Output/run_$timestamp/OrfProteins.bed ./Output/run_$timestamp/Unique_Protein_Regions.bed
python ./Uniqueness_scripts/Merge_Bedfile.py $output/Unique_Protein_Regions.bed 10 $output/Unique_Protein_Regions_merged.bed
bedtools subtract -a $output/ORFProteins.bed -b $output/Unique_Protein_Regions_merged.bed > $output/ProteinRegionsforBlast.bed
bedtools getfasta -fi $output/ORFProteins.fa -bed $output/ProteinRegionsforBlast.bed -fo $output/ProteinsforBlast.fa

conda deactivate
makeblastdb -in $proteins -out $output/ProteinDatabase -dbtype prot

#use BlastP to align the translated ORFs to the proteins (currently using 20 threads, change -num_threads otherwise
#-outfmt keyword standard results in file format:
#qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore

blastp -query $output/ProteinsforBlast.fa -db $output/ProteinDatabase -out $output/OrfsAlign.txt -evalue 10 -num_threads 8 -outfmt "6 std"

#sort the blastp output by the second column to group by proteins
sort -k2 $output/OrfsAlign.txt > $output/OrfsAlign_sorted.txt
#rm ${output}/OrfsAlign.txt

#run the detection script to parse
python ./SplitOrfs-master/DetectValidSplitOrfMatches.py $output/OrfsAlign_sorted.txt > $output/ValidProteinORFPairs.txt

#sort file per Orf-transcript ID on column 3. Here it is important to omit the head while sorting
cat $output/ValidProteinORFPairs.txt | awk 'NR<2{print ;next}{print | "sort -k3"}'  > $output/ValidProteinORFPairs_sortCol3.txt
python ./SplitOrfs-master/getLongestOrfMatches.py $output/ValidProteinORFPairs_sortCol3.txt > $output/UniqueProteinORFPairs.txt

python ./SplitOrfs-master/makeBed.py $output/UniqueProteinORFPairs.txt > $output/UniqueProteinMatches.bed

bedtools intersect -a $output/UniqueProteinMatches.bed -b $annotation -wa -F 1  -wb > $output/intersectResults.txt
python ./SplitOrfs-master/addFunctionalOverlap.py $output/UniqueProteinORFPairs.txt $output/intersectResults.txt > $output/UniqueProteinORFPairs_annotated.txt

R -e "rmarkdown::render('Split-ORF_Report.Rmd',output_file='./Output/run_$timestamp/Split-ORF_Report.html',params=list(args = c('/Output/run_$timestamp/ValidProteinORFPairs.txt','/Output/run_$timestamp/UniqueProteinORFPairs_annotated.txt')))"

#activate py3_7 environment in order to execute Select_validORF_DNA_sequences.py and Find_Unique_Regions.py --> Dependency BioSeqIO which needs python3
conda activate py3_7

#Select_validORF_DNA_sequences.py extracts the Valid ORF sequences from the given transcripts.fa
echo 'Select SplitORF DNA Sequences'
python ./Uniqueness_scripts/SelectValidOrfSequences.py ./Output/run_$timestamp/ValidProteinORFPairs.txt ./Output/run_$timestamp/OrfProteins.bed ./Output/run_$timestamp/Valid_ORF_Proteins.bed
python ./Uniqueness_scripts/Select_validORF_DNA_sequences.py $transcripts ./Output/run_$timestamp/Valid_ORF_Proteins.bed ./Output/run_$timestamp/ValidORF_DNA_Sequences.fa 

#Determine Unique DNA regions by calling mummer maxmatch with a minimum length of 20, annotating the matches (non unique regions) in a bedfile and using bedtools subtract to get the non matching regions
#which are then annotated as the unique regions in another bedfile
echo "Align ORF-transcripts(DNA) to protein coding transcripts with mummer -maxmatch"
mummer -maxmatch $proteinCodingTranscripts ./Output/run_$timestamp/ValidORF_DNA_Sequences.fa > ./Output/run_$timestamp/DNA_maxmatch.mums
echo "Select the non matching regions as unique regions and save as bedfile"
python ./Uniqueness_scripts/Find_Unique_Regions.py ./Output/run_$timestamp/DNA_maxmatch.mums ./Output/run_$timestamp/DNA_non_unique.bed ./Output/run_$timestamp/ValidORF_DNA_Sequences.fa ./Output/run_$timestamp/ValidORF_DNA_Sequences.bed ./Output/run_$timestamp/Unique_DNA_Regions.bed

#Merge the bedfile entries if the start and end positions for the same transcript only differ by the length parameter of MUMmer or less
#For more details see Merge_Bedfile.py
python ./Uniqueness_scripts/Merge_Bedfile.py ./Output/run_$timestamp/Unique_DNA_Regions.bed 20 ./Output/run_$timestamp/Unique_DNA_Regions_merged.bed

#Select the valid ORF-Proteins annotated in ValidProteinORFPairs_sortCol3.txt
python ./Uniqueness_scripts/SelectValidOrfSequences.py ./Output/run_$timestamp/ValidProteinORFPairs_sortCol3.txt ./Output/run_$timestamp/Unique_Protein_Regions_merged.bed ./Output/run_$timestamp/Unique_Protein_Regions_merged_valid.bed

#Get the genomic positions for the unique DNA and protein regions
source ./GetGenomicRegions.sh $genomicpositions ./Output/run_$timestamp/Unique_DNA_Regions_merged.bed ./Output/run_$timestamp/Unique_DNA_Regions_merged_genomic_with_chr.bed ./Output/run_$timestamp/Unique_DNA_Regions_merged_genomic_with_trID.bed ./Output/run_$timestamp/Unique_Protein_Regions_merged_valid.bed ./Output/run_$timestamp/Unique_Protein_Regions_merged_valid_genomic_with_chr.bed ./Output/run_$timestamp/Unique_Protein_Regions_merged_valid_genomic_with_trID.bed

#calculate overlap between unique DNA and protein regions
bedtools intersect -a ./Output/run_$timestamp/Unique_DNA_Regions_merged_genomic_with_trID.bed -b ./Output/run_$timestamp/Unique_Protein_Regions_merged_valid_genomic_with_trID.bed > ./Output/run_$timestamp/Unique_Regions_Overlap.bed

#Use bedtools getfasta to extract the fasta sequences of the unique regions annotated in the produced bedfiles
bedtools getfasta -fi ./Output/run_$timestamp/ValidORF_DNA_Sequences.fa -fo ./Output/run_$timestamp/Unique_DNA_regions.fa -bed ./Output/run_$timestamp/Unique_DNA_Regions_merged.bed
bedtools getfasta -fi ./Output/run_$timestamp/ORFProteins.fa -fo ./Output/run_$timestamp/Unique_Protein_regions.fa -bed ./Output/run_$timestamp/Unique_Protein_Regions_merged_valid.bed

#Reorganize Unique_DNA_Regions_merged.bed for later intersection with riboseq Alignment
python ./Uniqueness_scripts/Bedreorganize.py ./Output/run_$timestamp/Unique_DNA_Regions_merged.bed ./Output/run_$timestamp/Unique_DNA_Regions_for_comparison.bed

#Create a report with basics statistics of the uniqueness scripts
R -e "rmarkdown::render('Extended_Pipeline.Rmd',output_file='./Output/run_$timestamp/Uniqueness_Report.html',params=list(args = c('/Output/run_$timestamp/Unique_DNA_Regions.fa','/Output/run_$timestamp/Unique_Protein_Regions.fa')))"




