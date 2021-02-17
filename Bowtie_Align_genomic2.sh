#!/bin/bash

#USE THIS FILE FOR RI UNIQUE REGIONS AND GENOMIC ALIGNMENT. Bowtie_Align_genomic NEEDS TO BE PERFORMED BEFORE THIS SCRIPT CAN BE USED

#Help message:
usage="
Usage: ./Bowtie_Align.sh [-options] bedfile out unique_regions.bed transcripts.fa

bedfile				bedfile of the alignment created in Bowtie_Align_genomic.sh
Reads.fastq			Reads or transcripts that are to be aligned to the unique regions (can be fastq.gz)
out					Base name of the output files
transcripts.fa	fasta file with reference transcripts (NMD, RI etc.)

where:
-h			show this help"

#Check that all arguments are provided and are in the correct file format. Give error messages according to errors in the call and exit the programm if something goes wrong
RED='\033[0;31m' #Red colour for error messages
NC='\033[0m'	 #No colour for normal messages


#available options for the programm
while getopts ':h' option; do
  case "$option" in
    h) echo "$usage"
       exit 1
       ;;
   \?) printf "illegal option: -%s\n" "$OPTARG" >&2
       echo "$usage" >&2
       exit 1
       ;;
  esac
  shift 2
done

if [[ $# -ne 4 ]]; then #check for right number of arguments
  echo -e "${RED}
ERROR while executing the script!
Wrong number of arguments.${NC}"
  echo "$usage" >&2
  exit 1
fi

bedfile=$1
out=$2.sam
unique_regions=$3
referencetranscripts=$4

#determine the number of reads intersecting at least 33% with a unique region and sort the file accordingly
intersectbedfile=$(echo $bedfile | rev | cut -f 2- -d '.' | rev)_RI_intersect_counts.bed
echo "intersecting with unique regions"
bedtools intersect -c -F 0.33 -a $unique_regions -b $bedfile > $intersectbedfile
sortedBedfile=$(echo $intersectbedfile | rev | cut -f 2- -d '.' | rev)_sorted.bed
sort -k 4 -r -n $intersectbedfile > $sortedBedfile

#calculate the relative read count (number of aligning reads / length of unique region) and sort according to relative read count
intersectbedfilerelative=$(echo $intersectbedfile | rev | cut -f 2- -d '.' | rev)_relative.bed
cat $intersectbedfile | awk -v OFS='\t' '{print $1,$2,$3,$4,$4/($3-$2)}' > $intersectbedfilerelative
intersectbedfilerelativesorted=$(echo $intersectbedfilerelative | rev | cut -f 2- -d '.' | rev)_sorted.bed
sort -n -r -k 5 $intersectbedfilerelative > $intersectbedfilerelativesorted

#create a file with random regions of the same length distribution as original unique regions and again determine read count and relative read count
randomfile=$(echo $out | rev | cut -f 2- -d '.' | rev)_random_background_regions.bed
python ./PipeTest/Pipeline/Uniqueness_scripts/BackgroundRegions.py $unique_regions $referencetranscripts $randomfile
randomintersectfile=$(echo $out | rev | cut -f 2- -d '.' | rev)_RI_random_intersect_counts.bed
bedtools intersect -c -F 0.33 -a $randomfile -b $bedfile > $randomintersectfile
randomintersectfilerelative=$(echo $randomintersectfile | rev | cut -f 2- -d '.' | rev)_RI_relative.bed
cat $randomintersectfile | awk -v OFS='\t' '{print $1,$2,$3,$4,$4/($3-$2)}' > $randomintersectfilerelative
randomintersectfilesorted=$(echo $out | rev | cut -f 2- -d '.' | rev)_RI_random_intersect_counts_relative_sorted.bed
sort -n -r -k 5 $randomintersectfilerelative > $randomintersectfilesorted

