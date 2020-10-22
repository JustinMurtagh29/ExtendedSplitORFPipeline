#!/bin/bash

#Help message:
usage="
Usage: ./Bowtie_Align.sh [-options] numberOfThreads BowtieBaseName Reads.fastq out unique_regions.bed

numberOfThreads 		Int setting the number of threads to use with BOWTIE2
BowtieBaseName 			BaseName to use for creating BOWTIE2 Index files, or if allready created, BaseName of the existing index files
Reads.fastq			Reads or transcripts that are to be aligned to the unique regions (can be fastq.gz)
out					Base name of the output files
unique_regions.bed	Bedfile with the annotated unique regions	

where:
-h			show this help
-i unique_regions.fa	create new index files for the provided unique_regions.fa"

#Check that all arguments are provided and are in the correct file format. Give error messages according to errors in the call and exit the programm if something goes wrong
RED='\033[0;31m' #Red colour for error messages
NC='\033[0m'	 #No colour for normal messages


#available options for the programm
while getopts ':hi' option; do
  case "$option" in
    h) echo "$usage"
       exit 1
       ;;
	i) echo "creating new index files"
	   if [[ $# -ne 7 ]]; then #check for right number of arguments
		echo -e "${RED}
ERROR while executing the script!
Wrong number of arguments.${NC}"
		echo "$usage" >&2
		exit 1
	   fi
       bowtie2-build --threads $3 $2 $4
	   #echo $3
	   #echo $2
	   #echo $4
	   ;;
   \?) printf "illegal option: -%s\n" "$OPTARG" >&2
       echo "$usage" >&2
       exit 1
       ;;
  esac
  shift 2
done

if [[ $# -ne 5 ]]; then #check for right number of arguments
  echo -e "${RED}
ERROR while executing the script!
Wrong number of arguments.${NC}"
  echo "$usage" >&2
  exit 1
fi

numberOfThreads=$1
bowtieBaseName=$2
reads=$3
out=$4.sam
unique_regions=$5

#echo $numberOfThreads
#echo $bowtieBaseName
#echo $reads
#echo $out
#echo $unique_regions

bowtie2 --threads $numberOfThreads -x $BowtieBaseName -U $reads -S $out
bamfile=$(echo $out | rev | cut -f 2- -d '.' | rev).bam
samtools view -@ $numberOfThreads -S -b $out > $bamfile
bedfile=$(echo $out | rev | cut -f 2- -d '.' | rev).bed
bedtools bamtobed -i $bamfile > $bedfile
intersectbedfile=$(echo $bedfile | rev | cut -f 2- -d '.' | rev)_intersect_counts.bed
bedtools intersect -c -F 0.33 -a $unique_regions -b $bedfile > $intersectbedfile
sortedBedfile=$(echo $intersectbedfile | rev | cut -f 2- -d '.' | rev)_sorted.bed
sort -k 4 -r -n $intersectbedfile > $sortedBedfile

#echo $bamfile
#echo $bedfile
#echo $intersectbedfile
#echo $sortedBedfile
