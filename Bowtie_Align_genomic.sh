#!/bin/bash

#Help message:
usage="
Usage: ./Bowtie_Align.sh [-options] numberOfThreads BowtieBaseName Reads.fastq out unique_regions.bed transcripts.fa

numberOfThreads 		Int setting the number of threads to use with BOWTIE2
BowtieBaseName 			BaseName to use for creating BOWTIE2 Index files, or if allready created, BaseName of the existing index files
Reads.fastq			Reads or transcripts that are to be aligned to the unique regions (can be fastq.gz)
out					Base name of the output files
transcripts.fa	fasta file with reference transcripts (NMD, RI etc.)

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
	   if [[ $# -ne 8 ]]; then #check for right number of arguments
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

if [[ $# -ne 6 ]]; then #check for right number of arguments
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
referencetranscripts=$6

#echo $numberOfThreads
#echo $bowtieBaseName
#echo $reads
#echo $out
#echo $unique_regions

bowtie2 --threads $numberOfThreads -x $bowtieBaseName -U $reads -S $out
bamfile=$(echo $out | rev | cut -f 2- -d '.' | rev).bam
samtools view -q 42 -@ $numberOfThreads -S -b $out > $bamfile
#bedfile=$(echo $out | rev | cut -f 2- -d '.' | rev).bed
#echo "converting bam to bed"
#bedtools bamtobed -i $bamfile > $bedfile
#intersectbedfile=$(echo $bedfile | rev | cut -f 2- -d '.' | rev)_NMD_intersect_counts.bed
#echo "intersecting with unique regions"
#bedtools intersect -c -F 0.33 -a $unique_regions -b $bedfile > $intersectbedfile
#sortedBedfile=$(echo $intersectbedfile | rev | cut -f 2- -d '.' | rev)_sorted.bed
#sort -k 4 -r -n $intersectbedfile > $sortedBedfile
echo "sorting bamfile and indexing"
sortedbamfile=$(echo $bamfile | rev | cut -f 2- -d '.' | rev)_sorted.bam
samtools sort -@ 10 $bamfile > $sortedbamfile
samtools index -@ 10 $sortedbamfile
echo "calculating coverage"
bigwigfile=$(echo $bamfile | rev | cut -f 2- -d '.' | rev)_coverage.bigwig
bamCoverage -p max/2 -b $sortedbamfile -o $bigwigfile -of bigwig
bedgraphfile=$(echo $bamfile | rev | cut -f 2- -d '.' | rev)_coverage.bedgraph
bamCoverage -p max/2 -b $sortedbamfile -o $bedgraphfile -of bedgraph
#bamCoverage -p max/2 -b ./hs_iPScm_01_Ri_vs_NMD_sorted.bam -o hs_iPScm_01_Ri_vs_NMD_coverage.bedgraph -of bedgraph
#echo "adding relative intersection"
#coveragebedfile=$(echo $bedfile | rev | cut -f 2- -d '.' | rev)_coverage_counts.bed
#bedtools coverage -F 0.33 -a $unique_regions -b $bedfile > $coveragebedfile
#coveragebedfilesorted=$(echo $bedfile | rev | cut -f 2- -d '.' | rev)_coverage_counts_sorted.bed
#sort -k 7 -r -n $coveragebedfile > $coveragebedfilesorted
#intersectbedfilerelative=$(echo $intersectbedfile | rev | cut -f 2- -d '.' | rev)_relative.bed
#cat $intersectbedfile | awk -v OFS='\t' '{print $1,$2,$3,$4,$4/($3-$2)}' > $intersectbedfilerelative
#intersectbedfilerelativesorted=$(echo $intersectbedfilerelative | rev | cut -f 2- -d '.' | rev)_sorted.bed
#sort -n -r -k 5 $intersectbedfilerelative > $intersectbedfilerelativesorted
#randomfile=$(echo $out | rev | cut -f 2- -d '.' | rev)_NMD_random_background_regions.bed
#python ./PipeTest/Pipeline/Uniqueness_scripts/BackgroundRegions.py $unique_regions $referencetranscripts $randomfile
#randomintersectfile=$(echo $out | rev | cut -f 2- -d '.' | rev)_NMD_random_intersect_counts.bed
#bedtools intersect -c -F 0.33 -a $randomfile -b $bedfile > $randomintersectfile
#randomintersectfilerelative=$(echo $randomintersectfile | rev | cut -f 2- -d '.' | rev)_relative.bed
#cat $randomintersectfile | awk -v OFS='\t' '{print $1,$2,$3,$4,$4/($3-$2)}' > $randomintersectfilerelative
#randomintersectfilesorted=$(echo $out | rev | cut -f 2- -d '.' | rev)_NMD_random_intersect_counts_relative_sorted.bed
#sort -n -r -k 5 $randomintersectfilerelative > $randomintersectfilesorted

#echo $bamfile
#echo $bedfile
#echo $intersectbedfile
#echo $sortedBedfile
