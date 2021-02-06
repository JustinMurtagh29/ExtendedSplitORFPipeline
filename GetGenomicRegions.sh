#!/bin/bash

#Help message:
usage="
Usage: ./GetGenomicRegions.sh [-h] GenomicAnnotation.bed UniqueDNA.bed UniqueDNAOut.bed UniqueProtein.bed UniqueProteinOut.bed

GenomicAnnotation.bed 			a BED file containing the genomic positions for all used transcripts
UniqueDNA.bed 				the unique_DNA_regions_merged.bed file from the main pipeline
UniqueDNAOut.bed 			the name for the unique DNA regions file with genomic positions
UniqueProtein.bed			the unique_protein_regions_merged_valid.bed file from the main pipeline
UniqueProteinOut.bed 			the name for the unique protein regions file with genomic positions

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
if [ "$#" -ne 5 ]; then #check for right number of arguments
  echo -e "${RED}
ERROR while executing the Pipeline!
Wrong number of arguments.${NC}"
  echo "$usage" >&2
  exit 1
#check if any argument is a directory
elif  [ -d "$1" ] || [ -d "$2" ] || [ -d "$3" ] || [ -d "$4" ] || [ -d "$5" ]; then
  echo -e "${RED}
ERROR while executing the Pipeline!
One or more of the arguments are directories.${NC}"
  echo "$usage" >&2
  exit 1
#check if every argument has the correct file extension
elif ! [[ $1 =~ \.bed$ ]] || ! [[ $2 =~ \.bed$ ]] || ! [[ $3 =~ \.bed$ ]] || ! [[ $4 =~ \.bed$ ]] || ! [[ $5 =~ \.bed$ ]]; then
  echo -e "${RED}
ERROR while executing!
One or more of the arguments are not in the specified file format.${NC}"
  echo "$usage" >&2
  exit 1
fi
GenomicAnnotation=$1
UniqueDNA=$2
UniqueDNAOut=$3
UniqueProtein=$4
UniqueProteinOut=$5
python ./Uniqueness_scripts/genomicpositions.py $GenomicAnnotation $UniqueDNA $UniqueDNAOut
python ./Uniqueness_scripts/genomicpositions_protein.py $GenomicAnnotation $UniqueProtein $UniqueProteinOut