#!/bin/bash

# Get commandline options and arguments and initialize the starting variables  
cores=all
new_Genome_Graph=false
while getopts ":n:g:h" opt; do
  case $opt in
    # Help function showing the usage and possible options
    h)
      echo -e "Usage:\n"
      echo -e "\$ run_pipeline.sh [options] Proteins.fa Transcripts.fa Annotation.bed Genome.gtf ProteinCoding.fa\n\n"
      echo -e "Options:\n"
      echo -e "-n \t\t set the number of cores to be used in the snakemake step (Default = all)"
      echo -e "-g Genome.fa \t Create a new Graph file using the provided Genome fasta"
      ;;
    # The option -n sets the number of cores for the snakemake step
    # and will only accept integers as arguments. Defaults to all. 
    n)
      re='^[0-9]+$'
	if ! [[ $OPTARG =~ $re ]] ; then
  	  echo "Error: Argument of -$opt is not a number" >&2; exit 1
	else
	  cores=$OPTARG
	fi
      #echo "Numbber of cores was set to $cores" >&2
      ;;
    # The option -g will add the creation of a new genome graph file using the provided genome fasta
    g)
      if [[ $OPTARG == *.fa || $OPTARG == *.fasta ]]
      then
      	echo "A new graph file will be created using $OPTARG"
      	Genome=$OPTARG
      	new_Genome_Graph=true
      else
      	echo "$OPTARG provided to generate graph seems to be in the wrong format."
      fi
      ;;
    :)
      echo "Error: Option -$OPTARG needs Argument" >&2; exit 1
      ;;
    \?)
      echo "Error: Invalid option: -$OPTARG. Use -h flag for help." >&2; exit 1
      ;;
  esac
done
echo "Number of cores was set to $cores" >&2
shift $((OPTIND -1))

if ! [[ $# == 5 ]]
then
	echo "Invalid set of Arguments. Use -h flag for help" >&2; exit 1
fi
for var in "$@"
do
    if [[ $var == *[_\.\;\:]*.* ]]
    then
    	echo "Please make sure that your filenames only use letters in their names" >&2; exit 1
    fi
done


Proteins=$1
Transcripts=$2
Annotation=$3
GenomeAnnotation=$4
ProteinCoding=$5

# create run specific output folder
timestamp=$(date "+%Y.%m.%d-%H.%M.%S")
if [ -d "./Output/run_$timestamp" ]
then 
	echo "Directory ./Output/run_$timestamp already exists." >&2; exit 1
else
	mkdir ./Output/run_$timestamp
fi

#activate conda
source $(conda info --base)/etc/profile.d/conda.sh

#activate py2 environment in order to execute runSplitOrfs.sh
conda activate py2
echo 'Run SplitORFs'
source ./SplitOrfs-master/runSplitOrfs.sh ./Output/run_$timestamp $Proteins $Transcripts $Annotation
conda activate py3_7
echo 'Select SplitORF Sequences'
python ./Select_SplitORF_Sequences.py $Transcripts ./Output/run_$timestamp/UniqueProteinORFPairs_annotated.txt ./Output/run_$timestamp/UniqueProteinORFPairsSequences.fa 

# copy all input files necessary for Aeron to run into the Aeron/input folder
echo 'Copy all input files necessary for Aeron to run into the Aeron/input folder'
cp ./Output/run_/UniqueProteinORFPairsSequences.fa ./Aeron/input
echo -e '#####                     (33%)\r\c'
cp $GenomeAnnotation ./Aeron/input
echo -e '#############             (66%)\r\c'
cp $ProteinCoding ./Aeron/input
echo -e '#######################   (100%)\'

# If the provided gfa for hg38 is not applicable for your pipeline set -g and a new grap file will be created
if [ $new_Genome_Graph == true ]
then
	echo 'Copying Genome.fa to Aeron/input in order to create new graph file'
	cp $Genome ./Aeron/input
	GenomeBase="$(basename -- $Genome)"
	GenomeAnnBase="$(basename -- $Genome)"
	GenomeGraph="${GenomeBase%.*}""Graph.gfa"
	echo 'Creating new Graphfile'
	python ./Aeron/AeronScripts/GraphBuilder.py -e ./Aeron/input/$GenomeBase -g ./Aeron/input/$GenomeAnnBase -o ./Aeron/input/$GenomeGraph
fi


cd Aeron
# activate the snakemake environment and run snakemake with all available cores --> set cores using "-n" if you do not want to use all available cores
# make sure /usr/bin/python directs to a python3 executable
conda activate snakemake
snakemake --cores=3 #set via option
conda deactivate
cd ..

cp ./Aeron/output/*.json ./Output/run_$timestamp
# activate the py3_7 environment to run the python scripts for  unique regions --> dependencies need a python3 version < 3.8
conda activate py3_7

# delete linebreaks in fasta file of the transcripts for String comparison
TranscriptsBase=$(basename -- $Transcripts)
TranscriptsNoline="${TranscriptsBase%.*}""Nolinebreak.fa"
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' $Transcripts > ./Output/run_$timestamp/$TranscriptNoline

python ./Uniqueness_scripts/Find_Unique_DNA_Regions.py ./Output/run_$timestamp/*.json ./Output/run_$timestamp/Unique_DNA_Regions.fa ./Output/run_$timestamp/$TranscriptsNoline $ProteinCoding ./Output/run_$timestamp/Unique_DNA_Regions.bed


















