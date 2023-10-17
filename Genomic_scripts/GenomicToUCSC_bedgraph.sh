in=$1
ucOut=$(echo $in | rev | cut -f 2- -d '.' | rev)_UCSC.bed
scafOut=$(echo $ucOut | rev | cut -f 2- -d '.' | rev)_no_scaffold.bed

source ./PipeTest/Pipeline/Genomic_scripts/ScaffoldUmrechnung.sh $in
./PipeTest/Pipeline/Genomic_scripts/chromToUcsc -a ./PipeTest/Pipeline/Genomic_scripts/hg38.chromAlias.tsv -i $in -o $ucOut
cat $ucOut | awk '{if(index($1,"_")>0)print substr($1,1,index($1,"_")-1),$2,$3,$4; else print $1,$2,$3,$4}' > $scafOut