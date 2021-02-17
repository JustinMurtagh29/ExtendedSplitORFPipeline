#This script reorganizes the Unique_DNA_Regions.bed file as the ORF IDs and positions can not be used when comparing with the alignment of
#the Ribo-Seq reads to the transcripts
#The ORF-ID:Start:End is extracted from the first column and added as fourth column
import sys
uniquefile = open(sys.argv[1],'r')
uniquelist=[]
for line in uniquefile:
    elements = line.strip()
    elements = elements.split("\t")
    ids = elements[0].split(":")
    temp=[ids[0],elements[1],elements[2],ids[1],ids[2],ids[3]]
    uniquelist.append(temp)
with open(sys.argv[2],'w') as f:
    for i in uniquelist:
        f.write(i[0] + "\t" + i[1] + "\t" + i[2] + "\t" + i[3] + ":" + i[4] + ":" + i[5] + "\n")