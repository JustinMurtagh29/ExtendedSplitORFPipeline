# This script fetches the genomic positions of all transcripts and adds them to the positions of the found unique regions and adds the start position
# of the corresponding ORF.
# It takes as first Input a tsv with the genomic transcript positions in the following format:
# Gene stable ID	Transcript stable ID	Chromosome/scaffold name	Transcript start (bp)	Transcript end (bp)
# The second Input file is the unique DNA regions file produced by the main pipeline
# The third argument needed for this script is the name of the output file with chromosome/scaffold name
# The fourth argument needed for this script is the name of the output file with geneID|transcriptID
import sys
file = open(sys.argv[1],'r')
list={}
for line in file:
    elems = line.split("\t")
    ids=elems[0]+"|"+elems[1]
    temp=[elems[2],elems[3],elems[4]]
    list[ids]=temp
uniquefile = open(sys.argv[2],'r')
uniquelist=[]
for line in uniquefile:
    elements = line.strip()
    elements = elements.split("\t")
    uniIDs = elements[0].split(":")
    start = int(uniIDs[2]) + int(elements[1])
    end = int(uniIDs[2]) + int(elements[2])
    temp=[uniIDs[0],start,end,uniIDs[1],uniIDs[2],uniIDs[3]]
    uniquelist.append(temp)
with open(sys.argv[3],'w') as f:
    with open(sys.argv[4],'w') as f2:
        for i in uniquelist:
            i[1] = int(i[1]) + int(list[i[0]][1])
            i[2] = int(i[2]) + int(list[i[0]][1])
            f.write(list[i[0]][0] + "\t" + str(i[1]) + "\t" + str(i[2]) + "\n")
            f2.write(i[0] + ":" + i[3] + ":" + i[4] + ":" + i[5] +"\t" + str(i[1]) + "\t" + str(i[2]) + "\n")

