#This script is used to create a bedfile containing random regions of the given transcripts with the same length distribution as the given unique 
#region file.

from Bio import SeqIO
import sys
import random

#read unique region file (argument 1) and determine lengthdistribution
file = open(sys.argv[1],'r')
lengthdistribution=[]
for line in file:
    elems = line.split("\t")
    temp=int(elems[2])-int(elems[1])
    lengthdistribution.append(temp)

#read transcript.fa file into dict
source = SeqIO.parse(sys.argv[2], "fasta")
list={}
for line in source:
    list[line.id]=line.seq
randomlist=[]

#create a list with random line numbers to be chosen from the transcript.fa
for i in range(len(lengthdistribution)):
    randomlist.append(random.choice(range(len(list.keys()))))

#Go through the transcripts.fa file again. Check if the line is in randomlist and use this id with a random start postion and add the length of
#the next entry of lengthdistribution. If the end position should be higher than the last position of the transcript use the last position of the
#transcript.
source = SeqIO.parse(sys.argv[2], "fasta")
with open(sys.argv[3], "w") as out:
    i=0
    for line in source:
        while i in randomlist:
            length=lengthdistribution.pop()
            if length < len(line.seq):
                start = random.choice(range(len(line.seq)-length))
                end = start + length
                out.write(line.id + "\t" + str(start) + "\t" + str(end) + "\n")
            else:
                out.write(line.id + "\t" + "0" + "\t" + str(len(line.seq)) + "\n")
            randomlist.remove(i)
        i+=1