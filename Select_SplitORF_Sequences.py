import sys
from Bio import SeqIO
import csv

transcriptfile = sys.argv[1]
transcript_fastas = SeqIO.parse(open(transcriptfile),'fasta')
transcriptsequences = {}

for fasta in transcript_fastas:
        name, sequence = fasta.id, str(fasta.seq)
        transcriptsequences[name] = sequence

#print(transcriptsequences.keys())

with open(sys.argv[2]) as tsvfile:
    reader = csv.DictReader(tsvfile, dialect='excel-tab')
    TransID_list = []
    for row in reader:
        TransID_list.append(row['geneID']+'|'+row['OrfTransID'])

f = open(sys.argv[3], "w")
for id in TransID_list:
    if id in transcriptsequences.keys():
        f.write(">"+id+"\n")
        f.write(transcriptsequences[id]+"\n")
f.close()

