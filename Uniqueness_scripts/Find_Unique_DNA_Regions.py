import sys
import json
import re
import os
from Bio import SeqIO

print(sys.argv[1]) #json
print(sys.argv[2]) #unique_regions.fa
print(sys.argv[3]) #retained_nolinebreak.fa
print(sys.argv[4]) #gencode_transcripts.fa


# Function to check if a sequence occurs exactly once in a reference file
#
# Input:
# a reference file (The fasta file that is being analyzed for split-ORFs)
# and a sequence (the sequence of a node of the AERON output)
#
# Output:
# True if the sequence occurs exactly once (exact match)
# False if it occurs more than once or not at all
def check(sequencefile, sequence):
    f = open(sequencefile)
    result = re.findall(sequence, f.read())
    if len(result) == 1:
        return True   # Sequence occurs exactly once in the file
    else:
        return False  # Sequence occurs more than once or not at all

# Function to check if a sequence does not occur in a reference file
#
# Input:
# a reference file (Fasta file containing the transcripts of all protein coding regions)
# and a sequence (the concatenated Sequence which was unique according to 'check')
#
# Output:
# True if the sequence does not occur in the reference (sequence is unique to split-ORF transcript)
# False if it occurs
def unique_region(sequencefile, sequence):
    print("Checking if \'" + sequence + "\' does not occur in protein coding")
    f = open(sequencefile)
    if re.search(sequence, f.read()) == None:
        return True
    else:
        return False


# For easier accessibility the json file is read into a list
graph = []
for line in open(sys.argv[1], "r"):
    graph.append(json.loads(line))

# Initialize:
# previous_unique --> Boolean showing if the previous sequence was unique according to 'check'
# seq --> String containing the concatenated 'unique sequences'
# count --> Int showing how many nodes have been checked to get an idea of the progress made
previous_unique = False
seq = ""
count = 0

# Compute the uniqueness of the sequences of the nodes from AERON and concatenation of succeeding unique nodes
with open(sys.argv[2], "w") as f:
    for j in graph:
        for i in j['path']['mapping']:
            count += 1
            if count % 100 == 0:
                print("PROGRESS: " + str(count) + "nodes checked")
            if previous_unique == False:
                unique = check(sys.argv[3], i['edit'][0]['sequence'])
                if unique:
                    seq += i['edit'][0]['sequence']
                previous_unique = unique
            else:
                unique = check(sys.argv[3], i['edit'][0]['sequence'])
                if unique:
                    seq += i['edit'][0]['sequence']
                else:
                    is_unique = unique_region(sys.argv[4],seq)
                    if is_unique:
                        ref = open(sys.argv[3])
                        positions = re.search(seq,ref.read())
                        f.write(">" + j['name'] + "|" + str(positions.start()) + ":" + str(positions.end()) + " in " + os.path.basename(ref.name) + "\n")
                        f.write(seq + "\n")
                    seq = ""
                previous_unique = unique
        # If the last node is also unique call unique_regions with the momentary string seq and set seq to "" afterwards
        if not "".__eq__(seq):
            is_unique = unique_region(sys.argv[4], seq)
            if is_unique:
                ref = open(sys.argv[3])
                positions = re.search(seq, ref.read())
                f.write(">" + j['name'])
                # add + "|" + str(positions.start()) + ":" + str(positions.end()) + "|in:" + os.path.basename(ref.name) + "\n")
                # for addition of start and end position in the reference data --> total position not position in line
                f.write(seq + "\n")
            seq = ""

# Create a bedfile (sys.argv[5]) for the unique regions by reading in the unique_regions.fa created above (sys.argv[2])
# and checking where the sequences occur within the corresponding transcript sequence ((sys.argv[3]))
# The bedfile has the columns: Transcript-ID    Start of the unique region      End of the unique region
queryfile = SeqIO.parse(sys.argv[2], "fasta")
out = open(sys.argv[5], "w")
for qline in queryfile:
    # because SeqIO returns an iterator that can only be parsed once it needs to be renewed in every loop
    reffile = SeqIO.parse(sys.argv[3], "fasta")
    for rline in reffile:
        result = re.search(str(qline.seq), str(rline.seq))
        if result != None:
            out.write(qline.id + "\t" + str(result.start()) + "\t" + str(result.end()) + "\n")
