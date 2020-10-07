import sys

file = open(sys.argv[1],'r')
list=[]
for line in file:
    elems = line.split("\t")
    IDs=elems[4].split(",")
    list.extend(IDs)


file2 = open(sys.argv[2],'r')
with open(sys.argv[3],'w') as f:
    for line in file2:
        elems = line.split("\t")
        foo = ([pos for pos, char in enumerate(elems[0]) if char == ":"])
        ID = elems[0][foo[0]+1:foo[1]]
        if ID in list:
            f.write(line)
