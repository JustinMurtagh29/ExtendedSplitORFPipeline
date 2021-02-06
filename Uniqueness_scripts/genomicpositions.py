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
    temp=[elements[0],elements[1],elements[2]]
    uniquelist.append(temp)
with open(sys.argv[3],'w') as f:
    for i in uniquelist:
        i[1] = int(i[1]) + int(list[i[0]][1])
        i[2] = int(i[2]) + int(list[i[0]][1])
        f.write(list[i[0]][0] + "\t" + str(i[1]) + "\t" + str(i[2]) + "\n") #i[0] + "\t" +  set at start to get transcript and gene id

