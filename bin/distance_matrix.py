#!/bin/usr/env python3
import math, sys
import numpy as np
import pandas as pd
import scipy.spatial.distance as spd

df = pd.read_csv(sys.argv[1], sep='\t', header=0, index_col=0)
header = list(df)
genes = df.index.tolist()
ar = np.array(df)

# def euclidean(vect1, vect2):
#     return ((vect1-vect2)**2)

euclid = []
count = 1

##
c = len(genes)
# c = 100
##

reftype = str(sys.argv[2])
logbook = open((str(sys.argv[1])+'.log'), 'w')
eu_output = open("../Data/{}_Euclidean.tsv".format(reftype), 'w')
manhattan_output = open("../Data/{}_Manhattan.tsv".format(reftype), 'w')

eu_output.write('DM\t'+'\t'.join(genes[:c])+'\n')
manhattan_output.write('DM\t'+'\t'.join(genes[:c])+'\n')

for i in range(c):
    if count % 50 == 0:
        print ("{}/{} genes...".format(count, c-1))
    eu = []
    manhattan = []
    eu_output.write(str(genes[i])+'\t')
    manhattan_output.write(str(genes[i])+'\t')
    for j in range(c):
        location = "Position [{}, {}]".format(i,j)
        print (location)
        logbook.write(location+'\n')
        vect1 = ar[i,]
        vect2 = ar[j,]
        eu.append(spd.euclidean(vect1, vect2))
        manhattan.append(spd.cityblock(vect1, vect2))
        # output.write(str(math.sqrt(sum(euclidean(vect1, vect2))))+'\n')
        # eu.append(math.sqrt(sum(euclidean(vect1, vect2))))
    eu_output.write('\t'.join([str(x) for x in eu])+'\n')
    manhattan_output.write('\t'.join([str(x) for x in manhattan])+'\n')
    # euclid.append(eu)
    count += 1

print ("Done!")

# df = pd.DataFrame(euclid, index=genes, columns=genes)
# df.to_csv('euclid.tsv', index=True, header=True, sep='\t')

# print (df.iloc[0,])
# print (df.iloc[0,15:18])

# sys.exit()

# with open('ex.txt') as f:
#     print zip(*[line.split() for line in f])[1]

# filename = open(sys.argv[1], 'r')
# filename = filename.read().splitlines()
# header = filename[0].split('\t')
# del header[0]
# c = len(header)
# d = len(filename[1:])
# genes = []
# values = [[]] * c
# # print (len(values))
#
# print (filename[1])
# for z in filename[1]:
#     print (z)
#     z = z.split()
#     # genes.append(z[0])
#     # del z[0]
#     # print (z)
#
# sys.exit()

# count = 1
# for z in range(1,len(filename)):
#
#     row = filename[z].split('\t')
#     genes.append(row[0])
#     del row[0]
#     # sys.exit()
#     for a, col in enumerate(row):
#         # print (col)
#
#         values[a].append(col)
#     print (values[z])
#     count += 1
#     if count == 2:
#         break
        # values[col].append(row[col])

    # print (len(row))
    # for col in range(len(row)):
    #     if col == 0:
    #         genes.append(row[col])
    #     else:
    #         # print (col)
    #         values[col-1].append(float(row[col]))


    # row = filename[z].split('\t')
    # # genes.append(row[0])
    # for m in range(len(row)):
    #     if m == 0:
    #         # print (row[m])
    #         genes.append(row[m])
    #     else:
    #         # continue
    #         values[m-1].append(row[m])
    #         # print(m)#, row[m])
    #     # print (row[m])
    #     # values[m].extend(row[m+1])
# print ("\n\n\n\n\n")
# print (values, len(values))
# print (values[0], len(values), len(values[1]))

# sys.exit()

# for z in filename[1:]:
#     z = z.split('\t')



#
# c = len(filename[1:])
# dimension = len(header)
# genes = []
# values = []#[None] * c

# values = zip(*[line.split()[1:] for line in filename[1:]])
# for line in filename[1:]:
#     line = line.split()[1:]

# print (values)
# sys.exit()
#
# for lines in filename[1:]:
#     lines = lines.split()
#     genes.append(lines[0])
#     values.append(zip(*[float(x) for x in lines[1:]]))
# tissue_tpm = dict(zip(header, values))
# print (tissue_tpm)
# print (len(tissue_tpm))
# sys.exit()
#
# def euclidean(tissue_tpm):
#     eucvalues = []
#     for t1 in tissue_tpm:
#         print ("Calculating %s set..." % t1)
#         for t2 in tissue_tpm:
#             theset = [(a-b)**2 for a, b in
#                         zip(tissue_tpm[t1], tissue_tpm[t2])]
#             theset = math.sqrt(sum(theset))
#             eucvalues.append(theset)
#     eucvalues = np.array(eucvalues).reshape((c,c))
#     return eucvalues
#
# def manhattan(tissue_tpm):
#     eucvalues = []
#     for t1 in tissue_tpm:
#         print ("Calculating %s set..." % t1)
#         for t2 in tissue_tpm:
#             theset = [(abs(a-b)) for a, b in
#                         zip(tissue_tpm[t1], tissue_tpm[t2])]
#             theset = sum(theset)
#             eucvalues.append(theset)
#     eucvalues = np.array(eucvalues).reshape((c,c))
#     print (eucvalues)
#
# manhattan(tissue_tpm)

# if __name__ == '__main__':
#     euclidean("../Data/HPA_processed.tsv")
