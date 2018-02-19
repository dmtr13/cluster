#!/usr/bin/env python3
import math, sys, time
import numpy as np
import pandas as pd
import scipy.spatial.distance as spd

"""
A modular script that reads in the processed dataset (via process_raw.py)
and generates the distance matrices using Euclidean and Manhattan metrics.
"""
start = time.time()

df = pd.read_csv(sys.argv[1], sep='\t', header=0, index_col=0)
header = list(df)
genes = df.index.tolist()
ar = np.array(df)

##
c = len(genes)
##

reftype = str(sys.argv[2])
# logbook = open((str(sys.argv[1])+'.log'), 'w')

eu_output = open("../Data/{}_Euclidean.tsv".format(reftype), 'w')
eu_output.write('\t'.join(genes[:c])+'\n')

manhattan_output = open("../Data/{}_Manhattan.tsv".format(reftype), 'w')
manhattan_output.write('\t'.join(genes[:c])+'\n')

count = 1
for i in range(c):
    if count % 50 == 0:
        print ("{}/{} genes...".format(count, c-1))

    eu = []
    eu_output.write(str(genes[i])+'\t')

    manhattan = []
    manhattan_output.write(str(genes[i])+'\t')

    vect1 = ar[i,]
    for j in range(c):
        location = "[{}, {}]".format(i+1,j+1)
        print ("Position: {}".format(location))
        # logbook.write(location+'\n')
        vect2 = ar[j,]
        eu.append(spd.euclidean(vect1, vect2))
        manhattan.append(spd.cityblock(vect1, vect2))
    eu_output.write('\t'.join([str(x) for x in eu])+'\n')
    manhattan_output.write('\t'.join([str(x) for x in manhattan])+'\n')
    count += 1

eu_output.close()
manhattan_output.close()

end = time.gmtime(time.time()-start)
print ("Completed in {}.".format(time.strftime("%Hh %Mm %Ss", end)))
