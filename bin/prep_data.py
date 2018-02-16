#!/usr/bin/env python3
import math, sys, time, os
import numpy as np
import pandas as pd
import scipy.spatial.distance as spd

"""
A modular script that reads in the processed dataset (via process_raw.py)
and generates a normalised matrix for Euclidean and Manhattan metrics and also
the three-column format accepted by MCL.
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


eu_output = open("../Data/{}_Euclidean.tsv".format(reftype), 'w')
eu_output.write('\t'.join(genes[:c])+'\n')
eu_mcl = open("../Data/{}_EuclideanMCL.txt".format(reftype), 'w')


manhattan_output = open("../Data/{}_Manhattan.tsv".format(reftype), 'w')
manhattan_output.write('\t'.join(genes[:c])+'\n')
manhattan_mcl = open("../Data/{}_ManhattanMCL.txt".format(reftype), 'w')

def normalise(inlist):
    """
    Normalisation done by dividing each datapoint by the datapoint with the
    highest value in one particular set/list.
    """
    maxi = max(inlist)
    newlist = [i/maxi for i in inlist]
    return newlist

count = 1
try:
    for i in range(c):
        if count % 50 == 0:
            print ("{}/{} genes...".format(count, c))

        eu = []
        eu_output.write(str(genes[i])+'\t')

        manhattan = []
        manhattan_output.write(str(genes[i])+'\t')

        vect1 = ar[i,]
        for j in range(c):
            location = "[{}, {}]".format(i+1,j+1)
            print ("Position: {}".format(location))

            vect2 = ar[j,]
            eu.append(spd.euclidean(vect1, vect2))
            manhattan.append(spd.cityblock(vect1, vect2))

        eu = normalise(eu)
        eu_output.write('\t'.join([str(x) for x in eu])+'\n')
        for k in range(i, c):
            eu_mcl.write(genes[i]+'\t'+genes[k]+'\t'+str(eu[k])+'\n')

        manhattan = normalise(manhattan)
        manhattan_output.write('\t'.join([str(x) for x in manhattan])+'\n')
        for k in range(i, c):
            manhattan_mcl.write(genes[i]+'\t'+genes[k]+'\t'+str(manhattan[k])+'\n')

        count += 1
except:
    print ("ERROR!")
    logbook = open((str(sys.argv[1])+'.log'), 'w')
    logbook.write(location+'\n')

eu_output.close()
eu_mcl.close()

manhattan_output.close()
manhattan_mcl.close()

end = time.gmtime(time.time()-start)
print ("Completed in {}.".format(time.strftime("%Hh %Mm %Ss", end)))
