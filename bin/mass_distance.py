#!/usr/bin/env python3
import math, sys, time, os
import numpy as np
import pandas as pd
import multiprocessing as mp
from joblib import Parallel, delayed
import scipy.spatial.distance as spd

"""
A modular script that reads in the processed dataset (via process_raw.py)
and generates a normalised matrix for Euclidean and Manhattan metrics and also
the three-column format accepted by MCL. Multiprocessing enabled.
"""
start = time.time()
print ("Creating distance matrices using Euclidean metric...")

df = pd.read_csv(sys.argv[1], sep='\t', header=0, index_col=0)
header = list(df)
genes = df.index.tolist()
ar = np.array(df)
c = len(genes)
reftype = str(sys.argv[2])

def normalise(inlist):
    return [(i-min(inlist))/(max(inlist)-min(inlist)) for i in inlist]

eu_output_name = "../Data/MD/{}_{}_Euclidean_Norm".format(reftype, c)
eu_output = open(eu_output_name+'.tsv', 'w')
eu_output.write('\t'.join(genes[:c])+'\n')

eu_mcl_name = "../Data/MD/{}_{}_EuclideanMCL_Norm".format(reftype, c)
eu_mcl = open(eu_mcl_name+'.tsv', 'w')

def mass_distance(v1, v2):
    vec_len = len(v1)
    MDscore = 1
    assert vec_len == len(v2), "Unequal vector length!"

    for x in range(vec_len):
        minimum = min(v1[x], v2[x])
        maximum = max(v1[x], v2[x])
        v1v2 = np.concatenate((v1, v2))
        freq_x = len((np.where((v1v2 >= minimum) &
                                (v1v2 <= maximum))[0]))
        ## logx(A)+logx(B)=logx(A*B)
        MDscore *= freq_x
    ## MDScore = - (sum from i=1 to d log(MASS))
    MDscore = math.log10(MDscore)
    return MDscore

def calc_matrix(enum, array):
    print ("Processing {}/{}...".format(enum, c))
    eu = []
    for j in range(c):
        vect2 = ar[j, ]
        eu.append(mass_distance(array, vect2))
    if (enum+1) % 25 == 0:
        print ("Checkpoint: {}/{}".format(enum+1, c))
    eu = normalise(eu)
    eu = [1-x for x in eu]
    return [enum, genes[enum], eu]

results = Parallel(n_jobs=-1)(delayed(calc_matrix) \
                    (z, vect1) for z, vect1 in enumerate(ar) \
                    if z < c)


print ("Writing output...")
for r in results:
    eu_output.write(str(r[1])+'\t'+'\t'.join([str(k) for k in r[2]])+'\n')
    for h in range(int(r[0]), c):
        eu_mcl.write(str(r[1])+'\t'+genes[h]+'\t'+str(r[2][h])+'\n')

end = time.gmtime(time.time()-start)
print ("Mass-Distance similarity matrix of {} calculated in {}.\nDone!"
        .format(sys.argv[1], time.strftime("%Hh %Mm %Ss", end)))
