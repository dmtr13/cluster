#!/usr/bin/env python3
### (C) 2018 - Dimitri Wirjowerdojo #######
### https://github.com/dmtr13/cluster/ ####
import math, sys, time, random, tempfile, shutil, os
import numpy as np
import pandas as pd
import scipy.spatial.distance as spd
from joblib import Parallel, delayed

"""
Generates a random list of #n gene expression values across different tissues.
"""

start = time.time()

df = pd.read_csv(sys.argv[1], sep='\t', header=0, index_col=0)
header = list(df)
genes = df.index.tolist()
ar = np.array(df)
c = len(ar)

if "HPA" in str(sys.argv[1]):
    reftype = "HPA"
elif "GTEx" in str(sys.argv[1]):
    reftype = "GTEx"
else:
    namn = os.path.basename(sys.argv[1])
    reftype = namn.split('.')[0]

ns = int(sys.argv[2]) ## INT - number to sample
temp_no = []
for i in range(c):
    no = random.randint(0, c)
    # print (no)
    if len(temp_no) == ns:
        break
    if not no in temp_no and no != c:
        temp_no.append(no)

assert ns == len(temp_no), "WRONG LENGTH"

n_ar = np.array([ar[z, ] for z in temp_no]) # new array
n_gen = [genes[z] for z in temp_no]
df_R = pd.DataFrame(n_ar, columns=header, index=n_gen)

df_R.to_csv("../Data/{}_{}.tsv".format(sys.argv[2], reftype), sep='\t')
print ("DONE")
