#!/usr/bin/env python3
import math, sys, time, random, tempfile, shutil
import numpy as np
import pandas as pd
import scipy.spatial.distance as spd
from joblib import Parallel, delayed


start = time.time()

df = pd.read_csv(sys.argv[1], sep='\t', header=0, index_col=0)
header = list(df)
genes = df.index.tolist()
ar = np.array(df)
c = len(ar)

ns = int(sys.argv[2]) ## INT - number to sample
temp_no = []
for i in range(c):
    no = random.randint(0, c)
    # print (no)
    if len(temp_no) == ns:
        break
    if not no in temp_no and no != c:
        temp_no.append(no)

print (len(temp_no))

n_ar = np.array([ar[z, ] for z in temp_no]) # new array
n_gen = [genes[z] for z in temp_no]
df_R = pd.DataFrame(n_ar, columns=header, index=n_gen)

df_R.to_csv("../Data/2kR.tsv", sep='\t')
print ("DONE")
