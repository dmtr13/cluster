#!/usr/bin/env python3
### (C) 2018 - Dimitri Wirjowerdojo #######
### https://github.com/dmtr13/cluster/ ####
import sys, os, time, argparse, math
import numpy as np
import pandas as pd
from joblib import Parallel, delayed

"""
Takes a distance matrix, apply thresholding, and reshape into a 3-column format
acceptable for MCL. Multiprocessing enabled.
"""

### Defining arguments for input and thresholding cut-off
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', required=True, help="Input file")
parser.add_argument('-t', '--threshold', help="Cut off threshold. Default: prune 0.9 lowest. Input between 0-1.",
                    type=float, default=0.9)
args = parser.parse_args()

def thresholding(enum, vect):
    """ Takes the values of the top Y percent, everything else will be 0.
    """ ### Y <-- threshold value
    if (enum+1) % 25 == 0:
        print ("Checkpoint: {}/{}".format(enum+1, c))
    from scipy.stats import rankdata as rd
    vect = [max(vect) - el for el in vect]
    ranking = [int(x-1) for x in rd(vect)]
    ranked_vect = []
    for x in ranking:
        if x >= top:
            ranked_vect.extend([vect[x]])
        else:
            ranked_vect.extend([0])
    return [enum, genes[enum], ranked_vect]

start = time.time()
print ("Reading {}...".format(args.input))
df = pd.read_csv(args.input, sep='\t', header=0, index_col=0)
genes = list(df)
ar = np.array(df)
c = len(genes)

filename_extension = os.path.basename(args.input)
directory = os.path.dirname(args.input)
fn, ext = os.path.splitext(filename_extension)
top = math.ceil(c*args.threshold)


print ("Creating similarity matrix and applying thresholding at {}%...".format(
        100*args.threshold))
results = Parallel(n_jobs=-1)(delayed(thresholding) \
                    (enum, vect) for enum, vect in enumerate(ar)    )

print ("Writing file...")
outfile = open(directory+'/'+fn+'_'+str(100*args.threshold)+"_MCL"+ext, 'w')
# outfile.write('\t'.join(genes[:c])+'\n')
for r in results:
    for h in range(r[0], c):
        outfile.write(str(r[1])+'\t'+genes[h]+'\t'+str(r[2][h])+'\n')
outfile.close()

end = time.gmtime(time.time()-start)
print ("Created similarity matrix of {} with {}% thresholding in {}."
        .format(args.input, 100*args.threshold,
                time.strftime("%Hh %Mm %Ss", end))
                )
