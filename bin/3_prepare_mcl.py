#!/usr/bin/env python3
### (C) 2018 - Dimitri Wirjowerdojo #######
### https://github.com/dmtr13/cluster/ ####
import sys, os, time, argparse, math, csv
import numpy as np
import pandas as pd
from joblib import Parallel, delayed

"""
Takes a distance matrix, apply thresholding, and reshape into a 3-column format
acceptable for MCL. Multiprocessing enabled.
"""

### Defining arguments for input and thresholding cut-off
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input',
                    required=True,
                    help="Input file is output of create_matrix.py")
parser.add_argument('-t', '--threshold',
                    help="""Cut off threshold. Default: prune 0.9 lowest.
                            Input as float between 0-1.""",
                    nargs="*", type=float, default=[0.9])
args = parser.parse_args()

def thresholding(enum, vect):
    """ Takes the values of the top Y percent, everything else will be 0.
    """ ### Y <-- threshold value
    if (enum+1) % 100 == 0:
        print ("Thresholding: {}/{}".format(enum+1, c))
    from scipy.stats import rankdata as rd
    vect = np.subtract(max(vect), vect)
    ranking = np.subtract(rd(vect), np.ones(len(vect))).astype(int)
    ranked_vect = []
    for x in ranking:
        if x >= top:
            ranked_vect.extend([vect.item(x)])
        else:
            ranked_vect.extend([0])
    writelist = []
    for g in range(enum, len(ranked_vect)):
        writelist.append([genes[enum], genes[g], str(ranked_vect[g])])
    return writelist

start = time.time()
print ("Reading {}...".format(args.input))
df = pd.read_csv(args.input, sep='\t', header=0, index_col=0, memory_map=True)
genes = list(df)
ar = np.array(df)
c = len(genes)
csum = sum(range(1, c+1))

filename_extension = os.path.basename(args.input)
directory = os.path.dirname(args.input)
fn, ext = os.path.splitext(filename_extension)


for t in args.threshold:
    top = math.ceil(c*t)

    print ("Creating similarity matrix and applying thresholding at {}%..."
            .format(100*t))
    results = Parallel(n_jobs=-1)(delayed(thresholding) \
                        (enum, vect) for enum, vect in enumerate(ar)    )
    results = [j for i in results for j in i]

    print ("Writing file...")
    outfile = directory+'/'+fn+'_'+str(100*t)+"_MCL"+ext
    with open(outfile, 'w') as csvfile:
        writer = csv.writer(csvfile, delimiter='\t', lineterminator='\n',
                            quoting=csv.QUOTE_NONE, escapechar="\\")
        for z, line in enumerate(results):
            if (z+1) % 100000 == 0:
                print ("Writing: {}k/{}k".format((z+1)//1000, csum//1000))
            writer.writerow(line)
    results = None
    print ("Finished {}% thresholding.\n".format(100*t))


end = time.gmtime(time.time()-start)
print ("Created similarity matrix of {} in {}.\n"
        .format(args.input,
                time.strftime("%Hh %Mm %Ss", end))
                )
