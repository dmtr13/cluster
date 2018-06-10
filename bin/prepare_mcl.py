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
    if (enum+1) % 25 == 0:
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
    return [enum, genes[enum], ranked_vect]

def prep_write(vect):
    enum = vect[0]
    if (enum+1) % 500 == 0:
        print ("Writing: {}/{}".format(enum+1, c))
    towrite = []
    for j in range(vect[0], len(vect[2])):
        towrite.append([vect[1], genes[j], str(vect[2][j])])
    return towrite

start = time.time()
print ("Reading {}...".format(args.input))
df = pd.read_csv(args.input, sep='\t', header=0, index_col=0, memory_map=True)
genes = list(df)
ar = np.array(df)
c = len(genes)

filename_extension = os.path.basename(args.input)
directory = os.path.dirname(args.input)
fn, ext = os.path.splitext(filename_extension)


for t in args.threshold:
    top = math.ceil(c*t)

    print ("Creating similarity matrix and applying thresholding at {}%...".format(
            100*t))
    results = Parallel(n_jobs=-1)(delayed(thresholding) \
                        (enum, vect) for enum, vect in enumerate(ar)    )

    print ("Writing file...")
    outfile = directory+'/'+fn+'_'+str(100*t)+"_MCL"+ext
    re_write = Parallel(n_jobs=-1)(delayed(prep_write) \
                (vect) for vect in results)
    re_write = [j for i in re_write for j in i]


    with open(outfile, 'w') as csvfile:
        writer = csv.writer(csvfile, delimiter='\t', lineterminator='\n',
                            quoting=csv.QUOTE_NONE, escapechar="\\")
        for line in re_write:
            writer.writerow(line)

    print ("Finished {}% thresholding.\n".format(100*t))


end = time.gmtime(time.time()-start)
print ("Created similarity matrix of {} in {}.\n"
        .format(args.input,
                time.strftime("%Hh %Mm %Ss", end))
                )
