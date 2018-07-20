#!/usr/bin/env python3
### (C) 2018 - Dimitri Wirjowerdojo #######
### https://github.com/dmtr13/cluster/ ####
import math, sys, time, os, argparse
import numpy as np
import pandas as pd

from joblib import Parallel, delayed
import scipy.spatial.distance as spd

"""
A modular script that reads in the processed dataset (via process_raw.py)
and generates a distance matrix for choosable metrics and also the three-column
format accepted by MCL. Multiprocessing enabled.
"""

### List of Functions
lof = """\t[1] Relative-Euclidean \n
\t[2] Euclidean \n
\t[3] Mass-Distance \n
\t[4] Manhattan"""

### Input arguments
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', required=True, help="Input file")
parser.add_argument('-ref', '--reftype', help="Source of data. Default = HPA.",
                    type=str, default="HPA")
parser.add_argument('-f', '--function', help="""Function to create the similarity matrix for MCL.
                    Options:{lof}""".format(lof=lof),
                    type=int, default=1)
parser.add_argument('-n', '--null', help="Removes genes with null expression", type=bool, default=0)
args = parser.parse_args()

lof = {1:"Relative-Euclidean", 2:"Euclidean", 3:"Mass-Distance", 4:"Manhattan"}
if args.null == 1:
    print ("Creating a distance matrix using {} metric (no Null)...".format(lof[args.function]))
else:
    print ("Creating a distance matrix using {} metric...".format(lof[args.function]))

### Loading the basic stuff...
start = time.time()
df = pd.read_csv(args.input, sep='\t', header=0, index_col=0)
if args.null == 1:
    df = df.T.drop([col for col, val in df.T.sum().iteritems() if val==0], axis=1).T
header = list(df)
genes = df.index.tolist()
ar = np.array(df)
c = len(genes)
filename_extention = os.path.basename(args.input)
directory = os.path.dirname(args.input)
fn, ext = os.path.splitext(filename_extention)
reftype = str(args.reftype)


def normalise(inlist):
    return [(i-min(inlist))/(max(inlist)-min(inlist)) for i in inlist]

def rel_euclidean(vect1, vect2):
    summa = 0
    for v1, v2 in zip(vect1, vect2):
        ## The error is because it's dividing by 0. E.g. v1=0, v2=0, which
        ## then makes it nan.
        distance = (v1-v2)**2
        divider = (math.fabs(v1) + math.fabs(v2))/2
        if divider == 0:
            summa += 0
        else:
            summa += distance/divider
    return math.sqrt(summa)

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

def prune_top10pc(vect):
    """ Takes the values of the top Y percent, everything else will be 0.
    """ ### Y <-- threshold value
    from scipy.stats import rankdata as rd
    ranking = [int(x-1) for x in rd(vect)]
    ranked_vect = []
    for x in ranking:
        if x >= top10:
            ranked_vect.extend([vect[x]])
        else:
            ranked_vect.extend([0])
    return ranked_vect

def normalise_by_max(inlist):
    return [i-max(inlist) for i in inlist]


def calc_matrix(enum, array):
    # print ("Processing {}/{}...".format(enum, c))
    eu = []
    for j in range(c):
        vect2 = ar[j, ]
        if args.function == 1:
            eu.append(rel_euclidean(array, vect2))
        elif args.function == 2:
            eu.append(spd.euclidean(array, vect2))
        elif args.function == 3:
            eu.append(mass_distance(array, vect2))
        elif args.function == 4:
            eu.append(spd.cityblock(array, vect2))
        else:
            print ("Unspecified function.")
            sys.exit()
    if (enum+1) % 50 == 0:
        print ("Checkpoint: {}/{}".format(enum+1, c))
    return [enum, genes[enum], eu]

results = Parallel(n_jobs=-1)(delayed(calc_matrix) \
                    (z, vect1) for z, vect1 in enumerate(ar) \
                    if z < c)

if args.null == True:
    eu_output_name = "../Data/{}_{}_noNull".format(fn, lof[args.function])
else:
    eu_output_name = "../Data/{}_{}".format(fn, lof[args.function])
eu_output = open(eu_output_name+'.tsv', 'w')
eu_output.write('\t'.join(genes[:c])+'\n')

print ("Writing output...")
for r in results:
    eu_output.write(str(r[1])+'\t'+'\t'.join([str(k) for k in r[2]])+'\n')


end = time.gmtime(time.time()-start)
if args.null == 1:
    print ("Distance matrix (no Null) of {} for {} metric completed in {}."
            .format(args.input, lof[args.function], time.strftime("%Hh %Mm %Ss", end)))
else:
    print ("Distance matrix of {} for {} metric completed in {}."
        .format(args.input, lof[args.function], time.strftime("%Hh %Mm %Ss", end)))
