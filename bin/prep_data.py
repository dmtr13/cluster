#!/usr/bin/env python3
import math, sys, time, os, tempfile, shutil
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
import scipy.spatial.distance as spd

"""
A modular script that reads in the processed dataset (via process_raw.py)
and generates a normalised matrix for Euclidean and Manhattan metrics and also
the three-column format accepted by MCL.
"""
print ("Creating distance matrices using Euclidean metric...")
df = pd.read_csv(sys.argv[1], sep='\t', header=0, index_col=0)
start = time.time()
header = list(df)
genes = df.index.tolist()
ar = np.array(df)

##
c = len(genes)
##

if "HPA" in str(sys.argv[1]):
    reftype = "HPA"
elif "GTEx" in str(sys.argv[1]):
    reftype = "GTEx"
else:
    namn = os.path.basename(sys.argv[1])
    reftype = namn.split('.')[0]

# eu_output_name = "../Data/{}_Euclidean".format(reftype)
# eu_output = open(eu_output_name+'.tsv', 'w')
# eu_output.write('\t'.join(genes[:c])+'\n')
# eu_mcl_name = "../Data/{}_EuclideanMCL".format(reftype)
# eu_mcl = open(eu_mcl_name+'.tsv', 'w')

path = tempfile.mkdtemp()
eu_out_path = os.path.join(path, 'eu_out.mmap')
eu_out = np.memmap(eu_out_path, dtype=float, shape=(c,c), mode='w+')
# manhattan_output_name = "../Data/{}_Manhattan".format(reftype)
# manhattan_output = open(manhattan_output_name+'.tsv', 'w')
# manhattan_output.write('\t'.join(genes[:c])+'\n')
# manhattan_mcl_name = "../Data/{}_ManhattanMCL".format(reftype)
# manhattan_mcl = open(manhattan_mcl_name+'.tsv', 'w')

def normalise(inlist):
    return [(i-min(inlist))/(max(inlist)-min(inlist)) for i in inlist]

def calc_matrix(enum, array):
    if (enum+1) % 100 == 0:
        print ("CHECKPOINT {}/{}...".format(enum+1, c))
    else:
        print ("Processing {}/{}...".format(enum+1, c))
    eu = []
    for j in range(c):
        vect2 = ar[j, ]
        eu.append(spd.euclidean(array, vect2))
    eu_out[enum, ] = eu
    return True

# count = 1
# try:
#     for i in range(c):
#         eu = []
#         eu_output.write(str(genes[i])+'\t')
#
#         manhattan = []
#         manhattan_output.write(str(genes[i])+'\t')
#
#         vect1 = ar[i,]
#         for j in range(c):
#             location = "[{}, {}]".format(i+1,j+1)
#             print ("Position: {}".format(location))
#             vect2 = ar[j,]
#             eu.append(spd.euclidean(vect1, vect2))
#             # manhattan.append(spd.cityblock(vect1, vect2))
#
#         print ("Normalising Euclidean")
#         eu = normalise(eu)
#         eu_output.write('\t'.join([str(x) for x in eu])+'\n')
#
#         print ("Normalising Manhattan")
#         manhattan = normalise(manhattan)
#         manhattan_output.write('\t'.join([str(x) for x in manhattan])+'\n')
#
#         print ("Preparing MCL input...")
#         for k in range(i, c):
#             if k % 3000 == 0:
#                 print ("MCL Input Set: {}/{}".format(k, c))
#             eu_mcl.write(genes[i]+'\t'+genes[k]+'\t'+str(eu[k])+'\n')
#             manhattan_mcl.write(genes[i]+'\t'+genes[k]+'\t'+str(manhattan[k])+'\n')
#
#         if count % 5000 == 0:
#             print ("{}/{} genes...".format(count, c))
#         count += 1
# except:
#     print ("\nERROR!")
#     logbook = open((str(sys.argv[1])+'.log'), 'w')
#     print ("No. of genes: {}".format(c))
#     logbook.write("No. of genes{}\n".format(c))
#     print ("Last location: {}".format(location))
#     logbook.write("Last location: {}\n".format(location))
#     print ("Last MCL: i{}, k{}".format(i,k))
#     logbook.write("Last MCL: i{}, k{}\n".format(i,k))

results = Parallel(n_jobs=-2)(delayed(calc_matrix) \
                    (z, vect1) for z, vect1 in enumerate(ar) \
                    if z < c)

df_eu = pd.DataFrame(eu_out, columns=genes[:c], index=genes[:c])
df_eu.to_csv(path_or_buf="../Data/"+reftype+"_Euclidean.tsv", sep='\t')

try:
    shutil.rmtree(path)
    print ("Temporary file deleted.")
except:
    pass

print (df_eu.iloc[:3,:3])

# eu_output.close()
# eu_mcl.close()
# manhattan_output.close()
# manhattan_mcl.close()

end = time.gmtime(time.time()-start)
print ("Distance matrix of {} for Euclidean metric completed in {}."
        .format(sys.argv[1], time.strftime("%Hh %Mm %Ss", end)))

################################################################################
#### PART2: Inverse Covariance Matrix aka Partial Correlation Matrix ###########
################################################################################
# def incov(filehandle):
#     start = time.time()
#     print ("Processing {} MCL...".format(filehandle))
#     filename = filehandle+'.tsv'
#     df = pd.read_csv(filename, sep='\t', header=0, index_col=0)
#     ## Calculating covariance
#     df = df.cov()
#     ## Taking the inverse
#     df = pd.DataFrame(np.linalg.pinv(df.values), df.columns, df.index)
#     ## Normalising
#     x_scaled = np.array([normalise(z) for z in df.values])
#     df = pd.DataFrame(x_scaled, columns=df.columns, index=df.index)
#     ## Writing into file
#     outfile = "../Data/{}_InCov".format(reftype)
#     df.to_csv(path_or_buf=outfile+'.tsv', sep='\t')
#
#     outfile = open(outfile+'.tsv', 'r')
#     outname = open("../Data/{}_InCov_MCL.tsv".format(reftype), 'w')
#     header = list()
#     for z, line in enumerate(outfile):
#         ## this or readline is apparently better
#         ## at parsing a file with loads of lines
#         if z == 0: ## Collect information on header aka genes on the columns
#             line = line.rstrip('\n').split('\t')
#             if line[0] == '' or line[0] == 'DM':
#                 del line[0] ## Removing an empty tab char or 'DM'
#                 header.extend(line) ##
#                 n_gene = len(header)
#             else:
#                 header.extend(line)
#                 n_gene = len(header)
#         else:
#             if z % 250 == 0:
#                 print ("{}/{} genes processed.".format(z, n_gene))
#             line = line.rstrip('\n').split('\t')
#             gene, vals = line[0], [float(x) for x in line[1:]]
#             assert len(header) == len(vals), "Length not equal: {}-{}" \
#                     .format(len(header), len(vals))
#             for n in range(z-1, n_gene):
#                 outname.write(gene+'\t'+header[n]+'\t'+str(vals[n])+'\n')
#     outname.close()
#
#     end = time.gmtime(time.time()-start)
#     print ("Processing inverse covariance of {} completed in {}."
#             .format(filename, time.strftime("%Hh %Mm %Ss", end)))
#     return True
#
# incov(eu_output_name)
# incov(manhattan_output_name)

print ("Done!")
