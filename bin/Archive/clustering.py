#!/usr/bin/env python3
import sys
import numpy as np
from scipy.spatial.distance import squareform as ssds
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.cluster import hierarchy
import matplotlib.pyplot as plt

print ("Clustering {}...".format(sys.argv[1]))

def strip_first_col(fname, delimiter='\t'):
    with open(fname, 'r') as fin:
        for line in fin:
            try:
               yield line.split(delimiter, 1)[1]
            except IndexError:
               continue
arr = np.loadtxt(strip_first_col(sys.argv[1]), skiprows=1)
condensed = ssds(arr)

datatitle = sys.argv[1].split('/')[-1]
try:
    if ".tsv" in datatitle or ".csv" in datatitle:
        datatitle = datatitle.replace(".tsv","").replace(".csv","")
except:
    pass

with open(sys.argv[1], 'r') as g:
    genes = (g.readline()).strip('\n').split('\t')[1:]

Z = linkage(condensed, 'ward')
dn = dendrogram(Z, labels=genes, above_threshold_color='0.5', orientation='top')

hierarchy.set_link_color_palette(['m','c','y','k'])
# fig, axes = plt.subplots(1,2, figsize=(90,40))
# dn1 = dendrogram(Z, ax=axes[0], above_threshold_color='0.50', orientation='top',
                # labels=genes)
# dn2 = dendrogram(Z, ax=axes[1], above_threshold_color="0.50",
#                 orientation='right') #bcbddc
hierarchy.set_link_color_palette(None)
plt.suptitle(datatitle)

output = "../Figures/"+datatitle+'.'
format = 'pdf'
plt.savefig(output+format, format=format)

print ("Done!")

### LINKS
# scipy hierarchical clustering https://scipy.github.io/devdocs/generated/scipy.cluster.hierarchy.dendrogram.html#scipy.cluster.hierarchy.dendrogram
# ^^ tutorial https://joernhees.de/blog/2015/08/26/scipy-hierarchical-clustering-and-dendrogram-tutorial/
# https://docs.scipy.org/doc/scipy/reference/cluster.hierarchy.html#module-scipy.cluster.hierarchy
# https://stackoverflow.com/questions/20624137/numpy-loadtxt-skip-first-column
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.squareform.html
# https://stackoverflow.com/questions/11917779/how-to-plot-and-annotate-hierarchical-clustering-dendrograms-in-scipy-matplotlib?noredirect=1&lq=1
# https://stackoverflow.com/questions/16883412/how-do-i-get-the-subtrees-of-dendrogram-made-by-scipy-cluster-hierarchy?noredirect=1&lq=1
