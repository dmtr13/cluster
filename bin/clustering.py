#!/usr/bin/env python3
import numpy as np
from scipy.spatial.distance import squareform as ssds
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.cluster import hierarchy
import matplotlib.pyplot as plt

def strip_first_col(fname, delimiter='\t'):
    with open(fname, 'r') as fin:
        for line in fin:
            try:
               yield line.split(delimiter, 1)[1]
            except IndexError:
               continue
arr = np.loadtxt(strip_first_col("HPA_Euclidean.tsv"), skiprows=1)

clear
condensed = ssds(arr)
condensed

Z = linkage(condensed, 'ward')
Z

hierarchy.dendrogram(Z)
dn = hierarchy.dendrogram(Z)

hierarchy.set_link_color_palette(['m','c','y','k'])
fig, axes = plt.subplots(1,2, figsize=(8,3))
dn1 = hierarchy.dendrogram(Z, ax=axes[0], above_threshold_color='y', orientation='top')
dn2 = hierarchy.dendrogram(Z, ax=axes[1], above_threshold_color="#bcbddc", orientation='right')
hierarchy.set_link_color_palette(None)
plt.show()

### LINKS
# scipy hierarchical clustering https://scipy.github.io/devdocs/generated/scipy.cluster.hierarchy.dendrogram.html#scipy.cluster.hierarchy.dendrogram
# ^^ tutorial https://joernhees.de/blog/2015/08/26/scipy-hierarchical-clustering-and-dendrogram-tutorial/
# https://docs.scipy.org/doc/scipy/reference/cluster.hierarchy.html#module-scipy.cluster.hierarchy
# https://stackoverflow.com/questions/20624137/numpy-loadtxt-skip-first-column
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.squareform.html
# https://stackoverflow.com/questions/11917779/how-to-plot-and-annotate-hierarchical-clustering-dendrograms-in-scipy-matplotlib?noredirect=1&lq=1
# https://stackoverflow.com/questions/16883412/how-do-i-get-the-subtrees-of-dendrogram-made-by-scipy-cluster-hierarchy?noredirect=1&lq=1
