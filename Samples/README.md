# Workflow   

Step | Action | [Script](https://github.com/dmtr13/cluster/tree/master/bin) | Input | Output  
:---: | --- | :---: | --- | ---   
1 | Prepare RNAseq raw count | process_raw.py | raw count from tissues |  matrix of n_genes x n_tissues  
1.5 | Generate a small test set of M genes | gen_random.py | matrix of n_genes x n_tissues  | matrix of M_genes x n_tissues
2 | Calculate distance | prep_data.py | matrix of n/M_genes x n_tissues | Euclidean dist. (n/M_genes x n/M_genes)  
3 | Create similarity matrix by 1-normalised values | normalise_matrix.py | Euclidean dist. (n/M_genes x n/M_genes) | Similarity matrix (n/M_genes x n/M_genes)
4 | Condense matrix for MCL | prepare_mcl.py | Normalised matrix | ABC, 3-column format
5 | Run MCL | mcl | ABC, 3-column format | MCL results (out.) & .log  
