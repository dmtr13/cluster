#!/usr/bin/env python3
### (C) 2018 - Dimitri Wirjowerdojo #######
### https://github.com/dmtr13/cluster/ ####
import sys, os, time, argparse


"""
Takes a normalised matrix and converts into a 3-column format acceptable
for MCL.
"""
start = time.time()
parser.add_argument('-i', '--input', required=True, help="Input file")
parser.add_argument('-t', '--threshold', help="Cut off threshold. Default = 10\% highest.",
                    type=float, default=0.9)

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

print ("Creating similarity matrix and applying thresholding...")
filename_extention = os.path.basename(args.input)
directory = os.path.dirname(args.input)
fn, ext = os.path.splitext(filename_extention)

matrix = open(args.input, 'r')
header = list()

# filename, file_ext = os.path.splitext(sys.argv[1])
# outname = filename+"_MCLin"+file_ext
# outname = open(outname, 'w')

print ("Reshaping...")
for z, line in enumerate(matrix):
    ## this or readline is apparently better
    ## at parsing a file with loads of lines
    if z == 0: ## Collect information on header aka genes on the columns
        line = line.rstrip('\n').split('\t')
        if line[0] == '' or line[0] == 'DM':
            del line[0] ## Removing an empty tab char or 'DM'
            header.extend(line) ##
            n_gene = len(header)
        else:
            header.extend(line)
            n_gene = len(header)
    else:
        if z % 250 == 0:
            print ("{}/{} genes processed.".format(z, n_gene))
        line = line.rstrip('\n').split('\t')
        gene, vals = line[0], [float(x) for x in line[1:]]
        assert len(header) == len(vals), "Length not equal"
        for n in range(z-1, n_gene):
        # for n in range(n_gene):
            outname.write(gene+'\t'+header[n]+'\t'+str(vals[n])+'\n')
matrix.close()
outname.close()

end = time.gmtime(time.time()-start)
print ("Completed in {}.".format(time.strftime("%Hh %Mm %Ss", end)))
