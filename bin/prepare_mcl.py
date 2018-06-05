#!/usr/bin/env python3
### (C) 2018 - Dimitri Wirjowerdojo #######
### https://github.com/dmtr13/cluster/ ####
import sys, os, time

print ("Preparing {} for MCL...".format(sys.argv[1]))
"""
Takes a normalised matrix and converts into a 3-column format acceptable
for MCL.
"""
start = time.time()
matrix = open(sys.argv[1], 'r')
header = list()

filename, file_ext = os.path.splitext(sys.argv[1])
outname = filename+"_MCLin"+file_ext
outname = open(outname, 'w')

m = 1
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
