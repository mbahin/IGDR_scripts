#!/local/python/2.7/bin/python

# Mathieu Bahin, 25/07/14

# Script to get a file with the number of fusion for intervals among the genome.
# The inputs are an interval size (default is 100000), a karyotype, a breakpoint and an output filenames.
# The output is a png file.

import sys, re, argparse

# Getting options back
parser = argparse.ArgumentParser()
parser.add_argument('-k', dest='karyotype', required=True)
parser.add_argument('-b', dest='bkpt', required=True)
parser.add_argument('-i', dest='interval_size', default = 100000)
parser.add_argument('-o', dest='output', default = 'histo.txt')
options = parser.parse_args()

# Indexing the chromosome lengths
karyotype_file = open(options.karyotype, 'r')
chr_len = {}
for line in karyotype_file:
    chr = re.match(r'cfa(.*)',line.split('\t')[2]).group(1)
    chr_len[chr] = int(line.split('\t')[5])
karyotype_file.close()

# Initializing the interval table
interval = {}
for chr in range(1,39):
    chr = str(chr)
    interval[chr] = {}
    for i in range(0, (chr_len[chr] / options.interval_size) + 1):
        interval[chr][i] = 0
interval['X'] = {}
for i in range(0, (chr_len['X'] / options.interval_size) + 1):
        interval['X'][i] = 0

# Counting the breakpoints by interval
bkpt_file = open(options.bkpt, 'r')
for line in bkpt_file:
    interval[line.split('\t')[0]][int(line.split('\t')[1]) / options.interval_size] += 1
bkpt_file.close()

# Outputting results
output_file = open(options.output, 'w')
for chr in interval:
    for i in interval[chr]:
        line = 'cfa'+chr+'\t'+str(i * options.interval_size)+'\t'
        if ((i +1) * options.interval_size - 1) > chr_len[chr]:
            line += str(chr_len[chr])
        elif i == 0:
            line += str(options.interval_size - 1)
        else:
            line += str((i + 1) * options.interval_size - 1)
        output_file.write(line+'\t'+str(interval[chr][i])+'\n')
output_file.close()