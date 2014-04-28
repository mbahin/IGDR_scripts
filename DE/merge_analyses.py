# Mathieu Bahin, 04/04/14

import argparse

# Getting options back
parser = argparse.ArgumentParser()
parser.add_argument('-t', dest='threshold')
parser.add_argument('-d', dest='DESeq2_file')
parser.add_argument('-e', dest='edgeR_file')
options = parser.parse_args()

index = {}

# Indexing the DESeq2 results
DESeq2_file = open(options.DESeq2_file,'r')
DESeq2_file.readline()
for line in DESeq2_file:
    xloc = line.split(',')[0].split('"')[1]
    if not index.has_key(xloc):
        index[xloc] = {}
    index[xloc]['DESeq2'] = line.rstrip().split(',')[-1]
DESeq2_file.close()

# Indexing the edgeR results
edgeR_file = open(options.edgeR_file,'r')
edgeR_file.readline()
for line in edgeR_file:
    xloc = line.split(',')[0].split('"')[1]
    if not index.has_key(xloc):
        index[xloc] = {}
    index[xloc]['edgeR'] = line.rstrip().split(',')[-1]
edgeR_file.close()

# Producing the scores and venn files
scores_file = open('file.scores.csv','w')
venn_file = open('file.venn.csv','w')
gene_id_file = open('file.gene_id.list','w')
scores_file.write('XLOC,DESeq2,edgeR\n')
venn_file.write('XLOC,DESeq2,edgeR\n')
threshold = float(options.threshold)
for xloc in index:
    scores_file.write(xloc)
    venn_file.write(xloc)
    if index[xloc].has_key('DESeq2') and (index[xloc]['DESeq2'] != 'NA'):
        scores_file.write(','+index[xloc]['DESeq2'])
        venn_file.write(','+str(int(float(index[xloc]['DESeq2']) <= threshold)))
    else:
        scores_file.write(',NA')
        venn_file.write(',0')
    if index[xloc].has_key('edgeR'):
        scores_file.write(','+index[xloc]['edgeR']+'\n')
        venn_file.write(','+str(int(float(index[xloc]['edgeR']) <= threshold))+'\n')
        # Producing the DE gene list
        if (float(index[xloc]['edgeR']) <= threshold) and (index[xloc].has_key('DESeq2')) and (index[xloc]['DESeq2'] != 'NA') and (float(index[xloc]['DESeq2']) <= threshold):
            gene_id_file.write(xloc+'\n')
    else:
        scores_file.write(',NA\n')
        venn_file.write(',0\n')
scores_file.close()
venn_file.close()