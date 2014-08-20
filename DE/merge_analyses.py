#!/local/python/2.7/bin/python

# Mathieu Bahin, 04/04/14
# Last update: 16/07/14

import argparse

# Getting options back
parser = argparse.ArgumentParser()
parser.add_argument('-t', dest='threshold')
parser.add_argument('-d', dest='DESeq2_file')
parser.add_argument('-e', dest='edgeR_file')
options = parser.parse_args()

# Indexing the RLOCs index
RLOCs_file = open('/home/genouest/genouest/mbahin/Annotations/RLOCs_index.txt','r')
RLOCs_index = {}
for line in RLOCs_file:
    rloc = line.split('\t')[0]
    if not RLOCs_index.has_key(rloc):
        RLOCs_index[rloc] = {}
    RLOCs_index[rloc]['enscafgs'] = line.split('\t')[1].split(',')

# Indexing the cancer gene list
cancer_mutation_file = open('/home/genouest/genouest/mbahin/Annotations/mutation_gene_list.txt','r')

cancer_mutations = []
for line in cancer_mutation_file:
    cancer_mutations.append(line.split('\t')[0])

cancer_mutation_file.close()

# Indexing the DESeq2 results
DESeq2_file = open(options.DESeq2_file,'r')
results = {}
DESeq2_file.readline()
for line in DESeq2_file:
    xloc = line.split(',')[0].split('"')[1]
    if not results.has_key(xloc):
        results[xloc] = {}
    results[xloc]['DESeq2_score'] = line.rstrip().split(',')[-1]
    if line.split(',')[2] == 'NA':
        results[xloc]['DESeq2_fc'] = 'NA'
    elif float(line.split(',')[2]) >= 0:
        results[xloc]['DESeq2_fc'] = 'up'
    else:
        results[xloc]['DESeq2_fc'] = 'down'
DESeq2_file.close()

# Indexing the edgeR results
edgeR_file = open(options.edgeR_file,'r')
edgeR_file.readline()
for line in edgeR_file:
    xloc = line.split(',')[0].split('"')[1]
    if not results.has_key(xloc):
        results[xloc] = {}
    results[xloc]['edgeR_score'] = line.rstrip().split(',')[-1]
    if float(line.split(',')[1]) >= 0:
        results[xloc]['edgeR_fc'] = 'up'
    else:
        results[xloc]['edgeR_fc'] = 'down'
edgeR_file.close()

# Producing the scores and venn files
scores_file = open('file.scores.csv','w')
venn_file = open('file.venn.csv','w')
gene_id_file = open('file.gene_id.list','w')
scores_file.write('XLOC,DESeq2,edgeR\n')
venn_file.write('XLOC,DESeq2,edgeR\n')
threshold = float(options.threshold)

DE_genes = 0
upreg = 0
downreg = 0
ambiguous = 0
cancer_known = 0
#biotypes

for xloc in results:
    scores_file.write(xloc)
    venn_file.write(xloc)
    #if results[xloc].has_key('DESeq2_score') and (results[xloc]['DESeq2_score'] != 'NA') and (results[xloc]['DESeq2_fc'] != 'NA'):
    if results[xloc].has_key('DESeq2_score') and (results[xloc]['DESeq2_score'] != 'NA'):
        scores_file.write(','+results[xloc]['DESeq2_score'])
        venn_file.write(','+str(int(float(results[xloc]['DESeq2_score']) <= threshold)))
    else:
        scores_file.write(',NA')
        venn_file.write(',0')
    if results[xloc].has_key('edgeR_score'):
        scores_file.write(','+results[xloc]['edgeR_score']+'\n')
        venn_file.write(','+str(int(float(results[xloc]['edgeR_score']) <= threshold))+'\n')
        # Producing the DE gene list
        if (float(results[xloc]['edgeR_score']) <= threshold) and (results[xloc].has_key('DESeq2_score')) and (results[xloc]['DESeq2_score'] != 'NA') and (float(results[xloc]['DESeq2_score']) <= threshold):
            DE_genes += 1
            if RLOCs_index.has_key(xloc):
                enscafgs = ','.join((RLOCs_index[xloc]['enscafgs']))
            else:
                enscafgs = 'No_ENSCAFG'
            if (results[xloc]['DESeq2_fc'] == results[xloc]['edgeR_fc']) and (results[xloc]['DESeq2_fc'] == 'up'):
                upreg += 1
                gene_id_file.write(xloc+'\t'+str(enscafgs)+'\tUp-regulated\n')
            elif (results[xloc]['DESeq2_fc'] == results[xloc]['edgeR_fc']) and (results[xloc]['DESeq2_fc'] == 'down'):
                downreg += 1
                gene_id_file.write(xloc+'\t'+str(enscafgs)+'\tDown-regulated\n')
            else:
                ambiguous += 1
                gene_id_file.write(xloc+'\t'+str(enscafgs)+'\tambiguously regulated\n')
            """
            for enscafg in intersect[xloc]:
                if enscafg in cancer_mutations:
                    cancer_known += 1
                    break
            """
    else:
        scores_file.write(',NA\n')
        venn_file.write(',0\n')
scores_file.close()
venn_file.close()

# Printing some statistics
print 'DE gene(s) found:', DE_genes
print '\tUp-regulated gene(s):', upreg
print '\tDown-regulated gene(s):', downreg
print '\tAmbiguously regulated gene(s):', ambiguous
print 'Gene(s) known to be involved in cancer mutations: '+str(cancer_known)+' / '+str(len(cancer_mutations))