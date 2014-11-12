#!/local/python/2.7/bin/python

# Mathieu Bahin, 04/04/14
# Last update: 15/09/14

import argparse, shutil, sys

# Getting options back
parser = argparse.ArgumentParser()
parser.add_argument('-t', dest='threshold')
parser.add_argument('-d', dest='DESeq2_file')
parser.add_argument('-e', dest='edgeR_file')
options = parser.parse_args()

# Getting the functions from 'classics.py' and the RLOCs and ENSCAFGs indexes
sys.path.insert(1, '/home/genouest/genouest/mbahin/Scripts')
import get_indexes
RLOCs = get_indexes.RLOCs
ENSCAFGs = get_indexes.ENSCAFGs
sys.path.remove('/home/genouest/genouest/mbahin/Scripts')

# Indexing the cancer gene list
cancer_mutations = []
with open('/home/genouest/genouest/mbahin/Annotations/mutation_gene_list.txt','r') as cancer_mutation_file:
    for line in cancer_mutation_file:
        cancer_mutations.append(line.split('\t')[0])

# Indexing the DESeq2 results
results = {}
with open(options.DESeq2_file,'r') as DESeq2_file:
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

# Indexing the edgeR results
with open(options.edgeR_file,'r') as edgeR_file:
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
protein_coding = 0
lncRNA = 0
other = 0
cancer_known_nb = 0

for xloc in results:
    scores_file.write(xloc)
    venn_file.write(xloc)
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

            # Determining if the gene is in the candidates list
            if RLOCs.has_key(xloc):
                cancer_known = False
                for enscafg in RLOCs[xloc]['enscafgs']:
                    if enscafg in cancer_mutations:
                        cancer_known = True
                        cancer_known_nb += 1
                        break
                if cancer_known:
                    cancer_known = 'In_cancer_mutation_list'
                else:
                    cancer_known = 'Not_listed'

            # Writing the output file 'file.gene_id.list'
            gene_id_file.write(xloc+'\t'+','.join(RLOCs[xloc]['enscafgs'])+'\t')
            if RLOCs[xloc]['enscafgs'] == ['NA']:
                gene_id_file.write(RLOCs[xloc]['orthologous'])
            else:
                names = []
                for enscafg in RLOCs[xloc]['enscafgs']:
                    if (names == []) or (ENSCAFGs[enscafg]['consensus_names'] != 'NA'):
                        names.append('|'.join(ENSCAFGs[enscafg]['consensus_names']))
                gene_id_file.write('\t'+','.join(names))
            gene_id_file.write('\t'+','.join(RLOCs[xloc]['biotypes']))

            # Updating biotype statistics
            if RLOCs[xloc]['biotypes'] == ['protein_coding']:
                protein_coding += 1
            elif RLOCs[xloc]['biotypes'] == ['lncRNA']:
                lncRNA += 1
            else:
                other += 1

            # Determining the regulation and cancer list belonging
            if (results[xloc]['DESeq2_fc'] == results[xloc]['edgeR_fc']) and (results[xloc]['DESeq2_fc'] == 'up'):
                upreg += 1
                gene_id_file.write('\tUp-regulated\t'+cancer_known+'\n')
            elif (results[xloc]['DESeq2_fc'] == results[xloc]['edgeR_fc']) and (results[xloc]['DESeq2_fc'] == 'down'):
                downreg += 1
                gene_id_file.write('\tDown-regulated\t'+cancer_known+'\n')
            else:
                ambiguous += 1
                gene_id_file.write('\tambiguously regulated\t'+cancer_known+'\n')
    else:
        scores_file.write(',NA\n')
        venn_file.write(',0\n')
scores_file.close()
venn_file.close()

# Logging some statistics and printing it
if DE_genes == 0:
    print 'No DE gene found for this analysis.'
else:
    with open('file.statistics.txt','w') as stat_file:
        stat_file.write('DE gene(s) found: '+str(DE_genes)+'\n')
        stat_file.write('\tRegulation\n')
        stat_file.write('\t\tUp-regulated gene(s): '+str(upreg)+'\n')
        stat_file.write('\t\tDown-regulated gene(s): '+str(downreg)+'\n')
        stat_file.write('\t\tAmbiguously regulated gene(s): '+str(ambiguous)+'\n')
        stat_file.write('\tBiotypes\n')
        stat_file.write('\t\tProtein coding: '+str(protein_coding)+'\n')
        stat_file.write('\t\tLncRNA: '+str(lncRNA)+'\n')
        stat_file.write('\t\tAmbiguous / Other: '+str(other)+'\n')
        stat_file.write('\tGene(s) known to be involved in cancer mutations: '+str(cancer_known_nb)+' / '+str(len(cancer_mutations))+'\n')

	# Copying stats to the standard output
    with open ('file.statistics.txt','r') as stat_file:
        shutil.copyfileobj(stat_file, sys.stdout)