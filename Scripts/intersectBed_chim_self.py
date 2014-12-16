#!/local/python/2.7/bin/python

# Mathieu Bahin, 22/07/14

import argparse, os, subprocess

# Getting options back
parser = argparse.ArgumentParser()
parser.add_argument('-c', dest='chimCT', required=True)
options = parser.parse_args()

###### Functions

def process_intersectBed_line(line):
    ##### WARNING
    # Function to get the fasta file of the reads with one mate in the first interest region and the other in the second one.
    # The parameter indicate wether the operation is done in the sense or the reverse sense.
    #####

    if line.split('\t')[-1] != '0':
        rloc_new = line.split('\t')[14].split('"')[1]
        #tcpt = line.split('\t')[14].split('"')[3]
        tcpt_id = re.match(r'(.*)_([2-9]|[1-9][0-9]+)$',line.split('\t')[14].split('"')[3]
            if tcpt_id:
                tcpt_new = tcpt_id.group(1)
            else:
                tcpt_new = line.split('\t')[14].split('"')[3]

##### Functions end

#GTF = '/home/genouest/umr6061/recomgen/tderrien/DATA/canFam3/annotation/MasterAnnotation/BROADmRNA_lncRNA_antis.Ens75.gtfclean.mEns75.07-24-2014.gtf'
#GFF = '/home/genouest/umr6061/recomgen/dog/data/canFam3/annotation/MasterAnnotation/canfam3_cons_annot.gff'
GFF = '/home/genouest/umr6061/recomgen/dog/data/canFam3/annotation/MasterAnnotation/canfam3_cons_annot.gtf'

# Processing the input file lines
chimCT_file = open(options.chimCT,'r')
index = {}
for line in chimCT_file:
    if line.startswith('#'):
        continue

    # Getting the line information
    (chimID,chr1,end1,strand1,chr2,start2,strand2) = (line.split('\t')[1:8])
    if ('N/A' in chimID.split('---')) or (chimID.split('---')[0] == chimID.split('---')[1]):
        continue
    #print chr1,end1,strand1,chr2,start2,strand2
    if line.split('\t')[11] == '4':
        continue

    # Creating a BED file for the first part
    bkpt_file = open('bkpt.bed','w')
    bkpt_file.write('\t'.join((chr1,str(int(end1) - 1),str(int(end1) + 1),'part1','0','.'))+'\n')
    bkpt_file.close()

    # Intersecting with the reference file
    enscafg_list = []
    print 'First part'
    for line in subprocess.Popen('intersectBed -a bkpt.bed -b '+GFF+' -wao', shell=True, stdout=subprocess.PIPE).stdout:
        print line.rstrip()
        if (line.rstrip().split('\t')[-1] != '0') or (line.split('\t')[8] != 'mRNA'):
            enscafg_list.append(line.split('\t')[14].split(';')[0].split('=')[1])
            #if line.split('\t')[14].split(';')[0].split('=')[1] != 'NA':
                #enscafg_list.append(line.split('\t')[14].split('"')[5])
            #else:
                #enscafg_list.append(line.split('\t')[14].split('"')[3])
    if not index.has_key(chimID):
        index[chimID] = {}
    index[chimID]['first'] = enscafg_list

    # Creating a BED file for the second part
    bkpt_file = open('bkpt.bed','w')
    bkpt_file.write('\t'.join((chr2,str(int(start2) - 1),str(int(start2) + 1),'part2','0','.'))+'\n')
    bkpt_file.close()

    # Intersecting with the reference file
    enscafg_list = []
    print 'Second part'
    for line in subprocess.Popen('intersectBed -a bkpt.bed -b '+GFF+' -wao', shell=True, stdout=subprocess.PIPE).stdout:
        print line.rstrip()
        if line.rstrip().split('\t')[-1] != '0':
            if line.split('\t')[14].split('"')[5] != 'NA':
                enscafg_list.append(line.split('\t')[14].split('"')[5])
            else:
                enscafg_list.append(line.split('\t')[14].split('"')[3])
    index[chimID]['second'] = enscafg_list
    
chimCT_file.close()
#os.remove('bkpt.bed')

# Outputting the results
output_file = open('results.xls','w')
for chimID in index:
    output_file.write('\t'.join((chimID,','.join((index[chimID]['first'])),','.join((index[chimID]['second']))))+'\n')
output_file.close()