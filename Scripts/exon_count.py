#!/local/python/2.7/bin/python

# Mathieu Bahin, 11/09/14

# Script to get the exon counts for 2 features (usually involved in a fusion).

import argparse, os, sys, subprocess, glob
#from subprocess import Popen, PIPE

# Setting the environment
GTF = '/home/genouest/umr6061/recomgen/dog/data/canFam3/annotation/MasterAnnotation/canfam3_cons_annot.gtf'

# Getting options back
parser = argparse.ArgumentParser()
parser.add_argument('-d', dest='dir')
parser.add_argument('-f', dest='feat1', required=True)
parser.add_argument('-g', dest='feat2', required=True)
parser.add_argument('-s', dest='samples', required=True)
parser.add_argument('-p', dest='param', default='afn')
options = parser.parse_args()

# Getting the RLOCs and ENSCAFGs indexes from module 'get_indexes'
sys.path.insert(1, '/home/genouest/genouest/mbahin/Scripts')
import classics, get_indexes
RLOCs = get_indexes.RLOCs
ENSCAFGs = get_indexes.ENSCAFGs
sys.path.remove('/home/genouest/genouest/mbahin/Scripts')

# Creating a directory for the job
if options.dir == None:
    dir_name = options.feat1+'-'+options.feat2+'_exon_count'

classics.create_wd(dir_name)

# Getting the RLOC corresponding to the 2 requested features from the RLOCs/ENSCAFGs indexes.
rloc1 = classics.get_RLOC_from_feat(options.feat1)
rloc2 = classics.get_RLOC_from_feat(options.feat2)

# Creating the small GTF file corresponding to the two features of interest
print 'Creating the reduced GTF file...'
small_GTF = rloc1+'-'+rloc2+'.gtf'
with open(small_GTF, 'w') as small_GTF_file:
    exons = []
    i1 = 1
    i2 = 1
    with open(GTF, 'r') as GTF_file:
        for line in GTF_file:
            if line.split('\t')[2] == "exon":
                if line.split('\t')[8].split('"')[1]  == rloc1:
                    exons.append('\t'.join(line.split('\t')[0:8])+'\tgene_id "'+options.feat1+'_'+str(i1)+'"; transcript_id "'+options.feat1+'_'+str(i1)+'";')
                    i1 += 1
                elif line.split('\t')[8].split('"')[1]  == rloc2 :
                    exons.append('\t'.join(line.split('\t')[0:8])+'\tgene_id "'+options.feat2+'_'+str(i2)+'"; transcript_id "'+options.feat2+'_'+str(i2)+'";')
                    i2 += 1
        exons = sorted(set(exons))
        for line in exons:
            small_GTF_file.write(line+'\n')
            
# Creating the fasta file corresponding to the exons from the GTF file
os.system('~tderrien/bin/perl/script/FASTA/gtf2fasta.pl -infile '+small_GTF+' > '+rloc1+'-'+rloc2+'.fasta')

# Launching HTSeq-count
print 'Parallel executions of HTSeq-count...'
log_file = open('file.log','w')
log_file.write('File(s) processed:\n')
to_wait = []
for f in glob.glob(options.samples):
    log_file.write(f+'\n')
    jobID = subprocess.Popen(['qsub', '-terse', '/home/genouest/genouest/mbahin/HTSeq-count/HTSeq-count_launch.sh', '-i', f, '-g', os.getcwd()+'/'+small_GTF, '-'+options.param], stdout=subprocess.PIPE).communicate()[0].rstrip()
    #jobID = Popen(['qsub', '-terse', '/home/genouest/genouest/mbahin/HTSeq-count/HTSeq-count_launch.sh', '-i', f, '-g', os.getcwd()+'/'+small_GTF, '-'+options.param], stdout=PIPE).communicate()[0].rstrip()
    to_wait.append(jobID)
log_file.close()

# Waiting for the HTSeq-count to finish
subprocess.call('. /local/env/envpython-2.6.4.sh; export DRMAA_LIBRARY_PATH=/usr/local/sge/lib/lx24-amd64/libdrmaa.so; ~/Scripts/waiter.py '+' '.join(to_wait), stdout=subprocess.PIPE, shell=True)
#subprocess.call('. /local/env/envpython-2.6.4.sh; which python; export DRMAA_LIBRARY_PATH=/usr/local/sge/lib/lx24-amd64/libdrmaa.so; ~/Scripts/waiter.py '+' '.join(to_wait), stdout=PIPE, shell=True)

# Join count files
print 'Merging and formatting the results...'
#os.chdir('PDGFB-COL3A1_GIGA-HYS-2012_T04_BEMD_exon_count')
filenames = glob.glob('*_count.dir/*count')
#print filenames
with open('merged.count', 'w') as counts_file:
    file2_file = open(filenames[1],'r')
    with open(filenames[0]) as file1_file:
        for line in file1_file:
            counts_file.write(line.rstrip()+'\t'+file2_file.readline().split('\t')[1])
    file2_file.close()

# Cleaning directory
for f in glob.glob('HTSeq-count_launch.sh.[eo]*'):
    os.remove(f)