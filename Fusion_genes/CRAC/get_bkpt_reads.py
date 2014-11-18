#!/local/python/2.7/bin/python

# Mathieu Bahin, 12/06/14
# Last update: 15/09/14

# Script to get the split reads involved in and the paired-end reads around a breakpoint from files produced by CRAC.
# Inputs:
#   - a directory name where the results will be stored (option '-d', default is 'Spanning_reads')
#   - a directory with the results of the script 'process_chimCT_output.py' (option '-p', required)
#   - a breakpoint coordinates (e.g. "(15:60179225,1 / X:35772912,1 # 50,62)") (option '-b')
#   - a first feature (option '-f', required)
#   - a second feature (option '-g', required)
#   - a nucleotidic distance (option '-n')
#   - a flag for full feature mode (option '-t')
# If the mode 'full feature' is chosen, then all the paired-end reads with one mate on the first feature and the second mate on the second feature are retrieved. Otherwise, only the paired-end reads with one mate around the breakpoint are retrieved. If a feature is unknown, then, a kpb distance should be provided to search around the breakpoint.
# The outputs are:
#   - 'spanning_split_reads.fasta': a fasta file with the split reads (if 'full feature' mode is not chosen)
#   - 'spanning_PE_reads.fasta': a fasta file with the paired-end reads around the breakpoint or the paired-end reads for the 2 features (if 'full feature' mode is chosen)
#   - 'Contigs' (directory):
#      - 'file.log': the log file from CAP3
#      - 'PE_contigs.fasta': a fasta file with the contigs created by CAP3 from 'spanning_PE_reads.fasta' sequences
#      - 'PE_singlets.fasta': a fasta file with the sequences from 'spanning_PE_reads.fasta' that have not been contiged

import os, argparse, sys, re, shutil

# Getting the functions from 'classics.py' and the RLOCs and ENSCAFGs indexes
sys.path.insert(1, '/home/genouest/genouest/mbahin/Scripts')
import classics, get_indexes
RLOCs = get_indexes.RLOCs
ENSCAFGs = get_indexes.ENSCAFGs
gene_names = get_indexes.gene_names
sys.path.remove('/home/genouest/genouest/mbahin/Scripts')

###### Functions

def get_fasta(sense):
    #####
    # Function to get the fasta file of the reads with one mate in the first interest region and the other in the second one.
    # The parameter indicate wether the operation is done in the sense or the reverse sense.
    #####

    if options.feat2:
         # Getting first mates in the first area of interest
        if sense:
            os.system('samtools view -b -f 0x40 '+options.processed_dir+'/link_to_BAM_files_directory.ln/pairs.primary_alignment.sort.bam '+chr1+':'+feat1_beg+'-'+feat1_end+' > file.first_mates.bam')
        else:
            os.system('samtools view -b -f 0x80 '+options.processed_dir+'/link_to_BAM_files_directory.ln/pairs.primary_alignment.sort.bam '+chr1+':'+feat1_beg+'-'+feat1_end+' > file.first_mates.bam')

        # Getting second mates in the second area of interest
        if sense:
            os.system('samtools view -b -f 0x80 '+options.processed_dir+'/link_to_BAM_files_directory.ln/pairs.primary_alignment.sort.bam '+chr2+':'+feat2_beg+'-'+feat2_end+' > file.second_mates.bam')
        else:
            os.system('samtools view -b -f 0x40 '+options.processed_dir+'/link_to_BAM_files_directory.ln/pairs.primary_alignment.sort.bam '+chr2+':'+feat2_beg+'-'+feat2_end+' > file.second_mates.bam')
    else:
        os.system('samtools view '+options.processed_dir+'/link_to_BAM_files_directory.ln/pairs.primary_alignment.sort.bam '+chr1+':'+feat1_beg+'-'+feat1_end+' > file.merged_mates.sam')

    # Getting paired-end reads with the two mates in the two area of interest
    if options.feat2:
        os.system('samtools merge -nu file.merged_mates.bam file.first_mates.bam file.second_mates.bam')
        os.system('samtools view file.merged_mates.bam > file.merged_mates.sam')

    SAM_file = open('file.merged_mates.sam','r')
    if sense:
        output = open ('output.fasta','w')
    else:
        output = open ('output.rev.fasta','w')
    readID_list = {}
    for line in SAM_file:
        ID = line.split('\t')[0]
        if readID_list.has_key(ID):
            flag = line.split('\t')[1]
            seq = line.split('\t')[9]
            if int(flag) & 0x40:
                # First mate line found
                output.write('>'+ID+'/1\n'+seq+'\n>'+ID+'/2\n'+readID_list[ID]['seq']+'\n')
            else:
                # Second mate line found
                output.write('>'+ID+'/1\n'+readID_list[ID]['seq']+'\n>'+ID+'/2\n'+seq+'\n')
        else:
            readID_list[ID] = {}
            readID_list[ID]['seq'] = line.split('\t')[9]
    output.close()
    SAM_file.close()
    
    # Cleaning
    if options.feat2:
        os.remove('file.first_mates.bam')
        os.remove('file.second_mates.bam')
        os.remove('file.merged_mates.bam')
    os.remove('file.merged_mates.sam')

##### Functions end

# Setting the environment
os.system('. /local/env/envsamtools.sh')
GFF = '/home/genouest/umr6061/recomgen/dog/data/canFam3/annotation/MasterAnnotation/canfam3_cons_annot.gff'

# Getting options back
parser = argparse.ArgumentParser()
parser.add_argument('-d', dest='dir')
parser.add_argument('-p', dest='processed_dir', required=True)
parser.add_argument('-b', dest='bkpt')
parser.add_argument('-f', dest='feat1', required=True)
parser.add_argument('-g', dest='feat2')
parser.add_argument('-n', dest='nt', type=int, default=1)
parser.add_argument('-t', dest='total', action='store_true')
options = parser.parse_args()

if options.bkpt and (options.total == True):
    print "Breakpoint provided (option '-b') and full gene mode chosen (option '-t'). Aborting."
    sys.exit()
if not options.processed_dir.startswith('/'):
    print 'The "processed.dir" directory path specified must be an absolute path. Aborting.'
    sys.exit()

# Creating and going to a directory for the job
if options.dir != None:
    dir_name = options.dir
else:
    dir_name = 'Spanning_reads'
if not options.total:
    dir_name = dir_name+'.bkpt.dir'
else:
    dir_name = dir_name+'.full_feat.dir'

classics.create_wd(dir_name)

print '#####'

# Parsing the breakpoint
if options.bkpt:
    chr1 = options.bkpt.split('(')[1].split(':')[0]
    end1 = options.bkpt.split(':')[1].split(',')[0]
    strand1 = options.bkpt.split(',')[1].split(' ')[0]
    chr2 = options.bkpt.split(':')[1].split(' ')[2]
    start2 = options.bkpt.split(':')[2].split(',')[0]
    strand2 = options.bkpt.split(',')[2].split(' ')[0]

if not options.total:
#####
# Getting the spanning reads
#####

    # Building the spanning read file
    fasta_file = open(options.processed_dir+'/link_to_spanning_split_reads.ln.fasta','r')
    spanning_split_output = open ('spanning_split_reads.fasta','w')
    to_keep = False
    seq_nb = 1
    for line in fasta_file:
        if to_keep:
            spanning_split_output.write(line)
            to_keep = False
            continue
        if not line.startswith('>'):
            continue
        match = re.match(r'>(.*)@-?1@(.*)@(.*)@-?1@(.*)', line)
        if (match.group(1) in (chr1,chr2)) and (match.group(2) in (end1,start2)) and (match.group(3) in (chr1,chr2)) and (match.group(4) in (end1,start2)):
            to_keep = True
            spanning_split_output.write('>Spanning_read'+str(seq_nb)+'\n')
            seq_nb += 1

    spanning_split_output.close()
    fasta_file.close()

#####
# Indexing RLOCs coordinates
#####
with open(GFF,'r') as GFF_file:
    for line in GFF_file:
        if (not line.startswith('#')) and (line.split('\t')[2] == 'gene'):
            rloc = line.split('\t')[8].split(';')[0].split('=')[1]
            RLOCs[rloc]['chr'] = line.split('\t')[0]
            RLOCs[rloc]['start'] = line.split('\t')[3]
            RLOCs[rloc]['end'] = line.split('\t')[4]

#####
# Getting the paired-end reads around the breakpoint
#####
"""
# Checking the features
if options.feat1 not in feat_index:
    print 'Warning: the first feature you provided is a problem for now.\nSorry :). Aborting.\n#####'
    os.chdir('..')
    shutil.rmtree(dir_name)
    sys.exit()
if options.feat2 and (options.feat2 not in feat_index):
    print 'Warning: the second feature you provided is a problem for now.\nSorry :). Aborting.\n#####'
    os.chdir('..')
    shutil.rmtree(dir_name)
    sys.exit()
"""

# Defining the RLOCs requested
rloc1 = classics.get_RLOC_from_feat(options.feat1)
if options.feat2:
    rloc2 = classics.get_RLOC_from_feat(options.feat2)

# Getting the interval terminals
if options.total:
    if options.feat1 != 'No_match':
        feat1_beg = RLOCs[rloc1]['start']
        feat1_end = RLOCs[rloc1]['end']
    else:
        feat1_beg = str(int(end1) - 1000 * options.nt)
        feat1_end = str(int(end1) + 1000 * options.nt)
    if options.feat2:
        if options.feat2 != 'No_match':
            feat2_beg = RLOCs[rloc2]['start']
            feat2_end = RLOCs[rloc2]['end']
        else:
            feat2_beg = str(int(start2) - 1000 * options.nt)
            feat2_end = str(int(start2) + 1000 * options.nt)
else:
    if options.feat1 != 'No_match':
        if strand1 == '1':
            feat1_beg = RLOCs[rloc1]['start']
            feat1_end = end1
        else:
            feat1_beg = end1
            feat1_end = RLOCs[rloc1]['end']
    else:
        if strand1 == '1':
            feat1_beg = str(int(end1) - 1000 * options.nt)
            feat1_end = end1
        else:
            feat1_beg = end1
            feat1_end = str(int(end1) + 1000 * options.nt)
    if options.feat2:
        if options.feat2 != 'No_match':
            if strand2 == '1':
                feat2_beg = start2
                feat2_end = RLOCs[rloc2]['end']
            else:
                feat2_beg = RLOCs[rloc2]['start']
                feat2_end = start2
        else:
            if strand2 == '1':
                feat2_beg = start2
                feat2_end = str(int(start2) + 1000 * options.nt)
            else:
                feat2_beg = str(int(start2) - 1000 * options.nt)
                feat2_end = start2

# Finding chromosomes if full feature mode chosen (no breakpoint provided)
if options.total:
    chr1 = RLOCs[rloc1]['chr']
    if options.feat2:
        chr2 = RLOCs[rloc2]['chr']

# Getting a BAM file with paired-end reads in the area of interest
if options.feat2:
    print 'Searching reads within '+chr1+':'+feat1_beg+'-'+feat1_end+' and '+chr2+':'+feat2_beg+'-'+feat2_end
else:
    print 'Searching reads involving '+chr1+':'+feat1_beg+'-'+feat1_end

# Getting fasta files for sense and reverse sense fusions
if options.feat2:
    get_fasta(True)
    get_fasta(False)
else:
    get_fasta(True)

# Merging the sense and reverse fasta files
spanning_PE_output = open('spanning_PE_reads.fasta','w')
shutil.copyfileobj(open('output.fasta','r'), spanning_PE_output)
if options.feat2:
    shutil.copyfileobj(open('output.rev.fasta','r'), spanning_PE_output)
    ### TO MODIFY AND TEST add line 272 here ###
spanning_PE_output.close()
os.remove('output.fasta')
if options.feat2:
    os.remove('output.rev.fasta')

#####
# Making the contigs with CAP3
#####

# Checking if the spanning paired-end reads file is empty
if os.stat('spanning_PE_reads.fasta')[6] == 0:
    print 'Warning: "spanning_PE_reads.fasta" file is empty. No contig file can be created.\n#####'
    sys.exit()

# Creating a directory for contigs
os.makedirs('Contigs')

# Executing CAP3 command
os.system('/local/cap3/bin/cap3 spanning_PE_reads.fasta > Contigs/file.log')

# Checking if contig(s) have been created
if os.stat('spanning_PE_reads.fasta.cap.contigs')[6] == 0:
    print 'Warning: No contig could be created.'
else:
    # Renaming and reformatting files
    contigs_file = open('spanning_PE_reads.fasta.cap.contigs', 'r')
    new_contigs_file = open('Contigs/PE_contigs.fasta','w')
    seq = ''
    for line in contigs_file:
        if line.startswith('>'):
            if seq != '':
                new_contigs_file.write(seq+'\n')
                seq = ''
            new_contigs_file.write(line)
        else:
            seq = seq + line.rstrip()
    new_contigs_file.write(seq+'\n')
    new_contigs_file.close()
    contigs_file.close()

# Cleaning directory
os.rename('spanning_PE_reads.fasta.cap.singlets','Contigs/PE_singlets.fasta')
os.remove('spanning_PE_reads.fasta.cap.contigs')
os.remove('spanning_PE_reads.fasta.cap.ace')
os.remove('spanning_PE_reads.fasta.cap.contigs.links')
os.remove('spanning_PE_reads.fasta.cap.contigs.qual')
os.remove('spanning_PE_reads.fasta.cap.info')

print '#####'