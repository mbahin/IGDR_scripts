#!/local/python/2.7/bin/python

# Mathieu Bahin, 15/09/14

# Module gathering functions found in a lot of my scripts.

import os, sys

# Getting the RLOCs and ENSCAFGs indexes
sys.path.insert(1, '/home/genouest/genouest/mbahin/Scripts')
import get_indexes
RLOCs = get_indexes.RLOCs
ENSCAFGs = get_indexes.ENSCAFGs
gene_names = get_indexes.gene_names
sys.path.remove('/home/genouest/genouest/mbahin/Scripts')

def create_wd(dir_name):
    #####
    # Function to create and got to a new working directory.
    #####
    
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)
        os.chdir(dir_name)
    else:
        print 'The directory '+dir_name+' already exists. Aborting.'
        sys.exit()

def get_RLOC_from_feat(feat):
    #####
    # Getting the RLOC corresponding to a requested feature (ENSCAFG or gene_name) from the RLOCs/ENSCAFGs indexes.
    #####

    if feat.startswith('RLOC'):
        return feat
    elif feat.startswith('ENSCAFG'):
        return ENSCAFGs[feat]['RLOC']
    else:
        rloc = ''
        for enscafg in gene_names[feat]:
            if rloc and (ENSCAFGs[enscafg]['RLOC'] != rloc):
                print 'The gene name you requested can refer to several ENSCAFGs, please provide the ENSCAFG you are interested in. Aborting.'
                sys.exit()
            else:
                rloc = ENSCAFGs[enscafg]['RLOC']
        return rloc