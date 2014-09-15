#!/local/python/2.7/bin/python

# Mathieu Bahin, 12/09/14

# Module to index the RLOCs and ENSCAFGs index files into dictionaries.

RLOCs = {}
ENSCAFGs = {}
gene_names = {}

def index_RLOCs():
    #####
    # 1) Indexing the RLOCs from the file 'RLOCs_index.txt' into a dictionary, 'RLOCs'. There are 3 fields :
    #   - 'enscafgs' (list) : the ENSCAFG(s) list associated to the RLOC.
    #   - 'orthologous' : a human orthologous gene name if one can be found and if there is no ENSCAFG associated to the RLOC.
    #   - 'biotypes' (list) : the biotypes() (uniqued) list corresponding to the ENSCAFG(s) list.
    #
    # 2) Indexing the ENSCAFGs from the file 'ENSCAFGs_index.txt' into a dictionary, 'ENSCAFGs'. There are 4s fields :
    #   - 'Ensembl_name' : the Ensembl gene name.
    #   - 'BROAD_name' : the BROAD gene name.
    #   - 'consensus_names' (list) : the (uniqued) gene name(s) list from Ensembl and the BROAD.
    #   - 'RLOC' : the RLOC corresponding to the ENSCAFG.
    # The gene names are indexed into the dictionary 'gene_names' to link to the ENSCAFGs.
    #####

    with open('/home/genouest/umr6061/recomgen/dog/data/canFam3/annotation/Correspondence_Indexes/RLOCs_index.txt','r') as RLOCs_file:
        # Ignoring the header
        RLOCs_file.readline()
        # Processing each line
        for line in RLOCs_file:
            rloc = line.split('\t')[0]
            if not RLOCs.has_key(rloc):
                RLOCs[rloc] = {}
            RLOCs[rloc]['enscafgs'] = line.split('\t')[1].split(',')
            RLOCs[rloc]['orthologous'] = line.split('\t')[2]
            RLOCs[rloc]['biotypes'] = line.rstrip().split('\t')[3].split(',')

    with open('/home/genouest/umr6061/recomgen/dog/data/canFam3/annotation/Correspondence_Indexes/ENSCAFGs_index.txt','r') as ENSCAFGs_file:
        # Ignoring the header
        ENSCAFGs_file.readline()
        # Processing each line
        for line in ENSCAFGs_file:
            enscafg = line.split('\t')[0]
            if not ENSCAFGs.has_key(enscafg):
                ENSCAFGs[enscafg] = {}
            ENSCAFGs[enscafg]['Ensembl_name'] = line.split('\t')[1]
            ENSCAFGs[enscafg]['BROAD_name'] = line.split('\t')[2]
            ENSCAFGs[enscafg]['consensus_names'] = line.split('\t')[3].split('|')
            if ENSCAFGs[enscafg]['consensus_names'] != ['NA']:
                for name in line.split('\t')[3].split('|'):
                    if not gene_names.has_key(name):
                        gene_names[name] = []
                    gene_names[name].append(enscafg)
            ENSCAFGs[enscafg]['RLOC'] = line.rstrip().split('\t')[4]

index_RLOCs()