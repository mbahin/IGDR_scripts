#!/local/python/2.7/bin/python

# Mathieu Bahin, 12/09/14

# Module to index the RLOCs and ENSCAFGs index files into dictionaries.

RLOCs = {}
ENSCAFGs = {}

def index_RLOCs():
    global RLOCs
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
        ENSCAFGs[enscafg]['consensus_name'] = line.rstrip().split('\t')[3]