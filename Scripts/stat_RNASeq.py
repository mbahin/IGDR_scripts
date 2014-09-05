#!/local/python/2.7/bin/python

# Mathieu Bahin, 05/09/14

# Script to produce statistics on the Google spredsheet about RNA-Seqs. It has to be download in the 'tsv' format first.
# For each level of the structure there are counts for the total, the affected and the healthy.

import sys

# Indexing file information
projects = {}
breeds = {}
with open(sys.argv[1],'r') as RNASeq_file:
    RNASeq_file.readline()
    for line in RNASeq_file:
        breed = line.split('\t')[3]
        status = line.split('\t')[7]
        project = line.split('\t')[9]
        center = line.split('\t')[16]
        if project == 'Domestication':
            continue

        # Indexing project
        if project not in projects:
            projects[project] = {}
            projects[project]['count'] = 0
            projects[project]['Healthy'] = 0
            projects[project]['Affected'] = 0
        projects[project]['count'] += 1
        projects[project][status] += 1
        # Indexing breed
        if breed not in projects[project]:
            projects[project][breed] = {}
            projects[project][breed]['count'] = 0
            projects[project][breed]['Healthy'] = 0
            projects[project][breed]['Affected'] = 0
        projects[project][breed]['count'] += 1
        projects[project][breed][status] += 1
        if center not in projects[project][breed]:
            projects[project][breed][center] = {}
            projects[project][breed][center]['count'] = 0
            projects[project][breed][center]['Healthy'] = 0
            projects[project][breed][center]['Affected'] = 0
        projects[project][breed][center]['count'] += 1
        projects[project][breed][center][status] += 1

# Printing out results
for project in projects:
    print '#####\n# '+project+'\n#####'
    print str(projects[project]['count'])+' samples (Tumoral: '+str(projects[project]['Affected'])+' / Healthy: '+str(projects[project]['Healthy'])+')'
    for b in projects[project]:
        if b in (('count','Healthy','Affected')):
            continue
        print '\t'+str(projects[project][b]['count'])+' '+b+' (Tumoral: '+str(projects[project][b]['Affected'])+' / Healthy: '+str(projects[project][b]['Healthy'])+')'
        for c in projects[project][b]:
            if c in (('count','Healthy','Affected')):
                continue
            print '\t\t'+str(projects[project][b][c]['count'])+' '+c+' (Tumoral: '+str(projects[project][b][c]['Affected'])+' / Healthy: '+str(projects[project][b][c]['Healthy'])+')'
    print '\n'