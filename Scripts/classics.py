#!/local/python/2.7/bin/python

# Mathieu Bahin, 15/09/14

# Module gathering functions found in a lot of my scripts.

import os

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