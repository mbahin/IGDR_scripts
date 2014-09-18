#!/usr/bin/env python

import os, sys
import drmaa
import time

def main():
    """
    Wait for a list of jobs to finish.
    usage (to wait for 3 job ids): sge_waiter.py 415852 1452174 258852
    """
    s = drmaa.Session()
    s.initialize()
    
    jobids = []
    
    for arg in sys.argv[1:]:
        jobid = 0
        try:
            jobid = int(arg)
        except:
            print arg + " is not a valid job id (not an integer)"
        if jobid != 0:
            jobids.append(str(jobid))
    
    print "Waiting for jobs: %s" % ', '.join(jobids)
    
    while len(jobids) > 0:
        removed = False
        pos = len(jobids)-1
        try:
            print "Asking status for job " + jobids[pos]
            status = s.jobStatus(jobids[pos])
        except:
            print "Job " + jobids[pos] + " finished"
            jobids.pop()
            removed = True
        if not removed:
           print "waiting 10s"
           time.sleep(10)
        

    print 'All job finished'
    s.exit()
    
if __name__=='__main__':
    main()