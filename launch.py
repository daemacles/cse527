#!/usr/bin/env python
import sys
import os
from subprocess import Popen
import parallel_run

if __name__ == '__main__':
    try:
        iterations = int(sys.argv[1])
        num_bp = int(sys.argv[2])
        num_saved = int(sys.argv[3])
    except:
        print 'Usage: %s <Num iterations> <Num basepairs> <Num saved>'%sys.argv[0]
        sys.exit(1)

    clients = [('diglett.cs.washington.edu', 8),
               ('charliebrown.cs.washington.edu', 8),
               ('hobbes.cs.washington.edu', 7)]
    procs = []
    # Start the worker processes
    for client,num_jobs in clients:
        procs.append(Popen(['ssh', client, 'cd ~/Dropbox/cse527/project; nice -n 19 python parallel_run.py w %d %d'%(num_jobs, num_bp)]))
    print 'here'

    # Start the analyzing process
    parallel_run.runServer(iterations, num_saved)

    for proc in procs:
        proc.wait()
    if any(proc.returncode != 0 for proc in procs):
        print 'Something failed'
