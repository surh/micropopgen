#!/usr/bin/env python
import fyrd

commands = ['sleep 180','sleep 180','sleep 180',
            'sleep 180', 'sleep 180', 'sleep 180']
JOBS = []
for cmd in commands:
    JOBS.append(fyrd.Job(cmd, runpath = "run", outpath = "out",
                         clean_files = False, clean_outputs = False,
                         mem = '10MB', nodes = 1, time ='00:30:00',
                         cores = 1, partition = 'hbfraser',
                         scriptpath = 'sub', outfile = 'test.log',
                         errfile = 'test.err', max_jobs = 2))
    
    (JOBS[-1]).submit
    print(cmd)
    print(JOBS[-1].)
