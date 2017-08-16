#!/usr/bin/env python
import fyrd

commands = ['sleep 180','sleep 180','sleep 180',
            'sleep 180', 'sleep 180', 'sleep 180']
out = []
for cmd in commands:
    out.append(fyrd.submit(cmd, runpath = "run", outpath = "out",
                           clean_files = False, clean_outputs = False,
                           mem = '10MB', nodes = 1, time ='00:30:00',
                           cores = 1, partition = 'hbfraser',
                           scriptpath = 'sub', outfile = 'test.log',
                           errfile = 'test.err'))
    print(cmd) 
