#!/usr/bin/env python
import fyrd

commands = ['sleep 30','sleep 30','sleep 30',
            'sleep 30', 'sleep 30', 'sleep 30']
out = []
for cmd in commands:
    out.append(fyrd.submit(cmd, runpath = "run", outpath = "out",
                           clean_files = False, clean_outputs = False,
                           mem = '10MB', nodes = 1, time ='00:30:00',
                           threads = 1, partition = 'hbfraser'))
    print(cmd) 
