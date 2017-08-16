#!/usr/bin/env python
import fyrd

# commands = ['sleep 180','sleep 180','sleep 180',
#             'sleep 180', 'sleep 180', 'sleep 180']
# JOBS = []
# for cmd in commands:
#     JOBS.append(fyrd.Job(cmd, runpath = "run", outpath = "out",
#                          clean_files = False, clean_outputs = False,
#                          mem = '10MB', nodes = 1, time ='00:30:00',
#                          cores = 1, partition = 'hbfraser',
#                          scriptpath = 'sub', outfile = 'test.log',
#                          errfile = 'test.err', max_jobs = 2))
#     
#     (JOBS[-1]).submit
#     print(cmd)
#     print((JOBS[-1]).command)


cmd = 'sleep 180'
job1 = fyrd.Job(cmd,partition = "owners", outpath = 'out', scriptpath = 'sub')
job2 = fyrd.Job(cmd,partition = "owners", outpath = 'out', scriptpath = 'sub')
job3 = fyrd.Job(cmd,partition = "owners", outpath = 'out', scriptpath = 'sub')
job4 = fyrd.Job(cmd,partition = "owners", outpath = 'out', scriptpath = 'sub')
job5 = fyrd.Job(cmd,partition = "owners", outpath = 'out', scriptpath = 'sub')
job6 = fyrd.Job(cmd,partition = "owners", outpath = 'out', scriptpath = 'sub')

job1.submit(max_jobs=1)
job2.submit(max_jobs=2)
job3.submit(max_jobs=3)
job4.submit(max_jobs=4)
job5.submit(max_jobs=5)
job6.submit(max_jobs=6)