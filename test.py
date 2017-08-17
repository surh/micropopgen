#!/usr/bin/env python
import fyrd
import time

commands = ['sleep 60','sleep 60','sleep 60',
            'sleep 60', 'sleep 60', 'sleep 60']

def print_number(i):
    print("== Start {}".format(str(i)))
    time.sleep(30)
    print("== End {}".format(str(i)))
    return(i)

JOBS = []
iter = 0
for cmd in commands:
#     JOBS.append(fyrd.Job(cmd, runpath = "run", outpath = "out",
#                          clean_files = False, clean_outputs = False,
#                          mem = '10MB', nodes = 1, time ='00:30:00',
#                          cores = 1, partition = 'hbfraser',
#                          scriptpath = 'sub', outfile = 'test.log',
#                          errfile = 'test.err'))
    job = fyrd.Job(print_number(iter), runpath = "run", outpath = "out",
                     clean_files = False, clean_outputs = False,
                     mem = '10MB', nodes = 1, time ='00:30:00',
                     cores = 1, partition = 'hbfraser',
                     scriptpath = 'sub', outfile = 'test.log',
                     errfile = 'test.err')
    job.submit(max_jobs = 3)
    iter = iter + 1
         
#     (JOBS[-1]).submit(max_jobs = 3)
#     print(cmd)
    #print((JOBS[-1]).command)


# cmd = 'sleep 180'
# job1 = fyrd.Job(cmd,partition = "owners", outpath = 'out', scriptpath = 'sub')
# job2 = fyrd.Job(cmd,partition = "owners", outpath = 'out', scriptpath = 'sub')
# job3 = fyrd.Job(cmd,partition = "owners", outpath = 'out', scriptpath = 'sub')
# job4 = fyrd.Job(cmd,partition = "owners", outpath = 'out', scriptpath = 'sub')
# job5 = fyrd.Job(cmd,partition = "owners", outpath = 'out', scriptpath = 'sub')
# job6 = fyrd.Job(cmd,partition = "owners", outpath = 'out', scriptpath = 'sub')
# 
# job1.submit(max_jobs=1)
# job2.submit(max_jobs=2)
# job3.submit(max_jobs=3)
# job4.submit(max_jobs=4)
# job5.submit(max_jobs=5)
# job6.submit(max_jobs=6)