#!/usr/bin/env python
import fyrd

command1 = "ls -la"
command2 = "sleep 10"
 
# job1  = fyrd.submit(command1, runpath = "run", outpath = "out",
#                     clean_files = False, clean_outputs = False)
# print(job1)
# print(job1.get())
# job2 = fyrd.submit(command2, depends=job1)
# out  = job2.get()  # Will block until job completes



commands = ['ls -la /', 'ls -la /home/', 'pwd',
            'ls $SCRATHC', 'ls $PI_HOME','sleep 30',
            'sleep 30', 'sleep 30']
out = []
for cmd in commands:
    out.append(fyrd.submit(cmd, runpath = "run", outpath = "out",
                           clean_files = False, clean_outputs = False,
                           mem = '10MB', nodes = 1, time ='00:30:00',
                           threads = 1, partition = hbfraser))
    print(cmd) 

# for i in out:
#     print(i.get())

    



# def raise_me(something, power=2):
#     print(something)
#     return something**power
# outs = []
# if __name__ == '__main__':
#     for i in range(80):
#         outs.append(fyrd.submit(raise_me, (i,), {'power': 2},
#                                 mem='10MB', time='00:00:30'))
#     final_sum = 0
#     for i in outs:
#         final_sum += i.get()
#     print(final_sum)