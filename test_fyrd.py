import fyrd
import os
import time


job1 = fyrd.Job('sleep 10; echo "Hello" > myfile.test', runpath=os.getcwd())
job1.submit()
job2 = fyrd.Job('ls -l myfile.test', runpath=os.getcwd())
