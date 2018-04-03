#!/usr/bin/env python

import fyrd
import os
import time


def testfun(s,suf = '.txt'):
	print('.'.join([s,suf]))


#job1 = fyrd.Job('sleep 15; echo "Hello" > myfile.test', runpath=os.getcwd())
#job1.submit()

#j = job1

#job2 = fyrd.Job('ls -l myfile.test', runpath=os.getcwd(), depends=j)
#job2.submit()

s = "Hello"
job2 = fyrd.Job(testfun,s, runpath=os.getcwd())
job2.submit()
