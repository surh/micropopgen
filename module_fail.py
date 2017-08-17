#!/usr/bin/env python
import fyrd

command = 'blastp -h'

job = fyrd.Job(command, modules = "ncbi-blast+",
               partition = owners, time = "00:30:00",
               mem = "100MB")
job.submit()