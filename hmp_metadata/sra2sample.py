#!/usr/bin/env python
# Copyright (C) 2017 Sur Herrera Paredes


import download_runs
import subprocess
import os

def check_set_of_runs(runs, dir):
    for run in runs:
        run_sra = dir + "/" + run + ".sra"
        check = 0
        if os.path.exists(run_sra):
            command = 'vdb-validate ' + run_sra + " &"
            check = download_runs.run_command(command)
            #check = subprocess.run('vdb-validate ' + run_sra + " &", shell = True)
            if check != 0:
                raise ERROR("\rRun {} did not pass the validation".format(run))
        else:
            raise ERROR("\tRun {} file does not exist in {}".format(run,outdir))
        
        return(check)

def fastq_dump_runs(runs,indir,outdir):
    for run in runs:
        run_sra = indir + "/" + run + ".sra"
        
        if os.path.exists(run_sra):
            command = 'fastq-dump -O ' + outdir + ' --split-files ' + run_sra + " &"
            check = download_runs.run_command(command)
            #check = subprocess.run('fastq-dump -O ' + outdir + ' --split-files ' + run_sra + " &", shell = True)
            if check != 0:
                raise ERROR("\rRun {} could not be processed by fastq-dump".format(run))
        else:
            raise ERROR("\tRun {} file does not exist in {}".format(run,outdir))
        
        return(check)
        
def process_sample(sample,runs,indir,outdir):
    
    # Validate files
    try:
        check_set_of_runs(runs, indir)
    except (ERROR):
        raise ERROR("\tSample didn't pass check")
    
    # Proceed to fastq-dump
    try:
        fastq_dump_runs(runs,indir,outdir)
    except (ERROR):
        raise ERROR("\tSample could not be processed by fastq dump")
    
    # Proceed to concatenate
        
        

if __name__ == "__main__":
    import argparse
    
    indir = "./runs/"
    outdit = "./samples/"
    runs_file = "/home/sur/micropopgen/exp/2017/today8/runs_to_download.txt"
    
    runs_per_sample = download_runs.process_run_list(runs_file, 0, 1, True)
    
    for sample in runs_per_sample.keys():
        print(sample)
        #process_sample(sample, runs_per_sample[sample], indir, outdir)

        
    