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
    FILES = [[], []]
    for run in runs:
        run_sra = indir + "/" + run + ".sra"
        
        if os.path.exists(run_sra):
            command = 'fastq-dump -O ' + outdir + ' --split-files ' + run_sra + " &"
            check = download_runs.run_command(command)
            #check = subprocess.run('fastq-dump -O ' + outdir + ' --split-files ' + run_sra + " &", shell = True)
            if check != 0:
                raise ERROR("\rRun {} could not be processed by fastq-dump".format(run))
            else:
                read1 = outdir + "/" + run + "_1.fastq"                
                FILES[0].append(read1)
                read1 = outdir + "/" + run + "_2.fastq"                
                FILES[1].append(read2)
        else:
            raise ERROR("\tRun {} file does not exist in {}".format(run,outdir))
        
        return(FILES)

def concatenate_files(infiles, outfile):
    command = " ".join(infiles) 
    command = "cat " + command + " > " + outfile
    check = download_runs.run_command(command)
    
    if check != 0:
        raise ERROR("Could not concatenate files")
    
    return(check)

def concatenate_run(file_sets,outdir,name_prefix, extension = ".fastq"):
    i = 1
    FILES = []
    for files in file_sets:
        newfile = outdir + "/" + name_prefix + "_read" + i + extension
        try:
            concatenate_files(files, newfile)
            i += 1
            FILES.append(newfile)
        except (ERROR):
            raise ERROR("Could not concatenate files from read {}".format(i))
    
    return(FILES)
         
    
def process_sample(sample,runs,indir,outdir):
    
    # Validate files
    try:
        check_set_of_runs(runs, indir)
    except (ERROR):
        raise ERROR("\tSample didn't pass check")
    
    # Proceed to fastq-dump
    try:
        run_fastq = astq_dump_runs(runs,indir,outdir)
    except (ERROR):
        raise ERROR("\tSample {} could not be processed by fastq dump".format(sample))
    
    # Proceed to concatenate
    try:
        concatenated_files = concatenate_run(run_fastq, outdir, sample, ".fastq")
    except:
        raise ERROR("Could not concatenate files from sample {}".format(sample))
    
    return(concatenated_files)
        
        

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser()
    
    # Required arguments
    required = parser.add_argument_group("Required arguments")
    required.add_argument("--indir","-i", help = "Directory containing .sra files from all the runs to be processed",
                          type = str, required = True)
    required.add_argument("--fastq_dir","-f", help = "Directory where to place the fastq files that are produced from sra files",
                          type = str, required = True)
    required.add_argument("--outdir","-o", help = "Directory where to place the final fastq files, one per sample",
                          type = str, required = True)    
    required.add_argument("--map","-m", help = "Input tab delimited file that maps runs (SRR) and samples (SRS)",
                        type = str, required = True)
    
    
    
    
    
    indir = "./runs/"
    outdit = "./samples/"
    runs_file = "/home/sur/micropopgen/exp/2017/today8/runs_to_download.txt"
    
    runs_per_sample = download_runs.process_run_list(runs_file, 0, 1, True)
    
    for sample in runs_per_sample.keys():
        print(sample)
        #process_sample(sample, runs_per_sample[sample], indir, outdir)

        
    