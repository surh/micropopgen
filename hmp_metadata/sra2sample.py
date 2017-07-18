#!/usr/bin/env python
# Copyright (C) 2017 Sur Herrera Paredes


import download_runs
import subprocess
import os

class Error(Exception):
    """Base class for exceptions in this module."""
    pass

class IntegrityError(Error):
    """Exception raised for failure in the vdb-vvalidate."""
    pass

class MissingFileError(Error):
    """Exception raised for missing files."""
    pass

class ProcessError(Error):
    """Exception raised for failure in the some sra-tools call."""
    pass


def check_set_of_runs(runs, dir):
    print("\t Checking runs..,")
    for run in runs:
        run_sra = dir + "/" + run + ".sra"
        check = 0
        if os.path.exists(run_sra):
            command = 'vdb-validate ' + run_sra
            check = download_runs.run_command(command)
            #check = subprocess.run('vdb-validate ' + run_sra + " &", shell = True)
            if check.returncode != 0:
                raise IntegrityError("\rRun {} did not pass the validation".format(run))
        else:
            raise MissingFileError("\tRun {} file does not exist in {}".format(run,outdir))
        
    return(check)

def fastq_dump_runs(runs,indir,outdir,keep):
    if not os.path.isdir(indir):
        raise ERROR("Directory {} does not exisst".format(indir))
    if not os.path.isdir(outdir):
        print("\tCreating output directory {}".format(outdir))
        os.mkdir(outdir)
        
    FILES = [[], []]
    for run in runs:
        run_sra = indir + "/" + run + ".sra"
        
        if os.path.exists(run_sra):
            command = 'fastq-dump -I -O ' + outdir + ' --split-files --bzip2 ' + run_sra
            check = download_runs.run_command(command)
            #check = subprocess.run('fastq-dump -O ' + outdir + ' --split-files ' + run_sra + " &", shell = True)
            if check.returncode != 0:
                raise ProcessError("\rRun {} could not be processed by fastq-dump".format(run))
            else:
                read1 = outdir + "/" + run + "_1.fastq.bz2"                
                FILES[0].append(read1)
                read2 = outdir + "/" + run + "_2.fastq.bz2"                
                FILES[1].append(read2)
        else:
            raise MissingFileError("\tRun {} file does not exist in {}".format(run,outdir))
        
    return(FILES)

def concatenate_files(infiles, outfile):
    command = " ".join(infiles) 
    command = "cat " + command + " > " + outfile
    check = download_runs.run_command(command)
    
    if check != 0:
        raise ProcessError("Could not concatenate files")
    
    return(check)

def concatenate_run(file_sets,outdir,name_prefix, extension = ".fastq"):
    if not os.path.isdir(outdir):
        print("\tCreating output directory {}".format(outdir))
        os.mkdir(outdir)
    
    i = 1
    FILES = []
    for files in file_sets:
        newfile = outdir + "/" + name_prefix + "_read" + str(i) + extension
        try:
            concatenate_files(files, newfile)
            i += 1
            FILES.append(newfile)
        except (ProcessError):
            raise ProcessError("Could not concatenate files from read {}".format(i))
    
    return(FILES)
         
    
def process_sample(sample,runs,indir,fastqdir,outdir,keep = False):
    
    # Validate files
    try:
        check_set_of_runs(runs,indir)
    except (IntegrityError, MissingFileError):
        raise ProcessError("\tSample didn't pass check")
    
    # Proceed to fastq-dump
    try:
        run_fastq = fastq_dump_runs(runs,indir,fastqdir,keep)
    except (ProcessError, MissingFileError):
        raise ProcessError("\tSample {} could not be processed by fastq dump".format(sample))
    
    # Proceed to concatenate
    try:
        concatenated_files = concatenate_run(run_fastq, outdir, sample, ".fastq.bz2")
    except:
        raise ProcessError("Could not concatenate files from sample {}".format(sample))
    
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
    
    # Optional arguments
    parser.add_argument("--sample_col", help = "Column where the sample name is located in map", type = int,
                        default = 1)
    parser.add_argument("--run_col", help = "Column where the run accession is located in map. Run names must match names of sra files", type = int,
                        default = 2)
    parser.add_argument("--keep_intermediate", help = "Flag indicating whether to keep the intermediate fastq files.", action = "store_true")
    parser.add_argument("--header", help = "Flag indicating whether table has headers in the first row",
                        action = "store_true")
    
    args = parser.parse_args()
    args.sample_col -= 1
    args.run_col -= 1
    
#     
#     indir = "./runs/"
#     outdit = "./samples/"
#     runs_file = "/home/sur/micropopgen/exp/2017/today8/runs_to_download.txt"
#     
    runs_per_sample = download_runs.process_run_list(args.map, args.sample_col, args.run_col, args.header)
    
    for sample in runs_per_sample.keys():
        print("== Processing sample {}".format(sample))
        print(" ".join(runs_per_sample[sample]))
        process_sample(sample, runs_per_sample[sample], args.indir, args.fastq_dir, args.outdir, args.keep_intermediate)

        
    