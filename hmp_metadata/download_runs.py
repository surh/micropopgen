#!/usr/bin/env python
# Copyright (C) 2017 Sur Herrera Paredes

# Reads a table with sample and run columns, and dowloads all the run
# fastq files from SRA.

# Currently it keeps all temporary files

# It only uses qsub for submitting to cluster.
# It can use direct command line

# Creates a number of submissions by splitting the List
# of runs into n groups

# Adapted from download_all_runs.pl script

import csv
import os
import tempfile
from math import ceil

def process_run_list(file,sample_col,run_col,header = True):
    print("\n=============================================")
    with open(file,'r') as fh:
        print("> Processing map of runs")
        if header is True:
            colnames = fh.readline()
        reader = csv.reader(fh, delimiter = "\t")
        
        RUNS = dict()
        nruns = 0
        nsamples = 0
        for line in reader:
            sample = line[sample_col]
            run = line[run_col]
            #print([sample,run])
            if sample in RUNS:
                RUNS[sample].append(run)
            else:
                RUNS[sample] = [run]
                nsamples += 1
            nruns += 1
    fh.close()
    print("\tProcessed {} runs in {} samples".format(nruns,nsamples))
    print("=============================================")

    return(RUNS)

def create_submission_sets(runs_per_sample,outdir,split_by,ngroups):
    print("\n=============================================")
    # Create output directory
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    else:
        print("Outdir ({}) already exists. Using it.".format(outdir))
    
    submission_dir = tempfile.mkdtemp(suffix = None, prefix = 'submissions',
                                      dir = outdir)
    
    SUBMISSIONS = []
    if split_by == 'sample':
        print("== Entering splity by sample")
        for sample, runs in runs_per_sample.items():
            print("NOTHING")
    elif split_by == 'groups':
        samples = runs_per_sample.keys()
        total_samples = len(samples)
        samples_per_submission = ceil(total_samples / ngroups)
        
        print("\tSplitting {} samples into {} submissions".format(total_samples,ngroups))
        GROUPS = []
        i = 0
        group_i = 0
        for sample, runs in runs_per_sample.items():
            if (i % samples_per_submission) == 0:
                GROUPS.append([])
                group_i += 1
            GROUPS[group_i - 1].extend(runs)
            i += 1
        
        for group in GROUPS:
            #newfile = create_single_submission_file(group,submission_dir.name,outdir)
            #SUBMISSIONS.append(newfile)
            print(group)
    else:
        raise ValueError("Unrecognized split_by value")
    
    return(GROUPS)


# print(__name__)
# Run if called as script
if __name__ == "__main__":
    import argparse
        
    # Define arguments
    parser = argparse.ArgumentParser()
    required = parser.add_argument_group("Required arguments")
    required.add_argument("--infile","-i", help = "Input tab delimited file that maps runs (SRR) and samples (SRS)",
                        type = str, required = True)
    required.add_argument("--outdir","-o", help = "Output directory", type = str,
                          required = True)
    parser.add_argument("--method","-m", help = "method to submit jobs to the server", type = str,
                        default = 'qsub', choices = ['qsub','bash'])
    parser.add_argument("--sample_col", help = "Column number where the sample ID (SRS) is stored", type = int,
                        default = 1)
    parser.add_argument("--run_col", help = "Column number where the run ID (SRR) is stored", type = int,
                        default = 2)
    parser.add_argument("--header", help = "Flag indicating whether table has headers in the first row",
                        action = "store_true")
    parser.add_argument("--split_by", help = "Method to split runs for download", default = 'sample',
                        choices = ['sample','groups'])
    parser.add_argument('--logdir','-l', help = "Directory where to store log of run", default = './logs/',
                        type = str)
    parser.add_argument('--ngroups', help = "Number of groups to divide the runs into. Equal to the number of jobs that will be submitted",
                        default = 2, type = int)
    args = parser.parse_args()
    
    # Process column numbers
    args.run_col -= 1
    args.sample_col -= 1
    #print(args.infile)
    #print(args.run_col)
    #print(args.sample_col)
    
    runs_per_sample = process_run_list(args.infile, args.sample_col,
                                       args.run_col, args.header)
    create_submission_sets(runs_per_sample, args.outdir,
                           args.split_by, args.ngroups)


