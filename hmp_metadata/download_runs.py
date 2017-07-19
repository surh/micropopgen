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
import subprocess
from subprocess import CalledProcessError
from multiprocessing.dummy import Pool
import itertools

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

def create_submission_sets(runs_per_sample, split_by, ngroups):
    print("\n=============================================")
    GROUPS = dict()
    if split_by == 'sample':
        print("== Entering splity by sample")
        GROUPS = runs_per_sample
    elif split_by == 'groups':
        samples = runs_per_sample.keys()
        total_samples = len(samples)
        samples_per_submission = ceil(total_samples / ngroups)
        
        print("== Splitting {} samples into {} submissions".format(total_samples,ngroups))
        i = 0
        group_i = 0
        for sample, runs in runs_per_sample.items():
            if (i % samples_per_submission) == 0:
                #GROUPS.append([])
                id = 'group' + str(group_i)
                GROUPS[id] = []
                group_i += 1
            GROUPS[id].extend(runs)
            i += 1
    else:
        raise ValueError("Unrecognized split_by value")
    
    print("\tSplitted runs")
    print("=============================================")
    return(GROUPS)

def create_submission_files(groups, outdir, logdir):
    print("\n=============================================")
    # Create output directory
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    else:
        print("== Outdir ({}) already exists. Using it.".format(outdir))
    
    submission_dir = tempfile.mkdtemp(suffix = None, prefix = 'submissions',
                                      dir = outdir)
    
    SUBMISSIONS = []
    #i = 1
    for name, runs in groups:
        #name = str(i)
        #name = "group" + name + ".bash"
        #print(name)
        newfile = create_single_submission(name,runs,
                                           submission_dir,
                                           outdir,logdir)
        SUBMISSIONS.append(newfile)
        #print(newfile)
        #print(group)
        #i += 1
    
    print("\tCreated submision files")
    print("=============================================")
    return(SUBMISSIONS)
    

# def create_submission_sets(runs_per_sample,outdir,split_by,ngroups, logdir = "logs/"):
#     print("\n=============================================")
#     # Create output directory
#     if not os.path.exists(outdir):
#         os.mkdir(outdir)
#     else:
#         print("Outdir ({}) already exists. Using it.".format(outdir))
#     
#     submission_dir = tempfile.mkdtemp(suffix = None, prefix = 'submissions',
#                                       dir = outdir)
#     #print(submission_dir)
#     
#     SUBMISSIONS = []
#     if split_by == 'sample':
#         raise ValueError("Not implemented submission by sample, try perl script")
#         print("== Entering splity by sample")
#         for sample, runs in runs_per_sample.items():
#             print("NOTHING")
#     elif split_by == 'groups':
#         samples = runs_per_sample.keys()
#         total_samples = len(samples)
#         samples_per_submission = ceil(total_samples / ngroups)
#         
#         print("\tSplitting {} samples into {} submissions".format(total_samples,ngroups))
#         GROUPS = []
#         i = 0
#         group_i = 0
#         for sample, runs in runs_per_sample.items():
#             if (i % samples_per_submission) == 0:
#                 GROUPS.append([])
#                 group_i += 1
#             GROUPS[group_i - 1].extend(runs)
#             i += 1
#         
#         i = 1
#         for group in GROUPS:
#             name = str(i)
#             name = "group" + name + ".bash"
#             #print(name)
#             newfile = create_single_submission(name,group,
#                                                submission_dir,
#                                                outdir,logdir)
#             SUBMISSIONS.append(newfile)
#             #print(newfile)
#             #print(group)
#             i += 1
#     else:
#         raise ValueError("Unrecognized split_by value")
#     
#     print("\tSplitted runs")
#     print("=============================================")
#     return(SUBMISSIONS)

def create_single_submission(name, group,submission_dir,outdir,logdir):
    submission_file = submission_dir + "/" + name
    with open(submission_file,'w') as fh:
        fh.write("#!/bin/bash\n")
        fh.write("#PBS -N download." + name + "\n")
        fh.write("#PBS -d " + outdir + "\n")
        fh.write("#PBS -o " + logdir + "/download." + name + ".log\n")
        fh.write("#PBS -e " + logdir + "/download." + name + ".err\n")
        fh.write("#PBS -l mem=1000mb\n")
        
        # Add lines for every run in sample
        ascp_command = 'ascp -i /godot/hmp/aspera/asperaweb_id_dsa.openssh -k 1 -T -l200m'
        sra_prefix = 'anonftp\@ftp.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/SRR/'
        for run in group:
            run_location = run[0:6] + "/" + run + "/" + run + ".sra"
            command = " ".join([ascp_command, sra_prefix + "/" + run_location, outdir + "\n"])
            fh.write(command)
    fh.close()
    os.chmod(submission_file, 0o744)
    
    return(submission_file)

def run_command(command):
    status = 0;
    print("Executing:\n>{}".format(command))
    status = subprocess.run(command, shell = True)
    print("Status={}\n\n".format(status));

    return(status)

def qsub_submissions(submissions,logdir):
    
    if os.path.exists(logdir):
        print("Logdir exists")
    else:
        os.mkdir(logdir)
        
    for file in submissions:
        run_command("qsub " + file)
        
    print("==========SUBMISSIONS DONE==========\n\n")

def aspera_download(groups,outdir):
    
    ascp_command = 'ascp -i /godot/hmp/aspera/asperaweb_id_dsa.openssh -k 1 -T -l200m'
    sra_prefix = 'anonftp\@ftp.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/SRR/'
    
    FAILED = []
    for name, runs in groups.items():
        for run in runs:
            run_location = run[0:6] + "/" + run + "/" + run + ".sra"
            command = " ".join([ascp_command, sra_prefix + "/" + run_location, outdir + "\n"])
            print(command)
#             try:
#                 check = run_command(command)
#                 if check.returncode != 0:
#                     raise CalledProcessError("Aspera download failed"                                         )
#             except (CalledProcessError):
#                 print("\tWARNING: Failed downloading run {}".format(run))
#                 FAILED.append([run])
    FAILED = [[1],[2]]
    return(FAILED)

def write_table(outfile,rows, header = None, delimiter = "\t", verbose = False):
    with open(outfile,'w') as out_fh:
        writer = csv.writer(out_fh,delimiter = '\t')
        if verbose:
            print("\tWriting {}".format(outfile))
            
        nlines = 0
        if header is not None:
            writer.writerow(header)
            nlines += 1
        for row in rows:
            writer.writerow(row)
            nlines += 1
    out_fh.close()

    if verbose:
        print("\t\tWrote {} lines".format(nlines))
    
    return(nlines)    

  
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
                        default = 'qsub', choices = ['qsub','bash','python'])
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
    parser.add_argument('--threads', help = "Number of threads to use", default = 1,
                        type = int, choices = [1, 2, 3, 4])
    args = parser.parse_args()
    
    # Process column numbers
    args.run_col -= 1
    args.sample_col -= 1
    #print(args.infile)
    #print(args.run_col)
    #print(args.sample_col)
    
    runs_per_sample = process_run_list(args.infile, args.sample_col,
                                       args.run_col, args.header)
    submissions = create_submission_sets(runs_per_sample,
                                         args.split_by,
                                         args.ngroups)
    #print(submissions)
    # Submit files
    if args.method == 'qsub':
        submission_files = create_submission_files(submissions, args.outdir, args.logdir)
        qsub_submissions(submission_files,args.logdir)
    elif args.method == 'bash':
        submission_files = create_submission_files(submissions, args.outdir, args.logdir)
        for sub in submission_files:
            run_command(sub + " &")
    elif args.method == 'python':
        print("python")
        if args.threads > 1:
            # Convert dictionary into list of dictionaries
            submissions_threading = ()
            for sample in submissions:
                submissions_threading.append({sample : submissions[sample]})
            
            #print(submissions_threading)
            #for i in zip(submissions_threading, itertools.repeat(args.outdir)):
            #    print(i)
            
            # Then call threading
            pool = Pool(args.threads)
            failed = pool.starmap(aspera_download,
                                  submissions_threading)
            pool.close()
            pool.join()
        else:
            failed = aspera_download(submissions, args.outdir)
        
        print(failed)
        if len(failed) > 0:
            write_table('failed.txt', failed)
    else:
        raise ValueError("Method ($method) not recognized")

    
 
    

