#!/usr/bin/env python
# Copyright (C) 2017 Sur Herrera Paredes


#import download_runs
import subprocess
import os
import sutilspy
import itertools

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
            check = sutilspy.io.run_command(command)
            #check = subprocess.run('vdb-validate ' + run_sra + " &", shell = True)
            if check.returncode != 0:
                raise IntegrityError("\rRun {} did not pass the validation".format(run))
        else:
            raise MissingFileError("\tRun {} file does not exist in {}".format(run,dir))
        
    return(check)

def fastq_dump_runs(runs,indir,outdir,keep):
    if not os.path.isdir(indir):
        raise FileNotFoundError("Input directory {} does not exist".format(indir))
    if not os.path.isdir(outdir):
        print("\tCreating output directory {}".format(outdir))
        os.mkdir(outdir)
        
    FILES = [[], []]
    for run in runs:
        run_sra = indir + "/" + run + ".sra"
        
        if os.path.exists(run_sra):
            command = 'fastq-dump -I -O ' + outdir + ' --split-files --bzip2 ' + run_sra
            check = sutilspy.io.run_command(command)
            #check = subprocess.run('fastq-dump -O ' + outdir + ' --split-files ' + run_sra + " &", shell = True)
            if check.returncode != 0:
                raise ProcessError("\tRun {} could not be processed by fastq-dump".format(run))
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
    check = sutilspy.io.run_command(command)
    
    if check.returncode != 0:
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
    except IntegrityError as error:
        print("\tWARNING: Run(s) in sample {} did not pass integrity check. SKIPPING".format(sample))
        raise ProcessError("\tSample didn't pass check")
    except MissingFileError as error:
        print("\tWARNING: Missing file(s) for run(s) in sample {}. SKIPPING".format(sample))
        raise ProcessError("\tSample didn't pass check")    
    
    # Proceed to fastq-dump
    try:
        run_fastq = fastq_dump_runs(runs,indir,fastqdir,keep)
    except FileNotFoundError as error:
        print("\tERROR: Input directory {} does not exist. TERMINATING".format(indir))
        raise FileNotFoundError("Input directory {} does not exist".format(indir))
    except ProcessError as error:
        print("\tWARNING: Run(s) in sample {} could not be processed by fastq-dump. SKIPPING".format(sample))
        raise ProcessError("\tSample could not be processed  with fastq-dump")
    except MissingFileError as error:
        print("\tWARNING:Run(s) file(s) for sample {} missing".format(sample))
        raise ProcessError("Run(s) file(s) for sample {} missing".format(sample))
    
    # Proceed to concatenate
    try:
        concatenated_files = concatenate_run(run_fastq, outdir, sample, ".fastq.bz2")
    except ProcessError as error:
        print("\tWARNING. Could not concatenate files from sample {}. SKIPPING")
        raise ProcessError("Could not concatenate files from sample {}".format(sample))
    
    return(concatenated_files)

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

def qsub_sample(sample,runs,indir,fastqdir,outdir,logdir,submissionsdir,failedir,keep):
    # Create mappting ffile
    #map = [['Sample','Run']]
    #map.extend(zip(itertools.repeat(sample),runs))
    mapfile = "map." + sample + ".txt"
    sutilspy.io.write_table(outfile = mapfile,
                            rows = zip(itertools.repeat(sample),runs),
                            header = ['Sample','Run'],
                            delimiter = "\t",
                            verbose = True)
    
    # Create qsub file
    submission_file = submissionsdir + "/sra2sample." + sample + ".bash"
    
    with open(submission_file,'w') as fh:
        fh.write("#!/bin/bash\n")
        fh.write("#PBS -N sra2sample." + sample + "\n")
        #fh.write("#PBS -d " + outdir + "\n")
        fh.write("#PBS -o " + logdir + "/sra2sample." + sample + ".log\n")
        fh.write("#PBS -e " + logdir + "/sra2sample." + sample + ".err\n")
        fh.write("#PBS -l mem=1000mb\n")
        
        # Add lines for every run in sample
        failedfile = failedir + "/failed." + sample + ".txt"
        bin = "/home/sur/micropopgen/src/micropopgen/hmp_metadata/sra2sample.py"
        option = ["--indir", indir,
                  "--outdir",outdir,
                  "--fastq_dir", fastqdir,
                  "--map", mapfile, "--run_col", '2',
                  "--sample_col", '1', "--keep_intermediate",
                  "--header", "--failed", failedfile,
                  "--method", "serial"]
        option = " ".join(option)
        
        command = bin + " " + option
        fh.write("module load anaconda\n")
        fh.write(command)
    fh.close()
    os.chmod(submission_file, 0o744)
    
    # submit qsub
    #sutilspy.io.qsub_submissions(submission_file,logdir)

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
    parser.add_argument("--failed", help = "File to store the failed samples", type = str, default = 'failed.txt')
    parser.add_argument("--method", help = "Whether to process samples serially, or use a parallelized approach",
                        type = str, default = "serial", choices = ['serial','qsub'])
    parser.add_argument("--failed_dir", help = "If method is not serial, where to store files of failed samples",
                        type = str, default = "failed")
    parser.add_argument("--submissions_dir", help = "If method is cluster based. Where to store the submission bash files",
                        type = str, default = "submissions")
    parser.add_argument("--logdir", help = "If method is cluster-based, where to store the logfiles",
                         type = str, default = "logs")
    
    args = parser.parse_args()
    args.sample_col -= 1
    args.run_col -= 1
    
#     
#     indir = "./runs/"
#     outdit = "./samples/"
#     runs_file = "/home/sur/micropopgen/exp/2017/today8/runs_to_download.txt"
#     
    runs_per_sample = sutilspy.io.process_run_list(args.map, args.sample_col, args.run_col, args.header)
    
    # Make cluster output directories
    if args.method == 'qsub':
        if not os.path.isdir(args.failed_dir):
            os.mkdir(args.failed_dir)
        if not os.path.isdir(args.submissions_dir):
            os.mkdir(args.submissions_dir)
        if not os.path.isdir(args.logdir):
            os.mkdir(args.logdir)
    
    failed = []
    for sample in runs_per_sample.keys():
        print("== Processing sample {}".format(sample))
        print(" ".join(runs_per_sample[sample]))
        
        if args.method == 'qsub':
            qsub_sample(sample, runs_per_sample[sample], args.indir, 
                        args.fastq_dir, args.outdir,
                        args.logdir, args.submissions_dir,
                        args.failed_dir, args.keep_intermediate)
        else:
            try:
                files = process_sample(sample, runs_per_sample[sample], args.indir,args.fastq_dir, args.outdir, args.keep_intermediate)
            except FileNotFoundError as error:
                print("==Input directory {} does not exist==".format(args.indir))
                raise FileNotFoundError("ERROR:Input directory {} does not exist".format(args.indir))
            except (MissingFileError,ProcessError, IntegrityError) as error:
                print("\tSkipping sample {}".format(sample))
                failed.append([sample])
    if len(failed) > 0:
        write_table(args.failed,failed)
#/home/sur/micropopgen/src/micropopgen/hmp_metadata/sra2sample.py --indir runs --fastq_dir fastq --outdir samples --map map2.txt --failed failed.txt --header --keep_intermediate --run_col 2 --sample_col 1 --failed failed.txt
#/home/sur/micropopgen/src/micropopgen/hmp_metadata/sra2sample.py --indir /godot/hmp/WGS/runs/ --fastq_dir /godot/hmp/WGS/fastq --outdir /godot/hmp/WGS/samples --map /home/sur/micropopgen/data/hmp_download_records/2017-07-17.runs_to_download.txt --failed failed.txt --header --keep_intermediate --run_col 2 --sample_col 1 --method qsub --failed_dir failed --submissions_dir submissions --logdir logs
        
    