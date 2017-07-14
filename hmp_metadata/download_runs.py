#!/usr/bin/env perl
# Copyright (C) 2017 Sur Herrera Paredes

# Reads a table with sample and run columns, and dowloads all the run
# fastq files from SRA.

# Currently it keeps all temporary files

# It only uses qsub for submitting to cluster.
# It can use direct command line

# Creates a number of submissions by splitting the List
# of runs into n groups

# Adapted from download_all_runs.pl script

# Run if called as script
if __name__ == '__main__':
    import argparse
    
    # Set argument defaults
    infile = ''
    outdir = ''
    method = 'qsub'
    run_col = 5
    sample_col = 1
    #skip = 1
    header = True
    split_by = 'sample'
    logdir = 'logs/'
    ngroups = 2
    
    # Define arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--infile","-i", help = "Input tab delimited file that maps runs (SRR) and samples (SRS)",
                        type = str)
    
