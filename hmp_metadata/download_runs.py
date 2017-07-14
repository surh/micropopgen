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


# print(__name__)
# Run if called as script
if __name__ == "__main__":
    import argparse
        
    # Define arguments
    parser = argparse.ArgumentParser()
    required = parser.add_argument_group("Required arguments")
    required.add_argument("--infile","-i", help = "Input tab delimited file that maps runs (SRR) and samples (SRS)",
                        type = str)
    parser.add_argument("--outdir","-o", help = "Output directory", type = str, default = "./", required = True)
    parser.add_argument("--method","-m", help = "method to submit jobs to the server", type = str,
                        default = 'qsub', choices = ['qsub','bash'])
    parser.add_argument("--sample_col", help = "Column number where the sample ID (SRS) is stored", type = int,
                        default = 5)
    parser.add_argument("--header", help = "Flag indicating wether table has headers in the first row",
                        action = "store_true")
    parser.add_argument("--split_by", help = "Method to split runs for download", default = 'sample',
                        choices = ['sample'])
    parser.add_argument('--logdir','-l', help = "Directory where to store log of run", default = 'logs',
                        type = str)
    parser.add_argument('--ngroups', help = "Number of groups to divide the runs into. Equal to the number of jobs that will be submitted",
                        default = 2, type = int)
    args = parser.parse_args()




