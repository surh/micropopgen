#!/usr/bin/env python
# Copyright (C) 2017 Sur Herrera Paredes
import sutilspy
import os

if __name__ == "__main__":
    import argparse
    
    # Parse arguments
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    # Required arguments
    required = parser.add_argument_group("Required arguments")
    required.add_argument("--samples", help = "File with samples (one per line) to be processed",
                          type = str, required = True)
    required.add_argument("--indir", help ="Directory where sample files are located",
                          type = str, required = True)
    required.add_argument("--outdir", help = "Directory where to save output",
                          type = str, required = True)
    
    # Optional arguments
    parser.add_argument("--sample_col",help = "Column where sample id is located in --samples",
                        default = 1, type = int)
    parser.add_argument("--method", help = "Method to use for submissions", type = str,
                        default = 'qsub', choices = ['qsub'])
    parser.add_argument("--logdir", help = "If method is cluster-based, where to store the logfiles",
                         type = str, default = "logs")
    
    args = parser.parse_args()
    #args.sample_col -= 1
    
    
    # Read samples
    samples = sutilspy.io.return_column(infile = args.samples,
                                        col = args.sample_col,
                                        separator = '\t',
                                        header = False)
    #print(samples)
    
    pre_commands = []
    
    # Add module dependencies
    pre_commands.append("module load MIDAS/1.2.1")
    pre_commands.append("echo MIDAS database is $MIDAS_DB")
    bin = "run_midas.py"
    
    for sample in samples:
        sample_file_base = args.indir + "/" + sample
        read1 = sample_file_base + "_read1.fastq.bz2"
        read2 = sample_file_base + "_read2.fastq.bz2"
        
        if not os.path.isfile(read1):
            raise FileNotFoundError("File {} not found".format(read1))
        if not os.path.isfile(read2):
            raise FileNotFoundError("File {} not found".format(read2))
        
        midas_command = [bin,"species",args.outdir + "/" + sample,
                         "-1", read1, "-2", read2,"--remove_temp"]
        midas_command = " ".join(midas_command)
        print(midas_command)
        
    
    
    
    