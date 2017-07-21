#!/usr/bin/env python
# Copyright (C) 2017 Sur Herrera Paredes
import sutilspy

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
    print(samples)
    
    