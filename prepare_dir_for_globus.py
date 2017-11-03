#!/usr/bin/env python
# Copyright (C) 2017 Sur Herrera Paredes
import argparse
import os
import sutilspy

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    required = parser.add_argument_group("Required arguments")
    required.add_argument("--indir", help = "Path to directory where fastq files are located",
                          required = True)
    required.add_argument("--outdir", help = "Path to directory where to create links",
                          required = True)
    required.add_argument("--samples", help = 'File containing sample IDs to be linked. Must be one entry per line',
                          required = True)
    parser.add_argument("--overwrite", helpt = "Flag indicating that any existing file must be overwritten",
                        action = "store_true")
    
    args = parser.parse_args()
    
    samples = sutilspy.io.return_column(args.sample,1,header = False)
    print(samples)
    
    