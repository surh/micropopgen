#!/usr/bin/env python
# Copyright (C) 2017 Sur Herrera Paredes
import argparse
import os
import sutilspy
import shutil


def check_sample(r1_file, r2_file, path, extension='.fastq.bz2'):
    """Check whether the files for a given sample are present,
    and if the output directory exists"""

    # Check that files exist
    try:
        if not os.path.isfile(r1_file):
            raise FileNotFoundError
    except:
        print("""ERROR: File {} is missing""".format(r1_file))
        raise

    try:
        if not os.path.isfile(r2_file):
            raise FileNotFoundError
    except:
        print("""File {} is missing""".format(r2_file))
        raise

    # Check that output directory is present
    try:
        print("ISDIR")
        print(path)
        print(os.path.isdir(path))
        if not os.path.isdir(path):
            raise ValueError
    except:
        print("""ERROR: Output directory ({}) missing""".format(path))
        raise


def create_links(sample, indir, outdir, extension='.fastq.bz2',
                 checks=True):
    """Create symbolic links for files from sample"""
    # Create file names
    r1_file = ''.join([indir, '/', sample, '_read1', extension])
    r2_file = ''.join([indir, '/', sample, '_read2', extension])

    # Run checks
    try:
        if checks:
            check_sample(r1_file=r1_file, r2_file=r2_file, path=outdir)
    except:
        print(1)
        print("\tSample {} failed checks".format(sample))
        raise

    # Create links
    try:
        print(2)
        os.symlink(src=r1_file, dst=outdir)
        os.symlink(src=r1_file, dst=outdir)
    except:
        print(("ERROR:Could not create symbolic links "
               "for sample {}").format(sample))
        raise


if __name__ == '__main__':
    print("==============================================================")
    # Read arguments
    parser_format = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(formatter_class=parser_format)
    required = parser.add_argument_group("Required arguments")
    required.add_argument("--indir", help="""Path to directory where fastq
                          files are located""",
                          required=True)
    required.add_argument("--outdir", help="""Path to directory where to
                          create links""",
                          required=True)
    required.add_argument("--samples", help="""File containing sample IDs
                          to be linked. Must be one entry per line""",
                          required=True)
    parser.add_argument("--overwrite", help="""Flag indicating that any
                        existing file must be overwritten""",
                        action="store_true")
    parser.add_argument("--failure", help="""What to do upon failure of
                        a sample. Remove everythin (clear).
                        Keep going (continue)""", default='clear',
                        choices=['clear', 'continue'])
    parser.add_argument("--failed", help="""Filename where to store failed
                        samples""", default="failed.txt")

    print("Reading arguments")
    args = parser.parse_args()

    # Read list of samples
    try:
        print("Reading sample list...")
        samples = sutilspy.io.return_column(args.samples, 1, header=False)
        n = str(len(samples))
        print("{} samples read".format(n))
    except:
        print("ERROR: Failed reading sample list")
        raise

    # Create output directory
    try:
        print("Creating output directory")
        os.mkdir(args.outdir)
    except:
        print(("ERROR: Could not create output "
               "directory ({})").format(args.outdir))
        raise

    # Check every sample
    print("\n\n\t====Processing samples===")
    failed = []
    for s in samples:
        try:
            create_links(sample=s, indir=args.indir,
                         outdir=args.outdir, checks=True)
        except:
            print("Sample {} failed".format(s))
            if args.failure == 'clear':
                print("Cleaning and aborting")
                shutil.rmtree(path=args.outdir)
            elif args.failure == 'continue':
                failed.append([s])
                print("\tContinuing")

    # Write failed samples
    if len(failed) > 0:
        # print(failed)
        sutilspy.io.write_table(outfile=args.failed, rows=failed,
                                header=None, delimiter="\t",
                                verbose=False)
