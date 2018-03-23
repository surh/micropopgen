#!/usr/bin/env python
# Copyright (C) 2018 Sur Herrera Paredes

import fyrd
import argparse
import os


def process_arguments():
    # Read arguments
    parser_format = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(formatter_class=parser_format)
    required = parser.add_argument_group("Required arguments")

    # Define description
    parser.description = ("Uses EMBOSS' transeq to translate all fasta "
                          "sequences in all fasta files in a directory.\n"
                          "Currently only handles standard genetic code.")

    # Define required arguments
    required.add_argument("--indir", help=("Directory where fasta files are"),
                          required=True, type=str)
    required.add_argument("--outdir", help=("Directory where to write output "
                                            "files"),
                          required=True, type=str)

    # Define other arguments
    parser.add_argument("--fasta_suffix", help=("Suffix of fasta files in "
                                                "indir"),
                        type=str, default='.fna')
    parser.add_argument("--out_suffix", help=("Suffix of output files"),
                        type=str, default='.faa')
    parser.add_argument("--transeq", help=("Binary executable of transeq "
                                           "tool of the EMBOSS package"),
                        type=str, default='transeq')
    parser.add_argument("--logs", help=("Directory where to write logfiles "
                                        "from fyrd."),
                        type=str, default='logs/')
    parser.add_argument("--scripts", help=("Directory where to write scripts "
                                           "from fyrd."),
                        type=str, default='scripts/')
    parser.add_argument("--maxjobs", help=("Maximum number of fyrd jobs "
                                           "running simultaneously"),
                        type=int, default='500')
    parser.add_argument("--queue", help=("Queue (partition) to use for "
                                         "submitting jobs"),
                        type=str)
    parser.add_argument("--memory", help=("Amunt of memory to reserve per "
                                          "job"),
                        type=str, default="300mb")
    parser.add_argument("--time", help=("Amunt of time to reserve per "
                                        "job"),
                        type=str, default="1:00:00")

    # Read arguments
    print("Reading arguments")
    args = parser.parse_args()

    # Processing goes here if needed
    # Check if exectuable exists
    if which(args.transeq) is None:
        raise FileNotFoundError("Executable for transeq not found")
    else:
        args.transeq = which(args.transeq)

    # print(args.transeq)

    return args


def which(program):
    """Check if executable exists"""

    # From:
    # https://stackoverflow.com/questions/377017/test-if-executable-exists-in-python
    # Under MIT license

    def is_exe(fpath):
        """Check if path is executable"""

        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None


def strip_right(text, suffix):
    # tip from http://stackoverflow.com/questions/1038824
    if not text.endswith(suffix):
        return text
    # else
    return text[:len(text)-len(suffix)]


def transeq_file(filename, args, transeq='transeq',
                 indir='', outdir=''):
    """Use fyrd to run transeq on a given file"""

    # Get basename
    basename = strip_right(filename, args.fasta_suffix)

    # Build transeq filenames
    infile = '/'.join([indir, filename])
    outfile = '/'.join([outdir, basename])
    outfile = ''.join([outfile, args.out_suffix])

    # Build transeq command
    command = ' '.join([transeq,
                        "-sequence", infile,
                        '-outseq', outfile])

    # Build fyrd filenames
    job_name = '.'.join(['transeq', basename])
    print(job_name)

    print("\tCreating fyrd.Job")
    fyrd_job = fyrd.Job(command,
                        runpath=os.getcwd(), outpath=args.logs,
                        scriptpath=args.scripts,
                        clean_files=False, clean_outputs=False,
                        mem=args.memory, name=job_name,
                        outfile=job_name + ".log",
                        errfile=job_name + ".err",
                        partition=args.queue,
                        nodes=1, cores=1, time=args.time)

    # Submit joobs
    print("\tSubmitting job")
    fyrd_job.submit(max_jobs=args.maxjobs)


if __name__ == "__main__":
    args = process_arguments()

    # Get list of fasta files from indir
    fasta_files = os.listdir(args.indir)
    fasta_files = list(filter(lambda f: f.endswith(args.fasta_suffix),
                              fasta_files))
    # print(fasta_files)

    # Make directories
    if os.path.isdir(args.logs):
        raise FileExistsError("Directory for fyrd logs ({}) "
                              "already exists".format([args.logs]))
    else:
        os.mkdir(args.logs)
    if os.path.isdir(args.scripts):
        raise FileExistsError("Directory for fyrd scripts ({}) "
                              "already exists".format([args.scripts]))
    else:
        os.mkdir(args.scripts)
    if not os.path.isdir(args.outdir):
        os.mkdir(args.outdir)

    # Submit jobs
    for f in fasta_files:
        print(f)
        transeq_file(filename=f,  args=args,
                     transeq=args.transeq,
                     indir=args.indir,
                     outdir=args.outdir)
