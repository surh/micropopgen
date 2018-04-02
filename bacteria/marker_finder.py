#!/usr/bin/env python
# Copyright (C) 2018 Sur Herrera Paredes
# This file is based on the MarkerScanner.pl script of the AMPHORA
# pipeline by Martin Wu.

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import fyrd
import argparse
import os
from Bio import SearchIO


def process_arguments():
    # Read arguments
    parser_format = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(formatter_class=parser_format)
    required = parser.add_argument_group("Required arguments")

    # Define description
    parser.description = ("Uses HMMER hmmscan search all protein fasta "
                          "sequences in all fasta files in a directory, "
                          "against a database of profiles o marker genes")

    # Define required arguments
    required.add_argument("--indir", help=("Directory where fasta files are"),
                          required=True, type=str)
    required.add_argument("--outdir", help=("Directory where to write output "
                                            "files"),
                          required=True, type=str)
    required.add_argument("--db", help=("Database of hmm profiles for maker "
                                        "sequences"),
                          required=True, type=str)

    # Define other arguments
    parser.add_argument("--fasta_suffix", help=("Suffix of fasta files in "
                                                "indir"),
                        type=str, default='.faa')
    parser.add_argument("--out_suffix", help=("Suffix of output files"),
                        type=str, default='.hmms')
    parser.add_argument("--hmmscan", help=("Binary executable of hmmscan "
                                           "tool of the HMMER package"),
                        type=str, default='hmmscan')
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
    if which(args.hmmscan) is None:
        raise FileNotFoundError("Executable for hmmscan not found")
    else:
        args.hmmscan = which(args.hmmscan)

    return args


def which(program):
    """Check if executable exists. Returns path of executable."""

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
    # MIT License
    if not text.endswith(suffix):
        return text
    # else
    return text[:len(text)-len(suffix)]


def hmmscan_file(filename, db, args, hmmscan='hmmscan',
                 indir='', outdir=''):
    """Use fyrd to run hmmscan on a given file"""

    # Get basename
    basename = strip_right(filename, args.fasta_suffix)

    # Build hmmscan filenames
    infile = '/'.join([indir, filename])
    outfile = '/'.join([outdir, basename])
    outfile = ''.join([outfile, args.out_suffix])

    # Build hmmscan command
    command = ' '.join([hmmscan,
                        "-Z", str(5000),
                        "-E", str(1e-3),
                        db,
                        infile,
                        ">", outfile])

    # Build fyrd filenames
    job_name = '.'.join(['hmmscan', basename])
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

    return outfile, fyrd_job


def get_hmm_hits(hmmfile):
    """Read HMMER files and get hits"""

    hmmsearch = SearchIO.read(hmmfile, 'hmmer3-text')
    print("==Read==")


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

    # Submit hmmscan jobs
    hmm_files = dict()
    for f in fasta_files:
        hmmfile, job = hmmscan_file(filename=f, db=args.db, args=args,
                                    hmmscan=args.hmmscan,
                                    indir=args.indir,
                                    outdir=args.outdir)
        hmm_files['hmmscan_file'] = job
        print(f)

    # Submit hits_job
    for f, j in hmm_files.items():
        job = fyrd.Job(get_hmm_hits(f), depends=j)
        job.submit(max_jobs=args.maxjobs)
