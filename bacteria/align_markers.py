#!/usr/bin/env python
# Copyright (C) 2018 Sur Herrera Paredes
# This file is based on the AMPHORA
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

import argparse
import os
import fyrd

# from Bio import AlignIO
from Bio.Alphabet import generic_protein, single_letter_alphabet
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
import numpy as np


def filter_alignment(aln, gap_prop=0.99, remove_singletons=True,
                     alphabet=single_letter_alphabet):
    """Function to filter a numpy array that represents an alignment,
    where rows are records and columns are positions. Assumes gaps are
    given by '-'"""

    # Get sequence records
    nseqs = len(aln)
    rec_names = [r.id for r in aln]

    # Convert to numpu array
    a_array = align2array(aln)

    # Prepare index. Positions to be removed will be changed
    index = np.ones(a_array.shape[1], dtype=int)

    # Iterate over columns
    for i in range(a_array.shape[1]):
        c = a_array[:, i]
        # print(c)
        counts = np.unique(c, return_counts=True)

        # Remove constant columns
        if counts[0].shape == (1,):
            index[i] = 0
            continue

        # Count gaps
        ngaps = counts[1][b'-' == counts[0]]
        if ngaps.shape[0] == 0:
            # print("hello")
            ngaps = np.zeros(1, dtype=int)
        # print("ngaps")
        # print(type(ngaps))
        # print("ngaps:", ngaps)
        # print("ngaps/nseqs", ngaps/nseqs)
        # print("ngaps/nseqs > gap_prop", ngaps/nseqs > gap_prop)

        if ngaps / nseqs > gap_prop:
            index[i] = 0
            continue

        if remove_singletons:
            # DO SOMETHING
            # print(counts[1].size)
            # print(counts[1])
            # print(counts[1].min)
            if counts[1].size == 2 and counts[1].min() == 1:
                index[i] = 0
                continue

        # print("===")

    # Use index to slice array
    index = np.array(index, dtype=bool)
    # print(index.sum())
    filtered = a_array[:, index]

    # Convert back to alignent
    new_aln = array2align(arr=filtered, names=rec_names, alphabet=alphabet)

    return new_aln


def align2array(aln):
    """Convert multiple sequence alignment object to numpy array.
    Taken from tutorial."""

    a_array = np.array([list(rec) for rec in aln], np.character)

    return(a_array)


def array2align(arr, names, alphabet):
    """Convert numpy array to multiple sequence alignment.
    Adapted from documentation"""

    records = []

    # Iterate over array rows (i.e. records)
    for i in range(arr.shape[0]):
        seq = ''.join(np.array(arr[i], dtype=str))
        name = names[i]

        # Concatenate sequence records
        records.append(SeqRecord(Seq(seq, alphabet), id=name))

    # Convert to MSA
    new_aln = MultipleSeqAlignment(records)

    return(new_aln)


def process_arguments():
    # Read arguments
    parser_format = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(formatter_class=parser_format)
    required = parser.add_argument_group("Required arguments")

    # Define description
    parser.description = ("Script that takes a directory of marker files "
                          "and makes an alignment per marker file and "
                          "filters the alignment. Based on the AMPHORA "
                          "pipeline")

    # Define required arguments
    required.add_argument("--indir", help=("Directory with the individual "
                                           "marker files"),
                          required=True, type=str)
    required.add_argument("--outdir", help=("Directory where to store the "
                                            "output"),
                          required=True, type=str)

    # Define other arguments
    parser.add_argument("--muscle", help=("Executable of muscle"),
                        type=str,
                        default="muscle")
    parser.add_argument("--marker_suffix", help=("Suffix of files with "
                                                 "sequences"),
                        type=str, default='.faa')
    parser.add_argument("--ignore_markers", help=("A file of markers to "
                                                  "ignore."),
                        type=str, default='')
    parser.add_argument("--logs", help=("Directory for log files"),
                        type=str, default='logs/')
    parser.add_argument("--scripts", help=("Directory for script files"),
                        type=str, default='scripts/')
    parser.add_argument("--maxjobs", help=("Maximum number of jobs at a "
                                           "given time by fyrd."),
                        type=str, default=1000)
    parser.add_argument("--aln_mem", help=("Memory for alignment jobs"),
                        type=str, default="2000mb")
    parser.add_argument("--aln_time", help=("Time for alignment jobs"),
                        type=str, default="2:00:00")
    parser.add_argument("--aln_queue", help=("Qeueue for alignment jobs"),
                        type=str, default='')
    parser.add_argument("--filter_mem", help=("Memory for filter jobs"),
                        type=str, default="1000mb")
    parser.add_argument("--filter_time", help=("Time for filter jobs"),
                        type=str, default="1:00:00")
    parser.add_argument("--filter_queue", help=("Qeueue for filter jobs"),
                        type=str, default='')
    parser.add_argument("--gap_cov", help=("Maximum allowed propotion of "
                                           "gaps when filtering"),
                        type=float, default=0.99)
    parser.add_argument("--remove_singletons", help=("If included singletons "
                                                     "will also be filtered"),
                        action="store_true")

    # Read arguments
    print("Reading arguments")
    args = parser.parse_args()

    # Processing goes here if needed
    # Check if exectuable exists
    if which(args.muscle) is None:
        raise FileNotFoundError("Executable for hmmscan not found")
    else:
        args.muscle = which(args.muscle)

    return args


def concatenate_marker_files(indir, suffix, outdir='./', ignore=[]):
    # Get list of fasta files from indir
    fasta_files = os.listdir(indir)
    fasta_files = list(filter(lambda f: f.endswith(suffix),
                              fasta_files))

    # Get set of markers
    names = [strip_right(f, suffix) for f in fasta_files]
    markers = set([n.split('.').pop() for n in names])
    # print(markers)
    # markers = set(markers)
    # print(markers)

    print("Concatenating files per marker")
    outfiles = []
    for m in markers:
        # Check if it is in list to ignore
        if m in ignore:
            continue

        # Get files from marker
        marker_suffix = '.' + m + suffix
        files_from_marker = list(filter(lambda f: f.endswith(marker_suffix),
                                        fasta_files))
        files_from_marker = [''.join([indir, '/', f])
                             for f in files_from_marker]

        print("====", m, "====")
        # print(files_from_marker)
        # Build command
        outfile = ''.join([m, '.faa'])
        outfiles.append(outfile)
        outfile = ''.join([outdir, '/', outfile])
        command = ' '.join(['cat'] + files_from_marker + ['>', outfile])

        # Run command
        print(command)
        os.system(command)

    return outfiles


def strip_right(text, suffix):
    # tip from http://stackoverflow.com/questions/1038824
    # MIT License
    if not text.endswith(suffix):
        return text
    # else
    return text[:len(text)-len(suffix)]


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


def submit_align_markers(markersdir, args):
    """Use fyrd to submit alignment jobs"""

    # Create directory for concatenated aligned files
    alndir = ''.join([args.outdir, '/aln/'])
    if os.path.isdir(alndir):
        raise FileExistsError("Alignment dir already exists")
    else:
        os.mkdir(alndir)

    catfiles = os.listdir(markersdir)

    res = dict()
    for f in catfiles:
        infile = ''.join([markersdir, '/', f])
        marker = strip_right(f, args.marker_suffix)
        alnfile = ''.join([alndir, '/', marker, '.aln'])
        job_name = ''.join([marker, '.aln'])
        # print(infile)
        # print(alnfile)
        n, o, j = muscle_file(infile=infile, outfile=alnfile,
                              job_name=job_name, outpath=args.logs,
                              scriptpath=args.scripts,
                              partition=args.aln_queue,
                              time=args.aln_time, muscle=args.muscle,
                              memory=args.aln_mem, maxjobs=args.maxjobs)

        res[n] = [o, j]

    return(res)


def muscle_file(infile, outfile, job_name=None,
                outpath='./logs/', scriptpath='./scripts/',
                partition='', time='01:00:00', muscle='muscle',
                memory='2000mb', maxjobs=1000):
    """Use fyrd to call muscle for aligning each set"""

    # Build muscle command
    command = ' '.join([muscle,
                        "-in", infile,
                        "-out", outfile])

    # Build fyrd filenames
    if job_name is None:
        basename = os.path.basename(infile)
        job_name = '.'.join(['muscle', basename])
    print(job_name)

    print("\tCreating fyrd.Job")
    fyrd_job = fyrd.Job(command,
                        runpath=os.getcwd(), outpath=outpath,
                        scriptpath=scriptpath,
                        clean_files=False, clean_outputs=False,
                        mem=memory, name=job_name,
                        outfile=job_name + ".log",
                        errfile=job_name + ".err",
                        partition=partition,
                        nodes=1, cores=1, time=time)

    # Submit joobs
    print("\tSubmitting job")
    fyrd_job.submit(max_jobs=maxjobs)

    return job_name, outfile, fyrd_job


def submit_filter_alignments(alns, args):
    """Use fyrd to submit jobs for filtering alignments"""

    for n, o in alns.items():
        print("n", n)
        print("n", o[0])
        print("n", o[1])

    return 1


if __name__ == "__main__":
    args = process_arguments()

    # Check input files
    if not os.path.isdir(args.indir):
        raise FileNotFoundError("Indir not found")

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

    # Create output directory
    if os.path.isdir(args.outdir):
        raise FileExistsError("Outdir already exists")
    else:
        os.mkdir(args.outdir)

    # Read ignore list
    ignore = []

    # Concatenate fasta per marker
    # Create directory for concatenated files_from_marker
    markersdir = ''.join([args.outdir, '/cat/'])
    if os.path.isdir(markersdir):
        raise FileExistsError("Outdir already exists")
    else:
        os.mkdir(markersdir)
    concatenate_marker_files(indir=args.indir, suffix=args.marker_suffix,
                             outdir=markersdir, ignore=ignore)

    # Align fasta files per marker
    # print(args.muscle)
    alns = submit_align_markers(markersdir=markersdir, args=args)

    # Filter alignment per marker
    filtered = submit_filter_alignments(alns=alns, args=args)

    # Concatenate overall alignment
