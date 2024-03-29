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
# import subprocess

from Bio import AlignIO
from Bio.Alphabet import generic_protein, single_letter_alphabet
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
import numpy as np


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


def concatenate_alignments(alns, alphabet=single_letter_alphabet, gap='-'):
    """Take a list of multiple sequence alignments and
    concatenate them, fill with gaps where missing sequences."""

    # Get list of species from alignments
    species = []
    species_per_aln = []
    for a in alns:
        specs = [r.id for r in a]
        species.extend(specs)
        species_per_aln.append(specs)

    species = list(set(species))

    # Create empty alignmet
    new_aln = MultipleSeqAlignment([SeqRecord(Seq('', alphabet),
                                              id=s) for s in species])

    # Iterate over each species, re-ordering when neccessary
    for i in range(len(alns)):
        # print("alginment", i)
        specs = species_per_aln[i]
        if specs != species:
            # print("\treordering")
            # new_alns.append(reorder_alignment(aln=alns[i],
            #                                   specs=specs, species=species))
            new_aln = new_aln + reorder_alignment(aln=alns[i], specs=specs,
                                                  species=species,
                                                  alphabet=alphabet,
                                                  gap=gap)
        else:
            # print("matched")
            new_aln = new_aln + alns[i]

    return new_aln


def create_pipeline_outdirs(dir):
    """Create directories for output"""

    # Create directory for concatenated aligned files
    alndir = ''.join([dir, '/aln/'])
    if os.path.isdir(alndir):
        raise FileExistsError("Alignment dir already exists")
    else:
        os.mkdir(alndir)

    # Create output for filtered files
    fildir = ''.join([dir, '/filtered/'])
    if os.path.isdir(fildir):
        raise FileExistsError("Filtered dir already exists")
    else:
        os.mkdir(fildir)

    return alndir, fildir


def concatenate_and_align_markers(indir, suffix, args, outdir='./',
                                  ignore=[]):
    """Single function that calls function that performs
    concatenation and alignment per marker"""

    markers, fasta_files = get_marker_names(indir, suffix, ignore)
    alndir, fildir = create_pipeline_outdirs(outdir)
    catdir = ''.join([args.outdir, '/cat/'])

    print("Processing markers")
    # outfiles = []
    jobs = []
    for m in markers:
        # Check if it is in list to ignore. Redundant but safe
        if m in ignore:
            print("===Skipping", m)
            continue

        job_name = m + '.catalnfil'
        print(job_name)
        print("\tCreating fyrd job")
        job = fyrd.Job(concatenate_and_align_marker, m,
                       {'indir': indir,
                        'catdir': catdir,
                        'alndir': alndir,
                        'fildir': fildir,
                        'suffix': suffix,
                        'args': args},
                       runpath=os.getcwd(),
                       outpath=args.logs,
                       syspaths=[os.path.dirname(__file__)],
                       imports=[('from align_markers import '
                                 'filter_alignment_file, '
                                 'filter_alignment, '
                                 'align2array, '
                                 'array2align, '
                                 'muscle_file, '
                                 'get_marker_names, '
                                 'create_pipeline_outdirs')],
                       scriptpath=args.scripts,
                       clean_files=False,
                       clean_outputs=False,
                       mem=args.mem,
                       name=job_name,
                       outfile=job_name + ".log",
                       errfile=job_name + ".err",
                       partition=args.queue,
                       nodes=1, cores=1,
                       time=args.time)

        print("\tSubmitting job")
        job.submit(max_jobs=args.maxjobs)
        jobs.append(job)

    return jobs, fildir


def concatenate_and_align_marker(m, indir, catdir, alndir,
                                 fildir, suffix, args):
    """Perform all the alignment and concatenation for a single
    marker"""

    print("====", m, "====")
    # Get list of fasta files from indir
    fasta_files = os.listdir(indir)
    fasta_files = list(filter(lambda f: f.endswith(suffix),
                              fasta_files))

    # Get files from marker
    marker_suffix = '.' + m + suffix
    files_from_marker = list(filter(lambda f: f.endswith(marker_suffix),
                                    fasta_files))
    files_from_marker = [''.join([indir, '/', f])
                         for f in files_from_marker]

    print("Concatenating all sequences from marker")

    outfile = ''.join([m, '.faa'])
    catfile = outfile
    outfile = ''.join([catdir, '/', outfile])
    step = 100
    # Need to iterate over subsets of file to prevent error 32512
    # when command is too long
    for i in range(0, len(files_from_marker), step):
        # Build command
        end = min(i + step, len(files_from_marker))
        command = ' '.join(['cat'] + files_from_marker[i:end] + ['>>', outfile])

        # Run command
        print("Executing:\n>{}".format(command))
        # status = subprocess.run(command)
        status = os.system(command)
        print("Status=", status)

    print("Aligning all sequences from marker")
    infile = ''.join([catdir, '/', catfile])
    alnfile = ''.join([alndir, '/', m, '.aln'])
    job_name = ''.join([m, '.aln'])
    print("muscle", infile, alnfile)

    # Build muscle command
    command = ' '.join([args.muscle,
                        "-in", infile,
                        "-out", alnfile])
    basename = os.path.basename(infile)
    job_name = '.'.join(['muscle', basename])
    print(job_name)
    print("Executing:\n>{}".format(command))
    status = os.system(command)
    print("Status=", status)
    if not os.path.isfile(alnfile):
        raise FileNotFoundError("Output alignment file not found")

    # Filter alignment
    print("Filter alignment")
    filfile = ''.join([fildir, '/', m, '.filtered.aln'])
    job_name = m + '.filter'
    print(job_name)
    print("filter", alnfile, filfile)
    filter_alignment_file(alnfile, filfile, gap_prop=args.gap_prop)
    # print("hola")

    return filfile


def concatenate_marker_files(indir, suffix, outdir='./', ignore=[]):
    """Submit system job to concatenate all files from
    one marker"""

    markers, fasta_files = get_marker_names(indir, suffix, ignore)

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
    # print("Iterating")
    for i in range(a_array.shape[1]):
        c = a_array[:, i]
        # print("\t", c)
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


def filter_alignment_file(infile, outfile, gap_prop=0.99,
                          remove_singletons=True,
                          alphabet=single_letter_alphabet,
                          input_format='fasta',
                          output_format='fasta'):
    """Take file with a single alignment, filter it and
    write a new file"""

    print("\treading")
    aln = AlignIO.read(infile, input_format)
    print("\tfiltering")
    filtered = filter_alignment(aln=aln, gap_prop=gap_prop,
                                remove_singletons=remove_singletons,
                                alphabet=alphabet)
    print("\twriting")
    AlignIO.write(filtered, outfile, output_format)

    return(outfile)


def get_marker_names(indir, suffix, ignore=[]):
    """Get list of markers"""

    # Get list of fasta files from indir
    fasta_files = os.listdir(indir)
    fasta_files = list(filter(lambda f: f.endswith(suffix),
                              fasta_files))

    # Get set of markers
    names = [strip_right(f, suffix) for f in fasta_files]
    markers = set([n.split('.').pop() for n in names])

    # Remove markers to ignore
    markers = list(filter(lambda i: i not in ignore, markers))

    return markers, fasta_files


def muscle_file(infile, outfile, mode='fyrd', job_name=None,
                outpath='./logs/', scriptpath='./scripts/',
                partition='', time='01:00:00', muscle='muscle',
                memory='2000mb', maxjobs=1000):
    """Use fyrd to call muscle for aligning each set"""

    # Build muscle command
    command = ' '.join([muscle,
                        "-in", infile,
                        "-out", outfile,
                        "-quiet"])

    # Build fyrd filenames
    if job_name is None:
        basename = os.path.basename(infile)
        job_name = '.'.join(['muscle', basename])
    print(job_name)

    if mode == 'fyrd':
        print("\tCreating fyrd.Job")
        job = fyrd.Job(command,
                       runpath=os.getcwd(),
                       outpath=outpath,
                       scriptpath=scriptpath,
                       clean_files=False,
                       clean_outputs=False,
                       mem=memory,
                       name=job_name,
                       outfile=job_name + ".log",
                       errfile=job_name + ".err",
                       partition=partition,
                       nodes=1, cores=1,
                       time=time)

        # Submit joobs
        print("\tSubmitting job")
        job.submit(max_jobs=maxjobs)
    elif mode == 'bash':
        print("Executing:\n>{}".format(command))
        status = os.system(command)
        print("Status=", status)
        job = status
        if not os.path.isfile(outfile):
            raise FileNotFoundError("Output alignment file not found")
    else:
        raise ValueError("bash or fyrd")

    return job_name, outfile, job


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
    parser.add_argument("--mem", help=("Memory per fyrd job"),
                        type=str, default="10000mb")
    parser.add_argument("--time", help=("Time per fyrd job"),
                        type=str, default="10:00:00")
    parser.add_argument("--queue", help=("Qeueue for fyrd jobs"),
                        type=str, default='')
    parser.add_argument("--cat_mem", help=("Memory per fyrd job"),
                        type=str, default="10000mb")
    parser.add_argument("--cat_time", help=("Time per fyrd job"),
                        type=str, default="10:00:00")
    parser.add_argument("--cat_queue", help=("Qeueue for fyrd jobs"),
                        type=str, default='')
    parser.add_argument("--gap_prop", help=("Maximum allowed propotion of "
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


def reorder_alignment(aln, specs, species,
                      alphabet=single_letter_alphabet, gap='-'):
    """Take an alignment and reorder it acording to species list.
    Add records as gapped if missing"""

    new_aln = []
    missing_seq = ''.join([gap] * aln.get_alignment_length())
    for s in species:
        # Check if species exist in alignment
        try:
            i = specs.index(s)
        except ValueError:
            i = -1
        except:
            raise

        if i >= 0:
            new_aln.append(aln[i])
        elif i == -1:
            new_aln.append(SeqRecord(Seq(missing_seq, alphabet), id=s))

    new_aln = MultipleSeqAlignment(new_aln)

    return(new_aln)


def strip_right(text, suffix):
    # tip from http://stackoverflow.com/questions/1038824
    # MIT License
    if not text.endswith(suffix):
        return text
    # else
    return text[:len(text)-len(suffix)]


def submit_concatenate_alignments(indir, args):
    """Takes a directory of alignment files, reads them
    and creates a fyrd job to concatenate them all"""

    # Create output directory
    catalndir = ''.join([args.outdir, '/cataln/'])
    if os.path.isdir(catalndir):
        raise FileExistsError("Concatenated alignment dir already exists")
    else:
        os.mkdir(catalndir)

    # Read alignments
    alnfiles = os.listdir(indir)
    alns = []
    for f in alnfiles:
        alns.append(AlignIO.read(indir + '/' + f, 'fasta'))

    print("=============CONCATENATING ALIGNMENTS===============")
    job_name = 'cat.alns'
    print(job_name)
    print("\tCreating fyrd.Job")
    job = fyrd.Job(concatenate_alignments, [alns],
                   {'alphabet': single_letter_alphabet,
                    'gap': '-'},
                   runpath=os.getcwd(),
                   outpath=args.logs,
                   syspaths=[os.path.dirname(__file__)],
                   imports=[('from align_markers import '
                             'reorder_alignment')],
                   scriptpath=args.scripts,
                   clean_files=False, clean_outputs=False,
                   mem=args.cat_mem,
                   name=job_name,
                   outfile=job_name + ".log",
                   errfile=job_name + ".err",
                   partition=args.cat_queue,
                   nodes=1, cores=1,
                   time=args.cat_time)

    print("\tSubmitting job")

    res = job.submit(max_jobs=args.maxjobs)
    aln = res.get()

    print("\t Writing output file")
    outfile = catalndir + '/' + 'maker_concatenated_alignment.aln'
    AlignIO.write(aln, outfile, 'fasta')

    print("=============CONCATENATING ALIGNMENTS===============")

    return outfile


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
    ignore = ['mutM']

    # Concatenate fasta per marker
    # Create directory for concatenated files_from_marker
    markersdir = ''.join([args.outdir, '/cat/'])
    if os.path.isdir(markersdir):
        raise FileExistsError("Outdir already exists")
    else:
        os.mkdir(markersdir)

    jobs, fildir = concatenate_and_align_markers(indir=args.indir,
                                                 suffix=args.marker_suffix,
                                                 args=args,
                                                 outdir=args.outdir,
                                                 ignore=ignore)

    # Waiting for jobs to complete
    print("Waiting for jobs to complete")
    for j in jobs:
        j.wait()
    print("Finished alignment and filtering")

    # Concatenate overall alignment
    submit_concatenate_alignments(fildir, args)
