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
    parser.add_argument("--which_markers", help=("A file of markers to use. "
                                                 "if nothing is passed, then "
                                                 "all markers will be used"),
                        type=str, default='')

    # Read arguments
    print("Reading arguments")
    args = parser.parse_args()

    # Processing goes here if needed
    # Check if exectuable exists
    if which(args.muscle) is None:
        raise FileNotFoundError("Executable for hmmscan not found")
    else:
        args.hmmscan = which(args.muscle)

    return args


def concatenate_marker_files(indir, suffix, outdir='./'):
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
        # Get files from marker
        marker_suffix = '.' + m + suffix
        files_from_marker = list(filter(lambda f: f.endswith(marker_suffix),
                                        fasta_files))
        files_from_marker = [''.join([indir,'/',f]) for f in files_from_marker]

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


if __name__ == "__main__":
    args = process_arguments()

    # Check input files
    if not os.path.isdir(args.indir):
        raise FileNotFoundError("Indir not found")

    # Create output directory
    if os.path.isdir(args.outdir):
        raise FileExistsError("Outdir already exists")
    else:
        os.mkdir(args.outdir)

    # Concatenate fasta per marker
    concatenate_marker_files(indir=args.indir, suffix=args.marker_suffix)


    # Align fasta files per marker


    # Filter alignment per marker


    # Concatenate overall alignment
