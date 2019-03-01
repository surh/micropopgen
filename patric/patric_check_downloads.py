#!/usr/bin/env python
# Copyright (C) 2019 Sur Herrera Paredes

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
import pandas as pd
import shutil
import glob
from Bio import SeqIO
import re


def process_arguments():
    # Read arguments
    parser_format = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(formatter_class=parser_format)
    required = parser.add_argument_group("Required arguments")

    # Define description
    parser.description = ("Function that checks a folder that contains "
                          "genomes downloaded from PATRIC ftp site.")

    # Define required arguments
    required.add_argument("indir", help=("Directory containing patric "
                                         "downloads. Must contain a two-level "
                                         "directory structure where the first "
                                         "level is some arbitrary grouping "
                                         "and the second level contains one "
                                         "genome per directory."),
                          type=str)

    # Define other arguments
    parser.add_argument("--outfile", help=("Name of file for output"),
                        type=str,
                        default="results.txt")
    parser.add_argument("--clean", help=("If included, remove diecrories that "
                                         "do not have an .fna file"),
                        action="store_true",
                        default=False)
    parser.add_argument("--features", help=("If included, check fetures.tab "
                                            "file in genome dirs, and confirm "
                                            "that features are consistent "
                                            "with .fna"),
                        action="store_true",
                        default=False)

    parser.add_argument("--gff", help=("If included, check .gff file in "
                                       "genome dirs, and confirm "
                                       "that features are consistent with "
                                       ".fna"),
                        action="store_true",
                        default=False)

    # Read arguments
    print("Reading arguments")
    args = parser.parse_args()

    return args


def check_patric_genome(fna_file, features_file=[], gff_file=[]):
    """Checks if features defined in .features.tab and gff file are
    consistent with sequences in fna file"""

    if os.path.isfile(fna_file):
        fna_file_exists = True
    else:
        fna_file_exists = False

    # If needed calculate contig sizes
    n_feats = len(features_file)
    n_gff = len(gff_file)
    if fna_file_exists and (n_feats > 0 or n_gff > 0):
        contig_sizes = get_record_lengths(fna_file, 'fasta')
        if len(contig_sizes) > 0:
            fna_file_has_contigs = True
        else:
            fna_file_has_contigs = False

    # Check features if needed
    checked_features = False
    feat_success = 'NA'
    if n_feats > 0 & fna_file_has_contigs:
        checked_features = True
        for f in features_file:
            try:
                feat_success = check_patric_features(f, contig_sizes)
            except:
                feat_success = "FAILED"

            if feat_success is False or feat_success == "FAILED":
                break

    # Check gff if needed
    checked_gffs = False
    gff_success = 'NA'
    if n_gff > 0 & fna_file_has_contigs:
        checked_gffs = True
        for f in gff_file:
            try:
                gff_success = check_patric_gff(f, contig_sizes)
            except:
                gff_success = "FAILED"

            if gff_success is False or gff_success == "FAILED":
                break

    return fna_file_exists, fna_file_has_contigs, n_feats, checked_features,
    feat_success, n_gff, checked_gffs, gff_success


def check_genomes_dirs(indir, features=False, gff=False):
    """Takes a directory that contains a number of genome subdirectories,
    and checks that every genome subdirtectory has a .fna file

    Taken from patric_download_genomes.py"""

    if os.path.isdir(indir):
        # Find all genome subdirs
        specdirs = os.listdir(indir)
        Res = []
        for spec in specdirs:
            # Check if an fna file with the same name as the directory exists
            print("\tChecking genome {}".format(spec))

            # Build filenames
            fna_filename = os.path.join(indir,
                                        spec,
                                        ''.join([spec, '.fna']))

            if features:
                feat_files = glob.glob(indir + "/" + spec + "/*.features.tab")
            else:
                feat_files = []

            if gff:
                gff_files = glob.glob(indir + "/" + spec + "/*.gff")
            else:
                gff_files = []

            # Check genome
            check = check_patric_genome(fna_filename, feat_files, gff_files)
            Res.append([spec] + check)
    else:
        raise FileNotFoundError("Directory doesn't exist ({})".format(indir))

    # Create data frame with results
    colnames = ['ID', 'fna_file_exists', 'fna_file_has_contigs',
                'n_feats', 'checked_features', 'feat_success',
                'n_gff', 'checked_gffs', 'gff_success']
    Res = pd.DataFrame(Res, columns=colnames)

    return Res


def check_patric_gff(file, contig_sizes):
    """Confirming that gffs downloaded from patric have
    only features that fit into the contig sizes and that
    all the features belong to contigs present"""

    gffs = pd.read_csv(file, sep="\t", header=None, comment="#",
                       names=['seqname', 'source', 'feature',
                              'start', 'end', 'score', 'strand',
                              'frame', 'attribute'],
                       dtype={'seqname': str, 'source': str,
                              'feature': str,
                              'start': int, 'end': int,
                              'score': str, 'strand': str,
                              'frame': int, 'attribute': str})

    if gffs.shape[0] == 0:
        # No features is consistent with genome annotations
        return True

    # Process accession name
    # print("\t=>{}".format(file))
    accession = pd.Series([re.sub('\w+\|', '', s) for s in gffs.seqname])
    # print(accession.unique())

    correct = True
    for contig in contig_sizes:
        # print("Checking contig {}".format(contig))
        starts = gffs.start[accession == contig]
        ends = gffs.end[accession == contig]
        if any(starts < 1) or any(ends > contig_sizes[contig]):
            correct = False
            break

    if len(set(accession.unique()) - contig_sizes.keys()):
        correct = False

    return correct


def check_patric_features(file, contig_sizes):
    """Confirming that features.tab files downloaded from
    patric have only features that fit into the contig sizes,
    and that all the features belong to contigs present"""

    feats = pd.read_csv(file, sep="\t",
                        dtype={'genome_id': str, 'genome_name': str,
                               'accession': str, 'annotation': str,
                               'feature_type': str, 'patric_id': str,
                               'refseq_locus_tag': str, 'alt_locus_tag': str,
                               'uniprotkb_accession': str, 'start': int,
                               'end': int, 'strand': str,
                               'na_length': 'float64',
                               'gene': str, 'product': str,
                               'figfam_id': str, 'plfam_id': str,
                               'pgfam_id': str, 'go': str, 'ec': str,
                               'pathway': str})

    if feats.shape[0] == 0:
        # No features is consistent with genome annotations
        return True

    correct = True
    for contig in contig_sizes:
        # print("Checking contig {}".format(contig))
        starts = feats.start[feats.accession == contig]
        ends = feats.end[feats.accession == contig]
        if any(starts < 1) or any(ends > contig_sizes[contig]):
            correct = False
            break

    if len(set(feats.accession.unique()) - contig_sizes.keys()):
        correct = False

    return correct


def get_record_lengths(file, file_type='fasta'):
    """Read sequence file and get the length of each sequence"""

    record_lengths = dict()
    for record in SeqIO.parse(file, file_type):
        record_lengths[record.id] = len(record.seq)

    return record_lengths


if __name__ == "__main__":
    args = process_arguments()

    specdirs = os.listdir(args.indir)

    Res = pd.DataFrame()
    for s in specdirs:
        print("Processing {}".format(s))
        s_dir = os.path.join(args.indir, s)
        res = check_genomes_dirs(s_dir, features=args.features, gff=args.gff)
        res['path'] = s_dir + '/' + res.ID
        Res = Res.append(res)

    # Write results
    print("Writing results table")
    Res.to_csv(args.outfile, sep="\t", index=False)

    if args.clean:
        print("Cleaning")
        Res = Res[Res.fna != True]
        if Res.shape[0] > 0:
            for d in Res.path:
                print("Removing {}".format(d))
                shutil.rmtree(d)
        else:
            print("No directories to clean")
