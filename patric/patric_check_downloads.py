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
    parser.add_argument("--clean", help=("If included, remove failed "
                                         "folders"),
                        action="store_true",
                        default=False)
    parser.add_argument("--features", help=("If included, check fetures.tab "
                                            "file in genome dirs, and confirm "
                                            "that features are consistent with "
                                            ".fna"),
                        action="store_true",
                        default=False)

    parser.add_argument("--gff", help=("If included, check .gff file in genome "
                                       "dirs, and confirm "
                                       "that features are consistent with "
                                       ".fna"),
                        action="store_true",
                        default=False)

    # Read arguments
    print("Reading arguments")
    args = parser.parse_args()

    return args


def check_genomes_dirs(indir, features=False, gff=False):
    """Takes a directory that contains a number of genome subdirectories,
    and checks that every genome subdirtectory has a .fna file

    Taken from patric_download_genomes.py"""

    res = []
    if os.path.isdir(indir):
        # Find all genome subdirs
        specdirs = os.listdir(indir)
        for spec in specdirs:
            # Check if an fna file with the same name as the directory exists
            print("\tChecking genome {}".format(spec))
            fna_filename = os.path.join(indir,
                                        spec,
                                        ''.join([spec, '.fna']))
            if os.path.isfile(fna_filename):
                success = True
            else:
                success = False

            r = [spec, success]
            colnames = ['ID', 'fna']

            # If needed get contig sizes
            if features or gff:
                contig_sizes = get_record_lengths(fna_filename, 'fasta')

            # Check features
            if features:
                feat_files = glob.glob(indir + "/" + spec + "/*.features.tab")
                n_feats = len(feat_files)
                if n_feats > 0:
                    s = []
                    for f in feat_files:
                        s.append(check_patric_features(f, contig_sizes))
                    feat_success = all(s)
                else:
                    feat_success = 'NA'

                r.extend([n_feats, feat_success])
                colnames.extend(["feat_files", 'feats'])

            if gff:
                gff_files = glob.glob(indir + "/" + spec + "/*.gff")
                n_gff = len(gff_files)
                if n_gff > 0:
                    s = []
                    for f in gff_files:
                        s.append(check_patric_gff(f, contig_sizes))
                    gff_success = all(s)
                else:
                    gff_success = 'NA'

                r.extend([n_gff, gff_success])
                colnames.extend(["gff_files", 'gff'])

            res.append(r)

        # colnames don't need to be calculated every time
        res = pd.DataFrame(res, columns=colnames)
    else:
        raise FileNotFoundError("Genome directory does not exist")

    return res


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
