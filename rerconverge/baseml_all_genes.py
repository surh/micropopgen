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

from Bio import AlignIO
from skbio import TreeNode
from Bio import Align
from Bio.Phylo.PAML import baseml
import os
import io
import pandas as pd
import numpy as np
import argparse


def process_arguments():
    # Read arguments
    parser_format = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(formatter_class=parser_format)
    required = parser.add_argument_group("Required arguments")

    # Define description
    parser.description = ("Takes alginments from set of genes and coverage "
                          "information and uses baseml to calculate branch "
                          "lengths for each gene")

    # Define required arguments
    required.add_argument("--cov_file",
                          help=("Coverage file. First column must be 'name' "
                                "and every other column corresponds to a "
                                "sample."),
                          required=True, type=str)
    required.add_argument("--aln_dir",
                          help=("Directory with alignments. Files must be "
                                "named <gene>.aln.fasta corresponding to "
                                "gene names in cov_file."),
                          required=True, type=str)
    required.add_argument("--tre_file",
                          help=("Tree file containing all samples desired for "
                                "analysis. Must be newick format."),
                          required=True, type=str)

    # Define other arguments
    parser.add_argument("--outdir", help=("Name of output directory"),
                        type=str,
                        default="output/")
    parser.add_argument("--cov_threshold",
                        help=("Minimum coverage for a gene in a sample to be "
                              "included in baseml analysis"),
                        type=float,
                        default=0.8)
    parser.add_argument("--n_threshold",
                        help=("Minimum number of samples for a gene to be "
                              "included in baseml analysis."),
                        type=int,
                        default=5)
    parser.add_argument("--baseml",
                        help=("Path to baseml binary"),
                        type=str,
                        default='baseml')

    # Read arguments
    print("Reading arguments")
    args = parser.parse_args()

    # Processing goes here if needed

    return args


def subset_aln(infile, outfile, to_keep={},
               informat='fasta', outformat='fasta'):
    """Creates alignment file that contains only sequences in to_keep set."""

    aln = AlignIO.read(infile, informat)

    # Filter aln
    new_aln = Align.MultipleSeqAlignment([])
    for s in aln:
        if s.name in to_keep:
            new_aln.append(s)

    AlignIO.write(alignments=new_aln, handle=outfile, format=outformat)

    return len(new_aln)


def run_baseml(aln_file, tre_file, outdir="output/",
               model=7,
               clock=0,
               tree_format='newick',
               aln_format='fasta',
               baseml_bin='baseml'):
    """Estimate rates from dna alignment and phylogenetic tree
    using PAML's baseml program."""

    # Read original tree and aln
    aln = AlignIO.read(aln_file, aln_format)
    tre = TreeNode.read(tre_file, tree_format)

    # Homogenize tree and aln files
    aln_seqs = {s.id for s in aln}
    tre_tips = {n.name for n in tre.tips()}
    to_keep = aln_seqs & tre_tips
    # Filter aln
    new_aln = Align.MultipleSeqAlignment([])
    for s in aln:
        if s.name in to_keep:
            new_aln.append(s)
    # Filter tre
    new_tre = tre.shear(to_keep)
    new_tre.prune()

    # Write new files (eg. baseml input)
    new_tre_file = os.path.join(outdir, 'tree.newick')
    new_aln_file = os.path.join(outdir, 'aln.phylip')
    TreeNode.write(new_tre, new_tre_file, 'newick')
    # Custom phylip writer
    with open(new_aln_file, 'w') as oh:
        header = ' '.join([str(len(new_aln)), str(len(new_aln[0, :]))]) + "\n"
        oh.write(header)
        for s in new_aln:
            # Need TWO spaces between sequence name and sequence
            line = '  '.join([s.id, str(s.seq)]) + "\n"
            oh.write(line)
    oh.close()

    # run basml
    bml = baseml.Baseml(alignment=new_aln_file, tree=new_tre_file,
                        out_file=os.path.join(outdir, "baseml.out"),
                        working_dir=outdir)
    bml.set_options(model=model, runmode=0, clock=clock)
    res = bml.run(verbose=True, parse=True,
                  command=baseml_bin)

    return(res)


def baseml_all_genes(cov_file, aln_dir, tre_file, outdir="./output/",
                     cov_thres=0.8, n_threshold=5, baseml_bin="baseml"):
    """Run baseml on all genes with only samples above
    certain coverage threshold."""

    # Prepare output directory structure
    os.mkdir(outdir)
    os.mkdir(os.path.join(outdir, "gene_alns"))
    os.mkdir(os.path.join(outdir, "baseml"))
    os.mkdir(os.path.join(outdir, "gene_trees"))

    # Read coverage information
    covs = pd.read_csv(cov_file, sep="\t", dtype={'gene': np.character})
    # covs = covs.set_index('gene').head()
    covs = covs.set_index('gene')

    # Run baseml on every gene
    for g, c in covs.iterrows():
        aln_file = os.path.join(aln_dir, g + '.aln.fasta')

        # Skip if alignment does not exist
        if not os.path.isfile(aln_file):
            continue

        # Find samples to keep
        to_keep = set(c.index[c >= cov_thres])
        # print(g, len(to_keep))
        subset_aln_file = os.path.join(outdir, "gene_alns", g + '.aln.fasta')
        n_samples = subset_aln(infile=aln_file, outfile=subset_aln_file,
                               to_keep=to_keep)

        if n_samples < n_threshold:
            continue

        baseml_dir = os.path.join(outdir, "baseml", g)
        os.mkdir(baseml_dir)
        res = run_baseml(aln_file=subset_aln_file, tre_file=tre_file,
                         outdir=baseml_dir,
                         baseml_bin=baseml_bin)

        tre = TreeNode.read(io.StringIO(res.get('tree')))
        TreeNode.write(tre, file=os.path.join(outdir,
                                              "gene_trees",
                                              g + ".baseml.tre"))


if __name__ == "__main__":
    args = process_arguments()

    baseml_all_genes(cov_file=args.cov_file,
                     aln_dir=args.aln_dir,
                     tre_file=args.tre_file,
                     outdir=args.outdir,
                     baseml_bin=args.baseml,
                     cov_thres=args.cov_threshold,
                     n_threshold=args.n_threshold)
