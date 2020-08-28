#!/usr/bin/env python
# Copyright (C) 2020 Sur Herrera Paredes

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

from Bio import SeqIO
import argparse
import os

def process_arguments():
    # Read arguments
    parser_format = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(formatter_class=parser_format)
    required = parser.add_argument_group("Required arguments")

    # Define description
    parser.description = ("Convert TSV files from the "
                          "SNV catalogue of the UHGG into "
                          "VCF files.\n"
                          "This version was developped with "
                          "v1.0 of the catalogue from mgnify "
                          "and VCF version 4.3.")

    # Define required arguments
    required.add_argument("--input",
                          help=("Tab-delimited file (TSV) from the SNV"
                                "catalogue. It must be uncompressed."),
                          required=True, type=str)
    required.add_argument("--genome_fasta",
                          help=("FASTA file of the genome assembly"),
                          required=True, type=str)

    # Define other arguments
    parser.add_argument("--output",
                        help=("Name of the output VCF file"),
                        type=str,
                        default="snvs.vcf")
    parser.add_argument("--include_genomes",
                        help=("Flag indicating to include genomes "
                              "(i.e. sample) columns in the VCF file. "
                              "Otherwise these columns are excluded."),
                        action="store_true",
                        default=False)

    # Read arguments
    print("Reading arguments")
    args = parser.parse_args()

    # Processing goes here if needed
    # Check files exist
    if not os.path.isfile(args.input):
        raise FileNotFoundError("Input tsv file does not exist.")
    if not os.path.isfile(args.genome_fasta):
        raise FileNotFoundError("Genome fasta file does not exist.")

    return args


def contig_lengths(fasta_file):
    """Get contig lengths. The function simply produces a dictionary
    with the length of each record in a fasta file.

    :param fasta_file: Path to fasta file, defaults to None.
    :type fasta_fiile: str

    :return: A dictionary index by record IDs and with the record
    lengths as values.
    :rtype: A dictionary.
    """

    Ctgs = dict()
    for record in SeqIO.parse(assembly_fasta_file, 'fasta'):
        Ctgs[record.id] = len(record.seq)

    return Ctgs


def contig_field_format_vcf(fasta_file):
    """Creates contig field format headers for VCF (v4.3) files by
    assuming that each record in a FASTA file is a contig. It also returns
    the ID of the reference genome based on UHGG conventions (v1.0)

    :param fasta_file: Path to fasta file, defaults to None.
    :type fasta_fiile: str

    :return: A list with
        - Ctr_strs (:py:class:`list`) - VCF contig field headers as :py:class:`str`.
        - ref_id  (:py:class:`str`) - Reference genome ID according to UHGG conventions (v1.0).
    :rtype: A list.

    """

    Ctgs = contig_lengths(fasta_file)
    Ctg_strs = []
    for ctg in Ctgs:
        ctg_str = ''.join(['##contig=<ID=', ctg,
                   ',length=', str(Ctgs[ctg]),
                   '>'])
        Ctg_strs.append(ctg_str)

    # Get reference genome ID
    ref_id = '_'.join(ctg.split("_")[0:2])

    return Ctg_strs, ref_id


def tsv2vcf(snv_file, genome_fasta, outfile='snvs.vcf', include_genomes=False):
    """Convert TSV file from UHGG SNV catalogue into a VCF file.
    Assumes conventions from v1.0 of the catalogue

    :param snv_file: Path to SNV catalogue TSV file, defaults to None.
    :type snv_file: str
    :param genome_fasta: Path to FASTA file with genome assembly.
    :type genome_fasta: str
    :param outfile: Path for output VCF file, defaults to 'snvs.vcf'.
    :type outfile: str, optional
    :param include_genomes: Should genomes (i.e. VCF samples) be included
    in the output file, defaults to False.
    :type include_genomes: bool, optional

    :return: Number of SNVs written to outfile.
    :rtype: int
    """

    with open(snv_file, 'r') as ih, open(outfile, 'w') as oh:
        # Add headers
        oh.write('##fileformat=VCFv4.3' + "\n")
        oh.write('##INFO=<ID=NG,Number=1,Type=Integer,Description="Number of genomes with position">' + "\n")
        oh.write('##INFO=<ID=AF,Number=A,Type=Float,Description="Minor allele frequency">' + "\n")
        oh.write('##FORMAT=<ID=GT,Number=1,Type=Integer,Description="Genotype">' + "\n")
        # Contig definitions
        ctg_strs, ref_id = contig_field_format_vcf(genome_fasta)
        for ctg_str in ctg_strs:
            oh.write(ctg_str + "\n")

        # Read headers (and genome_names) from tsv file
        header = ih.readline()
        genome_names = header.rstrip().split("\t")[4:]
        # Add reference genome to genome_names
        genome_names.append(ref_id)

        # Create VCF header
        vcf_header = ['#CHROM', 'POS',
                  'ID', 'REF', 'ALT',
                  'QUAL', 'FILTER', 'INFO',
                 'FORMAT'] + genome_names
        oh.write("\t".join(vcf_header) + "\n")

        # Line by line VCF record creation
        i = 0
        for line in ih:
            Line = line.rstrip().split("\t")

            # Get genotypes
            genotypes = Line[4:]
            genotypes.append('0')
            n_missing = genotypes.count('255')
            n_present = len(genotypes) - n_missing
            n_minor = genotypes.count('1')
            maf = n_minor / n_present

            # Prepare vcf record
            info = ';'.join(['NG=' + str(n_present), 'AF=' + str(maf)])
            vcf_genotypes = ['.' if g == '255' else g for g in genotypes]
            vcf_line = [Line[0], Line[1],
                       '.', Line[2], Line[3],
                       '.', 'PASS',
                       info, 'GT']
            vcf_line.extend(vcf_genotypes)

            # write vcf record
            oh.write("\t".join(vcf_line) + "\n")
            i = i + 1
    oh.close()
    ih.close()

    return i


if __name__ == "__main__":
    args = process_arguments()

    n = tsv2vcf(snv_file=args.input,
                genome_fasta=args.genome_fasta,
                outfile=args.outfile,
                include_genomes=args.include_genomes)
    print("Finished writing {} SNVs".format(str(n)))
