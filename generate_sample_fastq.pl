#!/usr/bin/env perl
# Copyright (C) 2017 Sur Herrera Paredes

# Read a directory looking for SRA files and a mapping file
# of runs per sample. Then will find the complete samples ( and
# maybe intersect with list of samples) and test for integriity
# using sra-tools, convert to fastq and concatenate all reads
# from a sample.

use strict;
use warnings;
use Getopt::Long;
use File::Temp qw/ tempfile tempdir/;
#use POSIX qw/ ceil floor/;

my $mapfile = '';
my $indir = '';
my $outdir = '';
my $run_col = 5;
my $sample_col = 1;
my $skip = 1;
#my $method = 'qsub';
#my $split_by = 'sample';
#my $logdir = 'logs/';
#my $ngroups = 2;

my $opts = GetOptions('map|m=s' => \$mapfile,
			'indir|i=s' => \$indir,
			'outdir|o=s' => \$outdir);


