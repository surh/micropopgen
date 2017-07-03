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

# Load local module
use File::Basename qw(dirname);
use Cwd  qw(abs_path);
use lib dirname(abs_path $0) . '/lib';
use sra qw(:run2sample);

my $mapfile = '';
my $indir = '';
my $outdir = '';
my $run_col = 5;
my $sample_col = 1;
my $skip = 1;
#my $method = 'qsub';
#my $split_by = 'sample';
my $logdir = 'logs/';
#my $ngroups = 2;
my $samples_file = '';

my $opts = GetOptions('map|m=s' => \$mapfile,
			'indir|i=s' => \$indir,
			'samples|s=s' => \$samples_file,
			'outdir|o=s' => \$outdir);

# Read map file
$sample_col--;
$run_col--;
$logdir =  "$outdir/$logdir/";
my $runs_ref = process_run_list($mapfile,$sample_col,$run_col,$skip);
my $sample_of_run_ref = sample_of_run($mapfile,$sample_col,$run_col,$skip);

#print ">" . $sample_of_run_ref->{SRR060370} . "\n";

my $sample_list_ref = read_table($samples_file,0,0,0);

my @samples = match_sra_files_in_dir($indir,$runs_ref,$sample_of_run_ref,$sample_list_ref);

#### SUBROUTINES ####









