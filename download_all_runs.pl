#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

my $infile = '';
my $outdir = '';
my $method = 'qsub';
my $run_col = 5;
my $sample_col = 1;
my $skip = 1;

my $opts = GetOptions('infile|i=s' => \$infile,
			'outdir|o=s' => \$outdir);

$sample_col--;
$run_col--;
process_run_list($infile,$sample_col,$run_col);

### Subroutines
sub process_run_list{
	my ($infile,$sample_col,$run_col) = @_;
	
	open(my $in,$infile) or die "Can't open $infile ($!)";
	my %runs_per_sample;
	my $i = 0;
	while(<$in>){
		next unless $i++ >= $skip;
		chomp;
		my @line = split(/\t/,$_);
		my $sample = $line[$sample_col];
		my $run = $line[$run_col];
		if(exists($runs_per_sample{$sample})){
			push(@{$runs_per_sample{$sample}}, $run);
		}else{
			$runs_per_sample{$sample} = [$run];
		}
	}
	my $nruns = $i - $skip;
	my $nsamples = scalar keys %runs_per_sample;
	print "Processed $nruns runs in $nsamples samples.\n";
	return(\%runs_per_sample);
}


