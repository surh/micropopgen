#!/usr/bin/env perl
# Copyright (C) 2017 Sur Herrera Paredes

# Reads a table with sample and run columns, and dowloads all the run
# fastq files from SRA.

# Currently it keeps all temporary files

# It only uses qsub for submitting to cluster

# Creates one submission per sample

use strict;
use warnings;
use Getopt::Long;
use File::Temp qw/ tempfile tempdir/;

my $infile = '';
my $outdir = '';
my $method = 'qsub';
my $run_col = 5;
my $sample_col = 1;
my $skip = 1;
my $split_by = 'sample';
my $logdir = 'logs/';

my $opts = GetOptions('infile|i=s' => \$infile,
			'outdir|o=s' => \$outdir);


# Main steps. Convert into function
$sample_col--;
$run_col--;
$logdir =  "$outdir/$logdir/";
my $runs_ref = process_run_list($infile,$sample_col,$run_col);
my $submission_list_ref = create_qsub_submission($runs_ref,$outdir);
if($method eq 'qsub'){
	qsub_submissions($submission_list_ref,$logdir);
}else{
	die "Only method qsub is allowed for the momment.\n";
}

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


sub create_qsub_submission{
	my($runs_ref,$outdir) = @_;
	
	# Set working space
	mkdir $outdir unless -d $outdir;

	my $submissions_dir = tempdir( "submissions.XXXX", DIR => $outdir);
	print "Saving submission files to $submissions_dir\n";

	# Extract list of samples and pass it to function that creates single submission file
	my @submission_list;
	if ($split_by eq 'sample'){
		my $sample;
		for $sample (keys %$runs_ref){
			my $run_list_ref = $runs_ref->{$sample};
			my $sample_file = create_single_submission_file($sample,$run_list_ref,
									$submissions_dir,
									$logdir);
			push(@submission_list,$sample_file);
		}
	}else{
		die "Only split by sample is allowed for now \n";
	}

	return(\@submission_list);
}

sub create_single_submission_file{
	my ($name,$run_list_ref,$submissions_dir,$logdir) = @_;

	# Create new submission file
	my $submission_file = "$submissions_dir/$name.submission.bash";
	open(SUB,'>',$submission_file) or die "Can't create $submission_file ($!)";
	print SUB "#!/bin/bash\n";
	print SUB "#PBS -N download.$name\n";
	print SUB "#PBS -d $outdir\n";
	print SUB "#PBS -o $logdir/download.$name.log\n";
	print SUB "#PBS -e $logdir/download.$name.err\n";
	print SUB "#PBS -l mem=1000mb\n";
	print SUB "module load sra-tools/2.8.0\n";

	# Add lines for every run in sample
	for (@$run_list_ref){
		print SUB "fastq-dump $outdir $_ --gzip\n";
	}
	close SUB;

	# Set to executable
	chmod 0744, $submission_file;

	return($submission_file);
}

sub qsub_submissions{
	my ($submission_list_refi,$logdir) = @_;

	mkdir $logdir unless -d $logdir;

	my $file;
	run_command("qsub $_") foreach @$submission_list_ref[0..9];

	print "==========SUBMISSIONS DONE==========\n";
}

sub run_command{
	my ($command) = @_;

	my $status = 0;
	print "Executing:\n>$command\n";
	#my $status = system($command);
	print "Status=$status\n\n";

	return($status);
}
