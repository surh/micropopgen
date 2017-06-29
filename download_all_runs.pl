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
use POSIX qw/ ceil floor/;

my $infile = '';
my $outdir = '';
my $method = 'qsub';
my $run_col = 5;
my $sample_col = 1;
my $skip = 1;
my $split_by = 'sample';
my $logdir = 'logs/';
my $ngroups = 2;

my $opts = GetOptions('infile|i=s' => \$infile,
			'outdir|o=s' => \$outdir,
			'method|m=s' => \$method,
			'split|s=s' => \$split_by);

# Main steps. Convert into function
$sample_col--;
$run_col--;
$logdir =  "$outdir/$logdir/";
my $runs_ref = process_run_list($infile,$sample_col,$run_col,$skip);
my $submission_list_ref = create_submission_sets($runs_ref,$outdir,$split_by,$ngroups);
if($method eq 'qsub'){
	qsub_submissions($submission_list_ref,$logdir);
}elsif($method eq 'bash'){
	run_command("$_ &") foreach @$submission_list_ref;
}else{
	die "Method ($method) not recognized.\n";
}

### Subroutines
sub process_run_list{
	my ($infile,$sample_col,$run_col,$skip) = @_;
	
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


sub create_submission_sets{
	my($runs_ref,$outdir,$split_by,$ngroups) = @_;
	
	# Set working space
	mkdir $outdir unless -d $outdir;

	my $submissions_dir = tempdir( "submissions.XXXX", DIR => $outdir);
	print "Saving submission files to $submissions_dir\n";

	# Extract list of samples and pass it to function that creates single submission file
	my @submission_list;
	if ($split_by eq 'sample'){
		print "== Entering split by sample\n";
		my $sample;
		for $sample (keys %$runs_ref){
			my $run_list_ref = $runs_ref->{$sample};
			my $sample_file = create_single_submission_file($sample,$run_list_ref,
									$submissions_dir,
									$logdir,$outdir);
			push(@submission_list,$sample_file);
		}
	}elsif($split_by eq 'groups'){
		#print "== Entering split by number of groups\n";
		my @samples = keys %$runs_ref;
		my $total_samples = scalar @samples;
		my $samples_per_submission = ceil($total_samples / $ngroups);
		my $sample;
		my $group_list_ref = {};
		my $i = 0;
		my $group_i = 0;
		for $sample (@samples){
			#print "====Sample:$sample\n";
			if (($i % $samples_per_submission) == 0){
				$group_i++;	
				$group_list_ref->{"group.$group_i"} = [@{$runs_ref->{$sample}}];
				#print "$group_i\n";
				my @temp = @{$group_list_ref->{"group.$group_i"}};
				#print "group:" . "@temp" . "\n";
				
			}else{
				push(@{$group_list_ref->{"group.$group_i"}},@{$runs_ref->{$sample}});	
			}
			#print "@{$runs_ref->{$sample}}\n";
			#print "===$i\n";
			$i++;
		}
			
		my $group;
		for $group (keys %$group_list_ref){
			#my @temp = @{$group_list_ref->{$group}};
			#print "$group " . scalar(@temp) . "\n";
		    #print "$group " . "@temp" . "\n";
		    my $sample_file = create_single_submission_file($group,$group_list_ref->{$group},
									$submissions_dir,
									$logdir,$outdir);
			push(@submission_list,$sample_file);
		}
	}else{
		die "Unrecognize split_by ($split_by) \n";
	}

	return(\@submission_list);
}

sub create_single_submission_file{
	my ($name,$run_list_ref,$submissions_dir,$logdir,$outdir) = @_;

	# Create new submission file
	my $submission_file = "$submissions_dir/$name.submission.bash";
	open(SUB,'>',$submission_file) or die "Can't create $submission_file ($!)";
	print SUB "#!/bin/bash\n";
	print SUB "#PBS -N download.$name\n";
	print SUB "#PBS -d $outdir\n";
	print SUB "#PBS -o $logdir/download.$name.log\n";
	print SUB "#PBS -e $logdir/download.$name.err\n";
	print SUB "#PBS -l mem=1000mb\n";
	#print SUB "module load sra-tools/2.8.0\n";

	# Add lines for every run in sample
	my $ascp_command = "ascp -i /godot/hmp/aspera/asperaweb_id_dsa.openssh -k 1 -T -l200m";
	my $sra_prefix = "anonftp\@ftp.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/SRR/";
	
	for (@$run_list_ref){
		#print SUB "fastq-dump $outdir $_ --gzip\n";
		my $run_location = substr($_,0,6) . "/$_/$_.sra";
		print SUB "$ascp_command $sra_prefix/$run_location $outdir\n";
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
	run_command("qsub $_") foreach @$submission_list_ref[0..4];

	print "==========SUBMISSIONS DONE==========\n";
}

sub run_command{
	my ($command) = @_;

	my $status = 0;
	print "Executing:\n>$command\n";
	$status = system($command);
	print "Status=$status\n\n";
	sleep 1;

	return($status);
}
