package sra;

use strict;
use Exporter;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

$VERSION     = '0.0-1';
@ISA         = qw(Exporter);
@EXPORT      = ();
@EXPORT_OK   = qw(process_run_list match_sra_files_in_dir sample_of_run);
%EXPORT_TAGS = ( run2sample => [qw(&process_run_list match_sra_files_in_dir sample_of_run)]);

sub match_sra_files_in_dir{
	my ($indir,$runs_ref,$sample_of_run_ref,$samples_ref) = @_;
	
	opendir(DIR, $indir) or die "Can't open $indir ($!)";
	my ($file.@sra_files);
	while($file = readdir DIR){
		next unless (-f "$indir/$file") && ($file =~ /\.sra$/);
		
		
		my @sra_files =  readdir DIR;
			
	}
	close DIR;
	
	print scalar @sra_files, ":@sra_files" . "\n";
	
	# NEXT STAGE FIND SAMPLES TO KEEP
	# DEPURATE SAMPLES TO KEEP ACCORDING TO COMPLETNESS IN RUNS_REF
	# RETURN SAMPLES TO KEEP
	
}

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


sub sample_of_run{
	my ($infile,$sample_col,$run_col,$skip) = @_;
	
	open(my $in,$infile) or die "Can't open $infile ($!)";
	my %sample_of_run;
	my $i = 0;
	while(<$in>){
		next unless $i++ >= $skip;
		chomp;
		my @line = split(/\t/,$_);
		my $sample = $line[$sample_col];
		my $run = $line[$run_col];
		if(exists($sample_of_run{$run})){
			die "ERROR: Repeated run ($run), with two samples ($sample_of_run{$run}, $sample)\n";
		}else{
			$sample_of_run{$sample} = $sample;
		}
	}
	my $nruns = $i - $skip;
	my $nsamples = scalar keys %sample_of_run;
	print "Processed $nruns runs in $nsamples samples.\n";
	return(\%sample_of_run);
}


1;