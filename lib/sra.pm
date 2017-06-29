package sra;

use strict;
use Exporter;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

$VERSION     = '0.0-1';
@ISA         = qw(Exporter);
@EXPORT      = ();
@EXPORT_OK   = qw(process_run_list);
%EXPORT_TAGS = ( all => [qw(&process_run_list)]);


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

1;