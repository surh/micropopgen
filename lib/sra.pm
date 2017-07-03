package sra;

use strict;
use Exporter;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

$VERSION     = '0.0-1';
@ISA         = qw(Exporter);
@EXPORT      = ();
@EXPORT_OK   = qw(process_run_list match_sra_files_in_dir sample_of_run read_table);
%EXPORT_TAGS = ( run2sample => [qw(&process_run_list match_sra_files_in_dir sample_of_run read_table)]);

sub match_sra_files_in_dir{
	my ($indir,$runs_ref,$sample_of_run_ref,$samples_ref) = @_;
	
	print "Preparing to match run files with samples...\n";
	
	opendir(DIR, $indir) or die "Can't open $indir ($!)";
	my ($file, @sra_files, $sample, $run);
	while($file = readdir DIR){
		# Remove non .sra files
		next unless (-f "$indir/$file") && ($file =~ /\.sra$/);
		
		$run = $file;
		$file =~ s/\.sra$//;
		print ">" . $sample_of_run_ref->{SRR060370} . "\n";
		$sample = $sample_of_run_ref->{$run};
		print "\t==$run==\t$file\t$sample\n";

		push(@sra_files,$file);
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


sub read_table{
	my($infile,$col1,$col2,$skip) = @_;
	open(IN,$infile) or die "Can't open $infile ($!)";
	my $i = 0;
	my %Table;
	while(<IN>){
		chomp;
		next if $i++ < $skip;
		my @line = split(/\t/,$_);
		$Table{$line[$col1]} = $line[$col2];
	}
	close IN;
	
	return (\%Table);
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
			$sample_of_run{$run} = $sample;
			#print "==$run=>$sample==\n";
		}
	}
	my $nruns = $i - $skip;
	my $nsamples = scalar keys %sample_of_run;
	print "Processed $nruns runs in $nsamples samples.\n";
	return(\%sample_of_run);
}


1;