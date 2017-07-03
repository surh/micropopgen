package sra;

use strict;
use Exporter;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

# Dependencies
use File::Which;
use File::Temp;

$VERSION     = '0.0-1';
@ISA         = qw(Exporter);
@EXPORT      = ();
@EXPORT_OK   = qw(process_run_list match_sra_files_in_dir sample_of_run read_table
                  run_command check_integrity read_vdb_validate);
%EXPORT_TAGS = ( run2sample => [qw(&process_run_list match_sra_files_in_dir sample_of_run
                                   read_table run_command check_integrity read_vdb_validate)]);

sub check_integrity{
	# USA sra-tools to check the md5sum of .sra files.
	# If passed, convert to fastq
	my ($samples_ref,$indir,$runs_ref) = @_;
	
	print "\n=============================================\n";
	print "> Validating integrity of .sra files.\n";
	my $bin = 'vdb-validate';
	$bin = which($bin);
	# check if command exists
	die "Can't find $bin\n." unless -X $bin;
	
	# Create temporary directory for vdb-validate results
	my $tmpdir = File::Temp::tempdir(TEMPLATE => 'vdbXXXXX', CLEANUP => 0, DIR => $indir);
	print ">>Created $tmpdir to save results from $bin\n";
	
	my ($sample);
	# For each sample
	for $sample (@$samples_ref){
		my ($run,@runs);
		my $status = 0;
		# For every run of that sample
		for $run (@{$runs_ref->{$sample}}){
			my $file = "$indir/" . $run . ".sra";
			
			# Create temporary file for validation results
			my $tmp = File::Temp->new( TEMPLATE => 'tempXXXXX', DIR => $tmpdir, SUFFIX => '.txt', CLEANUP => 0);
			
			my $command = "$bin $file > $tmp";
			print ">$command\n";
			#my $out = run_command($command);
			my $check = read_vdb_validate($tmp);
			last if $check;
			
			#print "Check:$check\n";
			
			#push(@runs,$file);
			
		}
		
		# Check if all files passed, and dump fastq files.
		if($check){
			print ">>WARNING: Sample $sample had at least one run faling validation ($run).\n";
		}else{
			for $run (@{$runs_ref->{$sample}}){
				my $file = "$indir/" . $run . ".sra";
				my $command = "fastq-dump -O $outdir $file"; # Missing split command
				print ">$command\n";
				#my $out = run_command($command);
			}	
		}
		
		
	}
	print "=============================================\n";
}

sub read_vdb_validate{
	# Checks if a file with output from vdb-validate returns a positive result.
	my ($infile) = @_;
	open(IN,$infile) or die "Can't open $infile ($!)";
	my @file = <IN>;
	close IN;
	
#	2017-07-03T21:25:19 vdb-validate.2.8.0 info: Table 'SRR059326.sra' metadata: md5 ok
#	2017-07-03T21:25:23 vdb-validate.2.8.0 info: Column 'QUALITY': checksums ok
#	2017-07-03T21:25:24 vdb-validate.2.8.0 info: Column 'RD_FILTER': checksums ok
#	2017-07-03T21:25:25 vdb-validate.2.8.0 info: Column 'READ': checksums ok
#	2017-07-03T21:25:25 vdb-validate.2.8.0 info: Column 'X': checksums ok
#	2017-07-03T21:25:25 vdb-validate.2.8.0 info: Column 'Y': checksums ok
#	2017-07-03T21:25:26 vdb-validate.2.8.0 info: Table 'SRR059326.sra' is consistent
	
	return 1 if ($file[0] !~ /metadata: md5 ok$/);
	return 1 if ($file[1] !~ /Column 'QUALITY': checksums ok$/);
	return 1 if ($file[2] !~ /Column 'RD_FILTER': checksums ok$/);
	return 1 if ($file[3] !~ /Column 'READ': checksums o$/);
	return 1 if ($file[4] !~ /Column 'X': checksums ok$/);
	return 1 if ($file[5] !~ /Column 'Y': checksums ok$/);
	return 1 if ($file[6] !~ /is consistent$/);
	
	return 0;
}

sub match_sra_files_in_dir{
	# Compare .sra files in directory with expected samples from
	# required samples, and decide which samples to keep.
	my ($indir,$runs_ref,$sample_of_run_ref,$samples_ref) = @_;
	
	print ">Preparing to match run files with samples...\n";
	
	opendir(DIR, $indir) or die "Can't open $indir ($!)";
	my ($file, @sra_files, $sample, $run);
	my (%samples_to_keep);
	while($file = readdir DIR){
		# Remove non .sra files
		next unless (-f "$indir/$file") && ($file =~ /\.sra$/);
		
		# Check if run belongs to requested sample
		$run = $file;
		$run =~ s/\.sra$//;
		$sample = $sample_of_run_ref->{$run};
		#print "\t==$run==\t$file\t$sample\n";
		next unless exists($samples_ref->{$sample});
		
		# Save files tthat we are keeping
		push(@sra_files,$file);
		
		# Save samples to keep
		if(exists($samples_to_keep{$sample})){
			$samples_to_keep{$sample}++;
		}else{
			$samples_to_keep{$sample} = 1;
		}
	}
	close DIR;
	
	print ">>Found " . scalar @sra_files . " .sra files that match required samples\n";
	
	# Check if samples wanted have all runs.
	# IMPORTANT: WE ARE ASSUMING THAT THE SAME NUMBER OF RUNS AS EXPECTED
	# MEANS THAT ALL THE RUNS ARE PRESENT.
	my @keep;
	for $sample (keys %samples_to_keep){
		my $expected = scalar @{$runs_ref->{$sample}};
		my $found = $samples_to_keep{$sample};
		#print ">$sample=\t=$found=\t=$expected=\n";
		if($found == $expected){
			push(@keep,$sample);
		}
	}
	print ">>" . scalar @keep . " samples out of " . scalar (keys %samples_to_keep) . " were found to have a complete set of run files.\n"; 
	
	return \@keep;	
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

sub run_command{
	my ($command) = @_;

	my $status = 0;
	print "Executing:\n>$command\n";
	$status = system($command);
	print "Status=$status\n\n";
	sleep 1;

	return($status);
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