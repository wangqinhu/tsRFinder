#!/usr/bin/env perl 
#===============================================================================
#
#         File: tsRFinder.pl
#
#  Description: A tool for tsRNA analysis with NGS data
#
#       Author: Qinhu Wang
#               Northwest A&F University
#               wangqinhu@nwafu.edu.cn
#               http://wangqinhu.com
#               https://github.com/wangqinhu/tsRFinder
#
#===============================================================================

use strict;
use warnings;
use Env;
use Getopt::Std;

# Enviroment
my $version = '0.9';
my $tsR_dir = $ENV{"tsR_dir"};
my %option;
my %config;

#-------------------------------------------------------------------------------
#  Main
#-------------------------------------------------------------------------------

init();

# Global variables
my $mode     = $option{m} || $config{"mode"} || "run";
my $label    = $option{l} || $config{"label"};
my $refseq   = $option{g} || $config{"reference_genome"};
my $trna     = $option{t} || $config{"reference_tRNA"};
my $srna     = $option{s} || $config{"sRNA"};
my $adaptor  = $option{a} || $config{"adaptor"};
my $minrl    = $option{n} || $config{"min_read_length"} || 18;
my $maxrl    = $option{x} || $config{"max_read_length"} || 45;
my $fam_thr  = $option{f} || $config{"family_threshold"} || 72;
my $lab_trna = $option{w} || $config{"tRNA_with_label"} || "no";

find_tsRNA();

stop();

#-------------------------------------------------------------------------------
#  Subroutine
#-------------------------------------------------------------------------------

# A protocol to tsRNA prediction
sub find_tsRNA {

	# Create output directory
	create_directory();
	
	# Build reference tRNA dataset
	tRNA_scan();

	# Processing the raw sequencing data
	format_reads();

	# Mapping
	map_srna();

	# tsRNA identification
	define_tsRNA();

	# sRNA/tRNA distribution
	draw_distribution();

	# Cleavage position
	find_cleavage_site();

	# Draw cleavage sites
	draw_cleavage_site();

	# Family
	srna_family();

	# Write report file
	write_report();

}

# Init tsRFinder
sub init {

	# Options
	getopts("m:c:l:g:t:s:a:n:x:f:w:hv", \%option) or die "$!\n" . usage();

	# Defualt mode: run
	unless ($option{m}) {
		$option{m} = "run";
	}

	# Check help
	if ($option{h}) {
		usage();
	}

	# Check version
	if ($option{v}) {
		version();
	}

	# Check enviroment
	if (!defined($tsR_dir)) {
		print "Exit: environment tsR_dir is not correctly set up!\n";
		print "Check doc/manual.pdf to see how to do this.\n";
		exit;
	}
	chomp $tsR_dir;
	$tsR_dir =~ s/\/\s{0,}$//;    # replace extra "/" (and spaces) someone may include in the enviroment
	unless ( -d $tsR_dir ) {
		print "Exit: environment tsR_dir is set, but the directory is not exist!\n";
		exit;
	}

	# Check inputs
	if ( defined($option{c}) ) {
		start_log();
	} else {
		if ( defined($option{s}) ) {
			if ( defined($option{t}) or defined($option{g})) {
				if ( defined($option{l}) ) {
					if ( defined($option{a}) ) {
						start_log();
					} else {
						print "No adaptor sequence specified!";
						usage();
					}
				} else {
					print "No label specified!\n";
					usage();
				}
			} else {
				print "No reference tRNA or genome sequence file specified!\n";
				usage();
			}
		} else {
			print "No sRNA file specified!\n";
			usage();
		}
	}

	# Check configs
	my $config_file = $option{c};
	if ( $config_file ) {
		if ( -e $config_file ) {
			print_log("Parsing configuration file $config_file");
			parse_config($config_file);
		} else {
			print_log("Configuration file $config_file does not exist!");
			exit;
		}
	}

	# Check critical dependencies
	check_install("bowtie-build");
	check_install("bowtie");
	check_install("R");

	# Indicating modes
	if ($option{m} eq "debug") {
		print_log("tsRFinder init success: Running in debug mode");
	} else {
		if ($option{m} eq "run") {
			print_log("tsRFinder init success!");
		} else {
			die "Unknown  mode: $option{m}!\n";
		}
	}

}

# Stop tsRFinder
sub stop {

	stop_log();
	if ($mode ne "debug") {
		clean_data();
	} else {
		print_log("Attention: you are in debug mode, temporary files were keeped!");
	}
	print "Finshed!\n";
	exit;

}

# Clean the temporary files
sub clean_data {

	system("rm -rf $label/_trna*");
	system("rm -rf $label/_raw");
	system("rm -rf $label/_srna_map");
	unless ($trna) {
		system("rm $label/trna.ss");
	}
	system("rm $label/tRNA.fas");
	system("rm *.Rout");
	system("rm BDI.txt infile.txt");
	system("rm *.len");
	system("rm tscs.*");
	system("mv ./tsRFinder.log $label/");

}

# Parse config
sub parse_config {

	my ($file) = @_;

	open (FILE, $file) or die "Cannot open file $file: $!\n";
	while (<FILE>) {
		next if /^\s{0,}\#/;
		if (/\s{0,}(\S+)\s{0,}\:\s{0,}(\S+)\s{0,}\#{0,}/) {
			$config{$1} = $2;
		} elsif (/\s{0,}(\S+)\s{0,}\:\s{0,}\#{0,}/) {
			$config{$1} = undef;
		} else {
			if ($option{m} eq "debug") {
				print_log("Unknown configuration line: $_");
			} else {
				next;
			}
		}
	}
	close FILE;

}

# Create the output directory
sub create_directory {

	if ( $label ~~ ["tsRFinder", "demo", "lib", "doc"] ) {
		print "Your label setting conflict with tsRFinder build-in directory, please change to another label!\n";
		print "Label tsRFinder, demo, lib and doc are not allowed!\n";
		exit;
	}

	if ( -e "$label" ) {
		print "$label exist, remove it? [N/y] ";
		my $answer = <>;
		chomp $answer;
		if ($answer ~~ ["Y", "y", "YES", "yes"] ) {
			print_log("Directory $label will be removed now!");
			system("rm -rf $label");
		} else {
			print_log("Exit: directory $label exists");
			exit;
		}
	}
	system("mkdir $label");

}

# Fetech tRNA from genome sequence or user input
sub tRNA_scan {

	unless (  defined($trna) ) {
	    print_log("No reference tRNA supplied, predicting by tRNAscan-SE ...");
		if ( -e $refseq ) {
			predict_tRNA();
		} else {
			print_log("Exit: No reference genome specified!");
			exit;
		}
	} else {
		print_log("Using tRNA $trna");
		my $status = check_valid_tRNA($trna);
		if ($status == 1) {
			print_log("The tRNA file supplied is valid");
		} else {
			print_log("The tRNA file supplied is NOT valid");
			exit;
		}
	}

}

# Check if the tRNA sequence supplied is valid
sub check_valid_tRNA {

	my ($file) = @_;

	open (FILE, "$file") or die "Cannot open file $file: $!\n";

	# A typical tRNA file should contain at least three lines:
	while (<FILE>) {
		# First, begin with ">" and followed by accession number
		if (/^>\S+/) {
			my $seq = <FILE>;
			# Second, have nucleotide sequences
			if (/[^ATCGatcg]+/) {
				my $str = <FILE>;
				# Third, have structure information
				# Please do NOT use ">" for matching becuase all fasta file have it
				if ($str =~ /[\<|\)]+/) {
					format_tRNA($trna);
					return 1;  # valid;
				} else {
					return 0;  # non-valid
				}
			} else {
				return 0;      # non-valid
			}
		} else {
			return 0;          # non-valid
		}
	}

}

# Format tRNA SS file if a valid tRNA file is supplied
sub format_tRNA {

	my ($file) = @_;

	unless ( -e $file ) {
		print_log("file $file does not exist!");
		exit;
	}

	print_log("Formating tRNA ...");
	open (IN, "$file") or die "Cannot open file $file: $!\n";
	open (TSQ, ">$label/tRNA.fa") or die "Cannot create file tRNA.fa: $!\n";
	open (TSS, ">$label/tRNA.fas") or die "Cannot create file tRNA.fas: $!\n";
	while (<IN>) {
		if (/\(+/) {
			s/\(/\>/g;
			s/\)/\</g;
			print TSS $_;
		} elsif (/\<+/) {
			print TSS $_;
		} else {
			print TSQ $_;
			print TSS $_;
		}
	}
	close TSQ;
	close TSS;
	close IN;

}

# A protocal for tRNA prediction
sub predict_tRNA {
	
	# To build a reference tRNA sequence and its sencondary structure file
	run_tRNAscanSE();

	# Extract tRNA sequence and structure information
	extract_tRNA();

	# To build a unique tRNA sequence and structure reference
	unique_tRNA();

}

# tRNAscanSE control
sub run_tRNAscanSE {

	my $ss_file = "trna.ss";

	print_log("Predicting tRNA ...");

	if ( check_install("tRNAscan-SE") ) {
		my $correct_install = `tRNAscan-SE --version 2>&1`;
		if ($correct_install !~ /^Option/) {
			print_log("tRNAscan-SE is not correctly installed, e.g. PERL5LIB is not correctly set.");
			exit;
		} else {
			if ($mode eq "debug") {
				system("tRNAscan-SE -f $label/$ss_file $refseq");
			} else {
				system("tRNAscan-SE -f $label/$ss_file $refseq 1>/dev/null 2>&1");
			}
		}
	}

}

# Extract reliable tRNA sequence and structure
sub extract_tRNA {

	print_log("Extracting tRNA sequences ... ");
	
	my $ss_file = "trna.ss";
	my $trna_file = "_trna.sq";

	print_log("Processing file $ss_file ...");

	# Hold the ss_file content in $trna_ss		
	open (SS, "$label/$ss_file") or die "Cannot open file $ss_file: $!\n";
	my $trna_ss = undef;
	while (<SS>) {
		$trna_ss .= $_;
	}
	close SS;
	
	# Split the sequence and secondary structure of each tRNA
	open (OUT, ">$label/$trna_file") or die "Cannot open file $trna_file: $!\n";
	my @trna = split /\n\n/, $trna_ss;
	my @trnas = ();
	foreach (@trna) {
		# Skip undet and pseudogene tRNA
		next if /Undet|pseudogene/;
		# Parsing SS to extract id, seq, and str
		if (/Type\:\s+(\S+)\s+Anticodon\:\s+(\S+)/) {
			my $id = $1 . $2;
			my $seq = undef;
			my $str = undef;
			if (/Seq\:\s+(\S+)\nStr\:\s+(\S+)/) {
				$seq = $1;
				$str = $2;
				# Discard the sequence with more than 5 'N'
				next if ($seq =~ /N{5,}/i);
				# Output tRNA sequence and its sencondary structure
				print OUT $id, "\t", $seq, "\n";
				print OUT $str, "\n";
			}
		}
	}
	close OUT;

}

# Generate a non-redundant tRNA sequence and structure reference
sub unique_tRNA {

	print_log("Generating non-redundant tRNA ...");

	system("mkdir $label/_trna");

	# Create triplet files
	my @base = ("A", "T", "C", "G");
	my @code = ();
	my $i = 0;
	for my $a (@base) {
		for my $b (@base) {
			for my $c (@base) {
				$code[$i] = $a.$b.$c;
				system("touch $label/_trna/$code[$i]");
				$i++;
			}
		}
	}

	# Generate summary file without redundancy
	my %count = ();
	print_log("The number of unique triplet:");
	my $trna_file = "_trna.sq"; 
	
	# Foreach triplet
	foreach my $ac (@code) {
		# Prepare to write the sencondary sequences
		open (UNI, ">$label/_trna/$ac") or die "Cannot open file $ac: $!\n";
		# For counting
		$count{$ac} = 0;
		# Read the tRNA file
		open (TRNA, "$label/$trna_file") or die "Cannot open file $trna_file: $!\n";
	    # Store the unique sequences
		my %unique = ();
		# id index of each triplet (ac)
		my $idi = 0;
		# Parse the tRNA file
		while (<TRNA>) {
			# Split id and sequence
			my ($id, $seq) = split /\t/;
			# The next line hold the secondary structure of tRNA
			my $str = <TRNA>;
			# If the id, seq, and str belongs to the current triplet
			if ($ac eq substr($id,3,3)) {
				# Use hash to remove redundancy
				if (!exists $unique{$seq}) {
					$unique{$seq} = 1;
					$count{$ac}++;
					$idi++;
					print UNI ">";
					if ($lab_trna eq "yes") {
						print UNI $label;
					}
					print UNI "tRNA-", $id, $idi, "\n", $seq;
					print UNI "$str";
				}
			}
		}

		close TRNA;
		close UNI;
		print_log("$ac\t$count{$ac}");
	}
	
	# Write tRNA sequence file with sencondary structure
	system("cat $label/_trna/* > $label/tRNA.fas");
	# Write tRNA sequence file
	system("awk 'NR%3' $label/tRNA.fas > $label/tRNA.fa");

}

# Format the input sRNA file to adapt tsRFinder
sub format_reads {

	my ( $file ) = ( $srna );

	unless ( -e $file ) {
		print_log("No such file: $file");
		exit;
	}

	# Detect the file type of sRNA supplied	
	my $file_type = fastq_or_fasta($file);

	# Formating the raw data or clean data
	# For fastq format
	if ($file_type eq "fq") {
		print_log("Processing the raw fastq file ...");
		# Build a clean sRNA reads from raw data
		process_raw($file);
		# Format fasta file required by tsRFinder
		format_fasta("$label/_raw/sRNA.fa", "$label/sRNA.fa");
	# For fasta format
	} elsif ($file_type eq "fa") {
		my $fasta = check_fasta($file);
		if ($fasta eq "valid") {
			system("cp $srna $label/sRNA.fa");
		} else {
			# Format fasta file required by tsRFinder
			format_fasta("$srna", "$label/sRNA.fa");
		}
	# If not a fastq or fasta file
	} else {
		print_log("Exit: unknown filetype detected! $file ");
		exit;
	}

}

# Detect file type: fastq or fasta
sub fastq_or_fasta {

	my ($file) = @_;

	my $head = `head -1 $file`;

	# fastq file begin with "@"
	if ($head =~ /^@/) {
		return "fq";
	# fasta file begin with ">"
	} elsif ($head =~ /^>/) {
		return "fa";
	} else {
		return "unkown";
	}

}

# Check fasta
sub check_fasta {

	my ($file) = @_;

	my $head = `head -1 $file`;

	# Typically, a fasta head carrying both id and num sections is required
	unless ($head =~ /^\>(\S+)[\-|\||\_|\t|\s+](\d+)/) {
		print_log("Unsupport fasta format dectected, see the demo file or using your raw data instead");
		exit;
	}

	if ($head =~ /^\>[A-Z]{1,20}\d{1,20}\_\d+/i) {
		print_log("Valid sRNA file found!");
		return "valid";
	}

	return "unknown";

}

# Define a special fasta file required by tsRFinder
sub format_fasta {

	my ($input, $output) = @_;

	print_log("Formatting fasta reads ...");

	my $head = `head -1 $input`;
	open (IN, $input) or die "Cannot open file $input: $!\n";
	open (OUT, ">$output") or die "Cannot open file $output: $!\n";
	if ($head =~ /^\>[A-Z]{1,20}\d{1,20}[\||\t|\-]\d+/i) {
		while (<IN>) {
			s/[\||\t|\-]/\_/;
			print OUT $_;
		}
	} else {
		while (<IN>){
			# Detect, split id and num by [-, |, _, \s, \t]
			if (/^\>(\S+)[\-|\||\_|\t|\s+](\d+)/) {
				my $id = $1;
				my $num = $2;
				# Usually, small RNA sequencing will generate ~ 10,000,000 reads
				# and 7 characters are enough to represent the unique reads
				my $zlen = 7 - length($id);
				my $str = '0' x $zlen;
				# The head is formated as follow:
				# lab0000001
				# ...
				# lab9999999
				$id = $label . $str . $id . "_" . $num;
				print OUT ">$id\n";
			} else {
				print OUT $_;
			}
		}
	}
	close IN;
	close OUT;

}

# Pipeline for processing the raw sequencing data
sub process_raw {

	my ($file) = @_;

	unless ( -e $file) {
		print_log("No such file: $file");
		exit;
	}
	
	system("mkdir $label/_raw");

	# fastq --> fasta
	print_log("Converting fastq to fasta ...");
	fastq_to_fasta($file);

	# Remove adaptor
	print_log("Clipping adaptor ...");
	unless ( $adaptor =~ /[^ATCGatcg]/) {
		if ( length($adaptor) >= 30 ) {
			print_log("Adaptor is too long, please check it!");
			exit;	
		} else {
			fasta_clipper();
		}
	} else {
		print_log("Unkown adaptor: $adaptor");
		exit;
	}

	# Control reads length, too short reads will have problem in mapping
	# And hard for experimental validation
	# A high quality sequencing also generate very few too short reads
	# We recommend the reads length is between 15-50 nt
	# 18-30 nt is OK but may lose some large tsRNA reads
	print_log("Remove reads length less than 15 or more than 50 ...");
	if ( $minrl >= 15 && $minrl <=50 ) {
		if ( $maxrl >= 15 or $maxrl <=50 ) {
			if ( $minrl <= $maxrl ) {
				fasta_len_filter();
			} else {
				print_log("min_read_length larger than max_read_length, please check your configuration file!\n");
				exit;
			}
		} else {
			print_log("max_read_length less than 15 or more than 50!");
		}
	} else {
		print_log("min_read_length less than 15 or more than 50!");
		exit;
	}
	
	# Collapse the fasta sequence to generate a non-redundant fasta file	
	print_log("Collapsing ...");
	collapse_fasta();

}

# Convert fastq to fasta
sub fastq_to_fasta {

	my ($fastq) = @_;
	my $fasta = "$label/_raw/srna.f";

	my $id = undef;
	my $seq = undef;
	my $discard = undef;
	open (IN, $fastq) or die "Cannot open file $fastq: $!\n";
	open (OUT, ">$fasta") or die "Cannot open file $fasta: $!\n";
	while ($id = <IN>) {
		$id =~ s/^\@/>/;
		$seq = <IN>;
		print OUT $id, $seq;
		$discard = <IN>;
		$discard = <IN>;
	}
	close IN;
	close OUT;

}

# Clip adaptor
sub fasta_clipper {

	my $fasta = "$label/_raw/srna.f";
	my $fastc = "$label/_raw/srna.fc";

	# Covert adaptor to capital
	$adaptor =~ tr/atcgn/ATCGN/;
	$adaptor = substr($adaptor,0,6);

	my $id = undef;
	my $seq = undef;
	open (IN, $fasta) or die "Cannot open file $fasta: $!\n";
	open (OUT, ">$fastc") or die "Cannot open file $fastc: $!\n";
	while ($id = <IN>) {
		$seq = <IN>;
		$seq =~ tr/atcgn/ATCGN/;
		chomp $seq;
		if ($seq =~ /(\S+)$adaptor/) {
			$seq = $1;
		}
		print OUT $id, $seq, "\n";
	}
	close IN;
	close OUT;

}

# Fasta sequence length filter
sub fasta_len_filter {

	my $fastc = "$label/_raw/srna.fc";
	my $fastl = "$label/_raw/srna.fl";
	my $len = $minrl;

	my $id  = undef;
	my $seq = undef;
	open (IN, $fastc) or die "Cannot open file $fastc: $!\n";
	open (OUT, ">$fastl") or die "Cannot open file $fastl: $!\n";
	while ($id = <IN>) {
		$seq = <IN>;
		chomp $seq;
		$seq =~ s/\s+//g;
		if (length($seq) >= $len) {
			print OUT $id, $seq, "\n";
		}
	}
	close IN;
	close OUT;

}

# Collapse fasta
sub collapse_fasta {

	my $fastl = "$label/_raw/srna.fl";
	my $fasta = "$label/_raw/sRNA.fa";

	my %num = ();
	my $id = undef;
	my $seq = undef;
	open (IN, $fastl) or die "Cannot open file $fastl: $!\n";
	while ($id = <IN>) {
		$seq = <IN>;
		$num{$seq}++;
	}
	close IN;

	my $index = 0;
	open (OUT, ">$fasta") or die "Cannot open file $fasta: $!\n";
	foreach my $seq (reverse sort {$num{$a} <=> $num{$b}} keys %num) {
		$index++;
		print OUT ">$index-$num{$seq}\n$seq";
	}
	close OUT;

}

# Begin log file
sub start_log {

	my $pwd =`pwd`;

	chomp $pwd;
	
	if ( -e "$pwd/tsRFinder.log") {
		unlink "tsRFinder.log";
	}

	my $time = `date`;

	print_log("Start Time: $time\nCurrent working directory: $pwd");

}

# Stop log file
sub stop_log {

	my $time = `date`;
	print_log("\nStop Time: $time");

}

# Log function: output in screen and save it to log file
sub print_log {

	my ($msg) = @_;

	# Standard output
	print $msg, "\n";

	# Save to log file
	open (LOG, ">>tsRFinder.log") or die "Open file tsRFinder.log failed: $!\n";
	print LOG $msg, "\n";
	close LOG;

}

# sRNA mapping
sub map_srna {

	my $pwd = `pwd`;
	chomp $pwd;
	my $refseq = "$pwd/$label/tRNA.fa";
	my $rnaseq = "$pwd/$label/sRNA.fa";

	# Directories and setting
	my $dir = "$pwd/$label/_srna_map";
	system("mkdir $dir");
	system("mkdir $dir/map");
	system("mkdir $dir/tmp");
	system("mkdir $dir/img");
	system("mkdir $dir/db");

	my $refseqdb = $dir . "/db/refseq.db";
	my $dir_map = $dir . "/map";
	my $dir_tmp = $dir . "/tmp";
	my $dir_img = $dir . "/img";

	# Building and Mapping
	print_log("Building reference ...");
	if ($mode eq "debug") {
		system("bowtie-build $refseq $refseqdb");
	} else {
		system("bowtie-build $refseq $refseqdb 1>/dev/null 2>&1");
	}
	print_log("Mapping ...");
	chdir "$dir_map";
	if ($mode eq "debug") {
		system("bowtie -v 0 -k 1 -a --best --strata $refseqdb -f $rnaseq --refout -t");
	} else {
		system("bowtie -v 0 -k 1 -a --best --strata $refseqdb -f $rnaseq --refout -t 1>/dev/null 2>&1");
	}
	chdir "$pwd";

	# Get reference size
	open (REF, "$refseq") or die "Cannot open $refseq: $!\n";
	my %refsize = ();
	my $refid = undef;
	while (<REF>) {
		chomp;
		if (/^>(\S+)/) {
			$refid = $1;
		} else {
			s/\s//g;
			$refsize{$refid} += length ($_);
		}
	} 
	close REF;

	# Process map file
	opendir (DM, $dir_map ) or die "Cannot open $dir_map: $!\n";
	print_log("Processing map file ...");
	foreach my $file (sort readdir DM) {
		
		# Process only *.map file
		next unless $file =~ /\.map$/;

		# Parse map
		open (MAP, "$dir_map/$file") or die "Cannot open file $file: $!\n";
		my %map = ();
		my $cid = undef;
		while (<MAP>) {
			my ($id, $strand, $refid, $start, $read, $discard) = split /\t/, $_, 6;
			my ($sid, $num) = split /_/, $id;
			# bdi use 0 - index, so the end = len - 1
			my $end = $start + length($read) - 1;
			$cid = $refid;
			$id = $sid . "|". $strand . "|" . $start;
			$map{$id}{"strand"} = $strand;
			$map{$id}{"start"} = $start;
			$map{$id}{"end"} = $end;
			$map{$id}{"num"} = $num;
		}
		close MAP;
		if (!defined $cid) {
			next;
		}

		# Counting reads
		my @coordplus = ();
		my @coordminus = ();
		for (my $i = 0; $i < $refsize{$cid}; $i++) {
			$coordplus[$i] = 0;
			$coordminus[$i] = 0;
			foreach my $id (keys %map) {
				if ($i >= $map{$id}{"start"} && $i <= $map{$id}{"end"}) {
					if ($map{$id}{"strand"} eq "+") {
						$coordplus[$i] += $map{$id}{"num"};
					} elsif ($map{$id}{"strand"} eq "-") {
						$coordminus[$i] += $map{$id}{"num"};
					} else {
						print "Opps, I do not recognize the map strand, check the map!\n";
					}
				}
			}
		}

		# Build block
		open (CUR, ">$dir_tmp/$cid.txt");
		print CUR "position\tplus\tminus\n";
		for (my $j = 0; $j < $refsize{$cid}; $j++) {
			print CUR $j, "\t", $coordplus[$j], "\t", $coordminus[$j], "\n";
		}

		close CUR;

	}
	closedir DM;

	# Plot pattern  	
	chdir "$dir";
	if ($mode eq "debug") {
		system("R CMD BATCH $tsR_dir/lib/draw_map.r");
	} else {
		system("R CMD BATCH $tsR_dir/lib/draw_map.r 1>/dev/null 2>&1");
	}
	system("mv $dir_img $pwd/$label/images");
	chdir "$pwd";

}

# Identify tsRNA from tRNA map
sub define_tsRNA {

	# File variables 
	my $outputmap  = "$label/tsRNA.tmap";        # Output text map file
	my $matureseq  = "$label/tsRNA.seq";         # Mature tsRNA file
	my $report     = "$label/tsRNA.report.xls";  # Report file
	my $trna_reads = "$label/tRNA.read.fa";      # tRNA reads
	my $refseq     = "$label/tRNA.fas";	         # Reference SS file 

	print_log("Finding tsRNA ...");

	# Get reference trna size, sequence and clover structure
	open (REF, $refseq) or die "Cannot open $refseq: $!\n";

	print_log("\tReading tRNA reference ...");

	# Hold map file information
	my %refsize = ();
	my %seq = ();
	my %clover = ();
	my $refid = undef;

	# Include tRNA sequnce, size and structure
	while (<REF>) {
		chomp;
		if (/^>(\S+)/) {
			$refid = $1;
		} else {
			s/\s//g;
			$seq{$refid} = $_;
			$refsize{$refid} = length ($seq{$refid});
			my $next = <REF>;
			chomp $next;
			$clover{$refid} = $next;
		}
	} 
	close REF;

	print_log("\tParsing tRNA map ...");

	# Parse tRNA map
	my $dir_map = "$label/_srna_map/map";
	opendir (DM, $dir_map) or die "Cannot open $dir_map: $!\n";
	open (OUT, ">$outputmap") or die "Cannot open $outputmap: $!\n";
	open (MAT, ">$matureseq") or die "Cannot open $matureseq: $!\n";
	open (RPT, ">$report") or die "Cannot open file $report:$!\n";

	# Write report file head
	print RPT "#tRNA\ttRNA_num\t";
	print RPT "tsR5_name\ttsR5_id\ttsR5_num\ttsR5_read\t";
	print RPT "tsR3_name\ttsR3_id\ttsR3_num\ttsR3_read\n";

	# Hold tRNA reads
	my %trna_reads = ();

	foreach my $file (sort readdir DM) {

		# Process only *.map file
		next unless $file =~ /\.map$/;

		print_log("\tParsing $file ...");

		# Read in mapping information
		open (MAP, "$dir_map/$file") or die "Cannot open file $file: $!\n";
		my @map = ();
		my $i = 0;
		my $refid = undef;
		while (<MAP>) {
			chomp;
			my ($id, $strand, $rid, $start, $read, $discard) = split /\t/, $_, 6;
			# Skip minus strand reads
			next if ($strand eq "-");
			my ($sid, $num) = split /_/, $id;
			my $len = length($read);
			# used for calculation the right "*" or "-" length, no need to minus 1
			my $end = $start + $len;
			$map[$i++] = {	"start"	=> $start,
							"end" 	=> $end,
							"rid" 	=> $rid,
							"num"	=> $num,
							"len"	=> $len,
							"sid" 	=> $sid,
							"read" 	=> $read };
			$refid = $rid;
		}
		close MAP;

		next unless (defined ($refid));
		
		# Extract tRNA reads / tRNA abundance
		my %trna_num = ();
		foreach my $map (@map) {
			my $id = $map->{"sid"} . "_" . $map->{"num"};
			$trna_reads{$id} = $map->{"read"};
			$trna_num{$map->{"rid"}} += $map->{"num"};
		}

		# Caculate BDI
		#
		#                   X-Xmin
		# BDI = floor ( 9 ----------- )
		#                  Xmax-Xmin
		#
		system("cp $dir_map/../tmp/$refid.txt infile.txt");
		if ($mode eq "debug") {
			system("R CMD BATCH $tsR_dir/lib/bdi.r");
		} else {
			system("R CMD BATCH $tsR_dir/lib/bdi.r 1>/dev/null 2>&1");
		}
		open (BDI, "BDI.txt") or die  "Cannot open file BDI.txt!\n";
		my @bdi = ();
		while (<BDI>) {
			chomp;
			push @bdi, $_;
		}
		close BDI;

		# Parse tRNA struture
		my ($arm5, $loop) = ();
		if ( $clover{$refid} =~ 
			m/	(^[\>|\.]+\>\.+\<[\<|\.]+		# first loop
				  [\>|\.]+\>)(\.+)\<[\<|\.]+	# second loop
				  [\>|\.]+\>\.+\<[\<|\.]+		# third loop, varloop or the last loop
			/x) {	
			$arm5 = $1;
			$loop = $2;
		} else {
			print_log("Warning: no clover structure found in $refid!");
		}

		# Find small RNA
		# Definition: top 1 reads in each arms
		# Divide reads
		my @map5 = ();
		my @map3 = ();
		foreach my $map (@map) {
			if ($map->{"end"} lt length ($arm5 . $loop)) {
				push @map5, $map;
			} elsif ($map->{"start"} gt length($arm5)) {
				push @map3, $map;
			} else {
				# loop reads not belong to tsRNA
				next;
			}
		}

		# Escape if no tsRNA reads found (loop reads only)
		next if ( @map5 + @map3 < 1);

		# Search 5' tsRNA
		my @mature5 = sort by_tmap_num @map5;
		my ($srna5) = $mature5[0];

		# Search 3' tsRNA
		my @mature3 = sort by_tmap_num @map3;
		my ($srna3) = $mature3[0];

		# Output basic information
		print OUT "<tmap", "\n";
		print OUT $seq{$refid}, "\t$refid\n";
		print OUT $clover{$refid}, "\n";

		# Output small RNA region
		# The 5' tsRNA
		unless (!defined $srna5) {
			print OUT format_tsRNA($srna5, \%refsize), "\n";
			print MAT $srna5->{"rid"}, "-5p\t", $srna5->{"sid"}, "_", $srna5->{"num"}, "\t", $srna5->{"read"}, "\n";
		}

		# The 3' tsRNA
		unless (!defined $srna3) {
			print OUT format_tsRNA($srna3, \%refsize), "\n";
			print MAT $srna3->{"rid"}, "-3p\t", $srna3->{"sid"}, "_", $srna3->{"num"}, "\t", $srna3->{"read"}, "\n";
		}

		# Output BDI
		print OUT @bdi, "\tBDI\n";

		# Output map
		foreach my $map (sort by_tmap @map) {
			my $left = $map->{"start"};
			my $right = $refsize{$refid} - $map->{"end"};
			my $lstr = "-" x $left;
			my $rstr = "-" x $right;
			print OUT $lstr, $map->{"read"}, $rstr, "\t", $map->{"sid"}, "\t", $map->{"num"}, "\n";
		}

		print OUT "/>\n";

		
		# Write report file
		# The tRNA and its expression
		print RPT $refid, "\t", $trna_num{$refid}, "\t";
		# The 5' tsRNA name, id, seq and expression
		unless (!defined $srna5) {
			print RPT $srna5->{"rid"}, "-5p\t", $srna5->{"sid"}, "\t", $srna5->{"num"}, "\t", $srna5->{"read"}, "\t";
		} else {
			print RPT "NA\tNA\t0\tNA\t";
		}
		# The 3' tsRNA name, id, seq and expression
		unless (!defined $srna3) {
			print RPT $srna3->{"rid"}, "-3p\t", $srna3->{"sid"}, "\t", $srna3->{"num"}, "\t", $srna3->{"read"}, "\n";
		} else {
			print RPT "NA\tNA\t0\tNA\n";
		}

	}

	close MAT;
	close OUT;
	close RPT;
	closedir DM;

	# Output tRNA reads
	open (READ, ">$trna_reads") or die "Cannot open file $trna_reads:$!\n";
	foreach my $id (sort keys %trna_reads) {
		print READ ">", $id, "\n", $trna_reads{$id}, "\n";
	}
	close READ;

}

# Sorting reads by start position on reference
sub by_tmap {

	# sort with start
	$a->{"start"} <=> $b->{"start"}
	# then by num
	or $b->{"num"} <=> $a->{"num"}
	# if cannot sort, by sid
	or $a->{"sid"} cmp $b->{"sid"}

}

# Sorting reads by abundance
sub by_tmap_num {

	# the most abundant read first
	$b->{"num"} <=> $a->{"num"}
	# with small start position
	or $a->{"start"} <=> $b->{"start"}
	# if cannot sort, by sid
	or $a->{"sid"} cmp $b->{"sid"}

}

# Format output of tRNA tmap
sub format_tsRNA {

	my ($map, $refsize) = @_;

	return if (!defined $map);
	
	my $left = $map->{"start"};
	my $right = $refsize->{$map->{"rid"}} - $map->{"end"};
	# produce a lot of stars as the read flanking
	my $lstr = "*" x $left;
	my $rstr = "*" x $right;
	my $matureline = $lstr . $map->{"read"} . $rstr . "\t" . $map->{"sid"} . "\t" . $map->{"num"};

	return $matureline;

}

# write report file
sub write_report {

	print_log("\n\n");
	print_log("---------");
	print_log(" SUMMARY ");
	print_log("---------\n");

	# tRNA summary
	my $ntrna = `grep '>' $label/tRNA.fa | wc -l`;
	$ntrna =~ s/\s//g;
	chomp $ntrna;
	print_log("    tRNA seq : $label/tRNA.fa");
	print_log("       Total : $ntrna\n");

	# sRNA summary
	my ($st, $su) = read_stat("$label/sRNA.fa");
	my ($tt, $tu) = read_stat("$label/tRNA.read.fa");
	print_log("  sRNA reads : $label/sRNA.fa");
	print_log("       Total : $st");
	print_log("      Unique : $su\n");

	# tRNA reads summary
	print_log("  tRNA reads : $label/tRNA.read.fa");
	print_log("       Total : $tt");
	print_log("      Unique : $tu\n");

	# tsRNA seq summary
	my $ntsrna = `cat $label/tsRNA.seq | wc -l`;
	$ntsrna =~ s/\s//g;
	chomp $ntsrna;
	my $nr = '!a[$2]++';
	my $uts = `awk '$nr' $label/tsRNA.seq | wc -l`;
	$uts =~ s/\s//g;
	chomp $uts;
	print_log("   tsRNA seq : $label/tsRNA.seq");
	print_log("       Total : $ntsrna");
	print_log("      Unique : $uts\n");

	# tsRNA report
	print_log("tsRNA report : $label/tsRNA.report.xls\n");

	# tsRNA map
	print_log("    text map : $label/tsRNA.tmap\n");
	print_log("  visual map : $label/images\n");

	# sRNA/tRNA distribution
	print_log("distribution : $label/distribution.pdf\n");

	# tRNA cleavage site
	print_log("    cleavage :");
	print_log("      detail : $label/cleavage.txt");
	print_log("     profile : $label/cleavage_profile.pdf\n");

	# tRNA family
	print_log("tsRNA family : $label/tsRNA.fam\n");

	# statistical measurement
	stat_index();

	print_log("---------");

}

# Small RNA reads statistic
sub read_stat {

	# Require fasta file for input
	my ($file) = @_;

	my $a = '$2';
	my $b = '$1';
	my $total  = `grep ">" $file | awk -F _ '{print $a}' | awk '{sum+=$b} END {print sum}'`;
	my $unique = `grep ">" $file | wc -l`;

	chomp $total;
	$unique =~ s/\s//g;

	return ($total, $unique);

}

# Generate small RNA length distribution
sub draw_distribution {

	print_log("Drawing sRNA/tRNA distribution ...");
	srna_len_stat("$label/sRNA.fa");
	system("mv read.len srna.len");
	srna_len_stat("$label/tRNA.read.fa");
	system("mv read.len trna.len");
	if ($mode eq "debug") {
		system("R CMD BATCH $tsR_dir/lib/draw_distribution.r");
	} else {
		system("R CMD BATCH $tsR_dir/lib/draw_distribution.r 1>/dev/null 2>&1");
	}
	system("mv distribution.pdf $label/");

}

# Small RNA length statistics
sub srna_len_stat {

	my ( $file ) = @_;

	open (IN, $file) or die "Cannot open file $file: $!\n";
	my @read = ();

	foreach ($minrl .. $maxrl) {
		$read[$_] = 0;
	}

	my $read = undef;
	while (<IN>) {
		if (/^\S+[\_|\||\-](\S+)/) {
			$read = <IN>;
			chomp $read;
			$read[length($read)] += $1;
		}
	}
	close IN;

	open (DAT, ">read.len") or die "Cannot open file read.len: $!\n";
	foreach ($minrl .. $maxrl) {
		print DAT $_, "\t", $read[$_], "\n";
	}
	close DAT;

}

sub stat_index {

	my ($sm_hash) = statistical_measure("$label/tsRNA.tmap");

	my @sen = ();
	my @spe = ();
    my @acc = ();
	my $i = 0;

	foreach (keys %$sm_hash) {
		my $stat = $sm_hash->{$_};
		my ($tp, $fn, $fp, $tn) = split /\t/, $stat;
		$sen[$i] = sensitivity($tp, $fn);
		$spe[$i] = specificity($tn, $fp);
		$acc[$i] = accuracy($tp, $fn, $fp, $tn);
		$i++;
	}

	my $sen = mean(@sen);
	my $spe = mean(@spe);
	my $acc = mean(@acc);

	$sen = format_percent($sen,4);
	$spe = format_percent($spe,4);
	$acc = format_percent($acc,4);

	print_log("stat. by BDI :");
	print_log(" Sensitivity : $sen");
	print_log(" Specificity : $spe");
	print_log("    Accuracy : $acc\n");

}

# Format stat index
sub format_percent {

	my($num, $dig) = @_;
	my $percent = undef;

	if ( $num > 0 ) {
		$percent = 	int($num * 10**$dig + 0.5 ) / 10**$dig;
	} else {
		$percent = 	int($num * 10**$dig - 0.4 ) / 10**$dig;
	}

	return $percent;

}
# Math: mean
sub mean {

	my @num = @_;

	my $mean = 0;
	foreach (@num) {
		$mean += $_;
	}
	$mean /= @num;

	return $mean;

}

# True postive rate
sub sensitivity {

	my ($tp, $fn) = @_;
	my $sen = $tp / ($tp + $fn);
	return $sen;

}

# True nagtive rate
sub specificity {

	my ($tn, $fp) = @_;
	my $spe = $tn / ($fp + $tn);
	return $spe;

}

# Prediction accuracy
sub accuracy {

	my ($tp, $fn, $fp, $tn) = @_;
	my $acc = ($tp + $tn) / ($tp + $fn + $fp + $tn);
	return $acc;

}

# Caculate ture/false postive/negative
sub statistical_measure {

	my ( $tmap ) = @_;

	my %tmap = ();
	my $id = undef;
	my ($tp, $fn, $fp, $tn) = (0, 0, 0, 0);

	open (IN, $tmap) or die "Cannot open file $tmap: $!\n";
	while (<IN>) {
		if ( /^\<tmap/ ) {
			$id = <IN>;
			chomp $id;
			my ($discard, $id) = split /\t/, $id;
			my $dicard = <IN>;
			my $seq1 = <IN>;
			my $seq2 = undef;
			my $bdi  = undef;
			my $next = <IN>;
			if ( $next =~ /\*/) {
				$seq2 = $next;
				$bdi  = <IN>;
			} else {
				$bdi = $next;
			}
			($seq1, $discard) = split /\t/, $seq1;
			if (defined($seq2)) {
				($seq2, $discard) = split /\t/, $seq2;
			}
			($bdi, $discard) = split /\t/, $bdi;
			for (my $i = 0; $i< length($seq1); $i++) {
				my $cbdi = substr($bdi,$i,1);
				# if bdi > 1
				if ( $cbdi >= 1 ) {
					# if has tsRNA base => true postive
					if (substr($seq1,$i,1) ne "*" ) {
						$tp++;
					} elsif (defined($seq2) && substr($seq2,$i,1) ne "*") {
						$tp++;
					# if has no tsRNA base => false negative
					} else {
						$fn++;
					}
				# if bdi = 0
				} else {
					# if has tsRNA base => false postive
					if (substr($seq1,$i,1) ne "*" ) {
						$fp++;
					} elsif (defined($seq2) && substr($seq2,$i,1) ne "*") {
						$fp++;
					# if has no tsRNA base => true negative
					} else {
						$tn++;
					}
				}

				#reset current bdi
				$cbdi = 0;
			}
			$tmap{$id} = "$tp\t$fn\t$fp\t$tn";
			# reset vars
			($tp, $fn, $fp, $tn) = (0, 0, 0, 0);
		}
	}

	return (\%tmap);

}

# Find cleavage sites
sub find_cleavage_site {

	my ($cp_hash) = cleavage_pattern("$label/tsRNA.tmap");

	print_log("Gathering cleavage information ...");

	open (CLVG, ">$label/cleavage.txt");
	print CLVG "tRNA\ttsR5\ttsR3\n";
	foreach (keys %$cp_hash) {
		print CLVG $_, "\t", $cp_hash->{$_},"\n";
	}
	close CLVG;

}

# Cleavage information
# Definition:
#                Anticodon
#                ---------
#                    |
#                  -----
# NNNNNNNNNNN--N-N-X-X-X-N-N-NNNNNNNNNNN
#             -1 0 1 2 3 4 5   tRNA site
#         >>.  . . . . . . < <
#            |  | | | | | | |
#           -2 -1 0 1 2 3 4 5  clevage site
#
# N - nucleotide
# X - Anticodon nucleotide
sub cleavage_pattern {

	my ( $tmap ) = @_;

	# Distance from anticodon
	my %dac = ();
	my $id  = undef;

	open (IN, $tmap) or die "Cannot open file $tmap: $!\n";
	while (<IN>) {
		if ( /^\<tmap/ ) {

			# read basic information
			my $trnaseq = <IN>;
			chomp $trnaseq;
			($trnaseq, $id) = split /\t/, $trnaseq;
			my $str  = <IN>;
			chomp $str;
			my $seq1 = <IN>;
			my $discard = undef;
			($seq1, $discard) = split /\t/, $seq1;
			my $seq2 = undef;
			my $next = <IN>;
			if ( $next =~ /\*/) {
				$seq2 = $next;
				($seq2, $discard) = split /\t/, $seq2;
			}

			# Begin stat
			# 1   Determine anticodon position
			# 1.1 Capture anticodon
			my $anticodon = undef;
			($discard, $anticodon) = split /\-/, $id;
			$anticodon = substr($anticodon, 3, 3);
			# 1.2 Parse tRNA struture
			my ($arm5, $loop) = ();
			if ( $str =~
				m/	(^[\>|\.]+\>\.+\<[\<|\.]+		# first loop
					  [\>|\.]+\>)(\.+)\<[\<|\.]+	# second loop
					  [\>|\.]+\>\.+\<[\<|\.]+		# third loop, varloop or the last loop
				/x) {
				$arm5 = length($1);
				$loop = length($2);
			}
			# 1.3 Determine the start and end of anticodon position
			my $anticodon_region = substr($trnaseq, $arm5, $loop);
			$anticodon_region =~ /(\S+)$anticodon\S+/;
			my $anticodon_start = $arm5 + length($1);
			my $anticodon_end = $anticodon_start + 3;
			# 2   Caculate 5' tsRNA (end) distance from anticodon (start)
			if ($seq1 =~ /(\*{0,20}[ATCGNatcgn]+)\*/) {
				my $seq1_end = length($1);
				$dac{$id} = $seq1_end - $anticodon_start;
			# Note: if 5' tsRNA not found, seq1 belongs to 3' tsRNA
			} else {
				$dac{$id} = "NA";
				$seq2 = $seq1;
			}
			# To keep consistence, use the same index as used in 5' tsRNA
			# 3   Caculate 3' tsRNA (start) distance from anticodon (start)
			if (defined($seq2) && $seq2 =~ /(^\*{20,})[ATCGNatcgn]+/) {
				my $seq2_start = length($1);
				my $dac = $seq2_start - $anticodon_end + 3;
				$dac{$id} = $dac{$id} . "\t" . $dac;
			} else {
				$dac{$id} = $dac{$id} . "\t" . "NA";
			}
		}
	}

	return \%dac;

}

# Class tsRNA
sub srna_family {

	my $file = "$label/tsRNA.seq";
	my @nwalign = ("nwalign", "nwalign_linux32", "nwalign_linux64");
	my $nwalign = $nwalign[ os_index() ];

	# Family threshold is specified by user, with a default value of 72
	# Although it can be lower, but we strongly recommend you to use a higher value, such as 100.
	# The family threshlod value can be estimated:
	# Score = 4 X Number of matches
	$fam_thr =~ s/\s//g;

	print_log("Classing sRNA ...");

	my %seq = ();
	my @fam = ();

	open (IN, $file) or die "Cannot open file $file: $!\n";
	while (<IN>) {
		chomp;
		my @w = split /\t/;
		$seq{$w[0]} = $w[2];
	}
	close IN;

	# Pairwise distance
	my %pairwise = ();
	my @ids = sort keys %seq;

	# Class status:
	# 0 - not classed
	# 1 - classed
	my %class = ();
	foreach (@ids) {
		$class{$_} = 0;
	}

	my $fam_index = 0;
	for (my $i = 0; $i < @ids; $i++) {
		my $id1 = $ids[$i];
		if ($class{$id1} == 0) {
			$class{$id1} = 1;
			$fam[$fam_index] = $id1;
			for (my $j = $i+1; $j < @ids; $j++ ) {
				my $id2 = $ids[$j];
				my $id = $id1 . "\t" . $ids[$j];
				$pairwise{$id}= `$tsR_dir/lib/$nwalign -s -a $seq{$id1} -b $seq{$id2}`;
				if ($pairwise{$id} > $fam_thr) {
					$class{$id2} = 1;
					$fam[$fam_index] .= "\t$id2";
				}
			}
			$fam_index++;
		} else {
			next;
		}
	}

	open (CLS, ">$label/tsRNA.fam") or die "Cannot open file tsRNA.fam: $!\n";
	foreach (reverse sort by_len @fam) {
		print CLS $_, "\n";
	}
	close CLS;

}

# Sort by length
sub by_len {

	length($a) <=> length($b);

}

# Dectect Operating System
sub os_index {

	my $os = $^O;
	my $index = undef;
	# To deal with the ones who complile nwalign themselves
	my $file = `file -b $tsR_dir/lib/nwalign`;

	if($os  eq  "darwin") {
		# In OS X, check if nwalign had been compiled on a different system
		if ($file =~ /^Mach/) {
			$index = 0;
		} else {
			print_log("Incorrect nwalign detected:");
			print_log("lib/nwalign: $file");
			print_log("Goto $tsR_dir/lib/src and type make before running");
		}
	} elsif($os  eq  "linux") {
		# my $bit=`getconf LONG_BIT`;
		# chomp $bit;
		# some libs used by getconf may not correctly installed
		# use uname instead
		my $bit = `uname -m`;

		# Unless nwalign is complied by user
		unless ($file =~ /^Mach/ ) {
			$index = 0;
		} elsif ($bit =~ "64") {
			$index = 2;
		} else {
			$index = 1;
		}
	} elsif($os  eq  "MSWin32") {
		die "tsRFinder does not support Windows!\n";
	} else {
		die "Unknown Operating System! Try to compile nwalign in $tsR_dir/lib/src yourself.\n";
	}
	return $index;

}

# Draw cleavage site of tRNA
sub draw_cleavage_site {

	print_log("Building cleavage profile ...");

	open (IN, "$label/cleavage.txt") or die "Cannot open file $label/cleavage.txt: $!\n";

	my %tsR5 = ();
	my %tsR3 = ();

	while (<IN>) {
		chomp;
		my @w = split /\t/;
		if ($w[1] =~ /^(\-?\d+)$/) {
			$tsR5{$1}++;
		}
		if ($w[2] =~ /^(\-?\d+)$/) {
			$tsR3{$1}++;
		}
	}

	open (TC5, ">tscs.5") or die "Cannot create file tscs.5: $!\n";
	foreach (sort by_num keys %tsR5) {
		print TC5 $_, "\t", $tsR5{$_}, "\n";
	}
	close TC5;

	open (TC3, ">tscs.3") or die "Cannot create file tscs.3: $!\n";
	foreach (sort by_num keys %tsR3) {
		print TC3 $_, "\t", "-", $tsR3{$_}, "\n";
	}
	close TC3;

	if ($mode eq "debug") {
		system("R CMD BATCH $tsR_dir/lib/draw_cleavage_site.r");
	} else {
		system("R CMD BATCH $tsR_dir/lib/draw_cleavage_site.r 1>/dev/null 2>&1");
	}
	system("mv cleavage_profile.pdf $label/");

}

# Sort by numbers
sub by_num {

	$a <=> $b;

}

# Check dependency
sub check_install {

	my ($app) = @_;

	my $which = `which $app 2>&1`;
	chomp $which;
	my @w = split /\//, $which;
	if (!defined($w[-1]) or $w[-1] ne "$app") {
		print_log("\n\n\nWARNING:\n$app IS NOT INSTALLED!!\n\n\n");
		exit;
	} else {
		return "1";
	}

}

# usage
sub usage {

print <<USAGE;

tsRFinder usage:

    tsRFinder.pl <option>

    -c  Configuration file
    -l  Label
    -g  Reference genomic sequence
    -t  Reference tRNA sequence
    -s  Small RNA sequence
    -a  Adaptor sequence
    -n  Min read length            [defalut 18]
    -x  Max read length            [default 45]
    -f  Small RNA family threshold [default 72]
    -w  tRNA with/without label    [defualt no]
    -m  mode, run/debug            [defualt run]
    -h  Help
    -v  Version

Example:

    tsRFinder.pl -c demo/tsR.conf

USAGE
exit;

}

# version
sub version {

print <<VERSION;

tsRFinder: version $version

To check if this is the latest version, visit:

https://github.com/wangqinhu/tsRFinder

VERSION
exit;

}


=head1 NAME

  tsRFinder - A tool for tsRNA analysis with NGS data

=head1 SYNOPSIS

  tsRFinder.pl <option>
  For example: ./tsRFinder.pl -c demo/tsR.conf
  Type "./tsRFinder.pl -h" to see all the options.
 
=head1 AUTHOR

  Qinhu Wang
  Northwest A&F University
  E-mail: wangqinhu@nwafu.edu.cn
  
=head1 HOME

  https://github.com/wangqinhu/tsRFinder
  
=head1 LICENSE

  The MIT License.
  
=cut
