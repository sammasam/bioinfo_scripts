#!usr/bin/perl
use warnings;
use strict;

###########################################################################
##### get split fastq file and make split fasta file 
### header id = $ARGV[2]

open(INFILE, "$ARGV[0]") or die "error opening fastq file!\n";
open(OUTFILE, ">>$ARGV[1]") or die "error opening output!\n";

#############################################################################

###########################################################################
##### and call global variables

##### variables for putting each read ( id, seq, and qual) into an array
my $header = $ARGV[2];
my $id;
my $seq;
my $space;
my $qual;
my @seq_array = ();

my @read_array = ();

##### variables for hash of each @read_array

my %read_hash = ();
my $key_num = 1;

##### variables for accessing hash

my @array_to_print = ();

my $array_key = 1;

my $read_id;
my $read_seq;

######################################################################################
##### iteratively put each sequence and qual in an array and populate hash with arrays

while(my $line = <INFILE>) {

	#### match sequence header

	if($line =~ m/^$header/) {
		chomp($line);
		$id = $line;
	
	##### move to seq

		my $nextline1 = <INFILE>;
		chomp($nextline1);
		$seq = $nextline1;
		
	##### move to +

		my $nextline2 = <INFILE>;
		chomp($nextline2);
		$space = $nextline2;
		
	##### move to qual

		my $nextline3 = <INFILE>;
		chomp($nextline3);
		$qual = $nextline3;	

	##### put $id, $seq, and $qual into an array

		push(@seq_array, $id, $seq, $qual);
	
	##### move @seq_array to @read_array

		@read_array = @seq_array;
	
	##### clear variables and @seq_array
	
		@seq_array = ();
		$id = "";
		$seq = "";
		$space = "";
		$qual = "";

		}

	##### put array into %read_hash and empty working array @read_array

	if (@read_array) {
		$read_hash{$key_num} = [ @read_array ];
		$key_num ++;
		@read_array = ();
		}
	}

##### hash is populated so close infile

close(INFILE);

##########(test)########## are id, seq, and qual in the hash with numerical keys?
#    foreach $key_num (sort keys %read_hash) {
#    print "$key_num: @{$read_hash{$key_num}}\n";
#	}
################################################################################

##### iterate over each read set in hash and print to fasta file: split.fasta


foreach $key_num (keys %read_hash) {
	@array_to_print = @{ $read_hash{$array_key}};
	$read_id = $array_to_print[0];
	$read_seq = $array_to_print[1];
	
	while($read_id) {
		print OUTFILE ">$read_id\n";				
		print OUTFILE "$read_seq\n";		
		$read_id = "";
		$read_seq = "";
		$array_key ++;
		@array_to_print = ();
		}
	}

close(OUTFILE);

exit;
















