#!usr/bin/perl
use warnings;
use strict;

###########################################################################
##### extract GIs from BLAST hits file and output total hits ##############
##### for each GID as well as list of query seqs for each GID #############
###########################################################################
##### 1st ARGV[0] = BLAST File ############################################
##### 2rd ARGV[1] = gene name for output files ############################
###########################################################################
###########################################################################
###########################################################################


###########################################################################
##### VARIABLES FOR TWO HASHES ############################################
###########################################################################
##### 1ST HASH FOR COUNT TOTALS FOR EACH GID ##############################
##### 2ND HASH FOR ARRAY OF QUERY IDs FOR EACH GID ########################
###########################################################################

my $gene = $ARGV[1];
my $countlist = "$gene\.blast_counts\.txt";
my $hitlist = "$gene\.all_hits_list\.txt";
my $output = "$gene\.hits_by_gid\.txt";


my $gid = "";
my @line_array = ();
my $total_count = 0;
my $ref_count = 0;

my %count_hash = ();
my $geneid = "";	# # key
my $hit_count = 0;	# # value

my $hitid = "";
my %hit_hash = ();
my $subject = "";
my @hit_array = ();	# # value

my $next_output = 0;
my @count_array = ();
my $gnid = "";
my $gnid_count = 0;
my @hits_array = ();

###########################################################################

open(HITLIST, ">>$hitlist") or die "error opening all hits output!\n";
open(INFILE, "$ARGV[0]") or die "error opening fastq file!\n";

while(<INFILE>) {
	chomp($_);
	$total_count ++;
	@line_array = split(/\t/, $_);
	$gid = $line_array[4];
	$hitid = $line_array[7];
	$hit_count = 0;
	
	if(exists $count_hash{$gid}) {
		$hit_count = $count_hash{$gid};
		$hit_count ++;
		$count_hash{$gid} = $hit_count;
		}

	else {
		$count_hash{$gid} = 1;
		$ref_count ++;
		print HITLIST "$hitid\n";
		}
	
	if(exists $hit_hash{$gid}) {
		@hit_array = @{$hit_hash{$gid}};
		push (@hit_array, $hitid);
		$hit_hash{$gid} = [@hit_array];
		@hit_array = ();
		}

	else {
		push (@hit_array, $hitid);
		$hit_hash{$gid} = [@hit_array];
		@hit_array = ();
		}
	}
	
close(INFILE);
close(HITLIST);

###########################################################################

open (COUNTFILE, ">>$countlist") or die "error opening output file!\n";
if(%count_hash) {
	foreach $geneid (sort {$count_hash{$b} <=> $count_hash{$a}} keys %count_hash) {
		print COUNTFILE "$geneid"."\t"."$count_hash{$geneid}"."\n";
		}
	%count_hash = ();
	$next_output = 1;
	}
close(COUNTFILE);

###########################################################################

print "Total sequences: $total_count\n";
print "Total Ref Sequences: $ref_count\n";

###########################################################################

open (OUTFILE, ">>$output") or die "error opening output!\n";

if($next_output == 1) {
	open (REFFILE, "$countlist") or die "run error!\n";

	while (<REFFILE>) {
		chomp($_);
		@count_array = split(/\t/, $_);
		$gnid = $count_array[0];
		$gnid_count = $count_array[1];
		@hits_array = @{$hit_hash{$gnid}};
		
		print OUTFILE ">$gnid\n";

		foreach(@hits_array) {
			print OUTFILE "$_\t";
			}
		print OUTFILE "\n";
		}
	}
close (REFFILE);
close (OUTFILE);
%hit_hash = ();

###########################################################################

exit;