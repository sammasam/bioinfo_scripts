#!usr/bin/perl
use warnings;
use strict;

###########################################################################
##### PROTEIN ALIGNMENT TO NUCLEOTIDE CONVERTER ###########################
##### get nucleotide sequence for amino acid alignment ####################
##### replace residues with codons and keep gaps ##########################
###########################################################################
##### 1st ARGV = nt sequence file #########################################
##### 2nd ARGV = aa alignment file ########################################
##### 3rd ARGV = output nt alignment ######################################
###########################################################################



###########################################################################
##### GLOBAL VARS #########################################################
###########################################################################


##### variables for putting each nucleotide sequence  
##### into a hash with headers as keys and seq as value

my $seq_key = "";

my $seq_header= "";
my $seq = "";
my %seq_hash = ();
my $seq_counter = 0;

##### variables for hash for each alignment and header array

my @headers = ();

my $align_key = "";

my $align_header = "";
my $align_seq = "";
my %align_hash = ();
my $align_counter = 0;

##### variables for manipulating string

my $header = "";
my $alignment = "";
my $nucleotides = "";
my @aligned = ();
my $start = 0;
my $codon = "";

###########################################################################
##### HASH FOR SEQUENCES ##################################################
###########################################################################

open(SEQFILE, "$ARGV[0]") or die "error opening fastq file!\n";

while(<SEQFILE>) {
	chomp($_);
	#### match sequence header -> $seq_header
	if($_ =~ m/^\>/) {
		if($seq){	
		$seq_hash{$seq_header} = $seq;
		}
	$seq_header = $_;
	$seq_counter ++;
	### reset $seq
	$seq = "";
	}
	else {
	$seq .= $_;
	}
}

##### store last sequence in %seq_hash

$seq_hash{$seq_header} = $seq;

close(SEQFILE);

###########################################################################
##### TEST ##### print the hash (keys, values) ############################
###########################################################################
#
#foreach $seq_key (keys %seq_hash)
#{
#  print "Key: $seq_key\nValue: $seq_hash{$seq_key}\n";
#}	
###########################################################################
##### NEXT HASH FOR ALIGNMENT##############################################
###########################################################################
##### ADD ARRAY OF HEADERS ################################################
###########################################################################


open(ALIGNFILE, "$ARGV[1]") or die "error opening fastq file!\n";

while(<ALIGNFILE>) {
	chomp($_);
	### match  header -> $align_header
	if($_ =~ m/^\>/) {
	push(@headers, $_);

		if($align_seq) {	
		$align_hash{$align_header} = $align_seq;
		}

	$align_header = $_;
	$align_counter ++;
	### reset $align_seq
	$align_seq = "";
	}
	else {
	$align_seq .= $_;
	}
}

##### store last sequence in %align_hash

$align_hash{$align_header} = $align_seq;

close(ALIGNFILE);


###########################################################################
##### TEST ##### print the hash (keys, values)
#foreach $align_key (keys %align_hash)
#{
#  print "Key: $align_key\nValue: $align_hash{$align_key}\n";
#}	
#
###########################################################################
##### TEST ##### print array
#foreach(@headers)
#{
#print "$_\n";
#}
#
###########################################################################



###########################################################################
##### MAKE OUTFILE ########################################################
###########################################################################

open(OUTFILE, ">>$ARGV[2]") or die "error opening output!\n";

	foreach(@headers) {
	$header = $_;
	$alignment = $align_hash{$header};
	$nucleotides = $seq_hash{$header};
	@aligned = split(//, $alignment);
	$start = 0;
	$codon = "";

		if(@aligned) {
		print OUTFILE "$header\n";

			foreach(@aligned) {
				if($_ =~ m/-/) {
				print OUTFILE "---";
				}
				elsif($_ =~m/[A-Z]/) {
				$codon = substr($nucleotides, $start, 3);
				print OUTFILE "$codon";
				$start += 3;
				$codon = "";
				}
			}
		}
	print OUTFILE "\n";
	}

close(OUTFILE);

print "nt sequences: $seq_counter\n";
print "aligned seqs: $align_counter\n";

###########################################################################

exit;
















