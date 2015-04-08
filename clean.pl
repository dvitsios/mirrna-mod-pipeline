#!/usr/bin/perl

use warnings;
use strict;

$|=1;

my $READ_LENGTH_LOW_THRESHOLD = 10;
my $READ_LENGTH_HIGH_THRESHOLD = 32; #default=32, assuming a max. tolerance of 10nt long modifications in total (both from 5p and 3p)


my $have_fq = 0;
my $id = "";
my $sequence = "";
my $count = 1;
my $good=0;
my $total=0;


#open(my $print_fh, '>', 'clean_pl_out.txt');


while (<STDIN>) {
	$have_fq = 1 if $. == 1 && /^@/;
	chomp();
	if ($have_fq) {
	
	      #print $print_fh "FASTQ FILE\n";

	      if (/^@(\S+)/ && $. % 4 == 1) {
		 $total++;
		 $id=$1;
			 if (/(\d+)$/) {
				$count = $1;
			 } else {
				$count = 1;
			 }
	      } elsif ($. % 4 == 2) {
		 $sequence = $_;
		 $total++;
		 if ($sequence && $id) {
		 	print ">$id" . "_x$count" . "\n$sequence\n";
		 	$good++;
		 }
	      }
	      else {
		 $id = "";
	      }
	} 
	else {
		if (/^(>.*)/) {
		        #print $print_fh "NO FASTQ FILE\n";
			$total++;
			$id=$1;	
		}
	      	else {
			$sequence=$_;
	
			if ($sequence){

				if ((length($sequence) >= $READ_LENGTH_LOW_THRESHOLD) && (length($sequence) <= $READ_LENGTH_HIGH_THRESHOLD)){
					print "$id\n$sequence\n";
					#print $print_fh "$id\n$sequence\n";
					$good++;
				}
			}
		}
	}
}

print STDERR "cleaner.pl Passed $good of $total Total Sequences\n";
