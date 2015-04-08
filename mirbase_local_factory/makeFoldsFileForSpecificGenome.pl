#!/usr/bin/env perl

use strict;
use warnings;

my $species = $ARGV[0];


open(CUSTOM_SPECIES_FOLDS_FH, ">", "$species"."_hairpin_precursors.fa.folds");

open(HAIRPIN_FH, ,"<", "./hairpin.fa");


my $FOUND_VALID_ENTRY = 'FALSE';
my $hairpinSeq = "";


while(my $line = <HAIRPIN_FH>){

    chomp $line;    
    
    if($line =~ /^>/){

        if($hairpinSeq ne ""){
		    
            print TMP_RNAFOLD_OUTPUT_FH "$hairpinSeq\n";


	    my $tmpRNAfoldOutput = `RNAfold < tmpRNAfoldInput.fa`;
	    print CUSTOM_SPECIES_FOLDS_FH "$tmpRNAfoldOutput";
            
	    $FOUND_VALID_ENTRY = 'FALSE';
            $hairpinSeq = "";
        }
        
        
        
        #$line =~ s/>//;
        if($line =~ /^>$species/){

	    $line =~ s/stem-loop//;
	    open(TMP_RNAFOLD_OUTPUT_FH, ">", "tmpRNAfoldInput.fa");
            $FOUND_VALID_ENTRY = 'TRUE';
            #my @vals = split(/\s/, $line);
            print TMP_RNAFOLD_OUTPUT_FH "$line\n"; 
        } 
    } else{
        if($FOUND_VALID_ENTRY eq 'TRUE'){ 
            $hairpinSeq = $hairpinSeq.$line;
        }
    }

}

print TMP_RNAFOLD_OUTPUT_FH "$hairpinSeq\n";

my $tmpRNAfoldOutput = `RNAfold < tmpRNAfoldInput.fa`;
print CUSTOM_SPECIES_FOLDS_FH "$tmpRNAfoldOutput";

close CUSTOM_SPECIES_FOLDS_FH;
