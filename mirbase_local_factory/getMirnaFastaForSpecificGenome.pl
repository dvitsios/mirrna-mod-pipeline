#!/usr/bin/env perl

use strict;
use warnings;

my $species = $ARGV[0];


open(CUSTOM_SPECIES_FASTA_FH, ">", "$species"."_hairpin_precursors.fa");

open(HAIRPIN_FH, ,"<", "./hairpin.fa");

my $FOUND_VALID_ENTRY = 'FALSE';
my $hairpinSeq = "";

while(my $line = <HAIRPIN_FH>){
    
    
    if($line =~ /^>/){

        if($hairpinSeq ne ""){
            print CUSTOM_SPECIES_FASTA_FH "$hairpinSeq";
            $FOUND_VALID_ENTRY = 'FALSE';
            $hairpinSeq = "";
        }
        
        
        
        #$line =~ s/>//;
        if($line =~ /^>$species/){

            $FOUND_VALID_ENTRY = 'TRUE';
            #my @vals = split(/\s/, $line);
            #print CUSTOM_SPECIES_FASTA_FH "$vals[0]\t$vals[1]\t"; 
            print CUSTOM_SPECIES_FASTA_FH $line; 
        } 
    } else{
        if($FOUND_VALID_ENTRY eq 'TRUE'){
            
            $hairpinSeq = $hairpinSeq.$line;
            
        }
    }

}

print CUSTOM_SPECIES_FASTA_FH "$hairpinSeq";
