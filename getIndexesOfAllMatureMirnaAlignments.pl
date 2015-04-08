#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Std;

# -d: blast hits depth, number of blast hits to (initially) retain and then filter to keep only the valid ones
my %opts = ();
getopt('d', \%opts);


my $script_path=$0;
$script_path=~ s/\/[^\/]+$//;

my $blast_hits_depth = 10; #default value for blast_hits_depth: 10
if(defined($opts{'d'})){
	$blast_hits_depth = $opts{'d'};
} 

my $blastEval = 0.001;
my $wordsize = 7;
my $mirnaIdsAndSeqsFile = "matureToUnivIDMappings.txt";


=begin comment
my $fullSeqBlastOutput = "tmpCurMirFullSeqBlastOutput.txt";
my $species = "$script_path/data/hsa_hairpin_precursors.fa";


my $moption = 8;
my $eval = 0.001;
my $wordsize = 11;
system("echo -e \">testId\nCCCCTGGAGTGTGACAATGGTGTTTGCCCC\" | $script_path/blastall -p blastn -e $eval -d $species -m $moption -v 1 -b $blast_hits_depth -o $fullSeqBlastOutput -F F");
exit;
}
=end comment
=cut


my %matureseqs = ();

my $globalMatureMirAlignmentsFilename = "indexesOfAllMatureMirnaAlignments.txt";
open(MY_OUTPUT_FH, '>', $globalMatureMirAlignmentsFilename);
print MY_OUTPUT_FH "MIRNA_ID\tHAIRPIN_ID\tFULL_ALIGN_HAIRPIN_START_INDEX\tFULL_ALIGN_HAIRPIN_END_INDEX\n";	

getSequencesForAllMatureMirnas();

close(MY_OUTPUT_FH);



###########################################
############   SUBROUTINES   ##############
###########################################
sub getSequencesForAllMatureMirnas{

	open(MYFILE, $mirnaIdsAndSeqsFile);
	
	while(<MYFILE>){
		my ($mirId, $univId, $mirSeq) = split("\t",$_);	
		$mirId =~ s/R/r/;
		$matureseqs{$mirId} = $mirSeq;

		my $genomeStr = substr($mirId, 0, index($mirId, '-'));

		
		run_blast_for_cur_seq($mirId, $mirSeq, $genomeStr);

	}
	close(MYFILE);
}

sub run_blast_for_cur_seq{

	my $mirId = shift;
	my $mirSeq = shift;
	my $genomeStr = shift;


	my $fullSeqBlastOutput = "tmpCurMirFullSeqBlastOutput.txt";
	my $species = "$script_path/data/$genomeStr"."_hairpin_precursors.fa";
	

	system("echo -e \">$mirId\n$mirSeq\" | $script_path/blastall -p blastn -e $blastEval -d $species -m 8 -o $fullSeqBlastOutput -v 1 -b $blast_hits_depth -F F -W $wordsize");             
	my $valid_blast_hits_ref = get_valid_blast_hits_from_single_miRNA_BLAST($fullSeqBlastOutput);
	my @valid_blast_hits = @$valid_blast_hits_ref;


	# keep only one entry per with hit the same hairpin seq
	my %hairpins_checked = (); 

	foreach(@valid_blast_hits){
	
		my $cur_full_seq_hit = $_;
	
		#print "cur_full_seq_hit: $cur_full_seq_hit\n";

		# 'fs' stands for full sequence: it refers to the blast output when we align the full original sequence of the
		# currently found $id against the whole mirbase.	
		my ($fs_id,$fs_match,$fs_perc,$fs_mm,$fs_mismatch,$fs_gaps,$fs_qs,$fs_qe,$fs_ss,$fs_se,$fs_ev)=split("\t",$cur_full_seq_hit);	

		if(!defined($hairpins_checked{$fs_match})){
			print MY_OUTPUT_FH "$mirId\t$fs_match\t$fs_ss\t$fs_se\n";

			# define it so that it doesn't write another entry with the same hairpin hit
			$hairpins_checked{$fs_match} = 1;	
		} 



	}	
}



sub get_valid_blast_hits_from_single_miRNA_BLAST{

        # get the first best hit from a single miRNA BLAST that complies with the condition: $ss > $se

        my $tmpBlastOutputFile = shift;

	my @valid_blast_hits = ();

        my $cnt = 0;
        my $prev_blast_hit_id = "";
        my $cur_blast_hit_id = "";
        my $CUR_HIT_IS_CHECKED = "FALSE";

        open(FULL_SEQ_PROC, '<', "$tmpBlastOutputFile");
        while(<FULL_SEQ_PROC>){
                $cnt++;
                chomp();

                my ($id,$match,$perc,$mm,$mismatch,$gaps,$qs,$qe,$ss,$se,$ev)=split("\t",$_);

                if($cnt == 1){
                        $prev_blast_hit_id = $id;
                        $cur_blast_hit_id = $id;

                        if($ss > $se){
                                next;
                        } else{
                                push(@valid_blast_hits, $_);
				#return($_);
                        }
                } else{

                        if($ss > $se){
                                next;
                        } else{
                                push(@valid_blast_hits, $_);
				#return($_);
                        }
                }

        }
        close(FULL_SEQ_PROC);
        system("rm $tmpBlastOutputFile");

	return(\@valid_blast_hits);
}

