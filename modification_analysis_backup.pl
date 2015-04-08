#!/usr/bin/env perl
$|=1;

use strict; 
#use warnings;
use Data::Dumper;
use Cwd;
use JSON::XS qw(encode_json decode_json);
use File::Slurp qw(read_file write_file);
use Switch;


my $base_script_dir = getcwd;
chdir $base_script_dir;

print "base_script_dir: $base_script_dir\n";


# basic parameters
my $DEBUG_MODE = 'FALSE'; #when set to 'TRUE' it prints more info about my data content and conditions getting triggered
my $blast_hits_depth = 3;
my $blastEval = 0.001;
my $wordsize = 7;

my $VERBOSE = 'FALSE';
my $CONSIDER_READS_WITH_BOTH_MODS_TWICE_IN_TOTAL_READS = 'FALSE'; # FALSE is more correct to get the right ratios for SNPs etc.
my $PRINT_ALL_ALIGNMENTS_AT_HTML = 'TRUE';
my $SNP_RATIO_THRESHOLD = 0.7;
my $ADD_POSITION_OF_SNP_OR_ADAR_INFO_IN_FINAL_IDENTIFIER = 'TRUE';
my $KEEP_INTERNAL_MODS = 'FALSE';
my $STRETCH_LENGTH_THRESHOLD_FOR_SNPS = 2;

my $mod_types_separator = '=>';
my %all_mature_mirnas_align_indexes = ();
get_all_all_mature_mirnas_align_indexes();
#print Dumper(\%all_mature_mirnas_align_indexes);




# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
# [FOR EXPERIMENTING / DEBUGGINF PUPROSE ONLY]:
# PARAMETERS in order to retain the reads that in the end map to a certain list of miRNAs that I want to study separately.
# Currently, the interesting cases are miRNAs with histo-mod profile expanded in both arms of the hairpin.
# This indicates miRNAs that do exist actually but have not been annotated yet, maybe due to lack of convincing evidence.
my $RETAIN_SEQS_FOR_DEBUG_DATASET = 'FALSE';
my $ARM_INFO_IN_MIRS = 'FALSE';

# Two separate arms, no arm info in the match id
my @mir_ids_to_retain = ('hsa-mir-320a', 'hsa-mir-448','hsa-mir-764','hsa-mir-544a','hsa-mir-1264','hsa-mir-217');

# One arm with adjacent 'extra arm', with arm info in the match id
#my @mir_ids_to_retain = ('hsa-mir-370', 'hsa-mir-219a-2', 'hsa-mir-138-2', 'hsa-mir-181b', 'hsa-mir-124', 'hsa-mir-331');


# strong candidates to be non-actual miRNAs
#my @mir_ids_to_retain = ('hsa-mir-4792', 'hsa-mir-7641', 'hsa-mir-664a', 'hsa-mir-6087', 'hsa-mir-3607', 'hsa-mir-297', 'hsa-mir-4791', 'hsa-mir-1273a', 'hsa-mir-4485', 'hsa-mir-1469', 'hsa-mir-6516', 'hsa-mir-449b', 'hsa-mir-5010', 'hsa-mir-4508', 'hsa-mir-451b', 'hsa-mir-4324', 'hsa-mir-4492', 'hsa-mir-636', 'hsa-mir-4419a', 'hsa-mir-1178', 'hsa-mir-1268a');

my $debug_data_arm_info_identifier;
if($ARM_INFO_IN_MIRS eq 'TRUE'){
	$debug_data_arm_info_identifier = 'AND'
} else{
	$debug_data_arm_info_identifier = 'NO';
}

my $TMP_DEBUG_FH;
my $debug_dataset_filename;
if($RETAIN_SEQS_FOR_DEBUG_DATASET eq 'TRUE'){
	#$debug_dataset_filename = "mirs_with_2_arms_".$debug_data_arm_info_identifier."_arm_info.fa";
	$debug_dataset_filename = "non_actual_mir_candidates.fa";
	open($TMP_DEBUG_FH, '>', $debug_dataset_filename);
}
# =*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*


# ######################
# Current implementation requires explicit definition
# of the output folder in the arguments list!!!
# ######################


my %valid_mature_ids = getAllValidMatureIds();
#my %valid_mature_ids = %$valid_mature_ids_ref;  # if using json


my $upload_dir;

my @filenames;
for(my $i=0; $i<$#ARGV-2; $i++){
	$filenames[$i] = $ARGV[$i];
}
print "filenames:\n";
print @filenames;

my $outputFolder = $ARGV[$#ARGV-2];
print "\noutputFolder: $outputFolder\n";

my $RUN_FROM_PIPELINE = $ARGV[$#ARGV-1];

print "RUN_FROM_PIPELINE: $RUN_FROM_PIPELINE\n";

my $curTaxonName = $ARGV[$#ARGV];
print "curTaxonName: $curTaxonName\n";


print STDERR "RUN_FROM_PIPELINE: $RUN_FROM_PIPELINE\n";



if($RUN_FROM_PIPELINE == 1){
	$upload_dir = "";
} else {
	$upload_dir = ".";
}

foreach(@filenames){
	print "input: $_\n";
}

print "output: $outputFolder\n";


my $script_path=$0;
$script_path=~ s/\/[^\/]+$//;
print "$script_path\n";


my $species = "";

if($curTaxonName ne ''){
	$species = "$script_path/data/$curTaxonName.fasta";
} else {
	$species = "$script_path/data/hsa.fasta";
}



# ============================================== FILE HANDLING ==============================================
# Files to be created:
# 1. Create and open 'processed.counts' file.
my $output_processed_counts_file = $outputFolder."/processed.counts";
print "output_processed_counts_file: $output_processed_counts_file\n";
open (FILEOUT4,">$output_processed_counts_file");

# 2. Create new2arms.counts file: list with miRNAs that give a histo-mod profile with 2 arms, 
# one of which has not been annotated in mirbase. 
my $new_mirs_with_2_arms_fname = $outputFolder."/new2arms.counts";
open(NEW_2ARMS_FH, '>', $new_mirs_with_2_arms_fname);


# 3. Create and open 'non_actual_mirs.counts' file: list with potentially non actual miRNAs.
my $non_actual_mirs_fname = $outputFolder."/non_actual_mirs.counts";
open(NON_ACTUAL_MIRS_FH, '>', $non_actual_mirs_fname);
# ===========================================================================================================


print "miRBase Mapping Script, Enrightlab\n";
print "Will store final tabular counts in: $output_processed_counts_file\n";
print "-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n";
my $i=0;


my $edits=1;
my $depth_threshold=5;

my %mirs=();
my %reads=();
my %depths=();
my %rawdepth=();
my %maxdepth=();
my %totalreads=();
my %allreads=();
my %alldepths=();
my %threepreads=();
my %fivepreads=();
my %mmreads=();
my %mirseqs=();
my %totaldepths=();


print "\nProcessing files:\n\n";

print "Reading Hairpin lengths\n";

my $id;
my %hairpinseq=();
my %foldseq=();
my %hairpinlength=();
my %matureseqs = ();


# for each miRNA hit store its read ID and hairpin identifier that it has been mapped to  
my %id_to_hairpin_match_hashmap = ();
my %id_to_mature_match_hashmap = ();

# store all mature miRNA sequences to %matureseqs Hash:
getSequencesForAllMatureMirnas();


open(FILE,$species);

while(<FILE>){
	chomp;
	if (/^>(\S+)/){
		$id = $1;
	} else {
		$hairpinseq{$id}.=$_;
		$hairpinlength{$id}=length($hairpinseq{$id});
	}
}
close(FILE);
#print Dumper(\%hairpinseq);



open(FILE,"$species.folds");
while(<FILE>){
        chomp;
        if (/^>(\S+)/){
                $id=$1;
        } else {
                $foldseq{$id}.=$_;
        }
}
close(FILE);
#print Dumper(\%foldseq);



##########################################
############## MAIN LOOP #################
##########################################
foreach my $filename(@filenames){

	print NEW_2ARMS_FH "filename: $filename|\tmiRNA\treads\thisto\tmodif\n";

	my $lastSlashIndex = rindex($filename, '/');
	my $filename_id = substr($filename, $lastSlashIndex+1);


print "Species: $species Filename:$filename\n";

my $aln_file="$filename" . ".alignments.html";
my $mat_file="$filename" . ".mat";

print "Storing This Samples Alignments in: $aln_file\n";
print "Storing This Samples Matrices in: $mat_file\n";

open (FILEOUT1,">$aln_file");
open (FILEOUT2,">$mat_file");
print FILEOUT2 "\tA\tC\tG\tU\n";

my %alignments=();
undef(%rawdepth);

my $depth=0;

%mirseqs=();
my $mid="";


if($edits){
	print "Storing Sequences\n";
	open (FILE,"gunzip -c $upload_dir/$filename |");

	while(<FILE>){
		chomp;
		
		if(/^>(\S+)/){
			$mid=$1;
		} else {
			$mirseqs{$mid}.=$_;
			$mirseqs{$mid}=~ s/T/U/g;
		}

	}
	close(FILE);
	print "done\n";
}
#print Dumper(\%mirseqs);


my $cat_function="cat";
print "Processing $upload_dir/$filename\n";


my $blastOutputFile = "$filename.blast_out.txt";
print "gunzip -dcf $upload_dir/$filename | dos2unix | $script_path/clean.pl |$script_path/blastall -p blastn -e $blastEval -d $species -m 8 -o $blastOutputFile -v 1 -b $blast_hits_depth -F F -W $wordsize 2> $filename.blast.stderr";
system("gunzip -dcf $upload_dir/$filename | dos2unix | $script_path/clean.pl |$script_path/blastall -p blastn -e $blastEval -d $species -m 8 -o $blastOutputFile -v 1 -b $blast_hits_depth -F F -W $wordsize 2> $filename.blast.stderr");


my $SHOW_PLAIN_BLAST_OUTPUT = 'FALSE';
if($SHOW_PLAIN_BLAST_OUTPUT eq 'TRUE'){
	open(PROC, '<', "$blastOutputFile");
	while(<PROC>){
		print "$_\n";
	}
	close(PROC);
}

# keep the first valid hit for each miRNA hit at the blast output
keep_best_blast_hit_per_miRNA($blastOutputFile);



# Initially collect the information about the total number of reads
# for each idientified miRNA (with all types of any modifications) 
my $totalreads_per_mirna_ref = get_totalreads_per_mirna($blastOutputFile, $filename);
my %totalreads_per_mirna = %$totalreads_per_mirna_ref;

print "totalreads_pre_mirna:\n";
print Dumper(\%totalreads_per_mirna);


# Identify and classify stretches using the info from the %totalreads_per_mirna hash.
my $cnt = 0;
open(PROC, '<', $blastOutputFile);
while(<PROC>){

	$cnt++;
	chomp();
	# qs, qe: query start, query end
	# ss, se: subject start, subject end
	# mm: alignment length (?)
	my ($id,$match,$perc,$mm,$mismatch,$gaps,$qs,$qe,$ss,$se,$ev)=split("\t",$_);
	



	print "===*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*===\n";
	print "$_\n";

	if($RETAIN_SEQS_FOR_DEBUG_DATASET eq 'TRUE'){

		#if(index($match, 'hsa-mir-219a-2') != -1 || index($match, 'hsa-mir-370') != -1 || index($match, 'hsa-mir-138-2') != -1 ){
		#	print "[selected mir] id: $id, match: $match\n";			
		#}
	
		my @selected_read_matches = grep(/^$match$/, @mir_ids_to_retain);
		my $selected_read_matches_len = scalar @selected_read_matches;


		#if(index($match, $mir_id_to_retain) == -1){
		#	next;
		#} else{
		#	print $TMP_DEBUG_FH ">$id\n$mirseqs{$id}\n";
		#}
		
		if($selected_read_matches_len >= 1){
			print $TMP_DEBUG_FH ">$id\n$mirseqs{$id}\n";
		} else{
			next;
		}

	}
	



	print "\n";
	
	my $fail=0;
	my $sample_reads++;

	my $nice_id = $id;
	$nice_id =~ s/_x(\d+)/ Depth:$1 Modification:/g;
	my $this_depth=$1;

	if($VERBOSE eq 'TRUE'){
		print "$_\n";
		print "this_depth: $this_depth\n";
	}


	$allreads{$filename}++;
	$alldepths{$filename}+=$this_depth;

	my $mirlength=length($mirseqs{$id});

	# Remove 'edge effect' hits where the miR maps off the boundary of a precursor
	if (($ss-$qs) < 0){
		print "[FAIL]: ss < qs\n";
		$fail=1;
		print "[DBG]: critical error => (ss < qs)!\n";
	}

	if ($edits){

		if (($ss+($mirlength - $qs)) > $hairpinlength{$match}) {
			if($DEBUG_MODE eq 'TRUE'){
				print "mirlength: $mirlength, match: $match, hairpinlength{$match}: $hairpinlength{$match}\n";
				print "[FAIL]: ss+(mirlength - qs)) > hairpinlength{match}\n";
			}
			$fail=1;
		}

		#Exclude low depth hits
		if ($this_depth < $depth_threshold){
			if($DEBUG_MODE eq 'TRUE'){
				print "[FAIL]: this_depth < depth_threshold\n";
			}
			$fail=1;
		}
	}



	# after passing the initial validation process
	if (!$fail){


		# keep the coupling information between the current read id and the aligned hairpin id
		$id_to_hairpin_match_hashmap{$id} = $match;

		my $arm="";

		if ((($ss+$se)/2) <= ($hairpinlength{$match}/2)){
			$arm="5p";
		} else {
			$arm="3p";
		}

		my $aligned="";

		if ($edits){
			# add spaces for the 'missing' positions starting from the begining of the hairpin.		
			for ($i=0; $i<($ss-$qs); $i++){
				$aligned .= " ";
			} 	

			$aligned .= $mirseqs{$id} . "\t$nice_id";
		}


		my $initial_id=$match;
		
		$match = get_correct_id_with_arm_info($match, $arm);

		
		# keep the coupling information between the current read id and the aligned mature miRNA id
		$id_to_mature_match_hashmap{$id} = $match;

		#print ">>$mirseqs{$id} $mirlength\n";
		my $five_nt="";
		my $three_nt="";	
		if ($edits){
			if ($mm < $mirlength){
				if ($qs != '1'){
					$five_nt=substr($mirseqs{$id},0,$qs-1);
				}
				if ($qe < $mirlength){
					$three_nt=substr($mirseqs{$id},$qe,$mirlength-$qe);
				}
			}
		}


		# => <===============================================> <=

		my $curMirnaFullSeq = $matureseqs{$match};
		if($VERBOSE eq 'TRUE'){
			print "match: $match, curMirnaFullSeq: $curMirnaFullSeq\n";
		}



		my $cur_mir_ratio = $this_depth / $totalreads_per_mirna{$match};

	
		if($VERBOSE eq 'TRUE'){
			#print ">>> stretchSeq_5p: $stretchSeq_5p\n";
			#print ">>> stretchSeq_3p: $stretchSeq_3p\n";
			print "this depth: $this_depth\n";
			print "totalreads for this miRNA: $totalreads_per_mirna{$match}\n";
			print "ratio of this miRNA: $cur_mir_ratio\n";	
		}


		# get alignment indexes from full miRNA sequence.
		#my $fs_ss = $all_mature_mirnas_align_indexes{$match}[0]; 
		#my $fs_se = $all_mature_mirnas_align_indexes{$match}[1];
		my $cur_aligned_hairpin_id = $id_to_hairpin_match_hashmap{$id};

		my $fs_ss = $all_mature_mirnas_align_indexes{$match}{$cur_aligned_hairpin_id}[0]; 
		my $fs_se = $all_mature_mirnas_align_indexes{$match}{$cur_aligned_hairpin_id}[1]; 
	      	


		# check content of the 5p and 3p stretches so that it is classified eventually
	        # either as a modification or as ADAR-editing / SNP.
		
		

		# identify first 'internal' ADARS and SNPs in case we have a single mismatch
		my $internal_mismatch_identifier = '';
		if($mismatch == 1){

			my $mod_position = -1;

			my $querySeq = $mirseqs{$id};
			my $subjectSeq = $hairpinseq{$initial_id};

			my $aligned_query_seq = substr($querySeq, $qs-1, $qe-$qs+1);
			my $aligned_subject_seq = substr($subjectSeq, $ss-1, $se-$ss+1);

			#print "aligned_query_seq: $aligned_query_seq\n";
			#print "aligned_subjt_seq: $aligned_subject_seq\n";

			my $stretch_length = 0;
			if($five_nt){
				$stretch_length = length($five_nt);
			}
			my $original_to_cur_hairpin_align_offset = $ss - $fs_ss;

			#print "stretch_length: $stretch_length\n";
			#print "ss: $ss, fs_ss: $fs_ss\n";
			#print "original_to_cur_hairpin_align_offset: $original_to_cur_hairpin_align_offset\n";

			for(my $al=0; $al<length($aligned_query_seq); $al++){
				
				my $qchar = substr($aligned_query_seq, $al, 1);
				my $schar = substr($aligned_subject_seq, $al, 1);
		
				#print "qchar: $qchar, schar: $schar\n";
	
				if(($qchar eq 'G') & ($schar eq 'A')){
					$mod_position = $qs + $al;

					my $mod_position_relative_to_mir_start = $mod_position - ($stretch_length - $original_to_cur_hairpin_align_offset);
	
					$internal_mismatch_identifier = "internal_adar_G_$mod_position_relative_to_mir_start";
					last;
				} elsif($qchar ne $schar){
					if($cur_mir_ratio > $SNP_RATIO_THRESHOLD){
						$mod_position = $qs + $al;
						my $mod_position_relative_to_mir_start = $mod_position - ($stretch_length - $original_to_cur_hairpin_align_offset);
						
						$internal_mismatch_identifier = "internal_snp_$qchar"."_$mod_position_relative_to_mir_start";
						last;
					} else{
						if($KEEP_INTERNAL_MODS eq 'TRUE'){
							$mod_position = $qs + $al;
							my $mod_position_relative_to_mir_start = $mod_position - ($stretch_length - $original_to_cur_hairpin_align_offset);
							
							$internal_mismatch_identifier = "internal_mod_$mod_position_relative_to_mir_start";
							last;
						}
					}
				}

			}	
		}





		if ($edits == 0){
			$five_nt="";
			$three_nt="";
		}
		my $modification = "";
		if ($id=~/(\d+)$/){
			$depth=$1;
		} else {
			$depth=1;
		}
		$rawdepth{$match}+=$depth;
		my $FOUND_FIVE_OR_THREE_NT = 'FALSE';
	


	
		


		my $FIVEp_STRETCH_BEYOND_MIR_LIMITS = 'FALSE';
		my $five_prime_missing_stretch_length = $ss - $fs_ss;
		my $THREEp_STRETCH_BEYOND_MIR_LIMITS = 'FALSE';
		my $three_prime_missing_stretch_length = $fs_se - $se;


		# hits from multiple loci may cause minor issue, not correct length in the end..
		# taking care of that (sub)optimally.
		if($five_prime_missing_stretch_length < 0){
			$FIVEp_STRETCH_BEYOND_MIR_LIMITS = 'TRUE';
		} 
		if($three_prime_missing_stretch_length < 0){
			$THREEp_STRETCH_BEYOND_MIR_LIMITS = 'TRUE';
		}



		if ($five_nt){

			my $five_nt_offset = $ss - $fs_ss - length($five_nt);

			print "5p-stretch: $five_nt\n"; 

			my $armEnd = "5p";
			my $stretch_final_identifier = classify_stretch($five_nt, $match, $initial_id, $mismatch, $gaps, $qs, $ss, $qe, $se, $armEnd, $cur_mir_ratio, $FIVEp_STRETCH_BEYOND_MIR_LIMITS);

			$stretch_final_identifier .= "_$five_nt_offset";

			print "stretch_final_identifier: $stretch_final_identifier\n";	
			$modification = $stretch_final_identifier;
			if($internal_mismatch_identifier ne ''){
				$modification .= $mod_types_separator.$internal_mismatch_identifier;
			}


			#$depths{$filename}{$match} += $depth;
			#$mirs{$match} = 1;
			#$modification = $stretch_final_identifier;

			$depths{$filename}{$match.$modification} += $depth;
			$mirs{$match.$modification} = 1;
			
			#if ($maxdepth{$match} < $depth){
			#	$maxdepth{$match} = $depth;
			#}
			if ($maxdepth{$match.$modification} < $depth){
				$maxdepth{$match.$modification} = $depth;
			}

			$FOUND_FIVE_OR_THREE_NT = 'TRUE';
			
			$alignments{$match}.=$aligned . " " . $modification . "\n";
		}
		if ($three_nt){

			my $three_nt_offset = $se - $fs_se + 1;

			print "3p-stretch: $three_nt\n"; 
			
			my $armEnd = "3p";
			my $stretch_final_identifier = classify_stretch($three_nt, $match, $initial_id, $mismatch, $gaps, $qs, $ss, $qe, $se, $armEnd, $cur_mir_ratio, $THREEp_STRETCH_BEYOND_MIR_LIMITS);
			
			$stretch_final_identifier .= "_$three_nt_offset";
			print "stretch_final_identifier: $stretch_final_identifier\n";	
	
				
			$modification = $stretch_final_identifier;
			if(($internal_mismatch_identifier ne '') & (!$five_nt)){
				$modification .= $mod_types_separator.$internal_mismatch_identifier;
			}
	

			$FOUND_FIVE_OR_THREE_NT = 'TRUE';
			# if the depth for this miRNA has already been counted
			# due to the already identified _5p_ mods that occur too
			# create a sequence that has only the 3p strech and
			# empty strings in all the other positions of the hairpin.
			if($five_nt){
	
				my @aligned_vals = split('\t', $aligned);
				my $tmp_aligned_seq = $aligned_vals[0];
				my $rest_of_aligned = $aligned_vals[1];
				
				my $tmp_full_length = length($tmp_aligned_seq);
				my $trailing_spaces_length = $tmp_full_length - $mirlength;


				my $tmp_aligned = '';
				for(my $t_len=0; $t_len<$trailing_spaces_length; $t_len++){
					$tmp_aligned .= ' ';
				}
				for(my $t_len=0; $t_len<$mirlength; $t_len++){
					$tmp_aligned .= 'x'; # 'x' is used to denote empty position at the miRNA
				}
			        $aligned = $tmp_aligned."\t".$rest_of_aligned;
				
				print "new aligned with 5p&3p mods: ###$aligned###\n";
				$modification .= $mod_types_separator."doubled";

			}

			#$depths{$filename}{$match} += $depth;
			#$mirs{$match} = 1;
			#$modification = $stretch_final_identifier;
		
			$depths{$filename}{$match.$modification} += $depth;
                        $mirs{$match.$modification} = 1;

                        if ($maxdepth{$match.$modification} < $depth){
                                $maxdepth{$match.$modification} = $depth;
                        }

			$alignments{$match}.=$aligned . " " . $modification . "\n";
		} 

		if($FOUND_FIVE_OR_THREE_NT eq 'FALSE'){
			# add any 'internal' SNPs, ADARs or just mods.

			$modification = $mod_types_separator."no_mod";
			if($internal_mismatch_identifier ne ''){
				$modification .= $mod_types_separator.$internal_mismatch_identifier;
			}
			

                        $mirs{$match.$modification} = 1;
			$reads{$filename}{$match}++;
			$depths{$filename}{$match.$modification} += $depth;

	
			#$mirs{$match}=1;
		        #$reads{$filename}{$match}++;
		        #$depths{$filename}{$match}+=$depth;

			my $tmp_match_and_mod_id = $match.$modification;

			if ($depths{$filename}{$tmp_match_and_mod_id} > $maxdepth{$tmp_match_and_mod_id}){
			        $maxdepth{$tmp_match_and_mod_id}=$depths{$filename}{$tmp_match_and_mod_id};
			}
			$alignments{$match}.=$aligned . " " . $modification . "\n";
		}

		if ($edits){
			#$alignments{$match}.=$aligned . " " . $modification . "\n";
			# temporary solution to the hits from multiple-loci for the same miRNA
			if(!defined($hairpinseq{$match})){
				$hairpinseq{$match}=$hairpinseq{$initial_id};
			}	
			
			if($VERBOSE eq 'TRUE'){
				print "match: $match\n";
				print "hairpinseq{match}: $hairpinseq{$match}\n";
				print "initial_id: $initial_id\n";
				print "hairpinseq{initial_id}: $hairpinseq{$initial_id}\n";
			}

			$foldseq{$match}=$foldseq{$initial_id};
		}

	 	# check if there is any inconsitency with that
	 	# because when I have both 5p and 3p mods in a read
	 	# then I consider two separate reads of the same depth
	 	# but with only the one type mod each time (5p or 3p)
	 	# It doesn't interfere though with the algorithm's logic
	 	# since those two variables (totalreads, totaldepths)
	 	# are used only for printing or writing to file purposes. 	
		$totalreads{$filename}++;
		$totaldepths{$filename}+=$depth;

	} else {
		#FAIL
	}


	if($VERBOSE eq 'TRUE'){	
		print "==============================================\n\n\n";
	}

}
close(PROC);


print "mirBase\tUnique\tDepth\n";
print "TOTAL:\t$totalreads{$filename}\t$totaldepths{$filename}\n";
print "INITIAL:\t$allreads{$filename}\t$alldepths{$filename}\n";
print "---------------------------------\n";

print FILEOUT1 "<HTML>\n";
print FILEOUT1 "<BODY BGCOLOR=\"#ECE0F8\">\n";
print FILEOUT1 "<FONT FACE=\"Arial\">\n";
print FILEOUT1 "<H2>MicroRNA Modification Report for $filename</H2>\n";
print FILEOUT1 "<BR>Reads with less than depth $depth_threshold excluded from this analysis\n";
print FILEOUT1 "<BR><BR><BR><BR>\n";


print Dumper(\%alignments);


# for each (alined) miRNA, with all the possible modification patters it has been found with.
foreach my $th (sort special2(keys(%alignments))){

	print FILEOUT1 "<HR>\n";
	print FILEOUT1 "<H3>$th</H3>\n";
	print FILEOUT1 "Depth: <b>$rawdepth{$th}</b> reads\n";
	print FILEOUT1 "<H4>Aligned Reads on Precursor</H4>\n";
	print FILEOUT1 "<pre style=\"width: 1024px; background-color: FFFFFF; border: 1px solid #bebab0; -webkit-border-radius: 10px; -khtml-border-radius: 10px; -moz-border-radius: 10px; border-radius: 10px;\">\n";

	if($PRINT_ALL_ALIGNMENTS_AT_HTML eq 'TRUE'){
		print FILEOUT1 formatseq("$hairpinseq{$th}\n");
		#print FILEOUT1 "$foldseq{$th}\n";

		my $first_left_parenthesis_idx = index($foldseq{$th}, '(');
		my $first_dot_index = index($foldseq{$th}, '.');
		my $first_fold_char_idx = $first_left_parenthesis_idx < $first_dot_index? $first_left_parenthesis_idx: $first_dot_index;

		my $tmp_only_fold_seq = substr($foldseq{$th}, $first_fold_char_idx);
		print FILEOUT1 "$tmp_only_fold_seq\n";
		print FILEOUT1 formatseq("$alignments{$th}\n");
	}
	
	print FILEOUT1 "</pre>\n";	


	#print "th: $th\n";
	#print "hairpin{th}: $hairpinseq{$th}\n";


	my @matrix=();
	my @histo=();
	my @modif=();

	my $hlength=length($hairpinseq{$th});
	my @array=split("\n",$alignments{$th});


	for (my $j=0;$j<=$hlength;$j++){
		$histo[$j]=0;
		$modif[$j]=0;
	}
	my @nts = ('A', 'C', 'G', 'U');
	
	# foreach of the alignments, with information about the associated modifications
	# appended to their miRNA identifier, analyze 'errors' content in order to decide
	# which modifications are actually modifications and which ones are SNPs or ADAR-editing.
	# In the end, decide which miRNAs are not actually miRNAs.
	for (my $i=0; $i<=$#array; $i++){
		
		my $string=$array[$i];
		#$string=~ s/(\s+\S+)\s+.*/$1/g; #keep only the sequence in the 'string' variable ## obsolete, replaced by the split('\t', $string)

	
		my @string_vals = split('\t', $string);

		$string = $string_vals[0];

		my @rest_string_vals = split('\s', $string_vals[1]);
		my $tmp_cur_depth = $rest_string_vals[1];
		$tmp_cur_depth =~ s/Depth://;	

		#print "=*=*=*> tmp_cur_depth: $tmp_cur_depth\n";


		my $strLen = length($string);
		#print "string: $string, length: $strLen\n"; 		


		for (my $j=0;$j<=$hlength;$j++){
			
			my $matrix_col_index = -1;			
			my $cur_read_nt = substr($string,$j,1);
			my $cur_hairpin_nt = substr($hairpinseq{$th},$j,1);
 
	
			switch($cur_read_nt){
				case 'A' {$matrix_col_index = 1;}
				case 'C' {$matrix_col_index = 2;}
				case 'G' {$matrix_col_index = 3;}
				case 'U' {$matrix_col_index = 4;} 	

			}
			
			if($matrix_col_index != -1){
				$matrix[$j][$matrix_col_index]++;
				#$matrix[$j][$matrix_col_index] += $tmp_cur_depth;
			}

			
			my @read_nt_matches = grep(/^$cur_read_nt$/, @nts);
			my @hairpin_nt_matches = grep(/^$cur_hairpin_nt$/, @nts);
			
			my $read_match_len = scalar @read_nt_matches;
			my $hairpin_match_len = scalar @hairpin_nt_matches;

			if($read_match_len == 1 & $hairpin_match_len == 1){
				if($cur_read_nt ne $cur_hairpin_nt){
					#$modif[$j]++;
					$modif[$j] += $tmp_cur_depth;
				} else{
					#$histo[$j]++;
					$histo[$j] += $tmp_cur_depth;
				}	
			}

		} 
	}

	my $histo_prevalence_ratios_ref = get_hist_prevalence_ratios(\@histo, \@modif, $th);
	my @histo_prevalence_ratios = @$histo_prevalence_ratios_ref;
	

	my @histo_modif_union_arr = ();

	my $histo_len = scalar @histo;

	for(my $uIdx=0; $uIdx<$histo_len; $uIdx++){
		$histo_modif_union_arr[$uIdx] = $histo[$uIdx] + $modif[$uIdx];
	}


	if($VERBOSE eq 'TRUE'){
		print "\n----------------------------\nmiRNA: $th\n";
		
		print "==> histo_prevalence_ratios array:\n";
		print "@histo_prevalence_ratios\n";

		print "modif array:\n";
		print "@modif\n";

		print "histo array:\n";
		print "@histo\n";

		print "histo_modif_union_arr:\n";
		print "@histo_modif_union_arr\n\n";
		
		#print Dumper(\%totalreads_per_mirna);
	}


	my $HAS_TWO_ARMS = detect_mirs_with_two_distinct_arms_one_non_annotated(\@histo_modif_union_arr);
	if($HAS_TWO_ARMS eq 'TRUE'){
		print NEW_2ARMS_FH ">$th\t$totalreads_per_mirna{$th}\t@histo\t@modif\n";
	}	

	
	#detect_non_actual_mirnas(\@histo, \@modif, \@histo_prevalence_ratios);	




	my $max_height=0;
	for (my $i=0;$i<=$#histo;$i++){
		if (($histo[$i]+$modif[$i])>$max_height){
			$max_height=($histo[$i]+$modif[$i]);
		}
	}


	print FILEOUT1 "<H4>Distribution across Precursor</H4>\n";
	print FILEOUT1 "<TABLE style=\'padding:0px; background-color:white; border:1px solid black; border-spacing:0; border-collapse:collapse;\'><TR>\n";
	for (my $i=0;$i<=$#histo;$i++){
		

		my $height1=($histo[$i]/$max_height)*100;
		my $height2=($modif[$i]/$max_height)*100;
		print FILEOUT1 "<TD>\n";
		print FILEOUT1 "<div style=\'position:relative; height:100px; width:6px\'>\n";
    		print FILEOUT1 "<div style=\'background-color:blue; height:" . $height1 . "px; position:absolute; bottom:0px; width:6px \' />\n";
		print FILEOUT1 "<div style=\'background-color:red; height:" . $height2 . "px;  position:absolute; bottom:". $height1 ."px; width:6px \' />\n";
		print FILEOUT1 "</div>\n";
		print FILEOUT1 "</TD>\n";
	}
	print FILEOUT1 "</TR></TABLE>\n";


	for (my $i=0;$i<=$hlength;$i++){
		print FILEOUT2 "$th" ."_"."$i:\t";
		for (my $j=1;$j<=4;$j++){
			if (!$matrix[$i][$j]){
				print FILEOUT2 "0\t";
			}
			print FILEOUT2 "$matrix[$i][$j]\t";
		}
		print FILEOUT2 "\n";
	}


}
print NEW_2ARMS_FH "-total depth: $totaldepths{$filename}";
close(NEW_2ARMS_FH);


#foreach my $th (sort special2(keys(%alignments))){

#}


print FILEOUT1 "</HTML>\n";
}
close(NON_ACTUAL_MIRS_FH);


foreach my $fz(@filenames){
        print FILEOUT4 "\t$fz";
}
print FILEOUT4 "\n";

foreach my $mir(sort special(keys(%mirs))){
        print FILEOUT4 "$mir";
        foreach my $fz(@filenames){
                if (!$depths{$fz}{$mir}){
                        $depths{$fz}{$mir}=0;
                }
                print FILEOUT4 "\t$depths{$fz}{$mir}";
        }
        print FILEOUT4 "\n";
}

print FILEOUT4 "totaldepth";
foreach my $fz(@filenames){
        print FILEOUT4 "\t$totaldepths{$fz}";
}
print FILEOUT4 "\n";




sub special {
	return($maxdepth{$b} <=> $maxdepth{$a});
}

sub special2 {
	return($rawdepth{$b} <=> $rawdepth{$a});
}

sub formatseq{

	my $seq=$_[0];

	$seq=~ s/A/<font bgcolor=\"black\" color=\"red\">A<\/font>/g;
	$seq=~ s/U/<font bgcolor=\"black\" color=\"orange\">U<\/font>/g;
	$seq=~ s/C/<font bgcolor=\"black\" color=\"blue\">C<\/font>/g;
	$seq=~ s/G/<font bgcolor=\"black\" color=\"purple\">G<\/font>/g;
	$seq=~ s/N/<font bgcolor=\"black\" color=\"grey\">N<\/font>/g;

	return($seq);
}



############################################
############# MY SUBROUTINES ###############
############################################
sub detect_mirs_with_two_distinct_arms_one_non_annotated{

	my @histo_and_modif_arr = @{$_[0]};

	my $NON_ZERO_ELEMS_THRESHOLD_FOR_TWO_ARMS_DETECTION = 35; #emperical value
	my $SINGLE_ARM_LENGTH_THRESHOLD = 15; #emperical value
	my $DIST_BETWEEN_ARMS_THRESHOLD = 8;  #emperical value
	my $UNKNOWN_PASSENGER_STRAND_MIN_READS_RATIO = 0.1;

	my $HAS_TWO_ARMS = 'FALSE';
	#my @matches = grep { $histo_and_modif_arr[$_] ~~ 0 } 0 .. $#histo_and_modif_arr;
	#my $test_num_nzelems = (scalar @matches) - (scalar @histo_and_modif_arr);

	my $test_num_nzelems = 0;

	foreach(@histo_and_modif_arr){
		if($_ != 0){
			$test_num_nzelems++;
		}
	}


	# check only miRNAs with total profile length of at least 35nt.
	if($test_num_nzelems >= $NON_ZERO_ELEMS_THRESHOLD_FOR_TWO_ARMS_DETECTION){


		my $left_arm_length = 0;
		my $right_arm_length = 0;
		my $distance_between_arms = 0;

		my $left_arm_start = -1;
		my $left_arm_end = -1;
		my $right_arm_start = -1;
		my $right_arm_end = -1;


		my $FOUND_LEFT_ARM = 'FALSE';
		my $FOUND_RIGHT_ARM = 'FALSE';
		my $FINISHED_LEFT_ARM = 'FALSE';

	
		my $cnt_index = 0;

		foreach(@histo_and_modif_arr){
			my $el = $_;

			#if($FOUND_LEFT_ARM eq 'FALSE' & $FOUND_RIGHT_ARM eq 'FALSE' & $FINISHED_LEFT_ARM eq 'FALSE' & $el == 0){
			#	next; 
			#} 
			if($FOUND_LEFT_ARM eq 'FALSE' & $FOUND_RIGHT_ARM eq 'FALSE' & $FINISHED_LEFT_ARM eq 'FALSE' & $el != 0){
				$FOUND_LEFT_ARM = 'TRUE';
				$left_arm_start = $cnt_index;
				$left_arm_length++;
			} elsif($FOUND_LEFT_ARM eq 'TRUE' & $FOUND_RIGHT_ARM eq 'FALSE' & $FINISHED_LEFT_ARM eq 'FALSE' & $el != 0){
				$left_arm_length++;
			} elsif($FOUND_LEFT_ARM eq 'TRUE' & $FOUND_RIGHT_ARM eq 'FALSE' & $FINISHED_LEFT_ARM eq 'FALSE' & $el == 0){
				$distance_between_arms++;
				$FINISHED_LEFT_ARM = 'TRUE';
			} elsif($FOUND_LEFT_ARM eq 'TRUE' & $FOUND_RIGHT_ARM eq 'FALSE' & $FINISHED_LEFT_ARM eq 'TRUE' & $el == 0){
				$distance_between_arms++;
			} elsif($FOUND_LEFT_ARM eq 'TRUE' & $FOUND_RIGHT_ARM eq 'FALSE' & $FINISHED_LEFT_ARM eq 'TRUE' & $el != 0){
				$right_arm_start = $cnt_index;
				$right_arm_length++;
				$FOUND_RIGHT_ARM = 'TRUE';
			} elsif($FOUND_LEFT_ARM eq 'TRUE' & $FOUND_RIGHT_ARM eq 'TRUE' & $FINISHED_LEFT_ARM eq 'TRUE' & $el != 0){
				$right_arm_length++;
			} elsif($FOUND_LEFT_ARM eq 'TRUE' & $FOUND_RIGHT_ARM eq 'TRUE' & $FINISHED_LEFT_ARM eq 'TRUE' & $el == 0){
				last;
			}

			$cnt_index++;
		}

		if($left_arm_length > 0 & $right_arm_length > 0){
			$left_arm_end = $left_arm_start + $left_arm_length;
			$right_arm_end = $right_arm_start + $right_arm_length;


			my $left_arm_reads_sum = 0;
			for(my $i=$left_arm_start; $i<$left_arm_end; $i++){
				$left_arm_reads_sum += $histo_and_modif_arr[$i];
			}
			my $left_arm_reads_mean = $left_arm_reads_sum / $left_arm_length;

			
			my $right_arm_reads_sum = 0;
			for(my $i=$right_arm_start; $i<$right_arm_end; $i++){
				$right_arm_reads_sum += $histo_and_modif_arr[$i];
			}
			my $right_arm_reads_mean = $right_arm_reads_sum / $right_arm_length;

			my $unknown_arm_reads_mean = $left_arm_reads_mean < $right_arm_reads_mean? $left_arm_reads_mean: $right_arm_reads_mean;	
			my $unknown_arm_ratio = $unknown_arm_reads_mean/($left_arm_reads_mean + $right_arm_reads_mean);

			
			if($VERBOSE eq 'TRUE'){	
				print "left_arm_length: $left_arm_length\n";
				print "right_arm_length: $right_arm_length\n";
				print "distance_between_arms: $distance_between_arms\n";
				print "test_num_nzelems: $test_num_nzelems\n";
				print "-----\n";
				print "left_arm_reads_mean: $left_arm_reads_mean\n";
				print "right_arm_reads_mean: $right_arm_reads_mean\n";
				print "unknown_arm_ratio: $unknown_arm_ratio\n";
			}
				

			if($unknown_arm_ratio >= $UNKNOWN_PASSENGER_STRAND_MIN_READS_RATIO){
				if($left_arm_length>=$SINGLE_ARM_LENGTH_THRESHOLD & $right_arm_length>=$SINGLE_ARM_LENGTH_THRESHOLD & $distance_between_arms>=$DIST_BETWEEN_ARMS_THRESHOLD){
					$HAS_TWO_ARMS = 'TRUE';
				}
			}
		}	
	}
	
	#TRUE of FALSE
	return($HAS_TWO_ARMS); 
}



# ==> get_hist_prevalence_ratios <==
# This module is independent from the module that has to detect 
# any SNPs or ADAR-editing events.
# In this case we only worry about the final sequencing profile 
# acquired for each of the identified miRNAs and we want to infer
# only from this profile wether this miRNA a real miRNA or just an artefact.
sub get_hist_prevalence_ratios{	

	my @histo = @{$_[0]};
	my @modif = @{$_[1]};
	my $th = $_[2];

	#print "th: $th\n";

	my @histo_prevalence_ratios = ();
	for(my $i=0;$i<=$#histo;$i++){
	
		#print "i: $i , modif: $modif[$i], histo: $histo[$i]\n";

		if($histo[$i] != 0){
			$histo_prevalence_ratios[$i] = $histo[$i] / ($histo[$i]+$modif[$i]);
		} else{
			$histo_prevalence_ratios[$i] = 0;
		}
	}
	
	return(\@histo_prevalence_ratios);
}

sub checkIfMatureMirIdIsValid{

	my $subject_match = shift;
	my $FOUND_MID_HIT = 'FALSE';
	
	if(defined($valid_mature_ids{$subject_match})){
		$FOUND_MID_HIT = 'TRUE';
	}
	return $FOUND_MID_HIT;
}



# %hairpinseq: sequences from all miRNA sequences in mirbase
# %hairpinlength: lengths of all miRNA sequences in mirbase
# %mirseqs:	sequences from all miRNA sequences in the input file


sub getSequencesForAllMatureMirnas{

	my $mirnaIdsAndSeqsFile = "matureToUnivIDMappings.txt";
	open(MYFILE, $mirnaIdsAndSeqsFile);
	
	while(<MYFILE>){
		my ($mirId, $univId, $mirSeq) = split("\t",$_);	
		$mirId =~ s/R/r/;
		$matureseqs{$mirId} = $mirSeq;
	}
	close(MYFILE);
}


sub keep_best_blast_hit_per_miRNA{

        # keep for each miRNA the first blast (this will also have the best match score, e.g. 
        # no mismatches and/or gaps) hit that has $ss < $se and store only these entries in the initial blast output file.

        # blastOutputFile := $filename.blast_out.txt
	my $blastOutputFile = shift;

        my $cnt = 0;
        my $prev_blast_hit_id = "";
        my $cur_blast_hit_id = "";
        my $CUR_HIT_IS_CHECKED = "FALSE";

        my $filteredBlastOutputFile = $blastOutputFile.".filtered";
        open(FILTERED_FH, '>', $filteredBlastOutputFile);


        open(BEST_HIT_PROC, '<', "$blastOutputFile");
        while(<BEST_HIT_PROC>){

                $cnt++;
                chomp();

                my ($id,$match,$perc,$mm,$mismatch,$gaps,$qs,$qe,$ss,$se,$ev)=split("\t",$_);

                if($cnt == 1){
                        $prev_blast_hit_id = $id;
                        $cur_blast_hit_id = $id;

                        if($ss > $se){
                                next;
                        } else{
                                $CUR_HIT_IS_CHECKED = 'TRUE';
                                print FILTERED_FH $_."\n";
                        }
                } else{

                        $prev_blast_hit_id = $cur_blast_hit_id;
                        $cur_blast_hit_id = $id;

                        if($prev_blast_hit_id ne $cur_blast_hit_id){

                                $CUR_HIT_IS_CHECKED = "FALSE";
                                if($ss > $se){
                                        next;
                                } else{
                                       $CUR_HIT_IS_CHECKED = 'TRUE';
                                        print FILTERED_FH $_."\n";
                                }

                        } else{
                                if($CUR_HIT_IS_CHECKED eq "TRUE"){
                                        next;
                                } else{
                                        if($ss > $se){
                                                next;
                                        } else{
                                                $CUR_HIT_IS_CHECKED = 'TRUE';
                                                print FILTERED_FH $_."\n";
                                        }
                                }
                        }
                }

        }
        close(BEST_HIT_PROC);
        close(FILTERED_FH);

        system("rm $blastOutputFile");
        system("mv $filteredBlastOutputFile $blastOutputFile");

}

sub get_best_blast_hit_from_single_miRNA_BLAST{

	# get the first best hit from a single miRNA BLAST that complies with the condition: $ss > $se

	my $tmpBlastOutputFile = shift;


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
				close(FULL_SEQ_PROC);
				system("rm $tmpBlastOutputFile");
                                return($_);
                        }
                } else{

			if($ss > $se){
				next;
			} else{
				close(FULL_SEQ_PROC);
				system("rm $tmpBlastOutputFile");
				return($_);
			}
                }

        }
        close(FULL_SEQ_PROC);
	system("rm $tmpBlastOutputFile");

}


sub get_correct_id_with_arm_info{

	my $match = shift;
	my $arm = shift;

	my $subject_match = $match."-$arm";

	# check if the $subject_match is a valid mature mir id
	my $FOUND_MID_HIT = checkIfMatureMirIdIsValid($subject_match);

	if($FOUND_MID_HIT eq 'TRUE'){
		$match=$subject_match;
	} else{ 

		#print "subject_match: $subject_match\n";				
		my @match_vals = split('-', $match);

		#print "match_vals length: $#match_vals\n";		

		if(@match_vals > 3){ #remove loci info from mir id	

			#print "in match_vals > 3\n";

			my $last_dash_index = rindex($match, '-');

			$match = substr($match, 0, $last_dash_index);
			#print "match: $match\n";
			
			my $new_subject_match=$match."-$arm";
			#print "match: $new_subject_match\n";

			my $VALID_ID_CHECK = checkIfMatureMirIdIsValid($new_subject_match);
			if($VALID_ID_CHECK eq 'TRUE'){
				$match=$new_subject_match;
			} else{ # final verification step
				my $NEW_VALID_ID_CHECK = checkIfMatureMirIdIsValid($match);
			
				if($NEW_VALID_ID_CHECK eq 'FALSE'){

					if($arm eq '5p') { $arm = '3p';}
					elsif($arm eq '3p') { $arm = '5p';}

					$match = $match."-$arm";
					
					# remove arm info instead of the loci info
					my $LOCICOMPLARM_VALID_ID_CHECK = checkIfMatureMirIdIsValid($match);
									
					if($LOCICOMPLARM_VALID_ID_CHECK eq 'FALSE'){
						print "[LOCICOMPLARM_VALID_ID_CHECK]: mir id $match cannot be identified\n";
						next;
						#exit;
					}
				}
			}
		} else{ # remove redundant-non valid 5p/3p info   
			#print "With arm match: $match\n";
			#$match = substr($match, 0, rindex($match, '-'));

			#print "No arm match: $match\n";
			# final verification step
			my $NOARM_VALID_ID_CHECK = checkIfMatureMirIdIsValid($match);
			
			if($NOARM_VALID_ID_CHECK eq 'FALSE'){
				#try if there is a mir with the complementary arm instead
				if($arm eq '5p') { $arm = '3p';}
				elsif($arm eq '3p') { $arm = '5p';}

				$match = $match."-$arm";

				my $COMPLARM_VALID_ID_CHECK = checkIfMatureMirIdIsValid($match);

				if($COMPLARM_VALID_ID_CHECK eq 'FALSE'){
					print "[COMPLARM_VALID_ID_CHECK]: mir id $match cannot be identified\n";
					next;
					#exit;
				}
			}
		}
	}
	return $match;
}


sub getAllValidMatureIds{

        chdir $base_script_dir;
	#chdir '.';

        my %valid_mature_ids = ();

	my $mature_ids_file = "./matureToUnivIDMappings.txt";
        open MATURE_IDS_FH, $mature_ids_file or die $!;

        my $lnCnt = 0;
        while(my $line = <MATURE_IDS_FH>){

                if($lnCnt !=0){
                        my @vals = split('\t', $line);
                        my $mid = $vals[0];
                        $mid =~ s/miR/mir/;

                        $valid_mature_ids{$mid} = 1;
                }

                $lnCnt++;
        }
        close(MATURE_IDS_FH);

        return %valid_mature_ids;
}


sub getAllValidMatureIdsFromJsonFile{
	my $json;

	my $BREAK_JSON_READ_WHILE_LOOP = 'TRUE';
	my $cnt = 0;
	while(1){
		# probably useless when this script is not running as the main thread
		try{
			$json = read_file('./valid_mature_ids_dump.json', { binmode => ':raw' });
			if(defined($json)){
				print "found json after $cnt unsuccessful trials!!\n";	
				$BREAK_JSON_READ_WHILE_LOOP = 'TRUE';
				last;
		 	} 
		} catch{	
			$cnt++;
			"Exception caught while trying to read .json file: $_\n";
			$BREAK_JSON_READ_WHILE_LOOP = 'FALSE';
			
			sleep 0.5;
			continue;	
		}
	}
    	%valid_mature_ids = %{ decode_json $json };

	return(\%valid_mature_ids);
}



# Collect the information about the total number of reads
# for each idientified miRNA (with all types of any modifications) 
# before starting the classification process for the 5p and 3p stretches
# (modifications | SNPs | ADAR-editing)
sub get_totalreads_per_mirna{

	my $blastOutputFile = shift;
	my $filename = shift;

	my %totalreads_per_mirna = ();

	my $cnt = 0;
	open(PROC, '<', $blastOutputFile);
	while(<PROC>){

		$cnt++;
		chomp();
		# qs, qe: query start, query end
		# ss, se: subject start, subject end
		# mm: alignment length (?)
		my ($id,$match,$perc,$mm,$mismatch,$gaps,$qs,$qe,$ss,$se,$ev)=split("\t",$_);

		print "$_\n";

		my $fail=0;
		my $sample_reads++;

		my $nice_id = $id;
		$nice_id =~ s/_x(\d+)/ Depth:$1 Modification:/g;
		my $this_depth=$1;

		print "this_depth: $this_depth\n";

		$allreads{$filename}++;
		$alldepths{$filename}+=$this_depth;

		my $mirlength=length($mirseqs{$id});

		# Remove 'edge effect' hits where the miR maps off the boundary of a precursor
		if (($ss-$qs) < 0){
			print "[FAIL]: ss < qs\n";
			$fail=1;
			print "[DBG]: critical error => (ss < qs)!\n";
		}

		if ($edits){

			if (($ss+($mirlength - $qs)) > $hairpinlength{$match}) {
				if($DEBUG_MODE eq 'TRUE'){
					print "mirlength: $mirlength, match: $match, hairpinlength{$match}: $hairpinlength{$match}\n";
					print "[FAIL]: ss+(mirlength - qs)) > hairpinlength{match}\n";
				}
				$fail=1;
			}

			#Exclude low depth hits
			if ($this_depth < $depth_threshold){
				if($DEBUG_MODE eq 'TRUE'){
					print "[FAIL]: this_depth < depth_threshold\n";
				}
				$fail=1;
			}
		}
			
		if (!$fail){

			my $cur_reads_multiplier = 1;			
			my $arm="";

			if ((($ss+$se)/2) <= ($hairpinlength{$match}/2)){
				$arm="5p";
			} else {
				$arm="3p";
			}

			$match = get_correct_id_with_arm_info($match, $arm);

			my $five_nt="";
			my $three_nt="";	
			if ($mm < $mirlength){
				if ($qs != '1'){
					$five_nt=substr($mirseqs{$id},0,$qs-1);
				}
				if ($qe < $mirlength){
					$three_nt=substr($mirseqs{$id},$qe,$mirlength-$qe);
				}
			}


			if($CONSIDER_READS_WITH_BOTH_MODS_TWICE_IN_TOTAL_READS eq 'TRUE'){
				# consider the reads with both 3p and 5p mods twice in the total reads number.	
				if($five_nt ne "" & $three_nt ne ""){
					$cur_reads_multiplier = 2;
				}
			}

	
			if(!defined($totalreads_per_mirna{$match})){
				$totalreads_per_mirna{$match} = $this_depth * $cur_reads_multiplier;
			} else{
				$totalreads_per_mirna{$match} += $this_depth * $cur_reads_multiplier; 
			}
		}
	}
	
	return(\%totalreads_per_mirna);
}


################################################################
############      Stretches classification module ##############
################################################################
sub classify_stretch{
	
	my $stretch = shift;
	my $match = shift;
	my $initial_id = shift;
	my $mismatch = shift;
	my $gaps = shift;
	my $qs = shift;
	my $ss = shift;
	my $qe = shift;
	my $se = shift;
	my $armEnd = shift;
	my $cur_mir_ratio = shift;
	my $STRETCH_BEYOND_MIR_LIMITS = shift;


	my $stretch_len = length($stretch);
	my $cur_hairpin_length = length($hairpinseq{$initial_id});
	my $hairpin_substr_start = -1;



	if($armEnd eq "5p"){
		$hairpin_substr_start = $ss - $stretch_len - 1;
		
		if($hairpin_substr_start < 0){
			$hairpin_substr_start = 0;
			$stretch_len = $ss;	
		}	
	} else{ #3p
		$hairpin_substr_start = $se;
		
		if( ($hairpin_substr_start + $stretch_len) > $cur_hairpin_length){ 
			$stretch_len = $cur_hairpin_length - $se + 1;	
		}
	}

	
	# for '3p' check if $hairpin_substr_start > length(hairpinseq{$initial_id} 

	my $hairpin_substr = substr($hairpinseq{$initial_id}, $hairpin_substr_start, $stretch_len);


	if($VERBOSE eq 'TRUE'){	
		print "###############################\n";	
		print "armEnd: $armEnd\n";
		print "stretch: $stretch\n";
		print "hairpin_substr: $hairpin_substr\n";
		print "full hairpin: $hairpinseq{$initial_id}\n";
		print "cur_mir_ratio: $cur_mir_ratio\n";						
		print "###############################\n\n";	
	}

	

	my $seqs_diff_array_ref = get_diff_array_from_two_short_sequences($stretch, $hairpin_substr, $armEnd);
	my @seqs_diff_array = @$seqs_diff_array_ref;

	print "seqs_diff_array_ref: \n";
	print "@seqs_diff_array\n";
	
	# characterization will be: 'MOD', 'ADAR' or 'SNP'
	# based on the ratio of the current miRNA read counts
	my ($stretch_characterization_type, $adar_indexes_ref, $potential_snp_indexes_ref) = assign_stretch_characterization(\@seqs_diff_array, $cur_mir_ratio, $id, $qs, $qe, $armEnd, $stretch_len, $STRETCH_BEYOND_MIR_LIMITS); 

	my @adar_indexes = @$adar_indexes_ref;
	my @potential_snp_indexes = @$potential_snp_indexes_ref;


	my $stretch_final_identifier = '';

	
	if($stretch_characterization_type eq 'MOD'){	
		$stretch_final_identifier = $mod_types_separator.'nont_'.$armEnd."_$stretch";

	} elsif($stretch_characterization_type eq 'ADAR'){

		my $adar_offset;
		if(@adar_indexes){
			print "adar indexes: @adar_indexes\n";
			$adar_offset = $adar_indexes[0];
		}

		my $adar_index;
		if($armEnd eq "5p"){
			$adar_index = $qs - $adar_offset;
		} elsif($armEnd eq "3p"){
			$adar_index = $qe + $adar_offset;
		}



		if($ADD_POSITION_OF_SNP_OR_ADAR_INFO_IN_FINAL_IDENTIFIER eq 'TRUE'){
			my $tmp_stretch = '';
				
			for(my $sIdx=0; $sIdx<=$#seqs_diff_array; $sIdx++){
				if($seqs_diff_array[$sIdx] == 1){
					$tmp_stretch .= 'G';
					# by pass the construction of the full stretch with 'x'
					last;
				} else{	
					$tmp_stretch .= 'x';
				}
			}

			$stretch = $tmp_stretch;		
		}

		#$stretch_final_identifier = $mod_types_separator.'adar_'.$armEnd."_$stretch";
		$stretch_final_identifier = $mod_types_separator.'adar'."_$armEnd"."_$stretch";
	
	} elsif($stretch_characterization_type eq 'SNP'){
	
		my $snp_offset;
		if(@potential_snp_indexes){
			print "potential snp indexes: @potential_snp_indexes\n";
			$snp_offset = $potential_snp_indexes[0];
		}
	
		my $snp_index;
		if($armEnd eq "5p"){
			$snp_index = $qs - $snp_offset;
		} elsif($armEnd eq "3p"){
			$snp_index = $qe + $snp_offset;
		}


		if($ADD_POSITION_OF_SNP_OR_ADAR_INFO_IN_FINAL_IDENTIFIER eq 'TRUE'){
			my $tmp_stretch = '';
			
			for(my $sIdx=0; $sIdx<=$#seqs_diff_array; $sIdx++){
				if($seqs_diff_array[$sIdx] == -1){
					$tmp_stretch .= substr($stretch, $sIdx, 1);
					# by pass the construction of the full stretch with 'x'
					last;
				} else{	
					$tmp_stretch .= 'x';
				}
			}

			$stretch = $tmp_stretch;
		}	

		#$stretch_final_identifier = $mod_types_separator.'snp_'.$armEnd."_$stretch"."_$snp_index";
		$stretch_final_identifier = $mod_types_separator.'snp'."_$armEnd"."_$stretch";
	}


	return($stretch_final_identifier);	
}

sub assign_stretch_characterization{
	
	my $seqs_diff_array_ref = shift;
	my $cur_mir_ratio = shift;
	my $id = shift;
	my $qs = shift; 
	my $qe = shift;
	my $armEnd = shift;
	my $stretch_len = shift;
	my $STRETCH_BEYOND_MIR_LIMITS = shift;
	


	my $cur_read_len = length($mirseqs{$id});

	my @seqs_diff_array = @$seqs_diff_array_ref;
	my $stretch_characterization_type = ''; # possible values: MOD, ADAR or SNP

	my $adar_favorable_score = 0;
	my $snp_or_mod_favorable_score = 0;

	my @adar_indexes = ();
	my @potential_snp_indexes = ();

	my $tmp_index = 1;
	foreach(@seqs_diff_array){
		my $curDiffValue = $_;

		if($curDiffValue == 1){
			@adar_indexes = push(@adar_indexes, $tmp_index);
			$adar_favorable_score++;

		} elsif($curDiffValue == -1){
			@potential_snp_indexes = push(@potential_snp_indexes, $tmp_index);
			$snp_or_mod_favorable_score++;
		}

		$tmp_index++;
	}

	if($snp_or_mod_favorable_score == 0){
		
		if($adar_favorable_score > 0){
			$stretch_characterization_type = "ADAR";
		} else {
			print "[error]: no ADAR or SNP/MOD found in current stretch?!";
			exit;		
		}

	} elsif($snp_or_mod_favorable_score == 1){ # in order to have SNP we need only one recorded change
		
		if($STRETCH_BEYOND_MIR_LIMITS ne 'TRUE' & $stretch_len > $STRETCH_LENGTH_THRESHOLD_FOR_SNPS){ #set stricter criteria for identifying SNPs in stretches
			if($cur_mir_ratio > $SNP_RATIO_THRESHOLD){
				$stretch_characterization_type = "SNP";
			} else {
				$stretch_characterization_type = "MOD";
			}
		} else{
			$stretch_characterization_type = "MOD";
		}
	} else {
		$stretch_characterization_type = "MOD";
	}
	

	return ($stretch_characterization_type, \@adar_indexes, \@potential_snp_indexes);


}


sub get_diff_array_from_two_short_sequences {

	# Values mapping:
	# 1: ADAR editing
	# -1: SNP or MOD
	# 0: no change
	
	my $stretch = shift;
	my $hairpin_substr = shift;
	my $armEnd = shift;

	print "stretch: $stretch\n";
	print "hairpin_substr: $hairpin_substr\n";

	# For a 5p stretch ==> we reverse the sequences to do a 3' -> 5' check
	if($armEnd eq "5p"){
		$stretch = reverse($stretch);
		$hairpin_substr = reverse($hairpin_substr);
	}

	my $total_num_of_mismatches = 0;	
	# use the shortest length in order to get the diff between the sequences
	my $length_to_check = (length($stretch) < length($hairpin_substr)) ? length($stretch) : length($hairpin_substr);	

	if($VERBOSE eq 'TRUE'){
		print "(reversed or not) stretch: $stretch\n";
		print "(reversed or not) hairpin_substr: $hairpin_substr\n";
		print "length_to_check: $length_to_check\n";
	}

	my @seqs_diff_array = ();
	for(my $offset=0; $offset<$length_to_check; $offset++){
		if(substr($stretch, $offset, 1) eq substr($hairpin_substr, $offset, 1)){
			$seqs_diff_array[$offset] = 0;  # no change
		} elsif(substr($stretch, $offset, 1) eq 'G' & substr($hairpin_substr, $offset, 1) eq 'A'){
			$seqs_diff_array[$offset] = 1; # ADAR-editing candidate
		} else{
			$seqs_diff_array[$offset] = -1; # SNP or MOD candidate
		}
	}

	# reverse seqs_diff_array for '5p' arm
	if($armEnd eq "5p"){
		@seqs_diff_array = reverse @seqs_diff_array;
	}

	return(\@seqs_diff_array);
}


sub get_all_all_mature_mirnas_align_indexes{

	my $all_indexes_filename = "indexesOfAllMatureMirnaAlignments.txt";
	open(INDEXES_FH, '<', $all_indexes_filename);

	my $line = 0;
	while(<INDEXES_FH>){
		if($line != 0){ 
			(my $mirId, my $hairpinId, my $sIdx, my $eIdx) = split('\t', $_);
			my @indexes = ($sIdx, $eIdx);
			$all_mature_mirnas_align_indexes{$mirId}{$hairpinId} = [@indexes];
		}
		$line++;
	}
}

if($RETAIN_SEQS_FOR_DEBUG_DATASET eq 'TRUE'){
	close(TMP_DEBUG_FH);
	
	my $debug_dataset_zipped_filename = $debug_dataset_filename."gz";
	if(-e $debug_dataset_zipped_filename){
		system("rm $debug_dataset_zipped_filename");
	} 
	system("gzip $debug_dataset_filename");
}

1;
