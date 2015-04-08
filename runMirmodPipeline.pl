#!/usr/bin/env perl
#use v5.16.0;
use strict;
use warnings;
use threads;
use threads::shared;
use Data::Dumper;
use Cwd;
use JSON::XS qw(encode_json decode_json);
use File::Slurp qw(read_file write_file);
use Fcntl qw(:flock);

#require "modification_analysis.pl"; 


my $root_input_dir = $ARGV[0];


my $base_run_dir = getcwd;

my $mature_ids_file = "matureToUnivIDMappings.txt";
my %valid_mature_ids = getAllValidMatureIds();
{
    my $json = encode_json \%valid_mature_ids;
    write_file('valid_mature_ids_dump.json', { binmode => ':raw' }, $json);
}



opendir my $dh, $root_input_dir
  or die "$0: opendir: $!";

my @datasets = grep {-d "$root_input_dir/$_" && ! /^\.{1,2}$/} readdir($dh);

my %taxas_for_each_dataset :shared= ();
# get taxon for each dataset:
my $datasets_taxa_file = "DATASET_TAXA_MAPPING.txt";
open my $datasets_taxa_fh, $datasets_taxa_file or die "Could not open $datasets_taxa_file: $!";	

while(my $d_tax_line = <$datasets_taxa_fh>){
		
	chomp $d_tax_line;

	my @d_tax_vals = split(/\t/, $d_tax_line);
	my $curAccessionNum = $d_tax_vals[0];
	my $curTaxonName = $d_tax_vals[1]; 

	$taxas_for_each_dataset{$curAccessionNum} = $curTaxonName;
}

close $datasets_taxa_fh;
#print Dumper(\%taxas_for_each_dataset);


my @datasets_threads = ();
# parallelize it in separate thread for each dataset!
# e.g. PRJNA19003/
foreach(@datasets){

	my $cur_dataset_shared :shared = $_;
	
	#push @datasets_threads, threads->create( sub{
	
	push @datasets_threads, threads->create(\&runMirmodForADatasetThread, $_);
	#sleep 1;	

	sub runMirmodForADatasetThread{

		chdir $base_run_dir;
		#lock($cur_dataset_shared);
		my $cur_dataset = shift;


		print "\n\n$cur_dataset:\n";
		my $cur_dataset_dir = $root_input_dir."/$cur_dataset";

		opendir(sep_results_dh, $cur_dataset_dir) 
			or die "$0: opendir:$! ($cur_dataset_dir)";

		my @samples_dirs = grep {-d "$cur_dataset_dir/$_" && ! /^\.{1,2}$/} readdir(sep_results_dh);
		close sep_results_dh;



		my $CUR_ACCESSION_NUM_IS_VALID_TO_RUN = 'FALSE';
		my $curTaxonName = "";
		lock(%taxas_for_each_dataset);
		if(defined($taxas_for_each_dataset{$cur_dataset})){
			$CUR_ACCESSION_NUM_IS_VALID_TO_RUN = 'TRUE';
			$curTaxonName = $taxas_for_each_dataset{$cur_dataset};
		}	
		#find the $curTaxonName from a pre-prepared list of dataset-taxa associations
	


		if($CUR_ACCESSION_NUM_IS_VALID_TO_RUN eq 'TRUE'){

			my @samples_dirs_threads = ();
	
			#e.g. 1/, 2/, etc
			foreach(@samples_dirs){
				
				my $cur_sample = $_;
				my $cur_sample_dir = $cur_dataset_dir."/$cur_sample/";
		
	
				#print "cur_sample_dir: $cur_sample_dir\n";
#				opendir my $cleaned_zipped_file_dh, $cur_sample_dir
#					 or print $cur_dataset and exit;
#					#or die "$0: opendir:$!";

				my @clean_zipped_files = "*.fa.gz";
				
				
#				while(my $file = readdir($cleaned_zipped_file_dh)) {

#					next unless (-f "$cur_sample_dir/$file");
#					next unless ($file =~ m/\.gz$/);

#					print "$file\n";

#					push(@clean_zipped_files, $file);
#				}
#				close $cleaned_zipped_file_dh;
			

#				if(@clean_zipped_files == 0){
#					print "clean_zipped_files == 0\n";
#					exit;				
#				}


				if(@clean_zipped_files > 1){
					print "more than one .fa.gz files have been detected in $cur_sample_dir\n";
					exit;
				} elsif(@clean_zipped_files == 1){

					my $modScriptInput =  $cur_sample_dir."/$clean_zipped_files[0]";
					my $modScriptOutput = $cur_sample_dir;
					
					push @samples_dirs_threads, threads->create(sub{
						print "Running modification_analysis for $modScriptInput\n";
						# chdir to base dir, where I should be anyway
						chdir $base_run_dir;
						my $system_call_return = 0;
						my $mod_analysis_perl_script = $base_run_dir."/modification_analysis.pl";
						
						my $call_mirmod_cnt = 0;
						while(1){
							$system_call_return = system("perl $mod_analysis_perl_script $modScriptInput $modScriptOutput 0 $curTaxonName");
							if($system_call_return == 0){
								if($call_mirmod_cnt > 0){
									print "Modification_analysis.pl exception was resolved after $call_mirmod_cnt failure(s)!\n";
								}
								#print "Successfully finished mirmod_pipeline after $call_mirmod_cnt failures\n";
								last;
							} else{
								$call_mirmod_cnt++;
								print "system_error: $system_call_return\n";
								sleep 1;
							}
						}
					});
				}
			}

			foreach(@samples_dirs_threads){
 			       $_->join();
			}	

			my $numberOfFilesForRAnalysis = $#samples_dirs + 1;
			my $mirmodDirWithCountsForCurAccessNumPath = $cur_dataset_dir;
			my $curAccessNum = $cur_dataset;

			mergeProcessedCountsFilesForAnAccessNum($numberOfFilesForRAnalysis, $mirmodDirWithCountsForCurAccessNumPath, $curAccessNum, $root_input_dir);
		}
	}
}


foreach(@datasets_threads){
       $_->join();
}


sub mergeProcessedCountsFilesForAnAccessNum{

        my $numberOfFilesForRAnalysis = shift;
        my $mirmodDirWithCountsForCurAccessNumPath = shift;
        my $curAccessNum = shift;
	my $root_input_dir = shift;


        my $mergedCountsFileSavePath =$root_input_dir;

        my $mergeCountsFiles_rScript= "mergeCountsFiles.R";

	# call R script for merging processed.counts files passing the parent dir of the dirs './1/processedc.counts\, './2/processed.counts' etc and the number of dirs to look at
        my $mergeFilesROutput = `Rscript $mergeCountsFiles_rScript $mirmodDirWithCountsForCurAccessNumPath $numberOfFilesForRAnalysis $mergedCountsFileSavePath $curAccessNum`;

        print $? if $?;
        print $mergeFilesROutput,"\n";
}


sub getAllValidMatureIds{

	chdir $base_run_dir;

	my %valid_mature_ids = ();

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
                         
