###
This work is Chimira's predecessor.
You can find the release version of Chimira here:

[http://wwwdev.ebi.ac.uk/enright-dev/chimira/index.php](http://wwwdev.ebi.ac.uk/enright-dev/chimira/index.php)

## Run example ##
./modification_analysis.pl input_test/hsa-mir-320a_dataset.fa.gz out 0 hsa
####
--------------------

## global indexing of miRNA positions:
     5p-mod...{    mature miRNA   }...3p-mod
         xxxxx|ooooooooooooooooooo|xxxxx     	 
            -1 0123....        -10 1234

## The available modification types are:
1. no_mod
2. nont_3p_PATTERN_INDEX     e.g. 'nont_3p_CACC_-3'
3. nont_5p_PATTERN_INDEX                - " - 
4. snp_PATTERN_SNP-INDEX     e.g. 'snp_3p_C_-2' 
5. adar_PATTERN_ADAR-INDEX   e.g  'adar_5p_G_-1'
6. internal_adar_G_INDEX     e.g. 'internal_adar_G_11'
7. internal_snp_NT_INDEX          'internal_snp_A_6'  # there can be either only internal_adar or internal_snp, otherwise the miRNA is not recognized currently as a valid hit. 

## extra post-identifiers:
- 'doubled'		e.g. mmu-mir-143-3p=>nont_3p_U_1=>doubled
Ignore the counts from these mirmod hits when quantifying for templated miRNAs expression.
The reason is that these counts will be considered from the 'sister' miRNA hit that will have recorded the 5p modification
for the same miRNA molecule.
#! Should take care of that in the initial parsing process in the R analysis.


Each modification type is separated by the '=>' separator from all the other modification types.
Different formations of modification complexes are allowed, e.g.:
'=>nont_3p_CACC_-3=>internal_snp_G_11'
'=>no_mod=>internal_snp_A_6'
'=>snp_3p_C_-2=>internal_adar_G_11'




```perl
print "MIRNA\t";
print "MODIFICATION_TYPE\t"; # - (for no_mod), mod (for nont), snp (for snp) or adar (for adar)
print "MODIFICATION_ARM\t"; # 5p or 3p
print "MODIFICATION_PATTERN\t"; # the modification sequence
print "MODIFICATION_POSITION\t"; # an integer index, positive or negative

print "INTERNAL_MOD_TYPE\t"; # - (if there is not), adar (for internal_adar) or snp (for internal_snp)
print "INTERNAL_MOD_PATTERN\t";  # always 'G' for adar or any of the other 3 nts for SNPs. 
print "INTERNAL_MOD_POSITION\t"; # an integer index, always positive.

print "DOUBLED\t"; # yes/-: if yes, then ignore the counts from that hit in all cases that
# don't refer to quantification between 5p and 3p 'mods' only. In global quantification events,
# # I\ll be using the counts from the 5p mod equivalent hit.
```


###
Mirmod_pipeline specifications:
```
 - Only 1 mismatch is permitted. If we have two, we don't get an alignment hit anyway.
 [-] If I have 0 mismatches:
 	I should look for potential mismatches at the 5p and 3p stretch, if there are any stretches.
 	> For the 5p stretch:
 		start comparing (3' to 5' direction) the stretch with the part of the full miRNA sequence
 		that was not aligned. 
 		- If I have ONLY 1 mismatch in that stretch then keep this information either as:
 			* ADAR editing in that position, if we have G instead of the original A
 			* SNP in that position
 		- Else If I have > 1 mismatches
 			* 5p-modification
	> For the 3p stretch:
		start comparing (5' to 3' direction) the stretch with the part of the full miRNA sequence
		that was not aligned.
		- If I have ONLY 1 mismatch in that stretch then keep this information either as:
                       * ADAR editing in that position, if we have G instead of the original A
                       * SNP in that position
		- Else If I have > 1 mismatches
                       * 3p-modification
 [-] If I have 1 mismatch:
 	(comment): approach for the 1st version.
 	> Consider both 5p and 3p stretches as modifications 
```


####
Run: 
 - 1/0 stands for run from my 'extended seqimp' pipeline / run from the command line
 - $curTaxonName can be either Homo_sapiens or Mus_musculus at this stage 
 - $scriptInput = "$curMirmodInputFolder/*.fa.gz";
 - $scriptOutput = "$curMirmodInputFolder";
./modification_analysis.pl $scriptInput $scriptOutput 0 $curTaxonName

#### create blastable database:
makeblastdb -in mmu_hairpin_precursors.fa -dbtype nucl -title mmu_hairpin_precursors -out ./mmu_hairpin_precursors

#### Takes as input a set of CLEAN UNIQ files from the output directory of Kraken
../../mirmod_pipeline/modification_analysis.pl *.fa.gz

- Generates a table of miRNA edits from all samples for downstream R analysis

#### dvitsios:
*** COMMENTS FOR THE PIPELINE ***

The final output dir contains all dataset dirs that have been processed 
with a processed counts file for each file in each dataset.
All the processed counts within a dataset are later on merged
using an R script, the 'mergeCountsFiles.R' which is called inside
the main getAndProcessDataFromENA.pl script.


