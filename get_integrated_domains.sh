#!/usr/bin/env bash
#Copy from https://github.com/peritob/Myrtaceae_NLR_workflow/blob/master/get_integrated_domains.sh

prefix_name=$1
NLR_type=$2

#First extract all Pfam codes from the .tsv - thereby excluding all coils and other database outputs.
grep 'Pfam' result/${prefix_name}_braker_aa.fasta.tsv > result/${prefix_name}_pfam.tsv
#Then extract all the genes based on the lists from NLR_clasification. For example,
grep -F -w -f result/${prefix_name}_braker_aa_${NLR_type}.list result/${prefix_name}_pfam.tsv > result/${prefix_name}_${NLR_type}_Pfams.tsv
#Then extract all domains that are not present in NLR classifications using the parameter -v with fgrep. Make a list of NLR Pfams and call the txt file "NLR_Pfams.list" (PF00931, PF18052, PF08263, PF12799, PF13306, PF13855, PF13516, PF00560, PF07725, PF05659, PF01582, PF17862, PF23598, PF23247) associated with NLR genes.
grep -F -wv -f NLR_Pfams.list result/${prefix_name}_${NLR_type}_Pfams.tsv > result/${prefix_name}_${NLR_type}_ID_Pfams.tsv
#Obtain meaningful column data 
gawk 'BEGIN {FS="\t"; OFS="\t"}; {print $1, $5, $6, $7, $8}' result/${prefix_name}_${NLR_type}_ID_Pfams.tsv > result/${prefix_name}_${NLR_type}_ID_Pfams.data
