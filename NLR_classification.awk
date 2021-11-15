#Adapted from awk script from https://github.com/peritob/Myrtaceae_NLR_workflow
#Use for NLR classification step in the snakemake pipeline
#This script can be used individually by running 'bash run_classification.awk {interproscan.tsv} {augustus.gff3}'

#Identify NBARC
# this is how it works: awk takes a tab delimited file and looks in colum 6 for "NB-ARC domain" and removes the ".t1" from gene id in column 1. Then sorts and prints unique gene id into 	temp_$1_0.
gawk 'BEGIN {FS="\t"} $6=="NB-ARC domain" { print $1 }' $1 | sort -k 1b,1 | uniq > ${1%_augustus_aa.fasta.tsv}_NBARC.list
gawk 'BEGIN {FS="\t"} $6=="NB-ARC domain" {split($1, a, "."); print a[1]}' $1 | sort -k 1b,1 | uniq > 	temp_$1_0
# then joins the gene id 	temp_$1_0 with the augustus.gff3 matched from column 9 and does some formatting. then outputs all mathed lines and adds the scaffold location to column 4 and 5. The strand notation for scaffold is transposed into column 7. The logic for this is: the hidden Markov model derived NBARC sequences were identified from the original genome as stranded + or -. This is likely to be the final true strandedness of the gene.  
join 	temp_$1_0 \
  <(gawk 'BEGIN {OFS="\t"} {split($9, a, "[=\\.;]"); print a[2], NR, $0}' $2 | sort -k 1b,1) | \
sort -k2,2 -g | \
awk 'BEGIN {OFS="\t"} {split($3, a, "[:\\-+]"); print  a[1], $4, $5, (a[2]+$6), (a[2]+$7), $8, substr($3,length($3)),  $10, $11}' >${1%_augustus_aa.fasta.tsv}_NBARC.gff3
for files in 	temp_$1_* ; do rm ${files} ; done

#Identify NLR
gawk 'BEGIN {FS="\t"} $5=="PF00931"  { print $1 }' $1 | sort -k 1b,1 | uniq > 	temp_$1_0
gawk 'BEGIN {FS="\t"} ($5=="G3DSA:3.80.10.10" || $5=="PF08263" || $5=="PF13516" || $5=="PF12799" || $5=="PF13306" || $5=="PF13855")  { print $1 }' $1 | sort -k 1b,1 | uniq > 	temp_$1_1
join 	temp_$1_0 	temp_$1_1 > ${1%_augustus_aa.fasta.tsv}_NLR.list
# this script outputs the gff3 file for all genes with NB-ARC and LRR domains in the augustus.gff3
# this is how it works: awk takes a tab delimited file and looks in colum 6 for "NB-ARC domain" OR in column 5 for all LRR Pfam identifiers and removes the ".t1" from gene id in column 1. Then sorts and prints unique gene id into 	temp_$1_0.
# then joins the gene id list 	temp_$1_2 with the augustus.gff3 matched from column 9 and does some formatting. then outputs all matched lines and adds the scaffold location to column 4 and 5. Note that the strand notation for scaffold is transposed into column 7. The logic for this is: the hidden Markov model derived NBARC sequences were identified from the original genome as stranded + or -.
gawk 'BEGIN {FS="\t"} $5=="PF00931" {split($1, a, "."); print a[1]}' $1 | sort -k 1b,1 | uniq > 	temp_$1_0
gawk 'BEGIN {FS="\t"} ($5=="G3DSA:1.10.8.430" || $5=="G3DSA:3.80.10.10" || $5=="PF08263" || $5=="PF13516"  || $5=="PF12799" || $5=="PF13306" || $5=="PF13855") {split($1, a, "."); print a[1]}' $1 | sort -k 1b,1 | uniq > 	temp_$1_1
join 	temp_$1_0 	temp_$1_1 > 	temp_$1_2
join 	temp_$1_2 \
  <(gawk 'BEGIN {OFS="\t"} {split($9, a, "[=\\.;]"); print a[2], NR, $0}' $2 | sort -k 1b,1) | \
sort -k2,2 -g | \
gawk 'BEGIN {OFS="\t"} {split($3, a, "[:\\-+]"); print  a[1], $4, $5, (a[2]+$6), (a[2]+$7), $8, substr($3,length($3)),  $10, $11}' >${1%_augustus_aa.fasta.tsv}_NLR.gff3
for files in 	temp_$1_* ; do rm ${files} ; done

#Identify RNL
# this script outputs the gff3 file for all genes with NB-ARC and RPW8 domains in the augustus.gff3
# this is how it works: awk takes a tab delimited file and looks in colum 6 for "NB-ARC domain" OR in column 5 for all RPW8 Pfam identifiers and removes the ".t1" from gene id in column 1. Then sorts and prints unique gene id into 	temp_$1_0.
# then joins the gene id list 	temp_$1_2 with the augustus.gff3 matched from column 9 and does some formatting. then outputs all matched lines and adds the scaffold location to column 4 and 5. Note that the strand notation for scaffold is transposed into column 7. The logic for this is: the hidden Markov model derived NBARC sequences were identified from the original genome as stranded + or -. 
gawk 'BEGIN {FS="\t"} $5=="PF00931" {split($1, a, "."); print a[1]}' $1 | sort -k 1b,1 | uniq > 	temp_$1_F
gawk 'BEGIN {FS="\t"} $5=="PF05659" {split($1, a, "."); print a[1]}' $1 | sort -k 1b,1 | uniq > 	temp_$1_G
gawk 'BEGIN {FS="\t"} ($5=="G3DSA:3.80.10.10" || $5=="PF08263" || $5=="PF13516"  || $5=="PF12799" || $5=="PF13306" || $5=="PF13855") {split($1, a, "."); print a[1]}' $1 | sort -k 1b,1 | uniq > 	temp_$1_H
join 	temp_$1_F 	temp_$1_G > 	temp_$1_I
join 	temp_$1_I 	temp_$1_H > 	temp_$1_J
cp 	temp_$1_J ${1%_augustus_aa.fasta.tsv}_RNL.list
join 	temp_$1_J \
  <(gawk 'BEGIN {OFS="\t"} {split($9, a, "[=\\.;]"); print a[2], NR, $0}' $2 | sort -k 1b,1) | \
sort -k2,2 -g | \
gawk 'BEGIN {OFS="\t"} {split($3, a, "[:\\-+]"); print  a[1], $4, $5, (a[2]+$6), (a[2]+$7), $8, substr($3,length($3)),  $10, $11}' >${1%_augustus_aa.fasta.tsv}_RNL.gff3
for files in 	temp_$1_* ; do rm ${files} ; done

#RxNL
# this script outputs the gff3 file for all genes with NB-ARC and Rx domains in the augustus.gff3
# this is how it works: awk takes a tab delimited file and looks in colum 6 for "NB-ARC domain" OR in column 5 for all Rx Pfam identifiers and removes the ".t1" from gene id in column 1. Then sorts and prints unique gene id into 	temp_$1_0.
# then joins the gene id list 	temp_$1_2 with the augustus.gff3 matched from column 9 and does some formatting. then outputs all matched lines and adds the scaffold location to column 4 and 5. Note that the strand notation for scaffold is transposed into column 7. The logic for this is: the hidden Markov model derived NBARC sequences were identified from the original genome as stranded + or -. 
gawk 'BEGIN {FS="\t"} $5=="PF00931" {split($1, a, "."); print a[1]}' $1 | sort -k 1b,1 | uniq > 	temp_$1_K
gawk 'BEGIN {FS="\t"} $5=="PF18052" {split($1, a, "."); print a[1]}' $1 | sort -k 1b,1 | uniq > 	temp_$1_L
gawk 'BEGIN {FS="\t"} ($5=="G3DSA:3.80.10.10" || $5=="PF08263" || $5=="PF13516"  || $5=="PF12799" || $5=="PF13306" || $5=="PF13855") {split($1, a, "."); print a[1]}' $1 | sort -k 1b,1 | uniq > 	temp_$1_M
join 	temp_$1_K 	temp_$1_L > 	temp_$1_N
join 	temp_$1_N 	temp_$1_M > 	temp_$1_O
cp 	temp_$1_O ${1%_augustus_aa.fasta.tsv}_RxNL.list
join 	temp_$1_O \
  <(gawk 'BEGIN {OFS="\t"} {split($9, a, "[=\\.;]"); print a[2], NR, $0}' $2 | sort -k 1b,1) | \
sort -k2,2 -g | \
gawk 'BEGIN {OFS="\t"} {split($3, a, "[:\\-+]"); print  a[1], $4, $5, (a[2]+$6), (a[2]+$7), $8, substr($3,length($3)),  $10, $11}' >${1%_augustus_aa.fasta.tsv}_RxNL.gff3
for files in 	temp_$1_* ; do rm ${files} ; done

#TNL
# this script outputs the gff3 file for all genes with NB-ARC and TIR domains in the augustus.gff3
# this is how it works: awk takes a tab delimited file and looks in colum 6 for "NB-ARC domain" OR in column 5 for all TIR Pfam identifiers and removes the ".t1" from gene id in column 1. Then sorts and prints unique gene id into 	temp_$1_0.
# then joins the gene id list 	temp_$1_2 with the augustus.gff3 matched from column 9 and does some formatting. then outputs all matched lines and adds the scaffold location to column 4 and 5. Note that the strand notation for scaffold is transposed into column 7. The logic for this is: the hidden Markov model derived NBARC sequences were identified from the original genome as stranded + or -. This is likely to be the final true strandedness of the gene rather than the braker output.
gawk 'BEGIN {FS="\t"} $5=="PF00931" {split($1, a, "."); print a[1]}' $1 | sort -k 1b,1 | uniq > 	temp_$1_0
gawk 'BEGIN {FS="\t"} $5=="PF01582" {split($1, a, "."); print a[1]}' $1 | sort -k 1b,1 | uniq > 	temp_$1_1
gawk 'BEGIN {FS="\t"} ($5=="G3DSA:1.10.8.430" || $5=="G3DSA:3.80.10.10" || $5=="PF08263" || $5=="PF13516"  || $5=="PF12799" || $5=="PF13306" || $5=="PF13855") {split($1, a, "."); print a[1]}' $1 | sort -k 1b,1 | uniq > 	temp_$1_2
join 	temp_$1_0 	temp_$1_1 > 	temp_$1_3
join 	temp_$1_3 	temp_$1_2 > 	temp_$1_4
cp 	temp_$1_4 ${1%_augustus_aa.fasta.tsv}_TNL.list
join 	temp_$1_4 \
  <(gawk 'BEGIN {OFS="\t"} {split($9, a, "[=\\.;]"); print a[2], NR, $0}' $2 | sort -k 1b,1) | \
sort -k2,2 -g | \
gawk 'BEGIN {OFS="\t"} {split($3, a, "[:\\-+]"); print  a[1], $4, $5, (a[2]+$6), (a[2]+$7), $8, substr($3,length($3)),  $10, $11}' >${1%_augustus_aa.fasta.tsv}_TNL.gff3
for files in 	temp_$1_* ; do rm ${files} ; done

#BNL
# this script outputs the gff3 file for all genes with NB-ARC and BED domains in the augustus.gff3
# this is how it works: awk takes a tab delimited file and looks in colum 6 for "NB-ARC domain" OR in column 5 for all BED Pfam identifiers and removes the ".t1" from gene id in column 1. Then sorts and prints unique gene id into 	temp_$1_0.
# then joins the gene id list 	temp_$1_2 with the augustus.gff3 matched from column 9 and does some formatting. then outputs all matched lines and adds the scaffold location to column 4 and 5. Note that the strand notation for scaffold is transposed into column 7. The logic for this is: the hidden Markov model derived NBARC sequences were identified from the original genome as stranded + or -. 
gawk 'BEGIN {FS="\t"} $5=="PF00931" {split($1, a, "."); print a[1]}' $1 | sort -k 1b,1 | uniq > 	temp_$1_P
gawk 'BEGIN {FS="\t"} $5=="PF02892" {split($1, a, "."); print a[1]}' $1 | sort -k 1b,1 | uniq > 	temp_$1_Q
gawk 'BEGIN {FS="\t"} ($5=="G3DSA:3.80.10.10" || $5=="PF08263" || $5=="PF13516"  || $5=="PF12799" || $5=="PF13306" || $5=="PF13855") {split($1, a, "."); print a[1]}' $1 | sort -k 1b,1 | uniq > 	temp_$1_R
join 	temp_$1_P 	temp_$1_Q > 	temp_$1_S
join 	temp_$1_R 	temp_$1_S > 	temp_$1_T
cp 	temp_$1_T ${1%_augustus_aa.fasta.tsv}_BNL.list
join 	temp_$1_T \
  <(gawk 'BEGIN {OFS="\t"} {split($9, a, "[=\\.;]"); print a[2], NR, $0}' $2 | sort -k 1b,1) | \
sort -k2,2 -g | \
gawk 'BEGIN {OFS="\t"} {split($3, a, "[:\\-+]"); print  a[1], $4, $5, (a[2]+$6), (a[2]+$7), $8, substr($3,length($3)),  $10, $11}' >${1%_augustus_aa.fasta.tsv}_BNL.gff3
for files in 	temp_$1_* ; do rm ${files} ; done

#CNL
# this script outputs the gff3 file for all genes with NB-ARC and Coils domains in the augustus.gff3
# this is how it works: awk takes a tab delimited file and looks in colum 6 for "NB-ARC domain" OR in column 5 for all Coil Pfam identifiers and removes the ".t1" from gene id in column 1. Then sorts and prints unique gene id into 	temp_$1_0.
# then joins the gene id list 	temp_$1_2 with the augustus.gff3 matched from column 9 and does some formatting. then outputs all matched lines and adds the scaffold location to column 4 and 5. Note that the strand notation for scaffold is transposed into column 7. The logic for this is: the hidden Markov model derived NBARC sequences were identified from the original genome as stranded + or -. 
gawk 'BEGIN {FS="\t"} $5=="PF00931" {split($1, a, "."); print a[1]}' $1 | sort -k 1b,1 | uniq > 	temp_$1_A
gawk 'BEGIN {FS="\t"} $5=="Coil" {split($1, a, "."); print a[1]}' $1 | sort -k 1b,1 | uniq > 	temp_$1_B
gawk 'BEGIN {FS="\t"} ($5=="G3DSA:3.80.10.10" || $5=="PF08263" || $5=="PF13516"  || $5=="PF12799" || $5=="PF13306" || $5=="PF13855") {split($1, a, "."); print a[1]}' $1 | sort -k 1b,1 | uniq > 	temp_$1_C
join 	temp_$1_A 	temp_$1_B > 	temp_$1_D
join 	temp_$1_D 	temp_$1_C > 	temp_$1_E
cp 	temp_$1_E ${1%_augustus_aa.fasta.tsv}_CNL.list
join 	temp_$1_E \
  <(gawk 'BEGIN {OFS="\t"} {split($9, a, "[=\\.;]"); print a[2], NR, $0}' $2 | sort -k 1b,1) | \
sort -k2,2 -g | \
gawk 'BEGIN {OFS="\t"} {split($3, a, "[:\\-+]"); print  a[1], $4, $5, (a[2]+$6), (a[2]+$7), $8, substr($3,length($3)),  $10, $11}' >${1%_augustus_aa.fasta.tsv}_CNL.gff3
for files in 	temp_$1_* ; do rm ${files} ; done

#JNL
# this script outputs the gff3 file for all genes with NB-ARC and Jacalin domains in the augustus.gff3
# this is how it works: awk takes a tab delimited file and looks in colum 6 for "NB-ARC domain" OR in column 5 for all Jacalin Pfam identifiers and removes the ".t1" from gene id in column 1. Then sorts and prints unique gene id into 	temp_$1_0.
# then joins the gene id list 	temp_$1_2 with the augustus.gff3 matched from column 9 and does some formatting. then outputs all matched lines and adds the scaffold location to column 4 and 5. Note that the strand notation for scaffold is transposed into column 7. The logic for this is: the hidden Markov model derived NBARC sequences were identified from the original genome as stranded + or -. 
gawk 'BEGIN {FS="\t"} $5=="PF00931" {split($1, a, "."); print a[1]}' $1 | sort -k 1b,1 | uniq > 	temp_$1_U
gawk 'BEGIN {FS="\t"} $5=="PF01419" {split($1, a, "."); print a[1]}' $1 | sort -k 1b,1 | uniq > 	temp_$1_V
gawk 'BEGIN {FS="\t"} ($5=="G3DSA:3.80.10.10" || $5=="PF08263" || $5=="PF13516"  || $5=="PF12799" || $5=="PF13306" || $5=="PF13855") {split($1, a, "."); print a[1]}' $1 | sort -k 1b,1 | uniq > 	temp_$1_W
join 	temp_$1_U 	temp_$1_V > 	temp_$1_X
join 	temp_$1_X 	temp_$1_W > 	temp_$1_Y
cp 	temp_$1_Y ${1%_augustus_aa.fasta.tsv}_JNL.list
join 	temp_$1_Y \
  <(gawk 'BEGIN {OFS="\t"} {split($9, a, "[=\\.;]"); print a[2], NR, $0}' $2 | sort -k 1b,1) | \
sort -k2,2 -g | \
gawk 'BEGIN {OFS="\t"} {split($3, a, "[:\\-+]"); print  a[1], $4, $5, (a[2]+$6), (a[2]+$7), $8, substr($3,length($3)),  $10, $11}' >${1%_augustus_aa.fasta.tsv}_JNL.gff3
for files in 	temp_$1_* ; do rm ${files} ; done
