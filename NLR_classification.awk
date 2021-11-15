#Adapted from awk script from https://github.com/peritob/Myrtaceae_NLR_workflow
#Use for NLR classification step in the snakemake pipeline
#This script can be used individually by running 'bash run_classification.awk {interproscan.tsv} {augustus.gff3}'

#Identify NBARC
# this is how it works: awk takes a tab delimited file and looks in colum 6 for "NB-ARC domain" and removes the ".t1" from gene id in column 1. Then sorts and prints unique gene id into temp0.
gawk 'BEGIN {FS="\t"} $6=="NB-ARC domain" { print $1 }' $1 | sort -k 1b,1 | uniq > ${1%_augustus_aa.fasta.tsv}_NBARC.list
gawk 'BEGIN {FS="\t"} $6=="NB-ARC domain" {split($1, a, "."); print a[1]}' $1 | sort -k 1b,1 | uniq > temp0
# then joins the gene id temp0 with the augustus.gff3 matched from column 9 and does some formatting. then outputs all mathed lines and adds the scaffold location to column 4 and 5. The strand notation for scaffold is transposed into column 7. The logic for this is: the hidden Markov model derived NBARC sequences were identified from the original genome as stranded + or -. This is likely to be the final true strandedness of the gene.  
join temp0 \
  <(gawk 'BEGIN {OFS="\t"} {split($9, a, "[=\\.;]"); print a[2], NR, $0}' $2 | sort -k 1b,1) | \
sort -k2,2 -g | \
awk 'BEGIN {OFS="\t"} {split($3, a, "[:\\-+]"); print  a[1], $4, $5, (a[2]+$6), (a[2]+$7), $8, substr($3,length($3)),  $10, $11}' >${1%_augustus_aa.fasta.tsv}_NBARC.gff3
for files in temp*
  do  rm ${files}
done

#Identify NLR
gawk 'BEGIN {FS="\t"} $5=="PF00931"  { print $1 }' $1 | sort -k 1b,1 | uniq > temp0
gawk 'BEGIN {FS="\t"} ($5=="G3DSA:3.80.10.10" || $5=="PF08263" || $5=="PF13516" || $5=="PF12799" || $5=="PF13306" || $5=="PF13855")  { print $1 }' $1 | sort -k 1b,1 | uniq > temp1
join temp0 temp1 > ${1%_augustus_aa.fasta.tsv}_NLR.list
# this script outputs the gff3 file for all genes with NB-ARC and LRR domains in the augustus.gff3
# this is how it works: awk takes a tab delimited file and looks in colum 6 for "NB-ARC domain" OR in column 5 for all LRR Pfam identifiers and removes the ".t1" from gene id in column 1. Then sorts and prints unique gene id into temp0.
# then joins the gene id list temp2 with the augustus.gff3 matched from column 9 and does some formatting. then outputs all matched lines and adds the scaffold location to column 4 and 5. Note that the strand notation for scaffold is transposed into column 7. The logic for this is: the hidden Markov model derived NBARC sequences were identified from the original genome as stranded + or -.
gawk 'BEGIN {FS="\t"} $5=="PF00931" {split($1, a, "."); print a[1]}' $1 | sort -k 1b,1 | uniq > temp0
gawk 'BEGIN {FS="\t"} ($5=="G3DSA:1.10.8.430" || $5=="G3DSA:3.80.10.10" || $5=="PF08263" || $5=="PF13516"  || $5=="PF12799" || $5=="PF13306" || $5=="PF13855") {split($1, a, "."); print a[1]}' $1 | sort -k 1b,1 | uniq > temp1
join temp0 temp1 > temp2
join temp2 \
  <(gawk 'BEGIN {OFS="\t"} {split($9, a, "[=\\.;]"); print a[2], NR, $0}' $2 | sort -k 1b,1) | \
sort -k2,2 -g | \
gawk 'BEGIN {OFS="\t"} {split($3, a, "[:\\-+]"); print  a[1], $4, $5, (a[2]+$6), (a[2]+$7), $8, substr($3,length($3)),  $10, $11}' >${1%_augustus_aa.fasta.tsv}_NLR.gff3
for files in temp*
  do  rm ${files}
done

#Identify RNL
# this script outputs the gff3 file for all genes with NB-ARC and RPW8 domains in the augustus.gff3
# this is how it works: awk takes a tab delimited file and looks in colum 6 for "NB-ARC domain" OR in column 5 for all RPW8 Pfam identifiers and removes the ".t1" from gene id in column 1. Then sorts and prints unique gene id into temp0.
# then joins the gene id list temp2 with the augustus.gff3 matched from column 9 and does some formatting. then outputs all matched lines and adds the scaffold location to column 4 and 5. Note that the strand notation for scaffold is transposed into column 7. The logic for this is: the hidden Markov model derived NBARC sequences were identified from the original genome as stranded + or -. 
gawk 'BEGIN {FS="\t"} $5=="PF00931" {split($1, a, "."); print a[1]}' $1 | sort -k 1b,1 | uniq > tempF
gawk 'BEGIN {FS="\t"} $5=="PF05659" {split($1, a, "."); print a[1]}' $1 | sort -k 1b,1 | uniq > tempG
gawk 'BEGIN {FS="\t"} ($5=="G3DSA:3.80.10.10" || $5=="PF08263" || $5=="PF13516"  || $5=="PF12799" || $5=="PF13306" || $5=="PF13855") {split($1, a, "."); print a[1]}' $1 | sort -k 1b,1 | uniq > tempH
join tempF tempG > tempI
join tempI tempH > tempJ
cp tempJ ${1%_augustus_aa.fasta.tsv}_RNL.list
join tempJ \
  <(gawk 'BEGIN {OFS="\t"} {split($9, a, "[=\\.;]"); print a[2], NR, $0}' $2 | sort -k 1b,1) | \
sort -k2,2 -g | \
gawk 'BEGIN {OFS="\t"} {split($3, a, "[:\\-+]"); print  a[1], $4, $5, (a[2]+$6), (a[2]+$7), $8, substr($3,length($3)),  $10, $11}' >${1%_augustus_aa.fasta.tsv}_RNL.gff3
for files in temp*
  do  rm ${files}
done

#RxNL
# this script outputs the gff3 file for all genes with NB-ARC and Rx domains in the augustus.gff3
# this is how it works: awk takes a tab delimited file and looks in colum 6 for "NB-ARC domain" OR in column 5 for all Rx Pfam identifiers and removes the ".t1" from gene id in column 1. Then sorts and prints unique gene id into temp0.
# then joins the gene id list temp2 with the augustus.gff3 matched from column 9 and does some formatting. then outputs all matched lines and adds the scaffold location to column 4 and 5. Note that the strand notation for scaffold is transposed into column 7. The logic for this is: the hidden Markov model derived NBARC sequences were identified from the original genome as stranded + or -. 
gawk 'BEGIN {FS="\t"} $5=="PF00931" {split($1, a, "."); print a[1]}' $1 | sort -k 1b,1 | uniq > tempK
gawk 'BEGIN {FS="\t"} $5=="PF18052" {split($1, a, "."); print a[1]}' $1 | sort -k 1b,1 | uniq > tempL
gawk 'BEGIN {FS="\t"} ($5=="G3DSA:3.80.10.10" || $5=="PF08263" || $5=="PF13516"  || $5=="PF12799" || $5=="PF13306" || $5=="PF13855") {split($1, a, "."); print a[1]}' $1 | sort -k 1b,1 | uniq > tempM
join tempK tempL > tempN
join tempN tempM > tempO
cp tempO ${1%_augustus_aa.fasta.tsv}_RxNL.list
join tempO \
  <(gawk 'BEGIN {OFS="\t"} {split($9, a, "[=\\.;]"); print a[2], NR, $0}' $2 | sort -k 1b,1) | \
sort -k2,2 -g | \
gawk 'BEGIN {OFS="\t"} {split($3, a, "[:\\-+]"); print  a[1], $4, $5, (a[2]+$6), (a[2]+$7), $8, substr($3,length($3)),  $10, $11}' >${1%_augustus_aa.fasta.tsv}_RxNL.gff3

#TNL
# this script outputs the gff3 file for all genes with NB-ARC and TIR domains in the augustus.gff3
# this is how it works: awk takes a tab delimited file and looks in colum 6 for "NB-ARC domain" OR in column 5 for all TIR Pfam identifiers and removes the ".t1" from gene id in column 1. Then sorts and prints unique gene id into temp0.
# then joins the gene id list temp2 with the augustus.gff3 matched from column 9 and does some formatting. then outputs all matched lines and adds the scaffold location to column 4 and 5. Note that the strand notation for scaffold is transposed into column 7. The logic for this is: the hidden Markov model derived NBARC sequences were identified from the original genome as stranded + or -. This is likely to be the final true strandedness of the gene rather than the braker output.
gawk 'BEGIN {FS="\t"} $5=="PF00931" {split($1, a, "."); print a[1]}' $1 | sort -k 1b,1 | uniq > temp0
gawk 'BEGIN {FS="\t"} $5=="PF01582" {split($1, a, "."); print a[1]}' $1 | sort -k 1b,1 | uniq > temp1
gawk 'BEGIN {FS="\t"} ($5=="G3DSA:1.10.8.430" || $5=="G3DSA:3.80.10.10" || $5=="PF08263" || $5=="PF13516"  || $5=="PF12799" || $5=="PF13306" || $5=="PF13855") {split($1, a, "."); print a[1]}' $1 | sort -k 1b,1 | uniq > temp2
join temp0 temp1 > temp3
join temp3 temp2 > temp4
cp temp4 ${1%_augustus_aa.fasta.tsv}_TNL.list
join temp4 \
  <(gawk 'BEGIN {OFS="\t"} {split($9, a, "[=\\.;]"); print a[2], NR, $0}' $2 | sort -k 1b,1) | \
sort -k2,2 -g | \
gawk 'BEGIN {OFS="\t"} {split($3, a, "[:\\-+]"); print  a[1], $4, $5, (a[2]+$6), (a[2]+$7), $8, substr($3,length($3)),  $10, $11}' >${1%_augustus_aa.fasta.tsv}_TNL.gff3
for files in temp*
  do  rm ${files}
done

#BNL
# this script outputs the gff3 file for all genes with NB-ARC and BED domains in the augustus.gff3
# this is how it works: awk takes a tab delimited file and looks in colum 6 for "NB-ARC domain" OR in column 5 for all BED Pfam identifiers and removes the ".t1" from gene id in column 1. Then sorts and prints unique gene id into temp0.
# then joins the gene id list temp2 with the augustus.gff3 matched from column 9 and does some formatting. then outputs all matched lines and adds the scaffold location to column 4 and 5. Note that the strand notation for scaffold is transposed into column 7. The logic for this is: the hidden Markov model derived NBARC sequences were identified from the original genome as stranded + or -. 
gawk 'BEGIN {FS="\t"} $5=="PF00931" {split($1, a, "."); print a[1]}' $1 | sort -k 1b,1 | uniq > tempP
gawk 'BEGIN {FS="\t"} $5=="PF02892" {split($1, a, "."); print a[1]}' $1 | sort -k 1b,1 | uniq > tempQ
gawk 'BEGIN {FS="\t"} ($5=="G3DSA:3.80.10.10" || $5=="PF08263" || $5=="PF13516"  || $5=="PF12799" || $5=="PF13306" || $5=="PF13855") {split($1, a, "."); print a[1]}' $1 | sort -k 1b,1 | uniq > tempR
join tempP tempQ > tempS
join tempR tempS > tempT
cp tempT ${1%_augustus_aa.fasta.tsv}_BNL.list
join tempT \
  <(gawk 'BEGIN {OFS="\t"} {split($9, a, "[=\\.;]"); print a[2], NR, $0}' $2 | sort -k 1b,1) | \
sort -k2,2 -g | \
gawk 'BEGIN {OFS="\t"} {split($3, a, "[:\\-+]"); print  a[1], $4, $5, (a[2]+$6), (a[2]+$7), $8, substr($3,length($3)),  $10, $11}' >${1%_augustus_aa.fasta.tsv}_BNL.gff3
for files in temp*
  do  rm ${files}
done

#CNL
# this script outputs the gff3 file for all genes with NB-ARC and Coils domains in the augustus.gff3
# this is how it works: awk takes a tab delimited file and looks in colum 6 for "NB-ARC domain" OR in column 5 for all Coil Pfam identifiers and removes the ".t1" from gene id in column 1. Then sorts and prints unique gene id into temp0.
# then joins the gene id list temp2 with the augustus.gff3 matched from column 9 and does some formatting. then outputs all matched lines and adds the scaffold location to column 4 and 5. Note that the strand notation for scaffold is transposed into column 7. The logic for this is: the hidden Markov model derived NBARC sequences were identified from the original genome as stranded + or -. 
gawk 'BEGIN {FS="\t"} $5=="PF00931" {split($1, a, "."); print a[1]}' $1 | sort -k 1b,1 | uniq > tempA
gawk 'BEGIN {FS="\t"} $5=="Coil" {split($1, a, "."); print a[1]}' $1 | sort -k 1b,1 | uniq > tempB
gawk 'BEGIN {FS="\t"} ($5=="G3DSA:3.80.10.10" || $5=="PF08263" || $5=="PF13516"  || $5=="PF12799" || $5=="PF13306" || $5=="PF13855") {split($1, a, "."); print a[1]}' $1 | sort -k 1b,1 | uniq > tempC
join tempA tempB > tempD
join tempD tempC > tempE
cp tempE ${1%_augustus_aa.fasta.tsv}_CNL.list
join tempE \
  <(gawk 'BEGIN {OFS="\t"} {split($9, a, "[=\\.;]"); print a[2], NR, $0}' $2 | sort -k 1b,1) | \
sort -k2,2 -g | \
gawk 'BEGIN {OFS="\t"} {split($3, a, "[:\\-+]"); print  a[1], $4, $5, (a[2]+$6), (a[2]+$7), $8, substr($3,length($3)),  $10, $11}' >${1%_augustus_aa.fasta.tsv}_CNL.gff3
for files in temp*
  do  rm ${files}
done

#JNL
# this script outputs the gff3 file for all genes with NB-ARC and Jacalin domains in the augustus.gff3
# this is how it works: awk takes a tab delimited file and looks in colum 6 for "NB-ARC domain" OR in column 5 for all Jacalin Pfam identifiers and removes the ".t1" from gene id in column 1. Then sorts and prints unique gene id into temp0.
# then joins the gene id list temp2 with the augustus.gff3 matched from column 9 and does some formatting. then outputs all matched lines and adds the scaffold location to column 4 and 5. Note that the strand notation for scaffold is transposed into column 7. The logic for this is: the hidden Markov model derived NBARC sequences were identified from the original genome as stranded + or -. 
gawk 'BEGIN {FS="\t"} $5=="PF00931" {split($1, a, "."); print a[1]}' $1 | sort -k 1b,1 | uniq > tempU
gawk 'BEGIN {FS="\t"} $5=="PF01419" {split($1, a, "."); print a[1]}' $1 | sort -k 1b,1 | uniq > tempV
gawk 'BEGIN {FS="\t"} ($5=="G3DSA:3.80.10.10" || $5=="PF08263" || $5=="PF13516"  || $5=="PF12799" || $5=="PF13306" || $5=="PF13855") {split($1, a, "."); print a[1]}' $1 | sort -k 1b,1 | uniq > tempW
join tempU tempV > tempX
join tempX tempW > tempY
cp tempY ${1%_augustus_aa.fasta.tsv}_JNL.list
join tempY \
  <(gawk 'BEGIN {OFS="\t"} {split($9, a, "[=\\.;]"); print a[2], NR, $0}' $2 | sort -k 1b,1) | \
sort -k2,2 -g | \
gawk 'BEGIN {OFS="\t"} {split($3, a, "[:\\-+]"); print  a[1], $4, $5, (a[2]+$6), (a[2]+$7), $8, substr($3,length($3)),  $10, $11}' >${1%_augustus_aa.fasta.tsv}_JNL.gff3
for files in temp*
  do  rm ${files}
done
