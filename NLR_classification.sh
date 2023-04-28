#!/usr/bin/env bash

#From a bash script from https://github.com/peritob/Myrtaceae_NLR_workflow

#Use for NLR classification step in the snakemake pipeline, need seqtk and gff3tools, can be installed with conda inastall -c bioconda seqtk -y and pip install gff3tool
#gff3tools is from https://github.com/NAL-i5K/GFF3toolkit
#This script can be used individually by running 'bash run_classification.sh {interproscan.tsv} {augustus.gff3} {corresponding fasta file}'
TSV_FILE=$1
GFF3_FILE=$2
PROTEIN_FILE=$3
GENOME_FILE=$4

# try to get a sensible prefix for the output files
OUT_PREFIX=${TSV_FILE%.tsv}
OUT_PREFIX=${OUT_PREFIX%.fasta}
OUT_PREFIX=${OUT_PREFIX%_augustus_aa}

echo tsv file: $TSV_FILE
echo gff3 file: $GFF3_FILE
echo protein file: $PROTEIN_FILE
echo genome file: $GENOME_FILE
echo prefix: $OUT_PREFIX

# Build a temporary working file which is the gff3 file but with two columns inserted on the left.
# The first column contains the gene ID without the extension (extracted from column 9).
# The second column contains the row number so we can return to the original gff3 order.
gawk 'BEGIN {FS="\t"; OFS="\t"} {split($9, a, "[=\\.;]"); print a[2], NR, $0}' $GFF3_FILE | sort -k 1b,1 > GFF3_temp


#
# Identify NBARC
#
gawk 'BEGIN {FS="\t"} $5=="PF00931" {print $1}' $TSV_FILE | sort -k 1b,1 | uniq > ${OUT_PREFIX}_NBARC.list
gawk '{split($1, a, "."); print a[1]}' ${OUT_PREFIX}_NBARC.list | sort -k 1b,1 | uniq > NBARC_temp_1
join NBARC_temp_1 GFF3_temp | sort -k2,2 -n > NBARC_temp_2
# Reformat rows to produce new gff3 file.  This outputs the annotated gene co-ordinates and orientation based on the original genome using braker outputs.
gawk 'BEGIN {OFS="\t"} {split($3, a, "[:\\-+]"); print  a[1], $4, $5, (a[2]+$6), (a[2]+$7), $8, $9, $10, $11}' NBARC_temp_2 > ${OUT_PREFIX}_NBARC.gff3
seqtk subseq ${PROTEIN_FILE} ${OUT_PREFIX}_NBARC.list >${OUT_PREFIX}_NBARC.fasta
gff3_to_fasta -g ${OUT_PREFIX}_NBARC.gff3 -f ${GENOME_FILE} -st cds -noQC -d simple -o ${OUT_PREFIX}_NBARC
for files in NBARC_temp_*; do rm ${files}; done


#
# Identify NLR
#
gawk 'BEGIN {FS="\t"} ($5=="G3DSA:3.80.10.10" || $5=="PF08263" || $5=="PF13516" || $5=="PF12799" || $5=="PF13306" || $5=="PF13855") {print $1}' $TSV_FILE | sort -k 1b,1 | uniq > NLR_temp_1
join ${OUT_PREFIX}_NBARC.list NLR_temp_1 > ${OUT_PREFIX}_NLR.list
gawk '{split($1, a, "."); print a[1]}' ${OUT_PREFIX}_NLR.list | sort -k 1b,1 | uniq > NLR_temp_2
join NLR_temp_2 GFF3_temp | sort -k2,2 -n > NLR_temp_3
# Reformat rows to produce new gff3 file
gawk 'BEGIN {OFS="\t"} {split($3, a, "[:\\-+]"); print  a[1], $4, $5, (a[2]+$6), (a[2]+$7), $8, $9, $10, $11}' NLR_temp_3 > ${OUT_PREFIX}_NLR.gff3
seqtk subseq ${PROTEIN_FILE} ${OUT_PREFIX}_NLR.list >${OUT_PREFIX}_NLR.fasta
gff3_to_fasta -g ${OUT_PREFIX}_NLR.gff3 -f ${GENOME_FILE} -st cds -noQC -d simple -o ${OUT_PREFIX}_NLR
for files in NLR_temp_*; do rm ${files}; done


#
# Identify RNB
#
gawk 'BEGIN {FS="\t"} $5=="PF05659" {print $1}' $TSV_FILE | sort -k 1b,1 | uniq > RNB_temp_1
join ${OUT_PREFIX}_NBARC.list RNB_temp_1 > ${OUT_PREFIX}_RNB.list
gawk '{split($1, a, "."); print a[1]}' ${OUT_PREFIX}_RNB.list | sort -k 1b,1 | uniq > RNB_temp_2
join RNB_temp_2 GFF3_temp | sort -k2,2 -n > RNB_temp_3
# Reformat rows to produce new gff3 file
gawk 'BEGIN {OFS="\t"} {split($3, a, "[:\\-+]"); print  a[1], $4, $5, (a[2]+$6), (a[2]+$7), $8, $9, $10, $11}' RNB_temp_3 > ${OUT_PREFIX}_RNB.gff3
seqtk subseq ${PROTEIN_FILE} ${OUT_PREFIX}_RNB.list >${OUT_PREFIX}_RNB.fasta
gff3_to_fasta -g ${OUT_PREFIX}_RNB.gff3 -f ${GENOME_FILE} -st cds -noQC -d simple -o ${OUT_PREFIX}_RNB
for files in RNB_temp_* ; do rm ${files} ; done

#
# Identify RNL
#
gawk 'BEGIN {FS="\t"} $5=="PF05659" {print $1}' $TSV_FILE | sort -k 1b,1 | uniq > RNL_temp_1
join ${OUT_PREFIX}_NLR.list RNL_temp_1 > ${OUT_PREFIX}_RNL.list
gawk '{split($1, a, "."); print a[1]}' ${OUT_PREFIX}_RNL.list | sort -k 1b,1 | uniq > RNL_temp_2
join RNL_temp_2 GFF3_temp | sort -k2,2 -n > RNL_temp_3
# Reformat rows to produce new gff3 file
gawk 'BEGIN {OFS="\t"} {split($3, a, "[:\\-+]"); print  a[1], $4, $5, (a[2]+$6), (a[2]+$7), $8, $9, $10, $11}' RNL_temp_3 > ${OUT_PREFIX}_RNL.gff3
seqtk subseq ${PROTEIN_FILE} ${OUT_PREFIX}_RNL.list >${OUT_PREFIX}_RNL.fasta
gff3_to_fasta -g ${OUT_PREFIX}_RNL.gff3 -f ${GENOME_FILE} -st cds -noQC -d simple -o ${OUT_PREFIX}_RNL
for files in RNL_temp_* ; do rm ${files} ; done


#
# Identify RxNB
#
gawk 'BEGIN {FS="\t"} $5=="PF18052" {print $1}' $TSV_FILE | sort -k 1b,1 | uniq > RxNB_temp_1
join ${OUT_PREFIX}_NBARC.list RxNB_temp_1 > ${OUT_PREFIX}_RxNB.list
gawk '{split($1, a, "."); print a[1]}' ${OUT_PREFIX}_RxNB.list | sort -k 1b,1 | uniq > RxNB_temp_2
join RxNB_temp_2 GFF3_temp | sort -k2,2 -n > RxNB_temp_3
# Reformat rows to produce new gff3 file
gawk 'BEGIN {OFS="\t"} {split($3, a, "[:\\-+]"); print  a[1], $4, $5, (a[2]+$6), (a[2]+$7), $8, $9, $10, $11}' RxNB_temp_3 > ${OUT_PREFIX}_RxNB.gff3
seqtk subseq ${PROTEIN_FILE} ${OUT_PREFIX}_RxNB.list >${OUT_PREFIX}_RxNB.fasta
gff3_to_fasta -g ${OUT_PREFIX}_RxNB.gff3 -f ${GENOME_FILE} -st cds -noQC -d simple -o ${OUT_PREFIX}_RxNB
for files in RxNB_temp_* ; do rm ${files} ; done


#
# Identify RxNL
#
gawk 'BEGIN {FS="\t"} $5=="PF18052" {print $1}' $TSV_FILE | sort -k 1b,1 | uniq > RxNL_temp_1
join ${OUT_PREFIX}_NLR.list RxNL_temp_1 > ${OUT_PREFIX}_RxNL.list
gawk '{split($1, a, "."); print a[1]}' ${OUT_PREFIX}_RxNL.list | sort -k 1b,1 | uniq > RxNL_temp_2
join RxNL_temp_2 GFF3_temp | sort -k2,2 -n > RxNL_temp_3
# Reformat rows to produce new gff3 file
gawk 'BEGIN {OFS="\t"} {split($3, a, "[:\\-+]"); print  a[1], $4, $5, (a[2]+$6), (a[2]+$7), $8, $9, $10, $11}' RxNL_temp_3 > ${OUT_PREFIX}_RxNL.gff3
seqtk subseq ${PROTEIN_FILE} ${OUT_PREFIX}_RxNL.list >${OUT_PREFIX}_RxNL.fasta
gff3_to_fasta -g ${OUT_PREFIX}_RxNL.gff3 -f ${GENOME_FILE} -st cds -noQC -d simple -o ${OUT_PREFIX}_RxNL
for files in RxNL_temp_* ; do rm ${files} ; done

#
# Identify TNB
#
gawk 'BEGIN {FS="\t"} $5=="PF01582" {print $1}' $TSV_FILE | sort -k 1b,1 | uniq > TNB_temp_1
join ${OUT_PREFIX}_NBARC.list TNB_temp_1 > ${OUT_PREFIX}_TNB.list
gawk '{split($1, a, "."); print a[1]}' ${OUT_PREFIX}_TNB.list | sort -k 1b,1 | uniq > TNB_temp_2
join TNB_temp_2 GFF3_temp | sort -k2,2 -n > TNB_temp_3
# Reformat rows to produce new gff3 file
gawk 'BEGIN {OFS="\t"} {split($3, a, "[:\\-+]"); print  a[1], $4, $5, (a[2]+$6), (a[2]+$7), $8, $9, $10, $11}' TNB_temp_3 > ${OUT_PREFIX}_TNB.gff3
seqtk subseq ${PROTEIN_FILE} ${OUT_PREFIX}_TNB.list >${OUT_PREFIX}_TNB.fasta
gff3_to_fasta -g ${OUT_PREFIX}_TNB.gff3 -f ${GENOME_FILE} -st cds -noQC -d simple -o ${OUT_PREFIX}_TNB
for files in TNB_temp_* ; do rm ${files} ; done


#
# Identify TNL
#
gawk 'BEGIN {FS="\t"} $5=="PF01582" {print $1}' $TSV_FILE | sort -k 1b,1 | uniq > TNL_temp_1
join ${OUT_PREFIX}_NLR.list TNL_temp_1 > ${OUT_PREFIX}_TNL.list
gawk '{split($1, a, "."); print a[1]}' ${OUT_PREFIX}_TNL.list | sort -k 1b,1 | uniq > TNL_temp_2
join TNL_temp_2 GFF3_temp | sort -k2,2 -n > TNL_temp_3
# Reformat rows to produce new gff3 file
gawk 'BEGIN {OFS="\t"} {split($3, a, "[:\\-+]"); print  a[1], $4, $5, (a[2]+$6), (a[2]+$7), $8, $9, $10, $11}' TNL_temp_3 > ${OUT_PREFIX}_TNL.gff3
seqtk subseq ${PROTEIN_FILE} ${OUT_PREFIX}_TNL.list >${OUT_PREFIX}_TNL.fasta
gff3_to_fasta -g ${OUT_PREFIX}_TNL.gff3 -f ${GENOME_FILE} -st cds -noQC -d simple -o ${OUT_PREFIX}_TNL
for files in TNL_temp_* ; do rm ${files} ; done


#
# Identify CNB
#
gawk 'BEGIN {FS="\t"} $5=="Coil" {print $1}' $TSV_FILE | sort -k 1b,1 | uniq > CNB_temp_1
join ${OUT_PREFIX}_NBARC.list CNB_temp_1 > CNB_temp_2

# Do not want genes that have any of these
gawk 'BEGIN {FS="\t"} ($5=="PF18052" || $5=="PF05659" || $5=="PF01419" || $5=="PF01582") {print $1}' $TSV_FILE | sort -k 1b,1 | uniq > CNB_temp_exclude
grep -Fxv -f CNB_temp_exclude CNB_temp_2 > ${OUT_PREFIX}_CNB.list

gawk '{split($1, a, "."); print a[1]}' ${OUT_PREFIX}_CNB.list | sort -k 1b,1 | uniq > CNB_temp_3
join CNB_temp_3 GFF3_temp | sort -k2,2 -n > CNB_temp_4
# Reformat rows to produce new gff3 file
gawk 'BEGIN {OFS="\t"} {split($3, a, "[:\\-+]"); print  a[1], $4, $5, (a[2]+$6), (a[2]+$7), $8, $9, $10, $11}' CNB_temp_4 > ${OUT_PREFIX}_CNB.gff3
seqtk subseq ${PROTEIN_FILE} ${OUT_PREFIX}_CNB.list >${OUT_PREFIX}_CNB.fasta
gff3_to_fasta -g ${OUT_PREFIX}_CNB.gff3 -f ${GENOME_FILE} -st cds -noQC -d simple -o ${OUT_PREFIX}_CNB
for files in CNB_temp_* ; do rm ${files} ; done


#
# Identify CNL
#

# Working on from CNB
join ${OUT_PREFIX}_CNB.list ${OUT_PREFIX}_NLR.list > ${OUT_PREFIX}_CNL.list 
gawk '{split($1, a, "."); print a[1]}' ${OUT_PREFIX}_CNL.list | sort -k 1b,1 | uniq > CNL_temp_2
join CNL_temp_2 GFF3_temp | sort -k2,2 -n > CNL_temp_3
# Reformat rows to produce new gff3 file
gawk 'BEGIN {OFS="\t"} {split($3, a, "[:\\-+]"); print  a[1], $4, $5, (a[2]+$6), (a[2]+$7), $8, $9, $10, $11}' CNL_temp_3 > ${OUT_PREFIX}_CNL.gff3
seqtk subseq ${PROTEIN_FILE} ${OUT_PREFIX}_CNL.list >${OUT_PREFIX}_CNL.fasta
gff3_to_fasta -g ${OUT_PREFIX}_CNL.gff3 -f ${GENOME_FILE} -st cds -noQC -d simple -o ${OUT_PREFIX}_CNL
for files in CNL_temp_* ; do rm ${files} ; done


# Identify JNB
#
gawk 'BEGIN {FS="\t"} $5=="PF01419" {print $1}' $TSV_FILE | sort -k 1b,1 | uniq > JNB_temp_1
join ${OUT_PREFIX}_NBARC.list JNB_temp_1 > ${OUT_PREFIX}_JNB.list
gawk '{split($1, a, "."); print a[1]}' ${OUT_PREFIX}_JNB.list | sort -k 1b,1 | uniq > JNB_temp_2
join JNB_temp_2 GFF3_temp | sort -k2,2 -n > JNB_temp_3
# Reformat rows to produce new gff3 file
gawk 'BEGIN {OFS="\t"} {split($3, a, "[:\\-+]"); print  a[1], $4, $5, (a[2]+$6), (a[2]+$7), $8, $9, $10, $11}' JNB_temp_3 > ${OUT_PREFIX}_JNB.gff3
seqtk subseq ${PROTEIN_FILE} ${OUT_PREFIX}_JNB.list >${OUT_PREFIX}_JNB.fasta
gff3_to_fasta -g ${OUT_PREFIX}_JNB.gff3 -f ${GENOME_FILE} -st cds -noQC -d simple -o ${OUT_PREFIX}_JNB
for files in JNB_temp_* ; do rm ${files} ; done

#
# Identify JNL
#
gawk 'BEGIN {FS="\t"} $5=="PF01419" {print $1}' $TSV_FILE | sort -k 1b,1 | uniq > JNL_temp_1
join ${OUT_PREFIX}_NLR.list JNL_temp_1 > ${OUT_PREFIX}_JNL.list
gawk '{split($1, a, "."); print a[1]}' ${OUT_PREFIX}_JNL.list | sort -k 1b,1 | uniq > JNL_temp_2
join JNL_temp_2 GFF3_temp | sort -k2,2 -n > JNL_temp_3
# Reformat rows to produce new gff3 file
gawk 'BEGIN {OFS="\t"} {split($3, a, "[:\\-+]"); print a[1], $4, $5, (a[2]+$6), (a[2]+$7), $8, $9, $10, $11}' JNL_temp_3 > ${OUT_PREFIX}_JNL.gff3
seqtk subseq ${PROTEIN_FILE} ${OUT_PREFIX}_JNL.list >${OUT_PREFIX}_JNL.fasta
gff3_to_fasta -g ${OUT_PREFIX}_JNL.gff3 -f ${GENOME_FILE} -st cds -noQC -d simple -o ${OUT_PREFIX}_JNL
for files in JNL_temp_* ; do rm ${files} ; done

rm GFF3_temp
