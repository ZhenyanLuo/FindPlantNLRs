#!/usr/bin/env bash

#From a bash script from https://github.com/peritob/Myrtaceae_NLR_workflow

#Use for NLR classification step in the snakemake pipeline, need seqtk and gff3tools, can be installed with conda inastall -c bioconda seqtk -y and pip install gff3tool
#gff3tools is from https://github.com/NAL-i5K/GFF3toolkit
#This script can be used individually by running 'bash run_classification.sh {interproscan.tsv} {augustus.gff3} {corresponding fasta file}'

set -euo pipefail

# Help/Usage function
usage() {
    echo "Usage: $0 --prefix <prefix> --output_dir <output directory> <interproscan.tsv> <augustus.gff3> <protein.fasta> <genome.fasta>"
    exit 1
}

# Initialise empty tracking variables
PREFIX=""
OUTPUT_DIR="."
TSV_FILE=""
GFF3_FILE=""
PROTEIN_FILE=""
GENOME_FILE=""
CODING_SEQ_FILE=""

# Parse command line flags/switches
while [[ $# -gt 0 ]]; do
    case "$1" in
        --prefix)
            if [[ -n "${2:-}" && ! "${2:-}" =~ ^- ]]; then
                PREFIX="$2"
                shift 2
            else
                echo "Error: Argument for $1 is missing." >&2
                usage
            fi
            ;;
        --output_dir)
            if [[ -n "${2:-}" && ! "${2:-}" =~ ^- ]]; then
                OUTPUT_DIR="$2"
                shift 2
            else
                echo "Error: Argument for $1 is missing." >&2
                usage
            fi
            ;;
        -h|--help)
            usage
            ;;
        -*)
            echo "Error: Unknown option $1" >&2
            usage
            ;;
        *)
            # Collect positional parameters (the files)
            if [[ -z "$TSV_FILE" ]]; then
                TSV_FILE="$1"
            elif [[ -z "$GFF3_FILE" ]]; then
                GFF3_FILE="$1"
            elif [[ -z "$PROTEIN_FILE" ]]; then
                PROTEIN_FILE="$1"
            elif [[ -z "$GENOME_FILE" ]]; then
                GENOME_FILE="$1"
            elif [[ -z "$CODING_SEQ_FILE" ]]; then
                CODING_SEQ_FILE="$1"
            else
                echo "Error: Too many positional arguments provided." >&2
                usage
            fi
            shift
            ;;
    esac
done

# Ensure all 4 mandatory input files were passed
if [[ -z "$TSV_FILE" || -z "$GFF3_FILE" || -z "$PROTEIN_FILE" || -z "$GENOME_FILE" ]]; then
    echo "Error: Missing required input files." >&2
    usage
fi

# Verify files exist before proceeding
for f in "$TSV_FILE" "$GFF3_FILE" "$PROTEIN_FILE" "$GENOME_FILE"; do
    if [[ ! -f "$f" ]]; then
        echo "Error: File '$f' does not exist." >&2
        exit 1
    fi
done

# Fallback automated naming logic if a prefix switch wasn't specified
if [[ -z "$PREFIX" ]]; then
    PREFIX="${TSV_FILE%.tsv}"
    PREFIX="${PREFIX%.fasta}"
    PREFIX="${PREFIX%_augustus_aa}"
fi

OUTPUT_DIR="${OUTPUT_DIR%/}/"

echo "tsv file:        $TSV_FILE"
echo "gff3 file:       $GFF3_FILE"
echo "protein file:    $PROTEIN_FILE"
echo "genome file:     $GENOME_FILE"
echo "coding seq file: $CODING_SEQ_FILE"
echo "prefix:          $PREFIX"

# Build a temporary working file which is the gff3 file but with two columns inserted on the left.
# The first column contains the gene ID without the extension (extracted from column 9).
# The second column contains the row number so we can return to the original gff3 order.
awk 'BEGIN {FS="\t"; OFS="\t"} {split($9, a, "[=\\.;]"); print a[2], NR, $0}' $GFF3_FILE | sort -k 1b,1 > GFF3_temp


#
# Identify NBARC
#
awk 'BEGIN {FS="\t"} $5=="PF00931" {print $1}' $TSV_FILE | sort -k 1b,1 | uniq > ${OUTPUT_DIR}${PREFIX}_NBARC.list
awk '{split($1, a, "."); print a[1]}' ${OUTPUT_DIR}${PREFIX}_NBARC.list | sort -k 1b,1 | uniq > NBARC_temp_1
join NBARC_temp_1 GFF3_temp | sort -k2,2 -n > NBARC_temp_2
# Reformat rows to produce new gff3 file.  This outputs the annotated gene co-ordinates and orientation based on the original genome using braker outputs.
awk 'BEGIN {OFS="\t"} {split($3, a, "[:\\-+]"); print  a[1], $4, $5, (a[2]+$6), (a[2]+$7), $8, $9, $10, $11}' NBARC_temp_2 > ${OUTPUT_DIR}${PREFIX}_NBARC.gff3
seqtk subseq ${PROTEIN_FILE} ${OUTPUT_DIR}${PREFIX}_NBARC.list > ${OUTPUT_DIR}${PREFIX}_NBARC_AA.fasta
seqtk subseq ${CODING_SEQ_FILE} ${OUTPUT_DIR}${PREFIX}_NBARC.list > ${OUTPUT_DIR}${PREFIX}_NBARC_CDS.fasta
for files in NBARC_temp_*; do rm ${files}; done


#
# Identify NLR
#
awk 'BEGIN {FS="\t"} ($5=="G3DSA:3.80.10.10" || $5=="PF08263" || $5=="PF13516" || $5=="PF12799" || $5=="PF13306" || $5=="PF13855") {print $1}' $TSV_FILE | sort -k 1b,1 | uniq > NLR_temp_1
join ${OUTPUT_DIR}${PREFIX}_NBARC.list NLR_temp_1 > ${OUTPUT_DIR}${PREFIX}_NLR.list
awk '{split($1, a, "."); print a[1]}' ${OUTPUT_DIR}${PREFIX}_NLR.list | sort -k 1b,1 | uniq > NLR_temp_2
join NLR_temp_2 GFF3_temp | sort -k2,2 -n > NLR_temp_3
# Reformat rows to produce new gff3 file
awk 'BEGIN {OFS="\t"} {split($3, a, "[:\\-+]"); print  a[1], $4, $5, (a[2]+$6), (a[2]+$7), $8, $9, $10, $11}' NLR_temp_3 > ${OUTPUT_DIR}${PREFIX}_NLR.gff3
seqtk subseq ${PROTEIN_FILE} ${OUTPUT_DIR}${PREFIX}_NLR.list > ${OUTPUT_DIR}${PREFIX}_NLR_AA.fasta
seqtk subseq ${CODING_SEQ_FILE} ${OUTPUT_DIR}${PREFIX}_NLR.list > ${OUTPUT_DIR}${PREFIX}_NLR_CDS.fasta
for files in NLR_temp_*; do rm ${files}; done


#
# Identify RNB
#
awk 'BEGIN {FS="\t"} $5=="PF05659" {print $1}' $TSV_FILE | sort -k 1b,1 | uniq > RNB_temp_1
join ${OUTPUT_DIR}${PREFIX}_NBARC.list RNB_temp_1 > ${OUTPUT_DIR}${PREFIX}_RNB.list
awk '{split($1, a, "."); print a[1]}' ${OUTPUT_DIR}${PREFIX}_RNB.list | sort -k 1b,1 | uniq > RNB_temp_2
join RNB_temp_2 GFF3_temp | sort -k2,2 -n > RNB_temp_3
# Reformat rows to produce new gff3 file
awk 'BEGIN {OFS="\t"} {split($3, a, "[:\\-+]"); print  a[1], $4, $5, (a[2]+$6), (a[2]+$7), $8, $9, $10, $11}' RNB_temp_3 > ${OUTPUT_DIR}${PREFIX}_RNB.gff3
seqtk subseq ${PROTEIN_FILE} ${OUTPUT_DIR}${PREFIX}_RNB.list > ${OUTPUT_DIR}${PREFIX}_RNB_AA.fasta
seqtk subseq ${CODING_SEQ_FILE} ${OUTPUT_DIR}${PREFIX}_RNB.list > ${OUTPUT_DIR}${PREFIX}_RNB_CDS.fasta
for files in RNB_temp_* ; do rm ${files} ; done

#
# Identify RNL
#
awk 'BEGIN {FS="\t"} $5=="PF05659" {print $1}' $TSV_FILE | sort -k 1b,1 | uniq > RNL_temp_1
join ${OUTPUT_DIR}${PREFIX}_NLR.list RNL_temp_1 > ${OUTPUT_DIR}${PREFIX}_RNL.list
awk '{split($1, a, "."); print a[1]}' ${OUTPUT_DIR}${PREFIX}_RNL.list | sort -k 1b,1 | uniq > RNL_temp_2
join RNL_temp_2 GFF3_temp | sort -k2,2 -n > RNL_temp_3
# Reformat rows to produce new gff3 file
awk 'BEGIN {OFS="\t"} {split($3, a, "[:\\-+]"); print  a[1], $4, $5, (a[2]+$6), (a[2]+$7), $8, $9, $10, $11}' RNL_temp_3 > ${OUTPUT_DIR}${PREFIX}_RNL.gff3
seqtk subseq ${PROTEIN_FILE} ${OUTPUT_DIR}${PREFIX}_RNL.list > ${OUTPUT_DIR}${PREFIX}_RNL_AA.fasta
seqtk subseq ${CODING_SEQ_FILE} ${OUTPUT_DIR}${PREFIX}_RNL.list > ${OUTPUT_DIR}${PREFIX}_RNL_CDS.fasta
for files in RNL_temp_* ; do rm ${files} ; done


#
# Identify RxNB
#
awk 'BEGIN {FS="\t"} $5=="PF18052" {print $1}' $TSV_FILE | sort -k 1b,1 | uniq > RxNB_temp_1
join ${OUTPUT_DIR}${PREFIX}_NBARC.list RxNB_temp_1 > ${OUTPUT_DIR}${PREFIX}_RxNB.list
awk '{split($1, a, "."); print a[1]}' ${OUTPUT_DIR}${PREFIX}_RxNB.list | sort -k 1b,1 | uniq > RxNB_temp_2
join RxNB_temp_2 GFF3_temp | sort -k2,2 -n > RxNB_temp_3
# Reformat rows to produce new gff3 file
awk 'BEGIN {OFS="\t"} {split($3, a, "[:\\-+]"); print  a[1], $4, $5, (a[2]+$6), (a[2]+$7), $8, $9, $10, $11}' RxNB_temp_3 > ${OUTPUT_DIR}${PREFIX}_RxNB.gff3
seqtk subseq ${PROTEIN_FILE} ${OUTPUT_DIR}${PREFIX}_RxNB.list > ${OUTPUT_DIR}${PREFIX}_RxNB_AA.fasta
seqtk subseq ${CODING_SEQ_FILE} ${OUTPUT_DIR}${PREFIX}_RxNB.list > ${OUTPUT_DIR}${PREFIX}_RxNB_CDS.fasta
for files in RxNB_temp_* ; do rm ${files} ; done


#
# Identify RxNL
#
awk 'BEGIN {FS="\t"} $5=="PF18052" {print $1}' $TSV_FILE | sort -k 1b,1 | uniq > RxNL_temp_1
join ${OUTPUT_DIR}${PREFIX}_NLR.list RxNL_temp_1 > ${OUTPUT_DIR}${PREFIX}_RxNL.list
awk '{split($1, a, "."); print a[1]}' ${OUTPUT_DIR}${PREFIX}_RxNL.list | sort -k 1b,1 | uniq > RxNL_temp_2
join RxNL_temp_2 GFF3_temp | sort -k2,2 -n > RxNL_temp_3
# Reformat rows to produce new gff3 file
awk 'BEGIN {OFS="\t"} {split($3, a, "[:\\-+]"); print  a[1], $4, $5, (a[2]+$6), (a[2]+$7), $8, $9, $10, $11}' RxNL_temp_3 > ${OUTPUT_DIR}${PREFIX}_RxNL.gff3
seqtk subseq ${PROTEIN_FILE} ${OUTPUT_DIR}${PREFIX}_RxNL.list > ${OUTPUT_DIR}${PREFIX}_RxNL_AA.fasta
seqtk subseq ${CODING_SEQ_FILE} ${OUTPUT_DIR}${PREFIX}_RxNL.list > ${OUTPUT_DIR}${PREFIX}_RxNL_CDS.fasta
for files in RxNL_temp_* ; do rm ${files} ; done

#
# Identify TNB
#
awk 'BEGIN {FS="\t"} $5=="PF01582" {print $1}' $TSV_FILE | sort -k 1b,1 | uniq > TNB_temp_1
join ${OUTPUT_DIR}${PREFIX}_NBARC.list TNB_temp_1 > ${OUTPUT_DIR}${PREFIX}_TNB.list
awk '{split($1, a, "."); print a[1]}' ${OUTPUT_DIR}${PREFIX}_TNB.list | sort -k 1b,1 | uniq > TNB_temp_2
join TNB_temp_2 GFF3_temp | sort -k2,2 -n > TNB_temp_3
# Reformat rows to produce new gff3 file
awk 'BEGIN {OFS="\t"} {split($3, a, "[:\\-+]"); print  a[1], $4, $5, (a[2]+$6), (a[2]+$7), $8, $9, $10, $11}' TNB_temp_3 > ${OUTPUT_DIR}${PREFIX}_TNB.gff3
seqtk subseq ${PROTEIN_FILE} ${OUTPUT_DIR}${PREFIX}_TNB.list > ${OUTPUT_DIR}${PREFIX}_TNB_AA.fasta
seqtk subseq ${CODING_SEQ_FILE} ${OUTPUT_DIR}${PREFIX}_TNB.list > ${OUTPUT_DIR}${PREFIX}_TNB_CDS.fasta
for files in TNB_temp_* ; do rm ${files} ; done


#
# Identify TNL
#
awk 'BEGIN {FS="\t"} $5=="PF01582" {print $1}' $TSV_FILE | sort -k 1b,1 | uniq > TNL_temp_1
join ${OUTPUT_DIR}${PREFIX}_NLR.list TNL_temp_1 > ${OUTPUT_DIR}${PREFIX}_TNL.list
awk '{split($1, a, "."); print a[1]}' ${OUTPUT_DIR}${PREFIX}_TNL.list | sort -k 1b,1 | uniq > TNL_temp_2
join TNL_temp_2 GFF3_temp | sort -k2,2 -n > TNL_temp_3
# Reformat rows to produce new gff3 file
awk 'BEGIN {OFS="\t"} {split($3, a, "[:\\-+]"); print  a[1], $4, $5, (a[2]+$6), (a[2]+$7), $8, $9, $10, $11}' TNL_temp_3 > ${OUTPUT_DIR}${PREFIX}_TNL.gff3
seqtk subseq ${PROTEIN_FILE} ${OUTPUT_DIR}${PREFIX}_TNL.list > ${OUTPUT_DIR}${PREFIX}_TNL_AA.fasta
seqtk subseq ${CODING_SEQ_FILE} ${OUTPUT_DIR}${PREFIX}_TNL.list > ${OUTPUT_DIR}${PREFIX}_TNL_CDS.fasta
for files in TNL_temp_* ; do rm ${files} ; done


#
# Identify CNB
#
awk 'BEGIN {FS="\t"} $5=="Coil" {print $1}' $TSV_FILE | sort -k 1b,1 | uniq > CNB_temp_1
join ${OUTPUT_DIR}${PREFIX}_NBARC.list CNB_temp_1 > CNB_temp_2

# Do not want genes that have any of these
awk 'BEGIN {FS="\t"} ($5=="PF18052" || $5=="PF05659" || $5=="PF01419" || $5=="PF01582") {print $1}' $TSV_FILE | sort -k 1b,1 | uniq > CNB_temp_exclude
grep -Fxv -f CNB_temp_exclude CNB_temp_2 > ${OUTPUT_DIR}${PREFIX}_CNB.list

awk '{split($1, a, "."); print a[1]}' ${OUTPUT_DIR}${PREFIX}_CNB.list | sort -k 1b,1 | uniq > CNB_temp_3
join CNB_temp_3 GFF3_temp | sort -k2,2 -n > CNB_temp_4
# Reformat rows to produce new gff3 file
awk 'BEGIN {OFS="\t"} {split($3, a, "[:\\-+]"); print  a[1], $4, $5, (a[2]+$6), (a[2]+$7), $8, $9, $10, $11}' CNB_temp_4 > ${OUTPUT_DIR}${PREFIX}_CNB.gff3
seqtk subseq ${PROTEIN_FILE} ${OUTPUT_DIR}${PREFIX}_CNB.list > ${OUTPUT_DIR}${PREFIX}_CNB_AA.fasta
seqtk subseq ${CODING_SEQ_FILE} ${OUTPUT_DIR}${PREFIX}_CNB.list > ${OUTPUT_DIR}${PREFIX}_CNB_CDS.fasta
for files in CNB_temp_* ; do rm ${files} ; done


#
# Identify CNL
#

# Working on from CNB
join ${OUTPUT_DIR}${PREFIX}_CNB.list ${OUTPUT_DIR}${PREFIX}_NLR.list > ${OUTPUT_DIR}${PREFIX}_CNL.list 
awk '{split($1, a, "."); print a[1]}' ${OUTPUT_DIR}${PREFIX}_CNL.list | sort -k 1b,1 | uniq > CNL_temp_2
join CNL_temp_2 GFF3_temp | sort -k2,2 -n > CNL_temp_3
# Reformat rows to produce new gff3 file
awk 'BEGIN {OFS="\t"} {split($3, a, "[:\\-+]"); print  a[1], $4, $5, (a[2]+$6), (a[2]+$7), $8, $9, $10, $11}' CNL_temp_3 > ${OUTPUT_DIR}${PREFIX}_CNL.gff3
seqtk subseq ${PROTEIN_FILE} ${OUTPUT_DIR}${PREFIX}_CNL.list > ${OUTPUT_DIR}${PREFIX}_CNL_AA.fasta
seqtk subseq ${CODING_SEQ_FILE} ${OUTPUT_DIR}${PREFIX}_CNL.list > ${OUTPUT_DIR}${PREFIX}_CNL_CDS.fasta
for files in CNL_temp_* ; do rm ${files} ; done


# Identify JNB
#
awk 'BEGIN {FS="\t"} $5=="PF01419" {print $1}' $TSV_FILE | sort -k 1b,1 | uniq > JNB_temp_1
join ${OUTPUT_DIR}${PREFIX}_NBARC.list JNB_temp_1 > ${OUTPUT_DIR}${PREFIX}_JNB.list
awk '{split($1, a, "."); print a[1]}' ${OUTPUT_DIR}${PREFIX}_JNB.list | sort -k 1b,1 | uniq > JNB_temp_2
join JNB_temp_2 GFF3_temp | sort -k2,2 -n > JNB_temp_3
# Reformat rows to produce new gff3 file
awk 'BEGIN {OFS="\t"} {split($3, a, "[:\\-+]"); print  a[1], $4, $5, (a[2]+$6), (a[2]+$7), $8, $9, $10, $11}' JNB_temp_3 > ${OUTPUT_DIR}${PREFIX}_JNB.gff3
seqtk subseq ${PROTEIN_FILE} ${OUTPUT_DIR}${PREFIX}_JNB.list > ${OUTPUT_DIR}${PREFIX}_JNB_AA.fasta
seqtk subseq ${CODING_SEQ_FILE} ${OUTPUT_DIR}${PREFIX}_JNB.list > ${OUTPUT_DIR}${PREFIX}_JNB_CDS.fasta
for files in JNB_temp_* ; do rm ${files} ; done

#
# Identify JNL
#
awk 'BEGIN {FS="\t"} $5=="PF01419" {print $1}' $TSV_FILE | sort -k 1b,1 | uniq > JNL_temp_1
join ${OUTPUT_DIR}${PREFIX}_NLR.list JNL_temp_1 > ${OUTPUT_DIR}${PREFIX}_JNL.list
awk '{split($1, a, "."); print a[1]}' ${OUTPUT_DIR}${PREFIX}_JNL.list | sort -k 1b,1 | uniq > JNL_temp_2
join JNL_temp_2 GFF3_temp | sort -k2,2 -n > JNL_temp_3
# Reformat rows to produce new gff3 file
awk 'BEGIN {OFS="\t"} {split($3, a, "[:\\-+]"); print a[1], $4, $5, (a[2]+$6), (a[2]+$7), $8, $9, $10, $11}' JNL_temp_3 > ${OUTPUT_DIR}${PREFIX}_JNL.gff3
seqtk subseq ${PROTEIN_FILE} ${OUTPUT_DIR}${PREFIX}_JNL.list > ${OUTPUT_DIR}${PREFIX}_JNL_AA.fasta
seqtk subseq ${CODING_SEQ_FILE} ${OUTPUT_DIR}${PREFIX}_JNL.list > ${OUTPUT_DIR}${PREFIX}_JNL_CDS.fasta
for files in JNL_temp_* ; do rm ${files} ; done

rm GFF3_temp
