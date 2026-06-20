#!/usr/bin/env bash

#Copy from https://github.com/peritob/Myrtaceae_NLR_workflow/blob/master/get_integrated_domains.sh

set -euo pipefail

usage() {
    echo "Usage: $0 --prefix <prefix> --output_dir <output directory> --type <NLR type>"
    exit 1
}

PREFIX=""
OUTPUT_DIR="."
NLR_TYPE=""

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
        --type)
            if [[ -n "${2:-}" && ! "${2:-}" =~ ^- ]]; then
                NLR_TYPE="$2"
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
    esac
done

OUTPUT_DIR="${OUTPUT_DIR%/}/"

if [[ -z "$PREFIX" || -z "$NLR_TYPE" ]]; then
    echo "Error: Missing required input." >&2
    usage
fi

#First extract all Pfam codes from the .tsv - thereby excluding all coils and other database outputs.
grep 'Pfam' ${OUTPUT_DIR}${PREFIX}.tsv > ${OUTPUT_DIR}${PREFIX}_pfam.tsv
#Then extract all the genes based on the lists from NLR_clasification. For example,
grep -F -w -f ${OUTPUT_DIR}${PREFIX}_${NLR_TYPE}.list ${OUTPUT_DIR}${PREFIX}_pfam.tsv > ${OUTPUT_DIR}${PREFIX}_${NLR_TYPE}_Pfams.tsv
#Then extract all domains that are not present in NLR classifications using the parameter -v with fgrep. Make a list of NLR Pfams and call the txt file "NLR_Pfams.list" (PF00931, PF18052, PF08263, PF12799, PF13306, PF13855, PF13516, PF00560, PF07725, PF05659, PF01582, PF17862, PF23598, PF23247) associated with NLR genes.
grep -F -wv -f NLR_Pfams.list ${OUTPUT_DIR}${PREFIX}_${NLR_TYPE}_Pfams.tsv > tmp/${PREFIX}_${NLR_TYPE}_ID_Pfams.tsv
#Obtain meaningful column data 
awk 'BEGIN {FS="\t"; OFS="\t"}; {print $1, $5, $6, $7, $8}' tmp/${PREFIX}_${NLR_TYPE}_ID_Pfams.tsv > ${OUTPUT_DIR}${PREFIX}_${NLR_TYPE}_ID_Pfams.tsv

