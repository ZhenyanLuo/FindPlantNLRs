import os


configfile: "FindPlantNLRs.config"


SAMPLES, = glob_wildcards("./genome/{sample}.fasta")

if not os.path.exists("tmp"):
    os.mkdir("tmp")
if not os.path.exists("result"):
    os.mkdir("result")

rule all:
    input:
        expand("result/{sample}_NLRclassification.done", sample=SAMPLES),
        expand("result/{sample}_pfam.tsv", sample=SAMPLES),
        expand("result/{sample}_NLR_Pfams.tsv", sample=SAMPLES),
        expand("result/{sample}_NLR_ID_Pfams.tsv", sample=SAMPLES)


# Chopping the genome sequence into overlapping subsequences
rule chop_sequence:
    input:
        "genome/{sample}.fasta"
    output:
        "tmp/{sample}.choppedseq.fasta"
    log:
        "logs/chop_sequence/{sample}.log"
    shell:
        r"""
        java -jar {config[NLR-Annotator]}/ChopSequence.jar -i {input} -o {output} -l 20000 -p 5000 > {log} 2>&1
        """


# Searching the chopped subsequences for pre-determined -associated motifs
rule search_motifs:
    input:
        "tmp/{sample}.choppedseq.fasta"
    output:
        "tmp/{sample}.parser.xml"
    log:
        "logs/search_motifs/{sample}.log"
    shell:
        r"""
        java -jar {config[NLR-Annotator]}/NLR-Parser3.jar -t 10 -y {config[meme]} -x {config[meme_xml]} -i {input} -c {output} > {log} 2>&1
        """


# Generate the GFF format of  loci for the searched motifs
rule run_NLR_annotator:
    input:
        "tmp/{sample}.parser.xml"
    output:
        "tmp/{sample}.parser.gff"
    log:
        "logs/run_NLR_annotator/{sample}.log"
    shell:
        r"""
        java -jar {config[NLR-Annotator]}/NLR-Annotator.jar -i {input} -g {output} > {log} 2>&1
        """


# Indexing genome sequence
rule index_genome:
    input:
        "genome/{sample}.fasta"
    output:
        "genome/{sample}.fasta.fai"
    log:
        "logs/index_genome/{sample}.log"
    shell:
        "samtools faidx {input} > {log} 2>&1"


# Creat a file which have contig length information.
rule get_contig_length:
    input:
        "genome/{sample}.fasta.fai"
    output:
        "genome/{sample}.genomefile"
    log:
        "logs/get_contig_length/{sample}.log"
    shell:
        "cut -d $'\t' -f1,2 {input} > {output} 2> {log}"


# Convert format of parser.gff in order to incorporate with hmm output later
rule convert_parser_format:
    input:
        "tmp/{sample}.parser.gff"
    output:
        "tmp/{sample}.parser.bed"
    log:
        "logs/convert_parser_format/{sample}.log"
    shell:
        r"""
        awk -v OFS='\t' '{{if ($7 == "+") {{print $1, $4, $5, $1, "forward", $7}} else if ($7 == "-") print $1, $4, $5, $1, "reverse", $7}}' {input} > {output} 2> {log}
        """


# Convert to 20kbflanking bed file with bedtools
rule get_parser_20kbflanking:
    input:
        bed="tmp/{sample}.parser.bed",
        genomefile="genome/{sample}.genomefile"
    output:
        "tmp/{sample}_parser.20kbflanking.bed"
    log:
        "logs/get_parser_20kbflanking/{sample}.log"
    shell:
        r"""
        bedtools slop -b 20000 -s -i {input.bed} -g {input.genomefile} | bedtools sort -i - | bedtools merge -s -d 1 -c 1,5,6 -o distinct,distinct,distinct > {output} 2> {log}
        """


# ----------------------------------------------Use blast to identify genes which cannot be detected by  annotator pipeline------------------------------------------
# -----------------------------------------------------------------part 2--------------------------------------------------------------------------------------------
# Make a genome database for detecting nucleotide or protein query sequence
# Build blastdb and dectect whether there are genes which cannot be captured by using -parser by using tblastn# remember to form a folder which include blastprotein
rule run_tblastn:
    input:
        genome="genome/{sample}.fasta",
    output:
        "tmp/{sample}.tblastnout.outfmt6"
    params:
        "tmp/{sample}.database"
    log:
        "logs/run_tblastn/{sample}.log"
    shell:
        r"""
        makeblastdb -in {input.genome} -dbtype nucl -parse_seqids -out {params} > {log} 2>&1
        tblastn -query {config[ref]} -db {params} -num_threads {config[blast_threads]} -evalue 0.001 -outfmt 6 -out {output} >> {log} 2>&1
        """


# Convert tblastn file into bed, and add two coloums for strand information
rule convert_tblastn_format:
    input:
        "tmp/{sample}.tblastnout.outfmt6",
    output:
        "tmp/{sample}.tblastnout.bed"
    log:
        "logs/convert_tblastn_format/{sample}.log"
    shell:
        r"""
        awk -v OFS='\t' '{{if ($10 - $9 > 0) {{print $2, $9, $10, $1, "forward", "+"}} else if ($10 - $9 < 0) print $2, $10, $9, $1, "reverse", "-"}}' {input} > {output} 2> {log}
        """


##Generate 20kb flanking bed file for blast
rule get_blast_20kb:
    input:
        bed="tmp/{sample}.tblastnout.bed",
        genomefile="genome/{sample}.genomefile"
    output:
        "tmp/{sample}.blast.20kbflanking.bed"
    log:
        "logs/get_blast_20kb/{sample}.log"
    shell:
        r"""
        bedtools slop -b 20000 -s -i {input.bed} -g {input.genomefile} | bedtools sort -i - | bedtools merge -s -d 1 -c 1,5,6 -o distinct,distinct,distinct > {output} 2> {log}
        """


# Adapted from Peri Tobias' s scripts------------------------------------------------------------------------------------------------
# Use nhmmer to search for conserved nucleotide binding domain shared by Apaf-1, Resistance proteins and CED4 from coiled-coil  and TIR  sequences
# -----------------------------------------------------------------part 3--------------------------------------------------------------------------------------------
# Peri has already prepared hmm profiles which are named as EG_nonTIRhmm and EG_TIRhmm respectively.
rule find_TIR:
    input:
        nonTIR="ref_db/EG_nonTIRhmm",
        TIR="ref_db/EG_TIRhmm",
        genome="genome/{sample}.fasta"
    output:
        TIR_out="tmp/{sample}.TIRout",
        nonTIR_out="tmp/{sample}.nonTIRout"
    log:
        "logs/find_TIR/{sample}.log"
    shell:
        r"""
        nhmmer --dna {input.nonTIR} {input.genome} > {output.nonTIR_out} 2> {log}
        nhmmer --dna {input.TIR} {input.genome} > {output.TIR_out} 2>> {log}
        """


# convert both nhmmer output into bed file by using awk script
rule awk_convert:
    input:
        TIR_out="tmp/{sample}.TIRout",
        nonTIR_out="tmp/{sample}.nonTIRout",
    output:
        TIRoutbed="tmp/{sample}.TIRout.bed",
        nonTIRoutbed="tmp/{sample}.nonTIRout.bed",
    log:
        "logs/awk_convert/{sample}.log"
    shell:
        r"""
        awk '/Scores for complete hits/{{flag=1;next}}/------ inclusion threshold ------/{{flag=0}}flag' {input.TIR_out} | awk '(NR >2)' | awk -v OFS='\t' '{{if ($6 - $5 > 0) {{print $4, $5, $6, $4, "forward", "+"}} else if ($6 - $5 < 0) print $4, $6, $5, $4, "reverse", "-"}}' > {output.TIRoutbed} 2> {log}
        awk '/Scores for complete hits/{{flag=1;next}}/------ inclusion threshold ------/{{flag=0}}flag' {input.nonTIR_out} | awk '(NR >2)' | awk -v OFS='\t' '{{if ($6 - $5 > 0) {{print $4, $5, $6, $4, "forward", "+"}} else if ($6 - $5 < 0) print $4, $6, $5, $4, "reverse", "-"}}' >{output.nonTIRoutbed} 2>> {log}
        """


# convert both bed file into fasta file by using bedtools
rule get_fasta_sequence:
    input:
        TIR_bed="tmp/{sample}.TIRout.bed",
        nonTIR_bed="tmp/{sample}.nonTIRout.bed",
        genome="genome/{sample}.fasta"
    output:
        TIR_fasta="tmp/{sample}.TIR.fasta",
        nonTIR_fasta="tmp/{sample}.nonTIR.fasta",
    log:
        "logs/get_fasta_sequence/{sample}.log"
    shell:
        r"""
        bedtools getfasta -s -fi {input.genome} -bed {input.TIR_bed} -fo {output.TIR_fasta} > {log} 2>&1
        bedtools getfasta -s -fi {input.genome} -bed {input.nonTIR_bed} -fo {output.nonTIR_fasta} >> {log} 2>&1
        """


# Extract first 200 sequences from previously output nonTIR and TIR fasta file, change 200 into other number when neccessary
rule awk:
    input:
        nonTIR="tmp/{sample}.nonTIR.fasta",
        TIR="tmp/{sample}.TIR.fasta",
    output:
        nonTIR_200="tmp/{sample}.nonTIR_200.fasta",
        TIR_200="tmp/{sample}.TIR_200.fasta",
        NBARC400="tmp/{sample}_NBARC_400.fasta",
    log:
        "logs/awk/{sample}.log"
    shell:
        r"""
        awk '/^>/ {{n++}} n>200 {{exit}} 1' {input.nonTIR} > {output.nonTIR_200} 2> {log}
        awk '/^>/ {{n++}} n>200 {{exit}} 1' {input.TIR} > {output.TIR_200} 2>> {log}
        cat {output.nonTIR_200} {output.TIR_200} > {output.NBARC400} 2>> {log}
        """


# Remove duplicated sequences
rule run_seqkit:
    input:
        "tmp/{sample}_NBARC_400.fasta",
    output:
        "tmp/{sample}_NBARC.fasta",
    log:
        "logs/run_seqkit/{sample}.log"
    shell:
        "seqkit rmdup -D duplicates -n {input} > {output} 2> {log}"


# Generate alignment by using clustalo
rule run_clustalo:
    input:
        "tmp/{sample}_NBARC.fasta",
    output:
        "tmp/{sample}_NBARC.sto",
    log:
        "logs/run_clustalo/{sample}.log"
    shell:
        "clustalo -i {input} -o {output} --outfmt=st > {log} 2>&1"


# Build a profile hmm from an alignment
rule run_hmmbuild:
    input:
        "tmp/{sample}_NBARC.sto"
    output:
        "tmp/{sample}.hmm",
    log:
        "logs/run_hmmbuild/{sample}.log"
    shell:
        "hmmbuild -nucleic {output} {input} > {log} 2>&1"


# Use hmm profile built to search queries against genome
rule run_nhmmer:
    input:
        hmm="tmp/{sample}.hmm",
        genome="genome/{sample}.fasta"
    output:
        "tmp/{sample}_NBARCout",
    log:
        "logs/run_nhmmer/{sample}.log"
    shell:
        "nhmmer --dna {input.hmm} {input.genome} > {output} 2> {log}"


# Convert nhmmer output into bed file
# Required double-check
rule make_bed_hmmout:
    input:
        NBARC="tmp/{sample}_NBARCout",
    output:
        "tmp/{sample}_NBARC.bed",
    log:
        "logs/make_bed_hmmout/{sample}.log"
    shell:
        r"""
        cat {input.NBARC} | awk '/Scores for complete hits/{{flag=1;next}}/------ inclusion threshold ------/{{flag=0}}flag' | awk '(NR >2)' | awk -v OFS='\t' '{{if ($6 - $5 > 0) {{print $4, $5, $6, $4, "forward", "+"}} else if ($6 - $5 < 0) print $4, $6, $5, $4, "reverse", "-"}}' > {output} 2> {log}
        """


# Get 20kb upstream and downstream
rule get_NBARC_20kb:
    input:
        bed="tmp/{sample}_NBARC.bed",
        genomefile="genome/{sample}.genomefile",
    output:
        "tmp/{sample}_NBARC.20kbflanking.bed",
    log:
        "logs/get_NBARC_20kb/{sample}.log"
    shell:
        "bedtools slop -b 20000 -s -i {input.bed} -g {input.genomefile} > {output} 2> {log}"


# ---------------------------------------Now we have output from hmm, blast and _annotator, combine them into one file--------------------------------------------
# -----------------------------------------------------------------part 4--------------------------------------------------------------------------------------------

# Combine the output together
# Include blast file after getting the query file
rule combine_all_bed:
    input:
        hmm="tmp/{sample}_NBARC.20kbflanking.bed",
        annotator="tmp/{sample}_parser.20kbflanking.bed",
        blast="tmp/{sample}.blast.20kbflanking.bed",
    output:
        "tmp/{sample}_all20kbflanking.bed",
    log:
        "logs/combine_all_bed/{sample}.log"
    shell:
        r"""
        cat {input.hmm} {input.annotator} {input.blast} | awk -v OFS='\t' '{{print $1,$2,$3}}' | sort -k1,1 -k2,2n > {output} 2> {log}
        """


# Merge the bed file now
rule merge_all_20kbflanking:
    input:
        "tmp/{sample}_all20kbflanking.bed",
    output:
        "tmp/{sample}.all_20kbflanking_merged.bed",
    log:
        "logs/merge_all_20kbflanking/{sample}.log"
    shell:
        #         bedtools merge -s -d 1 -c 1,5,6 -o distinct,distinct,distinct -i {input} | awk -v OFS='\t' '{{print $1,$2,$3}}' | sort -k1,2 | uniq  > {output}
        #         bedtools merge -s -d 1 -c 1,5,6 -o distinct,distinct,distinct -i {input} | sort -k1,2 |uniq  > {output}
        "bedtools merge -i {input} > {output} 2> {log}"


# Convert bedfile into fasta
rule convert_20kbflankingbedfile_fasta:
    input:
        genome="genome/{sample}.fasta",
        bed="tmp/{sample}.all_20kbflanking_merged.bed"
    output:
        "tmp/{sample}.all_20kbflanking_merged.fasta"
    log:
        "logs/convert_20kbflankingbedfile_fasta/{sample}.log"
    shell:
        "bedtools getfasta -fi {input.genome} -bed {input.bed} -fo {output} > {log} 2>&1"


rule prep_braker:
    input:
        "tmp/{sample}.all_20kbflanking_merged.fasta"
    output:
        "tmp/{sample}_all_20kbflanking_removed.fasta"
    shell:
        r"""
        sed 's/(//;s/)//' {input} > {output}
        """

rule braker:
    input:
        "tmp/{sample}_all_20kbflanking_removed.fasta"
    output:
        gff3="tmp/{sample}_braker/braker.gff3",
        aa="tmp/{sample}_braker/braker.aa",
        codingseq="tmp/{sample}_braker/braker.codingseq"
    params:
        wd="tmp/{sample}_braker"
    log:
        "logs/braker/{sample}.log"
    shell:
        r"""
        perl {config[braker]}/braker.pl --gff3 --threads=15 --genome={input} --prot_seq={config[ref]} --workingdir={params.wd} >> {log} 2>&1
        """

rule prep_interproscan:
    input:
        "tmp/{sample}_braker/braker.aa"
    output:
        "tmp/{sample}_braker.aa"
    shell:
        r"""
        sed -e 's/\*//g' {input} > {output}
        """

# Interproscan
rule Interproscan_classification:
    input:
        aa="tmp/{sample}_braker.aa"
    output:
        interpro="result/{sample}.tsv"
    log:
        "logs/Interproscan_classification/{sample}.log"
    shell:
        r"""
        bash {config[interproscan]} -t p -appl Pfam, COILS, Gene3D -i {input.aa} -cpu 16 -f tsv -o {output.interpro} > {log} 2>&1
        """
        
# NLR_classification
rule NLR_classification:
    input:
        interpro="result/{sample}.tsv",
        fasta="tmp/{sample}_braker.aa",
        gff3="tmp/{sample}_braker/braker.gff3",
        genome="genome/{sample}.fasta",
        codingseq="tmp/{sample}_braker/braker.codingseq"
    output:
        nlr_list="result/{sample}_NLR.list",
        touch_file="result/{sample}_NLRclassification.done"
    params:
        prefix="{sample}",
        output_dir="result"
    log:
        "logs/Interproscan_classification/{sample}.log"
    shell:
        r"""
        NLR_classification.sh {input.interpro} {input.gff3} {input.fasta} {input.genome} {input.codingseq} --prefix {params.prefix} --output_dir {params.output_dir} >> {log} 2>&1
        touch {output.touch_file}
        """


# Find Integrate domain
rule Integrate_domain:
    input:
        "result/{sample}_NLR.list",
        "result/{sample}.tsv"
    output:
        "result/{sample}_pfam.tsv",
        "result/{sample}_NLR_Pfams.tsv",
        "result/{sample}_NLR_ID_Pfams.tsv"
    params:
        prefix="{sample}",
        output_dir="result",
        nlr_type="NLR"
    log:
        "logs/Integrate_domain/{sample}.log"
    shell:
        r"""
        get_integrated_domains.sh --prefix {params.prefix} --output_dir {params.output_dir} --type {params.nlr_type} > {log} 2>&1
        """
