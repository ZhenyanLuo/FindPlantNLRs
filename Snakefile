#Necessary packages:NLR-parser, NLR-annotator, BRAKER, samtools, blast, interproscan, etc.#
#use bash conda_install.sh command to install neccessary modules#
#Adapted from Tamene Tolessa's script#
SAMPLES = ["E_globulus"]
#Making a temp folder#
import os
path = 'tmp'
if not os.path.exists(path):
       os.mkdir(path)
result = 'result'
if not os.path.exists(result):
       os.mkdir(result)
#Remember to correct this since this should be the final output#
rule all:
     input:
       expand('tmp/{sample}.NLRparser.20kbflanking.fa', sample=SAMPLES)
#----------------------------------------------------------Start from NLR_annotator---------------------------------------------------------------------------------
#-----------------------------------------------------------------part 1--------------------------------------------------------------------------------------------
#Chopping the genome sequence into overlapping subsequences#
rule chop_sequence:
     input:
         "genome/{sample}.fa"
     output:
         "tmp/{sample}.choppedseq.fa"
     shell:
         "java -jar NLR-Annotator/ChopSequence.jar -i {input} -o {output} -l 20000 -p 5000"
#Searching the chopped subsequences for pre-determined NLR-associated motifs#
rule search_NLR_motifs:
     input:
         "tmp/{sample}.choppedseq.fa"
     output:
         "tmp/{sample}.NLRparser.xml"
     shell:
         "java -jar NLR-Annotator/NLR-Parser3.jar -t 10 -y ../meme/bin/mast -x NLR-Annotator/meme.xml -i {input} -c {output}"
#Generate the GFF format of NLR loci for the searched motifs#
rule NLR_annotator:
     input:
         "tmp/{sample}.NLRparser.xml"
     output:
         "tmp/{sample}.NLRparser.gff"
     shell:
         "java -jar NLR-Annotator/NLR-Annotator.jar -i {input} -g {output}"
#Indexing reference sequence#
rule index_ref_seq:
     input:
         "genome/{sample}.fa"
     output:
         "genome/{sample}.fa.fai"
     shell:
         "samtools faidx {input}"
#Create genome file. (?)#
rule convert_genome_format:
     input:
         "genome/{sample}.fa.fai"
     output:
         "genome/{sample}.genomefile"
     shell:
         "cut -d $'\t' -f1,2 {input} > {output}"
#Convert format of NLRparser.gff in order to incorporate with hmm output later# 
rule convert_NLRpaser:
      input:
          "tmp/{sample}.NLRparser.gff"
      output:
          "tmp/{sample}.NLRparser.bed"
      shell:
          """awk -v OFS='\\t' '{{if ($7 == "+") {{print $1, $4, $5, $1, "forward", $7}} else if ($7 == "-") print $1, $4, $5, $1, "reverse", $7}}' {input} >{output}"""
#Convert to 20kbflanking bed file with bedtools#
rule NLRpaser_20kbflanking:
       input:
          bed="tmp/{sample}.NLRparser.bed",
          genomefile="genome/{sample}.genomefile"
       output:
          "tmp/{sample}_NLRparser.20kbflanking.bed"
       shell:
          "bedtools slop -b 20000 -s -i {input.bed} -g {input.genomefile} | bedtools sort -i - |bedtools merge -s -d 1 -c 1,5,6 -o distinct,distinct,distinct, > {output}"
#Part 1 already tested and passed#              
              
#-------------------------------------------Use blast to identify genes which cannot be detected by NLR annotator pipeline------------------------------------------
#-----------------------------------------------------------------part 2--------------------------------------------------------------------------------------------
#Make a genome database for detecting nucleotide or protein query sequence#
#Build blastdb#
rule build_blast_database:
     input:
         "genome/{sample}.fa"
     output:
         "tmp/{sample}.database"
     shell:
         "makeblastdb -in {input} -dbtype nucl -parse_seqids -out {output}"
#Dectect whether there are genes which cannot be captured by using NLR-parser by using tblastn# remember to form a folder which include blastprotein#
rule run_tblastn:
     input:
         "genome/RefPlant_235NLR_ADR1_ZAR1_protein.fa"
     output:
         outfmt6="tmp/{sample}.tblastnout.outfmt6"
     params:
         "tmp/{sample}.database"    
     run:
         shell("tblastn -query {input} -db {params} -evalue 0.001 -outfmt 6 -o {output.outfmt6}")
#Convert tblastn file into bed, and add two coloums for strand information#
rule tblastn_to_bed:
     input:
         "tmp/{sample}.tblastnout.outfmt6"
     output:
         "tmp/{sample}.tblastnout.bed"
     shell:
         """awk -v OFS='\\t' '{{if ($10 - $9 > 0) {{print $2, $9, $10, $1, "forward", "+"}} else if ($10 - $9 < 0) print $2, $10, $9, $1, "reverse", "-"}}' {input} > {output}"""
##Generate 20kb flanking bed file for blast 
rule blast_20kb:
      input:
         bed="tmp/{sample}.tblastnout.bed",
         genomefile="genome/{sample}.genomefile"
     output:
         "tmp/{sample}.blast.20kbflanking.bed"
     shell:
         """ bedtools slop -b 20000 -s -i {input.bed} -g {input.genomefile} | bedtools sort -i - | bedtools merge -s -d 1 -c 1,5,6 -o distinct,distinct,distinct,  >  {output}"""      
#Adapted from Peri Tobias' s scripts------------------------------------------------------------------------------------------------
#Use nhmmer to search for conserved nucleotide binding domain shared by Apaf-1, Resistance proteins and CED4 from coiled-coil NLR and TIR NLR sequences#
#-----------------------------------------------------------------part 3--------------------------------------------------------------------------------------------
#Peri has already prepared hmm profiles which are named as EG_nonTIRhmm and EG_TIRhmm respectively. 
rule find_TIR:
     input:
         nonTIR="genome/EG_nonTIRhmm",
         TIR="genome/EG_TIRhmm",
         genome="genome/{sample}.fa"
     output:
         TIR_out="tmp/{sample}.TIRout",
         nonTIR_out="tmp/{sample}.nonTIRout"
     run:
         shell("nhmmer {input.nonTIR} {input.genome} > {output.nonTIR_out}")
         shell("nhmmer {input.TIR} {input.genome} > {output.TIR_out}")
#convert both nhmmer output into bed file by using awk script#
rule awk_convert:
     input:
         TIR_out="tmp/{sample}.TIRout",
         nonTIR_out="tmp/{sample}.nonTIRout"
     output:
         TIRoutbed="tmp/{sample}.TIRout.bed",
         nonTIRoutbed="tmp/{sample}.nonTIRout.bed"
     run:
         shell("awk -f Peris_NLR/Myrtaceae_NLR_workflow/make_bed_hmmOut.awk {input.TIR_out} > {output.TIRoutbed}")
         shell("awk -f Peris_NLR/Myrtaceae_NLR_workflow/make_bed_hmmOut.awk {input.nonTIR_out} > {output.nonTIRoutbed}")
#convert both bed file into fasta file by using bedtools#
rule bedtools:
     input:
         TIR_bed="tmp/{sample}.TIRout.bed",
         nonTIR_bed="tmp/{sample}.nonTIRout.bed",
         genome="genome/{sample}.fa"
     output:
         TIR_fasta="tmp/{sample}.TIR.fasta",
         nonTIR_fasta="tmp/{sample}.nonTIR.fasta"
     run:
         shell("bedtools getfasta -s -fi {input.genome} -bed {input.TIR_bed} -fo {output.TIR_fasta}")
         shell("bedtools getfasta -s -fi {input.genome} -bed {input.nonTIR_bed} -fo {output.nonTIR_fasta}")
#Extract first 200 sequences from previously output nonTIR and TIR fasta file, change 200 into other number when neccessary (why?#
rule awk:
     input:
         nonTIR="tmp/{sample}.nonTIR.fasta",
         TIR="tmp/{sample}.TIR.fasta"
     output:
         nonTIR_200="tmp/{sample}.nonTIR_200.fasta",
         TIR_200="tmp/{sample}.TIR_200.fasta",
         NBARC400="tmp/{sample}_NBARC_400.fasta"
     run:
         shell("""awk "/^>/ {{n++}} n>200 {{exit}} 1" {input.nonTIR} > {output.nonTIR_200}""") 
         shell("""awk "/^>/ {{n++}} n>200 {{exit}} 1" {input.TIR} > {output.TIR_200}""")
         shell("cat {output.nonTIR_200} {output.TIR_200} > {output.NBARC400}")
#Remove duplicated sequences#
rule seqkit:
     input:
         "tmp/{sample}_NBARC_400.fasta"
     output:
         "tmp/{sample}_NBARC.fasta"
     shell:
         "seqkit rmdup -D duplicates -n {input} > {output}"
#Generate alignment by using clustalo#
rule clustalo:
     input:
         "tmp/{sample}_NBARC.fasta"
     output:
         "tmp/{sample}_NBARC.sto"
     shell:
         "clustalo -i {input} -o {output} --outfmt=st"
#Build a profile hmm from an alignment#
rule hmmbuild:
     input:
         "tmp/{sample}_NBARC.sto"
     output:
         "tmp/{sample}.hmm"
     shell:
         "hmmbuild -nucleic {output} {input}"
#Use hmm profile built to search queries against genome#
rule nhmmer:
     input:
         hmm="tmp/{sample}.hmm",
         genome="genome/{sample}.fa"
     output:
         "tmp/{sample}_NBARCout"
     shell:
         "nhmmer {input.hmm} {input.genome} > {output}"
#Convert nhmmer output into bed file#
#Required double-check#
rule make_bed_hmmout:
     input:
         NBARC="tmp/{sample}_NBARCout"
     output:
         "tmp/{sample}_NBARC.bed"
     shell:
         "awk -f Peris_NLR/Myrtaceae_NLR_workflow/make_bed_hmmOut.awk {input.NBARC} > {output}"     
#Get 20kb upstream and downstream# 
rule NBARC_20flanking:
       input:
          bed="tmp/{sample}_NBARC.bed",
          genomefile="genome/{sample}.genomefile"
       output:
          "tmp/{sample}_NBARC.20kbflanking.bed"
       shell:
          "bedtools slop -b 20000 -s -i {input.bed} -g {input.genomefile} >  {output}"
#---------------------------------------Now we have output from hmm, blast and NLR_annotator, combine them into one file--------------------------------------------
#-----------------------------------------------------------------part 4--------------------------------------------------------------------------------------------
##Combine the output together#
#Include blast file after getting the query file#
rule combine_all_bed:
         input:
             hmm="tmp/{sample}_NBARC.20kbflanking.bed",
             annotator="tmp/{sample}_NLRparser.20kbflanking.bed",
             blast="tmp/{sample}.blast.20kbflanking.bed"
         output:
             "tmp/{sample}_all20kbflanking.bed"
         shell:
             "cat {input.hmm} {input.annotator} {input.blast}|sort -k1,1 -k2,2n >{output}"
#Merge the bed file now#
rule merge_all_20kbflanking:
     input:
         "tmp/{sample}_all20kbflanking.bed"
     output:
         "tmp/{sample}.all_20kbflanking_merged.bed"
     shell:
         "bedtools merge -s -d 1 -c 1,5,6 -o distinct,distinct,distinct, -i {input} > {output}"               
#Convert bedfile into fasta#
rule convert_20kbflankingbedfile_fasta:
     input:
         genome="genome/{sample}.fa",
         bed="tmp/{sample}.all_20kbflanking_merged.bed"
     output:
          "tmp/{sample}.all_20kbflanking_merged.fasta"
     shell:
          "bedtools getfasta -s -fi {input.genome} -bed {input.bed} -fo {output}"
#Convert all the sequences in 20kb flanking fasta into uppercase (not sure)#
rule convert_format:
     input:
         "tmp/{sample}.all_20kbflanking_merged.fasta"
     output:
         "tmp/{sample}.all_20kbflanking_merged_upper.fasta"
     shell:
         "cat {input} | awk '/^>/ {{print($0)}}; /^[^>]/ {{print(toupper($0))}}' > {output}"   
#Maybe use 20kbflanking.fa instead of NBARC_nt.fasta#
#Translate nucleotide NBARC sequeces including extended sequences#
rule translate:
     input:
         "tmp/{sample}.all_20kbflanking_merged_upper.fasta"
     output: 
         "tmp/{sample}.all_20kbflanking_merged_upper.faa"
     shell:  
         "Peris_NLR/Myrtaceae_NLR_workflow/translate.py {input} {output}"
#----------------------------------To classify the output of annotator and hmm#
#Remove * in stop codon, otherwise interproscan will not work#
rule remove_stop_codon:
     input:
         "tmp/{sample}.all_20kbflanking_merged_upper.faa"
     shell:
         "sed -i 's/*//g' {input}"
#Run Interproscan, database options: Pfam, coils, gene3D #
rule Interproscan:
     input:
         "tmp/{sample}.all_20kbflanking_merged_upper.faa"
     params:
         "tmp/"
     shell:
          "./interproscan/interproscan-5.50-84.0/interproscan.sh -t p -appl Pfam,COILS,Gene3D -i {input} -f tsv,gff3 -d {params}"
#Predict genes by using braker, remove special header first##Remember to use extended one#
#RefPlantNLR_aa.fa is from https://www.biorxiv.org/content/10.1101/2020.07.08.193961v2#
#Remember to correct the path for braker.pl#
####Please double check#
#This step is modified from Peri's script: braker_nlr.pbs#
#Remove sample species from config file of braker if you stopped once#
rule braker:
     input:
         raw="tmp/{sample}.all_20kbflanking_merged_upper.fasta",
         ref="genome/prothint_sequences.fa"
     output:
         removed="tmp/{sample}_all_20kbflanking_removed.fasta",
         hints_gtf="BRAKER/scripts/braker/{sample)_augustus.hints.gtf",
         gff3="BRAKER/scripts/braker/{sample}_augustus_out.gff3",
         hints_aa="BRAKER/scripts/braker/{sample}_augustus.hints.aa"
     params:
         "{sample}"
     run:
         shell("sed 's/(//;s/)//' {input.raw} > {output.removed}")
         shell("./scripts/braker.pl --cores=15 --genome={output.removed} --prot_seq={input.ref} --epmode --species={params} --gff3")
         shell("./Augustus/scripts/gtf2gff.pl <{output.hints_gtf} --printExon --out={output.gff3} --gff3")
         shell("mv scripts/braker/augustus.hints.gtf {output.hints_gtf}")
         shell("mv scripts/braker/augustus.hints.aa {output.hints_aa}")
         shell("mv scripts/braker/braker.gff3 scripts/braker/{sample}_braker.gff3")
         shell("rm -r scripts/braker/GeneMark-EP")
         shell("rm -r scirpts/braker/GeneMark-ES")
#Get fasta data from the braker.gff file#
#May use getAnnoFastaFromJoingenes.py in augustus foler instead#
rule braker_gff_to_fasta:
       input:
         genome="tmp/{sample}.all_20kbflanking_merged_upper.fasta",
         bed="scripts/braker/{sample}_braker.gtf"
       params:
         "tmp/{sample}_braker"   
       shell:
          "./Augustus/scripts/getAnnoFastaFromJoingenes.py -g {input.genome} -f {input.bed} -o {params}     
#Remove special characters and rename the augustus output#
rule braker_step2:
     input:
          "tmp/{sample}_braker.aa"      
     output:
          "tmp/{sample}_braker.faa"
     shell:
          "sed s/\*//g {input} > {output}"
###Tamene's version, use PF00931 and Pfam-A.hmm----------------------------------------------------------------
#Let's start from using hmm profile built ({sample}.hmm) #
rule hmmsearch:
     input:
           hint="tmp/{sample}_braker.aa" 
     output:
           noali="tmp/{sample}.NB-ARC_tblout_noali.txt",
           tblout="tmp/{sample}.NB-ARC_hmmsearch_tblout.perseqhit.txt"
     params:
           "Pfam/PF00931.hmm"
     shell:
           hmmsearch -o {output.noali} --tblout {out.tblout} --noali --notextw {params} {input.hint}"
#Sorting out and cutting of sequence IDs of domains
rule grep_hmmsearch:
     input:
           "tmp/{sample}.NB-ARC_hmmsearch_tblout.perseqhit.txt"
     output:
           "tmp/{sample}.NB-ARC_hmmsearch_perseqhit_seqID.txt"
     shell:
           "grep -v '#' {input} | sort -k5,5n | cut -d ' ' -f1 | uniq > {output}"
#Using the sequence ID and seqtk to pull out protein sequences of NB-ARC domains
rule seqtk:
     input:
           seqID="tmp/{sample}.NB-ARC_hmmsearch_perseqhit_seqID.txt",
           hint="BRAKER/scripts/braker/{sample}_augustus.hints.aa"
     output:
           "tmp/{sample}.NB-ARC_hmmsearch_perseqhit_protein.fa"
     shell:   
           "seqtk subseq {input.hint} {input.sedID} > {output}"
#Removing asterisk (*) from the end of each sequences
rule sed:
     input:
           "tmp/{sample}.NB-ARC_hmmsearch_perseqhit_protein.fa"    
     shell:
           "sed -i 's/*//g' {input}"
#Generate empty directory to store Interproscan output
rule mkdir:
     output:
           "result/Interpro_{sample}"
     shell:
           "mkdir -p {output}"
#Using InterProScan to detect the TIR, LRR and COIL sub-classes of NB-ARC domains
#Remember to update the path of interproscan.sh
rule interproscan_NBARC:
     input:
           "tmp/{sample}.NB-ARC_hmmsearch_perseqhit_protein.fa"    
     output:
           "result/Interpro_{sample}" 
     shell:
           "./interproscan/interproscan-5.50-84.0/interproscan.sh -t p -appl Pfam, COILS, Gene3D -i {input} -cpu 16 -f tsv, gff3 -d {output}"
#Search NB-ARC domain against library of Pfam
rule pfam_scan:
     input:
           "tmp/{sample}.NB-ARC_hmmsearch_perseqhit_protein.fa"
     output:
           mkdir="./Pfam_{sample}",
           Pfamscan="tmp/{sample}.protein.fa_pfamscan.txt"
     run:
           shell("mkdir -p {output.mkdir}")
           shell("./script/pfam_scan.pl -fasta {input} -dir ~/ddatabase/PfamScan -as -cpu 16 -outfile {output.Pfamscan}")
#Parsing the output of PfamScan output parser using the script 
#K-parse_Pfam_domains_v3.1.pl from https://github.com/krasileva-group/plant_rgenes is used in this step, ref: https://bmcbiol.biomedcentral.com/articles/10.1186/s12915-016-0228-7
rule K_parse:
     input:
           "tmp/sample}.protein.fa_pfamscan.txt"   
     output:
           "tmp/{sample}.protein.fa_pfamscan_parsed.verbose"
     shell:
           "perl ~/scripts/plant_rgenes/processing_scripts/K-parse_Pfam_domains_v3.1.pl --pfam {input} --evalue 0.001 --output {output} --verbose T"
#Parsing the output of PfamScan output parser using the script "K-parse_Pfam_domains_NLR-fusions-v2.4.2.pl"#
#K-parse_Pfam_domains_NLR-fusions-v2.4.1.pl from https://github.com/krasileva-group/plant_rgenes is used in this step, ref: https://bmcbiol.biomedcentral.com/articles/10.1186/s12915-016-0228-7               
#Remember to make a db_descriptions.txt file in genome folder
rule K_parse_fusion:
      input:
           verbose="tmp/{sample}.protein.fa_pfamscan_parsed.verbose",
           db="genome/db_descriptions.txt"
      output:
           mkdir="./Pfam_{sample}_parser"
      run:
           shell("mkdir -p {output.mkdir}")
           shell("mv {input.verbose} {output.mkdir}")
           shell("perl ~/scripts/plant_rgenes/processing_scripts/K-parse_Pfam_domains_NLR-fusions-v2.4.1.pl --indir {output.mkdir} --evalue 0.001 -o {output.mkdir} -d {input.db}")

