#Necessary packages:NLR-parser, NLR-annotator, BRAKER, samtools, blast, interproscan, etc.#
#use bash conda_install.sh command to install neccessary modules#
#Adapted from Tamene Tolessa's script#
SAMPLES = ["name_of_genome"]
#Making a temp folder#
import os
path = 'tmp'
if not os.path.exists(path):
       os.mkdir(path)
os.mkdir(result)
rule all:
     input:
       expand('genome/{sample}.fa', sample=SAMPLES)
#---------------Start from NLR_annotator-------------------------------------
#Chopping the genome sequence into overlapping subsequences#
rule chop_sequence:
     input: 
         "genome/{sample}.fa"
     output:
         "tmp/{sample}.choppedseq.fa"
     shell: 
         "java -jar ~/NLR-parser/scripts/ChopSequence.jar -i {input} -o {output} -l 20000 -p 5000" 
#Searching the chopped subsequences for pre-determined NLR-associated motifs#
rule search_NLR_motifs:
     input:
         "tmp/{sample}.choppedseq.fa"
     output:
         "tmp/{sample}.NLRparser.xml"
     shell:
         "java -jar ~/NLR-parser/scripts/NLR-Parser.jar -t 10 \
        -y ~/anaconda3/envs/NLR_Annotator/bin/meme_4.9.1/bin/mast \
        -x ~/NLR-parser/scripts/meme.xml -i {input} \
        -c {output}"
#Generate the GFF format of NLR loci for the searched motifs#
rule NLR_annotator:
     input:
         "tmp/{sample}.NLRparser.xml"
     output:
         "tmp/{sample}.NLRparser.gff"
     shell:
         "java -jar ~/NLR-parser/scripts/NLR-Annotator.jar -i {input} \
        -g {output}"
#Convert tblastn file into bed, get coloumn 1 2 9 10#
rule tblastn_to_bed:
     input:
         "tmp/{sample}.tblastnout.outfmt6"
     output:
         "tmp/{sample}.tblastnout.bed"
     shell:
         """cat {input} | awk "{{print $1"\\t"$2"\\t"$9"\\t"$10}}" > {output}"""
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
#Generate 20kb flanking BED file for blastx file#
rule generate_20kb_flanking_bed_for_blastx:
     input:
         bed="tmp/{sample}.tblastnout.bed",
         genomefile="genome/{sample}.genomefile"
     output:
         "tmp/{sample}.tblastn.20kbflanking.bed"
     shell:
         """bedtools slop -b 20000 -s -i {input.bed} -g {input.genomefile} | bedtools sort -i - | bedtools merge -s -d 100 -i - > {output}"""
#Generate 20kb flanking BED file for NLR-parser file#
rule generate_20kb_flanking _bed_for_NLR-parser:
     input:
         gff="tmp/{sample}.NLRparser.gff",
         genomefile="genome/{sample}.genomefile"
     output:
         "tmp/{sample}.NLRparser.20kbflanking.bed"
     shell:
        """ bedtools slop -b 20000 -s -i {input.gff} -g {input.genomefile} | bedtools sort -i - | bedtools merge -s -d 100 -i - >  {output}"""
#Merge the two bed files (combine blastn.bed and NLRparser.bed into one BED file#
rule merge_bed:
     input:
         tblastn="tmp/{sample}.tblastn.20kbflanking.bed",
         NLRparser="tmp/{sample}.NLRparser.20kbflanking.bed"
     output:
         "tmp/{sample}.all.20kbflanking.bed"
     shell:
         "cat {input.tblastn} {input.NLRparser} | bedtools sort -i - | bedtools merge -d 100 -i - > {output}"
#Convert the merged bed file into fasta format (? required double check)#
rule bed_to_fasta:
     input:
         genome="genome/{sample}.fa",
         flankingbed="tmp/{sample}.all.20kbflanking.bed"
     output:
         "tmp/{sample}.all.20kbflanking.fa"  
     shell:
         "bedtools getfasta -fi {input.genome} -bed {input.flankingbed} > {output}"
#Convert all the sequences in 20kb flanking fasta into uppercase (not sure)#
rule convert_format:
     input:
         "tmp/{sample}.all.20kbflanking.fa",
     output:
         "tmp/{sample}.all.20kbflanking_upper.fa"
     shell:
         "cat {input} | awk '/^>/ {{print($0)}; /^[^>]/ {print(toupper($0))}}' > {output}"         
#Gene prediction by BRAKER using extended regions around NB-ARCs by 20kb up and downsream##############################################################
#rule predict_by_braker:
#    input:
#         genome="/genome/{sample}.all_20kbflanking_upper.fa",
#         prot="tmp/prothint_sequences.fa"
#    # output:
#     shell:
#         "braker.pl --genome={input.genome} \
#        --prot_seq={input.prot} \
#        --species={sample} --epmode --cores=15 --softmasking --prg=ph \
#        --ALIGNMENT_TOOL_PATH=~/anaconda3/envs/braker2/bin/spaln --gff3"
#
#Adapted from Peri Tobias' s scripts------------------------------------------------------------------------------------------------
#Use nhmmer to search for conserved nucleotide binding domain shared by Apaf-1, Resistance proteins and CED4 from coiled-coil NLR and TIR NLR sequences#
#Peri has already prepared hmm profiles which are named as EG_nonTIRhmm and EG_TIRhmm respectively. 
rule find_TIR:
     input:
         nonTIR="genome/EG_nonTIRhmm",
         TIR="genome/EG_TIRhmm",
         genome="genome/{sample}.fa"
     output:
         TIR_out="tmp/{sample}.TIRout"
         nonTIR_out="tmp/{sample}.nonTIRout"
     run:
         shell("nhmmer {input.nonTIR} {input.genome} > {output.nonTIR_out}")
         shell("nhmmer {input.TIR} {input.genome} > {output.TIR_out}")
#convert both nhmmer output into bed file by using awk script#
rule awk:
     input:
         TIRout="tmp/{sample}.TIRout",
         nonTIRout="tmp/{sample}.nonTIRout"
     output:
         TIRoutbed="tmp/{sample}.TIRout.bed",
         nonTIRoutbed="tmp/{sample}.TIRout.bed"
     run:
         shell("awk -f make_bed_hmmOut.awk {input.TIRout} > {output.TIRoutbed}")
         shell("awk -f make_bed_hmmOut.awk {input.nonTIRout} > {output.nonTIRoutbed}")
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
         shell("awk "/^>/ {n++} n>200 {exit} 1" {input.nonTIR} > {output.nonTIR_200}") 
         shell("awk "/^>/ {n++} n>200 {exit} 1" {input.TIR} > {output.TIR_200}")
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
         genome="/genome/{sample}.fa"
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
         "awk -F script/make_bed_hmmOut.awk {input.NBARC} > {output}"
#Extract sequences from bed file#
#rule bedtools_NBARC:
#     input:
#         genome="/genome/{sample}.fa",
#         bed="tmp/{sample}_NBARC.bed"
#     output:
#         "tmp/{sample}_NBARC_nt.fasta"
#     shell:
#         "bedtools getfasta -s -fi {input.genome} -bed {input.bed} > {output}"
####
#Extend 20kb upstream and downstream##Remember to double check whether bed file can be used to generate 20kb flanking bed
rule generate_20kb_flanking_bed_for_NBARC:
     input:
         bed="tmp/{sample}_NBARC.bed",
         genomefile="genome/{sample}.genomefile"
     output:
         "tmp/{sample}.NBARC.20kbflanking.bed"
     shell:
        """ bedtools slop -b 20000 -s -i {input.bed} -g {input.genomefile} | bedtools sort -i - | bedtools merge -s -d 100 -i - >  {output}"""
#Convert 20kb flanking bed into fasta file#
rule bed_to_fasta:
     input:
         genome="genome/{sample}.fa",
         flankingbed="tmp/{sample}.NBARC.20kbflanking.bed"
     output:
         "tmp/{sample}.NBARC.20kbflanking.fa"
     shell:
         "bedtools getfasta -fi {input.genome} -bed {input.flankingbed} > {output}"
#Convert all the sequences in 20kb flanking fasta into uppercase (not sure)#
rule convert_format:
     input:
         "tmp/{sample}.NBARC.20kbflanking.fa",
     output:
         "tmp/{sample}.NBARC.20kbflanking_upper.fa"
     shell:
         "cat {input} | awk '/^>/ {{print($0)}; /^[^>]/ {print(toupper($0))}}' > {output}"   
####Combine with NLR_annotator output, and remove duplicated sequences#
rule combine:
     input:
         fasta_annotator="tmp/{sample}.all.20kbflanking.fa",
         fasta_hmm="tmp/{sample}.NBARC.20kbflanking_upper.fa"
     output:
         "tmp/{sample}.20kbflanking.fa
     shell:
         "cat {input.fasta_annotator} {input.fasta_hmm} > {output}"
#Maybe use 20kbflanking.fa instead of NBARC_nt.fasta#
#Translate nucleotide NBARC sequeces including extended sequences#
rule translate:
     input:
         "tmp/{sample}_NBARC_nt.fasta"
     output: 
         "tmp/{sample}_NBARC_aa.fasta"
     shell:  
         "script/translate.py {input} {output}"
#----------------------------------Filt#
#Remove * in stop codon, otherwise interproscan will not work#
rule remove_stop_codon:
     input:
         "tmp/{sample}_NBARC_nt.fasta"
     output:
         "tmp/{sample}_NBARC.faa"
     shell:
         "sed 's/*//g' {input} > {output}"
#Run Interproscan, database options: Pfam, coils, gene3D #
rule Interproscan:
     input:
         "tmp/{sample}_NBARC.faa"
     output:
         "tmp/{sample}.tsv
     shell:
          "./interproscan.sh -t p -appl Pfam,COILS,Gene3D -i {input} -f tsv,gff3 -o {output}"
#Predict genes by using braker, remove special header first##Remember to use extended one#
#RefPlantNLR_aa.fa is from https://www.biorxiv.org/content/10.1101/2020.07.08.193961v2#
#Remember to correct the path for braker.pl#
rule braker:
     input:
         raw="tmp/{sample}_NBARC_20kb.fasta",
         genome="genome/{sample}.fa",
         ref="genome/RefPlantNLR_aa.fa"
     output:
         removed="tmp/{sample}_NBARC_20kb_removed.fasta"
     run:
         shell("sed 's/(//;s/)//' {input.raw} > {output.removed}")
         shell("script/braker.pl --cores=6 --genome={input.genome} --prot_seq={input.ref} --ALIGNMENT_TOOL_PATH=/usr/local/genemark-es/4.59/ProtHint/bin/ --prg=ph --epmode --species={sample}
/usr/local/augustus/3.3.3/scripts/gtf2gff.pl /braker/augustus.hints.gtf --printExon --out=/braker/augustus.gff3 --gff3")
#Remove special characters and rename the augustus output#
rule braker_step2:
     input:
          "tmp/augustus.hints.aa"      
     output:
          "tmp/{sample}_augustus_aa.fasta"
     shell:
          "sed s/\*//g {input} > {output}"
#Filt the output based on domain identified#---------------------------------------Peris' s version-----------------------------------------------------------------
#PF00931 = NB-ARC domain, G3DSA:3.40.50.300 = P-loop containing nucleoside triphosphate hydrolase, PF08263 = LRR, PF12799 = LRR, PF13306 = LRR, PF13855 = LRR, PF13516 = LRR,
#G3DSA:1.10.8.430 = Gene3D LRR. Also need to search for this Gene3D output in interproscan because Pfam does not recognise all LRR signals.
#Also look for TIR domains, Coils etc.
#PF01582 = TIR#Combine the result of interproscan and braker together to identify NBARC#
#Identify NB_ARC: 
rule combine_interproscan_braker:
     input:
          tsv="tmp/{sample}.tsv",
          gff3="tmp/augustus_out.gff3"
     output:
          "result/{sample}_NBARC.gff3"     
     shell:
          "run_NBARC.sh {input.tsv} {input.gff3} > {output}"
#Identify TIR_NB: 
rule combine_interproscan_braker:
     input:
          tsv="tmp/{sample}.tsv",
          gff3="tmp/augustus_out.gff3"
     output:
          "result/{sample}_TIRNB.gff3"     
     shell:
          "run_TIRNB.sh {input.tsv} {input.gff3} > {output}"
#Identify NB_LRR: 
rule combine_interproscan_braker:
     input:
          tsv="tmp/{sample}.tsv",
          gff3="tmp/augustus_out.gff3"
     output:
          "result/{sample}_NBLRR.gff3"     
     shell:
          "run_NBLRR.sh {input.tsv} {input.gff3} > {output}"
#Identify TIR: 
rule combine_interproscan_braker:
     input:
          tsv="tmp/{sample}.tsv",
          gff3="tmp/augustus_out.gff3"
     output:
          "result/{sample}_TIR.gff3"     
     shell:
          "run_TIR.sh {input.tsv} {input.gff3} > {output}"
###Tamene's version, original one used PF00931 and Pfam-A.hmm----------------------------------------------------------------
#Let's start from using hmm profile built ({sample}.hmm) #
#add augustus.hints.aa#
rule hmmsearch:
     input:
           hmm="tmp/{sample}.hmm",
           hint="braker/augustus.hints.aa"
     output:
           noali="tmp/{sample}.NB-ARC_tblout_noali.txt",
           tblout="tmp/{sample}.NB-ARC_hmmsearch_tblout.perseqhit.txt"
     shell:
           "hmmsearch -o {output.noali} --tblout {output.tblout} --noali --notextw {input.hmm} {input.hint}"
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
           "tmp/{sample}.NB-ARC_hmmsearch_perseqhit_seqID.txt"   
     output:
           "tmp/{sample}.NB-ARC_hmmsearch_perseqhit_protein.fa"
     shell:   
           "seqtk subseq ./braker/augustus.hints.aa {input} > {output}"
#Removing asterisk (*) from the end of each sequences
rule sed:
     input:
           "tmp/{sample}.NB-ARC_hmmsearch_perseqhit_protein.fa"    
     shell:
           "sed -i 's/*//g' {input}"
#Generate empty directory to store Interproscan output
rule mkdir:
     output:
           "Interpro_{sample}"
     shell:
           "mkdir -p {output}"
#Using InterProScan to detect the TIR, LRR and COIL sub-classes of NB-ARC domains
#Remember to update the path of interproscan.sh
rule interproscan_NBARC:
     input:
           "tmp/{sample}.NB-ARC_hmmsearch_perseqhit_protein.fa"    
     output:
           "./Interpro_{sample}" 
     shell:
           "interproscan/interproscan.sh -t p -appl Pfam, COILS, Gene3D -i {input} -cpu 16 -f tsv, gff3 -d {output}"
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
#Run the filt script#
rule filt_script:
      input:
           

