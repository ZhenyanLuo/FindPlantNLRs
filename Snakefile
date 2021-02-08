#Necessary packages:NLR-parser, NLR-annotator, BRAKER, samtools, blast, interproscan#

SAMPLES = ["name_of_genome"]
#Making a temp folder#
import os
path = 'tmp'
if not os.path.exists(path):
       os.mkdir(path)

rule all:
     input:
       expand('genome/{sample}.fa', sample=SAMPLES)
#Chopping the genome sequence into overlapping subsequences#
rule chop_sequence:
     input: 
         fa="genome/{sample}.fa"
     output:
         "tmp/{sample}.choppedseq.fa"
     shell: 
         "java -jar ~/NLR-parser/scripts/ChopSequence.jar -i {input.fa} -o {output} \
         -l 20000 -p 5000" 

#Searching the chopped subsequences for pre-determined NLR-associated motifs#
rule step2:
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
rule step3:
     input:
         "tmp/{sample}.NLRparser.xml"
     output:
         "tmp/{sample}.NLRparser.gff"
     shell:
         "java -jar ~/NLR-parser/scripts/NLR-Annotator.jar -i {input} \
        -g {output}"
#Make a genome database for detecting nucleotide or protein query sequence#
rule step4:
     input:
         fa="genome/{sample}.fa"
     output:
         "tmp/{sample}.genome_nucl_database"
     shell:
         "makeblastdb -in {input.fa} -dbtype nucl -parse_seqids \
        -out {output}"
#Dectect whether there are genes which cannot be captured by using NLR-parser by using tblastn#
#remember to form a folder which include blastprotein#
rule step5:
     input:
         blastprotein="blastprotein/blastprotein",
         genomebase="tmp/{sample}.genome_nucl_database"
     output:
         "tmp/{sample}.tblastnout.outfmt6"
     shell:
         "tblastn -query {input.blastprotein} -db {input.genomebase} -evalue 0.001 \
         -outfmt 6 > {output}"
#Convert tblastn file into bed, get coloumn 1 2 9 10#
rule step6:
     input:
         "tmp/{sample}.tblastnout.outfmt6"
     output:
         "tmp/{sample}.tblastnout.bed"
     shell:
         """cat {input} | awk "{{print $1"\\t"$2"\\t"$9"\\t"$10}}" > {output}"""
#Indexing reference sequence#
rule step7:
     input:
         "genome/{sample}.fa"
     output:
         "genome/{sample}.fa.fai"
     shell:
         "samtools faidx {input}" 
#Create genome file. (?)#
rule step8:
     input:
         "genome/{sample}.fa.fai"
     output:
         "genome/{sample}.genomefile"
     shell:
         "cut -d $'\t' -f1,2 {input} > {output}"
#Generate 20kb flanking BED file for blastx file#
rule step9:
     input:
         bed="tmp/{sample}.tblastnout.bed",
         genomefile="genome/{sample}.genomefile"
     output:
         "tmp/{sample}.tblastn.20kbflanking.bed"
     shell:
         """bedtools slop -b 20000 -s -i {input.bed} \
        -g {input.genomefile} | bedtools sort -i - | bedtools merge \
        -s -d 100 -i - > {output}"""
#Generate 20kb flanking BED file for NLR-parser file#
rule step10:
     input:
         gff="tmp/{sample}.NLRparser.gff",
         genomefile="genome/{sample}.genomefile"
     output:
         "tmp/{sample}.NLRparser.20kbflanking.bed"
     shell:
        """ bedtools slop -b 20000 -s -i {input.gff} -g {input.genomefile} \
        | bedtools sort -i - | bedtools merge -s -d 100 -i - \
        >  {output}"""
#Merge the two bed files (combine blastn.bed and NLRparser.bed into one BED file#
rule step11:
     input:
         tblastn="tmp/{sample}.tblastn.20kbflanking.bed",
         NLRparser="tmp/{sample}.NLRparser.20kbflanking.bed"
     output:
         "tmp/{sample}.all.20kbflanking.bed"
     shell:
         "cat {input.tblastn} {input.NLRparser} \
        | bedtools sort -i - \
        | bedtools merge -d 100 -i - > {output}"
#Convert the merged bed file into fasta format (? required double check)#
rule step12:
     input:
         genome="genome/{sample}.fa",
         flankingbed="tmp/{sample}.all.20kbflanking.bed"
     output:
         "tmp/{sample}.all.20kbflanking.fa"  
     shell:
         "bedtools getfasta -fi {input.genome} -bed {input.flankingbed} \
        > {output}"
#Convert all the sequences in 20kb flanking fasta into uppercase (not sure)#
rule step13:
     input:
         "tmp/{sample}.all.20kbflanking.fa",
     output:
         "tmp/{sample}.all.20kbflanking_upper.fa"
     shell:
         "awk '/^>/ {{print($0)}; /^[^>]/ {print(toupper($0))}}'{input}>{output}"         
#Gene prediction by BRAKER using extended regions around NB-ARCs by 20kb up and downsream#
rule step14:
     input:
         genome="/genome/{sample}.all_20kbflanking_upper.fa",
         prot="tmp/prothint_sequences.fa"
    # output:
     shell:
         "braker.pl --genome={input.genome} \
        --prot_seq={input.prot} \
        --species={sample} --epmode --cores=15 --softmasking --prg=ph \
        --ALIGNMENT_TOOL_PATH=~/anaconda3/envs/braker2/bin/spaln --gff3"

#------------------------------------------------------------------------------------------------
#use nhmmer to search for conserved nucleotide binding domain shared by Apaf-1, Resistance proteins and CED4 from coiled-coil NLR and TIR NLR sequences#
rule find_nonTIR:
     input:
         nonTIR="/genome/EG_nonTIRhmm",
         genome="/genome/{sample}.fa"
     output:
         "tmp/{sample}.nonTIRout"
     shell:
         "nhmmer {input.nonTIR} {input.genome} > {output}"
#use nhmmer to search for conserved nucleotide binding domain shared by Apaf-1, Resistance proteins and CED4 from coiled-coil NLR and TIR NLR sequences#
rule find_TIR:
     input：
         TIR="/genome/EG_TIRhmm",
         genome="/genome/{sample}.fa"
     output:
         "tmp/{sample}.TIRout"
     shell:
         "nhmmer {input.TIR} {input.genome} > {output}"
#convert both nhmmer output into bed file by using awk script#
rule awk1:
     input:
         "tmp/{sample}.TIRout"
     output:
         "tmp/{sample}.TIRout.bed"
     shell:
         "awk -f make_bed_hmmOut.awk {input} > {output}"
rule awk2:
     input:
         "tmp/{sample}.nonTIRout"
     output:
         "tmp/{sample}.nonTIRout.bed"
     shell:
         "awk -f make_bed_hmmOut.awk {input} > {output}"
#convert both bed file into fasta file by using bedtools#
rule bedtools_1:
     input:
         TIR_bed="tmp/{sample}.TIRout.bed",
         genome="/genome/{sample}.fa"
     output:
         "tmp/{sample}.TIR.fasta"
     shell:
         "bedtools getfasta -s -fi {input.genome} -bed {input.TIR_bed} -fo {output}"
rule bedtools_2:
     input:
         nonTIR_bed="tmp/{sample}.nonTIRout.bed",
         genome="/genome/{sample}.fa"
     output:
         "tmp/{sample}.nonTIR.fasta"
     shell:
         "bedtools getfasta -s -fi {input.genome} -bed {input.nonTIR_bed} -fo {output}"
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
rule culstalo:
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
         "hummbuild -nucleic {output} {input}"
#Search queries against genome#
rule nhmmer:
     input:
         hmm="tmp/{sample}.hmm",
         genome="/genome/{sample}.fa"
     output:
         "tmp/{sample}_NBARCout"
     shell:
         "nhmmer {input.hmm} {input.genome} > {output}"
#Convert nhmmer output into bed file#
rule make_bed_hmmout:
     input:
         NBARC="tmp/{sample}_NBARCout",
         awk_script="script/make_bed_hmmOut.awk"
     output:
         "tmp/{sample}_NBARC.bed"
     shell:
         "awk -F {input.awk_script} {input.NBARC} > {output}"
#Extract sequences from bed file#
rule bedtools_NBARC:
     input:
         genome="/genome/{sample}.fa",
         bed="tmp/{sample}_NBARC.bed"
     output:
         "tmp/{sample}_NBARC_nt.fasta"
     shell:
         "bedtools getfasta -s -fi {input.genome} -bed {input.bed} > {output}"
####
####Combine with NLR_annotator output, and remove duplicated sequences#
rule combine:
     input:
         fasta_annotator="tmp/{sample}.all.20kbflanking.fa",
         fasta_hmm="tmp/{sample}_NBARC_nt.fasta"
     output:
         "tmp/
         




#Translate nucleotide NBARC sequeces including extended sequences#
rule translate:
     input:
         "tmp/{sample}_NBARC_nt.fasta"
     output: 
         "tmp/{sample}_NBARC_aa.fasta"
     shell:  
         "translate.py {input} {output}"
#----------------------------------Filt#
#Remove * in stop codon#
rule remove_stop_codon:
     input:
         "tmp/{sample}_NBARC_nt.fasta"
     output:
         "tmp/{sample}_NBARC.faa"
     shell:
         "sed 's/*//g' {input} >{output}"
#Run Interproscan, database options: Pfam, coils, gene3D #
rule Interproscan:
     input:
         "tmp/{sample}_NBARC.faa"
     shell:
          "./interproscan.sh -t p -appl Pfam,COILS,Gene3D -i {input} -f tsv,gff3 -d /tmp/"
#extend 20kb upstream and downstream#
#Predict genes by using braker, remove special header first##Remember to use extended one#
rule braker:
     input:
         raw="tmp/{sample}_NBARC_20kb.fasta",
         genome="/genome/{sample}.fa"
     output:
         removed="tmp/{sample}_NBARC_20kb_removed.fasta",
     run:
         shell("sed 's/(//;s/)//' {input.raw} > {output.removed} ")
         shell("braker.pl --cores=6 --genome={input.genome} --prot_seq=$wdir/1.data/RefPlantNLR_aa.fa --ALIGNMENT_TOOL_PATH=/usr/local/genemark-es/4.59/ProtHint/bin/ --prg=ph --epmode --species={sample}
/usr/local/augustus/3.3.3/scripts/gtf2gff.pl <$wdir/0.scripts.logs/braker/augustus.hints.gtf --printExon --out=$wdir/0.scripts.logs/braker/augustus.gff3 --gff3")

