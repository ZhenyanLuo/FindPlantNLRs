#Adapted from https://github.com/ttolessa/NLR-IDs_filter/blob/main/NLR_IDs_filter.sh


import os
genome="genome/"
configfile: "Intergrate_domain.config"
SAMPLES = list(set([x.split(".")[0] for x in os.listdir(genome) if x.endswith("fa")]))
path = 'tmp'
if not os.path.exists(path):
       os.mkdir(path)
result = 'result'
if not os.path.exists(result):
       os.mkdir(result)
pfam_parser = 'pfam_parser'
if not os.path.exists(pfam_parser):
       os.mkdir(pfam_parser)
rule all:
   input:
         expand('tmp/{sample}.pfamscan', sample=SAMPLES)
#Run pfam scan#
rule pfam_scan:
     input:
         "genome/{sample}.fa"
     output:
         "tmp/{sample}.pfamscan"
     params:
         "tmp"
     shell:
         "pfam_scan.pl -fasta {input} -dir {params} -o {output}"         
#Run K-parse_Pfam_domains_v3.1.pl
rule K-parse_pfamscan_parsed:
     input:
         "tmp/{sample}.pfamscan"
     output:
         "tmp/{sample}.protein.fa_pfamscan_parsed.verbose"
     params:
         "pfam_parser/{sample}"
     shell:
         "perl {config[plant_rgenes]}/processing_scripts/K-parse_Pfam_domains_v3.1.pl --pfam {input} --evalue 0.001 --output {output} --verbose T"
         "mkdir -p {params}"
         "mv {output} {params}/"
#Run K-parse_Pfam_domains_NLR-fusions-v2.4.2.pl
rule K-parse_pfamscan_NLR-fusions:
     params:
         dir="pfam_parser/{sample}",
         db="db_descriptions.txt",
         NLR_ID="NLR_IDs/{sample}_NLR_IDs"
     output:
         "tmp/K-parser_pfamscan_NLR-fusions_{sample}.done"
     shell:
         "perl {config[plant_rgenes]}/processing_scripts/K-parse_Pfam_domains_NLR-fusions-v2.4.2.pl --indir {params.dir} --evalue 0.001 -o {params.dir} -d {params.db}"
         "mkdir -p {parmas.NLR_ID}"
#Filter integrated domain
rule 
         
