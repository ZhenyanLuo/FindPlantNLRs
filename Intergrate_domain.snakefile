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
rule K-parse:
     input:
         "tmp/{sample}.pfamscan"
     output:
         "tmp/{sample}.protein.fa_pfamscan_parsed.verbose"
     shell:
         "perl ~/analysis/nlr_annotation_scripts/K-parse_Pfam_domains_v3.1.pl --pfam {input} --evalue 0.001 --output {output} --verbose T"
    
