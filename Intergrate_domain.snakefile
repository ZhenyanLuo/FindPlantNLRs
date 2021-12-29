#Adapted from https://github.com/ttolessa/NLR-IDs_filter/blob/main/NLR_IDs_filter.sh


import os
genome="genome/"
configfile: "FindPlantNLRs.config"
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
          
     shell:
         "pfam_scan.pl -fasta {input} -dir {params} -o {output}"
