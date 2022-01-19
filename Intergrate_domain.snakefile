#Adapted from https://github.com/ttolessa/NLR-IDs_filter/blob/main/NLR_IDs_filter.sh
#Before running this script, you should already have braker output (augustus.hints.aa)

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
         expand('tmp/{sample}.pfamscan', sample=SAMPLE)
   
rule hmmsearch:
   input:
        "result/{sample}_augustus_aa.fasta"
   output:
        noali="tmp/{sample}.NB-ARC_tblout_noali.txt",
        tblout="tmp/{sample}.NB-ARC_hmmsearch_tblout.perseqhit.txt",
        seqid="tmp/{sample}.NB-ARC_hmmsearch_perseqhit_seqID.txt",
        fasta="tmp/{sample}.NB-ARC_hmmsearch_perseqhit_proteinseq.fa"
   params:
        "{config[hmmfile]}"
   shell:
        """
        hmmsearch -o {output.noali} --tblout {output.tblout} --noali --notextw {params} {input}
        grep -v '#' {output.tblout} {output.seqid}| sort -k5,5n | cut -d ' ' -f1 |cut -d ':' -f2| uniq > {output.seqid}
        seqtk subseq {input} {output.seqid}|sed 's/*//g' > {output.fasta}
        """
#Run pfam scan#
rule pfam_scan:
     input:
         "tmp/{sample}.NB-ARC_hmmsearch_perseqhit_proteinseq.fa"
     output:
         "tmp/{sample}.pfamscan"
     params:
         "tmp"
     shell:
         "pfam_scan.pl -fasta {input} -dir {params} -as -cpu 16 -outfile {output}"         
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
rule filter:
     input:
         pfam="tmp/{sample}.pfamscan",
         fasta="
     output:
         pfamscan="pfam_parser/{sample}/{sample}_augustusprotein_pfamscan.txt",
         nlrid="pfam_parser/{sample}/{sample}_protein_nlrid_pfamscan.txt",
         nlrid_u="pfam_parser/{sample}/{sample}_protein_nlrid_geneid_uniq.txt",
     params:
         "pfam_parser/{sample}"
     shell:
         """
         grep -v '#' {input.pfam} | awk '/^$/ && !f{{f=1;next}}1' > {output.pfamscan}
#test this line seperately#
         awk 'NR==FNR{{for (i=1;i<=NF;i++) a[$i];next}} FNR==1 || ($7 in a)' {params}/*_wordcloud_*.txt {output.pfamscan} | sort | uniq > {output.nlrid}
         cat {output.nlr_id} | cut -d ' ' -f1 | sort | uniq > {output.nlrid_u}
         seqtk subseq ../${genomebase}.NB-ARC_hmmsearch_perseqhit_proteinseq.fa {output.nlrid_u} > ../../NLR_ID_proteins/${genomebase}_NLR_uniqid_protein.fa
         awk '{{print $1 " "$6 " " $7}}' ./${genomebase}_NLR_IDs/${genomebase}_protein_nlrid_pfamscan.txt > ./${genomebase}_NLR_IDs/${genomebase}_protein_nlrid_pfamid_multoccurrence.txt
         cat ${genomebase}.NBLRR.aa | grep '^>' | sed 's/>Species basename//g' > ./${genomebase}_NLR_IDs/${genomebase}_NBLRR_geneid.txt
         cat ${genomebase}.NBTIRs.aa | grep '^>' | sed 's/>Species basename//g' > ./${genomebase}_NLR_IDs/${genomebase}_NBTIR_geneid.txt
         cat ${genomebase}.NBCoils.aa | grep '^>' | sed 's/>Species basename//g' > ./${genomebase}_NLR_IDs/${genomebase}_NBCoil_geneid.txt
         cat ${genomebase}.NBRNLs.aa | grep '^>' | sed 's/>Species basename//g' > ./${genomebase}_NLR_IDs/${genomebase}_NBRNL_geneid.txt
         awk 'NR==FNR{{for (i=1;i<=NF;i++) a[$i];next}} FNR==1 || ($1 in a)' ${genomebase}_NBLRR_geneid.txt ${genomebase}_augustusprotein_pfamscan.txt | sort | uniq | awk '{print $1, $6, $7}' > ${genomebase}_NBLRRid_geneid_pfid.txt
         awk 'NR==FNR{{for (i=1;i<=NF;i++) a[$i];next}} FNR==1 || ($1 in a)' ${genomebase}_NBTIR_geneid.txt ${genomebase}_augustusprotein_pfamscan.txt | sort | uniq | awk '{print $1, $6, $7}' > ${genomebase}_NBTIRid_geneid_pfid.txt
         awk 'NR==FNR{{for (i=1;i<=NF;i++) a[$i];next}} FNR==1 || ($1 in a)' ${genomebase}_NBCoil_geneid.txt ${genomebase}_augustusprotein_pfamscan.txt | sort | uniq | awk '{print $1, $6, $7}' > ${genomebase}_NBCOILid_geneid_pfid.txt
         awk 'NR==FNR{{for (i=1;i<=NF;i++) a[$i];next}} FNR==1 || ($1 in a)' ${genomebase}_NBRNL_geneid.txt ${genomebase}_augustusprotein_pfamscan.txt | sort | uniq | awk '{print $1, $6, $7}' > ${genomebase}_NBRNLid_geneid_pfid.txt
         awk 'NR==FNR{{arr[$0];next}} $0 in arr' ${genomebase}_NBLRRid_geneid_pfid.txt ${genomebase}_protein_nlrid_pfamid_multoccurrence.txt > ${genomebase}_NBLRR_id_number.txt
         awk 'NR==FNR{{arr[$0];next}} $0 in arr' ${genomebase}_NBTIRid_geneid_pfid.txt ${genomebase}_protein_nlrid_pfamid_multoccurrence.txt > ${genomebase}_NBTIR_id_number.txt
         awk 'NR==FNR{{arr[$0];next}} $0 in arr' ${genomebase}_NBCOILid_geneid_pfid.txt ${genomebase}_protein_nlrid_pfamid_multoccurrence.txt > ${genomebase}_NBCOIL_id_number.txt
         awk 'NR==FNR{{arr[$0];next}} $0 in arr' ${genomebase}_NBRNLid_geneid_pfid.txt ${genomebase}_protein_nlrid_pfamid_multoccurrence.txt > ${genomebase}_NBRNL_id_number.txt
