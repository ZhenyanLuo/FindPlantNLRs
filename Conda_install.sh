conda install -y seqkit
conda install -y bedtools
conda install -y samtools
conda install -c bioconda -y blast
conda install -c bioconda -y pfam_scan
conda install -c bioconda -y augustus
#download interoproscan#
wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.48-83.0/interproscan-5.48-83.0-64-bit.tar.gz
gunzip interproscan-5.48-83.0-64-bit.tar.gz
git clone https://github.com/Gaius-Augustus/BRAKER.git
git clone https://github.com/steuernb/NLR-Annotator.git
git clone https://github.com/peritob/Myrtaceae_NLR_workflow.git
git clone https://github.com/krasileva-group/plant_rgenes.git
