conda create -n NLR -y
conda install -c bioconda snakemake=6.10.0 -y
conda install python=3.6.5 -y
conda install -c bioconda bedtools=2.3.0 -y
conda install -c bioconda samtools=1.9 -y
conda install -c bioconda clustalo=1.2.4 -y
conda install -c bioconda hmmer=3.3.2 -y
conda install -c bioconda blast=2.7.1 -y
conda install -c bioconda seqkit=2.0.0 -y
#download interoproscan#
wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.48-83.0/interproscan-5.48-83.0-64-bit.tar.gz
gunzip interproscan-5.48-83.0-64-bit.tar.gz
git clone https://github.com/Gaius-Augustus/BRAKER.git
git clone https://github.com/steuernb/NLR-Annotator.git
