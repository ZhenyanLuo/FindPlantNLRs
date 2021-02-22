conda install -y seqkit
conda install -y bedtools
conda install -y samtools
conda install -c bioconda -y blast
conda install -c bioconda -y pfam_scan
conda install -c bioconda -y augustus
conda install -c anaconda -y perl
#download interoproscan#
wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.48-83.0/interproscan-5.48-83.0-64-bit.tar.gz
gunzip interproscan-5.48-83.0-64-bit.tar.gz
git clone https://github.com/Gaius-Augustus/BRAKER.git
git clone https://github.com/steuernb/NLR-Annotator.git
git clone https://github.com/peritob/Myrtaceae_NLR_workflow.git
conda install openmpi-mpicc -y
#Install meme which is required for NLR-annotator#
wget https://meme-suite.org/meme/meme-software/5.3.3/meme-5.3.3.tar.gz
gunzip meme-5.3.3.tar.gz
tar -xvf meme-5.3.3.tar
cd meme-5.3.3
make
make test
make install
