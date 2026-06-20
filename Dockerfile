## Build Stage ##

FROM ubuntu:24.04 AS build

# build tools
#
RUN apt-get update && \
    apt-get install -y \
    build-essential \
    autoconf \
    automake \
    cmake \
    pkg-config \
    libtool \
    git \
    man-db \
    wget \
    vim-tiny &&\
    apt-get clean && rm -rf /var/lib/apt/lists/*

    
# libs
#
RUN apt-get update && \
    apt-get install -y \
    libgsl-dev \
    libboost-all-dev \
    libsuitesparse-dev \
    liblpsolve55-dev \
    libsqlite3-dev \
    libmysql++-dev \
    libboost-iostreams-dev \
    zlib1g-dev \
    libbamtools-dev \
    samtools \
    libjsoncpp-dev \
    libargtable2-dev \
    libhts-dev && \
    apt-get clean && rm -rf /var/lib/apt/lists/*


# start builds
#

# seqstats
#
RUN cd /tmp && \ 
    git clone --recursive https://github.com/clwgg/seqstats && \
    cd seqstats && \
    make && \
    mkdir -p /tmp/staging/usr/local/bin && \
    mv seqstats /tmp/staging/usr/local/bin/

# cdbfasta
# (apt has version 1.0)
# 
RUN cd /opt && \ 
    git clone https://github.com/gpertea/cdbfasta.git && \
    cd cdbfasta && \	
    make 


# Note make install will by default install in /opt
RUN cd /tmp && \
    git clone https://github.com/Gaius-Augustus/Augustus.git && \
    cd Augustus && \
    make && \
	mkdir -p /tmp/staging/usr/local/bin && \
    make DESTDIR=/tmp/staging install

RUN cd /opt && \
    git clone https://github.com/Gaius-Augustus/BRAKER.git && \
    cd BRAKER && \
    cd example && \
    wget http://bioinf.uni-greifswald.de/augustus/datasets/RNAseq.bam

# include ETP
#
RUN cd /opt && \
    git clone https://github.com/KatharinaHoff/GeneMark-ETP.git && \
    mv GeneMark-ETP ETP && \
    chmod a+x /opt/ETP/bin/*py /opt/ETP/bin/*pl /opt/ETP/tools/*

# bamtools
# (apt has version 2.5.2)
#
RUN cd /tmp && \
    git clone https://github.com/pezmaster31/bamtools.git && \
    cd bamtools && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make && \
    make DESTDIR=/tmp/staging install

# diamond
# (apt has diamond-aligner version 2.1.9-1)
#
RUN cd /tmp && \
    wget https://github.com/bbuchfink/diamond/releases/download/v2.2.0/diamond-linux64.tar.gz && \
    tar xzf diamond-linux64.tar.gz && \
    mkdir -p /tmp/staging/usr/local/bin && \
    mv diamond /tmp/staging/usr/local/bin/

# stringtie
# (apt has version 2.2.1)
#
RUN cd /opt && \
    wget https://github.com/gpertea/stringtie/releases/download/v3.0.3/stringtie-3.0.3.Linux_x86_64.tar.gz && \
    tar xzf stringtie-3.0.3.Linux_x86_64.tar.gz && \
    rm -f stringtie-3.0.3.Linux_x86_64.tar.gz && \
    mv stringtie-3.0.3.Linux_x86_64 stringtie-3.0.3

# bedtools
# (apt has version 2.31.1)
#
RUN cd /opt && \
    wget https://github.com/arq5x/bedtools2/releases/download/v2.31.1/bedtools-2.31.1.tar.gz && \
    tar xzf bedtools-2.31.1.tar.gz && \
    rm -f bedtools-2.31.1.tar.gz

# gffread
# (apt has version 0.12.7)
#
RUN cd /opt && \
    wget https://github.com/gpertea/gffread/releases/download/v0.12.9/gffread-0.12.9.Linux_x86_64.tar.gz && \
    tar xzf gffread-0.12.9.Linux_x86_64.tar.gz && \
    rm -f gffread-0.12.9.Linux_x86_64.tar.gz && \
    mv gffread-0.12.9.Linux_x86_64 gffread-0.12.9

# seqkit
# (apt has version 2.3.1)
#
RUN cd /tmp && \
    wget https://github.com/shenwei356/seqkit/releases/download/v2.13.0/seqkit_linux_amd64.tar.gz && \
    tar xzf seqkit_linux_amd64.tar.gz && \
    mkdir -p /tmp/staging/usr/local/bin && \
    mv seqkit /tmp/staging/usr/local/bin/

# clustalo
# (apt has version 1.2.4)
#
RUN cd /tmp && \
    wget https://github.com/GSLBiotech/clustal-omega/archive/refs/tags/1.2.4-cmake.tar.gz && \
    tar xzf 1.2.4-cmake.tar.gz && \
    cd clustal-omega-1.2.4-cmake/ && \
    ./configure && \
    make && \
    make DESTDIR=/tmp/staging install
    
# meme
#
RUN cd /tmp && \
    wget https://meme-suite.org/meme/meme-software/5.5.9/meme-5.5.9.tar.gz && \
    tar -xvzf meme-5.5.9.tar.gz && \
    rm -f meme-5.5.9.tar.gz && \
    cd /tmp/meme-5.5.9 && \
    ./configure --prefix=/opt/meme-5.5.9 --enable-build-libxml2 --enable-build-libxslt && \
    make && \
    make install


# ProtHint
#
RUN cd /opt && \
    wget https://github.com/gatech-genemark/ProtHint/releases/download/v2.6.0/ProtHint-2.6.0.tar.gz && \
    tar -xvzf ProtHint-2.6.0.tar.gz && \
    rm -f ProtHint-2.6.0.tar.gz

    
# blast
# (apt has version 2.12)
#
RUN cd /opt && \
    wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.17.0+-x64-linux.tar.gz && \
    tar -xvzf ncbi-blast-2.17.0+-x64-linux.tar.gz && \
    rm -f ncbi-blast-2.17.0+-x64-linux.tar.gz


# hmmer
# (apt has version 3.4)
#
RUN cd /tmp && \
    wget http://eddylab.org/software/hmmer/hmmer-3.4.tar.gz && \
    tar -xvzf hmmer-3.4.tar.gz && \
    cd hmmer-3.4 && \
    ./configure --prefix=/opt/hmmer-3.4 && \
    make && \
    make install
 
 
# NLR-Annotator
#
RUN cd /opt && \
	git clone -b nlr_parser3 https://github.com/steuernb/NLR-Annotator.git && \
    cd NLR-Annotator && \
    wget https://github.com/steuernb/NLR-Annotator/releases/download/v0.7-beta/meme.xml
    
# tsebra
#
RUN cd /opt && \
    git clone https://github.com/Gaius-Augustus/TSEBRA


# seqtk
# (apt has version 1.4-2)
#
RUN cd /tmp && \
    wget https://github.com/lh3/seqtk/archive/refs/tags/v1.5.tar.gz && \
    tar xzf v1.5.tar.gz && \
	cd seqtk-1.5 && \
	make && \
    mkdir -p /tmp/staging/usr/local/bin && \
    mv seqtk /tmp/staging/usr/local/bin/


# interproscan
#
RUN cd /opt && \
    wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.77-108.0/interproscan-5.77-108.0-64-bit.tar.gz && \
    tar xzf interproscan-5.77-108.0-64-bit.tar.gz && \
    rm -f interproscan-5.77-108.0-64-bit.tar.gz && \
    cd interproscan-5.77-108.0/data && \
    rm -rf antifam cdd funfam hamap ncbifam panther phobius pirsf pirsr prints prosite sfld smart superfamily tmhmm

# FindPlantNLRs
#
RUN cd /opt && \
    wget https://github.com/burntbridge/FindPlantNLRs/archive/refs/tags/v2.1.tar.gz && \
    tar xzf v2.1.tar.gz && \
    rm -f v2.1.tar.gz && \
    mv FindPlantNLRs-2.1 FindPlantNLRs && \
    chmod +x /opt/FindPlantNLRs/scripts/*

# Get the reference fasta
#
RUN cd /opt/FindPlantNLRs/ref_db && \
    wget https://doi.org/10.1371/journal.pbio.3001124.s013 && \
    mv journal.pbio.3001124.s013 ref.fasta
    
# copy stuff to /work

## Final Stage ##
FROM ubuntu:24.04 AS final

COPY --from=build /opt /opt
COPY --from=build /tmp/staging/usr /usr

ENV PATH=${PATH}:/opt/cdbfasta:/opt/hmmer-3.4/bin:/opt/ncbi-blast-2.17.0+/bin/:/opt/TSEBRA/bin
ENV PATH=${PATH}:/opt/BRAKER/scripts:/opt/augustus-3.4.0/bin:/opt/augustus-3.4.0/scripts
ENV PATH=${PATH}:/opt/ETP/bin:/opt/ETP/tools:/opt/ETP/bin/gmes/ProtHint/bin:/opt/ETP/bin/gmes
ENV PATH=${PATH}:/opt/FindPlantNLRs/scripts

ENV AUGUSTUS_CONFIG_PATH=/opt/augustus-3.4.0/config/
ENV AUGUSTUS_BIN_PATH=/opt/augustus-3.4.0/bin/
ENV GENEMARK_PATH=/opt/ETP/bin/
ENV BAMTOOLS_PATH=/usr/local/bin/
ENV DIAMOND_PATH=/usr/local/bin/
ENV BLAST_PATH=/opt/ncbi-blast-2.17.0+/
ENV SAMTOOLS_PATH=/usr/bin/
ENV CDBTOOLS_PATH=/opt/cdbfasta/

RUN apt-get update && \
    apt-get install -y \
    pkg-config \
    vim-tiny && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# augustus
#
RUN apt-get update && \
    apt-get install -y \
	libbamtools2.5.2 \
	libboost-iostreams1.83.0 \
	libboost-serialization1.83.0 \
	libc6 \
	libcolamd3 \
	libgcc-s1 \
	libgsl27 \
	libhts3t64 \
	libmysql++3t64 \
	libmysqlclient21 \
	libsqlite3-0 \
	libstdc++6 \
	libargtable2-0 \
	samtools \
	libdbd-mysql-perl && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

    

# python java
RUN apt-get update && \
    apt-get install -y \
    default-jre \
    python3 \
	python3-biopython \
    python3-pip && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# --break-system-packages is to get around PEP 668
# other ways are to use pipx or python environments
#
RUN pip install snakemake --break-system-packages
RUN pip install gff3tool --break-system-packages

# perl modules
RUN apt-get update && \
    apt-get install -y \
    libyaml-perl \
    libhash-merge-perl \
    libparallel-forkmanager-perl \
    libscalar-util-numeric-perl \
    libclass-data-inheritable-perl \
    libexception-class-perl \
    libtest-pod-perl \
    libfile-which-perl \
    libmce-perl \
    libthread-queue-perl \
    libmath-utils-perl \
    libscalar-list-utils-perl && \
    apt-get clean && rm -rf /var/lib/apt/lists/*


WORKDIR /work