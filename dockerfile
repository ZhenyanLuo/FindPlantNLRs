From yanolo/findplantnlrs:latest

RUN source /opt/conda/bin/activate FindPlantNLRs
RUN conda create -n Annotate_NLR -y
RUN cd /home && \
    git clone https://github.com/ZhenyanLuo/FindPlantNLRs.git --branch docker_version
RUN source /opt/conda/bin/deactivate FindPlantNLRs
RUN source /opt/conda/bin/activate Annotate_NLR && mamba install braker3 -y && mamba install openjdk=11.0.13 -y && mamba install gxx_linux-64 -y
RUN apt-get install gawk -y
RUN pip3 install gff3tool
COPY gm_key_64.gz /root/gm_key_64.gz
COPY gmes_linux_64.tar.gz /root/
RUN cd /root/ && tar -xvf gmes_linux_64.tar.gz
ENV PATH="/root/gmes_linux_64:$PATH"
ENV PERL5LIB=/opt/conda/envs/Annotate_NLR/lib/perl5/site_perl
RUN cd /root/ && zcat gm_key_64.gz > .gm_key
RUN cd /home/ && \
    git clone -b nlr_parser3 https://github.com/steuernb/NLR-Annotator && \
    wget https://github.com/steuernb/NLR-Annotator/releases/download/v0.7-beta/meme.xml
RUN cd /root/ && git clone https://github.com/gatech-genemark/ProtHint.git
ENV PROTHINT_PATH=/root/ProtHint/bin
ENV PATH="/opt/conda/envs/Annotate_NLR/bin:$PATH"
RUN cd /home/ & mkdir -p reference
#Change the file name of references as required
COPY journal.pbio.3001124.s013  /home/FindPlantNLRs/ref_db/ref.fasta
WORKDIR /home/FindPlantNLRs
ENTRYPOINT /bin/bash

