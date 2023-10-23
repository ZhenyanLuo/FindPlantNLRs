From yanolo/findplantnlrs:latest

RUN source /opt/conda/bin/activate FindPlantsNLR
RUN conda create -n AnnotatePlantsNLR -y
RUN cd /home && \
    git clone https://github.com/ZhenyanLuo/FindPlantNLRs.git --branch docker-version
RUN source /opt/conda/bin/deactivate FindPlantsNLR 
RUN source /opt/conda/bin/activate AnnotatePlantsNLR && mamba install braker3 -y && mamba install openjdk=11.0.13 -y && mamba install gxx_linux-64 -y
RUN apt-get install gawk
RUN pip3 install gff3tool
COPY gm_key_64.gz /root/gm_key_64.gz
COPY gmes_linux_64.tar.gz /root/
COPY EG_chr3.fa /home/FindPlantNLRs/genome/
RUN cd /root/ && tar -xvf gmes_linux_64.tar.gz
ENV PATH="/root/gmes_linux_64:$PATH"
ENV PERL5LIB=/opt/conda/envs/AnnotatePlantsNLR/lib/perl5/site_perl
RUN cd /root/ && zcat gm_key_64.gz > .gm_key
RUN cd /home/ && \
    git clone -b nlr_parser3 https://github.com/steuernb/NLR-Annotator && \
    wget https://github.com/steuernb/NLR-Annotator/releases/download/v0.7-beta/meme.xml
RUN cd /root/ && git clone https://github.com/gatech-genemark/ProtHint.git
ENV PROTHINT_PATH=/root/ProtHint/bin
ENV PATH="/opt/conda/envs/AnnotatePlantsNLR/bin:$PATH"
RUN cd /home/ & mkdir -p reference 
#Change the file name of references as required
COPY journal.pbio.3001124.s013  /home/FindPlantNLRs/ref_db/ref.fasta
RUN echo "NLR-Annotator: "/home/NLR-Annotator"" >/home/FindPlantNLRs/FindPlantNLRs.config
RUN echo "ref: "/home/FindPlantNLRs/ref_db/ref.fasta"" >>/home/FindPlantNLRs/FindPlantNLRs.config
RUN echo "meme: "/opt/conda/envs/FindPlantsNLR/bin/mast"" >>/home/FindPlantNLRs/FindPlantNLRs.config
RUN echo "meme_xml: "/home/meme.xml"" >>/home/FindPlantNLRs/FindPlantNLRs.config
RUN echo "braker: "/opt/conda/envs/AnnotatePlantsNLR/bin"" >>/home/FindPlantNLRs/FindPlantNLRs.config
RUN echo "augustus: "/opt/conda/envs/AnnotatePlantsNLR/bin/gtf2gff.pl"" >>/home/FindPlantNLRs/FindPlantNLRs.config
RUN echo "interproscan: "/home/interproscan/interproscan.sh"" >>/home/FindPlantNLRs/FindPlantNLRs.config
WORKDIR /home/FindPlantNLRs
ENTRYPOINT /bin/bash
