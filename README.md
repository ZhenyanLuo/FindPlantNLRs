## About FindPlantNLRs
We developed a comprehensive pipeline for annotating predicted NLR genes from a non-masked genome fasta file input. We combine loci identified using (1) NLR-annotator software (Steuernagel et al. 2020) with (2) a basic local alignment search tool (tblastn) (Altschul et al. 1990) using recently compiled and functionally validated NLR amino acid sequences and (3) a nucleotide iterative Hidden Markov Model (HMM) (Eddy 2010) to locate NBARC domains in genomes (Thrimawithana et al. 2019; Christie et al. 2016). While the pipeline was developed to seek NLR genes within Myrtaceae genomes, the supplied NBARC HMMs are suitable for any plant genome search due to the iterative step that builds a         unique species-specific HMM combined with the use of two other steps that incorporate broader models. The unmasked loci identified through these methods, and including 20 kb flanking regions, are then annotated with Braker2 software (Hoff et al. 2019) using protein hints from experimentally validated resistance genes (Kourelis et al. 2021). Annotated amino acid fasta files are screened for domains using Interproscan (Jones et al. 2014) and the predicted coding and amino acid sequences containing both NB-ARC and LRR domains are located back to scaffolds/chromosomes and  extracted in fasta and gff3 format. The complete described protocol including software versions, dependencies, HMMs and additional scripts are available here, 


## Overview of the whole pipeline:






## Dependencies

NLR_annotator: https://doi.org/10.1104/pp.19.01273
Steuernagel, B., Witek, K., Krattinger, S.G., Ramirez-Gonzalez, R.H., Schoonbeek, H.J., Yu, G., Baggs, E., Witek, A.I., Yadav, I., Krasileva, K.V. and Jones, J.D., 2020. The NLR-Annotator tool enables annotation of the intracellular immune receptor repertoire. Plant Physiology, 183(2), pp.468-482.

BRAKER: https://github.com/Gaius-Augustus/BRAKER

Interproscan:doi: https://interproscan-docs.readthedocs.io/en/latest/UserDocs.html#obtaining-a-copy-of-interproscan
10.1093/bioinformatics/btu031
Jones, P., Binns, D., Chang, H.Y., Fraser, M., Li, W., McAnulla, C., McWilliam, H., Maslen, J., Mitchell, A., Nuka, G. and Pesseat, S., 2014. InterProScan 5: genome-scale protein function classification. Bioinformatics, 30(9), pp.1236-1240.

Blast: https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download

hmmsearch: http://hmmer.org/download.html

Meme: https://meme-suite.org/meme/meme-software/4.9.1/readme.html

BRAKER2: https://github.com/Gaius-Augustus/BRAKER

```
dependencies:
  - _libgcc_mutex=0.1=conda_forge
  - _openmp_mutex=4.5=1_gnu
  - aioeasywebdav=2.4.0=py36h5fab9bb_1001
  - aiohttp=3.7.4.post0=py36h8f6f2f9_0
  - amply=0.1.4=py_0
  - appdirs=1.4.4=pyh9f0ad1d_0
  - argtable2=2.13=h14c3975_1001
  - asn1crypto=1.4.0=pyh9f0ad1d_0
  - async-timeout=3.0.1=py_1000
  - atk-1.0=2.36.0=h3371d22_4
  - attmap=0.13.2=pyhd8ed1ab_0
  - attrs=21.2.0=pyhd8ed1ab_0
  - backports=1.0=py_2
  - backports.functools_lru_cache=1.6.4=pyhd8ed1ab_0
  - bcrypt=3.2.0=py36h8f6f2f9_1
  - bedtools=2.30.0=h7d7f7ad_2
  - binutils_impl_linux-64=2.36.1=h193b22a_2
  - blast=2.7.1=h4422958_6
  - boost=1.67.0=py36h3e44d54_0
  - boost-cpp=1.67.0=h3a22d5f_0
  - boto3=1.19.12=pyhd8ed1ab_0
  - botocore=1.22.12=pyhd8ed1ab_0
  - brotlipy=0.7.0=py36h8f6f2f9_1001
  - bzip2=1.0.8=h7f98852_4
  - c-ares=1.18.1=h7f98852_0
  - ca-certificates=2021.10.8=ha878542_0
  - cachetools=4.2.4=pyhd8ed1ab_0
  - cairo=1.16.0=h18b612c_1001
  - certifi=2021.5.30=py36h5fab9bb_0
  - cffi=1.14.6=py36hd8eec40_1
  - chardet=4.0.0=py36h5fab9bb_1
  - charset-normalizer=2.0.0=pyhd8ed1ab_0
  - clustalo=1.2.4=h1b792b2_4
  - coincbc=2.10.5=hcee13e7_1
  - configargparse=1.5.3=pyhd8ed1ab_0
  - connection_pool=0.0.3=pyhd3deb0d_0
  - cryptography=2.5=py36hb7f436b_1
  - curl=7.64.0=h646f8bb_0
  - cycler=0.11.0=pyhd8ed1ab_0
  - datrie=0.8.2=py36h8f6f2f9_2
  - decorator=5.1.0=pyhd8ed1ab_0
  - docutils=0.17.1=py36h5fab9bb_0
  - dropbox=11.23.0=pyhd8ed1ab_0
  - expat=2.4.1=h9c3ff4c_0
  - fftw=3.3.10=nompi_h74d3f13_101
  - filechunkio=1.8=py_2
  - filelock=3.3.2=pyhd8ed1ab_0
  - font-ttf-dejavu-sans-mono=2.37=hab24e00_0
  - font-ttf-inconsolata=3.000=h77eed37_0
  - font-ttf-source-code-pro=2.038=h77eed37_0
  - font-ttf-ubuntu=0.83=hab24e00_0
  - fontconfig=2.13.1=hba837de_1005
  - fonts-conda-ecosystem=1=0
  - fonts-conda-forge=1=0
  - freetype=2.10.4=h0708190_1
  - fribidi=1.0.10=h36c2ea0_0
  - ftputil=5.0.1=pyhd8ed1ab_0
  - gcc=11.2.0=h702ea55_2
  - gcc_bootstrap_linux-64=11.2.0=h6cd0796_2
  - gcc_bootstrap_linux-s390x=11.2.0=h15f7010_2
  - gcc_impl_linux-64=11.2.0=h82a94d6_11
  - gdk-pixbuf=2.42.6=h04a7f16_0
  - gettext=0.19.8.1=h73d1719_1008
  - ghostscript=9.18=1
  - giflib=5.2.1=h36c2ea0_2
  - gitdb=4.0.9=pyhd8ed1ab_0
  - gitpython=3.1.18=pyhd8ed1ab_0
  - glib=2.70.0=h780b84a_1
  - glib-tools=2.70.0=h780b84a_1
  - gmp=6.1.2=hf484d3e_1000
  - gnutls=3.5.19=h2a4e5f8_1
  - google-api-core=1.31.2=pyhd8ed1ab_1
  - google-api-python-client=2.28.0=pyhd8ed1ab_0
  - google-auth=1.28.0=pyhd3eb1b0_0
  - google-auth-httplib2=0.1.0=pyhd8ed1ab_0
  - google-cloud-core=1.7.2=pyh6c4a22f_0
  - google-cloud-storage=1.19.0=py_0
  - google-crc32c=1.1.2=py36h0208b43_0
  - google-resumable-media=2.1.0=pyh6c4a22f_0
  - googleapis-common-protos=1.53.0=py36h5fab9bb_0
  - graphite2=1.3.13=h58526e2_1001
  - graphviz=2.42.3=h0511662_0
  - grpcio=1.16.0=py36h4f00d22_1000
  - gts=0.7.6=h64030ff_2
  - harfbuzz=2.4.0=h37c48d4_1
  - hmmer=3.3.2=h1b792b2_1
  - httplib2=0.20.2=pyhd8ed1ab_1
  - icu=58.2=hf484d3e_1000
  - idna=3.1=pyhd3deb0d_0
  - idna_ssl=1.1.0=py36h9f0ad1d_1001
  - imagemagick=7.0.11_12=pl5320h0cb4662_0
  - importlib-metadata=4.8.1=py36h5fab9bb_0
  - importlib_resources=5.4.0=pyhd8ed1ab_0
  - iniconfig=1.1.1=pyh9f0ad1d_0
  - ipython_genutils=0.2.0=py_1
  - jbig=2.1=h7f98852_2003
  - jinja2=3.0.2=pyhd8ed1ab_0
  - jmespath=0.10.0=pyh9f0ad1d_0
  - jpeg=9d=h36c2ea0_0
  - jsonschema=2.6.0=py36_1002
  - jupyter_core=4.8.1=py36h5fab9bb_0
  - kernel-headers_linux-64=2.6.32=he073ed8_15
  - kiwisolver=1.3.1=py36h605e78d_1
  - krb5=1.16.3=hc83ff2d_1000
  - lcms2=2.12=hddcbb42_0
  - ld_impl_linux-64=2.36.1=hea4e1c9_2
  - lerc=3.0=h9c3ff4c_0
  - libblas=3.9.0=12_linux64_openblas
  - libcblas=3.9.0=12_linux64_openblas
  - libcrc32c=1.1.2=h9c3ff4c_0
  - libcurl=7.64.0=h01ee5af_0
  - libdeflate=1.8=h7f98852_0
  - libedit=3.1.20170329=0
  - libffi=3.4.2=h9c3ff4c_4
  - libgcc=7.2.0=h69d50b8_2
  - libgcc-devel_linux-64=11.2.0=h0952999_11
  - libgcc-ng=11.2.0=h1d223b6_11
  - libgfortran-ng=11.2.0=h69a702a_11
  - libgfortran5=11.2.0=h5c6108e_11
  - libglib=2.70.0=h174f98d_1
  - libgomp=11.2.0=h1d223b6_11
  - libiconv=1.16=h516909a_0
  - liblapack=3.9.0=12_linux64_openblas
  - libopenblas=0.3.18=pthreads_h8fe5266_0
  - libpng=1.6.37=h21135ba_2
  - libprotobuf=3.18.0=h780b84a_1
  - librsvg=2.50.3=hfa39831_1
  - libsanitizer=11.2.0=he4da1e4_11
  - libsodium=1.0.18=h36c2ea0_1
  - libssh2=1.8.0=h1ad7b7a_1003
  - libstdcxx-ng=11.2.0=he4da1e4_11
  - libtiff=4.3.0=h6f004c6_2
  - libtool=2.4.6=h9c3ff4c_1008
  - libuuid=2.32.1=h7f98852_1000
  - libwebp=1.2.1=h3452ae3_0
  - libwebp-base=1.2.1=h7f98852_0
  - libxcb=1.13=h7f98852_1003
  - libxml2=2.9.12=h03d6c58_0
  - libzlib=1.2.11=h36c2ea0_1013
  - logmuse=0.2.6=pyh8c360ce_0
  - lz4-c=1.9.3=h9c3ff4c_1
  - markupsafe=2.0.1=py36h8f6f2f9_0
  - matplotlib-base=3.3.4=py36hd391965_0
  - more-itertools=8.10.0=pyhd8ed1ab_0
  - multidict=5.2.0=py36h8f6f2f9_0
  - nbformat=5.1.3=pyhd8ed1ab_0
  - ncurses=5.9=10
  - nettle=3.3=0
  - networkx=2.6.3=pyhd8ed1ab_1
  - numpy=1.19.5=py36hfc0c790_2
  - oauth2client=4.1.3=py_0
  - olefile=0.46=pyh9f0ad1d_1
  - openjpeg=2.4.0=hb52868f_1
  - openssl=1.0.2u=h516909a_0
  - packaging=21.2=pyhd8ed1ab_1
  - pandas=1.1.5=py36h284efc9_0
  - pango=1.42.4=h7062337_4
  - paramiko=2.8.0=pyhd8ed1ab_0
  - pcre=8.45=h9c3ff4c_0
  - peppy=0.31.2=pyhd8ed1ab_0
  - perl=5.32.1=1_h7f98852_perl5
  - perl-archive-tar=2.18=1
  - perl-exporter-tiny=0.042=1
  - perl-list-moreutils=0.413=1
  - perl-threaded=5.26.0=0
  - pillow=8.3.2=py36h676a545_0
  - pip=21.3.1=pyhd8ed1ab_0
  - pixman=0.38.0=h516909a_1003
  - pkg-config=0.29.2=h36c2ea0_1008
  - ply=3.11=py_1
  - prettytable=2.4.0=pyhd8ed1ab_0
  - protobuf=3.18.0=py36hc4f0c31_0
  - psutil=5.8.0=py36h8f6f2f9_1
  - pthread-stubs=0.4=h36c2ea0_1001
  - pulp=2.5.1=py36h5fab9bb_0
  - py=1.11.0=pyh6c4a22f_0
  - pyasn1=0.4.8=py_0
  - pyasn1-modules=0.2.7=py_0
  - pycparser=2.21=pyhd8ed1ab_0
  - pygments=2.10.0=pyhd8ed1ab_0
  - pygraphviz=1.6=py36he9cec78_1
  - pynacl=1.4.0=py36h8f6f2f9_2
  - pyopenssl=19.0.0=py36_0
  - pyparsing=2.4.7=pyhd8ed1ab_1
  - pysftp=0.2.9=py_1
  - pysocks=1.7.1=py36h5fab9bb_3
  - pytest=3.2.5=py36_0
  - python=3.6.5=1
  - python-dateutil=2.8.2=pyhd8ed1ab_0
  - python-irodsclient=1.0.0=pyhd8ed1ab_0
  - python_abi=3.6=2_cp36m
  - pytz=2021.3=pyhd8ed1ab_0
  - pyu2f=0.1.5=pyhd8ed1ab_0
  - pyyaml=5.4.1=py36h8f6f2f9_1
  - ratelimiter=1.2.0=py_1002
  - readline=7.0=0
  - requests=2.26.0=pyhd8ed1ab_0
  - rsa=4.7.2=pyh44b312d_0
  - s3transfer=0.5.0=pyhd8ed1ab_0
  - samtools=1.9=h46bd0b3_0
  - scipy=1.5.3=py36h81d768a_1
  - seqkit=2.0.0=h9ee0642_0
  - setuptools=58.0.4=py36h5fab9bb_2
  - simplejson=3.17.5=py36h8f6f2f9_0
  - six=1.16.0=pyh6c4a22f_0
  - slacker=0.14.0=py_0
  - smart_open=5.2.1=pyhd8ed1ab_0
  - smmap=3.0.5=pyh44b312d_0
  - snakemake=6.10.0=hdfd78af_0
  - snakemake-minimal=6.10.0=pyhdfd78af_0
  - sqlite=3.20.1=2
  - stone=3.2.1=pyhd8ed1ab_0
  - stopit=1.1.2=py_0
  - sysroot_linux-64=2.12=he073ed8_15
  - tabulate=0.8.9=pyhd8ed1ab_0
  - tk=8.6.11=h27826a3_1
  - toml=0.10.2=pyhd8ed1ab_0
  - toposort=1.7=pyhd8ed1ab_0
  - tornado=6.1=py36h8f6f2f9_1
  - traitlets=4.3.3=pyhd8ed1ab_2
  - typing-extensions=3.10.0.2=hd8ed1ab_0
  - typing_extensions=3.10.0.2=pyha770c72_0
  - tzdata=2021e=he74cb21_0
  - ubiquerg=0.6.1=pyh9f0ad1d_0
  - uritemplate=3.0.1=py_0
  - urllib3=1.26.7=pyhd8ed1ab_0
  - veracitools=0.1.3=py_0
  - wcwidth=0.2.5=pyh9f0ad1d_2
  - wheel=0.37.0=pyhd8ed1ab_1
  - wrapt=1.13.1=py36h8f6f2f9_0
  - xmlrunner=1.7.7=py_0
  - xorg-kbproto=1.0.7=h7f98852_1002
  - xorg-libice=1.0.10=h7f98852_0
  - xorg-libsm=1.2.3=hd9c2040_1000
  - xorg-libx11=1.7.2=h7f98852_0
  - xorg-libxau=1.0.9=h7f98852_0
  - xorg-libxdmcp=1.1.3=h7f98852_0
  - xorg-libxext=1.3.4=h7f98852_1
  - xorg-libxpm=3.5.13=h7f98852_0
  - xorg-libxrender=0.9.10=h7f98852_1003
  - xorg-libxt=1.2.1=h7f98852_2
  - xorg-renderproto=0.11.1=h7f98852_1002
  - xorg-xextproto=7.3.0=h7f98852_1002
  - xorg-xproto=7.0.31=h7f98852_1007
  - xz=5.2.5=h516909a_1
  - yaml=0.2.5=h516909a_0
  - yarl=1.6.3=py36h8f6f2f9_2
  - zipp=3.6.0=pyhd8ed1ab_0
  - zlib=1.2.11=h36c2ea0_1013
  - zstd=1.5.0=ha95c52a_0
```




## Before running the pipeline

### Step 1: Install required dependencies before running this pipeline

```
git clone https://github.com/ZhenyanLuo/FindPlantNLRs
cd FindPlantNLRs
```
Some packages can be installed by using the environment.yml

Create a new environment using the environment.yml
```
conda env create --name NLR -f environment.yml 
```
or install/update your environment with the environment.yml file
```
conda activate {your environment}
conda env update --file environment.yml
```
Install open-jdk for NLR-Annotator and Interproscan
```
conda activate {your environment}
wget https://anaconda.org/conda-forge/openjdk/11.0.9.1/download/linux-64/openjdk-11.0.9.1-h5cc2fde_1.tar.bz2
conda install openjdk-11.0.9.1-h5cc2fde_1.tar.bz2
rm openjdk-11.0.9.1-h5cc2fde_1.tar.bz2
```
Other dependencies need to be install manually to the FindPlantNLRs folder:

NLR_annotator: https://github.com/steuernb/NLR-Annotator
**Use nlr_parser3 branch only and don't forget to download meme.xml**
```
git clone -b nlr_parser3 https://github.com/steuernb/NLR-Annotator
```
BRAKER: https://github.com/Gaius-Augustus/BRAKER

Interproscan: https://www.ebi.ac.uk/interpro/download/

**Download and install meme-4.9.1 manually**

meme-4.9.1: https://meme-suite.org/meme/meme-software/4.9.1/readme.html

**Install Perl modules for BRAKER2**
```
cpan install Scalar::Util::Numeric Parallel::ForkManager File::HomeDir List::MoreUtils 
```

### Step 2: Get reference database

Recommended reference database for tblastn and braker can be downloaded from supplementary file S1 of 'RefPlantNLR is a comprehensive collection of experimentally validated plant disease resistance proteins from the NLR family'
https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3001124
Place the reference sequence in ref_db folder and add **absolute path** of the reference in the FindPlantNLRs_config.yaml file

### Step 3: Edit the FindPlantNLRs_config.yaml configure file by adding path to dependencies
```
#Path to NLR-Annotator

NLR-Annotator: "{path/to/your}/NLR-Annotator"

#Path to reference database

ref: "{path/to/your/ref/fasta}"

#Add path to meme

meme: "{path/to/your/meme}/mast"

#Add path to meme.xml from NLR-Annotator

meme_xml: "{path/to/your}/meme.xml"

#Add path to braker.pl

braker: "{path/to/your}/BRAKER"

#Add path to augustus script

augustus: "{path/to/your/augustus}/scripts/gtf2gff.pl

#Add path to Interproscan

interproscan: "{path/to/your/interproscan}/interproscan.sh"
```

### Step 4: Place your samples in the 'genome' folder, make sure their file names only have one dot and ends with 'fa' (e.g. E_globulus.fa) 
Chromosome 3 of E_grandis is a good choice
https://www.ncbi.nlm.nih.gov/nuccore/NC_052614.1
Click **Send to**, then choose **complete record**, **file**, **fasta format**, **create file** save as E_grandis_chr3.fasta, put this file into genome folder as test data
```
cat E_grandis_chr3.fasta|cut -d' ' -f1 >E_grandis_chr3.fa
```
Make sure sequence headers are short, unique and only have numerics and characters.

### Step 5: Test Snakemake pipeline by using testing file in genome/ folder
```
conda activate NLR
snakemake -s FindPlantNLRs --cores 16
```
### Step 6: Check your main output

### tmp:

tmp/E_grandis_chr3.blast.20kbflanking.bed

tmp/E_grandis_chr3_NBARC.20kbflanking.bed

tmp/E_grandis_chr3_parser.20kbflanking.bed

tmp/E_grandis_chr3.all_20kbflanking_merged.fasta

tmp/E_grandis_chr3.all_20kbflanking_merged_upper.fasta

### Result:

result/E_grandis_chr3_augustus_aa.fasta

result/E_grandis_chr3_augustus_aa.fasta.tsv

result/E_grandis_chr3_augustus.gff3

result/E_grandis_chr3_BNL.list

result/E_grandis_chr3_BNL.gff3

result/E_grandis_chr3_CNL.list

result/E_grandis_chr3_CNL.gff3

result/E_grandis_chr3_JNL.list

result/E_grandis_chr3_JNL.gff3

result/E_grandis_chr3_NBARC.list

result/E_grandis_chr3_NBARC.gff3

result/E_grandis_chr3_NLR.list

result/E_grandis_chr3_NLR.gff3

result/E_grandis_chr3_RNL.list

result/E_grandis_chr3_RNL.gff3

result/E_grandis_chr3_RxNL.list

result/E_grandis_chr3_RxNL.gff3

result/E_grandis_chr3_TNL.list

result/E_grandis_chr3_TNL.gff3

result/E_grandis_chr3_NLRclassification.done

Now is ready to use this pipeline on your own samples.


## Contributor
Tamene Tolessa, Peri Tobias, Benjamin Schwessinger, Zhenyan Luo

## How to cite this pipeline

Please also remember to cite all dependencies used in this pipeline
