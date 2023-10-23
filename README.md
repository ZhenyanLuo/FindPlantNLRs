## About FindPlantNLRs
We developed a comprehensive pipeline for annotating predicted NLR genes from a non-masked genome fasta file input. We identify loci using NLR-annotator software (Steuernagel _et al_. 2020), tblastn (Altschul _et al_. 1990) and Hidden Markov Model (HMM) (Eddy 2010). The unmasked loci identified through these methods, and including 20 kb flanking regions, are then annotated with Braker2 software (Hoff _et al_. 2019) using experimentally validated resistance genes as reference (Kourelis _et al_. 2021). Annotated amino acid fasta files are screened for domains using Interproscan (Jones _et al_. 2014) and the predicted coding and amino acid sequences containing both NB-ARC and LRR domains are located back to scaffolds/chromosomes and extracted in fasta and gff3 format.

## Overview of the whole pipeline
![NLR_flowchart2](https://user-images.githubusercontent.com/53864342/232355873-299f26a7-4776-442e-842e-13737dd605d8.jpg)







## Dependencies

NLR_annotator: https://doi.org/10.1104/pp.19.01273
Steuernagel, B., Witek, K., Krattinger, S.G., Ramirez-Gonzalez, R.H., Schoonbeek, H.J., Yu, G., Baggs, E., Witek, A.I., Yadav, I., Krasileva, K.V. and Jones, J.D., 2020. The NLR-Annotator tool enables annotation of the intracellular immune receptor repertoire. Plant Physiology, 183(2), pp.468-482.

BRAKER2: https://github.com/Gaius-Augustus/BRAKER

Interproscan: https://interproscan-docs.readthedocs.io/en/latest/UserDocs.html#obtaining-a-copy-of-interproscan
Jones, P., Binns, D., Chang, H.Y., Fraser, M., Li, W., McAnulla, C., McWilliam, H., Maslen, J., Mitchell, A., Nuka, G. and Pesseat, S., 2014. InterProScan 5: genome-scale protein function classification. Bioinformatics, 30(9), pp.1236-1240.

Blast: https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download

hmmsearch: http://hmmer.org/download.html

Meme: https://meme-suite.org/meme/meme-software/4.9.1/readme.html

BRAKER2: https://github.com/Gaius-Augustus/BRAKER

```
dependencies:
name: NLR
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - _libgcc_mutex=0.1=conda_forge
  - _openmp_mutex=4.5=2_gnu
  - argtable2=2.13=h14c3975_1001
  - bedtools=2.30.0=h468198e_3
  - blast=2.5.0=hc0b0e79_3
  - boost=1.80.0=py37h48bf904_1
  - boost-cpp=1.80.0=h75c5d50_0
  - bzip2=1.0.8=h7f98852_4
  - c-ares=1.18.1=h7f98852_0
  - ca-certificates=2022.9.24=ha878542_0
  - clustalo=1.2.4=h87f3376_5
  - htslib=1.11=hd3b49d5_2
  - icu=70.1=h27087fc_0
  - keyutils=1.6.1=h166bdaf_0
  - krb5=1.19.3=h3790be6_0
  - ld_impl_linux-64=2.39=hcc3a1bd_1
  - libblas=3.9.0=16_linux64_openblas
  - libcblas=3.9.0=16_linux64_openblas
  - libcurl=7.86.0=h7bff187_1
  - libdeflate=1.7=h7f98852_5
  - libedit=3.1.20191231=he28a2e2_2
  - libev=4.33=h516909a_1
  - libffi=3.4.2=h7f98852_5
  - libgcc=7.2.0=h69d50b8_2
  - libgcc-ng=12.2.0=h65d4601_19
  - libgfortran-ng=12.2.0=h69a702a_19
  - libgfortran5=12.2.0=h337968e_19
  - libgomp=12.2.0=h65d4601_19
  - liblapack=3.9.0=16_linux64_openblas
  - libnghttp2=1.47.0=hdcd2b5c_1
  - libnsl=2.0.0=h7f98852_0
  - libopenblas=0.3.21=pthreads_h78a6416_3
  - libssh2=1.10.0=haa6b8db_3
  - libstdcxx-ng=12.2.0=h46fd767_19
  - libzlib=1.2.13=h166bdaf_4
  - ncurses=6.2=h58526e2_4
  - numpy=1.21.6=py37h976b520_0
  - openssl=1.1.1s=h166bdaf_0
  - perl=5.22.2.1=0
  - perl-bioperl=1.6.924=2
  - perl-threaded=5.32.1=hdfd78af_1
  - pip=22.3.1=pyhd8ed1ab_0
  - python=3.7.12=hb7a2778_100_cpython
  - python_abi=3.7=3_cp37m
  - readline=8.1=h46c0cb4_0
  - samtools=1.11=h6270b1f_0
  - seqkit=2.3.1=h9ee0642_0
  - seqtk=1.3=h5bf99c6_3
  - setuptools=65.5.1=pyhd8ed1ab_0
  - sqlite=3.37.0=h9cd32fc_0
  - tk=8.6.12=h27826a3_0
  - wheel=0.38.4=pyhd8ed1ab_0
  - xz=5.2.6=h166bdaf_0
  - zlib=1.2.13=h166bdaf_4
  - zstd=1.5.2=h6239696_4
  - pip:
    - appdirs==1.4.4
    - attrs==22.1.0
    - certifi==2022.9.24
    - charset-normalizer==2.1.1
    - configargparse==1.5.3
    - connection-pool==0.0.3
    - datrie==0.8.2
    - docutils==0.19
    - dpath==2.1.1
    - fastjsonschema==2.16.2
    - gff3tool==2.1.0
    - gitdb==4.0.10
    - gitpython==3.1.29
    - idna==3.4
    - importlib-metadata==5.1.0
    - importlib-resources==5.10.0
    - jinja2==3.1.2
    - jsonschema==4.17.3
    - jupyter-core==4.12.0
    - markupsafe==2.1.1
    - nbformat==5.7.0
    - pkgutil-resolve-name==1.3.10
    - plac==1.3.5
    - psutil==5.9.4
    - pulp==2.7.0
    - pyrsistent==0.19.2
    - pyyaml==6.0
    - requests==2.28.1
    - reretry==0.11.1
    - smart-open==6.2.0
    - smmap==5.0.0
    - snakemake==7.18.2
    - stopit==1.1.2
    - tabulate==0.9.0
    - throttler==1.2.2
    - toposort==1.7
    - traitlets==5.6.0
    - typing-extensions==4.4.0
    - urllib3==1.26.13
    - wrapt==1.14.1
    - yte==1.5.1
    - zipp==3.11.0
```

## Install with docker
See (https://github.com/ZhenyanLuo/FindPlantNLRs/tree/docker_version)
Please have gm_key_64.gz, gmes_linux_64.tar.gz, journal.pbio.3001124.s013 being prepared in the same directory before running the following command

Build image from yanolo/findplantsnlr
```
docker build -t findplantnlrs -f dockerfile
```
Run container from the image, replace $PWD with the path of your input genome to mount your input folder to /home/FindPlantNLRs/genome
```
docker run -v $PWD:/home/FindPlantNLRs/genome -v ${interproscan}:/home/interproscan -fi findplantnlrs bash 
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

wget https://github.com/steuernb/NLR-Annotator/releases/download/v0.7-beta/meme.xml
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
Click **Send to**, then choose **complete record**-> **file**-> **fasta format**-> **create file** save as E_grandis_chr3.fasta, put this file into genome folder as test data
```
cat E_grandis_chr3.fasta|cut -d' ' -f1 >E_grandis_chr3.fa
```
It's header should look like this

>NC_052614.1

For your own data, make sure sequence headers are short, unique and only have numerics and characters.

### Step 5: Test Snakemake pipeline by using testing file in genome/ folder
```
conda activate NLR
snakemake -s FindPlantNLRs --cores 16 --wait-for-files
snakemake -s Annotate_NLR --cores 1 --wait-for-files
```
### Step 6: Check your main output

### tmp:

#### Three bed files cover identified loci from tblastn, hmm search and NLR-Annotator:

tmp/E_grandis_chr3.blast.20kbflanking.bed

tmp/E_grandis_chr3_NBARC.20kbflanking.bed

tmp/E_grandis_chr3_parser.20kbflanking.bed

#### A merged bed file and corresponding fasta file: 

tmp/E_grandis_chr3.all_20kbflanking_merged.fasta

tmp/E_grandis_chr3.all_20kbflanking_merged_upper.fasta

### Result:

#### Gene prediction result from Augustus:

result/E_grandis_chr3_augustus_aa.fasta

result/E_grandis_chr3_augustus_aa.fasta.tsv

result/E_grandis_chr3_augustus.gff3

#### Lists of different classes of R genes and corresponding bed files:

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

#### A file indicates the completion of classification:

result/E_grandis_chr3_NLRclassification.done

#### Now is ready to use this pipeline on your own samples.


## Contributor
Tamene Tolessa, Peri Tobias, Benjamin Schwessinger, Zhenyan Luo

## How to cite this pipeline

Please also remember to cite all dependencies used in this pipeline
