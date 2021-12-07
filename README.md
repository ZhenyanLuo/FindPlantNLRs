## Dependencies

NLR_annotator: https://doi.org/10.1104/pp.19.01273
Steuernagel, B., Witek, K., Krattinger, S.G., Ramirez-Gonzalez, R.H., Schoonbeek, H.J., Yu, G., Baggs, E., Witek, A.I., Yadav, I., Krasileva, K.V. and Jones, J.D., 2020. The NLR-Annotator tool enables annotation of the intracellular immune receptor repertoire. Plant Physiology, 183(2), pp.468-482.

BRAKER: https://github.com/Gaius-Augustus/BRAKER

Interproscan:doi: https://interproscan-docs.readthedocs.io/en/latest/UserDocs.html#obtaining-a-copy-of-interproscan
10.1093/bioinformatics/btu031
Jones, P., Binns, D., Chang, H.Y., Fraser, M., Li, W., McAnulla, C., McWilliam, H., Maslen, J., Mitchell, A., Nuka, G. and Pesseat, S., 2014. InterProScan 5: genome-scale protein function classification. Bioinformatics, 30(9), pp.1236-1240.

Blast: https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download

hmmsearch: http://hmmer.org/download.html



## Overview of the whole pipeline:

There should be a flowchart here.




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
tmp:

tmp/E_grandis_chr3.blast.20kbflanking.bed

tmp/E_grandis_chr3_NBARC.20kbflanking.bed

tmp/E_grandis_chr3_parser.20kbflanking.bed

tmp/E_grandis_chr3.all_20kbflanking_merged.fasta

tmp/E_grandis_chr3.all_20kbflanking_merged_upper.fasta

Result:


Now is ready to use this pipeline on your own samples.



## How to cite this pipeline

Please also remember to cite all dependencies used in this pipeline
