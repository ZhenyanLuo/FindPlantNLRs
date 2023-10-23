## This branch is for docker version of FindPlantNLRs pipeline 

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

BRAKER3: https://github.com/Gaius-Augustus/BRAKER

## Before running the pipeline

### Step 1: Set up docker (root or rootless)

### Step 2: Get reference database

Recommended reference database for tblastn and braker can be downloaded from supplementary file S1 of 'RefPlantNLR is a comprehensive collection of experimentally validated plant disease resistance proteins from the NLR family'
https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3001124

### Step 3: Get a few dependencies being prepared
GeneMark: http://exon.gatech.edu/GeneMark/
Please have gm_key_64.gz, gmes_linux_64.tar.gz, journal.pbio.3001124.s013 prepared in the same folder
Since interproscan database is large, we also recommend you to have interproscan downloaded on host
https://www.ebi.ac.uk/interpro/download/InterProScan/
### Step 4: Create the docker image 
Download the dockerfile from this branch, run
```
docker build -t findplantsnlr:latest -f dockerfile .
```
Once the docker image has been created suscessfully, you can run the following command to bind your interproscan:
```
docker run -v ${your interproscan}:/home/interproscan -v ${your input data file}:/home/FindPlantNLRs/genome -v ${result folder on your host}:/home/FindPlantNLRs/result -ti findplantsnlr bash
```
For your own data, make sure sequence headers are short, unique and only have numerics and characters.
### Step 5: Edit the FindPlantNLRs_config.yaml configure file to change threads used in tblastn 
## You don't need to change any of the path in most of the case
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
### Step 6: Launch the pipeline
For extract sequences potentially contain NLR genes
```
conda activate FindPlantNLRs
cd /home/FindPlantNLRs/
snakemake -s FindPlantNLRs -c {threads as you like} --wait-for-files
```
For annotate NLRs
```
conda activate Annotate_NLR
snakemake -s Annotate_NLR --cores 1 --wait-for-files
```



## Contributor
Tamene Tolessa, Peri Tobias, Benjamin Schwessinger, Zhenyan Luo

## How to cite this pipeline

Please also remember to cite all dependencies used in this pipeline
