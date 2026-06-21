## This is a revised version of FindPlantNLRs 

## About FindPlantNLRs
We developed a comprehensive pipeline for annotating predicted NLR genes from a non-masked genome fasta file input. We identify loci using NLR-annotator software (Steuernagel _et al_. 2020), tblastn (Altschul _et al_. 1990) and Hidden Markov Model (HMM) (Eddy 2010). The unmasked loci identified through these methods, and including 20 kb flanking regions, are then annotated with Braker2 software (Hoff _et al_. 2019) using experimentally validated resistance genes as reference (Kourelis _et al_. 2021). Annotated amino acid fasta files are screened for domains using Interproscan (Jones _et al_. 2014) and the predicted coding and amino acid sequences containing both NB-ARC and LRR domains are located back to scaffolds/chromosomes and extracted in fasta and gff3 format.

## Overview of the whole pipeline
![NLR_flowchart2](https://user-images.githubusercontent.com/53864342/232355873-299f26a7-4776-442e-842e-13737dd605d8.jpg)







## Changes in this version
This new docker version simplifies the installation and running. The outputs are not changed.

## To run the FindPlantNLRs pipeline

### Step 1: 
Download the dockerfile and build the docker image.
```
wget https://github.com/peritob/FindPlantNLRs/raw/refs/heads/main/Dockerfile
sudo docker build --tag fpn .
```


### Step 2:
Create a working directory e.g.
```
mkdir FindPlantNLRs
cd FindPlantNLRs
mkdir genome
```
### Step 3:
Move your genome file/s into "genome" directory. Genomes must have extension .fasta. Genomes should be unmasked (uppercase) and have simple headers (short, unique and only have numerics and characters).
For general use with NCBI genomes.
```
awk '/^>/ {print} /^[^>]/ {print toupper($0)}' input.fasta > output_upper.fasta
```
Fix/simplify headers
```
awk '/^>/ {print $1; next} {print}' input.fasta > output.fasta
```

### Step 4:
From the working directory you created (FindPlantNLRs), run the docker image.
```
sudo docker run --volume $(pwd):/work --interactive --tty --rm fpn bash
```
### Step 5:
Now you should be within the docker container. Start the pipeline with this command.

```
FindPlantNLRs.sh
```

### Other information
The pipeline will run through several stages including identifying NLR regions within the genome, annotating the regions, functional predictions followed by a final step of sorting by domain classes. Some steps take a while to run such as tblastn, braker and interproscan. 

### Outputs
All results will be in the result directory and intermediate files will be in tmp. Results include lists of genes by NLR class, integrated domain containing NLRs, protein and coding sequence fasta files by class. Once you have completed your run, it is safe to remove the tmp directory.

## Attributions
Please cite all dependencies used in this pipeline.

- NLR_annotator: https://doi.org/10.1104/pp.19.01273
- Interproscan: https://academic.oup.com/bioinformatics/article/30/9/1236/237988
- hmmsearch: http://eddylab.org/software/hmmer3/3.1b2/Userguide.pdf
- Meme: https://academic.oup.com/nar/article/43/W1/W39/2467905
- BRAKER3: https://genome.cshlp.org/content/34/5/769
- BLAST: https://www.sciencedirect.com/science/article/abs/pii/S0022283605803602
- Reference database (RefNLR) for tblastn and braker is the supplementary file S1: 
https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3001124

## Contributors
Tamene Tolessa, Peri Tobias, Benjamin Schwessinger, Zhenyan Luo

## How to cite this pipeline
If you find this pipeline useful for your work, please cite:  https://doi.org/10.1093/gigascience/giad102
