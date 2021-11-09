## Dependencies

NLR_annotator: https://doi.org/10.1104/pp.19.01273
Steuernagel, B., Witek, K., Krattinger, S.G., Ramirez-Gonzalez, R.H., Schoonbeek, H.J., Yu, G., Baggs, E., Witek, A.I., Yadav, I., Krasileva, K.V. and Jones, J.D., 2020. The NLR-Annotator tool enables annotation of the intracellular immune receptor repertoire. Plant Physiology, 183(2), pp.468-482.

BRAKER: https://github.com/Gaius-Augustus/BRAKER

Interproscan:doi: 10.1093/bioinformatics/btu031
Jones, P., Binns, D., Chang, H.Y., Fraser, M., Li, W., McAnulla, C., McWilliam, H., Maslen, J., Mitchell, A., Nuka, G. and Pesseat, S., 2014. InterProScan 5: genome-scale protein function classification. Bioinformatics, 30(9), pp.1236-1240.

Blast:

hmmsearch:








## Manual

### Create a conda environment and install neccessary dependencies
``` conda create -n NLR -y
conda activate NLR
conda install -c bioconda snakemake=6.10.0 -y
conda install python=3.6.5 -y
conda install -c bioconda bedtools=2.3.0 -y
conda install -c bioconda samtools=1.9 -y
conda install -c bioconda clustalo=1.2.4 -y
conda install -c bioconda hmmer=3.3.2 -y
conda install -c bioconda blast=2.7.1 -y
conda install -c bioconda seqkit=2.0.0 -y
conda install -c conda-forge openjdk -y
```

[Install NLR-Annotator and dependencies](https://github.com/steuernb/NLR-Annotator)

[Install Braker2 for gene prediction](https://github.com/Gaius-Augustus/BRAKER)

[Install Interproscan.sh](https://www.ebi.ac.uk/interpro/download/)
