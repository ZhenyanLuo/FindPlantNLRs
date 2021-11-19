## Dependencies

NLR_annotator: https://doi.org/10.1104/pp.19.01273
Steuernagel, B., Witek, K., Krattinger, S.G., Ramirez-Gonzalez, R.H., Schoonbeek, H.J., Yu, G., Baggs, E., Witek, A.I., Yadav, I., Krasileva, K.V. and Jones, J.D., 2020. The NLR-Annotator tool enables annotation of the intracellular immune receptor repertoire. Plant Physiology, 183(2), pp.468-482.

BRAKER: https://github.com/Gaius-Augustus/BRAKER

Interproscan:doi: https://interproscan-docs.readthedocs.io/en/latest/UserDocs.html#obtaining-a-copy-of-interproscan
10.1093/bioinformatics/btu031
Jones, P., Binns, D., Chang, H.Y., Fraser, M., Li, W., McAnulla, C., McWilliam, H., Maslen, J., Mitchell, A., Nuka, G. and Pesseat, S., 2014. InterProScan 5: genome-scale protein function classification. Bioinformatics, 30(9), pp.1236-1240.

Blast: https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download

hmmsearch: http://hmmer.org/download.html



Overview of the whole pipeline:

There should be a flowchart here.



### Dependencies need to be installed before running this pipeline

Some packages can be installed by using the environment.yml

Create a new environment using the environment.yml
```
conda env create --name NLR -f environment.yml 
```
or install/update your environment with the environment.yml file
```
conda env update --prefix ./env --file environment.yml --prune
```

Other dependencies need to be install manually:

NLR_annotator: https://github.com/steuernb/NLR-Annotator
Steuernagel, B., Witek, K., Krattinger, S.G., Ramirez-Gonzalez, R.H., Schoonbeek, H.J., Yu, G., Baggs, E., Witek, A.I., Yadav, I., Krasileva, K.V. and Jones, J.D., 2020. The NLR-Annotator tool enables annotation of the intracellular immune receptor repertoire. Plant Physiology, 183(2), pp.468-482.

BRAKER: https://github.com/Gaius-Augustus/BRAKER

Interproscan:doi: https://www.ebi.ac.uk/interpro/download/
Jones, P., Binns, D., Chang, H.Y., Fraser, M., Li, W., McAnulla, C., McWilliam, H., Maslen, J., Mitchell, A., Nuka, G. and Pesseat, S., 2014. InterProScan 5: genome-scale protein function classification. Bioinformatics, 30(9), pp.1236-1240.










### Test Snakemake pipeline by using testing file in genome/ folder
```
Snakemake -s Get_NLR --core 16
```



