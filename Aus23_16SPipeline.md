# Australia NSF Course 2023 16S Data Analysis Pipeline

## Setup of environment
### Creating Working Directory
Variable WD_path is set to the path to our working directory and QIIME2 environment is activated. Make sure to run these 3 lines with every new session
```
WD_path=/anvil/scratch/x-jnash12/Aus23_16S
module load biocontainers
module load qiime2
```

Creates working directory (do this only once - not with every new session)
```
mkdir ${WD_path}
```

## Downloading sklearn taxonomic classifiers on GreenGenes2 and SILVA databses
Oftentimes, QIIME2 users will train their own sklearn taxonomic classifiers based on reference databases trimmed to the gene region they amplified. However, QIIME2 has pre-trained classifiers for the 515F and 806R primer pair (which we used), so we can use those and don't need to train our own (which is memory intensive and can take 24+ hours).

```
cd ${WD_path}
wget https://data.qiime2.org/2023.9/common/silva-138-99-515-806-nb-classifier.qza
wget https://data.qiime2.org/classifiers/greengenes/gg_2022_10_backbone.v4.nb.qza
```