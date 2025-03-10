# Australia NSF Course 2023 16S Data Analysis Pipeline
This is a set of scripts to analyze 16S metabarcoding data (amplified with modified primer pair 515F and 806R from Cregger et al 2018) using QIIME2. These scripts were written to run on the Purdue ANVIL cluster using QIIME2 with the q2cli version 2023.7.0 that is available as a module on the cluster.  

## Setup of environment and downloads
### Code to run at beginning of every session
Variable WD_path is set to the path to our working directory. Make sure to run this with each new session
```
WD_path=/anvil/scratch/x-jnash12/Aus23_16S
```

### Creates working directory (do this only once - not with every new session)
```
mkdir ${WD_path}
```

### Downloads demultiplexed sequence data from Vilgalys Lab Dropbox into working directory using rclone
```
sbatch -o slurm-%j-DataDL.out --partition=shared --account=BIO230020 --export=ALL -t 2:00:00 -c 1 --wrap="rclone copy remote:'Vilgalys Lab'/NGS_RAW/Nash_8732_23101202 ${WD_path}/Nash_8732_23101202"
```

### Downloads sklearn taxonomic classifiers for GreenGenes2 and SILVA databses
Oftentimes, QIIME2 users will train their own sklearn taxonomic classifiers based on reference databases trimmed to the gene region they amplified. However, QIIME2 has pre-trained classifiers for the 515F and 806R primer pair (which we used), so we can use those and don't need to train our own (which is memory intensive and can take 24+ hours).

```
cd ${WD_path}
wget https://data.qiime2.org/2023.9/common/silva-138-99-515-806-nb-classifier.qza
wget https://data.qiime2.org/classifiers/greengenes/gg_2022_10_backbone.v4.nb.qza
```

### Preparing Manifest File
QIIME2 uses a special data format with the file type qza to store sequence project data. Fastq files must be imported into QIIME2 first. Because our data are already demultiplexed, we need to generate a tsv file that links sample IDs to the filepaths of the forward and reverse reads to be used in the import command.

```
ls ${WD_path}/Nash_8732_23101202 | grep "R1_001.fastq.gz" | sed 's/_R1_001.fastq.gz//g' | grep "16S" > ${WD_path}/filelist_16S
printf "%s\t%s\t%s\n" "sample-id" "forward-absolute-filepath" "reverse-absolute-filepath" > ${WD_path}/QIIMEManifest.tsv
for i in $(cat ${WD_path}/filelist_16S)
do
Sample=$(echo "$i" | awk -F '_' '{print $1}')
SampleType=$(echo "$i" | awk -F '_' '{print $2}')
printf "%s\t%s\t%s\n" "${Sample}_${SampleType}" "${WD_path}/Nash_8732_23101202/${i}_R1_001.fastq.gz" "${WD_path}/Nash_8732_23101202/${i}_R2_001.fastq.gz" >> ${WD_path}/QIIMEManifest.tsv
done
```

# QIIME Time
All following commands will be done within QIIME2. QIIME2 is available from Purdue ANVIL as a module that can be loaded. To do this, we create batch jobs for each command and load the QIIME2 module at the beginning of each script. To run the following chunks of code, you will need to create files using the nano command with the script name as noted in bold at top of code chunk, and then paste each chunk into the file. Then scripts can be run by typing sbatch {scriptname}

## Imports data into QIIME 
Uses manifest to import paired end sequence data into QIIME2 then creates a sequencing summary visualization in qzv format that can be downloaded onto our local machine and viewed at this link https://view.qiime2.org/

Takes a couple minutes

### DataImport
```
#!/bin/bash
#SBATCH -o slurm-%j-QIIME_Import.out
#SBATCH -c 1
#SBATCH --partition=shared 
#SBATCH -A BIO230020
#SBATCH --export=ALL
#SBATCH -t 1:00:00

module load biocontainers
module load qiime2
WD_path=/anvil/scratch/x-jnash12/Aus23_16S

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path ${WD_path}/QIIMEManifest.tsv \
  --output-path ${WD_path}/Aus_16S_demux.qza \
  --input-format PairedEndFastqManifestPhred33V2

qiime demux summarize \
  --i-data ${WD_path}/Aus_16S_demux.qza \
  --o-visualization ${WD_path}/Aus_16S_demux.qzv
 ```

## Cutadapt to trim primers and adapter sequences
Cutadapt can trim primers/adapters from both the 5' and 3' ends of the forward and reverse reads. Primers/adapters may occur on the 3' end in the event of read through (i.e. the 300 bp read extends to the end of the DNA fragment).

We provide reverse complements of the adapters and primers under the p-adapter flags in case read through happened. Normal primer sequences are provided in the p-front argument to trim them from the 5' end

Takes \~30 minutes

### RunCutadapt
```
#!/bin/bash
#SBATCH -o slurm-%j-cutadapt.out
#SBATCH -c 1
#SBATCH --partition=shared 
#SBATCH -A BIO230020
#SBATCH --export=ALL
#SBATCH -t 2:00:00

module load biocontainers
module load qiime2
WD_path=/anvil/scratch/x-jnash12/Aus23_16S

qiime cutadapt trim-paired \
	--i-demultiplexed-sequences ${WD_path}/Aus_16S_demux.qza \
	--p-adapter-f ATTAGAWACCCBDGTAGTCC \
	--p-adapter-f AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
	--p-front-f GTGCCAGCMGCCGCGGTAA \
	--p-front-f GTGCCAGCMGCWGCGGTAA \
	--p-front-f GTGCCAGCMGCCGCGGTCA \
	--p-front-f GTGKCAGCMGCCGCGGTAA \
	--p-front-f GCCTCCCTCGCGCCATCAGAGATGTGTATAAGAGACAG \
	--p-adapter-r TTACCGCGGCKGCTGMCAC \
	--p-adapter-r TGACCGCGGCKGCTGGCAC \
	--p-adapter-r TTACCGCWGCKGCTGGCAC \
	--p-adapter-r TTACCGCGGCKGCTGGCAC \
	--p-adapter-r CTGTCTCTTATACACATCTCTGATGGCGCGAGGGAGGC \
	--p-front-r GGACTACHVGGGTWTCTAAT \
	--p-front-r GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT \
	--p-error-rate 0.15 \
	--output-dir ${WD_path}/Aus_16S_demux_cutadapt \
	--verbose

qiime demux summarize \
  --i-data ${WD_path}/Aus_16S_demux_cutadapt/trimmed_sequences.qza \
  --o-visualization ${WD_path}/Aus_16S_demux_cutadapt.qzv
```

## Dada2 for ASV calling + QC
Dada2 is a pipeline for sequence error correction, chimera detection, read merging, and ASV calling. Output is a "feature table" (analagous to an OTU table) and a set of representative DNA sequences

Dada2 has very stringent quality filtering standards and will discard reads that do not meet. The --p-trunc-len-f and --p-trunc-len-r parameters allow us to trim off low quality regions at the end of the forward and reverse read, respectively, which helps more reads pass QC. The danger is if we trim too much, then the reads will not be long enough to merge. By inspecting the sequence quality plot after cutadapt, we can determine optimal parameters for this step.

The V4 region of 16S that we amplified is \~250bp after trimming primers/adapters. Reverse reads are often lower quality than forward reads and usually require trimming at an earlier position.  Dada2 requires 20-50 bp of overlap between the reads, so this leaves plenty of overlap remaining to recover the full V4 region (and we could even trim more if we are concerned about low quality bases). Here, we tested 4 different primer trimming settings as indicated below to determine which resulted in the highest percentage of reads being retained.

Trimming parameters (modify ForTrim and RevTrim variables in script below to match these):
* For: 180 Rev: 120
* For: 200 Rev: 140
* For: 220 Rev: 160
* For: 240 Rev: 180

Takes \~1 hour

### RunDada2_F180_R120
```
#!/bin/bash
#SBATCH -o slurm-%j-dada2.out
#SBATCH -c 48
#SBATCH --partition=shared 
#SBATCH -A BIO230020
#SBATCH --export=ALL
#SBATCH -t 24:00:00
#SBATCH --mem=100G

module load biocontainers
module load qiime2
ForTrim=180
RevTrim=120
WD_path=/anvil/scratch/x-jnash12/Aus23_16S

time qiime dada2 denoise-paired \
	--i-demultiplexed-seqs ${WD_path}/Aus_16S_demux_cutadapt/trimmed_sequences.qza \
	--p-trunc-len-f ${ForTrim} \
	--p-trunc-len-r ${RevTrim} \
	--p-trim-left-f 0 \
	--p-trim-left-r 0 \
	--p-n-threads 0 \
	--output-dir ${WD_path}/16S_FeatureTable_${ForTrim}_${RevTrim}

qiime metadata tabulate \
	--m-input-file ${WD_path}/16S_FeatureTable_${ForTrim}_${RevTrim}/denoising_stats.qza \
	--o-visualization ${WD_path}/denoising_stats_${ForTrim}_${RevTrim}.qzv

qiime feature-table summarize \
	--i-table ${WD_path}/16S_FeatureTable_${ForTrim}_${RevTrim}/table.qza \
	--o-visualization ${WD_path}/table_summary_${ForTrim}_${RevTrim}.qzv
```

## Taxonomy assignment with Greengenes2 and SILVA
Here we use the sklearn taxonomic classifiers trained on Greengenes2 and SILVA that we previously downloaded to assign taxonomy to the representative sequences that were generated by DADA2. We can compare the taxonomic assignments from SILVA and Greengenes2 to decide which we would like to use

The taxonomic classification with the Silva classifier was returning out of memory errors on the standard 256 gb nodes so I had to run on a special high memory 1 TB node

Takes \~10 minutes

### SilvaClassify
```
#!/bin/bash
#SBATCH -o slurm-%j-SILVA_classify.out
#SBATCH -c 128
#SBATCH --partition=highmem 
#SBATCH -A BIO230020
#SBATCH --export=ALL
#SBATCH -t 24:00:00

module load biocontainers
module load qiime2
WD_path=/anvil/scratch/x-jnash12/Aus23_16S

qiime feature-classifier classify-sklearn \
	--i-reads ${WD_path}/16S_FeatureTable_180_120/representative_sequences.qza \
	--i-classifier ${WD_path}/silva-138-99-515-806-nb-classifier.qza \
	--o-classification ${WD_path}/Silva_16S_taxonomy.qza \
	--p-n-jobs 128
```

### GGClassify
```
#!/bin/bash
#SBATCH -o slurm-%j-GG_classify.out
#SBATCH -c 96
#SBATCH --partition=shared 
#SBATCH -A BIO230020
#SBATCH --export=ALL
#SBATCH -t 24:00:00
#SBATCH --mem=200G

module load biocontainers
module load qiime2
WD_path=/anvil/scratch/x-jnash12/Aus23_16S

qiime feature-classifier classify-sklearn \
	--i-reads ${WD_path}/16S_FeatureTable_180_120/representative_sequences.qza \
	--i-classifier ${WD_path}/gg_2022_10_backbone.v4.nb.qza \
	--o-classification ${WD_path}/GG_16S_taxonomy.qza \
	--p-n-jobs 96
```

## Discard mitochondrial and chloroplast reads
Mitochondrial and chloroplasts ASVs are filtered out based on SILVA taxonomy

### FilterTable
```
#!/bin/bash
#SBATCH -o slurm-%j-FilterTable.out
#SBATCH -c 1
#SBATCH --partition=shared 
#SBATCH -A BIO230020
#SBATCH --export=ALL
#SBATCH -t 2:00:00


module load biocontainers
module load qiime2
WD_path=/anvil/scratch/x-jnash12/Aus23_16S

qiime taxa filter-table \
--i-table ${WD_path}/16S_FeatureTable_180_120/table.qza  \
--i-taxonomy ${WD_path}/Silva_16S_taxonomy.qza \
--p-exclude mitochondria,chloroplast \
--o-filtered-table ${WD_path}/16S_FeatureTable_180_120/table-no-mitochondria-no-chloroplast.qza

qiime feature-table summarize \
	--i-table ${WD_path}/16S_FeatureTable_180_120/table-no-mitochondria-no-chloroplast.qza \
	--o-visualization ${WD_path}/16S_FeatureTable_180_120/table-no-mitochondria-no-chloroplast_summary.qzv
```

### RareCurves100
```
#!/bin/bash
#SBATCH -o slurm-%j-RareCurve100.out
#SBATCH -c 1
#SBATCH --partition=shared 
#SBATCH -A BIO230020
#SBATCH --export=ALL
#SBATCH -t 12:00:00


module load biocontainers
module load qiime2
WD_path=/anvil/scratch/x-jnash12/Aus23_16S

qiime diversity alpha-rarefaction \
	--i-table ${WD_path}/16S_FeatureTable_180_120/table-no-mitochondria-no-chloroplast.qza \
	--p-max-depth 137314 \
	--p-metrics observed_features \
	--p-steps 200 \
	--p-iterations 100 \
	--o-visualization ${WD_path}/table-no-mitochondria-no-chloroplast_rarecurve100.qzv
```

## Alpha and beta diversity
The core-metrics command is used to generate a number of alpha and beta diversity outputs. A rarefaction depth of 24860 was used, which results in 1 sample being dropped (and also the negative control)

### RunCoreMetrics
```
#!/bin/bash
#SBATCH -o slurm-%j-CoreMetrics.out
#SBATCH -c 1
#SBATCH --partition=shared 
#SBATCH -A BIO230020
#SBATCH --export=ALL
#SBATCH -t 2:00:00

module load biocontainers
module load qiime2
WD_path=/anvil/scratch/x-jnash12/Aus23_16S

qiime diversity core-metrics \
  --i-table ${WD_path}/16S_FeatureTable_180_120/table-no-mitochondria-no-chloroplast.qza \
  --p-sampling-depth 24860 \
  --m-metadata-file ${WD_path}/MetabarcodingMetadata.txt \
  --output-dir ${WD_path}/core-metrics-results
```

## Alpha diversity significance tests
Conducts significance tests for alpha diversity metrics (richness, shannon, evenness)

### AlphaSignificance
```
#!/bin/bash
#SBATCH -o slurm-%j-CoreMetrics.out
#SBATCH -c 1
#SBATCH --partition=shared 
#SBATCH -A BIO230020
#SBATCH --export=ALL
#SBATCH -t 2:00:00

module load biocontainers
module load qiime2
WD_path=/anvil/scratch/x-jnash12/Aus23_16S

qiime diversity alpha-group-significance \
  --i-alpha-diversity ${WD_path}/core-metrics-results/observed_features_vector.qza \
  --m-metadata-file ${WD_path}/MetabarcodingMetadata.txt \
  --o-visualization ${WD_path}/core-metrics-results/richness-group-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity ${WD_path}/core-metrics-results/shannon_vector.qza \
  --m-metadata-file ${WD_path}/MetabarcodingMetadata.txt \
  --o-visualization ${WD_path}/core-metrics-results/shannon-group-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity ${WD_path}/core-metrics-results/evenness_vector.qza \
  --m-metadata-file ${WD_path}/MetabarcodingMetadata.txt \
  --o-visualization ${WD_path}/core-metrics-results/evenness-group-significance.qzv
```

## Creates taxonomy barplots using the silva and greengenes taxonomies
Generates qzv visualizations that can be used to interactively explore taxonomic differences between communities using GreenGenes2 and SILVA taxonomy

### TaxaBarplot
```
#!/bin/bash
#SBATCH -o slurm-%j-TaxaBarplot.out
#SBATCH -c 1
#SBATCH --partition=shared 
#SBATCH -A BIO230020
#SBATCH --export=ALL
#SBATCH -t 2:00:00

module load biocontainers
module load qiime2
WD_path=/anvil/scratch/x-jnash12/Aus23_16S

qiime taxa barplot \
  --i-table ${WD_path}/16S_FeatureTable_180_120/table-no-mitochondria-no-chloroplast.qza \
  --i-taxonomy ${WD_path}/Silva_16S_taxonomy.qza \
  --m-metadata-file ${WD_path}/MetabarcodingMetadata.txt \
  --o-visualization ${WD_path}/Silva_taxa_barplot.qzv

qiime taxa barplot \
  --i-table ${WD_path}/16S_FeatureTable_180_120/table-no-mitochondria-no-chloroplast.qza \
  --i-taxonomy ${WD_path}/GG_16S_taxonomy.qza \
  --m-metadata-file ${WD_path}/MetabarcodingMetadata.txt \
  --o-visualization ${WD_path}/GG_taxa_barplot.qzv
```

### ANCOMBC-TreeSpecies
Runs ancombc at ASV and Genus level to identify taxa that differ between Eucalyptus and Pine roots (removes all soil samples). FYI, right now this code isn't working - I think it's an installation issue
```
#!/bin/bash
#SBATCH -o slurm-%j-ANCOMBC.out
#SBATCH -c 1
#SBATCH --partition=shared 
#SBATCH -A BIO230020
#SBATCH --export=ALL
#SBATCH -t 2:00:00

module load biocontainers
module load qiime2
WD_path=/anvil/scratch/x-jnash12/Aus23_16S

qiime feature-table filter-samples \
	--i-table ${WD_path}/16S_FeatureTable_180_120/table-no-mitochondria-no-chloroplast.qza \
	--m-metadata-file ${WD_path}/MetabarcodingMetadata.txt \
	--p-where "[SampleType]='Root'" \
	--o-filtered-table ${WD_path}/table-no-mitochondria-no-chloroplast_RootOnly.qza

qiime taxa collapse \
  --i-table ${WD_path}/table-no-mitochondria-no-chloroplast_RootOnly.qza \
  --i-taxonomy ${WD_path}/Silva_16S_taxonomy.qza \
  --p-level 6 \
  --o-collapsed-table ${WD_path}/table-no-mitochondria-no-chloroplast_RootOnly_genus.qza

qiime composition ancombc \
  --i-table ${WD_path}/table-no-mitochondria-no-chloroplast_RootOnly_genus.qza  \
  --m-metadata-file ${WD_path}/MetabarcodingMetadata.txt \
  --p-formula 'TreeSpecies' \
  --o-differentials ${WD_path}/ancombc_TreeSpecies_genus.qza

qiime composition da-barplot \
  --i-data ${WD_path}/ancombc_TreeSpecies_genus.qza \
  --p-significance-threshold 0.01 \
  --o-visualization ${WD_path}/ancombc_TreeSpecies_genus.qzv

qiime composition ancombc \
  --i-table ${WD_path}/table-no-mitochondria-no-chloroplast_RootOnly.qza  \
  --m-metadata-file ${WD_path}/MetabarcodingMetadata.txt \
  --p-formula 'TreeSpecies' \
  --o-differentials ${WD_path}/ancombc_TreeSpecies_ASV.qza

qiime composition da-barplot \
  --i-data ${WD_path}/ancombc_TreeSpecies_ASV.qza \
  --p-significance-threshold 0.01 \
  --o-visualization ${WD_path}/ancombc_TreeSpecies_ASV.qzv
```

## Data Export
Now, we have generated the three output files necessary for downstream analyses - ASV table, representative sequences fasta file, and taxonomy table. These files need to be extracted from QIIME2 format into user-readable formats that can be downloaded to our local machines and imported into R

This script makes use of the QIIME2 export command to convert these files. I also do some copying/renaming to make sure files have informative names. Note that QIIME2's export command converts the ASV table to biom format (another weird format, sigh...) by default and we need to use the biom convert command to get this to a user-readable tsv format

### ExportData
```
#!/bin/bash
#SBATCH -o slurm-%j-ExportQIIME.out
#SBATCH -c 1
#SBATCH --partition=shared 
#SBATCH -A BIO230020
#SBATCH --export=ALL
#SBATCH -t 2:00:00

module load biocontainers
module load qiime2
WD_path=/anvil/scratch/x-jnash12/Aus23_16S

qiime tools export \
  --input-path ${WD_path}/16S_FeatureTable_180_120/table.qza \
  --output-path ${WD_path}/QIIME_16S_files

mv ${WD_path}/QIIME_16S_files/feature-table.biom ${WD_path}/QIIME_16S_files/Aus23_16S_ASV_table.biom

biom convert -i ${WD_path}/QIIME_16S_files/Aus23_16S_ASV_table.biom -o ${WD_path}/QIIME_16S_files/Aus23_16S_ASV_table.tsv --to-tsv

qiime tools export \
  --input-path ${WD_path}/16S_FeatureTable_180_120/representative_sequences.qza \
  --output-path ${WD_path}/QIIME_16S_files

mv ${WD_path}/QIIME_16S_files/dna-sequences.fasta ${WD_path}/QIIME_16S_files/Aus23_16S_ASV_repseqs.fasta

qiime tools extract \
  --input-path ${WD_path}/Silva_16S_taxonomy.qza \
  --output-path ${WD_path}/SilvaTaxonomy

qiime tools extract \
  --input-path ${WD_path}/GG_16S_taxonomy.qza \
  --output-path ${WD_path}/GGTaxonomy

cp ${WD_path}/SilvaTaxonomy/*/data/taxonomy.tsv ${WD_path}/QIIME_16S_files/Aus23_16S_SilvaTaxonomy_16S.tsv
cp ${WD_path}/GGTaxonomy/*/data/taxonomy.tsv ${WD_path}/QIIME_16S_files/Aus23_16S_GGTaxonomy_16S.tsv
```

### ImportData
If data needs to be imported back into QIIME (ASV table in biom format, taxonomy files, representative seqs), it can be done with this script

```
#!/bin/bash
#SBATCH -o slurm-%j-ImportQIIME.out
#SBATCH -c 1
#SBATCH --partition=shared 
#SBATCH -A BIO230020
#SBATCH --export=ALL
#SBATCH -t 2:00:00

module load biocontainers
module load qiime2
WD_path=/anvil/scratch/x-jnash12/Aus23_16S

mkdir ${WD_path}/16S_FeatureTable_180_120/

qiime tools import \
  --input-path ${WD_path}/Aus23_16S_ASV_table.biom \
  --type 'FeatureTable[Frequency]' \
  --input-format BIOMV210Format \
  --output-path ${WD_path}/16S_FeatureTable_180_120/table.qza

qiime tools import \
  --input-path ${WD_path}/Aus23_16S_ASV_repseqs.fasta  \
  --output-path ${WD_path}/16S_FeatureTable_180_120/representative_sequences.qza \
  --type 'FeatureData[Sequence]'

qiime tools import \
	--input-path ${WD_path}/Aus23_16S_SilvaTaxonomy_16S.tsv \
	--output-path ${WD_path}/Silva_16S_taxonomy.qza \
	--type FeatureData[Taxonomy]

qiime tools import \
	--input-path ${WD_path}/Aus23_16S_GGTaxonomy_16S.tsv \
	--output-path ${WD_path}/GG_16S_taxonomy.qza \
	--type FeatureData[Taxonomy]
```