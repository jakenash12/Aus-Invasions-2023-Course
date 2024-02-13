# Australia NSF Course 2023 18S Data Analysis Pipeline (w/ V-Xtractor & dada2)
This is a set of scripts to analyze 18S metabarcoding data (amplified with modified primer pair WANDA (Dumbrell et al., 2011) and AML2 from Lee et al, 2008) using QIIME2. These scripts were written to run on the CU ALPINE cluster using QIIME2 with the q2cli version 2023.5.0 that is available as a module on the cluster.  

## Setup of environment and downloads
### Code to run at beginning of every session
Variable WD_path is set to the path to our working directory. Make sure to run this with each new session
```
WD_path=/scratch/alpine/jasigs@colostate.edu/Aus23_18S
```

### Creates working directory (do this only once - not with every new session)
```
mkdir /scratch/alpine/jasigs@colostate.edu/Aus23_18S
```
### Files transferred via Globus transfer following VXtractor pipeline to isolate V4 region of SSU 
### See https://github.com/jakenash12/Aus-Invasions-2023-Course/blob/main/Aus23_18S_Vxtractor.md for steps of this pipeline

### Downloads sklearn taxonomic classifiers for Marjaam 
Oftentimes, QIIME2 users will train their own sklearn taxonomic classifiers based on reference databases trimmed to the gene region they amplified. We will need to train the MaarjAM classifier

## MaarjAM database citation
Öpik, M., Vanatoa, A., Vanatoa, E., Moora, M., Davison, J., Kalwij, J.M., Reier, Ü., Zobel, M. 2010. The online database MaarjAM reveals global and ecosystemic distribution patterns in arbuscular mycorrhizal fungi (Glomeromycota). New Phytologist 188: 223-241.

# Pulling from MaarjAM website
```
cd $/scratch/alpine/jasigs@colostate.edu/Aus23_18S
wget https://maarjam.ut.ee/resources/maarjam_database_SSU.qiime.2019.zip
unzip maarjam_database_SSU.qiime.2019.zip
```

### Preparing Manifest File
QIIME2 uses a special data format with the file type qza to store sequence project data. Fastq files must be imported into QIIME2 first. Because our data are already demultiplexed, we need to generate a tsv file that links sample IDs to the filepaths of the forward and reverse reads to be used in the import command.

```
ls ${WD_path}/VXtractorReads | grep "R1_001.V4.fastq.gz" | sed 's/_R1_001.V4.fastq.gz//g' | grep "18S" > ${WD_path}/filelist_18S
printf "%s\t%s\t%s\n" "sample-id" "forward-absolute-filepath" "reverse-absolute-filepath" > ${WD_path}/QIIMEManifest.Vxtractor.tsv
for i in $(cat ${WD_path}/filelist_18S)
do
Sample=$(echo "$i" | awk -F '_' '{print $1}')
SampleType=$(echo "$i" | awk -F '_' '{print $2}')
printf "%s\t%s\t%s\n" "${Sample}_${SampleType}" "${WD_path}/VXtractorReads/${i}_R1_001.V4.fastq.gz" "${WD_path}/VXtractorReads/${i}_R2_001.V4.fastq.gz" >> ${WD_path}/QIIMEManifest.Vxtractor.tsv
done
```
# QIIME Time
All following commands will be done within QIIME2. QIIME2 is available from CU ALPINE as a module that can be loaded. To do this, we create batch jobs for each command and load the QIIME2 module at the beginning of each script. To run the following chunks of code, you will need to create files using the nano command with the script name as noted in bold at top of code chunk, and then paste each chunk into the file. Then scripts can be run by typing sbatch {scriptname}

## Imports data into QIIME 
Uses manifest to import paired end sequence data into QIIME2 then creates a sequencing summary visualization in qzv format that can be downloaded onto our local machine and viewed at this link https://view.qiime2.org/

Takes a couple minutes

### DataImport
```

#!/bin/bash
#SBATCH --job-name=import.vx
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --partition=amilan
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jasigs@colostate.edu

module purge
module load anaconda
conda activate qiime2-2023.5

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path QIIMEManifest.Vxtractor.tsv \
  --output-path Aus_18S_demux.vxtractor.qza \
  --input-format PairedEndFastqManifestPhred33V2

qiime demux summarize \
  --i-data Aus_18S_demux.vxtractor.qza \
  --o-visualization Aus_18S_demux.vxtractor.qzv
 ```

### Denoising with dada2 and skipping cutadapt- just trimming primers off
## Parameters for trim-left-f and trim-left-r correspond to lengths of WANDA and AML2, respectively
qiime dada2 denoise-paired \
--i-demultiplexed-seqs Aus_18S_demux.vxtractor.qza \
--p-trim-left-f 20 \
--p-trim-left-r 22 \
--p-trunc-len-f 300 \
--p-trunc-len-r 300 \
--o-table Aus_18S_demux.vxtractor.table.qza \
--o-representative-sequences Aus_18S_demux.vxtractor.dada2-rep-seqs.qza \
--o-denoising-stats vxtractor-DADA2-stats.qza \
--verbose

# Visualize run stats
# Visualizing outputs
qiime metadata tabulate \
--m-input-file vxtractor-DADA2-stats.qza  \
--o-visualization vxtractor-DADA2-stats.qzv

# Visualize feature table
qiime feature-table summarize \
--i-table Aus_18S_demux.vxtractor.table.qza \
--o-visualization Aus_18S_demux.vxtractor.table.qzv \
--m-sample-metadata-file MetabarcodingMetadata.txt

# Visualize rep seqs
qiime feature-table tabulate-seqs \
--i-data Aus_18S_demux.vxtractor.dada2-rep-seqs.qza \
--o-visualization Aus_18S_demux.vxtractor.dada2-rep-seqs.qzv

# Trying feature classifier on all sequences at 95%
qiime feature-classifier classify-consensus-blast \
--i-query Aus_18S_demux.vxtractor.dada2-rep-seqs.qza \
--i-reference-reads MaarjAM-ref-seqs.qza \
--i-reference-taxonomy maarjam-ref-tax.qza \
--p-maxaccepts 1 \
--p-perc-identity 0.95 \
--p-query-cov 0.90 \
--p-strand both \
--p-evalue 1e-50 \
--p-min-consensus 0.51 \
--output-dir maarjam-95

# Visualizing classified sequences
qiime metadata tabulate \
--m-input-file maarjam-95/classification.qza \
--o-visualization maarjam-95/classification.qzv

# Pulling unassigned sequences to reassign at 90%
qiime taxa filter-seqs \
--i-sequences Aus_18S_demux.vxtractor.dada2-rep-seqs.qza \
--i-taxonomy maarjam-95/maarj95_classification.qza \
--p-include Unassigned \
--o-filtered-sequences maarjam-95/unassigned-rep-seqs-95.qza

# Reassigning unclassified seqs at 90%
qiime feature-classifier classify-consensus-blast \
--i-query maarjam-95/unassigned-rep-seqs-95.qza \
--i-reference-reads MaarjAM-ref-seqs.qza \
--i-reference-taxonomy maarjam-ref-tax.qza \
--p-maxaccepts 1 \
--p-perc-identity 0.90 \
--p-query-cov 0.90 \
--p-strand both \
--p-evalue 1e-50 \
--p-min-consensus 0.51 \
--output-dir maarjam-90

# Visualizing classified sequences
qiime metadata tabulate \
--m-input-file maarjam-90/classification.qza \
--o-visualization maarjam-90/classification.qzv

# Pulling unassigned sequences to reassign at 80%
qiime taxa filter-seqs \
--i-sequences maarjam-95/unassigned-rep-seqs-95.qza \
--i-taxonomy maarjam-90/maarj90_classification.qza \
--p-include Unassigned \
--o-filtered-sequences maarjam-90/unassigned-rep-seqs-90.qza

# Reassigning unclassified seqs at 80%
qiime feature-classifier classify-consensus-blast \
--i-query maarjam-90/unassigned-rep-seqs-90.qza \
--i-reference-reads MaarjAM-ref-seqs.qza \
--i-reference-taxonomy maarjam-ref-tax.qza \
--p-maxaccepts 1 \
--p-perc-identity 0.80 \
--p-query-cov 0.90 \
--p-strand both \
--p-evalue 1e-50 \
--p-min-consensus 0.51 \
--output-dir maarjam-80

# Visualizing classified sequences
qiime metadata tabulate \
--m-input-file maarjam-80/classification.qza \
--o-visualization maarjam-80/classification.qzv

#### Filtering out non-AMF taxa from our tables
# 95% table
qiime taxa filter-table \
--i-table Aus_18S_demux.vxtractor.table.qza \
--i-taxonomy maarjam-95/maarj95_classification.qza \
--p-include Glomeromycetes,Paraglomeromycetes \
--o-filtered-table maarjam-95/glom-only-maarjam-table-95.qza

# Filtering the unassigned sequences from 95% table
qiime taxa filter-table \
--i-table Aus_18S_demux.vxtractor.table.qza \
--i-taxonomy maarjam-95/maarj95_classification.qza \
--p-include Unassigned \
--o-filtered-table maarjam-95/unassigned-table-95.qza

# Filtering 90% table
qiime taxa filter-table \
--i-table maarjam-95/unassigned-table-95.qza \
--i-taxonomy maarjam-90/maarj90_classification.qza \
--p-include Glomeromycetes,Paraglomeromycetes \
--o-filtered-table maarjam-90/glom-only-maarjam-table-90.qza

# Filtering unassigned sequences from 95% table
qiime taxa filter-table \
--i-table maarjam-95/unassigned-table-95.qza \
--i-taxonomy maarjam-90/maarj90_classification.qza \
--p-include Unassigned \
--o-filtered-table maarjam-90/unassigned-table-90.qza

# Filter 80% table
qiime taxa filter-table \
--i-table maarjam-90/unassigned-table-90.qza \
--i-taxonomy maarjam-80/maarj80_classification.qza \
--p-include Glomeromycetes,Paraglomeromycetes \
--o-filtered-table maarjam-80/glom-only-maarjam-table-80.qza

# Making individual feature tables from AMF only
qiime taxa barplot \
--i-table maarjam-95/glom-only-maarjam-table-95.qza \
--i-taxonomy maarjam-95/maarj95_classification.qza \
--m-metadata-file MetabarcodingMetadata.txt \
--o-visualization maarjam-95/glom-only-taxa-barplot-95.qzv

qiime taxa barplot \
--i-table maarjam-90/glom-only-maarjam-table-90.qza \
--i-taxonomy maarjam-90/maarj90_classification.qza \
--m-metadata-file MetabarcodingMetadata.txt \
--o-visualization maarjam-90/glom-only-taxa-barplot-90.qzv

qiime taxa barplot \
--i-table maarjam-80/glom-only-maarjam-table-80.qza \
--i-taxonomy maarjam-80/maarj80_classification.qza \
--m-metadata-file MetabarcodingMetadata.txt \
--o-visualization maarjam-80/glom-only-taxa-barplot-80.qzv

# Merging the 3 feature tables
qiime feature-table merge \
--i-tables maarjam-95/glom-only-maarjam-table-95.qza \
--i-tables maarjam-90/glom-only-maarjam-table-90.qza \
--i-tables maarjam-80/glom-only-maarjam-table-80.qza \
--o-merged-table glom-only-final-table.qza \
--p-overlap-method sum

# Visualize feature table
qiime feature-table summarize \
--i-table glom-only-final-table.qza \
--o-visualization glom-only-final-table.qzv \
--m-sample-metadata-file MetabarcodingMetadata.txt