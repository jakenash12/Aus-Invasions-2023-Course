# Australia NSF Course 2023 18S Data Analysis Pipeline (w/ V-Xtractor)
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

### Only V4 regions of SSU are present in sequence count

# Vsearch Pipeline
# Merging paired end reads to a useable downstream artifact
qiime vsearch merge-pairs \
  --i-demultiplexed-seqs Aus_18S_demux.vxtractor.qza \
  --o-merged-sequences vsearch_merged_reads.vx.qza

# Dereplicating sequences (does not include any quality filtering)
qiime vsearch dereplicate-sequences \
  --i-sequences vsearch_merged_reads.vx.qza \
  --o-dereplicated-table vsearch_table.vx.qza \
  --o-dereplicated-sequences vsearch_rep-seqs.vx.qza

# Visualize rep seqs
#!/bin/bash
#SBATCH --job-name=tabseqs.vx
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --partition=amilan
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jasigs@colostate.edu

module purge
module load anaconda
conda activate qiime2-2023.5

qiime feature-table tabulate-seqs \
--i-data vsearch_rep-seqs.vx.qza \
--o-visualization vsearch_rep-seqs.vx.qzv

# Visualize feature table (1,514,680 features)
qiime feature-table summarize \
--i-table vsearch_table.vx.qza \
--o-visualization vsearch_table.vx.qzv \
--m-sample-metadata-file MetabarcodingMetadata.txt

### CLUSTERING AND TAXONOMY ###

# Open reference clustering at 95% (with Maarjam database)- Ran
#!/bin/bash
#SBATCH --job-name=Maarj95
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --partition=amilan
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jasigs@colostate.edu

module purge
module load anaconda
conda activate qiime2-2023.5

qiime vsearch cluster-features-open-reference \
  --i-table vsearch_table.vx.qza \
  --i-sequences vsearch_rep-seqs.vx.qza \
  --i-reference-sequences MaarjAM-ref-seqs.qza \
  --p-perc-identity 0.95 \
  --o-clustered-table table-maarjam-95.qza \
  --o-clustered-sequences rep-seqs-maarjam-95.qza \
  --o-new-reference-sequences new-ref-seqs-maarjam-95.qza


# Visualize Maarj95-Vsearch clustered feature table 
qiime feature-table summarize \
--i-table table-maarjam-95.qza \
--o-visualization table-maarjam-95.qzv \
--m-sample-metadata-file MetabarcodingMetadata.txt

# Visualize Maarj95-Vsearch rep seqs
qiime feature-table tabulate-seqs \
--i-data rep-seqs-maarjam-95.qza \
--o-visualization rep-seqs-maarjam-95.qzv

# Open reference clustering of unassigned 95% seqs at 90% (with Maarjam database)- Didn't work
#!/bin/bash
#SBATCH --job-name=Maarj90-UA
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --partition=amilan
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jasigs@colostate.edu

module purge
module load anaconda
conda activate qiime2-2023.5

qiime vsearch cluster-features-open-reference \
  --i-table vsearch_table.vx.qza \
  --i-sequences new-ref-seqs-maarjam-95.qza \
  --i-reference-sequences MaarjAM-ref-seqs.qza \
  --p-perc-identity 0.90 \
  --o-clustered-table table-maarjam-90.qza \
  --o-clustered-sequences rep-seqs-maarjam-90.qza \
  --o-new-reference-sequences new-ref-seqs-maarjam-90.qza

# Trying another feature classifier command (timed out-trying again; this is what Kacie used)
#!/bin/bash
#SBATCH --job-name=Maarj95-fc
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --partition=amilan
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jasigs@colostate.edu

module purge
module load anaconda
conda activate qiime2-2023.5

qiime feature-classifier classify-consensus-blast \
--i-query vsearch_rep-seqs.vx.qza \
--i-reference-reads MaarjAM-ref-seqs.qza \
--i-reference-taxonomy maarjam-ref-tax.qza \
--p-maxaccepts 1 \
--p-perc-identity 0.95 \
--p-query-cov 0.90 \
--p-strand both \
--p-evalue 1e-50 \
--p-min-consensus 0.51 \
--output-dir maarjam-95-fc

# Trying the OG classification route
#!/bin/bash
#SBATCH --job-name=classifier
#SBATCH --nodes=1
#SBATCH --ntasks=50
#SBATCH --partition=amilan
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jasigs@colostate.edu

#Activate qiime

module purge

module load anaconda

conda activate qiime2-2023.5

#Command
qiime feature-classifier fit-classifier-naive-bayes \
--i-reference-reads MaarjAM-ref-seqs.qza \
--i-reference-taxonomy maarjam-ref-tax.qza \
--o-classifier maarjam-classifier.qza
