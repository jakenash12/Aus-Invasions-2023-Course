# Australia NSF Course 2023 16S Data Analysis Pipeline

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
QIIME2 uses a special data format with the file type qza to store sequence project data. Fastq files must be imported into QIIME2 first. Because our data are already demultiplexed, we need to generate a tsv file that links sample IDs to the filepaths of the forward and reverse reads.

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
All following commands will be done within QIIME2. QIIME2 is available from Purdue ANVIL as a module that can be loaded. To do this, we create batch jobs for each command and load the QIIME2 module at the beginning of each script 

## Imports data into QIIME 
Uses manifest to import paired end sequence data into QIIME2 then creates a sequencing summary visualization in qzv format that can be downloaded onto our local machine and viewed at this link https://view.qiime2.org/

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

Takes \~20 minutes

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


```
#!/bin/bash
#SBATCH -o slurm-%j-cutadapt.out
#SBATCH -c 16
#SBATCH --partition=shared 
#SBATCH -A BIO230020
#SBATCH --export=ALL
#SBATCH -t 12:00:00

qiime dada2 denoise-paired \
	--i-demultiplexed-seqs ${WD_path}/Aus_16S_demux_cutadapt/trimmed_sequences.qza \
	--p-trunc-len-f 220 \
	--p-trunc-len-r 120 \
	--p-trim-left-f 0 \
	--p-trim-left-r 0 \
	--p-n-threads 0 \
	--output-dir ${WD_path}/16S_FeatureTable_120R
```