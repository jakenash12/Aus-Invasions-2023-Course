# Australia NSF Course 2023 18S Data Analysis Pipeline
This is a set of scripts to analyze 18S metabarcoding data (amplified with modified primer pair WANDA (Dumbrell et al., 2011) and AML2 from Lee et al, 2008) using QIIME2. These scripts were written to run on the CU ALPINE cluster using QIIME2 with the q2cli version 2023.5.0 that is available as a module on the cluster.  

## Setup of environment and downloads
### Code to run at beginning of every session
Variable WD_path is set to the path to our working directory. Make sure to run this with each new session
```
WD_path=/scratch/alpine/jasigs@colostate.edu
```

### Creates working directory (do this only once - not with every new session)
```
mkdir /scratch/alpine/jasigs@colostate.edu/Aus23_18S
```

### Downloads demultiplexed sequence data from Vilgalys Lab Dropbox into working directory using rclone (downloaded directly from folder)
```
sbatch -o slurm-%j-DataDL.out --partition=shared --account=BIO230020 --export=ALL -t 2:00:00 -c 1 --wrap="rclone copy remote:'Vilgalys Lab'/NGS_RAW/Nash_8732_23101202 ${WD_path}/Nash_8732_23101202"
```

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
ls ${WD_path}/Nash_8732_23101202 | grep "R1_001.fastq.gz" | sed 's/_R1_001.fastq.gz//g' | grep "18S" > ${WD_path}/filelist_18S
printf "%s\t%s\t%s\n" "sample-id" "forward-absolute-filepath" "reverse-absolute-filepath" > ${WD_path}/QIIMEManifest.tsv
for i in $(cat ${WD_path}/filelist_18S)
do
Sample=$(echo "$i" | awk -F '_' '{print $1}')
SampleType=$(echo "$i" | awk -F '_' '{print $2}')
printf "%s\t%s\t%s\n" "${Sample}_${SampleType}" "${WD_path}/Nash_8732_23101202/${i}_R1_001.fastq.gz" "${WD_path}/Nash_8732_23101202/${i}_R2_001.fastq.gz" >> ${WD_path}/QIIMEManifest.tsv
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
#SBATCH --job-name=import
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --partition=amilan
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jasigs@colostate.edu

module purge
module load anaconda
conda activate qiime2-2023.5
WD_path=/scratch/alpine/jasigs@colostate.edu/Aus23_18S

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path ${WD_path}/QIIMEManifest.tsv \
  --output-path ${WD_path}/Aus_18S_demux.qza \
  --input-format PairedEndFastqManifestPhred33V2

qiime demux summarize \
  --i-data ${WD_path}/Aus_18S_demux.qza \
  --o-visualization ${WD_path}/Aus_18S_demux.qzv
 ```

## Cutadapt to trim primers and adapter sequences
Cutadapt can trim primers/adapters from both the 5' and 3' ends of the forward and reverse reads. Primers/adapters may occur on the 3' end in the event of read through (i.e. the 300 bp read extends to the end of the DNA fragment).

We provide reverse complements of the adapters and primers under the p-adapter flags in case read through happened. Normal primer sequences are provided in the p-front argument to trim them from the 5' end

Takes \~30 minutes

### RunCutadapt (P-adapter-f is the reverse complement of AML2, p-front-f is WANDA, p-adapter-r is the reverse complement of WANDA, p-front-r is AML2)
```
#!/bin/bash
#SBATCH --job-name=import
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --partition=amilan
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jasigs@colostate.edu

module purge
module load anaconda
conda activate qiime2-2023.5
WD_path=/scratch/alpine/jasigs@colostate.edu/Aus23_18S

qiime cutadapt trim-paired \
	--i-demultiplexed-sequences ${WD_path}/Aus_18S_demux.qza \
	--p-adapter-f GGAAACCAAAGTGTTTGGGTTC \
	--p-adapter-f AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
	--p-front-f CAGCCGCGGTAATTCCAGCT \
	--p-front-f GCCTCCCTCGCGCCATCAGAGATGTGTATAAGAGACAG \
	--p-adapter-r AGCTGGAATTACCGCGGCTG \
	--p-adapter-r CTGTCTCTTATACACATCTCTGATGGCGCGAGGGAGGC \
	--p-front-r GAACCCAAACACTTTGGTTTCC \
	--p-front-r GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT \
	--p-error-rate 0.15 \
	--output-dir ${WD_path}/Aus_18S_demux_cutadapt \
	--verbose

qiime demux summarize \
  --i-data ${WD_path}/Aus_18S_demux_cutadapt/trimmed_sequences.qza \
  --o-visualization ${WD_path}/Aus_18S_demux_cutadapt.qzv
```
## Need to attempt with differing error rates (e.g., 0.05, 0.10)

## Dada2 for ASV calling + QC
Dada2 is a pipeline for sequence error correction, chimera detection, read merging, and ASV calling. Output is a "feature table" (analagous to an OTU table) and a set of representative DNA sequences

Dada2 has very stringent quality filtering standards and will discard reads that do not meet. The --p-trunc-len-f and --p-trunc-len-r parameters allow us to trim off low quality regions at the end of the forward and reverse read, respectively, which helps more reads pass QC. The danger is if we trim too much, then the reads will not be long enough to merge. By inspecting the sequence quality plot after cutadapt, we can determine optimal parameters for this step.

The WANDA-AML2 amplicon region of 18S that we amplified is \~550bp after trimming primers/adapters. Reverse reads are often lower quality than forward reads and usually require trimming at an earlier position.  Dada2 requires 20-50 bp of overlap between the reads. Here, we tested 4 different primer trimming settings as indicated below to determine which resulted in the highest percentage of reads being retained.

Given the lack of overlap, we will attempt untrimmed, ~20bp, and tbd based on reverse read quality

######## STOPPING HERE ON 12/18/2023 #########

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

## Creates taxonomy barplot using the silva taxonomy
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
  --input-path ${WD_path}/${WD_path}/16S_FeatureTable_180_120/table.qza \
  --output-path ${WD_path}/QIIME_16S_files

mv ${WD_path}/QIIME_16S_files/feature-table.biom ${WD_path}/QIIME_16S_files/Aus23_16S_ASV_table.biom

biom convert -i ${WD_path}/QIIME_16S_files/Aus23_16S_ASV_table.biom -o ${WD_path}/QIIME_16S_files/Aus23_16S_ASV_table.tsv --to-tsv

qiime tools export \
  --input-path ${WD_path}/${WD_path}/16S_FeatureTable_180_120/representative_sequences.qza \
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


