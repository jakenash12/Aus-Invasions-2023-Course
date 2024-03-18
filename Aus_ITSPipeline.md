### Setup of working directory
```
WD_path=/anvil/scratch/x-jnash12/Aus_ITS
mkdir ${WD_path}
cd ${WD_path}
```

### Downloads demultiplexed sequence data from Vilgalys Lab Dropbox into working directory using rclone
```
sbatch -o slurm-%j-DataDL.out --partition=shared --account=BIO230020 --export=ALL -t 2:00:00 -c 1 --wrap="rclone copy remote:'Vilgalys Lab'/NGS_RAW/Nash_8732_23101202 ${WD_path}/Nash_8732_23101202"
```

### Generates filelist to loop through, filtering to include only ITS sequences
```
ls ${WD_path}/Nash_8732_23101202" | grep "R1_001.fastq.gz" | grep "ITS" | sed 's/_R1_001.fastq.gz//g' > ${WD_path}/filelist
```

### Uses PEAR to merge paired reads
```
mkdir ${WD_path}/PEARReads/ && cd ${WD_path}/PEARReads/
for i in $(cat ${WD_path}/filelist)
do
sbatch -o slurm-%j-PEAR.out --partition=shared --account=BIO230020 --export=ALL -t 2:00:00 -c 1 --wrap="pear \
	-f /anvil/scratch/x-jnash12/Nash_8732_23101202/${i}_R1_001.fastq.gz \
	-r /anvil/scratch/x-jnash12/Nash_8732_23101202/${i}_R2_001.fastq.gz \
	-o ${WD_path}/PEARReads/${i}"
done
```

### Uses ITSxpress to extract ITS2 region from merged reads
ITSxpress uses a pre-trained database to identify sequences that contain the ITS2 region, then trims those sequences to remove the adjacent 5.8S and LSU sequences, which are more conserved and less useful for clustering and taxonomic identification
```
conda activate ITSxpress
mkdir ${WD_path}/ITSxpressReads/ && cd ${WD_path}/ITSxpressReads/
for i in $(cat ${WD_path}/filelist)
do
sbatch -o slurm-%j-ITSxpress.out --partition=shared --account=BIO230020 --export=ALL -t 24:00:00 -c 128 --wrap="itsxpress --fastq ${WD_path}/PEARReads/${i}.assembled.fastq --single_end --region ITS2 --taxa Fungi \
--log ${WD_path}/ITSxpressReads/${i}.logfile.txt --outfile ${WD_path}/ITSxpressReads/${i}.ITS.fastq.gz --threads 128"
done
```

### Prepares QIIME2 manifest using custom script.
The QIIME2 manifest is a tab separated value file that points to the file locations for each of the demultiplexed sequence files and associates that file with a unique sample ID
```
printf "%s\t%s\n" "sample-id" "absolute-filepath" > ${WD_path}/QIIMEManifest.tsv
for i in $(cat ${WD_path}/filelist)
do
SampleID=$(echo "$i" | grep -oP 'Mengyi.')
printf "%s\t%s\n" "${SampleID}" "${WD_path}/ITSxpressReads/${i}.ITS.fastq.gz" >> ${WD_path}/QIIMEManifest.tsv
done
```

### Imports data into QIIME2 using manifest to point to file locations
```
#!/bin/bash
#SBATCH -o slurm-%j-QIIME_Import.out
#SBATCH -c 1
#SBATCH --partition=shared 
#SBATCH -A BIO230020
#SBATCH --export=ALL
#SBATCH -t 4:00:00

module load biocontainers
module load qiime2
WD_path=/anvil/scratch/x-jnash12/Aus_ITS

qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path ${WD_path}/QIIMEManifest.tsv \
  --output-path ${WD_path}/ITS2_demux.qza \
  --input-format SingleEndFastqManifestPhred33V2

qiime demux summarize \
  --i-data ${WD_path}/ITS2_demux.qza \
  --o-visualization ${WD_path}/ITS2_demux.qzv

 ```

### Uses Dada2 to denoise, quality filter, chimera filter, and ASV call
```
#!/bin/bash
#SBATCH -o slurm-%j-dada2.out
#SBATCH -c 128
#SBATCH --partition=shared 
#SBATCH -A BIO230020
#SBATCH --export=ALL
#SBATCH -t 48:00:00

module load biocontainers
module load qiime2
WD_path=/anvil/scratch/x-jnash12/Aus_ITS

qiime dada2 denoise-single \
	--i-demultiplexed-seqs ${WD_path}/ITS2_demux.qza \
	--p-trunc-len 0 \
	--output-dir ${WD_path}/dada2out \
	--p-n-threads 0 \
	--o-table ${WD_path}/ITS2_dada2table.qza \
	--o-representative-sequences ${WD_path}/ITS2_dada2seqs.qza \
	--o-denoising-stats ${WD_path}/ITS2_dada2denoising.qza

qiime metadata tabulate \
	--m-input-file ${WD_path}/ITS2_dada2denoising.qza \
	--o-visualization ${WD_path}/ITS2_dada2denoising.qzv

qiime feature-table summarize \
	--i-table ${WD_path}/ITS2_dada2table.qza \
	--o-visualization ${WD_path}/ITS2_dada2table.qzv

qiime feature-table tabulate-seqs \
	--i-data ${WD_path}/ITS2_dada2seqs.qza \
	--o-visualization ${WD_path}/ITS2_dada2seqs.qzv
```


### Runs Vsearch clustering to generate 97% OTUs
Output is our final OTU table and a set of associated representative sequences that are used for taxonomic classification
```
#!/bin/bash
#SBATCH -o slurm-%j-vsearch_cluster.out
#SBATCH -c 128
#SBATCH --partition=shared 
#SBATCH -A BIO230020
#SBATCH --export=ALL
#SBATCH -t 6:00:00

module load biocontainers
module load qiime2
WD_path=/anvil/scratch/x-jnash12/Aus_ITS

qiime vsearch cluster-features-de-novo \
	--i-sequences ${WD_path}/ITS2_dada2seqs.qza\
	--i-table ${WD_path}/ITS2_dada2table.qza \
	--p-perc-identity 0.97 \
	--p-threads 0 \
	--o-clustered-table ${WD_path}/ITS2_Dada2_table97.qza \
	--o-clustered-sequences ${WD_path}/ITS2_Dada2_repseqs97.qza

qiime feature-table summarize \
	--i-table ${WD_path}/ITS2_Dada2_table97.qza \
	--o-visualization ${WD_path}/ITS2_Dada2_table97.qzv

qiime feature-table tabulate-seqs \
	--i-data ${WD_path}/ITS2_Dada2_repseqs97.qza \
	--o-visualization ${WD_path}/ITS2_Dada2_repseqs97.qzv
```

### Taxonomic Classification using QIIME2's sklearn machine learning classifier
Uses a pre-trained taxonomic classifier for the UNITE database of full length ITS sequences
```
#!/bin/bash
#SBATCH -o slurm-%j-sklearn_classify.out
#SBATCH -c 128
#SBATCH --partition=shared 
#SBATCH -A BIO230020
#SBATCH --export=ALL
#SBATCH -t 6:00:00

module load biocontainers
module load qiime2
WD_path=/anvil/scratch/x-jnash12/Aus_ITS
cd ${WD_path}

wget https://github.com/colinbrislawn/unite-train/releases/download/v9.0-v25.07.2023-qiime2-2023.9/unite_ver9_dynamic_all_25.07.2023-Q2-2023.9.qza -O unite_ver9_dynamic_all_25.07.2023-Q2-2023.9.qza

qiime feature-classifier classify-sklearn \
  --i-classifier ${WD_path}/unite_ver9_dynamic_all_25.07.2023-Q2-2023.9.qza \
  --i-reads ${WD_path}/ITS2_Dada2_repseqs97.qza \
  --o-classification ${WD_path}/ITS2_Dada2_repseqs97_taxonomy.qza \
  --p-n-jobs -1
```

### ExportData
Converts QIIME2 format files into tsv files that can be downloaded and read into R for statistics
```
#!/bin/bash
#SBATCH -o slurm-%j-ExportQIIME.out
#SBATCH -c 1
#SBATCH --partition=shared 
#SBATCH -A BIO230020
#SBATCH --export=ALL
#SBATCH -t 1:00:00

module load biocontainers
module load qiime2
WD_path=/anvil/scratch/x-jnash12/Aus_ITS

qiime tools extract \
  --input-path ${WD_path}/ITS2_Dada2_repseqs97_taxonomy.qza \
  --output-path ${WD_path}/ITS2_Dada2_repseqs97_taxonomy

qiime tools export \
  --input-path ${WD_path}/ITS2_Dada2_table97.qza \
  --output-path ${WD_path}/QIIME_exported_files

biom convert -i ${WD_path}/QIIME_exported_files/feature-table.biom -o ${WD_path}/QIIME_exported_files/ITS2_OTUTable_97.tsv --to-tsv

cp ${WD_path}/ITS2_Dada2_repseqs97_taxonomy/*/data/taxonomy.tsv ${WD_path}/QIIME_exported_files/ITS2_Dada2_repseqs97_taxonomy.tsv

qiime tools export \
  --input-path ${WD_path}/ITS2_Dada2_repseqs97.qza \
  --output-path ${WD_path}/QIIME_exported_files

mv ${WD_path}/QIIME_exported_files/dna-sequences.fasta ${WD_path}/QIIME_exported_files/ITS2_Dada2_repseqs97.fasta

```