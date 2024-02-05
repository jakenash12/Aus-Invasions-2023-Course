# Australia NSF Course 2023 18S v-xtractor
These scripts use v-xtractor to identify reads that contain the 18S region of interest and discards off-target sequences. Because v-xtractor only works with single fasta reads, reads are first converted to fasta files, then merged using PEAR. The fasta headers from the merged reads that are identified as 18S can then be used to subset the forward and reverse read files

## Setup of environment and downloads
### Code to run at beginning of every session
Variable WD_path is set to the path to our working directory. Make sure to run this with each new session
```
WD_path=/anvil/scratch/x-jnash12/Aus23_18S
cd ${WD_path}
```

### Downloads demultiplexed sequence data from Vilgalys Lab Dropbox into working directory using rclone
```
sbatch -o slurm-%j-DataDL.out --partition=shared --account=BIO230020 --export=ALL -t 2:00:00 -c 1 --wrap="rclone copy remote:'Vilgalys Lab'/NGS_RAW/Nash_8732_23101202 ${WD_path}/Nash_8732_23101202"
```

### Creates filelist containing base names of 18S files only (not including 16S or ITS2 files)
```
ls ${WD_path}/Nash_8732_23101202 | grep "18S" | grep "R1_001.fastq.gz" | sed 's/_R1_001.fastq.gz//g' > ${WD_path}/filelist_18S
```

### uses PEAR to merge reads
```
mkdir ${WD_path}/PEARReads/

for i in $(cat ${WD_path}/filelist_18S)
do
sbatch -o slurm-%j-PEAR18S.out --partition=shared --account=BIO230020 --export=ALL -t 1:00:00 -c 1 --wrap="pear \
	-f ${WD_path}/Nash_8732_23101202/${i}_R1_001.fastq.gz \
	-r ${WD_path}/Nash_8732_23101202/${i}_R2_001.fastq.gz \
	-o ${WD_path}/PEARReads/${i}"
done
```

### converts merged fastq files to fasta using seqret
```
for i in $(cat ${WD_path}/filelist_18S)
do
	sbatch -o slurm-%j-fastaconvert.out --partition=shared --account=BIO230020 --export=ALL -t 1:00:00 -c 1 --wrap="seqtk seq -a ${WD_path}/PEARReads/${i}.assembled.fastq > ${WD_path}/PEARReads/${i}.assembled.fasta"
done
```

### downloads v-xtractor and uses it to extract 18S V4 region. FYI, when installing hmmer (a dependency of V-xtracter), make sure to install version 3.0 because I found that some newer versions (specifically v3.4) do not work
```
mkdir ${WD_path}/vxtractor/
cd ${WD_path}/vxtractor/
git clone https://github.com/carden24/V-Xtractor
unzip ${WD_path}/vxtractor/V-Xtractor/HMMs.zip

for i in $(cat ${WD_path}/filelist_18S)
do
	sbatch -o slurm-%j-VxtractorV4short.out --partition=shared --account=BIO230020 --export=ALL -t 6:00:00 -c 1 --wrap="perl ${WD_path}/vxtractor/V-Xtractor/vxtractor.pl -h ${WD_path}/vxtractor/V-Xtractor/HMMs/fungi/ -r V4 -i long -o ${WD_path}/vxtractor/${i}.assembled.vxtractorV4.fasta -c ${WD_path}/vxtractor/${i}.assembled.vxtractorV4.csv ${WD_path}/PEARReads/${i}.assembled.fasta"
done
```

### Pulls out sequence headers of sequences identified as 18S by v-xtractor, then uses those to subset raw fastq (unmerged) files
```
for i in $(cat ${WD_path}/filelist_18S)
do
grep ">" ${WD_path}/vxtractor/${i}.assembled.vxtractorV4.fasta | sed -E 's/^([^_]*).*$/\1/g' | sed 's/>//g' > ${WD_path}/vxtractor/${i}.assembled.vxtractorV4.headersR1
done

for i in $(cat ${WD_path}/filelist_18S)
do
grep ">" ${WD_path}/vxtractor/${i}.assembled.vxtractorV4.fasta | sed -E 's/^([^_]*).*$/\1/g' | sed 's/>//g' | sed 's/ 1/ 2/g'  > ${WD_path}/vxtractor/${i}.assembled.vxtractorV4.headersR2
done

mkdir ${WD_path}/UnmergedVxtractor
```

### SubsetV4
```
#!/bin/bash
#SBATCH -o slurm-%j-R1Subset.out
#SBATCH -c 1
#SBATCH --partition=shared 
#SBATCH -A BIO230020
#SBATCH --export=ALL
#SBATCH -t 4:00:00

module load biocontainers
module load bbmap

WD_path=/anvil/scratch/x-jnash12/Aus23_18S

for i in $(cat ${WD_path}/filelist_18S)
do
filterbyname.sh in=${WD_path}/Nash_8732_23101202/${i}_R1_001.fastq.gz out=${WD_path}/UnmergedVxtractor/${i}_R1_001.V4.fastq.gz names=${WD_path}/vxtractor/${i}.assembled.vxtractorV4.headersR1 include=t
done

for i in $(cat ${WD_path}/filelist_18S)
do
filterbyname.sh in=${WD_path}/Nash_8732_23101202/${i}_R2_001.fastq.gz out=${WD_path}/UnmergedVxtractor/${i}_R2_001.V4.fastq.gz names=${WD_path}/vxtractor/${i}.assembled.vxtractorV4.headersR2 include=t
done
```

### Verifies that forward and reverse read files contain same number of reads
```
T=$(printf '\t')
for i in $(cat ${WD_path}/filelist_18S)
do
Sample=$(echo "${i%%_*}")
For=$(gunzip -c ${WD_path}/UnmergedVxtractor/${i}_R1_001.V4.fastq.gz | wc -l)
Rev=$(gunzip -c ${WD_path}/UnmergedVxtractor/${i}_R2_001.V4.fastq.gz | wc -l)
echo "${Sample}${T}${For}${T}${Rev}"
done
```

### Generates tsv summary file displaying number of reads retained at each step
```
T=$(printf '\t')
echo "Sample${T}Raw${T}Merged${T}ExtractedV4" > ${WD_path}/ReadCounts18S.tsv

for i in $(cat ${WD_path}/filelist_18S)
do
Sample=$(echo "${i%%_*}")
Raw=$(gunzip -c ${WD_path}/Nash_8732_23101202/${i}_R1_001.fastq.gz | echo $((`wc -l`/4)))
Merged=$(grep ">" ${WD_path}/PEARReads/${i}.assembled.fasta | wc -l)
ExtractedV4=$(cat ${WD_path}/vxtractor/${i}.assembled.vxtractorV4.headers | wc -l)
echo "${Sample}${T}${Raw}${T}${Merged}${T}${ExtractedV4}" >> ${WD_path}/ReadCounts18S.tsv
done
```


