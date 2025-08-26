# Processing MAGs

### Overview
There are several things we do with our MAGs once we generate them, including calculating their quality, classifying them taxonomically, annotating the genes in them, identifying redundant MAGs, and calculating "abundance" or coverage of MAGs in different metagenome samples.

## Generate a non-redundant MAG dataset

I generally prefer to use an average nucleotide identity cutoff of 98% (see InStrain documentation for some discussion on cutoffs)

```bash
conda activate /PATH/TO/CONDAENV/drep

dRep dereplicate drep98_MAGs_ALL -g MAGs_all/*fa --ignoreGenomeQuality -sa 0.98
```

## Taxonomic classification

We use the GTDB-tk and install it using conda.


### MAG quality

Since we have refined the MAGs we have also generally renamed and consolidated the MAGs into 1 directory so the metawrap stats files are not always usable. Generally we use checkM to calculate the completeness, contamination, GC, length, etc. I have used checkM directly or also used metawrap because I like the output files (takes WAY longer but again I like the output format)

checkM example

metawrap example



### Lineage abundance

We generally refer to the representative non-redundant MAGs in our dataset as "lineages" but I have seen other terms in the literature such as species, operational taxonomic units (OTUs), 


##############



conda activate /PATH/TO/CONDAENV/gtdbtk-2.4.0

cat > GTDB.sh

#!/bin/bash
#################
#SBATCH --job-name=GTDB2
#SBATCH --out=logs/sl-gtdb-%j.out
#SBATCH --partition=PARTITION
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=8gb
##SBATCH --time=24:00:00
#SBATCH --time=4:00:00
#################
source /PATH/TO/USER/.bashrc
conda activate /PATH/TO/CONDAENV/gtdbtk-2.4.0
mags_dir='/PATH/TO/SAMPLE/MAGs_coas_final_renamed'
gtdbtk classify_wf --cpus 16 --skip_ani_screen --genome_dir $mags_dir -x fa --out_dir gtdb_MAGs_coas_final_renamed

mags_dir='/PATH/TO/SAMPLE/ERM2_2017Sep27_15-40/BIN_REFINEMENT_FINAL/metawrap_50_10_bins'
gtdbtk classify_wf --cpus 16 --skip_ani_screen --genome_dir $mags_dir -x fa --out_dir gtdb_22


cat > checkM.sh
#!/bin/bash 
################# 
#SBATCH --job-name=checkM
#SBATCH --out=logs/sl-checkM—%j.out
#SBATCH --partition=PARTITION
#SBATCH -n 8
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=16gb
#SBATCH --time=4-00:00:00
################# 
source /PATH/TO/USER/.bashrc
conda activate /PATH/TO/CONDAENV/metawrap-1.3.2

bin_dir='/PATH/TO/SAMPLE/MAGs_coas_final_renamed' 
metawrap bin_refinement -o checkM_MAGs_coas_final_renamed -t 8 -m 256 -A $bin_dir -c 50 -x 10 

bin_dir='/PATH/TO/SAMPLE/ERM2_2017Sep27_15-40/BIN_REFINEMENT_FINAL/metawrap_50_10_bins' 
metawrap bin_refinement -o checkM_22 -t 8 -m 256 -A $bin_dir -c 50 -x 10 


##########################################################
Getting all the MAGs organized for read mapping
##########################################################
 
#1 move the dereplicated MAGs into one folder
################################
mkdir MAGs_ER171819_drep98_contigsrenamed
scp drep98_MAGs_ALL_renamed/dereplicated_genomes/*fa MAGs_ER171819_drep98_contigsrenamed

#count number of files in a directory
ls -1|wc -l

#2 add the file name to the contig names
################################
cd MAGs_ER171819_drep98_contigsrenamed/
for f in *.fa; do sed -i "s/^>/>${f%.*}_/" "$f"; done
# the %.* after the f removes the file extension so names look nicer

#3 make a single fasta file
################################
cat MAGs_ER171819_drep98_contigsrenamed/*fa > MAGs_ER171819_drep98_contigsrenamed.fasta

#4 get scaffold to bin file for mapping
################################
# Next we need to create a scaffold-to-bin file. This can easily be done using the 
# script parse_stb.py that comes with the program dRep:

conda activate /PATH/TO/CONDAENV/drep-3.2.2

parse_stb.py --reverse -f MAGs_ER171819_drep98_contigsrenamed/* -o MAGs_ER171819_drep98_contigsrenamed_scaffold_to_bin_file.stb


#5 call genes (can be done on bins or combined file)
################################
conda activate /PATH/TO/CONDAENV/metawrap-1.3.2

cat > prodigal.sh

#!/bin/bash 
################# 
#SBATCH --job-name=prodigal
#SBATCH --out=logs/sl-prodigal-%j.out
#SBATCH --partition=PARTITION
#SBATCH --mem=8G
#SBATCH -n 16
#SBATCH --time=04:00:00
################# 

prodigal -i MAGs_ER171819_drep98_contigsrenamed.fasta -o MAGs_ER171819_drep98_contigsrenamed.gene.coord.gff -a MAGs_ER171819_drep98_contigsrenamed.proteins.faa -d MAGs_ER171819_drep98_contigsrenamed.gene.fna -m -p single -f gff

scp -r arasmuss@login.sherlock.stanford.edu:/PATH/TO/SAMPLE/MAGs_ER171819_drep98_contigsrenamed.proteins.faa /Users/annarasmussen/Documents/EastRiver/coassemblies_ANR_floodplainER

#6 bowtie2 to recruit reads!
################################

conda activate /PATH/TO/CONDAENV/metawrap-1.3.2


conda activate /PATH/TO/CONDAENV/metawrap-1.3.2

cat > bowtie2build.sh
#!/bin/bash 
################# 
#SBATCH --job-name=bowtiebuild
#SBATCH --out=logs/sl-bowtie2-build-%j.out
#SBATCH --partition=PARTITION
#SBATCH --mem=72G
#SBATCH -n 24
#SBATCH --time=12:00:00
################# 
#### Mapping to genome database
bowtie2-build MAGs_ER171819_drep98_contigsrenamed.fasta MAGs_ER171819_drep98_contigsrenamed.fasta.bt2 --large-index --threads 24 
echo done


################################
SLURM ARRAY
################################

cat > bowtie2_slurm.sh 
#!/bin/bash

###############################
#                             #
#  1) Job Submission Options  #
#                             #
###############################

# Name
#SBATCH --job-name=bowtie

# Resources
# For MPI, increase ntasks-per-node
# For multithreading, increase cpus-per-task
#SBATCH --nodes=1
#SBATCH -n 24
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=1gb
#SBATCH --time=0-06:00:00 #will use much less time probably
#SBATCH --partition=PARTITION

# Logs
#SBATCH --output=logs/sl-bowtie2-13-%j.out
#SBATCH --mail-user=EMAIL@EMAIL.EMAIL
#SBATCH --mail-type=BEGIN,END,FAIL

# Environment
##SBATCH --export=ALL
#SBATCH --array=0-23

source /PATH/TO/USER/.bashrc
conda activate /PATH/TO/CONDAENV/metawrap-1.3.2

#  Show list of CPUs you ran on
echo $SLURM_JOB_NODELIST
echo ${SLURM_ARRAY_TASK_ID}

echo start

bash /PATH/TO/SAMPLE/bowtie2.sh ${SLURM_ARRAY_TASK_ID}

echo done


###

sbatch --dependency=afterany:65546339 bowtie2_slurm.sh 

################################
SCRIPT
################################

cat > bowtie2.sh 
arrayid=$(($1 + 1))

sample=$(sed -n "${arrayid}p" /PATH/TO/SAMPLE/sample_array_all.txt)

echo $sample

bowtie2 -q -p 24 -x MAGs_ER171819_drep98_contigsrenamed.fasta.bt2 -1 fastq_coassembled/${sample}_1.fastq -2 fastq_coassembled/${sample}_2.fastq -S MAGs_ER171819_drep98_contigsrenamed_${sample}.sam

echo finished

####### OUTPUT for bowtie read mapping #######
cd logs
cat sl-bowtie2-13* > bowtie2_ER171819_output.txt
cat sl-bowtie2-13-655* > bowtie2_ER171819_output.txt  
# How to print just the lines I need
awk 'NR%21==4||NR%21==5' bowtie2_output_ER171819_drep98.txt > bowtie2_reads.txt
awk 'NR%23==6||NR%23==7||NR%23==21' bowtie2_output_ER171819_drep98.txt > bowtie2_reads.txt

#replace all \n with \t then replace sentence after read count with \n to get in nice format

awk 'NR % 5 == 0' input > output
NUM=5
awk -v NUM=$NUM 'NR % NUM == 0' input > output

sed -n '0~5p' oldfile > newfile

###################################################################
7. MAG abundance
###################################################################


###############################################
#### Move files to personal computer for R ####
###############################################

scp -r arasmuss@login.sherlock.stanford.edu:/PATH/TO/SAMPLE/Riverton_refined_bins_manual/Instrain/MAGs_ER171819_drep98_contigsrenamed_*.IS/output/*IS_genome_info.tsv /Users/annarasmussen/Documents/Riverton/data/Instrain/
scp -r arasmuss@login.sherlock.stanford.edu:/PATH/TO/SAMPLE/Riverton_refined_bins_manual/Instrain/MAGs_ER171819_drep98_contigsrenamed_*.IS/output/*RVT2*IS_genome_info.tsv /Users/annarasmussen/Documents/Riverton/data/Instrain/
scp -r arasmuss@login.sherlock.stanford.edu:/PATH/TO/SAMPLE/Riverton_refined_bins_manual/Instrain/MAGs_ER171819_drep98_contigsrenamed_*.IS/output/*RVTP2*IS_genome_info.tsv /Users/annarasmussen/Documents/Riverton/data/Instrain/



################
#### Redo coverM for 3 samples that have EOF errors ####
################

## SAMPLES TO RERUN
5
13
15

##example output file for one with error
start
PTT1_JUN19_120cm
[2024-09-09T19:33:27Z INFO  coverm] CoverM version 0.4.0
[2024-09-09T19:33:27Z INFO  coverm] Using min-covered-fraction 0%
[W::bam_hdr_read] EOF marker is absent. The input is probably truncated
[2024-09-09T19:36:29Z INFO  coverm::genome] Of 1868812 reference IDs, 1868812 were assigned to a genome and 0 were not
[2024-09-09T19:36:30Z INFO  coverm::genome] In sample 'MAGs_ER171819_drep98_contigsrenamed_PTT1_JUN19_120cm.sorted', found 1808575 reads mapped out of 1882078 total (96.09%)
finished
done

##example of what inStrain runs
samtools view -S -@ 8 -b Instrain/MAGs_ER171819_drep98_contigsrenamed_PTT1_OCT19_120cm.sam > Instrain/MAGs_ER171819_drep98_contigsrenamed_PTT1_OCT19_120cm.bam
samtools sort Instrain/MAGs_ER171819_drep98_contigsrenamed_PTT1_OCT19_120cm.bam -o Instrain/MAGs_ER171819_drep98_contigsrenamed_PTT1_OCT19_120cm.sorted.bam -@ 8


#### Rerun samtools to get sorted bam files ####

cat > samtools_slurm.sh 
#!/bin/bash

###############################
#                             #
#  1) Job Submission Options  #
#                             #
###############################

# Name
#SBATCH --job-name=samtools

# Resources
# For MPI, increase ntasks-per-node
# For multithreading, increase cpus-per-task
#SBATCH --nodes=1
#SBATCH -n 12
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2gb
#SBATCH --time=0-1:00:00
#SBATCH --partition=PARTITION

# Logs
#SBATCH --output=logs/sl-samtools-13-%j.out
#SBATCH --mail-user=EMAIL@EMAIL.EMAIL
#SBATCH --mail-type=BEGIN,END,FAIL

# Environment
##SBATCH --export=ALL

##SBATCH --array=0-67
#SBATCH --array=0-23
source /PATH/TO/USER/.bashrc
conda activate /PATH/TO/CONDAENV/instrain-env
ml python/3.6.1

#  Show list of CPUs you ran on
echo $SLURM_JOB_NODELIST
echo ${SLURM_ARRAY_TASK_ID}

echo start

bash /PATH/TO/SAMPLE/samtools.sh ${SLURM_ARRAY_TASK_ID}

echo done



cat > samtools.sh
arrayid=$(($1 + 1))

sample=$(sed -n "${arrayid}p" /PATH/TO/SAMPLE/sample_array_all.txt)

echo $sample

##samtools view -S -@ 12 -b Instrain/MAGs_ER171819_drep98_contigsrenamed_${sample}.sam > Instrain/MAGs_ER171819_drep98_contigsrenamed_${sample}.bam
##samtools sort Instrain/MAGs_ER171819_drep98_contigsrenamed_${sample}.bam -o Instrain/MAGs_ER171819_drep98_contigsrenamed_${sample}.sorted.bam -@ 12
samtools view -S -@ 12 -b MAGs_ER171819_drep98_contigsrenamed_${sample}.sam > MAGs_ER171819_drep98_contigsrenamed_${sample}.bam
samtools sort MAGs_ER171819_drep98_contigsrenamed_${sample}.bam -o MAGs_ER171819_drep98_contigsrenamed_${sample}.sorted.bam -@ 12

echo finished


sbatch --dependency=aftercorr:65546396 samtools_slurm.sh

scp -r arasmuss@login.sherlock.stanford.edu:/PATH/TO/SAMPLE/*coverm_output.txt /Users/annarasmussen/Documents/EastRiver/coassemblies_ANR_floodplainER/coverM/


#2 coverM 
################################
# works much faster than inStrain for getting abundance

cat > coverm_slurm.sh 
#!/bin/bash

###############################
#                             #
#  1) Job Submission Options  #
#                             #
###############################

# Name
#SBATCH --job-name=coverM

# Resources
# For MPI, increase ntasks-per-node
# For multithreading, increase cpus-per-task
#SBATCH --nodes=1
#SBATCH -n 12
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2gb
#SBATCH --time=0-00:15:00
#SBATCH --partition=PARTITION

# Logs
#SBATCH --output=logs/sl-coverM-13-%j.out
#SBATCH --mail-user=EMAIL@EMAIL.EMAIL
#SBATCH --mail-type=BEGIN,END,FAIL

# Environment
##SBATCH --export=ALL


#SBATCH --array=0-23
source /PATH/TO/USER/.bashrc
conda activate /PATH/TO/CONDAENV/coverm-0.4.0
ml python/3.6.1

#  Show list of CPUs you ran on
echo $SLURM_JOB_NODELIST
echo ${SLURM_ARRAY_TASK_ID}

echo start

bash /PATH/TO/SAMPLE/coverm.sh ${SLURM_ARRAY_TASK_ID}

echo done




cat > coverm.sh
arrayid=$(($1 + 1))

sample=$(sed -n "${arrayid}p" /PATH/TO/SAMPLE/sample_array_all.txt)

echo $sample

coverm genome --bam-files MAGs_ER171819_drep98_contigsrenamed_${sample}.sorted.bam --genome-fasta-directory MAGs_ER171819_drep98_contigsrenamed --genome-fasta-extension fa -m count mean covered_fraction length -t 12 --min-covered-fraction 0 > ${sample}_coverm_output.txt

echo finished

sbatch --dependency=aftercorr:65546443 coverm_slurm.sh


###################################################################
8. Move specific MAGs to a folder using while and txt document with the names of MAGs you want
###################################################################

while IFS= read -r file
do
    scp -i -- /PATH/TO/SAMPLE/Riverton_refined_bins_manual/MAGs_drep98/dereplicated_genomes/"$file".fa /PATH/TO/SAMPLE/Riverton_refined_bins_manual/Nitrifier_MAGs/
done </PATH/TO/SAMPLE/Riverton_refined_bins_manual/Nitrifier_MAGs.txt


while IFS= read -r file
do
    scp -i -- /PATH/TO/SAMPLE/Riverton_refined_bins_manual/MAGs_2015_ALL_refined_manual/"$file".fa /PATH/TO/SAMPLE/Riverton_refined_bins_manual/Nitrifier_MAGs_2015/
done </PATH/TO/SAMPLE/Riverton_refined_bins_manual/Nitrifier_MAGs_2015.txt


while IFS= read -r file
do
    scp -i -- /PATH/TO/SAMPLE/Riverton_refined_bins_manual/MAGs_2019_ALL_refined_manual/"$file".fa /PATH/TO/SAMPLE/Riverton_refined_bins_manual/Nitrifier_MAGs_2019/
done </PATH/TO/SAMPLE/Riverton_refined_bins_manual/Nitrifier_MAGs_2019.txt

while IFS= read -r file
do
    scp -i -- /PATH/TO/SAMPLE/Riverton_initial_binning/MAGs_2017_ALL_refined/"$file".fa /PATH/TO/SAMPLE/Riverton_initial_binning/Nitrifier_MAGs_2017/
done </PATH/TO/SAMPLE/Riverton_initial_binning/Nitrifier_MAGs_2017.txt



