# Processing MAGs

### Overview
Once I have generated a bunch of MAGs from metagenomes, I start generating metadata for the MAGs, including: genome quality, taxonomic classification, gene calls and annotation, abundance/coverage in different metagenome samples, and identifying representative genomes.

I try to do all data wrangling in R and minimally process output files from each program.

See Metagenome_assembly_and_binning.md for more details, but generally at this stage I have a "My_SAMPLING_SITE" directory that has the following:
+ fastq (directory with metagenome forward and reverse fastq files)
+ MAGs_all (directory of all the MAGs from all the samples, renamed with sample included)
+ sample_array.txt (list of sample names)
+ logs (directory with the outputs from all the slurm submissions)
+ directories corresponding to each sample (have all of the assemblies, initial binning, and bin refinement directories)

## Generate a non-redundant MAG dataset

I generally prefer to use an average nucleotide identity cutoff of 98%. The [inStrain](https://instrain.readthedocs.io/en/latest/important_concepts.html#picking-and-evaluating-representative-genomes) documentation has some very helpful discussion on cutoffs and creating a non-redundant dataset.

I use [dRep](https://github.com/MrOlm/drep) to pick representative genomes, which I usually refer to as lineages, but can be called species, species representatives, OTUs, etc. in the literature.

```bash
conda activate /PATH/TO/CONDAENV/drep #I use conda to install programs

dRep dereplicate drep98_MAGs_ALL -g MAGs_all/*fa --ignoreGenomeQuality -sa 0.98
```

I then move the representative MAGs into one directory and rename the contig names to include the MAG names. I highly recommend keeping both the original and the renamed-contig versions of MAGs, NCBI and some other programs do not like long contig names.

```bash
mkdir MAGs_all_drep98
scp drep98_MAGs_ALL/dereplicated_genomes/*fa MAGs_all_drep98

cd MAGs_all_drep98
ls -1|wc -l #count number of files in a directory
for f in *.fa; do sed -i "s/^>/>${f%.*}_/" "$f"; done #add the bin name to the contig name without the .fa

cd ..
cat MAGs_all_drep98/*fa > MAGs_all_drep98.fasta #make a concatenated fasta of all the MAG sequences
```

## Taxonomic classification

I use the GTDB-tk to classify MAGs. I generally classify all the MAGs I have generated because sometimes the representative MAG selected by drep is missing the functional genes I am most interested in and I want to dig into all the MAGs from a certain Genus or something.

Here is an example with the slurm commands included.

```bash
cat > GTDB.sh

#!/bin/bash
#################
#SBATCH --job-name=GTDB
#SBATCH --out=logs/sl-gtdb-%j.out
#SBATCH --partition=PARTITION
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=8gb
##SBATCH --time=24:00:00
#SBATCH --time=4:00:00
#################
source /PATH/TO/USER/.bashrc
conda activate /PATH/TO/CONDAENV/gtdbtk

mags_dir='/PATH/TO/MY_SAMPLING_SITE/MAGs_all'
gtdbtk classify_wf --cpus 16 --skip_ani_screen --genome_dir $mags_dir -x fa --out_dir gtdb_MAGs_all

```

I then concatenate the 2 output files (1 for bacteria and archaea). Then in R, I do some data wrangling like separate the taxonomy columns. 

```R
gtdb <- read.csv("gtdb_output.txt", sep = "\t")
gtdb <- gtdb %>%
  mutate(classification2 = classification) %>% #make extra classification column to split up
  separate(., classification2, c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";") 
```

### MAG quality

Since I have renamed and consolidated all of the MAGs into 1 directory, the lovely *metaWRAP* stats from the bin refinement step are more difficult to use for quickly getting quality information such as completeness, contamination, GC, length. I have taken the stat files from each sample and renamed the bins in the txt file before then concatenating them all but if you have a lot of samples that can be a lot to manage. Generally, I use [CheckM](https://ecogenomics.github.io/CheckM/) again on the consolidated MAG set for simplicity. I have used *checkM* directly or also used *metaWRAP* because I like the output files (takes WAY longer than just straight checkM but I like the output format)

Here is a *checkM* example:

```bash
conda activate /PATH/TO/CONDAENV/metawrap #checkM is installed with this program or you could install checkM on its own

bin_dir='/PATH/TO/MY_SAMPLING_SITE/MAGs_all' 
checkm lineage_wf -t 24 -x fa --pplacer_threads 2 --file checkM_output.txt --tab_table $bin_dir checkM_MAGs_all
echo done
```bash


Here is a *metaWRAP* example:

```bash
conda activate /PATH/TO/CONDAENV/metawrap

bin_dir='/PATH/TO/MY_SAMPLING_SITE/MAGs_all' 
metawrap bin_refinement -o checkM_MAGs_all -t 8 -m 256 -A $bin_dir -c 50 -x 10 
```

### Lineage abundance

Once I have a dereplicated or non-redundant MAG dataset, I like to calculate the abundance and coverage of MAGs in the metagenome samples I am working with. I have used both [CoverM](https://github.com/wwood/CoverM) and [InStrain](https://instrain.readthedocs.io/en/latest/index.html) at this step. There are many different tools/programs involved in this process, but generally the first step is to recruit the metagenome reads back to the non-redundant MAG dataset.

Using [bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml), I first map the reads back to the genome database (the concatenated fasta file of all the dereplicated MAGs). 

```bash
conda activate /PATH/TO/CONDAENV/metawrap #or a bowtie2 specific install

bowtie2-build MAGs_all_drep98.fasta MAGs_all_drep98.fasta.bt2 --large-index --threads 24 
echo done
```

Then I set up an array so I can run all of the read recruitment for all of the metagenomes as shown below:

```bash
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
#SBATCH --time=0-03:00:00
#SBATCH --partition=PARTITION

# Logs
#SBATCH --output=logs/sl-bowtie2-%j.out
#SBATCH --mail-user=EMAIL@EMAIL.EMAIL
#SBATCH --mail-type=BEGIN,END,FAIL

# Environment
#SBATCH --array=0-2

source /PATH/TO/USER/.bashrc
conda activate /PATH/TO/CONDAENV/metawrap

#  Show list of CPUs you ran on
echo $SLURM_JOB_NODELIST
echo ${SLURM_ARRAY_TASK_ID}

echo start

bash /PATH/TO/MY_SAMPLING_SITE/bowtie2.sh ${SLURM_ARRAY_TASK_ID}

echo done
```

###
```bash
sbatch --dependency=afterany:65546339 bowtie2_slurm.sh #can submit using a dependency on the genome database step finishing
```bash

And here is the script that has the actual bowtie2 commands that the slurm script points to.

```bash
cat > bowtie2.sh 
arrayid=$(($1 + 1))

sample=$(sed -n "${arrayid}p" /PATH/TO/MY_SAMPLING_SITE/sample_array_all.txt)

echo $sample

bowtie2 -q -p 24 -x MAGs_all_drep98.fasta.bt2 -1 fastq/${sample}_1.fastq -2 fastq/${sample}_2.fastq -S MAGs_all_drep98_${sample}.sam

echo finished
```

Then I often like to take the output from the *bowtie2* run to get read totals and % of reads recruited to the non-redundant MAG dataset. This takes some data wrangling.

```bash
cd logs
cat sl-bowtie2-* > bowtie2_output_raw.txt #concaetnate bowtie2 outputs 
# How to print just the lines with sample name, read count, and % reads mapped
awk 'NR%21==4||NR%21==5||NR%21==19' bowtie2_output_raw.txt > bowtie2_output.txt

#in a text editor, replace all \n with \t then replace sentence after read count with \n to get in nice format

```bash

If I am using [CoverM](https://github.com/wwood/CoverM) to calculate coverage and metagenome reads recruited I then use [Samtools](https://www.htslib.org/) directly.

To save space I have omitted the slurm submission scripts and only included the Samtools commands below:

```bash
conda activate /PATH/TO/CONDAENV/samtools

samtools view -S -@ 12 -b MAGs_all_drep98_${sample}.sam > MAGs_all_drep98_${sample}.bam
samtools sort MAGs_all_drep98_${sample}.bam -o MAGs_all_drep98_${sample}.sorted.bam -@ 12

```

Then I will use the *CoverM* to calculate the reads recruited to each MAG. Again, for simplicity I am only including the *CoverM* commands but I use arrays to submit slurm scripts.


```bash
source /PATH/TO/USER/.bashrc
conda activate /PATH/TO/CONDAENV/coverm
ml python/3.6.1 #load a python module

coverm genome --bam-files MAGs_all_drep98_${sample}.sorted.bam --genome-fasta-directory MAGs_all_drep98 --genome-fasta-extension fa -m count mean covered_fraction length -t 12 --min-covered-fraction 0 > ${sample}_coverm_output.txt

```

Again, the dependency function is super helpful for submitting all these consecutive steps for each metagenome after the corresponding previous step finishes (bowtie2 to samtools to coverM)

```bash
sbatch --dependency=aftercorr:65546396 coverm_slurm.sh
```

I have also used [InStrain](https://instrain.readthedocs.io/en/latest/index.html) and think it is an awesome for digging into the microdiversity of strains across samples and even looking at SNPs and coverage at the gene level. Here is an *InStrain* example.

First we need to make a scaffold-to-bin (stb) file.

```bash
conda activate /PATH/TO/CONDAENV/drep

parse_stb.py --reverse -f MAGs_all_drep98/* -o MAGs_all_drep98_scaffold_to_bin_file.stb
```

To run inStrain I use the .sam files from each metagenomes (bowtie2 step), the concatenated fasta file of all the MAGs, the gene calls from the concatenated fasta file, and the stb file.

Again, I use a slurm array but have only included the inStrain commands below:

```bash
conda activate /PATH/TO/CONDAENVS/instrain-env
ml python/3.6.1

inStrain profile/PATH/TO/MY_SAMPLING_SITE/MAGs_all_drep98_${sample}.sam /PATH/TO/MY_SAMPLING_SITE/MAGs_all_drep98.fasta -o /PATH/TO/MY_SAMPLING_SITE/Instrain/MAGs_all_drep98_${sample}.IS -p 24 -g /PATH/TO/MY_SAMPLING_SITE/MAGs_all_drep98.gene.fna -s /PATH/TO/MY_SAMPLING_SITE/MAGs_all_drep98_scaffold_to_bin_file.stb --database_mode --skip_plot_generation

```

There are a lot of great outputs and the plots are helpful. Can be very memory and compute intensive.

There is also an inStrain compare function that I recommend looking into as a next step.


## Gene annotation

I have used a variety of different methods for gene annotation. Generally, my first step is making the gene calls, using [Prodigal](https://github.com/hyattpd/Prodigal) (which is conveniently installed with MetaWRAP). You could also call a prodigal environemnt if you have installed it separately.

```bash
conda activate /PATH/TO/CONDAENV/metawrap

prodigal -i MAGs_all_drep98.fasta -o MAGs_all_drep98.gene.coord.gff -a MAGs_all_drep98.proteins.faa -d MAGs_all_drep98.gene.fna -m -p single -f gff

```

Then I annotate the MAG gene calls using [kofam_scan](https://www.genome.jp/ftp/tools/kofam_scan/). Strongly recommend reading the documentation as it requires ruby, hmmsearch, and parallel. The installs of those programs are specified in the config file.

```bash
cat > config.yml

profile: /PATH/TO/CONDAENVS/kofam_scan-1.3.0/db/profiles
ko_list: /PATH/TO/CONDAENVS/kofam_scan-1.3.0/db/ko_list
hmmsearch: /PATH/TO/CONDAENVS/kofam_scan-1.3.0/bin/hmmsearch
parallel: /PATH/TO/CONDAENVS/kofam_scan-1.3.0/bin/parallel
cpu: 32
```

The actual script will then point to the kofam_scan install and the fasta file to be annotated. I chose the mapper format so I could then join the output file with a KO gene index but there are also some other programs. like [kofamparse](https://github.com/michauorin/kofamparse), for making the normal outputs more usable. 

```bash
ml ruby/3.1.2 #module load ruby

for genome in $(ls *.faa); do 
/PATH/TO/CONDAENVS/kofam_scan-1.3.0/kofam_scan-1.3.0/exec_annotation /PATH/TO/MY_SAMPLING_SITE/MAGs_all_drep98.proteins.faa -o /PATH/TO/MY_SAMPLING_SITE/kofamscan_out_MAGs_all_drep98.txt --format=mapper
done
```

### Next steps

Once I have all of the metadata for the MAGs (taxonomy, coverage, completeness, contamination, length, gene calls, gene annotations, etc) I get all of that data into R and start doing some data exploration and analysis. See the Community_analysis.md for how I use libraries like [phyloseq](https://joey711.github.io/phyloseq/) or the Genomic_analysis.md for how I use [anvi'o](https://anvio.org/). Both also have more examples of how I wrangle lots of data or visualize data in R but below are a few brief examples of how I get started.


#### *Wrangling all of the data in R* 

Here are just a few examples of how I start putting together all of the data. I am usually working on the computing cluster for all of the data processing steps and then move the output files to my personal computer to work on.

```R
#read in gtdb and checkM output files first
MAG.QC <- left_join(gtdb, checkM, by = "user_genome") #combine the gtdb and checkM outputs based on MAG names
```

And import all of the coverM data from the *coverm_output.txt* files. I am often working with a lot of metagenomes

```R
sample.array <- read.csv("sample_array.txt", sep = "\t", header = FALSE)[-c(1),]#text file with all metagenome names
sample.list <- sample_array #get a vector in case the input is a dataframe

for (i in (c(1:length(sample.list)))) {
   assign(sample.list[i], read.csv(paste("coverM/", sample.list[i],"_coverm_output.txt", sep = ""), sep = "\t") %>% mutate(Sample = sample.list[i]))
   assign(sample.list[i], setNames(get(sample.list[i]),  c("user_genome", "read.count", "mean", "covered.fraction", "length", "Sample")))
  }
  
library(data.table) #to use rbind list
coverM <- rbindlist(mget(sample.list)) #make a list of the samples to rbind!

coverM <- coverM %>%
  mutate(genome = str_split_i(user_genome, ".fa", 1)) #coverm keeps the file name so just need to remove the .fa to get the MAG name

```

Then I start joining all of the different datasets together

```R
coverM.dat <- coverM %>%
  left_join(metagenome.metadat , by = "Sample") %>% #file that has metagenome size in Gb and other metadata
  left_join(MAG.QC, by = "user_genome") %>% #add MAG taxonomy and quality data
  mutate(Kb = length/1000, #get genome size in Kilobases
         RPKG = (read.count/Kb)/Gb, #get reads recruited to kilobase of genome per gigabase of metagenome
         name = user_genome) 
         
coverM.dat.present <- coverM.dat %>%
  filter(covered.fraction >= 0.4) #only include MAGs with >= 40% genome coverage as "present" in a sample
```

#### *Moving specific MAGs* 

Once I have taxonomic information or other quality data, I often want to move a subset of MAGs (usually the nitrifiers) to a directory all on their own. I am a big fan of using *while*.

```R
#get the list of MAGs I want to move from R

MAG.QC %>%
filter(Order == "o__Nitrososphaerales") %>%
select(user_genome) %>%
write.csv(file = "Subset_of_MAGs_names.txt", quote = F, row.names = F)
```

Then in the terminal:

```bash
while IFS= read -r file
do
    scp -i -- /PATH/TO/SAMPLE/MAGs_all_drep98/"$file".fa /PATH/TO/SAMPLE/Subset_of_MAGs/
done </PATH/TO/SAMPLE/Subset_of_MAGs_names.txt
```


