---
output:
  html_document: default
  pdf_document: default
---
# Metagenome assembly and binning

This is an example of a metagenome assembly and binning pipeline through to refinement of metagenome-assembled genomes (MAGs) of medium to high quality. This pipeline uses primarily *bash* scripts and my favorite tool [MetaWRAP](https://github.com/bxlab/metaWRAP/tree/master) which has a great Github page and helpful tutorial. The code used here is based on the aforementioned [MetaWRAP tutorial](https://github.com/bxlab/metaWRAP/blob/master/Usage_tutorial.md) and edited to be used on a computing cluster with a [slurm workload manager](https://slurm.schedmd.com/documentation.html). Some data wrangling also uses *R*.

Examples of publications using this (or some iteration of this) pipeline include:

+[Rasmussen, Anna N., and Christopher A. Francis. "Genome-resolved metagenomic insights into massive seasonal ammonia-oxidizing archaea blooms in San Francisco Bay." Msystems 7.1 (2022): e01270-21.](https://doi-org.stanford.idm.oclc.org/10.1128/msystems.01270-21)
+[Rasmussen, Anna N., and Christopher A. Francis. "Dynamics and activity of an ammonia-oxidizing archaea bloom in South San Francisco Bay." The ISME Journal 18.1 (2024): wrae148.]( https://doi-org.stanford.idm.oclc.org/10.1093/ismejo/wrae148)
+[Rasmussen, Anna N., et al. "Diverse and unconventional methanogens, methanotrophs, and methylotrophs in metagenome-assembled genomes from subsurface sediments of the Slate River floodplain, Crested Butte, CO, USA." Msystems 9.7 (2024): e00314-24.](https://doi-org.stanford.idm.oclc.org/10.1128/msystems.00314-24)
+[Rasmussen, Anna N., et al. "Metagenome‐Assembled Genomes for Oligotrophic Nitrifiers From a Mountainous Gravelbed Floodplain." Environmental Microbiology 27.3 (2025): e70060.](https://doi-org.stanford.idm.oclc.org/10.1111/1462-2920.70060)

## Step 1: Prepare metagenomic data

I generally work in the university computing cluster and have a site-specific or sampling campaign-specific directory that I work in.


Make a project directory.

```bash
mkdir MY_SAMPLING_SITE
```

Make a list of sample names for use in arrays.

```bash
cat > sample_array.txt
SAMPLE_1
SAMPLE_2
SAMPLE_3

#Use Ctrl + d to exit
```

#### *1.1: Download metagenomes*

I generally analyze metagenomes sequenced through JGI. I download the filtered metagenome reads that JGI has removed contamination from instead of 'Raw Data' metagenome files. I use Globus to download the metagenomes.I am usually working with many metagenomes and download the entire "Filtered_Raw_Date" directory for each metagnome. This leads to a separate folder for each metagenome containing several files, the metagenome fastq files generally look something like:

```bash
52690.2.420404.TTGCGAAG-TTGCGAAG.filter-METAGENOME.fastq.gz #example file name straight from JGI
```

I often rename the files to match the sample names instead of the long JGI assigned name which is a manual process. If I have my wits about me, I will first move all the fastq files into one directory together.

```bash
mkdir fastq
mv */*fastq.gz fastq/ #move the fastq.gz files from all the directories into one directory

cd fastq #move into the fastq directory
mv 52690.2.420404.TTGCGAAG-TTGCGAAG.filter-METAGENOME.fastq SAMPLE_1.fastq.gz #start renaming files using mv
```

If you come up with a better way to rename files other than mv, go for it!


#### *1.2: Separate metagenome files into forward and reverse reads*

Example of how we separate the forward and reverse reads from the JGI-generated metagenomes.

```bash
for g in *.gz; do gunzip $g; done #unzip files if need be

for i in *.fastq; do paste - - - - < "$i" | tee >(awk 'BEGIN{FS="\t"; OFS="\n"} {if (match($1, " 1:N")) print $1,$2,$3,$4}' > "$i"_1.fastq ) | awk 'BEGIN{FS="\t"; OFS="\n"} {if (match($1, " 2:N")) print $1,$2,$3,$4}' > "$i"_2.fastq; done
```

#### Side note: Renaming files

Sometimes I haven't renamed the JGI files and need to fix ugly files names. Always good to plan ahead but just in case, here is an example using rename

```bash
#example ugly file name 
rename .filter-METAGENOME.fastq_ _ *fastq #replace .filter-METAGENOME.fastq_ with just a single _
```

or R

```R
R

filez = list.files(pattern = ".fastq")
head(filez)

for (file in filez) {
rm = strsplit(file, '[_.]')[[1]] #split the file name using _ or . separaters
file.rename(file, paste(rm[1], "_",rm[2],"_",rm[3],"_",rm[6], ".fastq", sep= "")) #select the important parts of the file name to keep
}
```


## Step 2: Assemble and/or co-assemble

I typically assemble individual samples. However, depending on how successfully/deeply samples were sequenced, if we sequenced the same sample multiple times, or if we are targeting low-abundance organisms I will coassemble multiple samples together. See coassembly subsection.

#### *2.1: Prepare more directories to work in*

I like to make a directory for each sample.

```
while IFS= read -r dir; do mkdir -p "$dir"; done <sample_array.txt
```

#### *2.2: Assemble*
If you aren't using arrays just follow the [MetaWRAP tutorial](https://github.com/bxlab/metaWRAP/blob/master/Usage_tutorial.md)- it is excellent! Your assembly command will look something like:

```bash
metawrap assembly -1 /PATH/TO/MY_SAMPLING_SITE/fastq/SAMPLE_1_1.fastq -2 /PATH/TO/MY_SAMPLING_SITE/fastq/SAMPLE_1_2.fastq -m 384 -t 8 --megahit -o SAMPLE_1/assembly_SAMPLE_1
```

If you use a slurm based job manager and want to use arrays, this is an example of what my slurm submission script looks like:

```bash
cat > assembly_slurm.sh
#!/bin/bash

###############################
#                             #
#  1) Job Submission Options  #
#                             #
###############################

# Name
#SBATCH --job-name=mh_assembly

# Resources
#SBATCH -n 8
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=48gb
#SBATCH --time=1-00:00:00
#SBATCH --partition=PARTITION_NAME

# Logs
#SBATCH --output=logs/sl-assembly-%j.out ##make sure to make a log file
#SBATCH --mail-user=EMAIL@EMAIL.EMAIL ##use your email not mine ;)
#SBATCH --mail-type=END,FAIL ##set only to notify if ends or fails, can add BEGIN,

# Environment
#SBATCH --array=0-2
source /PATH/TO/USER/.bashrc
conda activate /PATH/TO/CONDAENV/metawrap-1.3.2

#  Show list of CPUs you ran on
echo $SLURM_JOB_NODELIST
echo ${SLURM_ARRAY_TASK_ID}

echo start

bash /PATH/TO/MY_SAMPLING_SITE/assembly.sh ${SLURM_ARRAY_TASK_ID}

echo done

```

And here is the actual script for running the metawrap assembly pipeline

```bash
cat > assembly.sh 
arrayid=$(($1 + 1))

sample=$(sed -n "${arrayid}p" /PATH/TO/MY_SAMPLING_SITE/sample_array.txt)

echo $sample

metawrap assembly -1 /PATH/TO/MY_SAMPLING_SITE/fastq/${sample}_1.fastq -2 /PATH/TO/MY_SAMPLING_SITE/fastq/${sample}_2.fastq -m 384 -t 8 --megahit -o ${sample}/assembly_${sample}

echo $sample
echo finished
```


### *Coassembly*

We tend to know enough about our samples to decide which samples are appropriate to coassemble. However, *sourmash* could be used if you need more input deciding which samples have similar enough communities to reasonably co-assemble.

First, I concatenate the appropriate *.fastq* files. I often do this by making symbolic links of the *.fastq* files I want to concatenate in a directory that I name using the same new concatenated file name I want to use. I also make a new sample array with the new names and edit my scripts as appropriate.

Then we coassemble which uses more compute. If *.fastq* files are too large, sometimes co-assemblies are not possible with our computing resources.


## Step 3: Bin

I bin using all 3 binners available in *metawrap binning* and generally use multiple fastq files to calculate coverage information. I have done both single sample binning and binning with multiple *.fastq* files and at times have done both.

Make symbolic links of appropriate *.fastq* files in every sample directory or main directory

```bash
for dir in SAMPLE*/; do mkdir -- "$dir/fastq"; done
for dir in SAMPLE*/; do ln -s /PATH/TO/MY_SAMPLING_SITE/fastq/ -- "$dir/fastq"; done
```

Make a binning slurm submission script similar to the *assembly_slurm.sh* above but pointing to the binning script like the one below:

```bash
cat > init_binning.sh 
arrayid=$(($1 + 1))

sample=$(sed -n "${arrayid}p" /PATH/TO/MY_SAMPLING_SITE/sample_array.txt)

metawrap binning -o ${sample}/INIT_BINNING -t 16 -m 128 -l 2000 --universal -a /PATH/TO/MY_SAMPLING_SITE/${sample}/assembly_${sample}/final_assembly.fasta --metabat2 --maxbin2 --concoct /PATH/TO/MY_SAMPLING_SITE/${sample}/fastq/*_1.fastq /PATH/TO/MY_SAMPLING_SITE/${sample}/fastq/*_2.fastq

echo $sample
echo finished

```

If you have a slurm workload manager you can also use dependencies to submit scripts to run as soon as your previous step is finished.

```bash
sbatch --dependency=aftercorr:65324457 init_binning_slurm.sh #number is job number for assembly script
```

### *Coassembly binning*
I usually use all of the individual sample *.fastq* files to calculate coverage for binning rather than the concatenated fastq file.

## Step 4: Refine bins

In the step, we refine the bins to a finalized MAG dataset. This involves using the metawrap *metawrap bin_refinement* tool to combine redundant bins across the 3 binning algorithms (*Maxbin2*, *Metabat2*, *CONCOCT*) and filtering out bins based on *checkM*-calculated completeness and contamination thresholds. We generally choose < 10 % contamination and > 50% completeness.

```bash
cat > bin_refinement.sh 
arrayid=$(($1 + 1))

sample=$(sed -n "${arrayid}p" /PATH/TO/MY_SAMPLING_SITE/sample_array.txt)

echo $sample

metawrap bin_refinement -o ${sample}/BIN_REFINEMENT -t 24 -m 72 -A /PATH/TO/MY_SAMPLING_SITE/${sample}/INIT_BINNING/metabat2_bins -B /PATH/TO/MY_SAMPLING_SITE/${sample}/INIT_BINNING/maxbin2_bins -C /PATH/TO/MY_SAMPLING_SITE/${sample}/INIT_BINNING/concoct_bins -c 50 -x 10 

echo $sample
echo finished
```

### *Refining separate binning efforts for the same sample*

I often try multiple binning efforts for the same sample depending on what we know about the sample or how well our first binning attempts went. For example, we have tried binning using just the sample fastq files and then decided to bin using multiple *.fastq* files or with different minimum contig lengths. Luckily, the *metabat bin_refinement* tool can be used to refine all the different bins for the same sample. 



```bash
cat > bin_refinement_second_round.sh 
arrayid=$(($1 + 1))

sample=$(sed -n "${arrayid}p" /PATH/TO/MY_SAMPLING_SITE/sample_array.txt)

echo $sample

metawrap bin_refinement -o ${sample}/BIN_REFINEMENT_FINAL -t 24 -m 96 -A /PATH/TO/MY_SAMPLING_SITE/${sample}/BIN_REFINEMENT/metawrap_50_10_bins -B /PATH/TO/MY_SAMPLING_SITE/${sample}/BIN_REFINEMENT_2/metawrap_50_10_bins -c 50 -x 10 

echo $sample
echo finished
```

Rename the refined bins to include the sample name instead of just bin.1.fa to SAMPLE_1_bin_1.fa. This is a bit of a painful manual way depending on how many samples you have.

```bash
cd SAMPLE_1/BIN_REFINEMENT/metawrap_50_10_bins/

ls -1|wc -l #get count of bins in the directory if that's of interest to you
 
rename bin. SAMPLE_1_ *
```

Make a copy of all the MAGs in one folder. 

```bash
mkdir MAGs_all
scp SAMPLE*/BIN_REFINEMENT/metawrap_50_10_bins/*fa MAGs_all/
```

Once I have all of the MAGs from all of the metagenomes together in one directory, I move onto analyzing the MAGs.
