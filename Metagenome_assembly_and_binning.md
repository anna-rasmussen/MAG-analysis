# Metagenome assembly and binning

### Overview
This is an example of a metagenome assembly and binning pipeline through to refinement of metagenome-assembled genomes (MAGs) of medium to high quality. This pipeline uses primarily *bash* scripts and my favorite tool [MetaWRAP](https://github.com/bxlab/metaWRAP/tree/master) which has a great Github page and helpful tutorial. The code used here is based on the aforementioned [MetaWRAP tutorial](https://github.com/bxlab/metaWRAP/blob/master/Usage_tutorial.md) and edited to be used on a computing cluster with a [slurm workload manager](https://slurm.schedmd.com/documentation.html). Some data wrangling also uses *R*.

Examples of publications using this (or some iteration of this) pipeline include:

+ [Rasmussen, Anna N., and Christopher A. Francis. "Genome-resolved metagenomic insights into massive seasonal ammonia-oxidizing archaea blooms in San Francisco Bay." Msystems 7.1 (2022): e01270-21.](https://doi-org.stanford.idm.oclc.org/10.1128/msystems.01270-21)
+ [Rasmussen, Anna N., and Christopher A. Francis. "Dynamics and activity of an ammonia-oxidizing archaea bloom in South San Francisco Bay." The ISME Journal 18.1 (2024): wrae148.]( https://doi-org.stanford.idm.oclc.org/10.1093/ismejo/wrae148)
+ [Rasmussen, Anna N., et al. "Diverse and unconventional methanogens, methanotrophs, and methylotrophs in metagenome-assembled genomes from subsurface sediments of the Slate River floodplain, Crested Butte, CO, USA." Msystems 9.7 (2024): e00314-24.](https://doi-org.stanford.idm.oclc.org/10.1128/msystems.00314-24)
+ [Rasmussen, Anna N., et al. "Metagenome‐Assembled Genomes for Oligotrophic Nitrifiers From a Mountainous Gravelbed Floodplain." Environmental Microbiology 27.3 (2025): e70060.](https://doi-org.stanford.idm.oclc.org/10.1111/1462-2920.70060)

## Step 1: Prepare metagenomic data

I generally work in the university computing cluster and have a site-specific or sampling campaign-specific directory that I work in.


Make a project directory.

```bash
mkdir MY_SAMPLING_SITE
cd MY_SAMPLING_SITE
```

Make some files and directories in your main directory such as a list of sample names for use in arrays and a logs directory.

```bash
cat > sample_array.txt
SAMPLE_1
SAMPLE_2
SAMPLE_3

#Use Ctrl + d to exit

mkdir logs
```

### 1.1: Download metagenomes

I generally analyze metagenomes sequenced through JGI. I download the filtered metagenome reads that JGI has removed contamination from instead of 'Raw Data' metagenome files. I use Globus to download the metagenomes. I usually work with many metagenomes and download the entire "Filtered_Raw_Date" directory for each metagnome that includes additional report files. This ulitmately leads to a separate folder for each metagenome containing several files, the metagenome fastq files generally are named something like:

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

If you have a better way to rename multiple files using something other than mv, go for it!

### 1.2: Separate metagenome files into forward and reverse reads

Example of how we separate the forward and reverse reads from the JGI-generated metagenomes.

```bash
for g in *.gz; do gunzip $g; done #unzip files if need be

for i in *.fastq; do paste - - - - < "$i" | tee >(awk 'BEGIN{FS="\t"; OFS="\n"} {if (match($1, " 1:N")) print $1,$2,$3,$4}' > "$i"_1.fastq ) | awk 'BEGIN{FS="\t"; OFS="\n"} {if (match($1, " 2:N")) print $1,$2,$3,$4}' > "$i"_2.fastq; done
```

### *Side note: Renaming files*

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

### 2.1: Prepare more directories to work in

I like to make a directory for each sample with the sample names I will use throughought the rest of the pipeline.

```
while IFS= read -r dir; do mkdir -p "$dir"; done <sample_array.txt
```

### 2.2: Assemble
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


### *Coassembling metagenomes*

We tend to know enough about our samples to decide which samples are appropriate to coassemble. However, [sourmash](https://sourmash.readthedocs.io/en/latest/) could be used if you need more input deciding which samples have similar enough communities to reasonably coassemble.

First, I concatenate the appropriate *.fastq* files. My general approach is to make a new set of directories that have the new, coassembled sample names and then move the *.fastq* files I want to concatenate in the appropriate new sample directory. I also make a new sample array with the new names and edit my scripts as appropriate. Once I have the concatenated *.fastq* files I will start my coassembly. Note, if *.fastq* files are too large, sometimes co-assemblies are not possible with our computing resources.

Below is some text with examples from moving *.fastq* files for samples sequenced in triplicate into new folders.

```bash
cat > sample_array_coassemblies.txt #new sample array
SAMPLE_co1
SAMPLE_co2


cat > forwardfastq.txt #make a guide for moving files where you want them, I usually make in R w/ metagenome metadata
file_name1;Sample_coassembly
52481.2.358586.ATGGAAGG-ATGGAAGG_1.fastq;SAMPLE_co1
52481.2.358586.TCAAGGAC-TCAAGGAC_1.fastq;SAMPLE_co1
52481.2.358586.GATTACCG-GATTACCG_1.fastq;SAMPLE_co1
52481.2.358586.GTCTGATC-GTCTGATC_1.fastq;SAMPLE_co2
52481.2.358586.AATACGCG-AATACGCG_1.fastq;SAMPLE_co2
52481.2.358586.CGACGTTA-CGACGTTA_1.fastq;SAMPLE_co2

cat > reversefastq.txt #make a guide for moving files where you want them, I usually make in R w/ metagenome metadata
file_name1;Sample_coassembly
52481.2.358586.ATGGAAGG-ATGGAAGG_2.fastq;SAMPLE_co1
52481.2.358586.TCAAGGAC-TCAAGGAC_2.fastq;SAMPLE_co1
52481.2.358586.GATTACCG-GATTACCG_2.fastq;SAMPLE_co1
52481.2.358586.GTCTGATC-GTCTGATC_2.fastq;SAMPLE_co2
52481.2.358586.AATACGCG-AATACGCG_2.fastq;SAMPLE_co2
52481.2.358586.CGACGTTA-CGACGTTA_2.fastq;SAMPLE_co2

#make the new directories
while IFS= read -r dir; do mkdir -p "$dir"; done <sample_array_coassemblies.txt

#move the forward reads into the appropriate sample directory
cat forwardfastq.txt |
while read line; do
  IFS=';'
  set - $line
  mv $1 $2/
done

#move the reverse reads into the appropriate sample directory
cat reversefastq.txt |
while read line; do
  IFS=';'
  set - $line
  mv $1 $2/
done

#concatenate the forward and reverse fastq files

cat Sample_co1/*_1.fastq > Sample_co1_1.fastq
cat Sample_co2/*_1.fastq > Sample_co2_1.fastq
cat Sample_co1/*_2.fastq > Sample_co1_2.fastq
cat Sample_co2/*_2.fastq > Sample_co2_2.fastq
```


## Step 3: Bin

I bin using all 3 binners available in *metawrap binning* and generally use multiple *.fastq* files to calculate coverage information. I have done both single sample binning and binning with multiple *.fastq* files and at times have done both.

You can set the path for the *.fastq* files directly to the fastq directory or make symbolic links in every sample directory to point to the fastq directory as shown below:

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
I usually use all of the individual sample *.fastq* files to calculate coverage for binning rather than the concatenated *.fastq* file. I also often increase the minimum contig length to 2500 or 3000 to help get better quality bins and speed up the binning process.

## Step 4: Refine bins

In the step, we refine the bins to a finalized MAG dataset. This involves using the *metawrap bin_refinement* tool to combine redundant bins across the 3 binning algorithms (*Maxbin2*, *Metabat2*, *CONCOCT*) and filtering out bins based on *checkM*-calculated completeness and contamination thresholds. We generally choose < 10 % contamination and > 50% completeness.

Here is the script called upon by the slurm submission script.

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

I often try multiple binning efforts for the same sample depending on what we know about the sample or how well our first binning attempts went. For example, we have tried binning using just the sample *.fastq* files and then decided to bin using multiple *.fastq* files or with different minimum contig lengths. Luckily, the *metabat bin_refinement* tool can be used to refine all the different bins for the same sample. 



```bash
cat > bin_refinement_second_round.sh 
arrayid=$(($1 + 1))

sample=$(sed -n "${arrayid}p" /PATH/TO/MY_SAMPLING_SITE/sample_array.txt)

echo $sample

metawrap bin_refinement -o ${sample}/BIN_REFINEMENT_FINAL -t 24 -m 96 -A /PATH/TO/MY_SAMPLING_SITE/${sample}/BIN_REFINEMENT/metawrap_50_10_bins -B /PATH/TO/MY_SAMPLING_SITE/${sample}/BIN_REFINEMENT_2/metawrap_50_10_bins -c 50 -x 10 

echo $sample
echo finished
```

Rename the refined bins to include the sample name instead of just bin.1.fa to SAMPLE_1_bin_1.fa. This is a bit of a painful, manual way depending on how many samples you have.

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

Once I have all of the MAGs from all of the metagenomes together in one directory, I move onto analyzing the MAGs. See the MAG_processing.md for next steps!

