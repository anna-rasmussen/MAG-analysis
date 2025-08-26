# Genome-resolved metagenomic analysis pipeline

### Overview

This repository contains descriptions and examples of code used to carry out different steps in genome-resolved metagenomics, including: 
+ metagenome assembly and binning
+ metagenome-assembled genome (MAG) processing
+ microbial community analysis
+ genomic analysis.

The code I use for processing the metagenomes to MAGs is generally written as slurm submission scripts and uses a variety of tools that are predominantly python-based. I generally use bash and R for some data wrangling. Once I have output files (any kind of csv or txt files with data), I do most analysis and visualizing in R. I use some specialized software for gene alignments and phylogenies.

### Tools
In the vast majority of cases, the code here is taken directly from tutorials on the various bioinformatics tools and programs and relies heavily on the great documentation accompanying said tools. I highly recommend looking at the tutorials and documentation for all of the following (and of course if you use these tools cite them!):

+ [MetaWRAP](https://github.com/bxlab/metaWRAP/tree/master)
+ [bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
+ [Samtools](https://www.htslib.org/)
+ [Seqtk](https://github.com/lh3/seqtk)
+ [GTDB-Tk](https://ecogenomics.github.io/GTDBTk/)
+ [dRep](https://github.com/MrOlm/drep)
+ [kofam_scan](https://www.genome.jp/ftp/tools/kofam_scan/)
+ [CoverM](https://github.com/wwood/CoverM)
+ [InStrain](https://instrain.readthedocs.io/en/latest/index.html)
+ [phyloseq](https://joey711.github.io/phyloseq/)
+ [vegan](https://github.com/vegandevs/vegan)
+ [anvi'o](https://anvio.org/)
+ [IQ-TREE](https://iqtree.github.io/)
+ [FigTree](https://github.com/rambaut/figtree)
+ [iTOL](https://itol.embl.de/)
+ [slurm](https://slurm.schedmd.com/documentation.html)

### Published examples
Examples of publications using this (or some iteration of this) pipeline include:

+ [Rasmussen, Anna N., and Christopher A. Francis. "Genome-resolved metagenomic insights into massive seasonal ammonia-oxidizing archaea blooms in San Francisco Bay." Msystems 7.1 (2022): e01270-21.](https://doi-org.stanford.idm.oclc.org/10.1128/msystems.01270-21)
+ [Reji, Linta, et al. "Diverse ecophysiological adaptations of subsurface Thaumarchaeota in floodplain sediments revealed through genome-resolved metagenomics." The ISME Journal 16.4 (2022): 1140-1152.](https://doi-org.stanford.idm.oclc.org/10.1038/s41396-021-01167-7)
+ [Rasmussen, Anna N., and Christopher A. Francis. "Dynamics and activity of an ammonia-oxidizing archaea bloom in South San Francisco Bay." The ISME Journal 18.1 (2024): wrae148.]( https://doi-org.stanford.idm.oclc.org/10.1093/ismejo/wrae148)
+ [Rasmussen, Anna N., et al. "Diverse and unconventional methanogens, methanotrophs, and methylotrophs in metagenome-assembled genomes from subsurface sediments of the Slate River floodplain, Crested Butte, CO, USA." Msystems 9.7 (2024): e00314-24.](https://doi-org.stanford.idm.oclc.org/10.1128/msystems.00314-24)
+ [Rasmussen, Anna N., et al. "Metagenome‐Assembled Genomes for Oligotrophic Nitrifiers From a Mountainous Gravelbed Floodplain." Environmental Microbiology 27.3 (2025): e70060.](https://doi-org.stanford.idm.oclc.org/10.1111/1462-2920.70060)