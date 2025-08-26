# MAG-analysis
Genome-resolved metagenomic analysis pipeline 

This repository contains descriptions and examples of code used to carry out different steps in genome-resolved metagenomics, including markdown files for: 
+ metagenome assembly and binning
+ metagenome-assembled genome (MAG) processing
+ microbial community analysis
+ genomic analysis

The code I use for processing the metagenomes to MAGs is generally written as slurm submission scripts and uses a variety of tools that are predominantly python-based. I generally use bash and R for some data wrangling. Once I have output files (any kind of csv or txt files with data), I do most analysis and visualizing in R. I use some specialized software for gene alignments and phylogenies.

In the vast majority of cases, the code here is taken directly from tutorials on the various bioinformatics tools and programs and relies heavily on the great documentation accompanying said tools, including :

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
