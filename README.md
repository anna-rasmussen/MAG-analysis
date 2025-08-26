# README 

This project contains descriptions and examples of code used to carry out different steps in genome-resolved metagenomics used in several publications including:


The code I use for processing the metagenomes to metagenome-assembled genomes (MAGs) is generally written as slurm submission scripts and uses a variety of tools that I assume are predominantly python-based. I generally use bash running these programs and occasionally R for some data wrangling. Once I have output files (any kind of csv or txt files with data) I do most analysis and visualizing in R. I have access to geneious for making gene/protein alignments and phylogenetic trees. There are some specialized tools for visualizing phylogenetic trees that I use including iToL and figtree.


In the vast majority of cases, the code here is taken directly from tutorials on the various bioinformatics tools and programs and relies heavily on the great documentation accompanying said tools. I highly recommend looking at the documentation for all of the following (and of course if you use these tools site them!):

metawrap
bowtie2
samtools
seqtk
GTDB
drep
kofamscan
phyloseq
vegan
anvi'o
