# Genomic analysis

### Overview

This includes examples of phylogenomic and pangenomic analysis with [anvi'o](https://anvio.org/), pulling out specific genes from genomes, and visualizing gene content in R.

## Pangenomics and Phylogenomics

Follow the [anvi'o pangenomics workflow](https://merenlab.org/2016/11/08/pangenomics-v2/) and [anvi'o phylogenomics workflow](https://merenlab.org/2017/06/07/phylogenomics/). The *anvi'o* documentation is phenomenal so spend some time exploring.

### Example: Concatenated ribosomal tree

Check out the [anvi'o phylogenomics workflow](https://merenlab.org/2017/06/07/phylogenomics/) to get the full scoop on using the functions. I like using the *anvi-get-sequences-for-hmm-hits* function to make alignments of ribosomal or housekeeping genes. This is carried out in *bash* and then I usually use [IQ-TREE](https://iqtree.github.io/) to make the phylogenies and visualize those trees in [iTOL](https://itol.embl.de/) or [FigTree](https://github.com/rambaut/figtree).

Generally, I like to use the GTDB species representatives (or some subset of them) in my trees. The first step I do is download a csv of the GTDB sp. reps from GTDB that includes the NCBI accession numbers. Then I save all those accession numbers in a text document and turn them into a list.

Using the [NCBI datasets tool](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/), download the genome assemblies:

```bash
conda create -n ncbi_datasets #First create a conda environment
conda activate ncbi_datasets #activate
conda install -c conda-forge ncbi-datasets-cli #install the datasets conda package

datasets download genome accession GCA_014190895.1,GCA_002721545.1,GCF_900110495.1,GCF_900113575.1 --include genome --filename Nitrosomonadaceae_GTDBspReps.zip
```

After doing all the proper unzipping and moving of files I then create anvi'o databases for all the genomes I want to include in the tree and then annotate the housekeeping genes and tRNAs.

```bash
conda activate anvio 
cd /PATH/TO/GENOMES/

for i in *fa
do
anvi-script-reformat-fasta $i \
-o $i.contig \
-l 500 \
--simplify-names
done

for i in *contig
do
anvi-gen-contigs-database -f $i \
-o $i.db \
-T 5 \
--project-name project-Nitrosomonadaceae
done

#annotate the housekeeping genes
for i in *.db
do
anvi-run-hmms -c $i -T 2 --also-scan-trnas
done

```

I find estimating the genome completeness is always a good way to find bugs. Have an external document with the genome name and path to the genome.

```bash
anvi-estimate-genome-completeness --external-genomes external_Nitrosomonadaceae.txt -o completion_Nitrosomonadaceae.txt
```

Sometimes genomes have non AGCT characters and you need to remake one of the contigs databases:

```bash
for i in *SAMPLE_bin_62.fa.contig #doesn't need to be a for loop but in case you want to do a few...
do
anvi-script-reformat-fasta $i --seq-type NT --overwrite
anvi-gen-contigs-database -f $i -o $i.db -T 1 --project-name project-Nitrospirales --force-overwrite        
anvi-run-hmms -c $i.db --also-scan-trnas        
done
```

Here is an example of the *ribosomal_genes_bacteria.txt* file:

```bash
Ribosom_S12_S23
Ribosomal_L1
Ribosomal_L13
Ribosomal_L14
Ribosomal_L16
Ribosomal_L17
Ribosomal_L18p
Ribosomal_L19
Ribosomal_L2
Ribosomal_L20
Ribosomal_L21p
Ribosomal_L22
Ribosomal_L23
Ribosomal_L27
Ribosomal_L27A
Ribosomal_L28
Ribosomal_L29
Ribosomal_L3
Ribosomal_L32p
Ribosomal_L35p
Ribosomal_L4
Ribosomal_L5
Ribosomal_L6
Ribosomal_L9_C
Ribosomal_S10
Ribosomal_S11
Ribosomal_S13
Ribosomal_S15
Ribosomal_S16
Ribosomal_S17
Ribosomal_S19
Ribosomal_S2
Ribosomal_S20p
Ribosomal_S3_C
Ribosomal_S6
Ribosomal_S7
Ribosomal_S8
Ribosomal_S9
```

And finally, an example of the code I use to make the concatenated ribosomal gene tree:

```bash
# TOTAL genomes: 10
# 36 genes in list
# use only genomes with at least 18/36 genes
# use only genes that occure in 9/10 genomes
# ideally I like to have at least 20 genes in the tree

anvi-get-sequences-for-hmm-hits --external-genomes external_Nitrosomonadaceae.txt \
                                -o concatenated-proteins-Nitrosomonadaceae.fa \
                                --hmm-source Bacteria_71 \
                                --gene-names ribosomal_genes_bacteria.txt \
                                --return-best-hit \
                                --get-aa-sequences \
                                --concatenate\
                                --max-num-genes-missing-from-bin 16 \
                                --min-num-bins-gene-occurs 9

#remove position if it is a gap in over half of genomes                               
trimal -in concatenated-proteins-Nitrosomonadaceae.fa \
       -out concatenated-proteins-Nitrosomonadaceae-clean.fa \
       -gt 0.50 
       
#make a quick tree to see how things look (do I need more references, are their long branches?)   
anvi-gen-phylogenomic-tree -f concatenated-proteins-Nitrosomonadaceae-clean.fa \
                           -o phylogenomic-tree-fastree-Nitrosomonadaceae.txt

```

Once I have an alignment I feel good about, I make a phylogeny using *iqtree* and the model finder program.

```bash
conda activate /home/groups/caf/iqtree-env

cat > iqtree.sh 

#!/bin/bash 
################# 
#SBATCH --job-name=iqtree
#SBATCH --nodes=1
#SBATCH -n 16
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2gb
#SBATCH --time=0-08:00:00
#SBATCH --partition=PARTITION
#SBATCH --output=sl-iqtree-%j.out
#SBATCH --mail-user=EMAIL@EMAIL.EMAIL
#SBATCH --mail-type=END,FAIL

source /PATH/TO/USER/.bashrc
conda activate /PATH/TO/CONDAENV/iqtree-env

################# 

iqtree -s concatenated-proteins-Nitrosomonadaceae-clean.fa\
       -nt 16\
       -m MFP \
       -bb 1000
```

### Example: Functional gene tree

## Pulling out genes of interest

I often want to get functional gene sequences to make alignments and phylogenies. I take several approaches depending on the data. Read all of the documentation on how to make a [pangenome](https://merenlab.org/2016/11/08/pangenomics-v2/).

#### using *anvi'o*

If all of my genomes have been annotated via anvi'o, are in a pangenome, and I TRUST the anvi'o annotations then I will sometimes use *anvi'o * to pull out my favorite gene. For example, ammonia monoxygenase subunit A (*amoA*).

The first step is to view the pangenome in web browser

```bash
anvi-display-pan -p Nitrosomonadaceae/Nitrosomonadaceae-PAN.db -g Nitrosomonadaceae-GENOMES.db

# click draw to see pangenome
# order genes under Layers > Order by: > gene_cluster_frequencies
# MUST click save state AND save bins Bins > store bin collection
# store both under default
# must store the states to get summary
# find genes Search > Search functions > ammonia monooxygenase
# jot down the gene cluser IDs
# ctrl + c in terminal to quit
```

Then I get the gene I want based on the gene cluster IDs

```bash
anvi-get-sequences-for-gene-clusters -p Nitrosomonadaceae/Nitrosomonadaceae-PAN.db \
-g Nitrosomonadaceae-GENOMES.db \
-o Nitrosomonadaceae_GC00001575_amoA.faa \
--gene-cluster-id GC_00001575 
```

Also, since you have genes annotated and a pangenome I really like the gene summaries for pangenomes:

```bash
anvi-summarize -p Nitrosomonadaceae/Nitrosomonadaceae-PAN.db \
               -g Nitrosomonadaceae-GENOMES.db \
               -C default \
               -o Nitrosomonadaceae_PAN_SUMMARY
               
```

Also the compute functional enrichement if you have genera or clusters you want to compare. You have to add the Genus data first using *anvi-import-misc-data*.

```bash
anvi-compute-functional-enrichment-in-pan -p Nitrosomonadaceae/Nitrosomonadaceae-PAN.db \
                                          -g Nitrosomonadaceae-ENOMES.db \
                                          --category Genus_enrich \
                                          --annotation-source KOfam \
                                          -o FUNC_enrichment_KOfam-Nitrosomonadaceae-Genus.txt \
                                          --functional-occurrence-table-output FUNC_OCCURENCE_KOfam-Nitrosomonadaceae-limited-Genus.txt                                         
```

#### using *blastp*

```bash
conda activate anvio

#first need amino acid sequences for genomes, I call using prodigal
cd GENOME_DIRECTORY
for i in *fa; do prodigal -i $i -o $i.genes -a $i.faa; done
mv *genes genes/
mv *faa genes/
cd genes/

#Then make blast databases
for f in *.faa; do makeblastdb -in $f -title $f -dbtype prot -out databases/$f; done

#Move copies of blast queries into the folder genes
#Blast through genes faa

for f in *.faa; do blastp -db databases/$f -query reference_amoA.faa -evalue 1e-20 -out blastp_amoA_results_$f.txt; done

#get just the subject amino acids
for f in *fa.faa.txt; do sed -n '/^>/p; /^\Sbjct/p' $f > $f.edittest.txt; done
for f in *.edittest.txt; do sed -e '/^Sbjct/s/^Sbjct.......//' -e 's/.....$//' $f > ${f%.*.*.*.*.*}.faa; done

```

#### based on KEGG annotations

I do this step in R after wrangling all of the kofam_scan, GhostKOALA, or eggnog gene annotation data.

Here is an example of how I wrangle data from GhostKOALA:

```R
#import GhostKOALA output from the downloaded "detail view" of the annotaitons

GK.annot.dat <- read.csv("data/MAG_GK_gene_annotation.txt", 
                        sep = "\t",
                        header = FALSE)
#get the user genome names from the contig names (this is why I rename them!)

GK.annot.dat.edit <- GK.annot.dat %>% #make a file that has the user genome
  mutate(query_id2 = V1, #no rownames in output from GhostKOALA
         user_genome = str_split_i(query_id2, "_c_", 1), #remove contig info to get just the MAG name
         user_genome = str_split_i(user_genome, "-k141", 1), #remove contig info to get just the MAG name
         user_genome = str_split_i(user_genome, "_k141", 1), #remove contig info to get just the MAG name
         user_genome = str_split_i(user_genome, "-c_", 1),
         contig = str_split_i(query_id2, "_length", 1),
         contig = ifelse(str_detect(contig, "_c_"), paste(contig, "_remove", sep = ""), contig),
         contig = str_split_i(contig, "_*[[:digit:]]*_remove", 1))#remove contig info to get just the MAG name

annot.dat <-  GK.annot.dat.edit %>%
  left_join(., MAG.qa.dat.all, by = "user_genome") %>% #combine gene data with taxonomy data

```

Once I have the data in R I start looking for my favorite genes.


```{r}
annot.dat %>%
  filter(user_genome %in% drep98.genomes$user_genome) %>%
  filter(str_detect(V3, "amoA")) %>% #look for gene name of amoA
  select(V1) %>% #select the gene IDs like: SAMPLE_BIN_5_k141_00083837_2
  write.csv("data/MAGs_amoA_gene_headers.txt", quote = FALSE, row.names = FALSE)
```

Then I use seqkit to grab the faa sequences I want. I generally have a concatenated file of all the MAG genes (the same one I annotate the genes of in fact!)

```bash
seqkit grep -i -f MAGs_amoA_gene_headers.txt MAGs.proteins.faa > amoA_MAGs.faa
```

Then I do my alignments using [Geneious](https://www.geneious.com/) which is very much NOT an open source software but I really like.


## Gene presence/absence plots in R

First I pick some of my favorite genes that I want to investigate.

Also recommend exploring [KEGG](https://www.kegg.jp/) to find the genes you like.


```R
annot.dat %>%
  filter(str_detect(V3, "ammonia mono|nitrite oxido|nirK|hao"))%>%
  mutate(product = V3,
         metabolism = ifelse(str_detect(product,  "ammonia mono|amoA|amoB|amoC"), "AMO",
                              ifelse(str_detect(product,  "nxr"), "NAR/NXR",
                              ifelse(str_detect(product,  "nirK"), "NirK",
                              ifelse(str_detect(product,  "|hao|hydroxylamine"), "HAO")))),
         pathway = ifelse(metabolism %in% c("AMO", "NAR/NXR", "HAO", "NirK"), "Nitrification", "Other")),
                   gene_short = str_split_i(product, ";", 1)) %>%
  ggplot(aes(x = user_genome,
             y = gene_short,
             fill = Genus )) +
  geom_point(size=3, 
             shape = 21)+
  scale_color_manual(values = col.div)+
  scale_fill_manual(values = col.div)+
  scale_x_discrete(position = "top")+
  facet_nested(cols = vars(Family, Genus),
             rows=vars(pathway, metabolism),
             scales = "free",
             space = "free",
             labeller = label_wrap_gen(width = 10),
             strip = strip_nested(size = "variable", 
                                  background_x = elem_list_rect(fill = c(col.div[1:2], rep("grey80", 15))),
                                  text_x = elem_list_text(angle = c(0,0,90,90,0,0,90,90,90,90,90,90,0,0,0)))
             )+
  theme(axis.text.x = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(angle=90),
        strip.text.y = element_text(angle=0),
        strip.placement = "outside",
        strip.background = element_rect(color = "grey85", fill = "grey90"),
        legend.position = "none")+
  labs(x = NULL,
       y = "MAG")

```

