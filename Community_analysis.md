# Microbial community analyis

This is a general overview of the types of community analyses I have done on non-redundant MAG datasets. I use R and rely heavily on the packages [phyloseq](https://joey711.github.io/phyloseq/), [vegan](https://github.com/vegandevs/vegan), [ggplot2](https://ggplot2.tidyverse.org/) and the rest of the  [tidyverse](https://www.tidyverse.org/).

### Considerations

1. The non-redundant MAG dataset only represents a subset of the microbial community (the subset that successfully binned) so I always consider how well (what % of) the metagenome reads map back to the non-redundant MAG dataset and take diversity analyses with a grain of salt. 

4. Some cutoff must be used to decide if a MAG is "present" in a sample. MAG A, representing organism A, could recruit a small number of reads from a closely related organism, organism B that did not bin, despite organism A not truly being present in a sample. I use a cutoff of >= 40% genome coverage (*covered.fraction* from *coverM*). I like to look at the distributions of the *covered.fraction* and *mean* coverage calculated by *coverM* and get a sense of what the dataset looks like. I have seen other cutoffs like 50% coverage breadth or a minimum mean coverage.

2. RPKG is inherently a relative abundance metric

3. MAGs are incomplete and represent a population of very closely related organisms so, again, grain of salt.

## Get all the data you need into R

Here are just a few examples of how I start putting together all of the MAG and metagenome metadata. I am usually working on the computing cluster for all of the data processing steps and then move the output files to my personal computer to work on.

```R
#read in gtdb and checkM output files first
MAG.QC <- left_join(gtdb, checkM, by = "user_genome") #combine the gtdb and checkM outputs based on MAG names
```

And import all of the *coverM* data from the *coverm_output.txt* files. I am often working with a lot of metagenomes so like using a for loop when possible to read in the data.

```R
sample.array <- read.csv("sample_array.txt", sep = "\t", header = FALSE)[-c(1),] #text file with metagenome names
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

Then I join all of the different datasets together

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

## Microbial diversity analysis

### Make a phyloseq object

+ [phyloseq](https://joey711.github.io/phyloseq/) is an awesome tool for microbial community analysis that I used for 16S rRNA amplicon sequence analysis. I follow the same principles when using MAG data and turn the MAG Reads Per Kilobase of genome per Gigabase of metagenome (RPKG) data into an *otu_table*, the GTDB file into a *tax_table*, and the metagenome metadata/environmental data into and *env_table*. This takes some wrangling!


```R
library(phyloseq)

MAG.abund <- coverM.dat.present

MAG.abund %>%
  select(user_genome) %>%
  unique(.)

abund.matrix <- MAG.abund %>% #select abundance dataframe
  select(user_genome, Sample, RPKG) %>% 
  mutate() %>% #select genome, sample that reads were recruited from, and the reads per kilobase per gigabase of genome
  spread(., Sample, RPKG) #make a long table that has sample as column

row.names(abund.matrix) <- abund.matrix$user_genome #make user_genome the row names

abund.matrix <- abund.matrix %>%
  select(-user_genome) #remove user_genome column as this info is no in the row names

abund.matrix[is.na(abund.matrix)] <- 0 #replace nas with 0s

gtdb <- gtdb151719.all %>% data.frame()

row.names(gtdb) <- gtdb$user_genome

gtdb <- gtdb %>% 
  data.frame(.) %>%
  select(-user_genome) %>% 
  as.matrix(.)

env <- coverM.dat.present %>%
  select(Sample, Site, Date, Depth) %>%
  group_by(Sample) %>%
  summarise(
            Site = unique(Site),
            Date = unique(Date),
            Depth_cm = unique(Depth)) %>%
  data.frame(.)

row.names(env) <- env$Sample


OTU <- otu_table(abund.matrix, taxa_are_rows = TRUE)
TAX <- tax_table(gtdb)
ENV <- sample_data(env)

MAG.phy <- phyloseq(OTU, TAX, ENV)

```

Once I have a *phyloseq* object I generally like to do some alpha diversity calculations or ordinations. 

### Alpha diversity

*Note*: I often have to make a second phyloseq object with rounded abundance counts since the richness function wants the *otu_table* to have only integers

```R
richness <- estimate_richness(MAG.phy, measures = c("Observed", "Chao1","se.chao1","Shannon", "Simpson", "InvSimpson")) %>%
  data.frame()
```  

### Beta diversity

```R

ord.pcoa <- ordinate(MAG.phy, method = "PCoA", distance = "bray")

plot_scree(ord.pcoa) #look at the scree plot for the ordination axes

#check out PCoA plot
plot_ordination(MAGs.phy, 
                ord.pcoa, 
                type = "site", 
                axes=1:2
                )+
  geom_point(aes(color = Depth, 
                 fill = Depth,
                 shape = Site
                 ),
             size = 3,
             color = "black")+
  labs(color = "Depth",  fill = "Depth")+
  scale_color_distiller(palette = "RdYlBu")+
  scale_fill_distiller(palette = "RdYlBu")+
  scale_shape_manual(values = c(21, 23, 22, 24))+
  #scale_size_manual(values=c(3,2))+
  theme(aspect.ratio = 16.2/27.6, #change this depending on the variance explained on each axis
        legend.position = "bottom")+
  guides(shape = guide_legend(nrow = 2), size = guide_legend(nrow = 2))
```
       
## Relative abundance plots

I like making RPKG-based abundance plots to get a sense of spatial patterns of microbes across different sampling time points or depths, sites, stations, and/or dates. Depending on the data, I agglomerate the MAGs at a certain taxonomic level to make prettier plots.

```{r}
#define colors to use
col <- c(rev(brewer.pal(6, "RdYlBu")),
             rev(brewer.pal(4, "PiYG")), 
             rev(brewer.pal(4, "PuOr"))[1:2],
             rev(brewer.pal(4, "BrBG")),
             "black",
             "grey70")

#make an aggolmerated abundance table
MAG.abund.dat.Phylum.glom <- coverM.dat.present %>% #agglomerate all MAGs at the phylum level
  select(Site, 
         Depth, 
         Year,
         Date,
         Domain, 
         Phylum,
         RPKG) %>%
  group_by(Site, Depth, Year, Date, Domain, Phylum) %>%
  summarize(abund = sum(RPKG))

#Make a dataframe for determing the abundance of different phyla
Phylum.abund.order <- coverM.dat.present %>% #agglomerate all MAGs at the phylum level
  select(Site, 
         Depth, 
         Year,
         Date,
         Domain, 
         Phylum,
         RPKG) %>%
  group_by(Phylum) %>%
  summarize(abund = sum(RPKG)) %>% #get the summed abundance of each Phylum for ordering plot
  as.data.frame() 

#Make a list of the phylum order
list.order <- Phylum.abund.order[order(Phylum.abund.order[,2], decreasing = TRUE),]

top <- c(list.order$Phylum[1:16], "Other" ) #select top 16 phyla

MAG.abund.dat.Phylum.glom <- MAG.abund.dat.Phylum.glom %>%
  mutate(Phylum2 = ifelse(Phylum %in% top, paste(Phylum), "Other")) #label all the less abundant phyla as "Other""

MAG.abund.dat.Phylum.glom$Phylum2 <- factor(MAG.abund.dat.Phylum.glom$Phylum2, levels = c(list.order$Phylum, "Other")) #order phyla by abundance

MAG.abund.dat.Phylum.glom %>%
  ggplot(aes(x=Depth,
             y= abund,
             fill= Phylum2,
             color= Phylum2)) +
  geom_col(width = 5)+
  scale_color_manual(values = col.div)+
  scale_fill_manual(values = col.div)+
  scale_x_reverse()+
  facet_nested(cols = vars(Site, Year, Date),
             scales = "free", #make axis scales in each facet variable
             space = "free", #make the width or height of each facet variable
             strip = strip_nested(size = "variable", 
                                  background_x = elem_list_rect(fill = c("olivedrab", 
                                  "skyblue3", "red2", rep("grey80", 8))))
             )+
  theme(axis.text.x = element_text(angle=60, 
                                   hjust=1, 
                                   vjust=1),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="white", color = "grey85"),
        panel.grid.major = element_line(color = "grey95"),
        legend.position = "right",
        legend.direction = "vertical")+
  labs(y = "Abundance (RPKG)",
       x = "Depth (cm)",
       title = "Top 16 phyla"
       color = "Phylum",
       fill = "Phylum"
       )+
  guides(fill = guide_legend(ncol= 1), color = guide_legend(ncol= 1))+
  coord_flip()

ggsave("figures/Fig_abundance_phylum_top16.jpg", units="in", width = 8, height=4.5)
```



