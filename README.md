# Genomic surveillance of Salmonella enterica Minnesota strains from poultry products imported into South Africa
##### Authors: Vishnu Raghuram, Thendo Mafuna, Vignesh Ramnath, Hadrien Gourle, Josefin Blom, Kudakwashe Magwedere,Laura Carroll, Itumeleng Matle

SEPI - Salmonella enterica poultry imports . Data analysis code for the study 'Genomic monitoring of Non-typhoidal Salmonella enterica serotype Minnesota associated with poultry meat imports from Brazil to South Africa'

Link to raw data and intermediate files - https://doi.org/10.5281/zenodo.17224261

Link to preprint - https://www.medrxiv.org/content/10.1101/2025.04.16.25325939v2 (Currently still the version prior to revisions)

## Load libraries

```{r, echo=F}

library(data.table)
library(dplyr)
library(tidyverse)
library(ggtree)
library(scales)
library(cowplot)
library(treeio)
library(caper)
library(ggnewscale)
library(phytools)
library(Rlsd2)
library(umap)
library(ggside)
library(treeio)
library(ggsignif)
library(ggh4x)

```

### Full data table

``` r
ST548_full_table<-fread("data_tables/TableS1.tsv",header=T,sep="\t")
ST548_full_table<-as.data.frame(ST548_full_table)
rownames(ST548_full_table)<-ST548_full_table$sample_name

sepi<-ST548_full_table$sample_name[ST548_full_table$`Data source`=="POE" | ST548_full_table$`Data source`=="Supermarket"]
```

## Diversity (in core SNPs) of the 36 Brazil-to-ZA genomes

``` r
only_ZA_snps<-fread("Supplemental datasets/Dataset_S7.tsv",header=F,sep="\t",col.names=c("acc1","acc2","snps"))

only_ZA_snps_matrix <- dcast(only_ZA_snps, acc1 ~ acc2, value.var = "snps")
only_ZA_snps_matrix<-as.data.frame(only_ZA_snps_matrix)

rownames(only_ZA_snps_matrix)<-only_ZA_snps_matrix$acc1
only_ZA_snps_matrix <- only_ZA_snps_matrix[,-1]

my_color<-list(Continent= c(Africa = "#4CC3FF",Asia = "#FFC34C", Europe = "#290AD8",`North America` = "#FF2A00",`South America` = "#FFFF00"
),`Data source`=c(Supermarket = "black",POE="grey70",Public="white"))

only_ZA_heatmap<-pheatmap::pheatmap(only_ZA_snps_matrix,
                   show_rownames = T, 
                   show_colnames = T,
                   annotation_row = ST548_full_table[5],
                   annotation_names_row = F,
                   annotation_col = ST548_full_table[5],
                   annotation_names_col = F,
                   annotation_colors = my_color,
                   color = (c("#555555","#5B3794","#6E3A89","#74599C","#8085B8FF","#8DA3CAFF" ,"#A0BDD8FF" ,"#B7D5E4FF","#F1F1F1FF") ),
                   breaks = c(0,10,20,30,40,50,60,70,80),
                   legend_breaks = c(0,10,20,30,40,50,60,70,80),
                  display_numbers = T ,
                  number_format = "%d",
                  number_color = "black",
                  fontsize_number=8,
                  border_color = "white",
                  clustering_method = "average",
                  silent = T
                   )

Fig1A<-ggplotify::as.ggplot(only_ZA_heatmap)

# Remove self pairs and duplicate pairs
only_ZA_snps <- only_ZA_snps[only_ZA_snps$acc1 != only_ZA_snps$acc2,] %>%
  mutate(pair = ifelse(acc1 < acc2, paste(acc1, acc2, sep = "_"), paste(acc2, acc1, sep = "_"))) %>%
  distinct(pair, .keep_all = TRUE) %>%
  dplyr::select(-pair)

only_ZA_snps<-merge(only_ZA_snps,ST548_full_table[,c(1,5)],by.x="acc1",by.y="sample_name")
only_ZA_snps<-merge(only_ZA_snps,ST548_full_table[,c(1,5)],by.x="acc2",by.y="sample_name")
only_ZA_snps<-only_ZA_snps %>% mutate(concordance = ifelse(`Data source.x` < `Data source.y`, paste(`Data source.x`, `Data source.y`, sep = " vs "), paste(`Data source.y`, `Data source.x`, sep = " vs ")))

Fig1B<-ggplot(only_ZA_snps,aes(x=snps))+
    geom_histogram(fill="black",binwidth = 1)+
    theme_bw()+
    facet_wrap(concordance~.,nrow=3,scales = "free_y")+
    theme(axis.text = element_text(color="black",size=12))+
    theme(strip.text = element_text(color="black",size=12,face="bold"))+
    theme(axis.title = element_text(color="black",face="bold",size=14))+
    theme(panel.grid.minor = element_blank())+
    scale_x_continuous(breaks=pretty_breaks(n=12))+
    scale_y_continuous(breaks=pretty_breaks(n=10))+
    labs(x="Core SNP distance",y="Number of pairs")+scale_y_continuous(breaks = pretty_breaks(n=4))

###### hclust 
hc <- hclust(as.dist(only_ZA_snps_matrix),method = "average")
hc$height <- round(hc$height, 6) 

# clustering at 5 SNP intervals from 50 to 5 in a loop

for (i in seq(50,5,-5)){
  int_table<-as.data.frame(cutree(hc,h=i))
  int_table$acc<-rownames(int_table)
  colnames(int_table)[1]=i
  int_table<-int_table[,c(2,1)]
  assign(paste0("int_table_", i), int_table)
}

table_names <- paste0("int_table_", seq(50,5,-5))

hc_table<-lapply(mget(table_names), as.data.frame) %>% reduce(left_join,by="acc")

rm(list=table_names)
rm(int_table)

x<-as.data.frame(apply(hc_table[,-1], 2, max))
x$cluster<-row.names(x)
colnames(x)<-c("Number of clusters","SNP threshold")
x$`Number of clusters`<-as.numeric(x$`Number of clusters`)
x$`SNP threshold`<-as.numeric(x$`SNP threshold`)

Fig1C<-ggplot(x,aes(x=`Number of clusters`,y=`SNP threshold`))+
  geom_point(size=2)+
  geom_line(linetype="dashed")+
  theme_bw()+
  theme(axis.text = element_text(color="black",size=12))+
  theme(axis.title = element_text(color="black",face="bold",size=14))+
  theme(panel.grid.minor = element_blank())+
  scale_x_continuous(breaks=pretty_breaks(n=12))+
  scale_y_continuous(breaks=pretty_breaks(n=10))

Fig1<-plot_grid(Fig1A,NULL,plot_grid(Fig1B,NULL,Fig1C,nrow=3,ncol=1,rel_heights = c(1,0.02,0.5),labels = c("B","","C")),nrow=1,ncol=3,rel_widths = c(1,0.02,0.5),labels=c("A","",""))

ggsave(plot = Fig1,filename = "manuscript_figs/tiff/Fig 1.tif",device="tiff",units="in",dpi=600,width = 18,height=9,compression = "lzw")
ggsave(plot = Fig1,filename = "manuscript_figs/png/Fig 1.png",device="png",units="in",dpi=600,width = 18,height=9)
```

<p align="center">

<img src="manuscript_figs/png/Fig 1.png" width="70%" height="70%"/>

</p>

## Acquisition of publicly available data from Enterobase and NCBI

### Get NCBI assembly accessions

``` bash

##### Download NCBI genome summary (Was done on 18th July 2024 - file provided in 'data_tables' folder)
wget "https://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt"
head -2 assembly_summary_genbank.txt | tail -n +2 | sed 's/#//' > Salmonella_assembly_summary_genbank.txt
grep "Salmonella" assembly_summary_genbank.txt >> Salmonella_assembly_summary_genbank.txt
```

### Read Enterobase assembly table and filter samples

``` r
enterobase_assembly<-fread("data_tables/enterobase_all_assemblystats_20250729.txt",header=T,sep="\t",na.strings = c("","NA","ND","N/D","N/A","NaN"),select=c(1,6,8,12,13,27,29,34,35,36,37,38))

# Set country as "Unknown" for samples with NA in country column.
enterobase_assembly$Country[is.na(enterobase_assembly$Country)] <- 'Unknown'

enterobase_assembly$`Source Type`[is.na(enterobase_assembly$`Source Type`)] <- 'Unknown'

enterobase_assembly$`Continent`[is.na(enterobase_assembly$`Continent`)] <- 'Unknown'

# Filter only samples with assemblies,year of isolation, biosample/bioproject IDs and assembly stats
enterobase_assembly<-enterobase_assembly[complete.cases(enterobase_assembly)]

# Read MLST table and get only ST548
enterobase_mlst<-fread("data_tables/enterobase_all_mlst_20250729.txt",header=T,sep="\t",na.strings = c("","NA","ND","-","NaN"),select=c(1,34))
enterobase_mlst<-enterobase_mlst[enterobase_mlst$ST=="548"]
enterobase_assembly<-merge(enterobase_assembly,enterobase_mlst,by="Uberstrain")
remove(enterobase_mlst)

# Read Enterobase serotyping table and get only "Minnesota"
enterobase_serotype<-fread("data_tables/enterobase_all_serotype_20250729.txt",header=T,sep="\t",na.strings = c("","NA","ND","-","NaN"),select=c(1,33,37))
enterobase_serotype<-enterobase_serotype[complete.cases(enterobase_serotype)]
enterobase_serotype<-enterobase_serotype[enterobase_serotype$`SISTR1 Serovar`=="Minnesota" & enterobase_serotype$`SeqSero2 Serovar`=="Minnesota"]

enterobase_assembly<-merge(enterobase_assembly,enterobase_serotype,by="Uberstrain")
remove(enterobase_serotype)

ncbi_assembly<-fread("data_tables/Salmonella_assembly_summary_genbank_20250728.txt",header=T,sep="\t",select=c(1,3,12,20,26,31))

# NOTE: Not all assemblies in Enterobase were found in the NCBI summary table as of 29/07/2025.
enterobase_assembly<-merge(enterobase_assembly,ncbi_assembly,by.x="Sample ID",by.y="biosample",all.x=T)
remove(ncbi_assembly)

enterobase_assembly$`Data source` <- ifelse(is.na(enterobase_assembly$assembly_accession), "Public (SRA)", "Public (GenBank)")

# Filter out only samples with N50\>20000, \<500 contigs, species % >=85 for genomes in Genbank only,for samples without ncbi assembly keep everything
enterobase_assembly<-enterobase_assembly %>% separate_wider_delim(Species, delim = ";", names = c("Species", "Species_percent"))
enterobase_assembly$Species_percent<-gsub("%","",enterobase_assembly$Species_percent)
enterobase_assembly$Species_percent<-as.numeric(enterobase_assembly$Species_percent)

enterobase_assembly<-enterobase_assembly[(enterobase_assembly$Species=="Salmonella enterica" & enterobase_assembly$Species_percent>85 & enterobase_assembly$N50>50000 & enterobase_assembly$`Contig Number(>=200 bp)`<500 & enterobase_assembly$`Data source`=="Public (GenBank)") | (enterobase_assembly$`Data source`=="Public (SRA)"),]

write.table(enterobase_assembly,file="Supplemental datasets/Dataset_S2.tsv",col.names=T,row.names=F,sep="\t",quote=F)
```

### Download reads for samples not found in NCBI and assemble

``` bash

# Run this commad to get a list of SRA accessions for assemblies only in enterobase but not in NCBI

grep -f <(grep "Public (SRA)" enterobase_minnesota_filtered.tsv | cut -f1) data_tables/enterobase_all_assemblystats_20250729.txt | cut -f3 | cut -f1 -d";"

# Pipe the list of SRA accessions from the above command to fasterq-dump
fasterq-dump $1 -e 4 -O assemblies_not_in_ncbi/raw_reads/

# Using the same list of SRA accessions, run fastp on all the raw reads
fastp -w 4 -i raw_reads/$1\_1.fastq -I raw_reads/$1\_2.fastq -o trimmed_reads/$1\_1_trimmed.fastq  -O trimmed_reads/$1\_2_trimmed.fastq -j trimmed_reads/$1.json -h $1.html

# Assemble with shovill (skesa)
shovill --assembler skesa --R1 trimmed_reads/$1\_1_trimmed.fastq --R2 trimmed_reads/$1\_2_trimmed.fastq --outdir shovill/$1 --cpus 128

# Run QUAST using ATCC S. Minnesota genome as reference (same as in manuscript preprint)
quast -p shared quast -o ./quast -r ../GCF_000486855.2_ASM48685v2_genomic.fna -t 8 shovill/*/*.fna

# Remove sample SAMN41918439 as it has > 1000 contigs

# Mean genome length according to seqkit stats = 4769836
cat paths_to_all_ST548_assemblies_filtered.txt | xargs -I {} seqkit stats {} -T | grep -v "sum_len" | datamash mean 5
```

### Download assemblies and inspect for outliers using mashtree

``` bash

# download assemblies using ncbi datasets
cut -f16 enterobase_minnesota_filtered_onlyncbi.tsv | xargs -I {} datasets download genome accession {} --filename {}.zip --assembly-source GenBank

# mashtree
# 4753709 is the avg length of the 428 public genomes
mashtree --file-of-files paths_to_all_public_assemblies.txt --genomesize 4769836 > mashtree.dnd
```

``` r
mash_tree<-ggtree(read.tree("misc_files/FigS1_mashtree.dnd"),size=0.5)+
  geom_tippoint(size=1)+
  geom_tiplab(align=F,geom="text",size=1)+
    geom_treescale(x = 0.01,y=400,offset = 5)+
    theme_bw()+
    theme(axis.text = element_text(color="black",size=14))

# Delete sample SAMN14167312, it seems too distant to be in the same ST

ggsave(mash_tree,file="manuscript_figs/png/FigS1.png",device = "png",dpi=600,width = 8,height = 8)
```

<p align="center">

<img src="manuscript_figs/png/FigS1.png" width="70%" height="70%"/>

</p>

## Phylogeny for ST548 genomes

### Run RLSD2 on IQ-TREE

``` r
library(Rlsd2)

avgseqlen<-4769836

tree<-read.tree("misc_files/ST548_SA_R1_iqtree.treefile")
x<-as.data.frame(tree$tip.label)
colnames(x)<-"sample_name"

# using preset rate
lsd2(inputTree = tree,inputDate = "misc_files/lsd2/ST548_public_south_africa_lsd2_dates.txt",seqLen=avgseqlen,estimateRoot = "as",confidenceInterval = 1000,outFile = "misc_files/lsd2/ST548_public_south_africa_lsd2_fixed",givenRate = 1.27e-7)

# LSD2 estimated rate
lsd2(inputTree = tree,inputDate = "misc_files/lsd2/ST548_public_south_africa_lsd2_dates.txt",seqLen=avgseqlen,estimateRoot = "as",confidenceInterval = 1000,outFile = "misc_files/lsd2/ST548_public_south_africa_lsd2_estimated")
```

### Read annotation data

``` r
# Read AMRFinderPlus data
amrfindersummary<-fread("Supplemental datasets/Dataset_S3.tsv",header=T,sep="\t",na.strings = "NA")
amrfindersummary$Class[amrfindersummary$Class %in% c("LINCOSAMIDE","LINCOSAMIDE/MACROLIDE","MACROLIDE")]<-"LINCOSAMIDE/MACROLIDE"


# Read plasmidfinder data
plasmid<-fread("Supplemental datasets/Dataset_S4.tsv",sep="\t",header=T)

plasmid_annot<-plasmid[plasmid$`%IDENTITY`>80 & plasmid$`%COVERAGE`>80,c(1,6)] %>% mutate(present=1) %>% group_by(`#FILE`,GENE) %>% summarise(Freq = sum(present)) %>% pivot_wider(names_from = GENE,values_from = Freq,values_fill = 0,)

plasmid_annot<-as.data.frame(plasmid_annot)
rownames(plasmid_annot)<-plasmid_annot$`#FILE`
plasmid_annot<-plasmid_annot[,-1]

plasmid_annot<-plasmid_annot[, order(colSums(plasmid_annot), decreasing = TRUE)]

plasmid_annot<-plasmid_annot %>% mutate(across(everything(), ~ ifelse(. != "0", "PRESENT", "ABSENT")))

# Set colour palettes
cont_pal <- c(
  Africa = "#4CC3FF",
  Asia = "#FFC34C",
  Europe = "#290AD8",
  `North America` = "#FF2A00",
  `South America` = "#FFFF00",
  Unknown = "grey90"
)

antibiotic_pal <- c(
  AMINOGLYCOSIDE = "#50FF50",
  `BETA-LACTAM` = "#FF86FF",
  COLISTIN = "mistyrose",
  FOSFOMYCIN = "#005000",
  `LINCOSAMIDE/MACROLIDE` = "yellowgreen",
  PHENICOL = "#FFB6DB",
  QUINOLONE = "#00BB00",
  SULFONAMIDE = "#BB00BB",
  TETRACYCLINE = "#BBFFBB",
  TRIMETHOPRIM = "#500050",
  ABSENT = "#FFFFFF00"
)

amrfinder_presabs<- as.data.frame(amrfindersummary[ amrfindersummary$`Element subtype`=="AMR",c(1,7,12)]) %>%
  distinct() %>%
  pivot_wider(names_from = `Gene symbol`, values_from = Class, values_fill = "ABSENT")
amrfinder_presabs<-as.data.frame(amrfinder_presabs)
rownames(amrfinder_presabs)<-amrfinder_presabs$sample_name

ordered_gene_symbols <- as.data.frame(amrfindersummary[amrfindersummary$sample_name %in% ST548_full_table$sample_name & amrfindersummary$`Element subtype`=="AMR" &  amrfindersummary$`Gene symbol`!="mdsA" & amrfindersummary$`Gene symbol` !="mdsB",]) %>% 
  dplyr::select(`Gene symbol`, Class) %>% 
  distinct() %>%
  arrange(Class) %>%
  pull(`Gene symbol`)

amrfinder_presabs <- amrfinder_presabs[, ordered_gene_symbols]

ctxm8<-unique(amrfindersummary$sample_name[amrfindersummary$`Gene symbol`=="blaCTX-M-8"])
```

### LSD2-rate-estimated ML tree

``` r
tree_lsd2<-read.nexus("misc_files/lsd2/ST548_public_south_africa_lsd2_estimated.date.nexus")

high_amr_clade_node<-507
za_snp_cluster_0_node<-577

# Get samples in 'High-AMR clade'
high_amr_clade<-clade.members(high_amr_clade_node,tree_lsd2,tip.labels = T)
za_snp_cluster_0<-clade.members(za_snp_cluster_0_node,tree_lsd2,tip.labels = T)

## To get node labels
#x<-ggtree(tree_lsd2,mrsd="2025-06-01")+geom_text2(aes(label=node))
#ggsave(filename = "node_labels_lsd2_estimated.png",plot=x,device="png",width=25,height=40,units="in",dpi=600)

t<-ggtree(tree_lsd2,mrsd="2025-06-01") %<+% ST548_full_table +
    geom_tippoint(size=1,aes(color=Continent))+
    scale_color_manual(values=cont_pal)+
    new_scale_color()+
    geom_tiplab(align = T,as_ylab = F,geom = "text",size=3,aes(label=`Data source`,color=`Data source`))+
    scale_color_manual(values = c("grey50","#FFFFFF00","#FFFFFF00","black"),guide="none")+
    theme_tree2()+
    theme(panel.grid.major.x = element_line(color="grey85"))+
    theme(axis.text=element_text(size=10,color="black"))+
    theme(legend.position = "right") + 
    geom_hilight(node=high_amr_clade_node, type='roundrect',fill="#bbcdd5")+
    geom_cladelab(node=za_snp_cluster_0_node, label="ZA SNP\ncluster 0", align=F,hjust=0.5,offset=15,fontsize=3,textcolor="#007FFF",barcolor='#007FFF',barsize=1,geom="label",fill="white")+
  scale_x_continuous(breaks = c(1700,1750,1800,1850,1900,1950,2000))+
  guides(fill = guide_legend(override.aes = list(size=7)))
  
# Add dummy column with blank values for visualization purposes
filtered_amrfinder_presabs<-amrfinder_presabs[,c(1,2,6,11,12,32,34,37,42,44)]
filtered_amrfinder_presabs<-add_column(filtered_amrfinder_presabs,` `="ABSENT",.after=3)
filtered_amrfinder_presabs<-add_column(filtered_amrfinder_presabs,`  `="ABSENT",.after=6)
filtered_amrfinder_presabs<-add_column(filtered_amrfinder_presabs,`   `="ABSENT",.after=8)
filtered_amrfinder_presabs<-add_column(filtered_amrfinder_presabs,`    `="ABSENT",.after=10)
filtered_amrfinder_presabs<-add_column(filtered_amrfinder_presabs,`     `="ABSENT",.after=12)

tt<-gheatmap(t,filtered_amrfinder_presabs,offset = 80,width=1.452,colnames=T,color = "grey80",font.size =7,colnames_offset_y = 1,hjust = 0,colnames_angle = 90,colnames_position = 'top')+
  scale_fill_manual("AMR Class",values=antibiotic_pal)+
  theme(axis.text.x=element_text(size=12,angle=90,hjust=1,vjust=0.5,color="black",face="bold"))+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  theme(legend.text=element_text(size=12))+
  theme(legend.title=element_text(size=14,face="bold"))+
  geom_vline(xintercept = c(2115,2180,2201.5,2245,2266.5,2288,2309.5,2331,2352.5,2374,2395,2438),color="black")+
  new_scale_fill()

Fig2<-gheatmap(tt,plasmid_annot[c(1:5)],offset=475,width=0.75,colnames=T,color="grey80",colnames_angle = 90,colnames_offset_y = 1,colnames_position = "top",font.size = 7,hjust = 0)+
    scale_fill_manual("Plasmid",values=c("white","black"),na.value = "white",na.translate=F)+
    theme(axis.text.x=element_text(size=12,angle=90,hjust=1,vjust=0.5,color="black",face="bold"))+
    geom_vline(xintercept = c(2517,2684),color="black")+
    ggtree::vexpand(.13, 1)+
    guides(color = "none")

ggsave(plot=Fig2,filename = "manuscript_figs/png/Fig 2.png",device="png",units="in",dpi=600,width = 28,height=35)

ggsave(plot=Fig2,filename = "manuscript_figs/tiff/Fig 2.tif",device="tiff",units="in",dpi=600,width = 20,height=23,compression = "lzw")
```

<p align="center">

<img src="manuscript_figs/png/Fig 2.png" width="70%" height="70%"/>

</p>

### High AMR subclade BEAST

#### Subsampling by SNP-dist/year/location

#### Subsampling stratifying by year and location

### Supplemental phylogenies

#### Core SNP distances across 496 ST548 genomes

``` r
snps<-fread("Supplemental datasets/Dataset_S6.tsv",header=F,sep="\t",col.names=c("Ref_file","Query_file","gubbins_snps")) #snippy-gubbins alignment

snp_matrix <- dcast(snps, Ref_file ~ Query_file, value.var = "gubbins_snps")
snp_matrix<-as.data.frame(snp_matrix)

rownames(snp_matrix)<-snp_matrix$Ref_file
snp_matrix <- snp_matrix[,-1]

#annot<-merge(ST548_full_table[,c(1,6)],za_samples_metadata[,c(1,7)],by.x="sample_name",by.y="acc",all.x=T) %>% replace_na(list(`Data source` = 'Public data'))

ST548_full_table$`AMR group`[ST548_full_table$sample_name %in% high_amr_clade]<-"High-AMR group"
ST548_full_table$`AMR group`[!ST548_full_table$sample_name %in% high_amr_clade]<-"Low-AMR group"
ST548_full_table<-as.data.frame(ST548_full_table)
rownames(ST548_full_table)<-ST548_full_table$sample_name

my_color<-list(Continent= c(Africa = "#4CC3FF",Asia = "#FFC34C", Europe = "#290AD8",`North America` = "#FF2A00",`South America` = "#FFFF00",Unknown = "grey" 
),`Data source`=c(Supermarket = "black",POE="grey70",`Public (GenBank)`="white",`Public (SRA)`="white"),`AMR group`=c(`High-AMR group`="cyan",`Low-AMR group`="white"))

FigS2A<-pheatmap::pheatmap(snp_matrix,
                   show_rownames = F, 
                   show_colnames = F,
                   annotation_row = ST548_full_table[,c(4,5,11)],
                   annotation_names_row = F,
                   annotation_col = ST548_full_table[,c(4,5,11)],
                   annotation_names_col = F,
                   annotation_colors = my_color,
                   color = (c("#000000","#6B0077FF","#6F418DFF","#8085B8FF","#8DA3CAFF" ,"#A0BDD8FF" ,"#B7D5E4FF","#F1F1F1FF","#FFFFFF") ),
                   breaks = c(0,25,50,75,100,150,200,250),
                   legend_breaks = c(0,25,50,75,100,150,200,250),
                   fontsize_row = 2,
                   fontsize_col = 2,
                   clustering_method = "average",
                   annotation_legend = F,
                   silent = T
                   )

# Dummy plot just to generate a fig legend for the heatmap
FigS2Aleg<-ggplot(data=(ST548_full_table %>% mutate(`Data source` = ifelse(str_detect(`Data source`, "Public"), "Public", `Data source`))))+
  geom_point(aes(x=Continent,y=sample_name,fill=Continent),shape=22,size=5)+
  scale_fill_manual(values=c("#4CC3FF","#FFC34C","#290AD8","#FF2A00","#FFFF00","grey"))+
  ggnewscale::new_scale_fill()+
  geom_point(aes(x=`Data source`,y=sample_name,fill=`Data source`),shape=22,size=7)+
  scale_fill_manual(values=c("grey70","white","black"))+
  ggnewscale::new_scale_fill()+
  geom_point(aes(x=`AMR group`,y=sample_name,fill=`AMR group`),shape=22,size=7)+
  scale_fill_manual(values=c("cyan","white"))+
  theme_bw()+
  theme(legend.title = element_text(face="bold"))

snps<-merge(snps,ST548_full_table,by.x="Query_file",by.y="sample_name")
snps<-merge(snps,ST548_full_table,by.x="Ref_file",by.y="sample_name")
snps$concordance[snps$`AMR group.x`==snps$`AMR group.y`]<-"same AMR group"
snps$concordance[snps$`AMR group.x`!=snps$`AMR group.y`]<-"different AMR group"

FigS2B<-ggplot(snps[!(snps$`AMR group.x`=="Low-AMR group" & snps$`AMR group.y`=="Low-AMR group"),],aes(x=gubbins_snps))+
  geom_histogram(color="black",aes(fill=concordance),bins=50)+
  scale_fill_manual("",values=c("grey","salmon"),labels=c("High-AMR vs Low-AMR","High-AMR vs High-AMR"))+
  theme_bw()+
  theme(axis.text=element_text(size=12,color="black"))+
  theme(axis.title=element_text(size=14,color="black",face="bold"))+
  labs(x="Core SNP distance",y="Number of pairs")

FigS2<-plot_grid(plot_grid(ggplotify::as.ggplot(FigS2A),NULL,get_legend(FigS2Aleg),NULL,ncol=4,rel_widths = c(1,0.03,0.1,0.03),labels=c("A","")),NULL,FigS2B,ncol=1,labels = c("","","B"),rel_heights = c(1,0.1,0.5))

ggsave(plot=FigS2,filename = "manuscript_figs/tiff/FigS2.tif",device="tiff",units="in",dpi=600,width = 10,height=12,compression = "lzw")
ggsave(plot=FigS2,filename = "manuscript_figs/png/FigS2.png",device="png",units="in",dpi=600,width = 10,height=12)
```

<p align="center">

<img src="manuscript_figs/png/FigS2.png" width="70%" height="70%"/>

</p>

##### Get closest public genome to SEPI genome

``` r
find_closest_sample <- function(pairwise_df, subset_samples) {
  pairwise_df %>%
    # Filter out rows where both Ref_file and Query_file are in the subset
    dplyr::filter(!(Ref_file %in% subset_samples & Query_file %in% subset_samples)) %>%
    # Filter rows where Ref_file is in the subset
    dplyr::filter(Ref_file %in% subset_samples) %>%
    # Group by Ref_file
    dplyr::group_by(Ref_file) %>%
    # Find the row with the minimum SNPs for each Ref_file
    dplyr::slice_min(gubbins_snps, with_ties = T) %>%
    # Select relevant columns
    dplyr::select(Ref_file, Query_file, gubbins_snps) %>%
    # Rename the output columns for clarity
    dplyr::rename(Closest_public_sample = Query_file, Distance = gubbins_snps) %>%
    dplyr::ungroup()
}

za_closest_public<-find_closest_sample(snps,ST548_full_table$sample_name[ST548_full_table$`Data source`!="Public (SRA)" & ST548_full_table$`Data source`!="Public (GenBank)" ])
za_closest_public<-merge(za_closest_public,ST548_full_table[,c(-5,-7)],by.x="Closest_public_sample",by.y="sample_name")

#write.table(za_closest_public,file="data_tables/TableS2.tsv",row.names=F,col.names=T,sep="\t",quote=F)
```

##### Get closest public ZA genome to SEPI genome

``` r
find_closest_sample <- function(pairwise_df, subset_samples) {
  pairwise_df %>%
    # Filter out rows where both Ref_file and Query_file are in the subset
    dplyr::filter(!(Ref_file %in% subset_samples & Query_file %in% subset_samples)) %>%
    # Filter rows where Ref_file is in the subset
    dplyr::filter(Ref_file %in% subset_samples) %>%
    # Group by Ref_file
    dplyr::group_by(Ref_file) %>%
    # Find the row with the minimum SNPs for each Ref_file
    dplyr::slice_min(gubbins_snps, with_ties = T) %>%
    # Select relevant columns
    dplyr::select(Ref_file, Query_file, gubbins_snps) %>%
    # Rename the output columns for clarity
    dplyr::rename(Closest_ZA_public_sample = Query_file, Distance = gubbins_snps) %>%
    dplyr::ungroup()
}

snps_ZA<-snps[(snps$Country.x == "South Africa" & snps$Country.y=="South Africa") | (snps$Ref_file %in% sepi & snps$Query_file %in% sepi),]

za_closest_public<-find_closest_sample(snps_ZA,sepi)

# Get metadata from enterobase
enterobase_assembly<-fread("data_tables/enterobase_all_assemblystats_20250729.txt",header=T,sep="\t",na.strings = c("","NA","ND","N/D","N/A","NaN"),select=c(1,6,8,13,27,29,34,35,36,37,38))

za_closest_public<-merge(za_closest_public,enterobase_assembly[enterobase_assembly$`Sample ID` %in% c(unique(za_closest_public$Closest_ZA_public_sample)),c(3,4,2,6,5)],by.x="Closest_ZA_public_sample",by.y="Sample ID")

write.table(za_closest_public,file="data_tables/TableS3.tsv",row.names=F,col.names=T,sep="\t",quote=F)
```

#### Midpoint rooted tree

``` r
tree<-read.tree("misc_files/ST548_SA_R1_iqtree.treefile")
tree<-midpoint_root(tree)

t<-ggtree(tree)  %<+% ST548_full_table +
    geom_tippoint(size=2,aes(color=Continent))+
    scale_color_manual(values=cont_pal)+
    new_scale_color()+
    geom_tiplab(align = T,as_ylab = F,geom = "text",size=3,aes(label=`Data source`,color=`Data source`))+
  scale_color_manual(values = c("grey50","#FFFFFF00","#FFFFFF00","black"),guide="none")+
  theme_tree2()+
  theme(panel.grid.major.x = element_line(color="grey85"))+
  theme(axis.text=element_text(size=10,color="black"))+
  theme(legend.position = "right") + 
  geom_hilight(node=558, type='roundrect',fill="#bbcdd5")+
  geom_cladelab(node=716, label="ZA SNP\ncluster 0", align=F,hjust=-0.1,offset=0.000014,fontsize=3,textcolor="#007FFF",barcolor='#007FFF',barsize=1,geom="label",fill="white")+
  scale_x_continuous(breaks = c(0,0.5e-05,1e-05,1.5e-05,2e-05,2.5e-05,3e-05))+
  guides(fill = guide_legend(override.aes = list(size=7)))

tt<-gheatmap(t,amrfinder_presabs,offset = 0.000014,width=1,colnames=T,color = "grey80",font.size =3,colnames_offset_y = 1,hjust = 0,colnames_angle = 90,colnames_position = 'top')+
  scale_fill_manual("AMR Class",values=antibiotic_pal)+
  theme(axis.text.x=element_text(size=12,angle=90,hjust=1,vjust=0.5,color="black",face="bold"))+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  theme(legend.text=element_text(size=12))+
  theme(legend.title=element_text(size=14,face="bold"))+
  new_scale_fill()

FigS3<-gheatmap(tt,plasmid_annot,offset=0.000054,width=1,colnames=T,color="grey80",colnames_angle = 90,colnames_offset_y = 1,colnames_position = "top",font.size = 3,hjust = 0)+
  scale_fill_manual("Plasmid",values=c("white","black"),na.value = "white",na.translate=F)+
    theme(axis.text.x=element_text(size=12,angle=90,hjust=1,vjust=0.5,color="black",face="bold"))+
  ggtree::vexpand(.1, 1)+
  guides(color = "none")


ggsave(plot=FigS3,filename = "manuscript_figs/png/FigS3.png",device="png",units="in",dpi=600,width = 20,height=23)
ggsave(plot=FigS3,filename = "manuscript_figs/tiff/FigS3.tif",device="tiff",units="in",dpi=600,width = 20,height=23,compression = "lzw")
```

<p align="center">

<img src="manuscript_figs/png/FigS3.png" width="70%" height="70%"/>

</p>

#### LSD2-fixed-rate ML tree

``` r
tree_lsd2<-read.nexus("misc_files/lsd2/ST548_public_south_africa_lsd2_fixed.date.nexus")

t<-ggtree(tree_lsd2,mrsd="2025-06-01") %<+% ST548_full_table +
    geom_tippoint(size=2,aes(color=Continent))+
    scale_color_manual(values=cont_pal)+
    new_scale_color()+
    geom_tiplab(align = T,as_ylab = F,geom = "text",size=3,aes(label=`Data source`,color=`Data source`))+
    scale_color_manual(values = c("grey50","#FFFFFF00","#FFFFFF00","black"),guide="none")+
    theme_tree2()+
    theme(panel.grid.major.x = element_line(color="grey85"))+
    theme(axis.text=element_text(size=10,color="black"))+
    theme(legend.position = "right") + 
    geom_hilight(node=high_amr_clade_node, type='roundrect',fill="#bbcdd5")+
    geom_cladelab(node=za_snp_cluster_0_node, label="ZA SNP\ncluster 0", align=F,hjust=0.7,offset=0,fontsize=3,textcolor="#007FFF",barcolor='#007FFF',barsize=1,geom="label",fill="white")+
  scale_x_continuous(breaks = c(1700,1750,1800,1850,1900,1950,2000))+
  guides(fill = guide_legend(override.aes = list(size=7)))

tt<-gheatmap(t,amrfinder_presabs,offset = 80,width=0.51,colnames=T,color = "grey80",font.size =3,colnames_offset_y = 1,hjust = 0,colnames_angle = 90,colnames_position = 'top')+
  scale_fill_manual("AMR Class",values=antibiotic_pal)+
  theme(axis.text.x=element_text(size=12,angle=90,hjust=1,vjust=0.5,color="black",face="bold"))+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  theme(legend.text=element_text(size=12))+
  theme(legend.title=element_text(size=14,face="bold"))+
  new_scale_fill()

FigS4<-gheatmap(tt,plasmid_annot,offset=200,width=0.51,colnames=T,color="grey80",colnames_angle = 90,colnames_offset_y = 1,colnames_position = "top",font.size = 3,hjust = 0)+
  scale_fill_manual("Plasmid",values=c("white","black"),na.value = "white",na.translate=F)+
    theme(axis.text.x=element_text(size=12,angle=90,hjust=1,vjust=0.5,color="black",face="bold"))+
  ggtree::vexpand(.1, 1)+
  guides(color = "none")

ggsave(plot=FigS4,filename = "manuscript_figs/png/FigS4.png",device="png",units="in",dpi=600,width = 20,height=23)
ggsave(plot=FigS4,filename = "manuscript_figs/tiff/FigS4.tif",device="tiff",units="in",dpi=600,width = 20,height=23,compression = "lzw")
```

<p align="center">

<img src="manuscript_figs/png/FigS4.png" width="70%" height="70%"/>

</p>

#### FULL LSD2-rate-estimated ML tree

``` r
tree_lsd2<-read.beast("misc_files/lsd2/ST548_public_south_africa_lsd2_estimated.date.nexus")

t<-ggtree(tree_lsd2,mrsd="2025-06-01") %<+% ST548_full_table +
    geom_tippoint(size=2,aes(color=Continent))+
    scale_color_manual(values=cont_pal)+
    new_scale_color()+
    geom_tiplab(align = T,as_ylab = F,geom = "text",size=3,aes(label=`Data source`,color=`Data source`))+
    scale_color_manual(values = c("grey50","#FFFFFF00","#FFFFFF00","black"),guide="none")+
    theme_tree2()+
    theme(panel.grid.major.x = element_line(color="grey85"))+
    theme(axis.text=element_text(size=10,color="black"))+
    theme(legend.position = "right") + 
    geom_hilight(node=high_amr_clade_node, type='roundrect',fill="#bbcdd5")+
    geom_cladelab(node=za_snp_cluster_0_node, label="ZA SNP\ncluster 0", align=F,hjust=0.7,offset=0,fontsize=3,textcolor="#007FFF",barcolor='#007FFF',barsize=1,geom="label",fill="white")+
  geom_range("CI_date", color='grey', size=2, alpha=.5,linewidth=0.7)+
  scale_x_continuous(breaks = c(1700,1750,1800,1850,1900,1950,2000))+
  guides(fill = guide_legend(override.aes = list(size=7)))

tt<-gheatmap(t,amrfinder_presabs,offset = 80,width=0.51,colnames=T,color = "grey80",font.size =3,colnames_offset_y = 1,hjust = 0,colnames_angle = 90,colnames_position = 'top')+
  scale_fill_manual("AMR Class",values=antibiotic_pal)+
  theme(axis.text.x=element_text(size=12,angle=90,hjust=1,vjust=0.5,color="black",face="bold"))+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  theme(legend.text=element_text(size=12))+
  theme(legend.title=element_text(size=14,face="bold"))+
  new_scale_fill()

FigS5<-gheatmap(tt,plasmid_annot,offset=200,width=0.51,colnames=T,color="grey80",colnames_angle = 90,colnames_offset_y = 1,colnames_position = "top",font.size = 3,hjust = 0)+
  scale_fill_manual("Plasmid",values=c("white","black"),na.value = "white",na.translate=F)+
    theme(axis.text.x=element_text(size=12,angle=90,hjust=1,vjust=0.5,color="black",face="bold"))+
  ggtree::vexpand(.1, 1)+
  guides(color = "none")

ggsave(plot=FigS5,filename = "manuscript_figs/png/FigS5.png",device="png",units="in",dpi=600,width = 20,height=23)
ggsave(plot=FigS5,filename = "manuscript_figs/tiff/FigS5.tif",device="tiff",units="in",dpi=600,width = 20,height=23,compression = "lzw")
```

<p align="center">

<img src="manuscript_figs/png/FigS5.png" width="70%" height="70%"/>

</p>

#### High AMR subclade tree

``` r
# Number of AMR classes in each sample in the high-AMR clade
amrfindersummary[amrfindersummary$`Element type`=="AMR" & amrfindersummary$sample_name %in% high_amr_clade & amrfindersummary$Class!="EFFLUX",c(1,12)] %>% distinct() %>% group_by(sample_name) %>% count()
```

    # A tibble: 241 × 2
    # Groups:   sample_name [241]
       sample_name     n
       <chr>       <int>
     1 100N            5
     2 101N            4
     3 23              5
     4 256             5
     5 263             5
     6 264             4
     7 282             5
     8 290             4
     9 291             5
    10 292             4
    # ℹ 231 more rows

``` r
tree_highamr<-read.tree("misc_files/lsd2/high_amr_clade_full_iqtree.treefile")

tree_highamr<-midpoint.root(tree_highamr)

# Using midpoint-rooted tree
t<-ggtree(tree_highamr) %<+% ST548_full_table +
    geom_tippoint(size=2,aes(color=Continent))+
    scale_color_manual(values=cont_pal)+
    new_scale_color()+
    geom_tiplab(align = T,as_ylab = F,geom = "text",size=2,aes(label=`Data source`,color=`Data source`))+
   scale_color_manual(values = c("grey50","#FFFFFF00","#FFFFFF00","black"),guide="none")+
    theme_tree2()+
    theme(panel.grid.major.x = element_line(color="grey85"))+
    theme(axis.text=element_text(size=10,color="black"))+
    theme(legend.position = "right") + 
  scale_x_continuous(breaks=c(0,2.5e-6,5e-6,7.5e-6,1e-5))+
  guides(fill = guide_legend(override.aes = list(size=7)))

tt<-gheatmap(t,amrfinder_presabs,offset = 3.5e-6,width=0.51,colnames=T,color = "grey80",font.size =3,colnames_offset_y = 1,hjust = 0,colnames_angle = 90,colnames_position = 'top')+
  scale_fill_manual("AMR Class",values=antibiotic_pal)+
  theme(axis.text.x=element_text(size=12,angle=90,hjust=1,vjust=0.5,color="black",face="bold"))+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  theme(legend.text=element_text(size=12))+
  theme(legend.title=element_text(size=14,face="bold"))+
  new_scale_fill()

FigS6A<-gheatmap(tt,plasmid_annot,offset=1e-5,width=0.51,colnames=T,color="grey80",colnames_angle = 90,colnames_offset_y = 1,colnames_position = "top",font.size = 2,hjust = 0)+
  scale_fill_manual("Plasmid",values=c("white","black"),na.value = "white",na.translate=F)+
    theme(axis.text.x=element_text(size=12,angle=90,hjust=1,vjust=0.5,color="black",face="bold"))+
  ggtree::vexpand(.1, 1)+
  guides(color = "none")+
  theme(legend.position = "none")


# using estimated rate
tree_highamr<-read.tree("misc_files/lsd2/high_amr_clade_full_iqtree.treefile")
lsd2(inputTree = tree_highamr,inputDate = "misc_files/lsd2/high_amr_clade_full_dates.tsv",seqLen=4769836,estimateRoot = "as",confidenceInterval = 1000,outFile = "misc_files/lsd2/ST548_high_AMR_clade_lsd2")
```

    Reading output trees ...
    Done.

    $rate
    [1] 4.917513e-07

    $tMRCA
    [1] 2008.769

    $newickTree

    Phylogenetic tree with 241 tips and 145 internal nodes.

    Tip labels:
      SAMN20679892, 100N, 30, 313, SAMN12053334, SAMN14207196, ...
    Node labels:
      , 100, 100, 99, 98, 99, ...

    Rooted; includes branch lengths.

    $dateNexusTree
    'treedata' S4 object that stored information
    of
        'misc_files/lsd2/ST548_high_AMR_clade_lsd2.date.nexus'.

    ...@ phylo:

    Phylogenetic tree with 241 tips and 145 internal nodes.

    Tip labels:
      SAMN20679892, 100N, 30, 313, SAMN12053334, SAMN14207196, ...
    Node labels:
      , 100, 100, 99, 98, 99, ...

    Rooted; includes branch lengths.

    with the following features available:
      'CI_date', 'CI_height', 'date'.

    # The associated data tibble abstraction: 386 × 6
    # The 'node', 'label' and 'isTip' are from the phylo tree.
        node label         isTip CI_date   CI_height  date
       <int> <chr>         <lgl> <list>    <list>    <dbl>
     1     1 SAMN20679892  TRUE  <dbl [1]> <dbl [1]>  2021
     2     2 100N          TRUE  <dbl [1]> <dbl [1]>  2022
     3     3 30            TRUE  <dbl [1]> <dbl [1]>  2020
     4     4 313           TRUE  <dbl [1]> <dbl [1]>  2022
     5     5 SAMN12053334  TRUE  <dbl [1]> <dbl [1]>  2019
     6     6 SAMN14207196  TRUE  <dbl [1]> <dbl [1]>  2020
     7     7 SAMEA14356617 TRUE  <dbl [1]> <dbl [1]>  2022
     8     8 SAMEA14453033 TRUE  <dbl [1]> <dbl [1]>  2022
     9     9 SAMN09634088  TRUE  <dbl [1]> <dbl [1]>  2018
    10    10 SAMN40371311  TRUE  <dbl [1]> <dbl [1]>  2023
    # ℹ 376 more rows

    $outResultFiles
    [1] "misc_files/lsd2/ST548_high_AMR_clade_lsd2"           
    [2] "misc_files/lsd2/ST548_high_AMR_clade_lsd2.nwk"       
    [3] "misc_files/lsd2/ST548_high_AMR_clade_lsd2.date.nexus"

``` r
tree_highamr_lsd2<-read.nexus("misc_files/lsd2/ST548_high_AMR_clade_lsd2.date.nexus")

t<-ggtree(tree_highamr_lsd2,mrsd="2025-06-01") %<+% ST548_full_table +
    geom_tippoint(size=2,aes(color=Continent))+
    scale_color_manual(values=cont_pal)+
    new_scale_color()+
    geom_tiplab(align = T,as_ylab = F,geom = "text",size=3,aes(label=`Data source`,color=`Data source`))+
    scale_color_manual(values = c("grey50","#FFFFFF00","#FFFFFF00","black"),guide="none")+
    theme_tree2()+
    theme(panel.grid.major.x = element_line(color="grey85"))+
    theme(axis.text=element_text(size=10,color="black"))+
    theme(legend.position = "right") + 
  scale_x_continuous(breaks = c(2000,2005,2010,2015,2020))+
  guides(fill = guide_legend(override.aes = list(size=7)))

tt<-gheatmap(t,amrfinder_presabs,offset = 5,width=0.51,colnames=T,color = "grey80",font.size =3,colnames_offset_y = 1,hjust = 0,colnames_angle = 90,colnames_position = 'top')+
  scale_fill_manual("AMR Class",values=antibiotic_pal)+
  theme(axis.text.x=element_text(size=12,angle=90,hjust=1,vjust=0.5,color="black",face="bold"))+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  theme(legend.text=element_text(size=12))+
  theme(legend.title=element_text(size=14,face="bold"))+
  new_scale_fill()

FigS6B<-gheatmap(tt,plasmid_annot,offset=15,width=0.51,colnames=T,color="grey80",colnames_angle = 90,colnames_offset_y = 1,colnames_position = "top",font.size = 3,hjust = 0)+
  scale_fill_manual("Plasmid",values=c("white","black"),na.value = "white",na.translate=F)+
    theme(axis.text.x=element_text(size=12,angle=90,hjust=1,vjust=0.5,color="black",face="bold"))+
  ggtree::vexpand(.1, 1)+
  guides(color = "none")

FigS6Bleg<-get_legend(FigS6B)

FigS6<-plot_grid(plot_grid((FigS6A+theme(legend.position="none")),NULL,(FigS6B+theme(legend.position="none")),labels=c("A","","B"),ncol=1,rel_heights = c(1,0.1,1)),NULL,ggplotify::as.ggplot(FigS6Bleg),labels=c("","",""),ncol=3,rel_widths = c(1,0.05,0.2))

ggsave(plot=FigS6,filename = "manuscript_figs/png/FigS6.png",device="png",units="in",dpi=600,width = 20,height=23)
ggsave(plot=FigS6,filename = "manuscript_figs/tiff/FigS6.tif",device="tiff",units="in",dpi=600,width = 20,height=23,compression = "lzw")
```

<p align="center">

<img src="manuscript_figs/png/FigS6.png" width="70%" height="70%"/>

</p>

## Accessory genome based UMAP

``` r
panaroo_presabs<-fread("Supplemental datasets/Dataset_S5.tsv",header=T,sep="\t")

gene_freq<-cbind(panaroo_presabs[,1],rowSums(panaroo_presabs[,-1]))
colnames(gene_freq)<-c("Gene","Freq")
gene_freq$perc<-gene_freq$Freq*100/max(gene_freq$Freq)
int_genes<-as.data.frame(gene_freq$Gene[gene_freq$perc < 95 & gene_freq$perc > 5])

genes_list<-as.data.frame(int_genes)
colnames(genes_list)<-"Gene"

genes_list<-merge(genes_list,panaroo_presabs,by="Gene")
genes_list<-as.data.frame(t(genes_list))

#Convert gene names to header and add accession colmumn
header.true <- function(df) {
  names(df) <- as.character(unlist(df[1,]))
  df[-1,]
}
genes_list<-header.true(genes_list)
genes_list[1:ncol(genes_list)]<-sapply(genes_list[1:ncol(genes_list)],as.numeric)

genes_list_umap<-as.data.frame(genes_list)
set.seed(1000);umap_fit<-umap(genes_list_umap)
umap_df <- as.data.frame(umap_fit$layout)
colnames(umap_df)<-c("umap1","umap2")

#Saving accessions as another columm
genes_list_umap$sample<-rownames(genes_list_umap)

#Converting rownames from accessions to indexes to match with tsne table index
rownames(genes_list_umap) <- 1:nrow(genes_list_umap)

#Reeformatting umap table

umap_df$index <- 1:nrow(umap_df)
umap_df<-merge(umap_df,genes_list_umap[,c(ncol(genes_list_umap))],by.x='index',by.y='row.names')

umap_df<-umap_df[,c(4,2,3)]
colnames(umap_df)<-c("acc","umap1","umap2")
umap_df<-merge(ST548_full_table,umap_df,by.x="sample_name",by.y="acc")

Fig3<-ggplot(data=umap_df,aes(x=umap1,y=umap2,fill=Continent))+
  geom_hline(yintercept = 0,color="grey90")+
  geom_vline(xintercept = 0,color="grey90")+
  geom_point(shape=21,size=2.5,color="black")+
  scale_fill_manual("Continent",values=c(cont_pal))+
  theme_bw()+
  theme(axis.text = element_blank())+
  theme(axis.title = element_text(size=14,face="bold"))+
  theme(axis.ticks = element_blank())+
  theme(legend.position = "right")+
  theme(legend.title.position = "top",legend.title = element_text(face="bold"))+
  theme(legend.text = element_text(size=8))+
  guides(fill = guide_legend(ncol =2 ,override.aes = list(size=7)))+
  theme(legend.text = element_text(size=14))+
  xlim(-11,11)+
  ylim(-15,15)+
  guides(color = "none")+
  ggforce::geom_mark_ellipse(aes(color=`AMR group`),fill=NA)+
  scale_color_manual(values=c("cyan","black"))+
  theme(panel.grid = element_blank())+
  annotate("text", x=-6, y=-15, label= "Low-AMR group")+
  annotate("text", x=5, y=15, label= "High-AMR group")+
  geom_segment(aes(x = -6, y = -9.4, xend = -6, yend = -13.5), arrow = arrow(length = unit(0.5, "cm")))+
  geom_segment(aes(x = 5, y = 10, xend = 5, yend = 13.5), arrow = arrow(length = unit(0.5, "cm")))+
  ggnewscale::new_scale_fill()+
  geom_xsidehistogram(position="stack",bins=50,aes(fill=`AMR group`,color=`AMR group`)) +
  geom_ysidehistogram(position="stack",bins=50,aes(fill=`AMR group`,color=`AMR group`)) +
  scale_fill_manual(values=c("cyan","black"))+
  theme(ggside.panel.scale = .3)

ggsave(plot=Fig3,filename = "manuscript_figs/tiff/Fig 3.tif",device="tiff",units="in",dpi=600,width = 7,height=5,compression = "lzw")
ggsave(plot=Fig3,filename = "manuscript_figs/png/Fig 3.png",device="png",units="in",dpi=600,width = 7,height=5)
```

<p align="center">

<img src="manuscript_figs/png/Fig 3.png" width="70%" height="70%"/>

</p>

## Genomad

``` r
genomad_summary<-fread("Supplemental datasets/Dataset_S8.tsv",header=T,sep="\t")
genomad_summary<-as.data.frame(genomad_summary)

genomad_amrfinder<-merge(amrfindersummary[amrfindersummary$`Element subtype`=="AMR" & amrfindersummary$`Gene symbol`!="mdsA" & amrfindersummary$`Gene symbol` != "mdsB",c(1,3,7,12)],genomad_summary,by=c("sample_name","Contig id")) %>% pivot_longer(cols=5:7,names_to = "genomad_id",values_to = "genomad_score")

genomad_amrfinder$genomad_id<-gsub("_score","",genomad_amrfinder$genomad_id)

Fig4<-ggplot(genomad_amrfinder[genomad_amrfinder$`Gene symbol` %in% c("aadA1","aph(3'')-Ib","aph(3')-Ia","aph(6)-Id","blaCMY-2","blaCTX-M-8","qnrB19","sul2","tet(A)","tet(B)"),],aes(x=genomad_id,y=genomad_score,colour = genomad_id))+
  geom_boxplot(width=0.5,size=0.5,outlier.size = 0.5)+
  facet_wrap2(vars(Class,`Gene symbol`),nrow=1,strip = strip_nested(background_x = elem_list_rect(fill=c("#50FF50","#FF86FF","#00BB00","#BB00BB","#BBFFBB","white","white","white","white","white","white","white","white","white","white")) ))+
  theme_bw()+
  theme(axis.title.x = element_blank())+
  theme(axis.title.y = element_text(color="black",size=14,face="bold"))+
  theme(axis.text.x = element_text(color="black",size=14,face="bold",angle=45,hjust=1,vjust=1))+
  theme(axis.text.y = element_text(color="black",size=12))+
  theme(strip.text = element_text(color="black",size=14,face="bold"))+
  theme(strip.background = element_rect(fill = NA))+
  theme(legend.position = "none")+
  theme(legend.text = element_text(size=8))+
  guides(fill = guide_legend(ncol =3 ,override.aes = list(size=7)))+
  theme(legend.text = element_text(size=14))+
  scale_color_manual(values=c("black","black","black"))+
  labs(y="Genomad score")

ggsave(Fig4,filename = "manuscript_figs/tiff/Fig 4.tif",device="tiff",units="in",dpi=600,width = 18,height=5,compression = "lzw")
ggsave(Fig4,filename = "manuscript_figs/png/Fig 4.png",device="png",units="in",dpi=600,width = 18,height=5)

### Genomad test

genomad_summary$amr[genomad_summary$`Contig id` %in% amrfindersummary$`Contig id`]<-"yes"
genomad_summary$amr[!genomad_summary$`Contig id` %in% amrfindersummary$`Contig id`]<-"no"

wilcox.test(plasmid_score ~ amr, data = genomad_summary)
```


        Wilcoxon rank sum test with continuity correction

    data:  plasmid_score by amr
    W = 106626412, p-value < 2.2e-16
    alternative hypothesis: true location shift is not equal to 0

``` r
genomad_S1<-ggplot(genomad_summary,aes(sample=plasmid_score))+
  stat_qq() + 
  stat_qq_line()+
  theme_bw()+
  theme(axis.title = element_text(color="black",size=14,face="bold"))+
  theme(axis.text = element_text(color="black",size=12))+
  labs(x="Theoretical quantiles",y="Sample quantiles")

genomad_S2<-ggplot(genomad_summary,aes(x=plasmid_score))+
  geom_histogram(fill="black",bins=100)+
  theme_bw()+
  theme(axis.title = element_text(color="black",size=14,face="bold"))+
  theme(axis.text = element_text(color="black",size=12))+
  labs(x="Genomad plasmid score",y="Count")

genomad_S3<-ggplot(genomad_summary,aes(x=amr,y=plasmid_score))+
  geom_violin()+
  geom_boxplot(width=0.1,size=1,outliers = F)+
  theme_bw()+
 theme(axis.title = element_text(color="black",size=14,face="bold"))+
  theme(axis.text.y = element_text(color="black",size=12))+
  theme(axis.text.x = element_text(color="black",size=14,face="bold"))+
  scale_x_discrete(labels=c("AMR-","AMR+"))+
  labs(y="Genomad plasmid score",x="")+
  scale_y_continuous(breaks=pretty_breaks(n=10))+
  geom_signif(comparisons = list(c("no", "yes")))


FigS7<-plot_grid(plot_grid(genomad_S1,NULL,genomad_S2,nrow=3,labels = c("A","","B"),rel_heights = c(1,0.1,1)),NULL,genomad_S3,nrow=1,labels = c("","","C"),rel_widths = c(1,0.1,1))


ggsave(filename = "manuscript_figs/tiff/FigS7.tif",device="tiff",units="in",dpi=600,width = 7,height=7,compression = "lzw",plot=FigS7)
ggsave(filename = "manuscript_figs/png/FigS7.png",device="png",units="in",dpi=600,width = 7,height=7,plot=FigS7)
```

<p align="center">

<img src="manuscript_figs/png/Fig 4.png" width="70%" height="70%"/>

<img src="manuscript_figs/png/FigS7.png" width="70%" height="70%"/>

</p>

## Comparing with ATB blactxm8+ isolates

#### Script for getting metadata from biosampleID

``` bash

esearch -db biosample -query $biosample_ID | efetch -format xml | xtract -pattern BioSample -element @accession OrganismName \
-NAME "(NA)" \
-block Id -if Id@db_label -equals "Sample name" -NAME Id \
-block Ids -element "&NAME" \
-DATE "(NA)" \
-block Attribute -if Attribute@attribute_name -equals "collection_date" -DATE Attribute \
-block Attributes -element "&DATE" \
-LOC "(NA)" \
-block Attribute -if Attribute@attribute_name -equals "geo_loc_name" -LOC Attribute \
-block Attributes -element "&LOC" \
-HOST "(NA)" \
-block Attribute -if Attribute@attribute_name -equals "host" -HOST Attribute \
-block Attributes -element "&HOST" \
-SOURCE "(NA)" \
-block Attribute -if Attribute@attribute_name -equals "isolation_source" -SOURCE Attribute \
-block Attributes -element "&SOURCE" 

# 
```

#### blaCTX-M-8 positive genomes from ATB vs SEPI genomes core SNP distances

``` r
ctxm8<-amrfindersummary$sample_name[amrfindersummary$`Gene symbol`=="blaCTX-M-8"]

s<-fread("Supplemental datasets/Dataset_S17.tsv",header=F,sep="\t",col.names = c("sample1","sample2","snps"))

# remove reference 
s_matrix <- dcast(s[s$sample1!="GCA_006209225.1_PDT000316796.1" & s$sample2!="GCA_006209225.1_PDT000316796.1",], sample1 ~ sample2, value.var = "snps")
s_matrix<-as.data.frame(s_matrix)

rownames(s_matrix)<-s_matrix$sample1
s_matrix <- s_matrix[,-1]

s18<-fread("Supplemental datasets/Dataset_S18.tsv",header=T,sep="\t",na.strings = c("NA","missing","not collected", "not provided","Not available","not applicable"))

# Reformatting country column
s18$Country<-gsub(":.*","",s18$Country)

#Reformatting dates column
s18$collection_date<-gsub("-.*","",s18$collection_date)

# Writing continents column
library(countrycode)
s18$Continent<-countrycode(sourcevar = s18$Country,origin = "country.name",destination = "continent",custom_match = c(Canada = "North America" , USA = "North America"))
s18$Continent<-gsub("Americas","South America",s18$Continent)

s18$`Data source`<-"ATB"

# Remove ATB samples that were already included in the full dataset
s18<-s18[!s18$biosample %in% ctxm8,]
# Keep only Salmonella enterica Minnesota (i.e samples the SNP-dist estimation was done for)
s18<-s18[s18$biosample %in% unique(s$sample1),]

b<-rbind((s18[,c(1,7,4,5,8,9)]),(ST548_full_table[ST548_full_table$sample_name %in% ctxm8,c(1,6,2,3,4,5)]),use.names=F)
b<-as.data.frame(b)
rownames(b)<-b$biosample
b$collection_date<-as.factor(b$collection_date)

my_color<-list(Continent= c(Africa = "#4CC3FF",Asia = "#FFC34C", Europe = "#290AD8",`North America` = "#FF2A00",`South America` = "#FFFF00",Unknown = "white"), `Data source`=c(Supermarket = "black",ATB="steelblue4",`Public (SRA)` = "white"),collection_date=c(`2016`="#FEBEB1",`2017`="#EB6F5B",`2018`="#CE352E",`2019`="#B52027",`2020`="#9C0824",`2021`="#830717",`2022`="#69000C",`2023`="#2C1B07",Unknown = "white"))

FigS8<-pheatmap::pheatmap(s_matrix,
                   show_rownames = T, 
                   show_colnames = T,
                   annotation_row = b[,c(3,5,6)],
                   annotation_names_row = F,
                   annotation_col = b[,c(3,5,6)],
                   annotation_names_col = F,
                   annotation_colors = my_color,
                   color = (c("#555555","#5B3794","#6E3A89","#74599C","#8085B8FF","#8DA3CAFF" ,"#A0BDD8FF" ,"#B7D5E4FF","#F1F1F1FF","#FFFFFF") ),
                   breaks = c(0,10,25,50,75,100,150,200,250,300),
                   legend_breaks = c(0,10,25,50,75,100,150,200,250,300),
                  display_numbers = T ,
                  number_format = "%d",
                  number_color = "black",
                  fontsize_number=12,
                  border_color = "white",
                  clustering_method = "average",
                  silent=T
                   )

ggsave(filename = "manuscript_figs/tiff/FigS8.tif",plot=ggplotify::as.ggplot(FigS8),device="tiff",units="in",dpi=600,width = 16,height=12,compression = "lzw")
ggsave(filename = "manuscript_figs/png/FigS8.png",plot=ggplotify::as.ggplot(FigS8),device="png",units="in",dpi=600,width = 16,height=12)
```

<p align="center">

<img src="manuscript_figs/png/FigS8.png" width="70%" height="70%"/>

</p>

## Plasmid-mediated AMR combinations

``` r
ds9<-fread("Supplemental datasets/Dataset_S9.tsv",header=F,sep="\t",col.names=c("sample_name","Contig id","plasmid id","qlen","slen","length","qstart","qend","sstart","send","pident","evalue","bitscore"))

plasmid_lengths<-distinct(ds9[,c("plasmid id","slen")])

ds9<-ds9[ds9$length >= 0.95*ds9$qlen & qlen > 1000,]

ds9<-merge(ds9,amrfindersummary[amrfindersummary$`Element type`=="AMR" & amrfindersummary$`Gene symbol` != "mdsA" & amrfindersummary$`Gene symbol` != "mdsB"],by=c("sample_name","Contig id"))[,c(1,2,3,18,23,24)]

amr_by_plasmid <- ds9 %>%
  group_by(sample_name, `plasmid id`) %>%
  summarise(amr_genes = list(unique(`Class`)), .groups = "drop")

# get plasmids with 3 or more AMR genes
plasmids_with_mdr <- amr_by_plasmid[sapply(amr_by_plasmid$amr_genes, length) > 2,] %>% distinct()
# collapse list of amr classes into sorted string
plasmids_with_mdr$amr_genes<-lapply(plasmids_with_mdr$amr_genes, function(x) sort(trimws(x))) %>% 
    sapply(function(x) paste(x, collapse = "|"))


#plasmids_with_mdr <- merge(plasmids_with_mdr, plasmid_lengths, by="plasmid id")
plasmids_with_mdr[,c(1,3)] %>% distinct() %>% dplyr::select(2) %>% table() %>% sort(decreasing = T)
```

    amr_genes
                         AMINOGLYCOSIDE|SULFONAMIDE|TETRACYCLINE 
                                                             161 
                            BETA-LACTAM|SULFONAMIDE|TETRACYCLINE 
                                                              90 
             AMINOGLYCOSIDE|BETA-LACTAM|SULFONAMIDE|TETRACYCLINE 
                                                              71 
                         AMINOGLYCOSIDE|BETA-LACTAM|TETRACYCLINE 
                                                               8 
                          AMINOGLYCOSIDE|BETA-LACTAM|SULFONAMIDE 
                                                               7 
                AMINOGLYCOSIDE|PHENICOL|SULFONAMIDE|TRIMETHOPRIM 
                                                               7 
    AMINOGLYCOSIDE|BETA-LACTAM|PHENICOL|SULFONAMIDE|TRIMETHOPRIM 
                                                               5 
                             AMINOGLYCOSIDE|PHENICOL|SULFONAMIDE 
                                                               5 
                            AMINOGLYCOSIDE|PHENICOL|TETRACYCLINE 
                                                               5 
                               BETA-LACTAM|PHENICOL|TETRACYCLINE 
                                                               4 
                                 BETA-LACTAM|FOSFOMYCIN|PHENICOL 
                                                               3 
                               PHENICOL|SULFONAMIDE|TETRACYCLINE 
                                                               3 
                             AMINOGLYCOSIDE|BETA-LACTAM|PHENICOL 
                                                               2 
                AMINOGLYCOSIDE|BETA-LACTAM|PHENICOL|TETRACYCLINE 
                                                               2 
               AMINOGLYCOSIDE|LINCOSAMIDE/MACROLIDE|TETRACYCLINE 
                                                               2 
               AMINOGLYCOSIDE|LINCOSAMIDE/MACROLIDE|TRIMETHOPRIM 
                                                               2 
                AMINOGLYCOSIDE|PHENICOL|SULFONAMIDE|TETRACYCLINE 
                                                               2 
                  BETA-LACTAM|LINCOSAMIDE/MACROLIDE|TRIMETHOPRIM 
                                                               2 
                           SULFONAMIDE|TETRACYCLINE|TRIMETHOPRIM 
                                                               2 
                           AMINOGLYCOSIDE|BETA-LACTAM|FOSFOMYCIN 
                                                               1 
                AMINOGLYCOSIDE|BETA-LACTAM|LINCOSAMIDE/MACROLIDE 
                                                               1 
    AMINOGLYCOSIDE|BETA-LACTAM|PHENICOL|SULFONAMIDE|TETRACYCLINE 
                                                               1 
                              AMINOGLYCOSIDE|FOSFOMYCIN|PHENICOL 
                                                               1 
                 AMINOGLYCOSIDE|FOSFOMYCIN|PHENICOL|TETRACYCLINE 
                                                               1 
                          AMINOGLYCOSIDE|FOSFOMYCIN|TETRACYCLINE 
                                                               1 
                    BETA-LACTAM|FOSFOMYCIN|PHENICOL|TETRACYCLINE 
                                                               1 
                             BETA-LACTAM|FOSFOMYCIN|TETRACYCLINE 
                                                               1 
                   BETA-LACTAM|LINCOSAMIDE/MACROLIDE|SULFONAMIDE 
                                                               1 
                                BETA-LACTAM|PHENICOL|SULFONAMIDE 
                                                               1 
                                 FOSFOMYCIN|PHENICOL|SULFONAMIDE 
                                                               1 
                                FOSFOMYCIN|PHENICOL|TETRACYCLINE 
                                                               1 
                     LINCOSAMIDE/MACROLIDE|PHENICOL|TETRACYCLINE 
                                                               1 
                  LINCOSAMIDE/MACROLIDE|SULFONAMIDE|TETRACYCLINE 
                                                               1 

### blaCTX plasmid read mapping

#### Using bbmap - samtools - bedtools workflow to get read mapping information

``` bash

# Loop through list of reference fastas (blaCTX-M-8 positive/negative plasmids from PLSDB & rpoD reference sequence)

# Create a `script of scripts` - each line having the bbmap-samtools-bedtools workflow for each genome - reference sequence combination 
for refname in $(ls reference_fnas/); do
    while IFS=$'\t' read -r R1 R2; do 
    prefix=$( basename $R1 | sed 's/_1_trimmed.fastq.gz//'); 
    echo "bbmap.sh nodisk=t ref=reference_fnas/$refname.fna in=$R1 in2=$R2 out=stdout | samtools view -bS | samtools sort | bedtools genomecov -ibam stdin -bga > bedgraph_files/$prefix\_$refname.bg" ; done  < <(cat paths_to_reads.txt)
    done > script_of_scripts.sh

# Run script of scripts
bash ./script_of_scripts.sh

# Generate summary files by concatenating .bg files
cat bedgraph_files/*_rpoD.bg > blactxm8_rpoD_align_summary.txt
cat bedgraph_files/*_plasmidref.bg > blactxm8_plsdb_plasmids_align_summary.txt
```

#### Calculate coverage and number of mapped reads (Done on dardel)

``` r
library(IRanges)

# All SEPI mapped to mmseqs representative set of all blaCTX-M-8 positive plasmids in plsDB
df_plasmid<-fread("misc_files/align_blactxm8_reads_to_plasmids/blactxm8_plsdb_plasmids_align_summary.txt",sep="\t",header=F,col.names=c("sample_name","plasmid_name","start","end","no_of_reads"))

ctxm8<-amrfindersummary$sample_name[amrfindersummary$`Gene symbol`=="blaCTX-M-8"]

df_plasmid<-df_plasmid[df_plasmid$no_of_reads>0,]

# All SEPI mapped to rpoD from panaroo pan_genome_reference.fasta
df_rpoD<-fread("misc_files/align_blactxm8_reads_to_plasmids/blactxm8_rpoD_align_summary.txt",sep="\t",header=F,col.names=c("sample_name","plasmid_name","start","end","mapped"))

## Number of reads per sample mapped to rpod
rpoD_reads <- df_rpoD %>%
  group_by(sample_name) %>%
  summarise(total_rpoD_reads = sum(mapped), .groups = "drop")

# Plasmid lengths
df_lengths<-fread("misc_files/align_blactxm8_reads_to_plasmids/blactxm8_plsdb_plasmids.fasta.fai",sep="\t",header=F,select = c(1,2),col.names = c("plasmid_name","length"))

## Merge lengths with df_plasmid
df_plasmid <- merge(df_plasmid,df_lengths,by=c("plasmid_name")) %>%
  group_by(sample_name, plasmid_name, length)


compute_coverage <- function(df_sub, plasmid_length) {
  ranges <- IRanges(start = df_sub$start, end = df_sub$end)
  reduced <- reduce(ranges)
  covered_bases <- sum(width(reduced))
  return(covered_bases / plasmid_length)
}


results <- list()
for (g in group_split(df_plasmid)) {
  sample <- unique(g$sample_name)
  plasmid <- unique(g$plasmid_name)
  len <- unique(g$length)
  
  # Coverage
  cov <- compute_coverage(g, len)
  
  # Read count
  total_reads <- sum(g$no_of_reads)
  
  # rpoD normalization
  rpoD_sample_reads <- rpoD_reads %>%
    filter(sample_name == sample) %>%
    pull(total_rpoD_reads)
  
  norm_reads <- ifelse(length(rpoD_sample_reads) == 0 || rpoD_sample_reads == 0,
                       NA,
                       total_reads / rpoD_sample_reads)
  
  results[[length(results)+1]] <- data.frame(
    sample_name = sample,
    plasmid_name = plasmid,
    coverage = cov,
    total_reads = total_reads,
    total_rpoD_reads = rpoD_sample_reads,
    normalized_reads = norm_reads
  )
}

df_coverage<-bind_rows(results)

# Add plasmid ctxm8 information
p_b<-fread("misc_files/align_blactxm8_reads_to_plasmids/plasmid_blactxm8.tsv",header=F,sep="\t",col.names = c("plasmid_name","plasmid_blactxm8"))
df_coverage<-merge(df_coverage,p_b,by="plasmid_name")

# Add genomic ctxm8 information
df_coverage$blactxm8 <- ifelse(df_coverage$sample_name %in% ctxm8, "yes", "no")

# Rename column
df_coverage$`Log10 normalized no. of reads`<-log10(df_coverage$normalized_reads)

write.table(df_coverage,file="Supplemental datasets/Dataset_S14.tsv",quote=F,row.names=F,col.names=T,sep="\t")
```

#### Plots

##### SEPI mapping to blaCTX-M-8 plasmids

``` r
df_coverage<-fread("Supplemental datasets/Dataset_S14.tsv",sep="\t",header=T)

# Pangenome

plasmidpan<-fread("Supplemental datasets/Dataset_S12.tsv",header=T,sep="\t")

plasmidpan<-as.data.frame(plasmidpan)
rownames(plasmidpan)<-plasmidpan$Gene
plasmidpan<-plasmidpan[plasmidpan$`CP134392 (E. coli)`==1 | plasmidpan$`CP146619 (S. Mbandaka)`==1,]

plasmidpanheat<-t(plasmidpan[,-1])

fig5A<-ggplotify::as.ggplot(pheatmap::pheatmap(plasmidpanheat,silent = T,legend = F,fontsize = 10,cutree_rows = 3,border_color = "black",color = c("white","black")))

# Read mapping

## stats

### Test for normality
ggplot(df_coverage[df_coverage$sample_name %in% sepi],aes(sample=coverage))+stat_qq() + stat_qq_line()
```

![](SEPI_R1_files/figure-commonmark/unnamed-chunk-28-1.png)

``` r
ggplot(df_coverage[df_coverage$sample_name %in% sepi],aes(sample=`Log10 normalized no. of reads`))+ stat_qq() + stat_qq_line()
```

![](SEPI_R1_files/figure-commonmark/unnamed-chunk-28-2.png)

``` r
### Unequal variance test
var.test(df_coverage$coverage[df_coverage$sample_name %in% sepi & df_coverage$plasmid_blactxm8 =="yes" & df_coverage$blactxm8=="yes"],df_coverage$coverage[df_coverage$sample_name %in% sepi & df_coverage$plasmid_blactxm8 =="yes" & df_coverage$blactxm8=="no"])
```


        F test to compare two variances

    data:  df_coverage$coverage[df_coverage$sample_name %in% sepi & df_coverage$plasmid_blactxm8 == "yes" & df_coverage$blactxm8 == "yes"] and df_coverage$coverage[df_coverage$sample_name %in% sepi & df_coverage$plasmid_blactxm8 == "yes" & df_coverage$blactxm8 == "no"]
    F = 1.1526, num df = 279, denom df = 439, p-value = 0.1851
    alternative hypothesis: true ratio of variances is not equal to 1
    95 percent confidence interval:
     0.9342404 1.4292626
    sample estimates:
    ratio of variances 
              1.152597 

``` r
var.test(df_coverage$coverage[df_coverage$sample_name %in% sepi & df_coverage$plasmid_blactxm8 =="no" & df_coverage$blactxm8=="yes"],df_coverage$coverage[df_coverage$sample_name %in% sepi & df_coverage$plasmid_blactxm8 =="no" & df_coverage$blactxm8=="no"])
```


        F test to compare two variances

    data:  df_coverage$coverage[df_coverage$sample_name %in% sepi & df_coverage$plasmid_blactxm8 == "no" & df_coverage$blactxm8 == "yes"] and df_coverage$coverage[df_coverage$sample_name %in% sepi & df_coverage$plasmid_blactxm8 == "no" & df_coverage$blactxm8 == "no"]
    F = 0.8267, num df = 139, denom df = 219, p-value = 0.2232
    alternative hypothesis: true ratio of variances is not equal to 1
    95 percent confidence interval:
     0.6145868 1.1235228
    sample estimates:
    ratio of variances 
             0.8266989 

``` r
var.test(df_coverage$`Log10 normalized no. of reads`[df_coverage$sample_name %in% sepi & df_coverage$plasmid_blactxm8 =="yes" & df_coverage$blactxm8=="yes"],df_coverage$coverage[df_coverage$sample_name %in% sepi & df_coverage$plasmid_blactxm8 =="yes"])
```


        F test to compare two variances

    data:  df_coverage$`Log10 normalized no. of reads`[df_coverage$sample_name %in% sepi & df_coverage$plasmid_blactxm8 == "yes" & df_coverage$blactxm8 == "yes"] and df_coverage$coverage[df_coverage$sample_name %in% sepi & df_coverage$plasmid_blactxm8 == "yes"]
    F = 6.3584, num df = 279, denom df = 719, p-value < 2.2e-16
    alternative hypothesis: true ratio of variances is not equal to 1
    95 percent confidence interval:
     5.250059 7.767016
    sample estimates:
    ratio of variances 
              6.358429 

``` r
var.test(df_coverage$`Log10 normalized no. of reads`[df_coverage$sample_name %in% sepi & df_coverage$plasmid_blactxm8 =="no" & df_coverage$blactxm8=="yes"],df_coverage$coverage[df_coverage$sample_name %in% sepi & df_coverage$plasmid_blactxm8 =="no"])
```


        F test to compare two variances

    data:  df_coverage$`Log10 normalized no. of reads`[df_coverage$sample_name %in% sepi & df_coverage$plasmid_blactxm8 == "no" & df_coverage$blactxm8 == "yes"] and df_coverage$coverage[df_coverage$sample_name %in% sepi & df_coverage$plasmid_blactxm8 == "no"]
    F = 3.5861, num df = 139, denom df = 359, p-value < 2.2e-16
    alternative hypothesis: true ratio of variances is not equal to 1
    95 percent confidence interval:
     2.739846 4.775525
    sample estimates:
    ratio of variances 
               3.58614 

``` r
# Both close enough to normal dist, do parametric t test (no significant difference in variance)

readssig<-df_coverage[df_coverage$sample_name %in% sepi] %>%
  group_by(plasmid_blactxm8) %>%
  summarize(
    p_value = t.test(`Log10 normalized no. of reads` ~ blactxm8,var.equal = T)$p.value,
    blactxm8_pos_mean = mean(`Log10 normalized no. of reads`[blactxm8 == "yes"]),
    blactxm8_neg_mean = mean(`Log10 normalized no. of reads`[blactxm8 == "no"]),
    .groups = "drop"
  ) %>%
  mutate(p_adj = p.adjust(p_value, method = "bonferroni"))

  readssig$sig <- ifelse(readssig$p_adj < 0.01, "*", " ")
  readssig$metric<-"Log10 normalized no. of reads"

covsig<-df_coverage[df_coverage$sample_name %in% sepi] %>%
    group_by(plasmid_blactxm8) %>%
    summarize(
        p_value = t.test(`coverage` ~ blactxm8,var.equal = T)$p.value,
        blactxm8_pos_mean = mean(`coverage`[blactxm8 == "yes"]),
        blactxm8_neg_mean = mean(`coverage`[blactxm8 == "no"]),
        .groups = "drop"
    ) %>%
    mutate(p_adj = p.adjust(p_value, method = "bonferroni"))
  
  covsig$sig <- ifelse(covsig$p_adj < 0.01, "*", " ")
  covsig$metric<-"coverage"

siginf<-rbind(readssig,covsig)

fig5c<-df_coverage[df_coverage$sample_name %in% sepi,c(1,2,3,7,9,8)] %>%
    pivot_longer(cols = c(3,4),names_to = "metric",values_to = "value") %>%
    ggplot(aes(y=value,x=blactxm8,color=blactxm8))+
    geom_violin(width=0.5)+
    geom_boxplot(width=0.1,outliers = F)+
    scale_color_manual(values=c("black","red"))+
    facet_grid(metric~plasmid_blactxm8,scales = "free_y",labeller=as_labeller(c(no = "blaCTX-M-8 negative plasmid",yes="blaCTX-M-8 positive plasmid",coverage="Coverage",`Log10 normalized no. of reads`="Log10 normalized no. of reads")))+
    theme_bw()+
    theme(axis.title.x = element_text(color="black",size=14,face="bold"))+
    theme(axis.title.y = element_text(color="black",size=14,face="bold"))+
    theme(axis.text.x = element_text(color="black",size=12))+
    theme(axis.text.y = element_text(color="black",size=12))+
    theme(strip.text = element_text(color="black",size=10,face="bold"))+
    theme(legend.position = 'none')+
    theme(strip.background = element_rect(fill = NA))+
    theme(plot.margin = margin(t = 40, r = 20,b = 40,l = 20))+
    labs(x="blaCTX-M-8 detected in assembly")+
    geom_text(data=siginf,aes(label=sig,x=1.5,y=0),inherit.aes = FALSE,size = 12,fontface="bold")

fig5ac<-plot_grid(fig5A,fig5c,labels=c("A","C"),rel_widths = c(1,0.7))

fig5b<-ggplotify::as.ggplot(grid::rasterGrob(png::readPNG("misc_files/clinker_SEPI290_ecoli_mbandaka.png"),interpolate=T))

Fig5<-plot_grid(fig5ac,NULL,fig5b,labels=c("","","B"),rel_heights = c(1,0.1,0.4),nrow=3,ncol=1)


ggsave(plot=Fig5,filename = "manuscript_figs/tiff/Fig 5.tif",device="tiff",units="in",dpi=600,width = 20,height=20,compression = "lzw")
ggsave(plot=Fig5,filename = "manuscript_figs/png/Fig 5.png",device="png",units="in",dpi=600,width = 20,height=20)
```

##### ST548 (SRA Set + SEPI) read mapping

``` r
p_b<-fread("misc_files/align_blactxm8_reads_to_plasmids/plasmid_blactxm8.tsv",header=F,sep="\t",col.names = c("plasmid_name","plasmid_blactxm8"))

# Stats

df_coverage$genome_plasmid <- interaction(df_coverage$blactxm8, df_coverage$plasmid_blactxm8)

readssig<-df_coverage %>%
  group_by(plasmid_name) %>%
  summarize(
    p_value = t.test(`Log10 normalized no. of reads` ~ blactxm8)$p.value,
    blactxm8_pos_mean = mean(`Log10 normalized no. of reads`[blactxm8 == "yes"]),
    blactxm8_neg_mean = mean(`Log10 normalized no. of reads`[blactxm8 == "no"]),
    .groups = "drop"
  ) %>%
  mutate(p_adj = p.adjust(p_value, method = "fdr")) %>%
  merge(p_b,by="plasmid_name")

  readssig$sig <- ifelse(readssig$p_adj < 0.01, "*", " ")
  readssig$metric<-"Log10 normalized no. of reads"

covsig<-df_coverage %>%
  group_by(plasmid_name) %>%
  summarize(
    p_value = t.test(coverage ~ blactxm8)$p.value,
    blactxm8_pos_mean = mean(coverage[blactxm8 == "yes"]),
    blactxm8_neg_mean = mean(coverage[blactxm8 == "no"]),
    .groups = "drop"
  ) %>%
  mutate(p_adj = p.adjust(p_value, method = "fdr")) %>%
  merge(p_b,by="plasmid_name")
  
  covsig$sig <- ifelse(covsig$p_adj < 0.01, "*", " ")
  covsig$metric<-"coverage"

siginf<-rbind(readssig,covsig)

pla<-df_coverage[df_coverage$plasmid_blactxm8=="yes",c(1,2,3,7,9,8)] %>%
    pivot_longer(cols = c(3,4),names_to = "metric",values_to = "value") %>%
    ggplot(aes(y=value,x=blactxm8,color=blactxm8))+
    geom_violin(width=0.5)+
    #geom_point(size=1,position = position_jitterdodge(),aes(color=blactxm8),alpha=0.4)+
    scale_color_manual(values=c("black","red"))+
    facet_grid(metric~plasmid_name,scales = "free_y",labeller = label_wrap_gen(width=20))+
    theme_bw()+
    theme(axis.title.x = element_text(color="black",size=14,face="bold"))+
    theme(axis.title.y = element_text(color="black",size=14,face="bold"))+
    theme(axis.text.x = element_text(color="black",size=12))+
    theme(axis.text.y = element_text(color="black",size=12))+
    theme(strip.text = element_text(color="black",size=8,face="bold"))+
    theme(legend.position = 'none')+
    theme(strip.background = element_rect(fill = NA))+
    labs(x="blaCTX-M-8 detected in assembly",subtitle="Alignment against blaCTX-M-8+ plasmids")+
    geom_text(data=siginf[siginf$plasmid_name %in% df_coverage$plasmid_name[df_coverage$plasmid_blactxm8=="yes"] ,],aes(label=sig,x=0.5,y=-Inf),inherit.aes = FALSE,size = 8,fontface="bold")

plb<-df_coverage[df_coverage$plasmid_blactxm8=="no" ,c(1,2,3,7,9,8)] %>%
    pivot_longer(cols = c(3,4),names_to = "metric",values_to = "value") %>% 
    ggplot(aes(y=value,x=blactxm8,color=blactxm8))+
    geom_violin(width=0.5)+
    #geom_point(size=1,position = position_jitterdodge(),aes(color=blactxm8),alpha=0.4)+
    scale_color_manual(values=c("black","red"))+
    facet_grid(metric~plasmid_name,scales = "free_y",labeller = label_wrap_gen(width=20))+
    theme_bw()+
    theme(axis.title.x = element_text(color="black",size=14,face="bold"))+
    theme(axis.title.y = element_text(color="black",size=14,face="bold"))+
    theme(axis.text.x = element_text(color="black",size=12))+
    theme(axis.text.y = element_text(color="black",size=12))+
    theme(strip.text = element_text(color="black",size=10,face="bold"))+
    theme(legend.position = 'none')+
    theme(strip.background = element_rect(fill = NA))+
    labs(x="blaCTX-M-8 detected in assembly", subtitle="Alignment against blaCTX-M-8- plasmids")+
  geom_text(data=siginf[siginf$plasmid_name %in% df_coverage$plasmid_name[df_coverage$plasmid_blactxm8=="no"] ,],aes(label=sig,x=0.5,y=-Inf),inherit.aes = FALSE,size = 8,fontface="bold")

FigS9<-plot_grid(pla,NULL,plb,nrow=3,ncol=1,labels=c("A","","B"),rel_heights = c(1,0.1,1))

ggsave(plot=FigS9,filename = "manuscript_figs/tiff/FigS9.tif",device="tiff",units="in",dpi=600,width = 24,height=8,compression = "lzw")
ggsave(plot=FigS9,filename = "manuscript_figs/png/FigS9.png",device="png",units="in",dpi=600,width = 24,height=8)
```

### IS elements vs AMR associations

``` r
isescan<-fread("Supplemental datasets/Dataset_S15.tsv",header=F,sep="\t",col.names=c("sample_name","Contig id","IS_id","nIS","Genome_perc","bpIS","dnaLEN"))
isescan<-isescan[isescan$sample_name %in% ST548_full_table$sample_name,]

isescan_presabs<- as.data.frame(isescan[,c(1,3)]) %>%
  mutate(present = 1) %>%
  distinct() %>%
  pivot_wider(names_from = IS_id, values_from = present, values_fill = 0)
isescan_presabs<-as.data.frame(isescan_presabs)
rownames(isescan_presabs)<-isescan_presabs$sample_name
isescan_presabs<-isescan_presabs[,-1]



amrfinder_presabs01<- as.data.frame(amrfindersummary[ amrfindersummary$`Element subtype`=="AMR" & amrfindersummary$Class!="EFFLUX",c(1,7)]) %>%
  distinct() %>%
  mutate(present = 1) %>%
  merge(ST548_full_table[1],by="sample_name",all.y=T) %>%
  pivot_wider(names_from = `Gene symbol`, values_from = present, values_fill = 0) %>%
  dplyr::select(-c("NA"))
  
amrfinder_presabs01<-as.data.frame(amrfinder_presabs01)
rownames(amrfinder_presabs01)<-amrfinder_presabs01$sample_name
amrfinder_presabs01<-amrfinder_presabs01[,-1]


isescan_presabs <- as.data.frame(isescan_presabs)
amrfinder_presabs01 <- as.data.frame(amrfinder_presabs01)

cor_matrix <- cor(isescan_presabs, amrfinder_presabs01, method = ,use = "pairwise.complete.obs")
cor_matrix[is.na(cor_matrix)]<-0
cor_matrix <- cor_matrix[, ordered_gene_symbols]

cor_annot<-distinct(amrfindersummary[amrfindersummary$`Gene symbol` %in% ordered_gene_symbols,c(7,12)])
cor_annot<-as.data.frame(cor_annot)
rownames(cor_annot)<-cor_annot$`Gene symbol`

cor_annot_pal <- list(Class = c(
    AMINOGLYCOSIDE = "#50FF50",
    `BETA-LACTAM` = "#FF86FF",
    COLISTIN = "mistyrose",
    FOSFOMYCIN = "#005000",
    `LINCOSAMIDE/MACROLIDE` = "yellowgreen",
    PHENICOL = "#FFB6DB",
    QUINOLONE = "#00BB00",
    SULFONAMIDE = "#BB00BB",
    TETRACYCLINE = "#BBFFBB",
    TRIMETHOPRIM = "#500050",
    ABSENT = "#FFFFFF00")
)

FigS10<-ggplotify::as.ggplot(pheatmap::pheatmap(t(cor_matrix),
                   cluster_cols = F,
                   annotation_row = cor_annot[2],
                   annotation_colors = cor_annot_pal,
                   color = colorRampPalette(c("steelblue", "white","salmon"))( 12 ),
                  border_color = "black",
                  display_numbers = T,
                  number_color = "black",
                  legend_breaks = c(-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1),
                  silent=T
                   ))



ggsave(plot=FigS10,filename = "manuscript_figs/tiff/FigS10.tif",device="tiff",units="in",dpi=600,width = 14,height=10,compression = "lzw")
ggsave(plot=FigS10,filename = "manuscript_figs/png/FigS10.png",device="png",units="in",dpi=600,width = 14,height=10)
```

<p align="center">

<img src="manuscript_figs/png/FigS10.png" width="70%" height="70%"/>

</p>
