#### Libraries  ---------------------------------

rm(list = ls())

set.seed(2712)


library("caper")
library("patchwork")
library("ape")
library(dplyr)
library(ggplot2)
library(tidyverse)
library("lattice")
library(reshape2)
library(adephylo)
library(phylobase)
library(data.table)
library(phytools)
library("corrplot")
library("geiger")
library("ggpubr")
library(phylolm)
library(ggtree)
library("castor")
library(plotly)
library(ggfortify)
library("erer")
library(ggtreeExtra)
library("MoreTreeTools")
library("ggstar")
library(MASS)
library(PCAtest)
library(treemapify)

#### Import genomic and ecological data  ---------------------------------


Cichlids_chemoreception_df <- 
  read.table("Cichlids_chemoreception_df.csv",
           sep=",",
           header=TRUE)

Tribe_Sp_df <- 
  Cichlids_chemoreception_df %>%
  dplyr::select(Sp, Tribe) %>%
  distinct()

#### Colors palettes  ---------------------------------


#Define gene family colors 

chemoreceptor_family_colors <- 
  c(Total_OR="lightcoral",
    Total_TAAR = "#377EB8", 
    Total_V1R = "#FFFF33",
    Total_V2R = "#FF7F00",
    Total_T1R = "#4DAF4A", 
    Total_T2R = "#984EA3")





#Define the tribe colors (same as Ronco et al. 2021)
tribes.colors <- c(Bathybatini = "#242626" , Benthochromini = "#AE262A" ,
                   Boulengerochromini = "#59595C" , Cyphotilapiini = "#FDDF13" ,
                   Cyprichromini = "#F04D29" , Ectodini = "#9AB9D9" , 
                   Eretmodini = "#682E7A" , Lamprologini = "#C588BB" ,
                   Limnochromini = "#535CA9" , Perissodini = "#FBAC43" , 
                   Trematocarini = "#959170" , Tropheini = "#86C773" , 
                   Haplochromini = "#274e13" , outgroup = "gray", Other="gray")


tribes.colors_sec <- c(Bathybatini = "#242626" , 
                       Benthochromini = "#AE262A" , 
                       Boulengerochromini = "#59595C" , 
                       Cyphotilapiini = "#FDDF13" , 
                       Cyprichromini = "#F04D29" , 
                       Ectodini = "#9AB9D9" , 
                       Eretmodini = "#682E7A" , 
                       Lamprologini = "#C588BB" , 
                       Limnochromini = "#535CA9" , 
                       Perissodini = "#FBAC43" , 
                       Trematocarini = "#959170" , 
                       Tropheini = "#86C773" , 
                       Haplochromini = "#274e13" , 
                       outgroup = "gray", 
                       Other="gray", 
                       Oreochromini="gray")


pacbio.sp.colors <- c(Batmin = "#242626" , Cypfro = "#FDDF13" , 
               Cyplep = "#F04D29" , Cunlon = "#9AB9D9" , 
               Neomul = "#C588BB" , Simdia = "#86C773")


#Initale subfam colors

olfactory_family_colors <- 
  c(OR="lightcoral",
    TAAR = "#377EB8", 
    V1R = "#FFFF33",
    V2R = "#FF7F00")

  
  
T1R_subfam_colors <- 
  c(
    'T1R3'="#D81B60",
    "T1R1A"="#1E88E5",
    "T1R1B"="#FFC107"
  )



V1R_subfam_colors <- 
  c(
    'ORA1'="black",
    "ORA2"="#E31A1C",
    "ORA3"="green4",
    "ORA4"="darkblue",
    "ORA5"="darkorange4",
    "ORA6"="orchid1")



V2R_subfam_colors <- 
  c(
    'V2RD2-1'="black",
    "V2RD3-1"="#E31A1C",
    "V2RD10-1"="maroon",
    "V2RD11-1"="green1",
    "V2RD11-2"="palegreen2",
    "V2RD11-3"="green4",
    "V2RD11-4"="darkgreen",
    "V2RD11-5"="darkblue",
    "V2RD11-6"="dodgerblue2",
    "V2RD11-7"="darkturquoise",
    "V2RD4-1"="darkorange4",
    "V2RD5-1"="gray70",
    "V2RD7-1"="orchid1",
    "V2RD8-1"="deeppink1",
    "V2RD9-1"="yellow3"
  )




TAAR_subfam_colors <- 
  c(
    'TAARL-1'="black",
    "TAARA2-1"="#E31A1C",
    "TAARA2-2"="orange2",
    "TAARB4-1"="green4",
    "TAARB4-2"="darkblue",
    "TAARB4-3"="darkorange4",
    "TAARB4-4"="orchid1",
    "TAARB4-5"="yellow3"
  )


#Define OR colors

OR_DoC_df_info <- read.table("OR_DoC_df_info.csv", header=TRUE,  sep=",")
OR_DoC_df_info <- OR_DoC_df_info %>% filter(! ID %in% c("LBE2", "LBE9"))
OR_DoC_df_info_mean <- read.table("OR_DoC_df_info_mean.csv", header=TRUE, sep=",")
OR_DoC_df_info_mean_total <- OR_DoC_df_info_mean %>% filter(Subfamily == "Total_OR") 
OR_DoC_df_info_mean_wo_total <- 
  OR_DoC_df_info_mean %>% filter(Subfamily != "Total_OR") 



OR_DoC_df_info_mean_wototal <- OR_DoC_df_info_mean %>% filter(Subfamily != "Total_OR")

library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector <- col_vector[1:47]
vector_subfam <- as.vector(OR_DoC_df_info_mean_wototal %>% pull(Subfamily) %>% unique())
names(col_vector) <- vector_subfam
OR_subfam_colors <- col_vector


#### Import species tree  ---------------------------------

#load phylo tree
b1_tree <- read.tree("b1_tree_filt.nwk")
b1_tree_wo_Neospl <- 
  drop.tip(b1_tree, "Neospl")
b1_list_sp <- b1_tree_wo_Neospl$tip.label

b1_tree_wo_Neospl_NodeLabel <- makeNodeLabel(b1_tree_wo_Neospl, method="number", prefix="Node")

#Radiation species

radiation_tree_wo_Neospl <- drop.tip(b1_tree_wo_Neospl_NodeLabel, c("Gobeth", "Hetbut", "Tilbre", "Steult", "Tilspa"))



#### Import RNA-seq data  ---------------------------------


#Import read count on olfactory receptors
read_count <- 
  read.table("Read_count_all_samples.olfactory_receptors.OGG.csv", 
             sep=",")

colnames(read_count) <- c("Sample", "gene_family", "gene_name", "read_count", "gene_type", "clade", "subfamily", "OGG")

read_count <-
  read_count %>%
  mutate(Sp_misc = str_sub(Sample, 5, -1))

read_count <-
  read_count %>%
  mutate(TissueTubeID = str_sub(Sample, 1, 4))



#rename one bad species name
read_count <-
  read_count %>%
  mutate(Sp = str_replace(Sp_misc, "Pcypre", "Pcybri"))

#Sample size of the RNA-seq
Sample_size <- 
  as.data.frame(
    read_count %>% 
      dplyr::select(Sp, Sample) %>% 
      distinct() %>% 
      group_by(Sp) %>% 
      summarise(n())
  )



read_count %>%
  dplyr::select(Sp, Sample) %>%
  distinct()


#### Import genes lengths  ---------------------------------

#Import gene lengths

gene_length_bp <- 
  read.table("gene_length.tsv", sep="\t")
colnames(gene_length_bp) <- c("gene_name", "bp")
gene_length_kbp <- gene_length_bp %>% mutate(kbp = bp/1000)



#### Import specimen informations  ---------------------------------

#Import specimen and tissues ID
Specimen_table <- 
  read.table("Specimen_RNA_table.tsv", 
             sep="\t",
             header=TRUE)
Tissue_table <- 
  read.table("Tissue_RNA_table.tsv", 
             sep="\t",
             header=TRUE)
Specimen_info_table <- left_join(Tissue_table, Specimen_table, by="SpecimenID")



#### Combine tables  ---------------------------------

read_count_specimen <- 
  left_join(read_count, Specimen_info_table, by='TissueTubeID')

Summary_RNAseq_table <- 
  left_join(read_count_specimen, gene_length_kbp, by='gene_name')


#### Compute RPK values  ---------------------------------

#Sum the read count for each exon per gene
Summary_RNAseq_table <- 
  as.data.frame(
    Summary_RNAseq_table %>%
      group_by(Sample, gene_family, gene_name, gene_type, clade, subfamily, OGG,
               Sp_misc, TissueTubeID, Sp, SpecimenID, SL, TL,
               Weight, Sex, bp, kbp) %>%
      mutate(read_count_gene = sum(read_count)) %>%
      dplyr::select(-read_count) %>%
      distinct()
  )


#Now compute the RPK values, and remove non olfactory receptors subfamilies
RPK_values <- 
  as.data.frame(
    Summary_RNAseq_table %>%
      filter(! subfamily %in% c("Kappa-1", "Lambda-1", "Theta-1", "V2RD2-1", "TAARL-1")) %>%
      group_by(Sample, gene_name) %>%
      summarise(RPK = sum(read_count_gene)/kbp))

Summary_RNAseq_table <- 
  Summary_RNAseq_table %>%
  filter(! subfamily %in% c("Kappa-1", "Lambda-1", "Theta-1", "V2RD2-1", "TAARL-1")) 


Summary_RNAseq_table <- 
  left_join(
    Summary_RNAseq_table,
    RPK_values,
    by=c("Sample", "gene_name"))




#### Compute TPM values  ---------------------------------

#If RPK is equal to NA, then it should be 0
Summary_RNAseq_table <- Summary_RNAseq_table %>% 
  mutate_at(c('RPK'), ~replace_na(.,0))

#Compute TPM value per gene
Summary_RNAseq_table <- 
  as.data.frame(
    Summary_RNAseq_table %>%
      group_by(Sample) %>%
      mutate(sum_RPK = sum(RPK)) %>%
      ungroup() %>%
      mutate(scaling_factor = sum_RPK/1000000) %>%
      mutate(TPM = RPK/scaling_factor)
  )


#Verify that the sum of TPM per sample is equal to 1000000
Summary_RNAseq_table %>%
  group_by(Sample) %>%
  summarise(sum_TPM = sum(TPM))


#### Statistics - Number of species/ind/tribes  ---------------------------------

length(Summary_RNAseq_table %>%
         pull(Sample) %>%
         unique())

length(Summary_RNAseq_table %>%
  pull(Sp) %>%
  unique())

Summary_RNAseq_table %>%
         pull(Sp) %>%
         unique()
#### Statistics - Number of reads mapped  ---------------------------------

#Import read count table on the whole genome
read_count_GW <- 
  read.table("Read_count_all_samples.all_genes.csv", sep=",")
colnames(read_count_GW) <- c("Sample", "gene", "read_count")
read_count_GW <-
  read_count_GW %>%
  mutate(Sp_misc = str_sub(Sample, 5, -1))
read_count_GW <-
  read_count_GW %>%
  mutate(TissueTubeID = str_sub(Sample, 1, 4))
read_count_GW <-
  read_count_GW %>%
  mutate(Sp = str_replace(Sp_misc, "Pcypre", "Pcybri"))


#Count the total number of reads mapped to olfactory receptors
read_count_summary_OLR <-
  as.data.frame(
    read_count %>%
      group_by(Sample) %>% 
      summarise(total_OLR_reads = sum(read_count))
  )


#Count the total number of reads mapped to non-OLR genes
read_count_summary_GW <-
  as.data.frame(
    read_count_GW %>%
      group_by(Sample) %>% 
      summarise(total_GW_reads = sum(read_count))
  )


#Combine the two tables
read_count_summary_combine <-
  left_join(read_count_summary_OLR,
            read_count_summary_GW,
            by="Sample")


#Compute the fraction of reads mapped to olfactory receptors
read_count_summary_combine <- 
  read_count_summary_combine %>%
  mutate(OLR_reads_per_million = total_OLR_reads * 1000000 / total_GW_reads)

#The fraction of reads mapping to OLR is similar to previous studies:
#https://bmcbiol.biomedcentral.com/articles/10.1186/s12915-023-01661-8
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7349279/



#### Compute RE of olfactory receptors (OLR)  ---------------------------------

#For each individual, compute the relative expression of each gene
Summary_RNAseq_table_OLR_relative_individual <- 
  as.data.frame(
    Summary_RNAseq_table %>%
      group_by(Sample) %>%
      mutate(total_reads = sum(read_count_gene)) %>%
      mutate(All_TPM = sum(TPM)) %>%
      mutate(relative_expression = TPM / All_TPM))


#For each species, compute the weighted mean of relative expressions
list_sp_RNA <- Summary_RNAseq_table %>% pull(Sp) %>% unique()
Summary_RNAseq_table_OLR_relative_species <- as.data.frame(NULL)
for(curr_sp in list_sp_RNA){
  
  curr_Summary_RNAseq_table_OLR_relative_species <- 
    as.data.frame(
      Summary_RNAseq_table %>%
        filter(Sp == curr_sp) %>%
        group_by(Sample) %>%
        mutate(total_reads = sum(read_count_gene)) %>%
        ungroup() %>%
        group_by(Sample) %>%
        mutate(All_TPM = sum(TPM)) %>%
        ungroup() %>%
        group_by(gene_family, gene_type, clade, subfamily, gene_name, OGG, kbp) %>%
        summarise(weighted_TPM = weighted.mean(TPM,total_reads))
    )
  
  
  curr_sum_weighted_TPM <- 
    curr_Summary_RNAseq_table_OLR_relative_species %>% pull(weighted_TPM) %>% sum()
  
  curr_Summary_RNAseq_table_OLR_relative_species <- 
    curr_Summary_RNAseq_table_OLR_relative_species %>% 
    mutate(relative_expression = weighted_TPM/curr_sum_weighted_TPM) %>%
    mutate(Sp = curr_sp)
  
  
  Summary_RNAseq_table_OLR_relative_species <- 
    rbind(Summary_RNAseq_table_OLR_relative_species,
          curr_Summary_RNAseq_table_OLR_relative_species)
  
}


#### Compute RE of OR genes  ---------------------------------

#At the individual level
Summary_RNAseq_table_OR_relative_individual <- 
  as.data.frame(
    Summary_RNAseq_table %>%
      filter(gene_family == "OR") %>%
      group_by(Sample) %>%
      mutate(total_reads = sum(read_count_gene)) %>%
      mutate(All_TPM = sum(TPM)) %>%
      mutate(relative_expression = TPM / All_TPM))


#At the species level
list_sp_RNA <- Summary_RNAseq_table %>% pull(Sp) %>% unique()
Summary_RNAseq_table_OR_relative_species <- as.data.frame(NULL)
for(curr_sp in list_sp_RNA){
  
  curr_Summary_RNAseq_table_OR_relative_species <- 
    as.data.frame(
      Summary_RNAseq_table %>%
        filter(gene_family == "OR") %>%
        filter(Sp == curr_sp) %>%
        group_by(Sample) %>%
        mutate(total_reads = sum(read_count_gene)) %>%
        ungroup() %>%
        group_by(Sample) %>%
        mutate(All_TPM = sum(TPM)) %>%
        ungroup() %>%
        group_by(gene_family, gene_type, clade, subfamily, gene_name, OGG, kbp) %>%
        summarise(weighted_TPM = weighted.mean(TPM,total_reads))
    )
  
  
  curr_sum_weighted_TPM <- 
    curr_Summary_RNAseq_table_OR_relative_species %>% pull(weighted_TPM) %>% sum()
  
  curr_Summary_RNAseq_table_OR_relative_species <- 
    curr_Summary_RNAseq_table_OR_relative_species %>% 
    mutate(relative_expression = weighted_TPM/curr_sum_weighted_TPM) %>%
    mutate(Sp = curr_sp)
  
  
  Summary_RNAseq_table_OR_relative_species <- 
    rbind(Summary_RNAseq_table_OR_relative_species,
          curr_Summary_RNAseq_table_OR_relative_species)
  
}


#### Compute RE of V1R genes  ---------------------------------

#At the individual level
Summary_RNAseq_table_V1R_relative_individual <- 
  as.data.frame(
    Summary_RNAseq_table %>%
      filter(gene_family == "V1R") %>%
      group_by(Sample) %>%
      mutate(total_reads = sum(read_count_gene)) %>%
      mutate(All_TPM = sum(TPM)) %>%
      mutate(relative_expression = TPM / All_TPM))



#At the species level
list_sp_RNA <- Summary_RNAseq_table %>% pull(Sp) %>% unique()
Summary_RNAseq_table_V1R_relative_species <- as.data.frame(NULL)
for(curr_sp in list_sp_RNA){
  
  curr_Summary_RNAseq_table_V1R_relative_species <- 
    as.data.frame(
      Summary_RNAseq_table %>%
        filter(gene_family == "V1R") %>%
        filter(Sp == curr_sp) %>%
        group_by(Sample) %>%
        mutate(total_reads = sum(read_count_gene)) %>%
        ungroup() %>%
        group_by(Sample) %>%
        mutate(All_TPM = sum(TPM)) %>%
        ungroup() %>%
        group_by(gene_family, gene_type, clade, subfamily, gene_name, OGG, kbp) %>%
        summarise(weighted_TPM = weighted.mean(TPM,total_reads))
    )
  
  
  curr_sum_weighted_TPM <- 
    curr_Summary_RNAseq_table_V1R_relative_species %>% pull(weighted_TPM) %>% sum()
  
  curr_Summary_RNAseq_table_V1R_relative_species <- 
    curr_Summary_RNAseq_table_V1R_relative_species %>% 
    mutate(relative_expression = weighted_TPM/curr_sum_weighted_TPM) %>%
    mutate(Sp = curr_sp)
  
  
  Summary_RNAseq_table_V1R_relative_species <- 
    rbind(Summary_RNAseq_table_V1R_relative_species,
          curr_Summary_RNAseq_table_V1R_relative_species)
  
}



#### Compute RE of TAAR genes  ---------------------------------

#At the individual level
Summary_RNAseq_table_TAAR_relative_individual <- 
  as.data.frame(
    Summary_RNAseq_table %>%
      filter(gene_family == "TAAR") %>%
      group_by(Sample) %>%
      mutate(total_reads = sum(read_count_gene)) %>%
      mutate(All_TPM = sum(TPM)) %>%
      mutate(relative_expression = TPM / All_TPM))



#At the species level
list_sp_RNA <- Summary_RNAseq_table %>% pull(Sp) %>% unique()
Summary_RNAseq_table_TAAR_relative_species <- as.data.frame(NULL)
for(curr_sp in list_sp_RNA){
  
  curr_Summary_RNAseq_table_TAAR_relative_species <- 
    as.data.frame(
      Summary_RNAseq_table %>%
        filter(gene_family == "TAAR") %>%
        filter(Sp == curr_sp) %>%
        group_by(Sample) %>%
        mutate(total_reads = sum(read_count_gene)) %>%
        ungroup() %>%
        group_by(Sample) %>%
        mutate(All_TPM = sum(TPM)) %>%
        ungroup() %>%
        group_by(gene_family, gene_type, clade, subfamily, gene_name, OGG, kbp) %>%
        summarise(weighted_TPM = weighted.mean(TPM,total_reads))
    )
  
  
  curr_sum_weighted_TPM <- 
    curr_Summary_RNAseq_table_TAAR_relative_species %>% pull(weighted_TPM) %>% sum()
  
  curr_Summary_RNAseq_table_TAAR_relative_species <- 
    curr_Summary_RNAseq_table_TAAR_relative_species %>% 
    mutate(relative_expression = weighted_TPM/curr_sum_weighted_TPM) %>%
    mutate(Sp = curr_sp)
  
  
  Summary_RNAseq_table_TAAR_relative_species <- 
    rbind(Summary_RNAseq_table_TAAR_relative_species,
          curr_Summary_RNAseq_table_TAAR_relative_species)
  
}


#### Compute RE of V2R genes  ---------------------------------

#At the individual level
Summary_RNAseq_table_V2R_relative_individual <- 
  as.data.frame(
    Summary_RNAseq_table %>%
      filter(gene_family == "V2R") %>%
      group_by(Sample) %>%
      mutate(total_reads = sum(read_count_gene)) %>%
      mutate(All_TPM = sum(TPM)) %>%
      mutate(relative_expression = TPM / All_TPM))

Summary_RNAseq_table_V2R_relative_individual %>%
  filter(Sample == "VGB5Neopul") %>%
  pull(relative_expression) %>% sum()


#At the species level
list_sp_RNA <- Summary_RNAseq_table %>% pull(Sp) %>% unique()
Summary_RNAseq_table_V2R_relative_species <- as.data.frame(NULL)
for(curr_sp in list_sp_RNA){
  
  curr_Summary_RNAseq_table_V2R_relative_species <- 
    as.data.frame(
      Summary_RNAseq_table %>%
        filter(gene_family == "V2R") %>%
        filter(Sp == curr_sp) %>%
        group_by(Sample) %>%
        mutate(total_reads = sum(read_count_gene)) %>%
        ungroup() %>%
        group_by(Sample) %>%
        mutate(All_TPM = sum(TPM)) %>%
        ungroup() %>%
        group_by(gene_family, gene_type, clade, subfamily, gene_name, OGG, kbp) %>%
        summarise(weighted_TPM = weighted.mean(TPM,total_reads))
    )
  
  
  curr_sum_weighted_TPM <- 
    curr_Summary_RNAseq_table_V2R_relative_species %>% pull(weighted_TPM) %>% sum()
  
  curr_Summary_RNAseq_table_V2R_relative_species <- 
    curr_Summary_RNAseq_table_V2R_relative_species %>% 
    mutate(relative_expression = weighted_TPM/curr_sum_weighted_TPM) %>%
    mutate(Sp = curr_sp)
  
  
  Summary_RNAseq_table_V2R_relative_species <- 
    rbind(Summary_RNAseq_table_V2R_relative_species,
          curr_Summary_RNAseq_table_V2R_relative_species)
  
}


Summary_RNAseq_table_V2R_relative_species %>%
  filter(Sp == "Neomul") %>%
  pull(relative_expression) %>%
  sum()


#### Correlation gene expression between individuals  ---------------------------------

#tout dabord faire un test

RE_per_gene_long <- 
  Summary_RNAseq_table_OLR_relative_individual %>%
  dplyr::select(Sample, gene_family, gene_name, relative_expression)
  
RE_per_gene_wide_df <-
  as.data.frame(RE_per_gene_long %>% 
  pivot_wider(names_from = c(Sample), values_from = relative_expression), values_fill = 0)


cor.test(RE_per_gene_wide_df %>% pull(VGC2Asplep), 
         RE_per_gene_wide_df %>% pull(VGC1Asplep), 
         method="pearson")



RE_per_gene_wide_matrix <- 
  RE_per_gene_wide_df %>%
  dplyr::select(-c(gene_family, gene_name))
RE_per_gene_wide_matrix <- as.matrix(RE_per_gene_wide_matrix)


heatmap(RE_per_gene_wide_matrix)


#Lets compute a correlation matrix

cormat_RE <- round(cor(RE_per_gene_wide_matrix, method="pearson"),2)
melted_cormat_RE <- melt(cormat_RE)

# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}


sample_sp_df <- 
  read_count %>% 
  dplyr::select(Sample, Sp) %>% 
  distinct() %>% 
  arrange(Sp)


order_samples <- sample_sp_df$Sample
  



upper_tri <- get_upper_tri(cormat_RE)


upper_tri <- upper_tri[order_samples,,drop=FALSE]
upper_tri <- upper_tri[ , order_samples]



melted_cormat <- melt(upper_tri, na.rm = TRUE)






ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()




#### pGLS RE vs PCG - OLR  --------------------

#Compute the proportion of genes represented by each family
Cichlids_chemoreception_df <- 
  Cichlids_chemoreception_df %>%
  mutate(Total_OLR = Total_OR + Total_TAAR + Total_V1R + Total_V2R) %>%
  mutate(PCG_OR = Total_OR/Total_OLR) %>%
  mutate(PCG_TAAR = Total_TAAR/Total_OLR) %>%
  mutate(PCG_V1R = Total_V1R/Total_OLR) %>%
  mutate(PCG_V2R = Total_V2R/Total_OLR) 


#Sum the RE of each family
Summary_RNAseq_table_OLR_relative_species_fam <- 
  as.data.frame(
    Summary_RNAseq_table_OLR_relative_species %>%
      group_by(Sp, gene_family) %>%
      summarise(fam_relative_expr = sum(relative_expression),
                Sp) %>%
      ungroup() %>%
      distinct())


#Transform the table from long to wide
Summary_RNAseq_table_OLR_relative_species_fam_wide <- 
  as.data.frame(
    pivot_wider(Summary_RNAseq_table_OLR_relative_species_fam,
                names_from = c(gene_family), 
                values_from = fam_relative_expr))
colnames(Summary_RNAseq_table_OLR_relative_species_fam_wide) <- 
  c("Sp", "relative_exp_OR", "relative_exp_TAAR", "relative_exp_V1R", "relative_exp_V2R")


#Combine PCG and RE tables
OLR_RE_PCG_df <- left_join(Summary_RNAseq_table_OLR_relative_species_fam_wide, 
                                          Cichlids_chemoreception_df, 
                                          by="Sp")


#Create a caper data
caper_OLR_RE_PCG <-
  comparative.data(phy = radiation_tree_wo_Neospl, 
                   data = OLR_RE_PCG_df,
                   names.col = Sp, vcv = TRUE,
                   na.omit = FALSE, warn.dropped = TRUE)


#First test the RE vs PCG of OR genes

fit_phylo_OR_expression <- pgls(relative_exp_OR ~ PCG_OR,
                                data = caper_OLR_RE_PCG, 
                                lambda = "ML")

sum_fit_phy <- summary(fit_phylo_OR_expression)
PGLS_pente =   formatC(sum_fit_phy$coefficients[2], digits = 3)
PGLS_r2 =      formatC( sum_fit_phy$r.squared, digits = 2)
PGLS_pvalue =  formatC(sum_fit_phy$coefficients[8], digits = 3)
if (PGLS_pvalue == "   0"){ PGLS_pvalue = "< 2.2e-16"}
PGLS_cc <- coef(fit_phylo_OR_expression)
PGLS_fit_function <- function(x) PGLS_cc[1] + PGLS_cc[2]*x
PGLS_lambda <- sum_fit_phy$param[2]
lm_fit <- lm(data = OLR_RE_PCG_df, formula = relative_exp_OR ~ PCG_OR)
sum_fit <- summary(lm_fit)
GLS_r2 =      formatC( sum_fit$r.squared, digits = 3)
GLS_pente =   formatC(sum_fit$coefficients[2], digits = 3)
GLS_pvalue =  formatC(sum_fit$coefficients[8], digits = 3)
GLS_cc <- coef(lm_fit)
GLS_fit_function <- function(x) GLS_cc[1] + GLS_cc[2]*x


OLR_RE_PCG_df %>% 
  filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>%
  ggplot(aes( x = PCG_OR, y= relative_exp_OR)) +
  geom_point(aes(color=Tribe)) +
  stat_function(fun = GLS_fit_function, color="black") +
  labs(subtitle = bquote(R^2 ~ "=" ~ .(PGLS_r2) ~ ";" ~ P ~ "=" ~.(PGLS_pvalue) ~ ";" ~ lambda ~ "=" ~ .(PGLS_lambda))) +
  theme_classic() +
  xlab("Relative OR repertoire size") +
  ylab("Relative OR expression") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none") + 
  scale_color_manual(values=tribes.colors)



#RE vs PCG of TAAR genes

fit_phylo_TAAR_expression <- pgls(relative_exp_TAAR ~ PCG_TAAR,
                                  data = caper_OLR_RE_PCG, 
                                  lambda = "ML")

sum_fit_phy <- summary(fit_phylo_TAAR_expression)
PGLS_pente =   formatC(sum_fit_phy$coefficients[2], digits = 3)
PGLS_r2 =      formatC( sum_fit_phy$r.squared, digits = 2)
PGLS_pvalue =  formatC(sum_fit_phy$coefficients[8], digits = 3)
if (PGLS_pvalue == "   0"){ PGLS_pvalue = "< 2.2e-16"}
PGLS_cc <- coef(fit_phylo_TAAR_expression)
PGLS_fit_function <- function(x) PGLS_cc[1] + PGLS_cc[2]*x
PGLS_lambda <- sum_fit_phy$param[2]
lm_fit <- lm(data = OLR_RE_PCG_df, formula = relative_exp_TAAR ~ PCG_TAAR)
sum_fit <- summary(lm_fit)
GLS_r2 =      formatC( sum_fit$r.squared, digits = 3)
GLS_pente =   formatC(sum_fit$coefficients[2], digits = 3)
GLS_pvalue =  formatC(sum_fit$coefficients[8], digits = 3)
GLS_cc <- coef(lm_fit)
GLS_fit_function <- function(x) GLS_cc[1] + GLS_cc[2]*x


OLR_RE_PCG_df %>% 
  filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>%
  ggplot(aes( x = PCG_TAAR, y= relative_exp_TAAR)) +
  geom_point(aes(color=Tribe)) +
  stat_function(fun = GLS_fit_function, color="black") +
  labs(subtitle = bquote(R^2 ~ "=" ~ .(PGLS_r2) ~ ";" ~ P ~ "=" ~.(PGLS_pvalue) ~ ";" ~ lambda ~ "=" ~ .(PGLS_lambda))) +
  theme_classic() +
  xlab("Relative TAAR repertoire size") +
  ylab("Relative TAAR expression") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none") + 
  scale_color_manual(values=tribes.colors)



#RE vs PCG of V1R genes

fit_phylo_V1R_expression <- pgls(relative_exp_V1R ~ PCG_V1R,
                                 data = caper_OLR_RE_PCG, 
                                 lambda = "ML")

sum_fit_phy <- summary(fit_phylo_V1R_expression)
PGLS_pente =   formatC(sum_fit_phy$coefficients[2], digits = 3)
PGLS_r2 =      formatC( sum_fit_phy$r.squared, digits = 2)
PGLS_pvalue =  formatC(sum_fit_phy$coefficients[8], digits = 3)
if (PGLS_pvalue == "   0"){ PGLS_pvalue = "< 2.2e-16"}
PGLS_cc <- coef(fit_phylo_V1R_expression)
PGLS_fit_function <- function(x) PGLS_cc[1] + PGLS_cc[2]*x
PGLS_lambda <- sum_fit_phy$param[2]
lm_fit <- lm(data = OLR_RE_PCG_df, formula = relative_exp_V1R ~ PCG_V1R)
sum_fit <- summary(lm_fit)
GLS_r2 =      formatC( sum_fit$r.squared, digits = 3)
GLS_pente =   formatC(sum_fit$coefficients[2], digits = 3)
GLS_pvalue =  formatC(sum_fit$coefficients[8], digits = 3)
GLS_cc <- coef(lm_fit)
GLS_fit_function <- function(x) GLS_cc[1] + GLS_cc[2]*x


OLR_RE_PCG_df %>% 
  filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>%
  ggplot(aes( x = PCG_V1R, y= relative_exp_V1R)) +
  geom_point(aes(color=Tribe)) +
  stat_function(fun = GLS_fit_function, color="black") +
  labs(subtitle = bquote(R^2 ~ "=" ~ .(PGLS_r2) ~ ";" ~ P ~ "=" ~.(PGLS_pvalue) ~ ";" ~ lambda ~ "=" ~ .(PGLS_lambda))) +
  theme_classic() +
  xlab("Relative V1R repertoire size") +
  ylab("Relative V1R expression") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none") + 
  scale_color_manual(values=tribes.colors)



#RE vs PCG of V2R genes

fit_phylo_V2R_expression <- pgls(relative_exp_V2R ~ PCG_V2R,
                                 data = caper_OLR_RE_PCG, 
                                 lambda = "ML")

sum_fit_phy <- summary(fit_phylo_V2R_expression)
PGLS_pente =   formatC(sum_fit_phy$coefficients[2], digits = 3)
PGLS_r2 =      formatC( sum_fit_phy$r.squared, digits = 2)
PGLS_pvalue =  formatC(sum_fit_phy$coefficients[8], digits = 3)
if (PGLS_pvalue == "   0"){ PGLS_pvalue = "< 2.2e-16"}
PGLS_cc <- coef(fit_phylo_V2R_expression)
PGLS_fit_function <- function(x) PGLS_cc[1] + PGLS_cc[2]*x
PGLS_lambda <- sum_fit_phy$param[2]
lm_fit <- lm(data = OLR_RE_PCG_df, formula = relative_exp_V2R ~ PCG_V2R)
sum_fit <- summary(lm_fit)
GLS_r2 =      formatC( sum_fit$r.squared, digits = 3)
GLS_pente =   formatC(sum_fit$coefficients[2], digits = 3)
GLS_pvalue =  formatC(sum_fit$coefficients[8], digits = 3)
GLS_cc <- coef(lm_fit)
GLS_fit_function <- function(x) GLS_cc[1] + GLS_cc[2]*x


OLR_RE_PCG_df %>% 
  filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>%
  ggplot(aes( x = PCG_V2R, y= relative_exp_V2R)) +
  geom_point(aes(color=Tribe)) +
  stat_function(fun = GLS_fit_function, color="black") +
  labs(subtitle = bquote(R^2 ~ "=" ~ .(PGLS_r2) ~ ";" ~ P ~ "=" ~.(PGLS_pvalue) ~ ";" ~ lambda ~ "=" ~ .(PGLS_lambda))) +
  theme_classic() +
  xlab("Relative V2R repertoire size") +
  ylab("Relative V2R expression") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none") + 
  scale_color_manual(values=tribes.colors)


#### pGLS RE vs PCG - OR --------------------

#First import the detailed OR table
OR_DoC_df_info <- read.table("OR_DoC_df_info.csv", header=TRUE,  sep=",")
OR_DoC_df_info <- OR_DoC_df_info %>% filter(! ID %in% c("LBE2", "LBE9"))
OR_DoC_df_info_mean <- read.table("OR_DoC_df_info_mean.csv", header=TRUE, sep=",")
OR_DoC_df_info_mean_total <- OR_DoC_df_info_mean %>% filter(Subfamily == "Total_OR") 
OR_DoC_df_info_mean_wo_total <- 
  OR_DoC_df_info_mean %>% filter(Subfamily != "Total_OR") 

#Compute the PCG of each subfamily
total_gene_df <- 
  as.data.frame(
    OR_DoC_df_info_mean_wo_total %>%
      group_by(Sp) %>%
      summarise(total_gene_nb = sum(mean_normalized_nb)))

OR_DoC_df_info_mean_wo_total <- 
  left_join(OR_DoC_df_info_mean_wo_total, total_gene_df, by="Sp")

OR_DoC_df_info_mean_wo_total <- 
  OR_DoC_df_info_mean_wo_total %>%
  mutate(prop_gene_nb = mean_normalized_nb/total_gene_nb)


#Transform to a wide dataframe
OR_DoC_df_info_mean_wo_total_wide <- 
  as.data.frame(
    pivot_wider(OR_DoC_df_info_mean_wo_total %>%
                  dplyr::select(Sp, Subfamily, prop_gene_nb),
                names_from = c(Subfamily), 
                values_from = prop_gene_nb))


#Sum the RE per OR subfamily
Summary_RNAseq_table_OR_relative_species_subfam <- 
  as.data.frame(
    Summary_RNAseq_table_OR_relative_species %>%
      group_by(Sp, subfamily) %>%
      summarise(subfam_relative_expr = sum(relative_expression),
                Sp) %>%
      ungroup() %>%
      distinct())

Summary_RNAseq_table_OR_relative_species_subfam_wide <- 
  as.data.frame(
    pivot_wider(Summary_RNAseq_table_OR_relative_species_subfam,
                names_from = c(subfamily), 
                values_from = subfam_relative_expr))


#Combine the PCG and RE dataframes
OR_RE_PCG_df <- 
  left_join(Summary_RNAseq_table_OR_relative_species_subfam_wide, 
            OR_DoC_df_info_mean_wo_total_wide,
            by="Sp", 
            suffix = c(".RE", ".PCG"))


#Replace "-" by "_" for convenience in the rest of the script
names(OR_RE_PCG_df) <- 
  gsub(x = names(OR_RE_PCG_df), 
       pattern = "-", replacement = "_") 



#Create a caper data for pGLS computations
caper_OR_RE_PCG <-
  comparative.data(phy = radiation_tree_wo_Neospl, 
                   data = OR_RE_PCG_df,
                   names.col = Sp, vcv = TRUE,
                   na.omit = FALSE, warn.dropped = TRUE)



## Lets perform PGLS between RE and PCG for each subfamily

list_OR_subfamilies <- 
  Summary_RNAseq_table_OR_relative_species_subfam %>% 
  pull(subfamily) %>%
  unique()
list_OR_subfamilies <-
  gsub(x = list_OR_subfamilies, pattern = "-", replacement = "_")



OR_RE_vs_PCG_df <- as.data.frame(NULL)


for(subfam in list_OR_subfamilies){
  
  expression_col <- paste(subfam, "RE", sep=".")
  gene_nb_col <- paste(subfam, "PCG", sep=".")
  
  curr_OR_RE_PCG_df <- 
    OR_RE_PCG_df %>%
    dplyr::select(Sp, expression_col, gene_nb_col)
  
  colnames(curr_OR_RE_PCG_df) <- c("Sp", "expression", "proportion")
  
  curr_caper_OR_RE_PCG <-
    comparative.data(phy = radiation_tree_wo_Neospl, 
                     data = curr_OR_RE_PCG_df,
                     names.col = Sp, vcv = TRUE,
                     na.omit = FALSE, warn.dropped = TRUE)
  
  fit_phylo_OR_expression <-
    pgls(expression ~ proportion, 
         data = curr_caper_OR_RE_PCG, 
         lambda = "ML", 
         bounds=list(lambda=c(0,1)))
  
  
  sum_fit_phy <- summary(fit_phylo_OR_expression)
  PGLS_pente =   formatC(sum_fit_phy$coefficients[2], digits = 3)
  PGLS_r2 =      formatC( sum_fit_phy$r.squared, digits = 2)
  PGLS_pvalue =  formatC(sum_fit_phy$coefficients[8], digits = 3)
  if (PGLS_pvalue == "   0"){ PGLS_pvalue = "< 2.2e-16"}
  PGLS_cc <- coef(fit_phylo_OR_expression)
  PGLS_fit_function <- function(x) PGLS_cc[1] + PGLS_cc[2]*x
  PGLS_lambda <- sum_fit_phy$param[2]
  PGLS_r2
  lm_fit <- lm(data = curr_OR_RE_PCG_df, 
               formula = expression ~ proportion)
  sum_fit <- summary(lm_fit)
  GLS_r2 =      formatC( sum_fit$r.squared, digits = 3)
  GLS_pente =   formatC(sum_fit$coefficients[2], digits = 3)
  GLS_pvalue =  formatC(sum_fit$coefficients[8], digits = 3)
  GLS_cc <- coef(lm_fit)
  GLS_fit_function <- function(x) GLS_cc[1] + GLS_cc[2]*x
  
  xlab_curr <- paste("Proportion of", subfam, "genes")
  ylab_curr <- paste(subfam, "relative expression")
  
  
  curr_df <- as.data.frame(cbind(subfam, PGLS_pente, PGLS_r2, PGLS_pvalue, PGLS_lambda))
  colnames(curr_df) <- c("Subfamily", "PGLS_pente","PGLS_R2", "PGLS_Pvalue", "PGLS_Lambda")
  OR_RE_vs_PCG_df <- rbind(OR_RE_vs_PCG_df, curr_df)
  
  curr_OR_RE_PCG_df <- 
    left_join(curr_OR_RE_PCG_df, Tribe_Sp_df, by="Sp")
  
  g1 <- 
    curr_OR_RE_PCG_df %>% 
    ggplot(aes( x = proportion, y= expression, color=Tribe)) +
    geom_point()+
    stat_function(fun = GLS_fit_function, color="black")+
    labs(subtitle = bquote(R^2 ~ "=" ~ .(PGLS_r2) ~ ";" ~ P ~ "=" ~.(PGLS_pvalue) ~ ";" ~ lambda ~ "=" ~ .(PGLS_lambda))) +
    theme_classic() +
    xlab(xlab_curr) +
    ylab(ylab_curr) +
    theme(axis.text=element_text(size=16),
          axis.title=element_text(size=18),
          plot.subtitle=element_text(size=16),
          legend.position="none") + 
    scale_color_manual(values=tribes.colors)
  
  
  #print(g1)
  
  
  
}


nrow(OR_RE_vs_PCG_df %>%
       filter(as.numeric(PGLS_Pvalue) > 0.05))

nrow(OR_RE_vs_PCG_df %>%
       filter(as.numeric(PGLS_Pvalue) <= 0.05))



#### pGLS RE vs PCG - TAAR --------------------

#First import the detailed TAAR table
TAAR_DoC_df_info <- read.table("TAAR_DoC_df_info.csv", header=TRUE,  sep=",")
TAAR_DoC_df_info <- TAAR_DoC_df_info %>% filter(! ID %in% c("LBE2", "LBE9"))
TAAR_DoC_df_info_mean <- read.table("TAAR_DoC_df_info_mean.csv", header=TRUE, sep=",")
TAAR_DoC_df_info_mean_total <- TAAR_DoC_df_info_mean %>% filter(Subfamily == "Total_OR") 
TAAR_DoC_df_info_mean_wo_total <- 
  TAAR_DoC_df_info_mean %>% filter(Subfamily != "Total_OR") 

#Compute the PCG of each subfamily
total_gene_df <- 
  as.data.frame(
    TAAR_DoC_df_info_mean_wo_total %>%
      group_by(Sp) %>%
      summarise(total_gene_nb = sum(mean_normalized_nb)))

TAAR_DoC_df_info_mean_wo_total <- 
  left_join(TAAR_DoC_df_info_mean_wo_total, total_gene_df, by="Sp")

TAAR_DoC_df_info_mean_wo_total <- 
  TAAR_DoC_df_info_mean_wo_total %>%
  mutate(prop_gene_nb = mean_normalized_nb/total_gene_nb)


#Transform to a wide dataframe
TAAR_DoC_df_info_mean_wo_total_wide <- 
  as.data.frame(
    pivot_wider(TAAR_DoC_df_info_mean_wo_total %>%
                  dplyr::select(Sp, Subfamily, prop_gene_nb),
                names_from = c(Subfamily), 
                values_from = prop_gene_nb))


#Sum the RE per TAAR subfamily
Summary_RNAseq_table_TAAR_relative_species_subfam <- 
  as.data.frame(
    Summary_RNAseq_table_TAAR_relative_species %>%
      group_by(Sp, subfamily) %>%
      summarise(subfam_relative_expr = sum(relative_expression),
                Sp) %>%
      ungroup() %>%
      distinct())

Summary_RNAseq_table_TAAR_relative_species_subfam_wide <- 
  as.data.frame(
    pivot_wider(Summary_RNAseq_table_TAAR_relative_species_subfam,
                names_from = c(subfamily), 
                values_from = subfam_relative_expr))


#Combine the PCG and RE dataframes
TAAR_RE_PCG_df <- 
  left_join(Summary_RNAseq_table_TAAR_relative_species_subfam_wide, 
            TAAR_DoC_df_info_mean_wo_total_wide,
            by="Sp", 
            suffix = c(".RE", ".PCG"))


#Replace "-" by "_" for convenience in the rest of the script
names(TAAR_RE_PCG_df) <- 
  gsub(x = names(TAAR_RE_PCG_df), 
       pattern = "-", replacement = "_") 



#Create a caper data for pGLS computations
caper_TAAR_RE_PCG <-
  comparative.data(phy = radiation_tree_wo_Neospl, 
                   data = TAAR_RE_PCG_df,
                   names.col = Sp, vcv = TRUE,
                   na.omit = FALSE, warn.dropped = TRUE)



## Lets perform PGLS between RE and PCG for each subfamily

list_TAAR_subfamilies <- 
  Summary_RNAseq_table_TAAR_relative_species_subfam %>% 
  pull(subfamily) %>%
  unique()
list_TAAR_subfamilies <-
  gsub(x = list_TAAR_subfamilies, pattern = "-", replacement = "_")



TAAR_RE_vs_PCG_df <- as.data.frame(NULL)


for(subfam in list_TAAR_subfamilies){
  
  expression_col <- paste(subfam, "RE", sep=".")
  gene_nb_col <- paste(subfam, "PCG", sep=".")
  
  curr_TAAR_RE_PCG_df <- 
    TAAR_RE_PCG_df %>%
    dplyr::select(Sp, expression_col, gene_nb_col)
  
  colnames(curr_TAAR_RE_PCG_df) <- c("Sp", "expression", "proportion")
  
  curr_caper_TAAR_RE_PCG <-
    comparative.data(phy = radiation_tree_wo_Neospl, 
                     data = curr_TAAR_RE_PCG_df,
                     names.col = Sp, vcv = TRUE,
                     na.omit = FALSE, warn.dropped = TRUE)
  
  fit_phylo_TAAR_expression <-
    pgls(expression ~ proportion, 
         data = curr_caper_TAAR_RE_PCG, 
         lambda = "ML", 
         bounds=list(lambda=c(0,1)))
  
  
  sum_fit_phy <- summary(fit_phylo_TAAR_expression)
  PGLS_pente =   formatC(sum_fit_phy$coefficients[2], digits = 3)
  PGLS_r2 =      formatC( sum_fit_phy$r.squared, digits = 2)
  PGLS_pvalue =  formatC(sum_fit_phy$coefficients[8], digits = 3)
  if (PGLS_pvalue == "   0"){ PGLS_pvalue = "< 2.2e-16"}
  PGLS_cc <- coef(fit_phylo_TAAR_expression)
  PGLS_fit_function <- function(x) PGLS_cc[1] + PGLS_cc[2]*x
  PGLS_lambda <- sum_fit_phy$param[2]
  PGLS_r2
  lm_fit <- lm(data = curr_TAAR_RE_PCG_df, 
               formula = expression ~ proportion)
  sum_fit <- summary(lm_fit)
  GLS_r2 =      formatC( sum_fit$r.squared, digits = 3)
  GLS_pente =   formatC(sum_fit$coefficients[2], digits = 3)
  GLS_pvalue =  formatC(sum_fit$coefficients[8], digits = 3)
  GLS_cc <- coef(lm_fit)
  GLS_fit_function <- function(x) GLS_cc[1] + GLS_cc[2]*x
  
  xlab_curr <- paste("Proportion of", subfam, "genes")
  ylab_curr <- paste(subfam, "relative expression")
  
  
  curr_df <- as.data.frame(cbind(subfam, PGLS_pente, PGLS_r2, PGLS_pvalue, PGLS_lambda))
  colnames(curr_df) <- c("Subfamily", "PGLS_pente","PGLS_R2", "PGLS_Pvalue", "PGLS_Lambda")
  TAAR_RE_vs_PCG_df <- rbind(TAAR_RE_vs_PCG_df, curr_df)
  
  curr_TAAR_RE_PCG_df <- 
    left_join(curr_TAAR_RE_PCG_df, Tribe_Sp_df, by="Sp")
  
  g1 <- 
    curr_TAAR_RE_PCG_df %>% 
    ggplot(aes( x = proportion, y= expression, color=Tribe)) +
    geom_point()+
    stat_function(fun = GLS_fit_function, color="black")+
    labs(subtitle = bquote(R^2 ~ "=" ~ .(PGLS_r2) ~ ";" ~ P ~ "=" ~.(PGLS_pvalue) ~ ";" ~ lambda ~ "=" ~ .(PGLS_lambda))) +
    theme_classic() +
    xlab(xlab_curr) +
    ylab(ylab_curr) +
    theme(axis.text=element_text(size=16),
          axis.title=element_text(size=18),
          plot.subtitle=element_text(size=16),
          legend.position="none") + 
    scale_color_manual(values=tribes.colors)
  
  
  #print(g1)
  
  
  
}


nrow(TAAR_RE_vs_PCG_df %>%
       filter(as.numeric(PGLS_Pvalue) > 0.05))

nrow(TAAR_RE_vs_PCG_df %>%
       filter(as.numeric(PGLS_Pvalue) <= 0.05))




#### pGLS RE vs PCG - V2R --------------------

#First import the detailed V2R table
V2R_DoC_df_info <- read.table("V2R_DoC_df_info.csv", header=TRUE,  sep=",")
V2R_DoC_df_info <- V2R_DoC_df_info %>% filter(! ID %in% c("LBE2", "LBE9"))
V2R_DoC_df_info_mean <- read.table("V2R_DoC_df_info_mean.csv", header=TRUE, sep=",")
V2R_DoC_df_info_mean_total <- V2R_DoC_df_info_mean %>% filter(Subfamily == "Total_OR") 
V2R_DoC_df_info_mean_wo_total <- 
  V2R_DoC_df_info_mean %>% filter(Subfamily != "Total_OR") 

#Compute the PCG of each subfamily
total_gene_df <- 
  as.data.frame(
    V2R_DoC_df_info_mean_wo_total %>%
      group_by(Sp) %>%
      summarise(total_gene_nb = sum(mean_normalized_nb)))

V2R_DoC_df_info_mean_wo_total <- 
  left_join(V2R_DoC_df_info_mean_wo_total, total_gene_df, by="Sp")

V2R_DoC_df_info_mean_wo_total <- 
  V2R_DoC_df_info_mean_wo_total %>%
  mutate(prop_gene_nb = mean_normalized_nb/total_gene_nb)


#Transform to a wide dataframe
V2R_DoC_df_info_mean_wo_total_wide <- 
  as.data.frame(
    pivot_wider(V2R_DoC_df_info_mean_wo_total %>%
                  dplyr::select(Sp, Subfamily, prop_gene_nb),
                names_from = c(Subfamily), 
                values_from = prop_gene_nb))


#Sum the RE per V2R subfamily
Summary_RNAseq_table_V2R_relative_species_subfam <- 
  as.data.frame(
    Summary_RNAseq_table_V2R_relative_species %>%
      group_by(Sp, subfamily) %>%
      summarise(subfam_relative_expr = sum(relative_expression),
                Sp) %>%
      ungroup() %>%
      distinct())

Summary_RNAseq_table_V2R_relative_species_subfam_wide <- 
  as.data.frame(
    pivot_wider(Summary_RNAseq_table_V2R_relative_species_subfam,
                names_from = c(subfamily), 
                values_from = subfam_relative_expr))


#Combine the PCG and RE dataframes
V2R_RE_PCG_df <- 
  left_join(Summary_RNAseq_table_V2R_relative_species_subfam_wide, 
            V2R_DoC_df_info_mean_wo_total_wide,
            by="Sp", 
            suffix = c(".RE", ".PCG"))


#Replace "-" by "_" for convenience in the rest of the script
names(V2R_RE_PCG_df) <- 
  gsub(x = names(V2R_RE_PCG_df), 
       pattern = "-", replacement = "_") 



#Create a caper data for pGLS computations
caper_V2R_RE_PCG <-
  comparative.data(phy = radiation_tree_wo_Neospl, 
                   data = V2R_RE_PCG_df,
                   names.col = Sp, vcv = TRUE,
                   na.omit = FALSE, warn.dropped = TRUE)



## Lets perform PGLS between RE and PCG for each subfamily

list_V2R_subfamilies <- 
  Summary_RNAseq_table_V2R_relative_species_subfam %>% 
  pull(subfamily) %>%
  unique()
list_V2R_subfamilies <-
  gsub(x = list_V2R_subfamilies, pattern = "-", replacement = "_")



V2R_RE_vs_PCG_df <- as.data.frame(NULL)


for(subfam in list_V2R_subfamilies){
  
  expression_col <- paste(subfam, "RE", sep=".")
  gene_nb_col <- paste(subfam, "PCG", sep=".")
  
  curr_V2R_RE_PCG_df <- 
    V2R_RE_PCG_df %>%
    dplyr::select(Sp, expression_col, gene_nb_col)
  
  colnames(curr_V2R_RE_PCG_df) <- c("Sp", "expression", "proportion")
  
  curr_caper_V2R_RE_PCG <-
    comparative.data(phy = radiation_tree_wo_Neospl, 
                     data = curr_V2R_RE_PCG_df,
                     names.col = Sp, vcv = TRUE,
                     na.omit = FALSE, warn.dropped = TRUE)
  
  fit_phylo_V2R_expression <-
    pgls(expression ~ proportion, 
         data = curr_caper_V2R_RE_PCG, 
         lambda = "ML", 
         bounds=list(lambda=c(0,1)))
  
  
  sum_fit_phy <- summary(fit_phylo_V2R_expression)
  PGLS_pente =   formatC(sum_fit_phy$coefficients[2], digits = 3)
  PGLS_r2 =      formatC( sum_fit_phy$r.squared, digits = 2)
  PGLS_pvalue =  formatC(sum_fit_phy$coefficients[8], digits = 3)
  if (PGLS_pvalue == "   0"){ PGLS_pvalue = "< 2.2e-16"}
  PGLS_cc <- coef(fit_phylo_V2R_expression)
  PGLS_fit_function <- function(x) PGLS_cc[1] + PGLS_cc[2]*x
  PGLS_lambda <- sum_fit_phy$param[2]
  PGLS_r2
  lm_fit <- lm(data = curr_V2R_RE_PCG_df, 
               formula = expression ~ proportion)
  sum_fit <- summary(lm_fit)
  GLS_r2 =      formatC( sum_fit$r.squared, digits = 3)
  GLS_pente =   formatC(sum_fit$coefficients[2], digits = 3)
  GLS_pvalue =  formatC(sum_fit$coefficients[8], digits = 3)
  GLS_cc <- coef(lm_fit)
  GLS_fit_function <- function(x) GLS_cc[1] + GLS_cc[2]*x
  
  xlab_curr <- paste("Proportion of", subfam, "genes")
  ylab_curr <- paste(subfam, "relative expression")
  
  
  curr_df <- as.data.frame(cbind(subfam, PGLS_pente, PGLS_r2, PGLS_pvalue, PGLS_lambda))
  colnames(curr_df) <- c("Subfamily", "PGLS_pente","PGLS_R2", "PGLS_Pvalue", "PGLS_Lambda")
  V2R_RE_vs_PCG_df <- rbind(V2R_RE_vs_PCG_df, curr_df)
  
  curr_V2R_RE_PCG_df <- 
    left_join(curr_V2R_RE_PCG_df, Tribe_Sp_df, by="Sp")
  
  g1 <- 
    curr_V2R_RE_PCG_df %>% 
    ggplot(aes( x = proportion, y= expression, color=Tribe)) +
    geom_point()+
    stat_function(fun = GLS_fit_function, color="black")+
    labs(subtitle = bquote(R^2 ~ "=" ~ .(PGLS_r2) ~ ";" ~ P ~ "=" ~.(PGLS_pvalue) ~ ";" ~ lambda ~ "=" ~ .(PGLS_lambda))) +
    theme_classic() +
    xlab(xlab_curr) +
    ylab(ylab_curr) +
    theme(axis.text=element_text(size=16),
          axis.title=element_text(size=18),
          plot.subtitle=element_text(size=16),
          legend.position="none") + 
    scale_color_manual(values=tribes.colors)
  
  
  #print(g1)
  
  
  
}


nrow(V2R_RE_vs_PCG_df %>%
       filter(as.numeric(PGLS_Pvalue) > 0.05))

nrow(V2R_RE_vs_PCG_df %>%
       filter(as.numeric(PGLS_Pvalue) <= 0.05))


#### Summary RE vs PCG results  --------------------

All_RE_vs_PCG_df <- 
  rbind(V2R_RE_vs_PCG_df,
        OR_RE_vs_PCG_df,
        TAAR_RE_vs_PCG_df)

All_RE_vs_PCG_df$PGLS_Pvalue <- as.numeric(All_RE_vs_PCG_df$PGLS_Pvalue)


nrow(All_RE_vs_PCG_df %>% filter(PGLS_Pvalue >= 0.05)) #non significant
nrow(All_RE_vs_PCG_df %>% filter(PGLS_Pvalue < 0.05)) #significant


write.table(All_RE_vs_PCG_df,
            "Raw_R_plots/All_RE_vs_PCG.tsv",
            sep="\t",
            row.names=FALSE,
            col.names=TRUE,
            quote=FALSE)


#### Create color palette for phylogeny  --------------------

Tribe_for_tree <- 
  OR_DoC_df_info_mean_total %>%
  dplyr::select(Sp, Tribe)
Tribe_for_tree[(Tribe_for_tree$Tribe == "Serranochromini"),"Tribe"] <- "Haplochromini"

#create a node number and node label correspondance table

radiation_tree_wo_Neospl <- makeNodeLabel(radiation_tree_wo_Neospl, method="number", prefix="Node")

species_test <- "Neomul"
misc_table <- as.data.frame(species_test) %>% mutate(value = 1)
colnames(misc_table) <- c("node","value")
misc_table$node <- as.numeric(misc_table$node)
node_label_corresp <- 
  left_join(radiation_tree_wo_Neospl, misc_table, by = 'node') %>%
  dplyr::select(node, label)

node_nb <- 
  node_label_corresp %>%
  filter(., grepl("Node", label)) %>%
  pull(node)
species_names <- 
  node_label_corresp %>%
  filter(., !grepl("Node", label)) %>%
  pull(label)
species_and_nodes <- c(species_names, node_nb)
node_label_corresp <- node_label_corresp %>% mutate(tree_names = species_and_nodes)
node_label_corresp <- as.data.frame(node_label_corresp)



list_tribes <- 
  Tribe_for_tree %>% 
  distinct() %>%
  filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>% 
  group_by(Tribe) %>%
  summarise(count = n()) %>%
  filter(count > 1) %>%
  pull(Tribe) %>%
  unique()

df_nodes_tribes <- as.data.frame(NULL)
for(curr_tribe in list_tribes){
  
  curr_nodes <- 
    getMRCA(radiation_tree_wo_Neospl, 
            Tribe_for_tree %>% 
              filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>% 
              filter(Tribe == curr_tribe) %>%
              pull(Sp) %>%
              unique())
  curr_df <-
    as.data.frame(
      getDescendants(radiation_tree_wo_Neospl, node=curr_nodes)) %>%
    mutate(Tribe_curr = curr_tribe)
  
  colnames(curr_df) <- c("node", "Tribe")
  
  df_nodes_tribes <- rbind(df_nodes_tribes, curr_df)
}

Tribe_for_tree <- 
  left_join(node_label_corresp, df_nodes_tribes, by="node") %>%
  dplyr::select(label, Tribe)

#Remove Haplochromnii branches that are shared with tropheini branches

labels_to_remove <- 
  as.data.frame(
    Tribe_for_tree %>%
      filter(Tribe %in% c("Haplochromini", "Tropheini")) %>%
      group_by(label) %>%
      summarise(count = n()) %>%
      filter(count > 1)) %>%
  pull(label)

Tribe_lines_to_remove <- 
  Tribe_for_tree %>%
  filter(label %in% labels_to_remove) %>%
  filter(Tribe == "Haplochromini")


Tribe_for_tree <- anti_join(Tribe_for_tree, Tribe_lines_to_remove, by = c("label", "Tribe"))
Tribe_for_tree[(Tribe_for_tree$label == "Node225"),"Tribe"] <- "Tropheini"

to_add_df <- rbind(
  as.data.frame(cbind("Node203", "Eretmodini")),
  as.data.frame(cbind("Node207", "Haplochromini")),
  as.data.frame(cbind("Node195", "Perissodini")),
  as.data.frame(cbind("Node192", "Benthochromini")),
  as.data.frame(cbind("Node180", "Cyprichromini")),
  as.data.frame(cbind("Node139", "Ectodini")),
  as.data.frame(cbind("Node126", "Cyphotilapiini")),
  as.data.frame(cbind("Node130", "Limnochromini")),
  as.data.frame(cbind("Node19", "Lamprologini")),
  as.data.frame(cbind("Node10", "Bathybatini")),
  as.data.frame(cbind("Node4", "Trematocarini")),
  as.data.frame(cbind("Boumic", "Boulengerochromini"))
)
colnames(to_add_df) <- c("label", "Tribe")

Tribe_for_tree <- anti_join(Tribe_for_tree, to_add_df, by = c("label"))

Tribe_for_tree <- rbind(Tribe_for_tree, to_add_df)




#### Draw a phylogeny + RE of each OLR family + subfamily  --------------------



RNA_seq_tree <- 
  keep.tip(radiation_tree_wo_Neospl, 
           Summary_RNAseq_table_OR_relative_species_subfam %>% 
             filter(Sp != "Oretan") %>% 
             pull(Sp) %>% unique()
  )


Summary_RNAseq_table_V1R_relative_species_subfam <- 
  Summary_RNAseq_table_V1R_relative_species %>%
  dplyr::select(Sp, subfamily, relative_expression)


Summary_RNAseq_table_OR_relative_species_subfam$subfamily <- 
  factor(Summary_RNAseq_table_OR_relative_species_subfam$subfamily, 
         levels=rev(c("Eta-1","Eta-2","Eta-3","Eta-4","Eta-6",
                      "Eta-7","Eta-8","Eta-9","Eta-10","Eta-11","Eta-12",
                      "Eta-13","Eta-14","Eta-15","Beta-1","Beta-2","Epsilon-1","Epsilon-2",
                      "Epsilon-3","Epsilon-4","Zeta-1","Zeta-2","Zeta-3","Zeta-4","Delta-1",
                      "Delta-2","Delta-3","Delta-4","Delta-5","Delta-6","Delta-7","Delta-8",
                      "Delta-9","Delta-10","Delta-11","Delta-12","Delta-13","Delta-14","Delta-15",
                      "Delta-16","Delta-17","Delta-18","Delta-19","Delta-20")))


Summary_RNAseq_table_TAAR_relative_species_subfam$subfamily <- 
  factor(Summary_RNAseq_table_TAAR_relative_species_subfam$subfamily, 
         levels=rev(c('TAARA2-1',"TAARA2-2",'TAARB4-1',"TAARB4-2",
                      "TAARB4-3",'TAARB4-4','TAARB4-5')))


Summary_RNAseq_table_V2R_relative_species_subfam$subfamily <- 
  factor(Summary_RNAseq_table_V2R_relative_species_subfam$subfamily, 
         levels=rev(c("V2RD3-1", "V2RD4-1", "V2RD5-1", "V2RD7-1", 
                      "V2RD8-1", "V2RD9-1", "V2RD10-1", "V2RD11-1", "V2RD11-2", 
                      "V2RD11-3", "V2RD11-4", "V2RD11-5", "V2RD11-6", "V2RD11-7")))

#Draw the backbone first
p <- 
  ggtree(RNA_seq_tree)  %<+% Tribe_for_tree +
  aes(color=Tribe) + 
  scale_color_manual(values = tribes.colors) +
  geom_tiplab(size=3, offset=0.1, fontface=3) + 
  theme(legend.position='none') +
  geom_treescale(fontsize=4, linesize=1, offset=0.2, x=0, y=25, width=2) +
  ggnewscale::new_scale_color() 


#Figure with everything
p1 <- 
  p + 
  geom_facet(panel = 'OLR', 
             data = Summary_RNAseq_table_OLR_relative_species_fam, geom = geom_bar, 
             mapping = aes(x = as.numeric(fam_relative_expr), fill = as.factor(gene_family)), 
             orientation = 'y', width = 0.8, stat='identity', 
             position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values = olfactory_family_colors) +
  ggnewscale::new_scale_fill() +
  geom_facet(panel = 'OR', 
             data = Summary_RNAseq_table_OR_relative_species_subfam, geom = geom_bar, 
             mapping = aes(x = as.numeric(subfam_relative_expr), fill = as.factor(subfamily)), 
             orientation = 'y', width = 0.8, stat='identity', 
             position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values = OR_subfam_colors) +
  ggnewscale::new_scale_fill() +
  geom_facet(panel = 'TAAR', 
             data = Summary_RNAseq_table_TAAR_relative_species_subfam, geom = geom_bar, 
             mapping = aes(x = as.numeric(subfam_relative_expr), fill = as.factor(subfamily)), 
             orientation = 'y', width = 0.8, stat='identity', 
             position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values = TAAR_subfam_colors) +
  ggnewscale::new_scale_fill() +
  geom_facet(panel = 'V1R', 
             data = Summary_RNAseq_table_V1R_relative_species_subfam, geom = geom_bar, 
             mapping = aes(x = as.numeric(relative_expression), fill = as.factor(subfamily)), 
             orientation = 'y', width = 0.8, stat='identity', 
             position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values = V1R_subfam_colors) +
  ggnewscale::new_scale_fill() +
  geom_facet(panel = 'V2R', 
             data = Summary_RNAseq_table_V2R_relative_species_subfam, geom = geom_bar, 
             mapping = aes(x = as.numeric(subfam_relative_expr), fill = as.factor(subfamily)), 
             orientation = 'y', width = 0.8, stat='identity', 
             position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values = V2R_subfam_colors) +
  ggnewscale::new_scale_fill() +
  theme_tree2(legend.position = 'none') +
  xlim_tree(15)


#Figure with only OLR
p2 <- 
  p + 
  geom_facet(panel = 'OLR', 
             data = Summary_RNAseq_table_OLR_relative_species_fam, geom = geom_bar, 
             mapping = aes(x = as.numeric(fam_relative_expr), fill = as.factor(gene_family)), 
             orientation = 'y', width = 0.8, stat='identity', 
             position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values = olfactory_family_colors) +
  ggnewscale::new_scale_fill() +
  theme_tree2(legend.position = 'none') +
  xlim_tree(15)


#Figure with only OR
p3 <- 
  p + 
  geom_facet(panel = 'OR', 
             data = Summary_RNAseq_table_OR_relative_species_subfam, geom = geom_bar, 
             mapping = aes(x = as.numeric(subfam_relative_expr), fill = as.factor(subfamily)), 
             orientation = 'y', width = 0.8, stat='identity', 
             position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values = OR_subfam_colors) +
  ggnewscale::new_scale_fill() +
  theme_tree2(legend.position = 'none') +
  xlim_tree(15)


#Figure with only TAAR
p4 <- 
  p + 
  geom_facet(panel = 'TAAR', 
             data = Summary_RNAseq_table_TAAR_relative_species_subfam, geom = geom_bar, 
             mapping = aes(x = as.numeric(subfam_relative_expr), fill = as.factor(subfamily)), 
             orientation = 'y', width = 0.8, stat='identity', 
             position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values = TAAR_subfam_colors) +
  ggnewscale::new_scale_fill() +
  theme_tree2(legend.position = 'none') +
  xlim_tree(15)


#Figure with only V2R
p5 <- 
  p + 
  geom_facet(panel = 'V2R', 
             data = Summary_RNAseq_table_V2R_relative_species_subfam, geom = geom_bar, 
             mapping = aes(x = as.numeric(subfam_relative_expr), fill = as.factor(subfamily)), 
             orientation = 'y', width = 0.8, stat='identity', 
             position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values = V2R_subfam_colors) +
  ggnewscale::new_scale_fill() +
  theme_tree2(legend.position = 'none') +
  xlim_tree(15)


#Figure with only V1R
p6 <- 
  p + 
  geom_facet(panel = 'V1R', 
             data = Summary_RNAseq_table_V1R_relative_species_subfam, geom = geom_bar, 
             mapping = aes(x = as.numeric(relative_expression), fill = as.factor(subfamily)), 
             orientation = 'y', width = 0.8, stat='identity', 
             position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values = V1R_subfam_colors) +
  ggnewscale::new_scale_fill() +
  theme_tree2(legend.position = 'none') +
  xlim_tree(15)


pdf(file = "Raw_R_plots/Species_tree_Figure_expression.pdf",width = 12.34,  height = 4.61)

facet_widths(p1, c(Tree = 1))
facet_widths(p2, c(Tree = 1))
facet_widths(p3, c(Tree = 1))
facet_widths(p4, c(Tree = 1))
facet_widths(p5, c(Tree = 1))
facet_widths(p6, c(Tree = 1))

dev.off()







#### PCA of olfactory receptors expression --------------------


#Compute the relative expression of each subfamily compared to the whole OLR repertoire
Summary_RNAseq_table_OLR_relative_species_subfam <- 
  as.data.frame(
    Summary_RNAseq_table_OLR_relative_species %>%
      group_by(Sp, subfamily) %>%
      summarise(subfam_relative_expr = sum(relative_expression),
                Sp) %>%
      ungroup() %>%
      distinct())


Summary_test <- 
  as.data.frame(
    Summary_RNAseq_table_OLR_relative_species %>%
      group_by(Sp, subfamily) %>%
      summarise(subfam_relative_expr = sum(weighted_TPM),
                Sp) %>%
      ungroup() %>%
      distinct())

#Transform the long df to a wide df
Summary_RNAseq_table_OLR_relative_species_wide <- 
  as.data.frame(
    pivot_wider(Summary_RNAseq_table_OLR_relative_species_subfam,
                names_from = c(subfamily), 
                values_from = subfam_relative_expr))


Summary_test_wide <- 
  as.data.frame(
    pivot_wider(Summary_test,
                names_from = c(subfamily), 
                values_from = subfam_relative_expr))


#Add tribe info to the table, but keep one without any info for the PCA computation
Summary_RNAseq_table_OLR_relative_species_wide_woinfo <- 
  Summary_RNAseq_table_OLR_relative_species_wide
rownames(Summary_RNAseq_table_OLR_relative_species_wide_woinfo) <- Summary_RNAseq_table_OLR_relative_species_wide_woinfo[,1]
Summary_RNAseq_table_OLR_relative_species_wide_woinfo <- Summary_RNAseq_table_OLR_relative_species_wide_woinfo[,-1]
Summary_RNAseq_table_OLR_relative_species_wide_info <- 
  left_join(Summary_RNAseq_table_OLR_relative_species_wide, Tribe_Sp_df, by="Sp")


rownames(Summary_test_wide) <- Summary_test_wide[,1]
Summary_test_wide <- Summary_test_wide[,-1]

#Compute the PCA ! 
PCA_OLR_expression <- prcomp(Summary_RNAseq_table_OLR_relative_species_wide_woinfo, scale = TRUE)
#PCA_OLR_expression <- prcomp(Summary_test_wide, scale = TRUE)  #this is just to control that the same results is given when using TPM instead of proportions

#autoplot to see the results
autoplot(PCA_OLR_expression, data = Summary_RNAseq_table_OLR_relative_species_wide_woinfo) +
  theme_classic()

PC1 <- PCA_OLR_expression$x[,1]
PC2 <- PCA_OLR_expression$x[,2]


pdf(file = "Raw_R_plots/PCA_olfactory_receptors.pdf",width = 8.34,  height = 4.61)


ggplot(Summary_RNAseq_table_OLR_relative_species_wide_info,
       aes(x = PC1,
           y = PC2,
           color=Tribe)) +
  geom_point(size = 3) +
  scale_color_manual(values = tribes.colors_sec) +
  theme_classic() +
  xlab("PC1 (15.94%)") +
  ylab("PC2 (11.77%)") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none") 


dev.off()


## Check if PC1 and PC2 have a phylogenetic signal


Summary_RNAseq_table_OLR_relative_species_wide_info <- 
  Summary_RNAseq_table_OLR_relative_species_wide_info %>%
  mutate(PC1 = PCA_OLR_expression$x[,1]) %>%
  mutate(PC2 = PCA_OLR_expression$x[,2])


Sp_pc1 <- 
  Summary_RNAseq_table_OLR_relative_species_wide_info %>%
  dplyr::select(Sp, PC1)
Sp_pc1 <- 
  Sp_pc1 %>%
  filter(Sp %in% radiation_tree_wo_Neospl$tip.label)
tip_order <- radiation_tree_wo_Neospl$tip.label
Sp_pc1_ordered <- 
  enframe(tip_order, name = NULL, value = "Sp") %>% 
  left_join(Sp_pc1) %>% 
  column_to_rownames("Sp")

Sp_pc1_ordered <- Sp_pc1_ordered %>% filter(! is.na(PC1))
Sp_pc1_tree <- keep.tip(radiation_tree_wo_Neospl, row.names(Sp_pc1_ordered))
phylosig(tree = Sp_pc1_tree, x = Sp_pc1_ordered$PC1, method = "lambda", test = T)




Sp_pc2 <- 
  Summary_RNAseq_table_OLR_relative_species_wide_info %>%
  dplyr::select(Sp, PC2)
Sp_pc2 <- 
  Sp_pc2 %>%
  filter(Sp %in% radiation_tree_wo_Neospl$tip.label)
tip_order <- radiation_tree_wo_Neospl$tip.label
Sp_pc2_ordered <- 
  enframe(tip_order, name = NULL, value = "Sp") %>% 
  left_join(Sp_pc2) %>% 
  column_to_rownames("Sp")
Sp_pc2_ordered <- Sp_pc2_ordered %>% filter(! is.na(PC2))
Sp_pc2_tree <- keep.tip(radiation_tree_wo_Neospl, row.names(Sp_pc2_ordered))
phylosig(tree = Sp_pc2_tree, x = Sp_pc2_ordered$PC2, method = "lambda", test = T)

#Finally,  test the significance of the PCA

#Run it only if necessary, takes a lot of time

##test_PCA_signif <- 
##  PCAtest(Summary_RNAseq_table_OLR_relative_species_wide_woinfo, varcorr=FALSE, counter=FALSE, plot=FALSE)
##
##PCAtest(Summary_RNAseq_table_OLR_relative_species_wide_woinfo, varcorr=FALSE, counter=FALSE, plot=TRUE)
##
##
##PCA1_signif_columns <- c(1, 3, 4, 5, 6, 8, 9, 10, 12, 13, 20, 21, 27, 32, 34, 38, 40, 41, 42, 46, 47, 62, 64, 65, 66, 70, 71)
##PCA2_signif_columns <- c(7, 15, 16, 20, 28, 41, 48, 51, 52, 55, 57, 58, 60, 61, 63, 64, 67, 68, 71)
##
##
##colnames(Summary_RNAseq_table_OLR_relative_species_wide_woinfo[,PCA1_signif_columns])
##colnames(Summary_RNAseq_table_OLR_relative_species_wide_woinfo[,PCA2_signif_columns])



#### pGLS PCA axis vs ecology  --------------------

PCA_axis_df <- 
  Summary_RNAseq_table_OLR_relative_species_wide_info %>%
  dplyr::select(Sp, PC1, PC2)


PCA_axis_df_eco <- 
  left_join(PCA_axis_df, 
            Cichlids_chemoreception_df %>%
              dplyr::select(Sp, Tribe, mean_d15N, mean_d13C, breeding_type, breeding_mode, Food, habitat,
                            median_exploration),
            by="Sp")



caper_data_PCA <- 
  comparative.data(phy = radiation_tree_wo_Neospl, 
                   data = PCA_axis_df_eco,
                   names.col = Sp, vcv = TRUE,
                   na.omit = FALSE, warn.dropped = TRUE)



fit_phylo_PC1_habitat <- pgls(PC1 ~ habitat,
                              data = caper_data_PCA, lambda = "ML")

fit_phylo_PC2_habitat <- pgls(PC2 ~ habitat,
                              data = caper_data_PCA, lambda = "ML")


summary(fit_phylo_PC1_habitat)
summary(fit_phylo_PC2_habitat)


fit_phylo_PC1_Food <- pgls(PC1 ~ Food,
                           data = caper_data_PCA, lambda = "ML")

fit_phylo_PC2_Food <- pgls(PC2 ~ Food,
                           data = caper_data_PCA, lambda = "ML")

summary(fit_phylo_PC1_Food)
summary(fit_phylo_PC2_Food)



fit_phylo_PC1_d13C <- pgls(PC1 ~ mean_d13C,
                           data = caper_data_PCA, lambda = "ML")

fit_phylo_PC2_d13C <- pgls(PC2 ~ mean_d13C,
                           data = caper_data_PCA, lambda = "ML")

summary(fit_phylo_PC1_d13C)
summary(fit_phylo_PC2_d13C)


fit_phylo_PC1_d15N <- pgls(PC1 ~ mean_d15N,
                           data = caper_data_PCA, lambda = "ML")

fit_phylo_PC2_d15N <- pgls(PC2 ~ mean_d15N,
                           data = caper_data_PCA, lambda = "ML")


summary(fit_phylo_PC1_d15N)
summary(fit_phylo_PC2_d15N)



sum_fit_phy <- summary(fit_phylo_PC2_d15N)
PGLS_pente =   formatC(sum_fit_phy$coefficients[2], digits = 3)
PGLS_r2 =      formatC( sum_fit_phy$r.squared, digits = 2)
PGLS_pvalue =  formatC(sum_fit_phy$coefficients[8], digits = 3)
if (PGLS_pvalue == "   0"){ PGLS_pvalue = "< 2.2e-16"}
PGLS_cc <- coef(fit_phylo_PC2_d15N)
PGLS_fit_function <- function(x) PGLS_cc[1] + PGLS_cc[2]*x
PGLS_lambda <- sum_fit_phy$param[2]
PGLS_r2


lm_fit <- lm(data = PCA_axis_df_eco,
             formula = PC2 ~ mean_d15N)
sum_fit <- summary(lm_fit)
GLS_r2 =      formatC(sum_fit$r.squared, digits = 3)
GLS_pente =   formatC(sum_fit$coefficients[2], digits = 3)
GLS_pvalue =  formatC(sum_fit$coefficients[8], digits = 3)
GLS_cc <- coef(lm_fit)
GLS_fit_function <- function(x) GLS_cc[1] + GLS_cc[2]*x

pdf(file = "Raw_R_plots/expression_PC2_vs_d15N.pdf",width = 8.34,  height = 4.61)

PCA_axis_df_eco %>%
  filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>%
  ggplot(aes(x= mean_d15N, y = PC2, color=Tribe)) +
  geom_point(size = 3) +
  stat_function(fun = GLS_fit_function, color="black")+
  labs(subtitle = bquote(R^2 ~ "=" ~ .(PGLS_r2) ~ ";" ~ P ~ "=" ~.(PGLS_pvalue) ~ ";" ~ lambda ~ "=" ~ .(PGLS_lambda))) +
  theme_classic() +
  xlab("mean d15N") +
  ylab("PC2") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none") + 
  scale_color_manual(values=tribes.colors)


dev.off()


sum_fit_phy <- summary(fit_phylo_PC1_d15N)
PGLS_pente =   formatC(sum_fit_phy$coefficients[2], digits = 3)
PGLS_r2 =      formatC( sum_fit_phy$r.squared, digits = 2)
PGLS_pvalue =  formatC(sum_fit_phy$coefficients[8], digits = 3)
if (PGLS_pvalue == "   0"){ PGLS_pvalue = "< 2.2e-16"}
PGLS_cc <- coef(fit_phylo_PC1_d15N)
PGLS_fit_function <- function(x) PGLS_cc[1] + PGLS_cc[2]*x
PGLS_lambda <- sum_fit_phy$param[2]
PGLS_r2


lm_fit <- lm(data = PCA_axis_df_eco,
             formula = PC1 ~ mean_d15N)
sum_fit <- summary(lm_fit)
GLS_r2 =      formatC( sum_fit$r.squared, digits = 3)
GLS_pente =   formatC(sum_fit$coefficients[2], digits = 3)
GLS_pvalue =  formatC(sum_fit$coefficients[8], digits = 3)
GLS_cc <- coef(lm_fit)
GLS_fit_function <- function(x) GLS_cc[1] + GLS_cc[2]*x

pdf(file = "Raw_R_plots/expression_PC1_vs_d15N.pdf",width = 8.34,  height = 4.61)

PCA_axis_df_eco %>%
  filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>%
  ggplot(aes(x= mean_d15N, y = PC1, color=Tribe)) +
  geom_point() +
  stat_function(fun = GLS_fit_function, color="black")+
  labs(subtitle = bquote(R^2 ~ "=" ~ .(PGLS_r2) ~ ";" ~ P ~ "=" ~.(PGLS_pvalue) ~ ";" ~ lambda ~ "=" ~ .(PGLS_lambda))) +
  theme_classic() +
  xlab("mean d15N") +
  ylab("PC1") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none") + 
  scale_color_manual(values=tribes.colors)


dev.off()


sum_fit_phy <- summary(fit_phylo_PC2_d13C)
PGLS_pente =   formatC(sum_fit_phy$coefficients[2], digits = 3)
PGLS_r2 =      formatC( sum_fit_phy$r.squared, digits = 2)
PGLS_pvalue =  formatC(sum_fit_phy$coefficients[8], digits = 3)
if (PGLS_pvalue == "   0"){ PGLS_pvalue = "< 2.2e-16"}
PGLS_cc <- coef(fit_phylo_PC2_d13C)
PGLS_fit_function <- function(x) PGLS_cc[1] + PGLS_cc[2]*x
PGLS_lambda <- sum_fit_phy$param[2]
PGLS_r2


lm_fit <- lm(data = PCA_axis_df_eco,
             formula = PC2 ~ mean_d13C)
sum_fit <- summary(lm_fit)
GLS_r2 =      formatC( sum_fit$r.squared, digits = 3)
GLS_pente =   formatC(sum_fit$coefficients[2], digits = 3)
GLS_pvalue =  formatC(sum_fit$coefficients[8], digits = 3)
GLS_cc <- coef(lm_fit)
GLS_fit_function <- function(x) GLS_cc[1] + GLS_cc[2]*x

pdf(file = "Raw_R_plots/expression_PC2_vs_d13C.pdf",width = 8.34,  height = 4.61)

PCA_axis_df_eco %>%
  filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>%
  ggplot(aes(x= mean_d13C, y = PC2, color=Tribe)) +
  geom_point() +
  stat_function(fun = GLS_fit_function, color="black")+
  labs(subtitle = bquote(R^2 ~ "=" ~ .(PGLS_r2) ~ ";" ~ P ~ "=" ~.(PGLS_pvalue) ~ ";" ~ lambda ~ "=" ~ .(PGLS_lambda))) +
  theme_classic() +
  xlab("mean d13C") +
  ylab("PC2") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none") + 
  scale_color_manual(values=tribes.colors)


dev.off()


sum_fit_phy <- summary(fit_phylo_PC1_d13C)
PGLS_pente =   formatC(sum_fit_phy$coefficients[2], digits = 3)
PGLS_r2 =      formatC( sum_fit_phy$r.squared, digits = 2)
PGLS_pvalue =  formatC(sum_fit_phy$coefficients[8], digits = 3)
if (PGLS_pvalue == "   0"){ PGLS_pvalue = "< 2.2e-16"}
PGLS_cc <- coef(fit_phylo_PC1_d13C)
PGLS_fit_function <- function(x) PGLS_cc[1] + PGLS_cc[2]*x
PGLS_lambda <- sum_fit_phy$param[2]
PGLS_r2


lm_fit <- lm(data = PCA_axis_df_eco,
             formula = PC1 ~ mean_d13C)
sum_fit <- summary(lm_fit)
GLS_r2 =      formatC( sum_fit$r.squared, digits = 3)
GLS_pente =   formatC(sum_fit$coefficients[2], digits = 3)
GLS_pvalue =  formatC(sum_fit$coefficients[8], digits = 3)
GLS_cc <- coef(lm_fit)
GLS_fit_function <- function(x) GLS_cc[1] + GLS_cc[2]*x

pdf(file = "Raw_R_plots/expression_PC1_vs_d13C.pdf",width = 8.34,  height = 4.61)

PCA_axis_df_eco %>%
  filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>%
  ggplot(aes(x= mean_d13C, y = PC1, color=Tribe)) +
  geom_point() +
  stat_function(fun = GLS_fit_function, color="black")+
  labs(subtitle = bquote(R^2 ~ "=" ~ .(PGLS_r2) ~ ";" ~ P ~ "=" ~.(PGLS_pvalue) ~ ";" ~ lambda ~ "=" ~ .(PGLS_lambda))) +
  theme_classic() +
  xlab("mean d13C") +
  ylab("PC1") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none") + 
  scale_color_manual(values=tribes.colors)


dev.off()







fit_phylo_PC1_breed <- pgls(PC1 ~ breeding_mode,
                            data = caper_data_PCA, lambda = "ML")

fit_phylo_PC2_breed <- pgls(PC2 ~ breeding_mode,
                            data = caper_data_PCA, lambda = "ML")


summary(fit_phylo_PC1_breed)
summary(fit_phylo_PC2_breed)


pdf(file = "Raw_R_plots/Ecology_PC1_breedingMode.pdf",width = 8.34,  height = 4.61)


PCA_axis_df_eco %>%
  filter(Sp != "Neospl") %>%
  filter(! is.na(breeding_mode)) %>%
  ggplot(aes(x= breeding_mode, y = PC1, color=breeding_mode)) +
  geom_boxplot() +
  theme_classic() +
  ylab("PC1") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none") 


dev.off()




#### RE vs d13C vs d15N -- AIC strategy   --------------------

#Create a d15N_d13C_df


d15N_d13C_df <- 
  Cichlids_chemoreception_df %>%
  dplyr::select(Sp, mean_d15N, mean_d13C, Tribe) %>%
  distinct()


#Now combine table of the different OLR families
ALL_RE_PCG_df <- 
  left_join(OR_RE_PCG_df,
            TAAR_RE_PCG_df,
            by="Sp")
ALL_RE_PCG_df <- 
  left_join(ALL_RE_PCG_df,
            V2R_RE_PCG_df,
            by="Sp")



Summary_RNAseq_table_V1R_relative_species_subfam_wide <- 
  as.data.frame(
    pivot_wider(Summary_RNAseq_table_V1R_relative_species_subfam,
                names_from = c(subfamily), 
                values_from = relative_expression))


ALL_RE_PCG_df <- 
  left_join(ALL_RE_PCG_df,
            Summary_RNAseq_table_V1R_relative_species_subfam_wide,
            by="Sp")

ALL_RE_PCG_df <- 
  ALL_RE_PCG_df %>%
  mutate(ORA_nb = 1)



#Now lets make a list of subfamilies
list_all_subfamilies <- 
  colnames(ALL_RE_PCG_df)[colnames(ALL_RE_PCG_df) %>% grepl(pattern = ".RE")]
list_all_subfamilies <- c(list_all_subfamilies, "ORA1", "ORA2", "ORA3", "ORA4", "ORA5", "ORA6")

list_all_subfamilies <- gsub(".RE", "", list_all_subfamilies)

#Now lets make correlations between RE and  d15N and d13C, while taking into account or not the PCG


abnorm_termination_1 <- c("Delta_19")
abnorm_termination_2 <-  c("Delta_1")
abnorm_termination_3 <- c("V2RD11_2")


Summary_best_models <- as.data.frame(NULL)
AIC_models_OR_TAAR_V2R <- as.data.frame(NULL)
AIC_models_ORA <- as.data.frame(NULL)
for(subfam in list_all_subfamilies){
  
  
  #Generate a table with the current subfamily
  if(subfam %in% c("ORA1", "ORA2", "ORA3", "ORA4", "ORA5", "ORA6")){
    RE_col <- subfam
    PCG_col <- "ORA_nb"
  } else {
    RE_col <- paste(subfam, "RE", sep=".")
    PCG_col <- paste(subfam, "PCG", sep=".")
  }
  
  
  curr_ALL_RE_PCG_df <- 
    ALL_RE_PCG_df %>%
    dplyr::select(Sp, RE_col, PCG_col)
  
  curr_ALL_RE_PCG_df <- 
    left_join(curr_ALL_RE_PCG_df,
              d15N_d13C_df,
              by="Sp")
  
  colnames(curr_ALL_RE_PCG_df) <- c("Sp", "RE", "PCG", "mean_d15N", "mean_d13C", "Tribe")
  
  
  #Generate a caper object 
  curr_caper_RE_PCG <-
    comparative.data(phy = radiation_tree_wo_Neospl, 
                     data = curr_ALL_RE_PCG_df,
                     names.col = Sp, vcv = TRUE,
                     na.omit = FALSE, warn.dropped = TRUE)
  
  
  
  
  #pGLS between the relative expression and the prop, or d13C * prop
  if(! subfam %in% c("ORA1", "ORA2", "ORA3", "ORA4", "ORA5", "ORA6")){
    
    #Test all possible pGLS models including PCG, d15N and d13C
    fit_phylo_subfam_prop <-
      pgls(RE ~ PCG, 
           data = curr_caper_RE_PCG, 
           lambda = "ML", 
           bounds=list(lambda=c(0,1)))
  
    
    if(subfam %in% abnorm_termination_2){ 
      fit_phylo_subfam_d13C_abs <-
        pgls(RE ~ mean_d13C, 
             data = curr_caper_RE_PCG, 
             lambda = "ML", 
             bounds=list(lambda=c(0.1,1)))
    } else {
      fit_phylo_subfam_d13C_abs <-
        pgls(RE ~ mean_d13C, 
             data = curr_caper_RE_PCG, 
             lambda = "ML", 
             bounds=list(lambda=c(0,1)))
    }
    
    
    
    fit_phylo_subfam_d15N_abs <-
      pgls(RE ~ mean_d15N, 
           data = curr_caper_RE_PCG, 
           lambda = "ML", 
           bounds=list(lambda=c(0,1)))
    
    fit_phylo_subfam_d15N_d13C_abs_add <-
      pgls(RE ~ mean_d15N + mean_d13C, 
           data = curr_caper_RE_PCG, 
           lambda = "ML", 
           bounds=list(lambda=c(0,1)))
    
    fit_phylo_subfam_d15N_d13C_abs_int <-
      pgls(RE ~ mean_d15N * mean_d13C, 
           data = curr_caper_RE_PCG, 
           lambda = "ML", 
           bounds=list(lambda=c(0,1)))
    
    
    fit_phylo_subfam_d13C_prop_add <-
      pgls(RE ~ mean_d13C + PCG, 
           data = curr_caper_RE_PCG, 
           lambda = "ML", 
           bounds=list(lambda=c(0,1)))
    
    
    fit_phylo_subfam_d13C_prop_int <-
      pgls(RE ~ mean_d13C * PCG, 
           data = curr_caper_RE_PCG, 
           lambda = "ML", 
           bounds=list(lambda=c(0,1)))
    
    
    if(subfam %in% abnorm_termination_1){ 
      
      fit_phylo_subfam_d15N_prop_add <-
        pgls(RE ~ mean_d15N + PCG, 
             data = curr_caper_RE_PCG, 
             lambda = "ML", 
             bounds=list(lambda=c(0.2,1)))
    } else {
      fit_phylo_subfam_d15N_prop_add <-
        pgls(RE ~ mean_d15N + PCG, 
             data = curr_caper_RE_PCG, 
             lambda = "ML", 
             bounds=list(lambda=c(0,1)))
    }
    

    
    
    fit_phylo_subfam_d15N_prop_int <-
      pgls(RE ~ mean_d15N * PCG, 
           data = curr_caper_RE_PCG, 
           lambda = "ML", 
           bounds=list(lambda=c(0,1)))
    
    
    if(subfam %in% abnorm_termination_3){ 
      fit_phylo_subfam_d15N_d13C_prop_add <-
        pgls(RE ~ mean_d15N + mean_d13C + PCG, 
             data = curr_caper_RE_PCG, 
             lambda = "ML", 
             bounds=list(lambda=c(0.1,1)))
    } else {
      fit_phylo_subfam_d15N_d13C_prop_add <-
        pgls(RE ~ mean_d15N + mean_d13C + PCG, 
             data = curr_caper_RE_PCG, 
             lambda = "ML", 
             bounds=list(lambda=c(0,1)))
    }
    
    

    
    
    fit_phylo_subfam_d15N_d13C_prop_int <-
      pgls(RE ~ mean_d15N * mean_d13C + PCG, 
           data = curr_caper_RE_PCG, 
           lambda = "ML", 
           bounds=list(lambda=c(0,1)))
    
    
    fit_phylo_subfam_d13C_d15N_prop_fullint <-
      pgls(RE ~ mean_d15N * mean_d13C * PCG, 
           data = curr_caper_RE_PCG, 
           lambda = "ML", 
           bounds=list(lambda=c(0,1)))
    
    
    
    #Extract the AIC of the best complex model
    best_AIC_complex <- 
        head(
          as.data.frame(
            AIC(fit_phylo_subfam_prop,
                fit_phylo_subfam_d13C_abs,
                fit_phylo_subfam_d15N_abs,
                fit_phylo_subfam_d15N_d13C_abs_add,
                fit_phylo_subfam_d15N_d13C_abs_int,
                fit_phylo_subfam_d13C_prop_add,
                fit_phylo_subfam_d13C_prop_int,
                fit_phylo_subfam_d15N_prop_add,
                fit_phylo_subfam_d15N_prop_int,
                fit_phylo_subfam_d15N_d13C_prop_add,
                fit_phylo_subfam_d15N_d13C_prop_int,
                fit_phylo_subfam_d13C_d15N_prop_fullint
            )) %>%
            arrange(AIC), 1)
    best_AIC_value <- best_AIC_complex$AIC
    best_model <- rownames(best_AIC_complex)
  
    #Extract the AIC of the null model (only PCG)
    AIC_null_model <- AIC(fit_phylo_subfam_prop)
    
    #Lets make the difference between AIC values
    delta_AIC <- AIC_null_model - best_AIC_value
    
    
    #Make a nice AIC table
    curr_AIC <- as.data.frame(
      as.data.frame(AIC(fit_phylo_subfam_prop,
                        fit_phylo_subfam_d13C_abs,
                        fit_phylo_subfam_d15N_abs,
                        fit_phylo_subfam_d15N_d13C_abs_add,
                        fit_phylo_subfam_d15N_d13C_abs_int,
                        fit_phylo_subfam_d13C_prop_add,
                        fit_phylo_subfam_d13C_prop_int,
                        fit_phylo_subfam_d15N_prop_add,
                        fit_phylo_subfam_d15N_prop_int,
                        fit_phylo_subfam_d15N_d13C_prop_add,
                        fit_phylo_subfam_d15N_d13C_prop_int,
                        fit_phylo_subfam_d13C_d15N_prop_fullint)) %>%
        dplyr::select(-df) %>%
        tibble::rownames_to_column(var = "model") %>%
        pivot_wider(names_from = model, values_from = AIC)) %>% 
      mutate(subfamily = subfam)
    AIC_models_OR_TAAR_V2R <- rbind(AIC_models_OR_TAAR_V2R, curr_AIC)
    
    sum_fit_PCG <- summary(fit_phylo_subfam_prop)
    simple_PCG_r2 = formatC(sum_fit_PCG$adj.r.squared, digits = 2)
    
    

    if(best_model == "fit_phylo_subfam_d15N_d13C_abs_add" & (delta_AIC >= 2)){
      
      
      sum_fit_phy <- summary(fit_phylo_subfam_d15N_d13C_abs_add)
      PGLS_pente =   formatC(sum_fit_phy$coefficients[2], digits = 3)
      PGLS_r2 =      formatC( sum_fit_phy$adj.r.squared, digits = 2)
      PGLS_lambda <- sum_fit_phy$param[2]
      
      setwd("PGLS_parsing_folder/")
      sink(file = "Current_PGLS.txt")
      print(sum_fit_phy)
      sink(NULL)
      system("grep '^F-statistic' Current_PGLS.txt  | sed 's/.* p-value: //g' > curr_pvalue")
      PGLS_pvalue =  
        as.numeric(scan("PGLS_parsing_folder/curr_pvalue", what="character"))
      
      
    } else if (best_model == "fit_phylo_subfam_d15N_d13C_abs_int" & (delta_AIC >= 2)){
      
      sum_fit_phy <- summary(fit_phylo_subfam_d15N_d13C_abs_int)
      PGLS_pente =   formatC(sum_fit_phy$coefficients[2], digits = 3)
      PGLS_r2 =      formatC( sum_fit_phy$adj.r.squared, digits = 2)
      PGLS_lambda <- sum_fit_phy$param[2]
      
      setwd("PGLS_parsing_folder/")
      sink(file = "Current_PGLS.txt")
      print(sum_fit_phy)
      sink(NULL)
      system("grep '^F-statistic' Current_PGLS.txt  | sed 's/.* p-value: //g' > curr_pvalue")
      PGLS_pvalue =  
        as.numeric(scan("PGLS_parsing_folder/curr_pvalue", what="character"))
      
      
    } else if (best_model == "fit_phylo_subfam_d13C_prop_add" & (delta_AIC >= 2)){
      
      sum_fit_phy <- summary(fit_phylo_subfam_d13C_prop_add)
      PGLS_pente =   formatC(sum_fit_phy$coefficients[2], digits = 3)
      PGLS_r2 =      formatC( sum_fit_phy$adj.r.squared, digits = 2)
      PGLS_lambda <- sum_fit_phy$param[2]
      
      setwd("PGLS_parsing_folder/")
      sink(file = "Current_PGLS.txt")
      print(sum_fit_phy)
      sink(NULL)
      system("grep '^F-statistic' Current_PGLS.txt  | sed 's/.* p-value: //g' > curr_pvalue")
      PGLS_pvalue =  
        as.numeric(scan("PGLS_parsing_folder/curr_pvalue", what="character"))
      
      
    } else if (best_model == "fit_phylo_subfam_d13C_prop_int" & (delta_AIC >= 2)){
      
      sum_fit_phy <- summary(fit_phylo_subfam_d13C_prop_int)
      PGLS_pente =   formatC(sum_fit_phy$coefficients[2], digits = 3)
      PGLS_r2 =      formatC( sum_fit_phy$adj.r.squared, digits = 2)
      PGLS_lambda <- sum_fit_phy$param[2]
      
      setwd("PGLS_parsing_folder/")
      sink(file = "Current_PGLS.txt")
      print(sum_fit_phy)
      sink(NULL)
      system("grep '^F-statistic' Current_PGLS.txt  | sed 's/.* p-value: //g' > curr_pvalue")
      PGLS_pvalue =  
        as.numeric(scan("PGLS_parsing_folder/curr_pvalue", what="character"))
      
      
    } else if (best_model == "fit_phylo_subfam_d15N_prop_add" & (delta_AIC >= 2)){
      
      sum_fit_phy <- summary(fit_phylo_subfam_d15N_prop_add)
      PGLS_pente =   formatC(sum_fit_phy$coefficients[2], digits = 3)
      PGLS_r2 =      formatC( sum_fit_phy$adj.r.squared, digits = 2)
      PGLS_lambda <- sum_fit_phy$param[2]
      
      setwd("PGLS_parsing_folder/")
      sink(file = "Current_PGLS.txt")
      print(sum_fit_phy)
      sink(NULL)
      system("grep '^F-statistic' Current_PGLS.txt  | sed 's/.* p-value: //g' > curr_pvalue")
      PGLS_pvalue =  
        as.numeric(scan("PGLS_parsing_folder/curr_pvalue", what="character"))
      
      
    } else if (best_model == "fit_phylo_subfam_d15N_prop_int" & (delta_AIC >= 2)){
      
      sum_fit_phy <- summary(fit_phylo_subfam_d15N_prop_int)
      PGLS_pente =   formatC(sum_fit_phy$coefficients[2], digits = 3)
      PGLS_r2 =      formatC( sum_fit_phy$adj.r.squared, digits = 2)
      PGLS_lambda <- sum_fit_phy$param[2]
      
      setwd("PGLS_parsing_folder/")
      sink(file = "Current_PGLS.txt")
      print(sum_fit_phy)
      sink(NULL)
      system("grep '^F-statistic' Current_PGLS.txt  | sed 's/.* p-value: //g' > curr_pvalue")
      PGLS_pvalue =  
        as.numeric(scan("PGLS_parsing_folder/curr_pvalue", what="character"))
      
      
    } else if (best_model == "fit_phylo_subfam_d15N_d13C_prop_add" & (delta_AIC >= 2)){
      
      sum_fit_phy <- summary(fit_phylo_subfam_d15N_d13C_prop_add)
      PGLS_pente =   formatC(sum_fit_phy$coefficients[2], digits = 3)
      PGLS_r2 =      formatC( sum_fit_phy$adj.r.squared, digits = 2)
      PGLS_lambda <- sum_fit_phy$param[2]
      
      setwd("PGLS_parsing_folder/")
      sink(file = "Current_PGLS.txt")
      print(sum_fit_phy)
      sink(NULL)
      system("grep '^F-statistic' Current_PGLS.txt  | sed 's/.* p-value: //g' > curr_pvalue")
      PGLS_pvalue =  
        as.numeric(scan("PGLS_parsing_folder/curr_pvalue", what="character"))
      
      
    } else if (best_model == "fit_phylo_subfam_d15N_d13C_prop_int" & (delta_AIC >= 2)){
      
      sum_fit_phy <- summary(fit_phylo_subfam_d15N_d13C_prop_int)
      PGLS_pente =   formatC(sum_fit_phy$coefficients[2], digits = 3)
      PGLS_r2 =      formatC( sum_fit_phy$adj.r.squared, digits = 2)
      PGLS_lambda <- sum_fit_phy$param[2]
      
      setwd("PGLS_parsing_folder/")
      sink(file = "Current_PGLS.txt")
      print(sum_fit_phy)
      sink(NULL)
      system("grep '^F-statistic' Current_PGLS.txt  | sed 's/.* p-value: //g' > curr_pvalue")
      PGLS_pvalue =  
        as.numeric(scan("PGLS_parsing_folder/curr_pvalue", what="character"))
      
      
    } else if (best_model == "fit_phylo_subfam_d13C_d15N_prop_fullint" & (delta_AIC >= 2)){
      
      sum_fit_phy <- summary(fit_phylo_subfam_d13C_d15N_prop_fullint)
      PGLS_pente =   formatC(sum_fit_phy$coefficients[2], digits = 3)
      PGLS_r2 =      formatC( sum_fit_phy$adj.r.squared, digits = 2)
      PGLS_lambda <- sum_fit_phy$param[2]
      
      setwd("PGLS_parsing_folder/")
      sink(file = "Current_PGLS.txt")
      print(sum_fit_phy)
      sink(NULL)
      system("grep '^F-statistic' Current_PGLS.txt  | sed 's/.* p-value: //g' > curr_pvalue")
      PGLS_pvalue =  
        as.numeric(scan("PGLS_parsing_folder/curr_pvalue", what="character"))
      
      
    } else if (best_model == "fit_phylo_subfam_d15N_abs"){
      
      sum_fit_phy <- summary(fit_phylo_subfam_d15N_abs)
      PGLS_pente =   formatC(sum_fit_phy$coefficients[2], digits = 3)
      PGLS_r2 =      formatC( sum_fit_phy$adj.r.squared, digits = 2)
      PGLS_lambda <- sum_fit_phy$param[2]
      PGLS_pvalue =  formatC(sum_fit_phy$coefficients[8], digits = 3)
      
    } else if (best_model == "fit_phylo_subfam_d13C_abs"){
      
      sum_fit_phy <- summary(fit_phylo_subfam_d13C_abs)
      PGLS_pente =   formatC(sum_fit_phy$coefficients[2], digits = 3)
      PGLS_r2 =      formatC( sum_fit_phy$adj.r.squared, digits = 2)
      PGLS_lambda <- sum_fit_phy$param[2]
      PGLS_pvalue =  formatC(sum_fit_phy$coefficients[8], digits = 3)
      
    } else if (best_model == "fit_phylo_subfam_prop"){
      
      sum_fit_phy <- summary(fit_phylo_subfam_prop)
      PGLS_pente =   formatC(sum_fit_phy$coefficients[2], digits = 3)
      PGLS_r2 =      formatC( sum_fit_phy$adj.r.squared, digits = 2)
      PGLS_lambda <- sum_fit_phy$param[2]
      PGLS_pvalue =  formatC(sum_fit_phy$coefficients[8], digits = 3)
      
    } else {
      
      sum_fit_phy <- summary(fit_phylo_subfam_prop)
      PGLS_pente =   formatC(sum_fit_phy$coefficients[2], digits = 3)
      PGLS_r2 =      formatC( sum_fit_phy$adj.r.squared, digits = 2)
      PGLS_lambda <- sum_fit_phy$param[2]
      PGLS_pvalue =  formatC(sum_fit_phy$coefficients[8], digits = 3)
      
    }
    
    
    #Now lets do ORA genes
    
  }  else {
    
    
    fit_phylo_subfam_d13C_abs <-
      pgls(RE ~ mean_d13C, 
           data = curr_caper_RE_PCG, 
           lambda = "ML", 
           bounds=list(lambda=c(0,1)))
    
    fit_phylo_subfam_d15N_abs <-
      pgls(RE ~ mean_d15N, 
           data = curr_caper_RE_PCG, 
           lambda = "ML", 
           bounds=list(lambda=c(0,1)))
    
    fit_phylo_subfam_d15N_d13C_abs_add <-
      pgls(RE ~ mean_d15N + mean_d13C, 
           data = curr_caper_RE_PCG, 
           lambda = "ML", 
           bounds=list(lambda=c(0,1)))
    
    fit_phylo_subfam_d15N_d13C_abs_int <-
      pgls(RE ~ mean_d15N * mean_d13C, 
           data = curr_caper_RE_PCG, 
           lambda = "ML", 
           bounds=list(lambda=c(0,1)))
    
    
    
    
    best_model <- 
      rownames(
        head(
          as.data.frame(
            AIC(
              fit_phylo_subfam_d13C_abs,
              fit_phylo_subfam_d15N_abs,
              fit_phylo_subfam_d15N_d13C_abs_add,
              fit_phylo_subfam_d15N_d13C_abs_int
            )) %>%
            arrange(AIC), 1)
      )
    
    
  
    curr_AIC <- as.data.frame(
      as.data.frame(AIC(
        fit_phylo_subfam_d13C_abs,
        fit_phylo_subfam_d15N_abs,
        fit_phylo_subfam_d15N_d13C_abs_add,
        fit_phylo_subfam_d15N_d13C_abs_int)) %>%
        dplyr::select(-df) %>%
        tibble::rownames_to_column(var = "model") %>%
        pivot_wider(names_from = model, values_from = AIC)) %>% 
      mutate(subfamily = subfam)
    AIC_models_ORA <- rbind(AIC_models_ORA, curr_AIC)
    
    
    if(best_model == "fit_phylo_subfam_d15N_d13C_abs_add"){
      
      
      sum_fit_phy <- summary(fit_phylo_subfam_d15N_d13C_abs_add)
      PGLS_pente =   formatC(sum_fit_phy$coefficients[2], digits = 3)
      PGLS_r2 =      formatC( sum_fit_phy$adj.r.squared, digits = 2)
      PGLS_lambda <- sum_fit_phy$param[2]
      
      setwd("PGLS_parsing_folder/")
      sink(file = "Current_PGLS.txt")
      print(sum_fit_phy)
      sink(NULL)
      system("grep '^F-statistic' Current_PGLS.txt  | sed 's/.* p-value: //g' > curr_pvalue")
      PGLS_pvalue =  
        as.numeric(scan("PGLS_parsing_folder/curr_pvalue", what="character"))
      
      
    } else if (best_model == "fit_phylo_subfam_d15N_d13C_abs_int"){
      
      sum_fit_phy <- summary(fit_phylo_subfam_d15N_d13C_abs_int)
      PGLS_pente =   formatC(sum_fit_phy$coefficients[2], digits = 3)
      PGLS_r2 =      formatC( sum_fit_phy$adj.r.squared, digits = 2)
      PGLS_lambda <- sum_fit_phy$param[2]
      
      setwd("PGLS_parsing_folder/")
      sink(file = "Current_PGLS.txt")
      print(sum_fit_phy)
      sink(NULL)
      system("grep '^F-statistic' Current_PGLS.txt  | sed 's/.* p-value: //g' > curr_pvalue")
      PGLS_pvalue =  
        as.numeric(scan("PGLS_parsing_folder/curr_pvalue", what="character"))
      
      
    } else if (best_model == "fit_phylo_subfam_d15N_abs"){
      
      sum_fit_phy <- summary(fit_phylo_subfam_d15N_abs)
      PGLS_pente =   formatC(sum_fit_phy$coefficients[2], digits = 3)
      PGLS_r2 =      formatC( sum_fit_phy$adj.r.squared, digits = 2)
      PGLS_lambda <- sum_fit_phy$param[2]
      PGLS_pvalue =  formatC(sum_fit_phy$coefficients[8], digits = 3)
      
    } else if (best_model == "fit_phylo_subfam_d13C_abs"){
      
      sum_fit_phy <- summary(fit_phylo_subfam_d13C_abs)
      PGLS_pente =   formatC(sum_fit_phy$coefficients[2], digits = 3)
      PGLS_r2 =      formatC( sum_fit_phy$adj.r.squared, digits = 2)
      PGLS_lambda <- sum_fit_phy$param[2]
      PGLS_pvalue =  formatC(sum_fit_phy$coefficients[8], digits = 3)
      
    }
    
  }
  
  curr_df <- as.data.frame(cbind(PGLS_r2, PGLS_pvalue, PGLS_pente, PGLS_lambda, simple_PCG_r2))
  colnames(curr_df) <- c("R2", "pvalue", "pente", "lambda", "PCG_R2")
  curr_df <- curr_df %>% mutate(subfamily = subfam) %>% mutate(model = best_model)
  
  
  Summary_best_models <- rbind(Summary_best_models, curr_df)
  
  
  
}



#Parse the dataframe 

Summary_best_models$R2 <- as.numeric(Summary_best_models$R2)
Summary_best_models$pvalue <- as.numeric(Summary_best_models$pvalue)
Summary_best_models$pente <- as.numeric(Summary_best_models$pente)
Summary_best_models$lambda <- as.numeric(Summary_best_models$lambda)
Summary_best_models$PCG_R2 <- as.numeric(Summary_best_models$PCG_R2)


Summary_best_models <- 
  Summary_best_models %>% 
  mutate(predictor_pGLS =  case_when(
    model == "fit_phylo_subfam_prop" ~ "PCG",
    model == "fit_phylo_subfam_d13C_abs" ~ "d13C",
    model == "fit_phylo_subfam_d15N_abs" ~ "d15N",
    model == "fit_phylo_subfam_d15N_d13C_abs_add" ~ "d15N + d13C",
    model == "fit_phylo_subfam_d15N_d13C_abs_int" ~ "d15N * d13C",
    model == 'fit_phylo_subfam_d13C_prop_add' ~ "d13C + PCG",
    model == 'fit_phylo_subfam_d13C_prop_int' ~ "d13C * PCG",
    model == 'fit_phylo_subfam_d15N_prop_add' ~ "d15N + PCG",
    model == 'fit_phylo_subfam_d15N_prop_int' ~ "d15N * PCG",
    model == 'fit_phylo_subfam_d15N_d13C_prop_add' ~ "d15N + d13C + PCG",
    model == 'fit_phylo_subfam_d15N_d13C_prop_int' ~ "d15N * d13C + PCG",
    model == 'fit_phylo_subfam_d13C_d15N_prop_fullint' ~ "d15N * d13C * PCG"
  ))


colnames(AIC_models_OR_TAAR_V2R) <- c("PCG", "d13C", "d15N", "d15N + d13C", "d15N * d13C", "d13C + PCG", "d13C * PCG",
                          "d15N + PCG", "d15N * PCG", "d15N + d13C + PCG", "d15N * d13C + PCG", 
                          "d15N * d13C * PCG", "subfamily")

colnames(AIC_models_ORA) <- c("d13C", "d15N", "d15N + d13C", "d15N * d13C", "subfamily")







#Determine if the best pGLS model has a significant p-value or not

Summary_best_models <- 
  Summary_best_models %>%
  mutate(significant = case_when(
    pvalue <= 0.05 ~ model,
    pvalue > 0.05 ~ "Not_signif",
  ))



Summary_best_models %>%
  group_by(significant) %>%
  summarise(n())


#Determine the unique predictors corresponding to each model 
Summary_best_models <-
  Summary_best_models %>%
  mutate(signif_predictors = case_when(
    significant == "Not_signif" ~ "No significant pGLS model",
    significant == "fit_phylo_subfam_prop" ~ "PCG",
    significant == "fit_phylo_subfam_d15N_prop_int" ~ "PCG ; d15N",
    significant == "fit_phylo_subfam_d15N_prop_add" ~ "PCG ; d15N",
    significant == "fit_phylo_subfam_d15N_d13C_prop_int" ~ "PCG ; d15N ; d13C",
    significant == "fit_phylo_subfam_d15N_d13C_abs_int" ~ "d15N ; d13C",
    significant == "fit_phylo_subfam_d13C_prop_int" ~ "PCG ; d13C",
    significant == "fit_phylo_subfam_d13C_prop_add" ~ "PCG ; d13C",
    significant == "fit_phylo_subfam_d13C_d15N_prop_fullint" ~ "PCG ; d15N ; d13C",
    significant == "fit_phylo_subfam_d13C_abs" ~ "d13C",
  ))


Summary_best_models %>%
  group_by(signif_predictors) %>%
  summarise(n())




#Define colors of each model
model_colors <- 
  c("PCG"="#009E73",
    "PCG ; d15N" = "#56B4E9", 
    "PCG ; d15N ; d13C" = "#0072B2",
    "d15N ; d13C" = "#E69F00",
    "PCG ; d13C" = "#F0E442", 
    "d13C" = "#D55E00",
    "Not_signif" = "#999999")




#Transform the dataframe to a long df for plot

uniq_subfam <- 
  Summary_best_models %>% 
  pull(subfamily) %>% unique()

Ecological_niche_vs_RE_df <- as.data.frame(NULL)
for(curr_subfam in uniq_subfam){
  
  
  curr_R2 <- 
    Summary_best_models %>%
    filter(subfamily == curr_subfam) %>%
    pull(R2)
  
  curr_signif_predictors <- 
    Summary_best_models %>%
    filter(subfamily == curr_subfam) %>%
    pull(signif_predictors)
  
  if(curr_signif_predictors == "PCG"){
    curr_PCG <- as.data.frame(cbind(curr_subfam, "PCG", "Significant", curr_R2))
    curr_d15N <- as.data.frame(cbind(curr_subfam, "d15N", "non-significant", curr_R2))
    curr_d13C <- as.data.frame(cbind(curr_subfam, "d13C", "non-significant", curr_R2))
    colnames(curr_PCG) <- c("Response", "predictor", "signif", "R2")
    colnames(curr_d15N) <- c("Response", "predictor", "signif", "R2")
    colnames(curr_d13C) <- c("Response", "predictor", "signif", "R2")
    curr_df <- rbind(curr_PCG, curr_d15N, curr_d13C)
  }
  
  if(curr_signif_predictors == "No significant pGLS model"){
    curr_PCG <- as.data.frame(cbind(curr_subfam, "PCG", "non-significant", curr_R2))
    curr_d15N <- as.data.frame(cbind(curr_subfam, "d15N", "non-significant", curr_R2))
    curr_d13C <- as.data.frame(cbind(curr_subfam, "d13C", "non-significant", curr_R2))
    colnames(curr_PCG) <- c("Response", "predictor", "signif", "R2")
    colnames(curr_d15N) <- c("Response", "predictor", "signif", "R2")
    colnames(curr_d13C) <- c("Response", "predictor", "signif", "R2")
    curr_df <- rbind(curr_PCG, curr_d15N, curr_d13C)
  }
  
  
  if(curr_signif_predictors == "PCG ; d15N ; d13C"){
    curr_PCG <- as.data.frame(cbind(curr_subfam, "PCG", "Significant", curr_R2))
    curr_d15N <- as.data.frame(cbind(curr_subfam, "d15N", "Significant", curr_R2))
    curr_d13C <- as.data.frame(cbind(curr_subfam, "d13C", "Significant", curr_R2))
    colnames(curr_PCG) <- c("Response", "predictor", "signif", "R2")
    colnames(curr_d15N) <- c("Response", "predictor", "signif", "R2")
    colnames(curr_d13C) <- c("Response", "predictor", "signif", "R2")
    curr_df <- rbind(curr_PCG, curr_d15N, curr_d13C)
  }
  
  if(curr_signif_predictors == "PCG ; d13C"){
    curr_PCG <- as.data.frame(cbind(curr_subfam, "PCG", "Significant", curr_R2))
    curr_d15N <- as.data.frame(cbind(curr_subfam, "d15N", "non-significant", curr_R2))
    curr_d13C <- as.data.frame(cbind(curr_subfam, "d13C", "Significant", curr_R2))
    colnames(curr_PCG) <- c("Response", "predictor", "signif", "R2")
    colnames(curr_d15N) <- c("Response", "predictor", "signif", "R2")
    colnames(curr_d13C) <- c("Response", "predictor", "signif", "R2")
    curr_df <- rbind(curr_PCG, curr_d15N, curr_d13C)
  }
  
  
  if(curr_signif_predictors == "PCG ; d15N"){
    curr_PCG <- as.data.frame(cbind(curr_subfam, "PCG", "Significant", curr_R2))
    curr_d15N <- as.data.frame(cbind(curr_subfam, "d15N", "Significant", curr_R2))
    curr_d13C <- as.data.frame(cbind(curr_subfam, "d13C", "non-significant", curr_R2))
    colnames(curr_PCG) <- c("Response", "predictor", "signif", "R2")
    colnames(curr_d15N) <- c("Response", "predictor", "signif", "R2")
    colnames(curr_d13C) <- c("Response", "predictor", "signif", "R2")
    curr_df <- rbind(curr_PCG, curr_d15N, curr_d13C)
  }
  
  
  if(curr_signif_predictors == "d13C"){
    curr_PCG <- as.data.frame(cbind(curr_subfam, "PCG", "non-significant", curr_R2))
    curr_d15N <- as.data.frame(cbind(curr_subfam, "d15N", "non-significant", curr_R2))
    curr_d13C <- as.data.frame(cbind(curr_subfam, "d13C", "Significant", curr_R2))
    colnames(curr_PCG) <- c("Response", "predictor", "signif", "R2")
    colnames(curr_d15N) <- c("Response", "predictor", "signif", "R2")
    colnames(curr_d13C) <- c("Response", "predictor", "signif", "R2")
    curr_df <- rbind(curr_PCG, curr_d15N, curr_d13C)
  }
  
  
  if(curr_signif_predictors == "d15N ; d13C"){
    curr_PCG <- as.data.frame(cbind(curr_subfam, "PCG", "non-significant", curr_R2))
    curr_d15N <- as.data.frame(cbind(curr_subfam, "d15N", "Significant", curr_R2))
    curr_d13C <- as.data.frame(cbind(curr_subfam, "d13C", "Significant", curr_R2))
    colnames(curr_PCG) <- c("Response", "predictor", "signif", "R2")
    colnames(curr_d15N) <- c("Response", "predictor", "signif", "R2")
    colnames(curr_d13C) <- c("Response", "predictor", "signif", "R2")
    curr_df <- rbind(curr_PCG, curr_d15N, curr_d13C)
  }
  
  Ecological_niche_vs_RE_df <- rbind(Ecological_niche_vs_RE_df, curr_df)
}

Ecological_niche_vs_RE_df$R2 <- as.numeric(Ecological_niche_vs_RE_df$R2)

Signif_colors <-
  c("non-significant" = 0,
    "Significant" = 1)


Ecological_niche_vs_RE_df$predictor <-
  factor(Ecological_niche_vs_RE_df$predictor ,
         levels=c("PCG","d15N","d13C"))

Ecological_niche_vs_RE_df$Response <-
  factor(Ecological_niche_vs_RE_df$Response ,
         levels=c("Beta_1","Beta_2","Delta_1","Delta_3","Delta_4","Delta_5",
                  "Delta_6","Delta_7","Delta_8","Delta_9","Delta_10","Delta_11",
                  "Delta_12","Delta_13","Delta_14","Delta_15","Delta_16","Delta_17",
                  "Delta_18","Delta_19","Delta_2","Delta_20","Epsilon_1","Epsilon_2",
                  "Epsilon_3","Epsilon_4","Eta_1","Eta_2","Eta_3","Eta_4",
                  "Eta_6","Eta_7","Eta_8","Eta_9","Eta_10",
                  "Eta_11","Eta_12","Eta_13","Eta_14","Eta_15",
                  "Zeta_1","Zeta_2","Zeta_3","Zeta_4",
                  
                  "TAARA2_1","TAARA2_2","TAARB4_1","TAARB4_2","TAARB4_3","TAARB4_4","TAARB4_5",
                  
                  "V2RD3_1","V2RD4_1","V2RD5_1","V2RD7_1","V2RD8_1","V2RD9_1",
                  "V2RD10_1","V2RD11_1","V2RD11_2","V2RD11_3",
                  "V2RD11_4","V2RD11_5","V2RD11_6","V2RD11_7",
                  
                  "ORA1","ORA2","ORA3","ORA4","ORA5","ORA6")
  )


#Write tables

write.table(Summary_best_models %>% dplyr::select(subfamily, R2, pvalue, pente, lambda, predictor_pGLS),
            "Raw_R_plots/Summary_best_models_RE.tsv",
            sep="\t",
            col.names=TRUE,
            row.names=FALSE,
            quote=FALSE)



write.table(AIC_models_OR_TAAR_V2R,
            "Raw_R_plots/AIC_models_OR_TAAR_V2R.tsv",
            sep="\t",
            col.names=TRUE,
            row.names=FALSE,
            quote=FALSE)

write.table(AIC_models_ORA,
            "Raw_R_plots/AIC_models_ORA.tsv",
            sep="\t",
            col.names=TRUE,
            row.names=FALSE,
            quote=FALSE)


#Lets perform plots

pdf(file = "Raw_R_plots/Expression_model_heatplot.pdf",width = 8.34,  height = 3.61)


ggplot(Ecological_niche_vs_RE_df, 
       aes(y=predictor, x=Response, fill= R2, alpha=signif, color="black")) + 
  geom_tile(color = "black") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_viridis_c(option = "plasma") +
  scale_alpha_manual(values = Signif_colors) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("") +
  ylab("") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


ggplot(Ecological_niche_vs_RE_df, 
       aes(y=predictor, x=Response, fill= R2, alpha=signif, color="black")) + 
  geom_tile(color = "black") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme_minimal() +
  scale_fill_viridis_c(option = "plasma") +
  scale_alpha_manual(values = Signif_colors) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("") +
  ylab("") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


dev.off()





#pdf(file = "Raw_R_plots/Expression_model_treemap.pdf",width = 8.34,  height = 4.61)
#
#Summary_best_models_uncorrected %>%
#  group_by(signif_predictors) %>%
#  summarise(value = n()) %>%
#  ggplot(., aes(area = value, fill = signif_predictors, label = signif_predictors)) +
#  geom_treemap() +
#  geom_treemap_text(colour = "white",
#                    place = "centre",
#                    size = 15) +
#  scale_fill_manual(values = model_colors) +
#  theme(legend.position = "none")
#
#dev.off()



## Lets see the improvment in the R2


Summary_best_models_signif <- 
  Summary_best_models %>%
  filter(significant != "Not_signif") %>%
  filter(significant != "fit_phylo_subfam_prop")


Summary_best_models_signif_subset <-
  Summary_best_models_signif %>%
  dplyr::select(subfamily, R2, PCG_R2)

best_models <-
  Summary_best_models_signif %>%
  dplyr::select(subfamily, signif_predictors)

Summary_best_models_signif_subset_long <-
  Summary_best_models_signif_subset %>%
  pivot_longer(!subfamily ,names_to = "model", values_to = "value")


Summary_best_models_signif_subset_long <-
  left_join(Summary_best_models_signif_subset_long,
            best_models,
            by="subfamily")


Summary_best_models_signif_subset_long <-
  Summary_best_models_signif_subset_long %>%
  mutate(curr_model = case_when(
    model == "PCG_R2" ~ "PCG",
    model == "R2" ~ signif_predictors,
  ))




Summary_best_models_signif_subset_long$subfamily <- 
  factor(Summary_best_models_signif_subset_long$subfamily ,
         levels=rev(c("Eta_1","Eta_2", "Eta_3", "Eta_4", "Eta_6", "Eta_7", "Eta_8", "Eta_9",
                      "Eta_10","Eta_11","Eta_12","Eta_13","Eta_14","Eta_15","Beta_1",
                      "Beta_2","Epsilon_1","Epsilon_2","Epsilon_3","Epsilon_4","Zeta_1","Zeta_2","Zeta_3",
                      "Zeta_4","Delta_1","Delta_2","Delta_3","Delta_4","Delta_5","Delta_6","Delta_7","Delta_8",
                      "Delta_9","Delta_10","Delta_11","Delta_12","Delta_13","Delta_14","Delta_15","Delta_16","Delta_17",
                      "Delta_18","Delta_19","Delta_20","TAARA2_1","TAARA2_2","TAARB4_1","TAARB4_2",
                      "TAARB4_3","TAARB4_4","TAARB4_5","V2RD3_1","V2RD4_1","V2RD5_1",
                      "V2RD7_1","V2RD8_1","V2RD9_1","V2RD10_1","V2RD11_1","V2RD11_2" ,"V2RD11_3",
                      "V2RD11_4","V2RD11_5","V2RD11_6","V2RD11_7","T1R1B","T1R3", "ORA1", "ORA2", 'ORA3',
                      'ORA4', 'ORA5')))



pdf(file = "Raw_R_plots/Expression_model_adjustedRsquared.pdf",width = 8.34,  height = 4.61)

Summary_best_models_signif_subset_long_nonull <- 
  Summary_best_models_signif_subset_long %>%
  mutate(value_nonull = if_else(
    value > 0,
    value,
    0
  ))

Summary_best_models_signif_subset_long_nonull %>%
  filter(subfamily != "ORA5") %>% #remove ORA5 because non significant with PCG only
  ggplot(aes(y = subfamily, x = value_nonull)) +
  geom_line() +
  geom_point(data=Summary_best_models_signif_subset_long_nonull %>% filter(subfamily != "ORA5"), aes(x=value_nonull, y=subfamily, color=curr_model)) +
  scale_color_manual(values = model_colors) + 
  theme_classic() +
  xlab("Adjusted R2") +
  ylab("Subfamily") +
  theme(legend.position='none')

Summary_best_models_signif_subset_long %>%
  filter(subfamily != "ORA5") %>% #remove ORA5 because non significant with PCG only
  ggplot(aes(y = subfamily, x = value)) +
  geom_line() +
  geom_point(data=Summary_best_models_signif_subset_long %>% filter(subfamily != "ORA5"), aes(x=value, y=subfamily, color=curr_model)) +
  scale_color_manual(values = model_colors) + 
  theme_classic() +
  xlab("Adjusted R2") +
  ylab("Subfamily") +
  theme()

dev.off()




#### Expression evolutionnary rate -- OR  ---------------------------------

#Lets define unique pairs of species

species_list <- 
  Summary_RNAseq_table %>% 
  filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>%
  pull(Sp) %>% unique()


all_comb <- apply(combn(species_list,2),2,paste,collapse=',')
pairs_df <- as.data.frame(t(sapply(all_comb, function(x) unlist(strsplit(x, ",")))))
row.names(pairs_df) <- NULL
colnames(pairs_df) <- c("species_1", "species_2")


#Lets perform a spearman correlations for each species pair


spearman_rho_list <- c()
distance_list <- c()
for(pair_nb in 1:nrow(pairs_df)){
  
  sp1 <- pairs_df[pair_nb,]$species_1
  sp2 <- pairs_df[pair_nb,]$species_2
    


  sp1_expression <-
    Summary_RNAseq_table_OR_relative_species_subfam %>%
    filter(Sp == sp1) %>%
    pull(subfam_relative_expr)
  
  sp2_expression <-
    Summary_RNAseq_table_OR_relative_species_subfam %>%
    filter(Sp == sp2) %>%
    pull(subfam_relative_expr)
  
    
  rho_spearman <- 
    as.numeric(cor.test(x=sp1_expression, y=sp2_expression, method="spearman")$estimate)
  
  spearman_rho_list <- 
    c(spearman_rho_list, rho_spearman)
  
  curr_dist <- 
    (get_pairwise_distances(radiation_tree_wo_Neospl, 
                            sp1, 
                              sp2, 
                              as_edge_counts=FALSE, 
                              check_input=TRUE) / 2)
    
  distance_list <- 
    c(distance_list, curr_dist)
  
  
  
}

Rho_vs_dist_OR <-
  as.data.frame(
    cbind(spearman_rho_list, distance_list)
  )

colnames(Rho_vs_dist_OR) <- c("rho", "dist")


lm_rho_dist_OR <- 
  lm(data=Rho_vs_dist_OR,
     formula= rho~dist)
sum_rho_dist_OR <- summary(lm_rho_dist_OR)
pente_rho_dist_OR <- formatC(sum_rho_dist_OR$coefficients[2], digits = 3)
rho_dist_cc_OR <- coef(lm_rho_dist_OR)
rho_dist_F_OR <- function(x) rho_dist_cc_OR[1] + rho_dist_cc_OR[2]*x


Rho_vs_dist_OR %>%
  ggplot(., aes(x=dist, y=rho)) +
  geom_point() +
  theme_classic() 


Rho_vs_dist_OR %>%
  ggplot(., aes(x=dist, y=rho)) +
  geom_line(color="white") +
  theme_classic() +
  stat_function(fun = rho_dist_F_OR, color="black")

#### Expression evolutionnary rate -- TAAR  ---------------------------------

#Lets define unique pairs of species

species_list <- 
  Summary_RNAseq_table %>% 
  filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>%
  pull(Sp) %>% unique()


all_comb <- apply(combn(species_list,2),2,paste,collapse=',')
pairs_df <- as.data.frame(t(sapply(all_comb, function(x) unlist(strsplit(x, ",")))))
row.names(pairs_df) <- NULL
colnames(pairs_df) <- c("species_1", "species_2")


#Lets perform a spearman correlations for each species pair


spearman_rho_list <- c()
distance_list <- c()
for(pair_nb in 1:nrow(pairs_df)){
  
  sp1 <- pairs_df[pair_nb,]$species_1
  sp2 <- pairs_df[pair_nb,]$species_2
  
  
  
  sp1_expression <-
    Summary_RNAseq_table_TAAR_relative_species_subfam %>%
    filter(Sp == sp1) %>%
    pull(subfam_relative_expr)
  
  sp2_expression <-
    Summary_RNAseq_table_TAAR_relative_species_subfam %>%
    filter(Sp == sp2) %>%
    pull(subfam_relative_expr)
  
  
  rho_spearman <- 
    as.numeric(cor.test(x=sp1_expression, y=sp2_expression, method="spearman")$estimate)
  
  spearman_rho_list <- 
    c(spearman_rho_list, rho_spearman)
  
  curr_dist <- 
    (get_pairwise_distances(radiation_tree_wo_Neospl, 
                            sp1, 
                            sp2, 
                            as_edge_counts=FALSE, 
                            check_input=TRUE) / 2)
  
  distance_list <- 
    c(distance_list, curr_dist)
  
  
  
}

Rho_vs_dist_TAAR <-
  as.data.frame(
    cbind(spearman_rho_list, distance_list)
  )

colnames(Rho_vs_dist_TAAR) <- c("rho", "dist")


lm_rho_dist_TAAR <- 
  lm(data=Rho_vs_dist_TAAR,
     formula= rho~dist)
sum_rho_dist_TAAR <- summary(lm_rho_dist_TAAR)
pente_rho_dist_TAAR <- formatC(sum_rho_dist_TAAR$coefficients[2], digits = 3)
rho_dist_cc_TAAR <- coef(lm_rho_dist_TAAR)
rho_dist_F_TAAR <- function(x) rho_dist_cc_TAAR[1] + rho_dist_cc_TAAR[2]*x


Rho_vs_dist_TAAR %>%
  ggplot(., aes(x=dist, y=rho)) +
  geom_point() +
  theme_classic() 


Rho_vs_dist_TAAR %>%
  ggplot(., aes(x=dist, y=rho)) +
  geom_line(color="white") +
  theme_classic() +
  stat_function(fun = rho_dist_F_TAAR, color="black")


#### Expression evolutionnary rate -- V1R  ---------------------------------

#Lets define unique pairs of species

species_list <- 
  Summary_RNAseq_table %>% 
  filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>%
  pull(Sp) %>% unique()


all_comb <- apply(combn(species_list,2),2,paste,collapse=',')
pairs_df <- as.data.frame(t(sapply(all_comb, function(x) unlist(strsplit(x, ",")))))
row.names(pairs_df) <- NULL
colnames(pairs_df) <- c("species_1", "species_2")


#Lets perform a spearman correlations for each species pair


spearman_rho_list <- c()
distance_list <- c()
for(pair_nb in 1:nrow(pairs_df)){
  
  sp1 <- pairs_df[pair_nb,]$species_1
  sp2 <- pairs_df[pair_nb,]$species_2
  
  
  
  sp1_expression <-
    Summary_RNAseq_table_V1R_relative_species_subfam %>%
    filter(Sp == sp1) %>%
    pull(subfam_relative_expr)
  
  sp2_expression <-
    Summary_RNAseq_table_V1R_relative_species_subfam %>%
    filter(Sp == sp2) %>%
    pull(subfam_relative_expr)
  
  
  rho_spearman <- 
    as.numeric(cor.test(x=sp1_expression, y=sp2_expression, method="spearman")$estimate)
  
  spearman_rho_list <- 
    c(spearman_rho_list, rho_spearman)
  
  curr_dist <- 
    (get_pairwise_distances(radiation_tree_wo_Neospl, 
                            sp1, 
                            sp2, 
                            as_edge_counts=FALSE, 
                            check_input=TRUE) / 2)
  
  distance_list <- 
    c(distance_list, curr_dist)
  
  
  
}

Rho_vs_dist_V1R <-
  as.data.frame(
    cbind(spearman_rho_list, distance_list)
  )

colnames(Rho_vs_dist_V1R) <- c("rho", "dist")


lm_rho_dist_V1R <- 
  lm(data=Rho_vs_dist_V1R,
     formula= rho~dist)
sum_rho_dist_V1R <- summary(lm_rho_dist_V1R)
pente_rho_dist_V1R <- formatC(sum_rho_dist_V1R$coefficients[2], digits = 3)
rho_dist_cc_V1R <- coef(lm_rho_dist_V1R)
rho_dist_F_V1R <- function(x) rho_dist_cc_V1R[1] + rho_dist_cc_V1R[2]*x


Rho_vs_dist_V1R %>%
  ggplot(., aes(x=dist, y=rho)) +
  geom_point() +
  theme_classic() 


Rho_vs_dist_V1R %>%
  ggplot(., aes(x=dist, y=rho)) +
  geom_line(color="white") +
  theme_classic() +
  stat_function(fun = rho_dist_F_V1R, color="black")

#### Expression evolutionnary rate -- V2R  ---------------------------------
#Lets define unique pairs of species

species_list <- 
  Summary_RNAseq_table %>% 
  filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>%
  pull(Sp) %>% unique()


all_comb <- apply(combn(species_list,2),2,paste,collapse=',')
pairs_df <- as.data.frame(t(sapply(all_comb, function(x) unlist(strsplit(x, ",")))))
row.names(pairs_df) <- NULL
colnames(pairs_df) <- c("species_1", "species_2")


#Lets perform a spearman correlations for each species pair


spearman_rho_list <- c()
distance_list <- c()
for(pair_nb in 1:nrow(pairs_df)){
  
  sp1 <- pairs_df[pair_nb,]$species_1
  sp2 <- pairs_df[pair_nb,]$species_2
  
  
  
  sp1_expression <-
    Summary_RNAseq_table_V2R_relative_species_subfam %>%
    filter(Sp == sp1) %>%
    pull(subfam_relative_expr)
  
  sp2_expression <-
    Summary_RNAseq_table_V2R_relative_species_subfam %>%
    filter(Sp == sp2) %>%
    pull(subfam_relative_expr)
  
  
  rho_spearman <- 
    as.numeric(cor.test(x=sp1_expression, y=sp2_expression, method="spearman")$estimate)
  
  spearman_rho_list <- 
    c(spearman_rho_list, rho_spearman)
  
  curr_dist <- 
    (get_pairwise_distances(radiation_tree_wo_Neospl, 
                            sp1, 
                            sp2, 
                            as_edge_counts=FALSE, 
                            check_input=TRUE) / 2)
  
  distance_list <- 
    c(distance_list, curr_dist)
  
  
  
}

Rho_vs_dist_V2R <-
  as.data.frame(
    cbind(spearman_rho_list, distance_list)
  )

colnames(Rho_vs_dist_V2R) <- c("rho", "dist")


lm_rho_dist_V2R <- 
  lm(data=Rho_vs_dist_V2R,
     formula= rho~dist)
sum_rho_dist_V2R <- summary(lm_rho_dist_V2R)
pente_rho_dist_V2R <- formatC(sum_rho_dist_V2R$coefficients[2], digits = 3)
rho_dist_cc_V2R <- coef(lm_rho_dist_V2R)
rho_dist_F_V2R <- function(x) rho_dist_cc_V2R[1] + rho_dist_cc_V2R[2]*x


Rho_vs_dist_V2R %>%
  ggplot(., aes(x=dist, y=rho)) +
  geom_point() +
  theme_classic() 


Rho_vs_dist_V2R %>%
  ggplot(., aes(x=dist, y=rho)) +
  geom_line(color="white") +
  theme_classic() +
  stat_function(fun = rho_dist_F_V2R, color="black")


#### Summary -- Expression evolutionnary rate - geomline  ---------------------------------

Rho_vs_dist_OR <- 
  Rho_vs_dist_OR %>%
  mutate(mystat = (1-rho)/(dist)) %>%
  mutate(family = "OR")

Rho_vs_dist_TAAR <- 
  Rho_vs_dist_TAAR %>%
  mutate(mystat = (1-rho)/(dist)) %>%
  mutate(family = "TAAR")

Rho_vs_dist_V1R <- 
  Rho_vs_dist_V1R %>%
  mutate(mystat = (1-rho)/(dist)) %>%
  mutate(family = "V1R")

Rho_vs_dist_V2R <- 
  Rho_vs_dist_V2R %>%
  mutate(mystat = (1-rho)/(dist)) %>%
  mutate(family = "V2R")



Rho_vs_dist_all <- 
  rbind(Rho_vs_dist_OR, Rho_vs_dist_TAAR, Rho_vs_dist_V1R, Rho_vs_dist_V2R)



family_colors <- 
  c(OR="lightcoral",
    TAAR = "#377EB8", 
    V1R = "#FFFF33",
    V2R = "#FF7F00")

Rho_vs_dist_all %>%
  ggplot(., aes(x=dist, y=rho, color=family)) +
  geom_point() +
  scale_color_manual(values = family_colors) +
  stat_function(fun = rho_dist_F_OR, color="lightcoral") +
  stat_function(fun = rho_dist_F_TAAR, color="#377EB8") +
  stat_function(fun = rho_dist_F_V1R, color="#FFFF33") +
  stat_function(fun = rho_dist_F_V2R, color="#FF7F00") 


pdf(file = "Raw_R_plots/Rho_vs_dist_all.pdf",width = 8.34,  height = 4.61)

Rho_vs_dist_OR %>%
  ggplot(., aes(x=dist, y=rho)) +
  geom_point() +
  theme_classic() +
  stat_function(fun = rho_dist_F_OR, color="black") +
  xlab("Time since divergence (Mya)") +
  ylab("Spearman's rho") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none") 


Rho_vs_dist_TAAR %>%
  ggplot(., aes(x=dist, y=rho)) +
  geom_point() +
  theme_classic() +
  stat_function(fun = rho_dist_F_TAAR, color="black") +
  xlab("Time since divergence (Mya)") +
  ylab("Spearman's rho") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none") 



Rho_vs_dist_V1R %>%
  ggplot(., aes(x=dist, y=rho)) +
  geom_point() +
  theme_classic() +
  stat_function(fun = rho_dist_F_V1R, color="black") +
  xlab("Time since divergence (Mya)") +
  ylab("Spearman's rho") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none") 



Rho_vs_dist_V2R %>%
  ggplot(., aes(x=dist, y=rho)) +
  geom_point() +
  theme_classic() +
  stat_function(fun = rho_dist_F_V2R, color="black") +
  xlab("Time since divergence (Mya)") +
  ylab("Spearman's rho") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none") 



dev.off()


pdf(file = "Raw_R_plots/Rho_vs_dist_all_sameplot.pdf",width = 8.34,  height = 4.61)

Rho_vs_dist_V2R %>%
  ggplot(., aes(x=dist, y=rho)) +
  geom_line(color="white") +
  theme_classic() +
  stat_function(fun = rho_dist_F_OR, color="lightcoral") +
  stat_function(fun = rho_dist_F_TAAR, color="#377EB8") +
  stat_function(fun = rho_dist_F_V1R, color="#FFFF33") +
  stat_function(fun = rho_dist_F_V2R, color="#FF7F00") +
  xlim(3, 12.3) +
  xlab("Time since divergence") +
  ylab("Spearman's rho") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none") 

dev.off()

#### Summary -- Expression evolutionnary rate - boxplot  ---------------------------------



pdf(file = "Raw_R_plots/one_minus_rho_over_dist.pdf",width = 8.34,  height = 4.61)

Rho_vs_dist_all %>%
  ggplot(., aes(x=family, y=mystat, color=family)) +
  geom_boxplot() +
  theme_classic() +
  ylim(0, 0.10) +
  xlab("") +
  ylab("1- rho / time since divergence") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none") +
  scale_color_manual(values = olfactory_family_colors)

dev.off()

#model <- aov(mystat~family, data=Rho_vs_dist_all)
#summary(model)
#TukeyHSD(model, conf.level=.95)


wilcox.test(Rho_vs_dist_all %>% filter(family == "OR") %>% pull(mystat),
            Rho_vs_dist_all %>% filter(family == "TAAR") %>% pull(mystat))$p.value

wilcox.test(Rho_vs_dist_all %>% filter(family == "OR") %>% pull(mystat),
            Rho_vs_dist_all %>% filter(family == "V1R") %>% pull(mystat))$p.value

wilcox.test(Rho_vs_dist_all %>% filter(family == "OR") %>% pull(mystat),
            Rho_vs_dist_all %>% filter(family == "V2R") %>% pull(mystat))$p.value

wilcox.test(Rho_vs_dist_all %>% filter(family == "TAAR") %>% pull(mystat),
            Rho_vs_dist_all %>% filter(family == "V1R") %>% pull(mystat))$p.value


wilcox.test(Rho_vs_dist_all %>% filter(family == "TAAR") %>% pull(mystat),
            Rho_vs_dist_all %>% filter(family == "V2R") %>% pull(mystat))$p.value


wilcox.test(Rho_vs_dist_all %>% filter(family == "V1R") %>% pull(mystat),
            Rho_vs_dist_all %>% filter(family == "V2R") %>% pull(mystat))$p.value



#### Expression evolutionnary rate --  Boostrap OR vs TAAR   ---------------------------------

length(Summary_RNAseq_table_OR_relative_species_subfam %>% pull(subfamily) %>% unique())
length(Summary_RNAseq_table_TAAR_relative_species_subfam %>% pull(subfamily) %>% unique())
length(Summary_RNAseq_table_V1R_relative_species_subfam %>% pull(subfamily) %>% unique())
length(Summary_RNAseq_table_V2R_relative_species_subfam %>% pull(subfamily) %>% unique())


#Lets define unique pairs of species

species_list <- 
  Summary_RNAseq_table %>% 
  filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>%
  pull(Sp) %>% unique()


all_comb <- apply(combn(species_list,2),2,paste,collapse=',')
pairs_df <- as.data.frame(t(sapply(all_comb, function(x) unlist(strsplit(x, ",")))))
row.names(pairs_df) <- NULL
colnames(pairs_df) <- c("species_1", "species_2")


#Perform the same test as before but taking 7 subfamilies at random to compare with TAAR genes



Rho_vs_dist_OR_boostrap_7 <- as.data.frame(NULL)

for(i in 1:2){ #put 1000 for the real run
  
  random_OR_subfams <- 
    sample(Summary_RNAseq_table_OR_relative_species_subfam %>% 
             pull(subfamily) %>% 
             unique(), 7)
  
  spearman_rho_list <- c()
  distance_list <- c()
  
  for(pair_nb in 1:nrow(pairs_df)){
    
    sp1 <- pairs_df[pair_nb,]$species_1
    sp2 <- pairs_df[pair_nb,]$species_2
    
    sp1_expression <-
      Summary_RNAseq_table_OR_relative_species_subfam %>%
      filter(subfamily %in% random_OR_subfams) %>%
      filter(Sp == sp1) %>%
      pull(subfam_relative_expr)
    
    sp2_expression <-
      Summary_RNAseq_table_OR_relative_species_subfam %>%
      filter(subfamily %in% random_OR_subfams) %>%
      filter(Sp == sp2) %>%
      pull(subfam_relative_expr)
    
    rho_spearman <- 
      as.numeric(cor.test(x=sp1_expression, y=sp2_expression, method="spearman")$estimate)
    
    spearman_rho_list <- 
      c(spearman_rho_list, rho_spearman)
    
    curr_dist <- 
      (get_pairwise_distances(radiation_tree_wo_Neospl, 
                              sp1, 
                              sp2, 
                              as_edge_counts=FALSE, 
                              check_input=TRUE) / 2)
    
    distance_list <- 
      c(distance_list, curr_dist)
    
    
  }
  
  
  curr_Rho_vs_dist_OR <-
    as.data.frame(
      cbind(spearman_rho_list, distance_list)
    )
  colnames(curr_Rho_vs_dist_OR) <- c("rho", "dist")
  curr_Rho_vs_dist_OR <- curr_Rho_vs_dist_OR %>% mutate(bootstrap_iter = i)
  
  Rho_vs_dist_OR_boostrap_7 <- rbind(Rho_vs_dist_OR_boostrap_7, curr_Rho_vs_dist_OR)
  
}



#Now lets see if it is significantly different than the evolution of TAAR genes


Rho_vs_dist_OR_boostrap_7 <- 
  Rho_vs_dist_OR_boostrap_7 %>%
  mutate(mystat = (1-rho)/(dist)) %>%
  mutate(family = "OR")


#write.table(Rho_vs_dist_OR_boostrap_7,
#            "Rho_vs_dist_OR_boostrap_7.csv",
#            sep=",",
#            col.names=TRUE,
#            row.names=FALSE,
#            quote=FALSE)


Rho_vs_dist_OR_boostrap_7 <- 
  read.table(
    "Rho_vs_dist_OR_boostrap_7.csv",
    sep=",",
    header=TRUE
  )


Boostraps_wilcox_OR_TAAR <- as.data.frame(NULL)
for(i in 1:1000){
  
  wilcox_test <- 
    wilcox.test(
      Rho_vs_dist_OR_boostrap_7 %>%
        filter(bootstrap_iter == i) %>%
        pull(mystat),
      Rho_vs_dist_all %>%
        filter(family == "TAAR") %>%
        pull(mystat)
    )
  curr_pvalue <- wilcox_test$p.value
  
  mean_OR <- Rho_vs_dist_OR_boostrap_7 %>% filter(bootstrap_iter == i) %>% pull(mystat) %>% mean()
  mean_TAAR <- Rho_vs_dist_all %>% filter(family == "TAAR") %>% pull(mystat) %>% mean()
  
  
  curr_df <- as.data.frame(cbind(mean_OR, mean_TAAR, curr_pvalue, i))
  
  Boostraps_wilcox_OR_TAAR <- rbind(Boostraps_wilcox_OR_TAAR, curr_df)
  
}

Boostraps_wilcox_OR_TAAR$curr_pvalue <- as.numeric(Boostraps_wilcox_OR_TAAR$curr_pvalue)
Boostraps_wilcox_OR_TAAR$mean_OR <- as.numeric(Boostraps_wilcox_OR_TAAR$mean_OR)
Boostraps_wilcox_OR_TAAR$mean_TAAR <- as.numeric(Boostraps_wilcox_OR_TAAR$mean_TAAR)




Boostraps_wilcox_OR_TAAR <- 
  Boostraps_wilcox_OR_TAAR %>%
  mutate(mean_OR_TAAR = if_else(
    mean_OR < mean_TAAR,
    "OR<TAAR",
    "OR>TAAR"
  )) %>%
  mutate(significant = if_else(
    curr_pvalue < 0.05,
    "Significant",
    "Non-significant"
  ))


Boostraps_wilcox_OR_TAAR %>% 
  group_by(mean_OR_TAAR, significant) %>%
  summarise(n())

#### Expression evolutionnary rate --  Boostrap OR vs V2R   ---------------------------------

#Lets define unique pairs of species

species_list <- 
  Summary_RNAseq_table %>% 
  filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>%
  pull(Sp) %>% unique()


all_comb <- apply(combn(species_list,2),2,paste,collapse=',')
pairs_df <- as.data.frame(t(sapply(all_comb, function(x) unlist(strsplit(x, ",")))))
row.names(pairs_df) <- NULL
colnames(pairs_df) <- c("species_1", "species_2")


#Perform the same test as before but taking 14 subfamilies at random to compare with V2R genes



Rho_vs_dist_OR_boostrap_14 <- as.data.frame(NULL)

for(i in 1:2){ #put 1000 for the real run
  
  random_OR_subfams <- 
    sample(Summary_RNAseq_table_OR_relative_species_subfam %>% 
             pull(subfamily) %>% 
             unique(), 14)
  
  spearman_rho_list <- c()
  distance_list <- c()
  
  for(pair_nb in 1:nrow(pairs_df)){
    
    sp1 <- pairs_df[pair_nb,]$species_1
    sp2 <- pairs_df[pair_nb,]$species_2
    
    sp1_expression <-
      Summary_RNAseq_table_OR_relative_species_subfam %>%
      filter(subfamily %in% random_OR_subfams) %>%
      filter(Sp == sp1) %>%
      pull(subfam_relative_expr)
    
    sp2_expression <-
      Summary_RNAseq_table_OR_relative_species_subfam %>%
      filter(subfamily %in% random_OR_subfams) %>%
      filter(Sp == sp2) %>%
      pull(subfam_relative_expr)
    
    rho_spearman <- 
      as.numeric(cor.test(x=sp1_expression, y=sp2_expression, method="spearman")$estimate)
    
    spearman_rho_list <- 
      c(spearman_rho_list, rho_spearman)
    
    curr_dist <- 
      (get_pairwise_distances(radiation_tree_wo_Neospl, 
                              sp1, 
                              sp2, 
                              as_edge_counts=FALSE, 
                              check_input=TRUE) / 2)
    
    distance_list <- 
      c(distance_list, curr_dist)
    
    
  }
  
  
  curr_Rho_vs_dist_OR <-
    as.data.frame(
      cbind(spearman_rho_list, distance_list)
    )
  colnames(curr_Rho_vs_dist_OR) <- c("rho", "dist")
  curr_Rho_vs_dist_OR <- curr_Rho_vs_dist_OR %>% mutate(bootstrap_iter = i)
  
  Rho_vs_dist_OR_boostrap_14 <- rbind(Rho_vs_dist_OR_boostrap_14, curr_Rho_vs_dist_OR)
  
}



#Now lets see if it is significantly different than the evolution of V2R genes


Rho_vs_dist_OR_boostrap_14 <- 
  Rho_vs_dist_OR_boostrap_14 %>%
  mutate(mystat = (1-rho)/(dist)) %>%
  mutate(family = "OR")


#write.table(Rho_vs_dist_OR_boostrap_14,
#            "Rho_vs_dist_OR_boostrap_14.csv",
#            sep=",",
#            col.names=TRUE,
#            row.names=FALSE,
#            quote=FALSE)


Rho_vs_dist_OR_boostrap_14 <- 
  read.table(
    "Rho_vs_dist_OR_boostrap_14.csv",
    sep=",",
    header=TRUE
  )


Boostraps_wilcox_OR_V2R <- as.data.frame(NULL)
for(i in 1:1000){
  
  wilcox_test <- 
    wilcox.test(
      Rho_vs_dist_OR_boostrap_14 %>%
        filter(bootstrap_iter == i) %>%
        pull(mystat),
      Rho_vs_dist_all %>%
        filter(family == "V2R") %>%
        pull(mystat)
    )
  curr_pvalue <- wilcox_test$p.value
  
  mean_OR <- Rho_vs_dist_OR_boostrap_14 %>% filter(bootstrap_iter == i) %>% pull(mystat) %>% mean()
  mean_V2R <- Rho_vs_dist_all %>% filter(family == "V2R") %>% pull(mystat) %>% mean()
  
  
  curr_df <- as.data.frame(cbind(mean_OR, mean_V2R, curr_pvalue, i))
  
  Boostraps_wilcox_OR_V2R <- rbind(Boostraps_wilcox_OR_V2R, curr_df)
  
}

Boostraps_wilcox_OR_V2R$curr_pvalue <- as.numeric(Boostraps_wilcox_OR_V2R$curr_pvalue)
Boostraps_wilcox_OR_V2R$mean_OR <- as.numeric(Boostraps_wilcox_OR_V2R$mean_OR)
Boostraps_wilcox_OR_V2R$mean_V2R <- as.numeric(Boostraps_wilcox_OR_V2R$mean_V2R)




Boostraps_wilcox_OR_V2R <- 
  Boostraps_wilcox_OR_V2R %>%
  mutate(mean_OR_V2R = if_else(
    mean_OR < mean_V2R,
    "OR<V2R",
    "OR>V2R"
  )) %>%
  mutate(significant = if_else(
    curr_pvalue < 0.05,
    "Significant",
    "Non-significant"
  ))


Boostraps_wilcox_OR_V2R %>% 
  group_by(mean_OR_V2R, significant) %>%
  summarise(n())


#### Expression evolutionnary rate --  Boostrap OR vs V1R   ---------------------------------

#Lets define unique pairs of species

species_list <- 
  Summary_RNAseq_table %>% 
  filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>%
  pull(Sp) %>% unique()


all_comb <- apply(combn(species_list,2),2,paste,collapse=',')
pairs_df <- as.data.frame(t(sapply(all_comb, function(x) unlist(strsplit(x, ",")))))
row.names(pairs_df) <- NULL
colnames(pairs_df) <- c("species_1", "species_2")


#Perform the same test as before but taking 6 subfamilies at random to compare with V1R genes



Rho_vs_dist_OR_boostrap_6 <- as.data.frame(NULL)

for(i in 1:2){ #put 1000 for the real run
  
  random_OR_subfams <- 
    sample(Summary_RNAseq_table_OR_relative_species_subfam %>% 
             pull(subfamily) %>% 
             unique(), 6)
  
  spearman_rho_list <- c()
  distance_list <- c()
  
  for(pair_nb in 1:nrow(pairs_df)){
    
    sp1 <- pairs_df[pair_nb,]$species_1
    sp2 <- pairs_df[pair_nb,]$species_2
    
    sp1_expression <-
      Summary_RNAseq_table_OR_relative_species_subfam %>%
      filter(subfamily %in% random_OR_subfams) %>%
      filter(Sp == sp1) %>%
      pull(subfam_relative_expr)
    
    sp2_expression <-
      Summary_RNAseq_table_OR_relative_species_subfam %>%
      filter(subfamily %in% random_OR_subfams) %>%
      filter(Sp == sp2) %>%
      pull(subfam_relative_expr)
    
    rho_spearman <- 
      as.numeric(cor.test(x=sp1_expression, y=sp2_expression, method="spearman")$estimate)
    
    spearman_rho_list <- 
      c(spearman_rho_list, rho_spearman)
    
    curr_dist <- 
      (get_pairwise_distances(radiation_tree_wo_Neospl, 
                              sp1, 
                              sp2, 
                              as_edge_counts=FALSE, 
                              check_input=TRUE) / 2)
    
    distance_list <- 
      c(distance_list, curr_dist)
    
    
  }
  
  
  curr_Rho_vs_dist_OR <-
    as.data.frame(
      cbind(spearman_rho_list, distance_list)
    )
  colnames(curr_Rho_vs_dist_OR) <- c("rho", "dist")
  curr_Rho_vs_dist_OR <- curr_Rho_vs_dist_OR %>% mutate(bootstrap_iter = i)
  
  Rho_vs_dist_OR_boostrap_6 <- rbind(Rho_vs_dist_OR_boostrap_6, curr_Rho_vs_dist_OR)
  
}


#Now lets see if it is significantly different than the evolution of V1R genes


Rho_vs_dist_OR_boostrap_6 <- 
  Rho_vs_dist_OR_boostrap_6 %>%
  mutate(mystat = (1-rho)/(dist)) %>%
  mutate(family = "OR")


#write.table(Rho_vs_dist_OR_boostrap_6,
#            "Rho_vs_dist_OR_boostrap_6.csv",
#            sep=",",
#            col.names=TRUE,
#            row.names=FALSE,
#            quote=FALSE)


Rho_vs_dist_OR_boostrap_6 <- 
  read.table(
    "Rho_vs_dist_OR_boostrap_6.csv",
    sep=",",
    header=TRUE
  )


Boostraps_wilcox_OR_V1R <- as.data.frame(NULL)
for(i in 1:1000){
  
  wilcox_test <- 
    wilcox.test(
      Rho_vs_dist_OR_boostrap_6 %>%
        filter(bootstrap_iter == i) %>%
        pull(mystat),
      Rho_vs_dist_all %>%
        filter(family == "V1R") %>%
        pull(mystat)
    )
  curr_pvalue <- wilcox_test$p.value
  
  mean_OR <- Rho_vs_dist_OR_boostrap_6 %>% filter(bootstrap_iter == i) %>% pull(mystat) %>% mean()
  mean_V1R <- Rho_vs_dist_all %>% filter(family == "V1R") %>% pull(mystat) %>% mean()
  
  
  curr_df <- as.data.frame(cbind(mean_OR, mean_V1R, curr_pvalue, i))
  
  Boostraps_wilcox_OR_V1R <- rbind(Boostraps_wilcox_OR_V1R, curr_df)
  
}

Boostraps_wilcox_OR_V1R$curr_pvalue <- as.numeric(Boostraps_wilcox_OR_V1R$curr_pvalue)
Boostraps_wilcox_OR_V1R$mean_OR <- as.numeric(Boostraps_wilcox_OR_V1R$mean_OR)
Boostraps_wilcox_OR_V1R$mean_V1R <- as.numeric(Boostraps_wilcox_OR_V1R$mean_V1R)




Boostraps_wilcox_OR_V1R <- 
  Boostraps_wilcox_OR_V1R %>%
  mutate(mean_OR_V1R = if_else(
    mean_OR < mean_V1R,
    "OR<V1R",
    "OR>V1R"
  )) %>%
  mutate(significant = if_else(
    curr_pvalue < 0.05,
    "Significant",
    "Non-significant"
  ))


Boostraps_wilcox_OR_V1R %>% 
  group_by(mean_OR_V1R, significant) %>%
  summarise(n())

#### Expression evolutionnary rate --  Boostrap V2R vs V1R   ---------------------------------

#Lets define unique pairs of species

species_list <- 
  Summary_RNAseq_table %>% 
  filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>%
  pull(Sp) %>% unique()


all_comb <- apply(combn(species_list,2),2,paste,collapse=',')
pairs_df <- as.data.frame(t(sapply(all_comb, function(x) unlist(strsplit(x, ",")))))
row.names(pairs_df) <- NULL
colnames(pairs_df) <- c("species_1", "species_2")


#Perform the same test as before but taking 6 subfamilies at random to compare with V1R genes



Rho_vs_dist_V2R_boostrap_6 <- as.data.frame(NULL)

for(i in 1:2){ #put 1000 for the real run
  
  random_V2R_subfams <- 
    sample(Summary_RNAseq_table_V2R_relative_species_subfam %>% 
             pull(subfamily) %>% 
             unique(), 6)
  
  spearman_rho_list <- c()
  distance_list <- c()
  
  for(pair_nb in 1:nrow(pairs_df)){
    
    sp1 <- pairs_df[pair_nb,]$species_1
    sp2 <- pairs_df[pair_nb,]$species_2
    
    sp1_expression <-
      Summary_RNAseq_table_V2R_relative_species_subfam %>%
      filter(subfamily %in% random_V2R_subfams) %>%
      filter(Sp == sp1) %>%
      pull(subfam_relative_expr)
    
    sp2_expression <-
      Summary_RNAseq_table_V2R_relative_species_subfam %>%
      filter(subfamily %in% random_V2R_subfams) %>%
      filter(Sp == sp2) %>%
      pull(subfam_relative_expr)
    
    rho_spearman <- 
      as.numeric(cor.test(x=sp1_expression, y=sp2_expression, method="spearman")$estimate)
    
    spearman_rho_list <- 
      c(spearman_rho_list, rho_spearman)
    
    curr_dist <- 
      (get_pairwise_distances(radiation_tree_wo_Neospl, 
                              sp1, 
                              sp2, 
                              as_edge_counts=FALSE, 
                              check_input=TRUE) / 2)
    
    distance_list <- 
      c(distance_list, curr_dist)
    
    
  }
  
  
  curr_Rho_vs_dist_V2R <-
    as.data.frame(
      cbind(spearman_rho_list, distance_list)
    )
  colnames(curr_Rho_vs_dist_V2R) <- c("rho", "dist")
  curr_Rho_vs_dist_V2R <- curr_Rho_vs_dist_V2R %>% mutate(bootstrap_iter = i)
  
  Rho_vs_dist_V2R_boostrap_6 <- rbind(Rho_vs_dist_V2R_boostrap_6, curr_Rho_vs_dist_V2R)
  
}


#Now lets see if it is significantly different than the evolution of V1R genes


Rho_vs_dist_V2R_boostrap_6 <- 
  Rho_vs_dist_V2R_boostrap_6 %>%
  mutate(mystat = (1-rho)/(dist)) %>%
  mutate(family = "V2R")


#write.table(Rho_vs_dist_V2R_boostrap_6,
#            "Rho_vs_dist_V2R_boostrap_6.csv",
#            sep=",",
#            col.names=TRUE,
#            row.names=FALSE,
#            quote=FALSE)
#

Rho_vs_dist_V2R_boostrap_6 <- 
  read.table(
    "Rho_vs_dist_V2R_boostrap_6.csv",
    sep=",",
    header=TRUE
  )


Boostraps_wilcox_V2R_V1R <- as.data.frame(NULL)
for(i in 1:1000){
  
  wilcox_test <- 
    wilcox.test(
      Rho_vs_dist_V2R_boostrap_6 %>%
        filter(bootstrap_iter == i) %>%
        pull(mystat),
      Rho_vs_dist_all %>%
        filter(family == "V1R") %>%
        pull(mystat)
    )
  curr_pvalue <- wilcox_test$p.value
  
  mean_V2R <- Rho_vs_dist_V2R_boostrap_6 %>% filter(bootstrap_iter == i) %>% pull(mystat) %>% mean()
  mean_V1R <- Rho_vs_dist_all %>% filter(family == "V1R") %>% pull(mystat) %>% mean()
  
  
  curr_df <- as.data.frame(cbind(mean_V2R, mean_V1R, curr_pvalue, i))
  
  Boostraps_wilcox_V2R_V1R <- rbind(Boostraps_wilcox_V2R_V1R, curr_df)
  
}

Boostraps_wilcox_V2R_V1R$curr_pvalue <- as.numeric(Boostraps_wilcox_V2R_V1R$curr_pvalue)
Boostraps_wilcox_V2R_V1R$mean_V2R <- as.numeric(Boostraps_wilcox_V2R_V1R$mean_V2R)
Boostraps_wilcox_V2R_V1R$mean_V1R <- as.numeric(Boostraps_wilcox_V2R_V1R$mean_V1R)




Boostraps_wilcox_V2R_V1R <- 
  Boostraps_wilcox_V2R_V1R %>%
  mutate(mean_V2R_V1R = if_else(
    mean_V2R < mean_V1R,
    "V2R<V1R",
    "V2R>V1R"
  )) %>%
  mutate(significant = if_else(
    curr_pvalue < 0.05,
    "Significant",
    "Non-significant"
  ))


Boostraps_wilcox_V2R_V1R %>% 
  group_by(mean_V2R_V1R, significant) %>%
  summarise(n())

#### Expression evolutionnary rate --  Boostrap V2R vs TAAR   ---------------------------------

#Lets define unique pairs of species

species_list <- 
  Summary_RNAseq_table %>% 
  filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>%
  pull(Sp) %>% unique()


all_comb <- apply(combn(species_list,2),2,paste,collapse=',')
pairs_df <- as.data.frame(t(sapply(all_comb, function(x) unlist(strsplit(x, ",")))))
row.names(pairs_df) <- NULL
colnames(pairs_df) <- c("species_1", "species_2")


#Perform the same test as before but taking 7 subfamilies at random to compare with TAAR genes



Rho_vs_dist_V2R_boostrap_7 <- as.data.frame(NULL)

for(i in 1:1000){ #put 1000 for the real run
  
  random_V2R_subfams <- 
    sample(Summary_RNAseq_table_V2R_relative_species_subfam %>% 
             pull(subfamily) %>% 
             unique(), 7)
  
  spearman_rho_list <- c()
  distance_list <- c()
  
  for(pair_nb in 1:nrow(pairs_df)){
    
    sp1 <- pairs_df[pair_nb,]$species_1
    sp2 <- pairs_df[pair_nb,]$species_2
    
    sp1_expression <-
      Summary_RNAseq_table_V2R_relative_species_subfam %>%
      filter(subfamily %in% random_V2R_subfams) %>%
      filter(Sp == sp1) %>%
      pull(subfam_relative_expr)
    
    sp2_expression <-
      Summary_RNAseq_table_V2R_relative_species_subfam %>%
      filter(subfamily %in% random_V2R_subfams) %>%
      filter(Sp == sp2) %>%
      pull(subfam_relative_expr)
    
    rho_spearman <- 
      as.numeric(cor.test(x=sp1_expression, y=sp2_expression, method="spearman")$estimate)
    
    spearman_rho_list <- 
      c(spearman_rho_list, rho_spearman)
    
    curr_dist <- 
      (get_pairwise_distances(radiation_tree_wo_Neospl, 
                              sp1, 
                              sp2, 
                              as_edge_counts=FALSE, 
                              check_input=TRUE) / 2)
    
    distance_list <- 
      c(distance_list, curr_dist)
    
    
  }
  
  
  curr_Rho_vs_dist_V2R <-
    as.data.frame(
      cbind(spearman_rho_list, distance_list)
    )
  colnames(curr_Rho_vs_dist_V2R) <- c("rho", "dist")
  curr_Rho_vs_dist_V2R <- curr_Rho_vs_dist_V2R %>% mutate(bootstrap_iter = i)
  
  Rho_vs_dist_V2R_boostrap_7 <- rbind(Rho_vs_dist_V2R_boostrap_7, curr_Rho_vs_dist_V2R)
  
}



#Now lets see if it is significantly different than the evolution of TAAR genes


Rho_vs_dist_V2R_boostrap_7 <- 
  Rho_vs_dist_V2R_boostrap_7 %>%
  mutate(mystat = (1-rho)/(dist)) %>%
  mutate(family = "V2R")


#write.table(Rho_vs_dist_V2R_boostrap_7,
#            "Rho_vs_dist_V2R_boostrap_7.csv",
#            sep=",",
#            col.names=TRUE,
#            row.names=FALSE,
#            quote=FALSE)


Rho_vs_dist_V2R_boostrap_7 <- 
  read.table(
    "Rho_vs_dist_V2R_boostrap_7.csv",
    sep=",",
    header=TRUE
  )


Boostraps_wilcox_V2R_TAAR <- as.data.frame(NULL)
for(i in 1:1000){
  
  wilcox_test <- 
    wilcox.test(
      Rho_vs_dist_V2R_boostrap_7 %>%
        filter(bootstrap_iter == i) %>%
        pull(mystat),
      Rho_vs_dist_all %>%
        filter(family == "TAAR") %>%
        pull(mystat)
    )
  curr_pvalue <- wilcox_test$p.value
  
  mean_V2R <- Rho_vs_dist_V2R_boostrap_7 %>% filter(bootstrap_iter == i) %>% pull(mystat) %>% mean()
  mean_TAAR <- Rho_vs_dist_all %>% filter(family == "TAAR") %>% pull(mystat) %>% mean()
  
  
  curr_df <- as.data.frame(cbind(mean_V2R, mean_TAAR, curr_pvalue, i))
  
  Boostraps_wilcox_V2R_TAAR <- rbind(Boostraps_wilcox_V2R_TAAR, curr_df)
  
}

Boostraps_wilcox_V2R_TAAR$curr_pvalue <- as.numeric(Boostraps_wilcox_V2R_TAAR$curr_pvalue)
Boostraps_wilcox_V2R_TAAR$mean_V2R <- as.numeric(Boostraps_wilcox_V2R_TAAR$mean_V2R)
Boostraps_wilcox_V2R_TAAR$mean_TAAR <- as.numeric(Boostraps_wilcox_V2R_TAAR$mean_TAAR)




Boostraps_wilcox_V2R_TAAR <- 
  Boostraps_wilcox_V2R_TAAR %>%
  mutate(mean_V2R_TAAR = if_else(
    mean_V2R < mean_TAAR,
    "V2R<TAAR",
    "V2R>TAAR"
  )) %>%
  mutate(significant = if_else(
    curr_pvalue < 0.05,
    "Significant",
    "Non-significant"
  ))


Boostraps_wilcox_V2R_TAAR %>% 
  group_by(mean_V2R_TAAR, significant) %>%
  summarise(n())
