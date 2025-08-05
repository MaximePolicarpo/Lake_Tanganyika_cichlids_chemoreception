#### Libraries  ---------------------------------

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

#### Data load  ---------------------------------


#load db infos on cichlids
db_infos <- read.table("DNATube_2021-03-04_16-15-42_Simplified_parsed.tsv", sep="\t", header=TRUE)
colnames(db_infos) <- c("Sp", "Species", "ID", "Sex", "Tribe")
sp.colors <- c(Batmin = "#242626" , Cypfro = "#FDDF13" , Cyplep = "#F04D29" , Cunlon = "#9AB9D9" , Neomul = "#C588BB" , Simdia = "#86C773")


#load phylo tree
b1_tree <- read.tree("b1_tree.nwk")
b1_list_sp <- b1_tree$tip.label

#load pacbio mining results
pacbio_olfactory_receptor <- read.table("ALL_Chemoreceptors_classification.csv", header=FALSE, sep=",")
colnames(pacbio_olfactory_receptor) <- 
  c("Gene_type", "Gene_name", "Subfamily", "Clade", "OGG")

sequences_names <- pacbio_olfactory_receptor %>% pull(Gene_name)
all_genome_names <- c("Bathybates_minor","Cunningtonia_longiventralis","Cyphotilapia_frontosa","Cyprichromis_leptosoma","Neolamprologus_multifasciatus","Oreochromis_niloticus","Simochromis_diagramma")

corr_tips_species <- as.data.frame(NULL)
for(my_sp in all_genome_names){
  
  curr_sp_df <- 
    as.data.frame(grep(paste("^",my_sp, "---", sep=""), 
                       sequences_names, value = TRUE)) %>%
    mutate(genome_name = my_sp)
  colnames(curr_sp_df) <- c("Gene_name", "Species")
  
  corr_tips_species <- rbind(corr_tips_species, curr_sp_df)
  
}
pacbio_olfactory_receptor_sp <- 
  left_join(pacbio_olfactory_receptor, corr_tips_species,  by="Gene_name")
pacbio_olfactory_receptor <- pacbio_olfactory_receptor_sp



#### Compute number of genes per ID from depth of coverage of OR genes  ---------------------------------

raw_coverage_OR <- 
  read.table("ALL.OR.unmasked.raw.depth")
colnames(raw_coverage_OR) <- c("ID", "scaffold", "pos", "coverage")
mean_coverage_OR <- 
  read.table("ALL.OR.unmasked.mean.depth")
colnames(mean_coverage_OR) <- c("Gene_type", "Gene_name", "Clade", "Subfamily", 
                                "OGG", "scaffold", "start", "end", "ID", "coverage")
coverage_BUSCO <- 
  read.table("ALL.AllBUSCO.unmasked.depth")
colnames(coverage_BUSCO) <- c("ID", "coverage")


#all_ID <- scan("/scicore/home/salzburg/polica0000/Reads_Depth_Of_Coverage/good_ID_list.txt", what="character")
all_ID <- scan("good_ID_list.txt", what="character")

#stats
left_join(
  all_ID


OR_DoC_df <- as.data.frame(NULL)
for(curr_ID in all_ID){
  mean_OR_cov <- raw_coverage_OR %>% filter(ID == curr_ID) %>% pull(coverage) %>% mean()
  mean_BUSCO_cov <- coverage_BUSCO %>% filter(ID == curr_ID) %>% pull(coverage)
  inf_nb_OR <-  mean_OR_cov/mean_BUSCO_cov * 336
  
  
  normalized_nb_genes <- 
    mean_coverage_OR %>%
    filter(ID == curr_ID) %>%
    group_by(Gene_name) %>%
    mutate(cov_per_gene = mean(coverage)) %>%
    dplyr::select(- c(coverage,scaffold,start,end)) %>%
    distinct() %>%
    ungroup() %>%
    group_by(Subfamily) %>%
    mutate(cov_per_subfamily = mean(cov_per_gene)) %>%
    summarise(nb_onil_gene = n(), cov=cov_per_subfamily) %>%
    mutate(nb_gene = cov/mean_BUSCO_cov*nb_onil_gene) %>%
    distinct() %>%
    pull(nb_gene) %>% sum()
  
  Subfamily_nb <- 
    as.data.frame(
      mean_coverage_OR %>%
        filter(ID == curr_ID) %>%
        group_by(Gene_name) %>%
        mutate(cov_per_gene = mean(coverage)) %>%
        dplyr::select(- c(coverage,scaffold,start,end)) %>%
        distinct() %>%
        ungroup() %>%
        group_by(Subfamily) %>%
        mutate(cov_per_subfamily = mean(cov_per_gene)) %>%
        summarise(nb_onil_gene = n(), cov=cov_per_subfamily) %>%
        mutate(nb_gene = cov/mean_BUSCO_cov*nb_onil_gene) %>%
        distinct() %>%
        mutate(normalized_nb = (nb_gene*inf_nb_OR)/normalized_nb_genes) %>%
        dplyr::select(Subfamily, normalized_nb)
    )
  
  
  
  Subfamily_nb <- 
    Subfamily_nb %>%
    mutate(ID = curr_ID)
  
  OR_DoC_df <- 
    rbind(OR_DoC_df,
          Subfamily_nb)
  
}


rm(raw_coverage_OR)

#Add the number of Total OR genes
for(curr_ID in all_ID){
  total <- OR_DoC_df %>% filter(ID == curr_ID) %>% pull(normalized_nb) %>% sum()
  temp_df <- as.data.frame("Total_OR")
  temp_df <- cbind(temp_df, as.data.frame(total), as.data.frame(curr_ID))
  colnames(temp_df) <- c("Subfamily", "normalized_nb", "ID")
  
  OR_DoC_df <- rbind(OR_DoC_df, temp_df)
}




#### Compute number of genes per species from depth of coverage of OR genes  ---------------------------------

OR_DoC_df_info <- left_join(OR_DoC_df, db_infos, by="ID")

OR_DoC_df_info_mean <- 
  as.data.frame(
    OR_DoC_df_info %>%
      group_by(Species, Subfamily) %>%
      mutate(mean_normalized_nb = mean(normalized_nb))
  ) %>%
  dplyr::select(Subfamily, mean_normalized_nb, Sp, Species, Tribe) %>%
  distinct()



write.table(OR_DoC_df_info, "OR_DoC_df_info.csv", quote=FALSE,
            col.names=TRUE, row.names=FALSE, sep=",")

write.table(OR_DoC_df_info_mean, "OR_DoC_df_info_mean.csv", quote=FALSE,
            col.names=TRUE, row.names=FALSE, sep=",")


#### Compute number of genes per ID from depth of coverage of V2R genes  ---------------------------------

raw_coverage_V2R <- 
  read.table("ALL.V2R.unmasked.raw.depth")
colnames(raw_coverage_V2R) <- c("ID", "scaffold", "pos", "coverage")
mean_coverage_V2R <- 
  read.table("ALL.V2R.unmasked.mean.depth")
colnames(mean_coverage_V2R) <- c("Gene_type", "Gene_name", "Clade", "Subfamily", 
                                 "OGG", "scaffold", "start", "end", "ID", "coverage")
coverage_BUSCO <- 
  read.table("ALL.AllBUSCO.unmasked.depth")
colnames(coverage_BUSCO) <- c("ID", "coverage")


#all_ID <- scan("/scicore/home/salzburg/polica0000/Reads_Depth_Of_Coverage/good_ID_list.txt", what="character")
all_ID <- scan("good_ID_list.txt", what="character")

V2R_DoC_df <- as.data.frame(NULL)
for(curr_ID in all_ID){
  mean_V2R_cov <- raw_coverage_V2R %>% filter(ID == curr_ID) %>% pull(coverage) %>% mean()
  mean_BUSCO_cov <- coverage_BUSCO %>% filter(ID == curr_ID) %>% pull(coverage)
  inf_nb_V2R <-  mean_V2R_cov/mean_BUSCO_cov * 85
  
  
  normalized_nb_genes <- 
    mean_coverage_V2R %>%
    filter(ID == curr_ID) %>%
    group_by(Gene_name) %>%
    mutate(cov_per_gene = mean(coverage)) %>%
    dplyr::select(- c(coverage,scaffold,start,end)) %>%
    distinct() %>%
    ungroup() %>%
    group_by(Subfamily) %>%
    mutate(cov_per_subfamily = mean(cov_per_gene)) %>%
    summarise(nb_onil_gene = n(), cov=cov_per_subfamily) %>%
    mutate(nb_gene = cov/mean_BUSCO_cov*nb_onil_gene) %>%
    distinct() %>%
    pull(nb_gene) %>% sum()
  
  Subfamily_nb <- 
    as.data.frame(
      mean_coverage_V2R %>%
        filter(ID == curr_ID) %>%
        group_by(Gene_name) %>%
        mutate(cov_per_gene = mean(coverage)) %>%
        dplyr::select(- c(coverage,scaffold,start,end)) %>%
        distinct() %>%
        ungroup() %>%
        group_by(Subfamily) %>%
        mutate(cov_per_subfamily = mean(cov_per_gene)) %>%
        summarise(nb_onil_gene = n(), cov=cov_per_subfamily) %>%
        mutate(nb_gene = cov/mean_BUSCO_cov*nb_onil_gene) %>%
        distinct() %>%
        mutate(normalized_nb = (nb_gene*inf_nb_V2R)/normalized_nb_genes) %>%
        dplyr::select(Subfamily, normalized_nb)
    )
  
  
  
  Subfamily_nb <- 
    Subfamily_nb %>%
    mutate(ID = curr_ID)
  
  V2R_DoC_df <- 
    rbind(V2R_DoC_df,
          Subfamily_nb)
  
}


rm(raw_coverage_V2R)



#Add the number of Total V2R genes
for(curr_ID in all_ID){
  total <- V2R_DoC_df %>% filter(ID == curr_ID) %>% pull(normalized_nb) %>% sum()
  temp_df <- as.data.frame("Total_V2R")
  temp_df <- cbind(temp_df, as.data.frame(total), as.data.frame(curr_ID))
  colnames(temp_df) <- c("Subfamily", "normalized_nb", "ID")
  
  V2R_DoC_df <- rbind(V2R_DoC_df, temp_df)
}





#### Compute number of genes per species from depth of coverage of V2R genes  ---------------------------------

V2R_DoC_df_info <- left_join(V2R_DoC_df, db_infos, by="ID")

V2R_DoC_df_info_mean <- 
  as.data.frame(
    V2R_DoC_df_info %>%
      group_by(Species, Subfamily) %>%
      mutate(mean_normalized_nb = mean(normalized_nb))
  ) %>%
  dplyr::select(Subfamily, mean_normalized_nb, Sp, Species, Tribe) %>%
  distinct()



write.table(V2R_DoC_df_info, "V2R_DoC_df_info.csv", quote=FALSE,
            col.names=TRUE, row.names=FALSE, sep=",")

write.table(V2R_DoC_df_info_mean, "V2R_DoC_df_info_mean.csv", quote=FALSE,
            col.names=TRUE, row.names=FALSE, sep=",")



#### Compute number of genes per ID from depth of coverage of TAAR genes  ---------------------------------

raw_coverage_TAAR <- 
  read.table("ALL.TAAR.unmasked.raw.depth")
colnames(raw_coverage_TAAR) <- c("ID", "scaffold", "pos", "coverage")
mean_coverage_TAAR <- 
  read.table("ALL.TAAR.unmasked.mean.depth")
colnames(mean_coverage_TAAR) <- c("Gene_type", "Gene_name", "Clade", "Subfamily", 
                                  "OGG", "scaffold", "start", "end", "ID", "coverage")
coverage_BUSCO <- 
  read.table("ALL.AllBUSCO.unmasked.depth")
colnames(coverage_BUSCO) <- c("ID", "coverage")


#all_ID <- scan("/scicore/home/salzburg/polica0000/Reads_Depth_Of_Coverage/good_ID_list.txt", what="character")
all_ID <- scan("good_ID_list.txt", what="character")

TAAR_DoC_df <- as.data.frame(NULL)
for(curr_ID in all_ID){
  mean_TAAR_cov <- raw_coverage_TAAR %>% filter(ID == curr_ID) %>% pull(coverage) %>% mean()
  mean_BUSCO_cov <- coverage_BUSCO %>% filter(ID == curr_ID) %>% pull(coverage)
  inf_nb_TAAR <-  mean_TAAR_cov/mean_BUSCO_cov * 83
  
  
  normalized_nb_genes <- 
    mean_coverage_TAAR %>%
    filter(ID == curr_ID) %>%
    group_by(Gene_name) %>%
    mutate(cov_per_gene = mean(coverage)) %>%
    dplyr::select(- c(coverage,scaffold,start,end)) %>%
    distinct() %>%
    ungroup() %>%
    group_by(Subfamily) %>%
    mutate(cov_per_subfamily = mean(cov_per_gene)) %>%
    summarise(nb_onil_gene = n(), cov=cov_per_subfamily) %>%
    mutate(nb_gene = cov/mean_BUSCO_cov*nb_onil_gene) %>%
    distinct() %>%
    pull(nb_gene) %>% sum()
  
  Subfamily_nb <- 
    as.data.frame(
      mean_coverage_TAAR %>%
        filter(ID == curr_ID) %>%
        group_by(Gene_name) %>%
        mutate(cov_per_gene = mean(coverage)) %>%
        dplyr::select(- c(coverage,scaffold,start,end)) %>%
        distinct() %>%
        ungroup() %>%
        group_by(Subfamily) %>%
        mutate(cov_per_subfamily = mean(cov_per_gene)) %>%
        summarise(nb_onil_gene = n(), cov=cov_per_subfamily) %>%
        mutate(nb_gene = cov/mean_BUSCO_cov*nb_onil_gene) %>%
        distinct() %>%
        mutate(normalized_nb = (nb_gene*inf_nb_TAAR)/normalized_nb_genes) %>%
        dplyr::select(Subfamily, normalized_nb)
    )
  
  
  
  Subfamily_nb <- 
    Subfamily_nb %>%
    mutate(ID = curr_ID)
  
  TAAR_DoC_df <- 
    rbind(TAAR_DoC_df,
          Subfamily_nb)
  
}


rm(raw_coverage_TAAR)



#Add the number of Total TAAR genes
for(curr_ID in all_ID){
  total <- TAAR_DoC_df %>% filter(ID == curr_ID) %>% pull(normalized_nb) %>% sum()
  temp_df <- as.data.frame("Total_TAAR")
  temp_df <- cbind(temp_df, as.data.frame(total), as.data.frame(curr_ID))
  colnames(temp_df) <- c("Subfamily", "normalized_nb", "ID")
  
  TAAR_DoC_df <- rbind(TAAR_DoC_df, temp_df)
}





#### Compute number of genes per species from depth of coverage of TAAR genes  ---------------------------------

TAAR_DoC_df_info <- left_join(TAAR_DoC_df, db_infos, by="ID")

TAAR_DoC_df_info_mean <- 
  as.data.frame(
    TAAR_DoC_df_info %>%
      group_by(Species, Subfamily) %>%
      mutate(mean_normalized_nb = mean(normalized_nb))
  ) %>%
  dplyr::select(Subfamily, mean_normalized_nb, Sp, Species, Tribe) %>%
  distinct()



write.table(TAAR_DoC_df_info, "TAAR_DoC_df_info.csv", quote=FALSE,
            col.names=TRUE, row.names=FALSE, sep=",")

write.table(TAAR_DoC_df_info_mean, "TAAR_DoC_df_info_mean.csv", quote=FALSE,
            col.names=TRUE, row.names=FALSE, sep=",")



#### Compute number of genes per ID from depth of coverage of V1R genes  ---------------------------------

raw_coverage_V1R <- 
  read.table("ALL.V1R.unmasked.raw.depth")
colnames(raw_coverage_V1R) <- c("ID", "scaffold", "pos", "coverage")
mean_coverage_V1R <- 
  read.table("ALL.V1R.unmasked.mean.depth")
colnames(mean_coverage_V1R) <- c("Gene_type", "Gene_name", "Clade", "Subfamily", 
                                 "OGG", "scaffold", "start", "end", "ID", "coverage")
coverage_BUSCO <- 
  read.table("ALL.AllBUSCO.unmasked.depth")
colnames(coverage_BUSCO) <- c("ID", "coverage")


#all_ID <- scan("/scicore/home/salzburg/polica0000/Reads_Depth_Of_Coverage/good_ID_list.txt", what="character")
all_ID <- scan("good_ID_list.txt", what="character")

V1R_DoC_df <- as.data.frame(NULL)
for(curr_ID in all_ID){
  mean_V1R_cov <- raw_coverage_V1R %>% filter(ID == curr_ID) %>% pull(coverage) %>% mean()
  mean_BUSCO_cov <- coverage_BUSCO %>% filter(ID == curr_ID) %>% pull(coverage)
  inf_nb_V1R <-  mean_V1R_cov/mean_BUSCO_cov * 6
  
  
  normalized_nb_genes <- 
    mean_coverage_V1R %>%
    filter(ID == curr_ID) %>%
    group_by(Gene_name) %>%
    mutate(cov_per_gene = mean(coverage)) %>%
    dplyr::select(- c(coverage,scaffold,start,end)) %>%
    distinct() %>%
    ungroup() %>%
    group_by(Subfamily) %>%
    mutate(cov_per_subfamily = mean(cov_per_gene)) %>%
    summarise(nb_onil_gene = n(), cov=cov_per_subfamily) %>%
    mutate(nb_gene = cov/mean_BUSCO_cov*nb_onil_gene) %>%
    distinct() %>%
    pull(nb_gene) %>% sum()
  
  Subfamily_nb <- 
    as.data.frame(
      mean_coverage_V1R %>%
        filter(ID == curr_ID) %>%
        group_by(Gene_name) %>%
        mutate(cov_per_gene = mean(coverage)) %>%
        dplyr::select(- c(coverage,scaffold,start,end)) %>%
        distinct() %>%
        ungroup() %>%
        group_by(Subfamily) %>%
        mutate(cov_per_subfamily = mean(cov_per_gene)) %>%
        summarise(nb_onil_gene = n(), cov=cov_per_subfamily) %>%
        mutate(nb_gene = cov/mean_BUSCO_cov*nb_onil_gene) %>%
        distinct() %>%
        mutate(normalized_nb = (nb_gene*inf_nb_V1R)/normalized_nb_genes) %>%
        dplyr::select(Subfamily, normalized_nb)
    )
  
  
  
  Subfamily_nb <- 
    Subfamily_nb %>%
    mutate(ID = curr_ID)
  
  V1R_DoC_df <- 
    rbind(V1R_DoC_df,
          Subfamily_nb)
  
}


rm(raw_coverage_V1R)



#Add the number of Total V1R genes
for(curr_ID in all_ID){
  total <- V1R_DoC_df %>% filter(ID == curr_ID) %>% pull(normalized_nb) %>% sum()
  temp_df <- as.data.frame("Total_V1R")
  temp_df <- cbind(temp_df, as.data.frame(total), as.data.frame(curr_ID))
  colnames(temp_df) <- c("Subfamily", "normalized_nb", "ID")
  
  V1R_DoC_df <- rbind(V1R_DoC_df, temp_df)
}





#### Compute number of genes per species from depth of coverage of V1R genes  ---------------------------------

V1R_DoC_df_info <- left_join(V1R_DoC_df, db_infos, by="ID")

V1R_DoC_df_info_mean <- 
  as.data.frame(
    V1R_DoC_df_info %>%
      group_by(Species, Subfamily) %>%
      mutate(mean_normalized_nb = mean(normalized_nb))
  ) %>%
  dplyr::select(Subfamily, mean_normalized_nb, Sp, Species, Tribe) %>%
  distinct()



write.table(V1R_DoC_df_info, "V1R_DoC_df_info.csv", quote=FALSE,
            col.names=TRUE, row.names=FALSE, sep=",")

write.table(V1R_DoC_df_info_mean, "V1R_DoC_df_info_mean.csv", quote=FALSE,
            col.names=TRUE, row.names=FALSE, sep=",")



#### Compute number of genes per ID from depth of coverage of T1R genes  ---------------------------------

raw_coverage_T1R <- 
  read.table("ALL.T1R.unmasked.raw.depth")
colnames(raw_coverage_T1R) <- c("ID", "scaffold", "pos", "coverage")
mean_coverage_T1R <- 
  read.table("ALL.T1R.unmasked.mean.depth")
colnames(mean_coverage_T1R) <- c("Gene_type", "Gene_name", "Clade", "Subfamily", 
                                 "OGG", "scaffold", "start", "end", "ID", "coverage")
coverage_BUSCO <- 
  read.table("ALL.AllBUSCO.unmasked.depth")
colnames(coverage_BUSCO) <- c("ID", "coverage")


#all_ID <- scan("/scicore/home/salzburg/polica0000/Reads_Depth_Of_Coverage/good_ID_list.txt", what="character")
all_ID <- scan("good_ID_list.txt", what="character")

T1R_DoC_df <- as.data.frame(NULL)
for(curr_ID in all_ID){
  mean_T1R_cov <- raw_coverage_T1R %>% filter(ID == curr_ID) %>% pull(coverage) %>% mean()
  mean_BUSCO_cov <- coverage_BUSCO %>% filter(ID == curr_ID) %>% pull(coverage)
  inf_nb_T1R <-  mean_T1R_cov/mean_BUSCO_cov * 13
  
  
  normalized_nb_genes <- 
    mean_coverage_T1R %>%
    filter(ID == curr_ID) %>%
    group_by(Gene_name) %>%
    mutate(cov_per_gene = mean(coverage)) %>%
    dplyr::select(- c(coverage,scaffold,start,end)) %>%
    distinct() %>%
    ungroup() %>%
    group_by(Subfamily) %>%
    mutate(cov_per_subfamily = mean(cov_per_gene)) %>%
    summarise(nb_onil_gene = n(), cov=cov_per_subfamily) %>%
    mutate(nb_gene = cov/mean_BUSCO_cov*nb_onil_gene) %>%
    distinct() %>%
    pull(nb_gene) %>% sum()
  
  Subfamily_nb <- 
    as.data.frame(
      mean_coverage_T1R %>%
        filter(ID == curr_ID) %>%
        group_by(Gene_name) %>%
        mutate(cov_per_gene = mean(coverage)) %>%
        dplyr::select(- c(coverage,scaffold,start,end)) %>%
        distinct() %>%
        ungroup() %>%
        group_by(Subfamily) %>%
        mutate(cov_per_subfamily = mean(cov_per_gene)) %>%
        summarise(nb_onil_gene = n(), cov=cov_per_subfamily) %>%
        mutate(nb_gene = cov/mean_BUSCO_cov*nb_onil_gene) %>%
        distinct() %>%
        mutate(normalized_nb = (nb_gene*inf_nb_T1R)/normalized_nb_genes) %>%
        dplyr::select(Subfamily, normalized_nb)
    )
  
  
  
  Subfamily_nb <- 
    Subfamily_nb %>%
    mutate(ID = curr_ID)
  
  T1R_DoC_df <- 
    rbind(T1R_DoC_df,
          Subfamily_nb)
  
}


rm(raw_coverage_T1R)



#Add the number of Total T1R genes
for(curr_ID in all_ID){
  total <- T1R_DoC_df %>% filter(ID == curr_ID) %>% pull(normalized_nb) %>% sum()
  temp_df <- as.data.frame("Total_T1R")
  temp_df <- cbind(temp_df, as.data.frame(total), as.data.frame(curr_ID))
  colnames(temp_df) <- c("Subfamily", "normalized_nb", "ID")
  
  T1R_DoC_df <- rbind(T1R_DoC_df, temp_df)
}





#### Compute number of genes per species from depth of coverage of T1R genes  ---------------------------------

T1R_DoC_df_info <- left_join(T1R_DoC_df, db_infos, by="ID")

T1R_DoC_df_info_mean <- 
  as.data.frame(
    T1R_DoC_df_info %>%
      group_by(Species, Subfamily) %>%
      mutate(mean_normalized_nb = mean(normalized_nb))
  ) %>%
  dplyr::select(Subfamily, mean_normalized_nb, Sp, Species, Tribe) %>%
  distinct()



write.table(T1R_DoC_df_info, "T1R_DoC_df_info.csv", quote=FALSE,
            col.names=TRUE, row.names=FALSE, sep=",")

write.table(T1R_DoC_df_info_mean, "T1R_DoC_df_info_mean.csv", quote=FALSE,
            col.names=TRUE, row.names=FALSE, sep=",")



#### Compute number of genes per ID from depth of coverage of T2R genes  ---------------------------------

raw_coverage_T2R <- 
  read.table("ALL.T2R.unmasked.raw.depth")
colnames(raw_coverage_T2R) <- c("ID", "scaffold", "pos", "coverage")
mean_coverage_T2R <- 
  read.table("ALL.T2R.unmasked.mean.depth")
colnames(mean_coverage_T2R) <- c("Gene_type", "Gene_name", "Clade", "Subfamily", 
                                 "OGG", "scaffold", "start", "end", "ID", "coverage")
coverage_BUSCO <- 
  read.table("ALL.AllBUSCO.unmasked.depth")
colnames(coverage_BUSCO) <- c("ID", "coverage")


#all_ID <- scan("/scicore/home/salzburg/polica0000/Reads_Depth_Of_Coverage/good_ID_list.txt", what="character")
all_ID <- scan("good_ID_list.txt", what="character")

T2R_DoC_df <- as.data.frame(NULL)
for(curr_ID in all_ID){
  mean_T2R_cov <- raw_coverage_T2R %>% filter(ID == curr_ID) %>% pull(coverage) %>% mean()
  mean_BUSCO_cov <- coverage_BUSCO %>% filter(ID == curr_ID) %>% pull(coverage)
  inf_nb_T2R <-  mean_T2R_cov/mean_BUSCO_cov * 1
  
  
  normalized_nb_genes <- 
    mean_coverage_T2R %>%
    filter(ID == curr_ID) %>%
    group_by(Gene_name) %>%
    mutate(cov_per_gene = mean(coverage)) %>%
    dplyr::select(- c(coverage,scaffold,start,end)) %>%
    distinct() %>%
    ungroup() %>%
    group_by(Subfamily) %>%
    mutate(cov_per_subfamily = mean(cov_per_gene)) %>%
    summarise(nb_onil_gene = n(), cov=cov_per_subfamily) %>%
    mutate(nb_gene = cov/mean_BUSCO_cov*nb_onil_gene) %>%
    distinct() %>%
    pull(nb_gene) %>% sum()
  
  Subfamily_nb <- 
    as.data.frame(
      mean_coverage_T2R %>%
        filter(ID == curr_ID) %>%
        group_by(Gene_name) %>%
        mutate(cov_per_gene = mean(coverage)) %>%
        dplyr::select(- c(coverage,scaffold,start,end)) %>%
        distinct() %>%
        ungroup() %>%
        group_by(Subfamily) %>%
        mutate(cov_per_subfamily = mean(cov_per_gene)) %>%
        summarise(nb_onil_gene = n(), cov=cov_per_subfamily) %>%
        mutate(nb_gene = cov/mean_BUSCO_cov*nb_onil_gene) %>%
        distinct() %>%
        mutate(normalized_nb = (nb_gene*inf_nb_T2R)/normalized_nb_genes) %>%
        dplyr::select(Subfamily, normalized_nb)
    )
  
  
  
  Subfamily_nb <- 
    Subfamily_nb %>%
    mutate(ID = curr_ID)
  
  T2R_DoC_df <- 
    rbind(T2R_DoC_df,
          Subfamily_nb)
  
}


rm(raw_coverage_T2R)



#Add the number of Total T2R genes
for(curr_ID in all_ID){
  total <- T2R_DoC_df %>% filter(ID == curr_ID) %>% pull(normalized_nb) %>% sum()
  temp_df <- as.data.frame("Total_T2R")
  temp_df <- cbind(temp_df, as.data.frame(total), as.data.frame(curr_ID))
  colnames(temp_df) <- c("Subfamily", "normalized_nb", "ID")
  
  T2R_DoC_df <- rbind(T2R_DoC_df, temp_df)
}





#### Compute number of genes per species from depth of coverage of T2R genes  ---------------------------------

T2R_DoC_df_info <- left_join(T2R_DoC_df, db_infos, by="ID")

T2R_DoC_df_info_mean <- 
  as.data.frame(
    T2R_DoC_df_info %>%
      group_by(Species, Subfamily) %>%
      mutate(mean_normalized_nb = mean(normalized_nb))
  ) %>%
  dplyr::select(Subfamily, mean_normalized_nb, Sp, Species, Tribe) %>%
  distinct()



write.table(T2R_DoC_df_info, "T2R_DoC_df_info.csv", quote=FALSE,
            col.names=TRUE, row.names=FALSE, sep=",")

write.table(T2R_DoC_df_info_mean, "T2R_DoC_df_info_mean.csv", quote=FALSE,
            col.names=TRUE, row.names=FALSE, sep=",")




