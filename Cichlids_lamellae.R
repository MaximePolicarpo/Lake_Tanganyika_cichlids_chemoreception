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
library("ggstar")
library(MASS)
library(ggeffects)
library("viridis")

#### Data load - Genomic -- Cichlids  ---------------------------------


#load phylo tree
b1_tree <- read.tree("b1_tree_filt.nwk")
b1_tree_wo_Neospl <- 
  drop.tip(b1_tree, "Neospl")
b1_list_sp <- b1_tree_wo_Neospl$tip.label

b1_tree_wo_Neospl_NodeLabel <- makeNodeLabel(b1_tree_wo_Neospl, method="number", prefix="Node")

#load pacbio mining results
t_pacbio_olfactory_receptor <- read.table("pacbio_olfactory_receptor_nb.tsv", header=FALSE, sep=",")
colnames(t_pacbio_olfactory_receptor) <- c("Species","OR_Complete","OR_Pseudogene","OR_Truncated","OR_Edge",
                                           "OR_Ambiguous","TAAR_Complete","TAAR_Pseudogene","TAAR_Truncated",
                                           "TAAR_Edge","TAAR_Ambiguous","V1R_Complete","V1R_Pseudogene",
                                           "V1R_Truncated","V1R_Edge","V1R_Ambiguous","V2R_Complete","V2R_Pseudogene",
                                           "V2R_Truncated","V2R_Edge","V2R_Ambiguous","T1R_Complete","T1R_Pseudogene",
                                           "T1R_Truncated","T1R_Edge","T1R_Ambiguous","T2R_Complete","T2R_Pseudogene",
                                           "T2R_Truncated","T2R_Edge","T2R_Ambiguous")

#load depth of coverage analysis results

OR_DoC_df_info <- read.table("OR_DoC_df_info.csv", header=TRUE,  sep=",")
OR_DoC_df_info <- OR_DoC_df_info %>% filter(! ID %in% c("LBE2", "LBE9"))
OR_DoC_df_info_mean <- read.table("OR_DoC_df_info_mean.csv", header=TRUE, sep=",")
OR_DoC_df_info_mean_total <- OR_DoC_df_info_mean %>% filter(Subfamily == "Total_OR") 

TAAR_DoC_df_info <- read.table("TAAR_DoC_df_info.csv", header=TRUE,  sep=",")
TAAR_DoC_df_info <- TAAR_DoC_df_info %>% filter(! ID %in% c("LBE2", "LBE9"))
TAAR_DoC_df_info_mean <- read.table("TAAR_DoC_df_info_mean.csv", header=TRUE, sep=",")
TAAR_DoC_df_info_mean_total <- TAAR_DoC_df_info_mean %>% filter(Subfamily == "Total_TAAR") 

V1R_DoC_df_info <- read.table("V1R_DoC_df_info.csv", header=TRUE,  sep=",")
V1R_DoC_df_info <- V1R_DoC_df_info %>% filter(! ID %in% c("LBE2", "LBE9"))
V1R_DoC_df_info_mean <- read.table("V1R_DoC_df_info_mean.csv", header=TRUE, sep=",")
V1R_DoC_df_info_mean_total <- V1R_DoC_df_info_mean %>% filter(Subfamily == "Total_V1R") 

V2R_DoC_df_info <- read.table("V2R_DoC_df_info.csv", header=TRUE,  sep=",")
V2R_DoC_df_info <- V2R_DoC_df_info %>% filter(! ID %in% c("LBE2", "LBE9"))
V2R_DoC_df_info_mean <- read.table("V2R_DoC_df_info_mean.csv", header=TRUE, sep=",")
V2R_DoC_df_info_mean_total <- V2R_DoC_df_info_mean %>% filter(Subfamily == "Total_V2R") 

T1R_DoC_df_info <- read.table("T1R_DoC_df_info.csv", header=TRUE,  sep=",")
T1R_DoC_df_info <- T1R_DoC_df_info %>% filter(! ID %in% c("LBE2", "LBE9"))
T1R_DoC_df_info_mean <- read.table("T1R_DoC_df_info_mean.csv", header=TRUE, sep=",")
T1R_DoC_df_info_mean_total <- T1R_DoC_df_info_mean %>% filter(Subfamily == "Total_T1R") 

T2R_DoC_df_info <- read.table("T2R_DoC_df_info.csv", header=TRUE,  sep=",")
T2R_DoC_df_info <- T2R_DoC_df_info %>% filter(! ID %in% c("LBE2", "LBE9"))
T2R_DoC_df_info_mean <- read.table("T2R_DoC_df_info_mean.csv", header=TRUE, sep=",")
T2R_DoC_df_info_mean_total <- T2R_DoC_df_info_mean %>% filter(Subfamily == "Total_T2R") 

#load db infos on cichlids
db_infos <- read.table("DNATube_2021-03-04_16-15-42_Simplified_parsed.tsv", sep="\t", header=TRUE)
colnames(db_infos) <- c("Sp", "Species", "ID", "Sex", "Tribe")
#db_infos <- db_infos %>% filter(! Sp == "Neospl") #remove Neosplen

#load Tribe info from b1 tree
tribes_infos <- read.table("Tribes_Sp.tsv", sep="\t", header=FALSE)
colnames(tribes_infos) <- c("Tribe_n", "Sp")
db_infos <- left_join(db_infos, tribes_infos, by="Sp")
db_infos <- 
  db_infos %>% 
  mutate(final_Tribe = if_else(
    is.na(Tribe_n),
    Tribe,
    Tribe_n
  )) 

db_infos <- 
  db_infos %>%
  mutate(Final_Tribe = if_else(
    final_Tribe %in% c("Bathybatini","Benthochromini",
                       "Boulengerochromini","Cyphotilapiini",
                       "Cyprichromini","Ectodini","Eretmodini",
                       "Lamprologini","Limnochromini","Perissodini",
                       "Trematocarini","Tropheini","Haplochromini","outgroup"),
    final_Tribe,
    "Other"
  ))


#Define gene family colors 

chemoreceptor_family_colors <- 
  c(Total_OR="lightcoral",
    Total_TAAR = "#377EB8", 
    Total_V1R = "#FFFF33",
    Total_V2R = "#FF7F00",
    Total_T1R = "#4DAF4A", 
    Total_T2R = "#984EA3")


#Define the tribe colors (same as Ronco et al. 2021)
tribes.colors <- c(Bathybatini = "#242626" , Benthochromini = "#AE262A" , Boulengerochromini = "#59595C" , Cyphotilapiini = "#FDDF13" , Cyprichromini = "#F04D29" , Ectodini = "#9AB9D9" , Eretmodini = "#682E7A" , Lamprologini = "#C588BB" , Limnochromini = "#535CA9" , Perissodini = "#FBAC43" , Trematocarini = "#959170" , Tropheini = "#86C773" , Haplochromini = "#274e13" , outgroup = "gray", Other="gray")
sp.colors <- c(Batmin = "#242626" , Cypfro = "#FDDF13" , Cyplep = "#F04D29" , Cunlon = "#9AB9D9" , Neomul = "#C588BB" , Simdia = "#86C773")


#Importe table with ecological variables per species 
ecoCat <- read.table("Species_ecoCat.csv", sep=",")
colnames(ecoCat) <- c("Sp", "Tribe", "Food", "habitat")
ecoCat <- ecoCat %>% filter(! Sp == "Neospl")



#SI table
SI_infos <- read.table("Stable_isotope_table_simplified.tsv", sep="\t", header=TRUE)
colnames(SI_infos) <- c("ID", "Nmg", "Cmg", "d15N", "d13C", "Sp")
SI_infos <- SI_infos %>% filter(! Sp == "Neospl")


SI_infos_mean <- as.data.frame(SI_infos %>% 
                                 group_by(Sp) %>% 
                                 mutate(mean_d15N = mean(d15N),
                                        mean_d13C = mean(d13C)) %>%
                                 ungroup() %>%
                                 dplyr::select(mean_d15N, mean_d13C, Sp) %>%
                                 distinct())

#Breeding table
breed_infos <- read.table("Species_tribe_breeding.tsv", sep="\t", header=TRUE)
colnames(breed_infos) <- c("Sp", "Species", "Tribe", "breeding_type", "breeding_mode")
breed_infos <- breed_infos %>% filter(! Sp == "Neospl")

#OR table detailed
OR_illumina <- read.table("All_illumina_OR_results.tsv", header=FALSE, sep="\t")
colnames(OR_illumina) <- c("ID", "Complete", "Pseudogene", "Truncated", "Edge")

#TAAR table detailed
TAAR_illumina <- read.table("All_illumina_TAAR_results.tsv", header=FALSE, sep="\t")
colnames(TAAR_illumina) <- c("ID", "Complete", "Pseudogene", "Truncated", "Edge")

#V2R table detailed
V2R_illumina <- read.table("All_illumina_V2R_results.tsv", header=FALSE, sep="\t")
colnames(V2R_illumina) <- c("ID", "Complete", "Pseudogene", "Truncated", "Edge")

#V1R table detailed
V1R_illumina <- read.table("All_illumina_V1R_results.tsv", header=FALSE, sep="\t")
colnames(V1R_illumina) <- c("ID", "Complete", "Pseudogene", "Truncated", "Edge")

#T2R table detailed
T2R_illumina <- read.table("All_illumina_T2R_results.tsv", header=FALSE, sep="\t")
colnames(T2R_illumina) <- c("ID", "Complete", "Pseudogene", "Truncated", "Edge")

#T1R table detailed
T1R_illumina <- read.table("All_illumina_T1R_results.tsv", header=FALSE, sep="\t")
colnames(T1R_illumina) <- c("ID", "Complete", "Pseudogene", "Truncated", "Edge")

#Radiation species

radiation_tree_wo_Neospl <- drop.tip(b1_tree_wo_Neospl, c("Gobeth", "Hetbut", "Tilbre", "Steult", "Tilspa"))


#parse the tables

OLR_Doc_mean_sp <- 
  rbind(OR_DoC_df_info_mean_total, TAAR_DoC_df_info_mean_total, 
        V1R_DoC_df_info_mean_total, V2R_DoC_df_info_mean_total)

TR_Doc_mean_sp <- 
  rbind(T1R_DoC_df_info_mean_total, T2R_DoC_df_info_mean_total)

All_Doc_mean_sp <- 
  rbind(OR_DoC_df_info_mean_total, TAAR_DoC_df_info_mean_total, 
        V1R_DoC_df_info_mean_total, V2R_DoC_df_info_mean_total,
        T1R_DoC_df_info_mean_total, T2R_DoC_df_info_mean_total)

Tribe_for_tree <- 
  OLR_Doc_mean_sp %>%
  dplyr::select(Sp, Tribe)


All_Doc_mean_sp_wide <- 
  as.data.frame(
    pivot_wider(All_Doc_mean_sp, 
                names_from = c(Subfamily), values_from = mean_normalized_nb)
  )



#### Data load - Genomic -- Rayfinned fishes  ---------------------------------

#data from Policarpo et al. 2024
chemoreceptor_tidy_df <- read.table("Table_count_chemoreceptors.csv", sep=",", header=FALSE)
colnames(chemoreceptor_tidy_df) <- c("species", "gene_clade", "gene_type", "gene_family","number")

chemoreceptor_tidy_df_total <- 
  chemoreceptor_tidy_df %>%
  filter(gene_clade == "Total")

chemoreceptor_tidy_df_total_resh <-
  chemoreceptor_tidy_df_total %>%
  mutate(type_family = paste(gene_type, gene_family, sep="_")) %>%
  dplyr::select(species,type_family,number)

chemoreceptor_tidy_df_total_resh <-
  reshape(chemoreceptor_tidy_df_total_resh,
          idvar = "species",
          timevar   = "type_family",
          direction = "wide")


species_taxa <- read.table("species_list.tsv", header=TRUE, sep="\t")

busco_results <- read.table("BUSCO_results_table.tsv", header=TRUE, sep="\t")

species_taxa <- left_join(species_taxa, busco_results, by="genome_name")

species_taxa <-
  species_taxa %>%
  mutate(busco_perc_complete = Complete/(Complete + Fragmented + Missing)) %>%
  filter(busco_perc_complete > 0.8)


sp_tips_class <-
  species_taxa %>%
  dplyr::select(species, class)

chemoreceptor_tidy_df_total_resh_spclass <-
  left_join(sp_tips_class, 
            chemoreceptor_tidy_df_total_resh, 
            by="species")

chemoreceptor_table_RayFinnedFishes <- 
  chemoreceptor_tidy_df_total_resh_spclass %>%
  filter(class == "Actinopterygii")


sp_tips_order <-
  species_taxa %>%
  dplyr::select(species, order_ncbi)


lamellae_actino_nb <- read.table("Actino_lamellae_nb.tsv", header=TRUE, sep="\t")
colnames(lamellae_actino_nb) <- c("Accepted_species_name", "Unaccepted_species_name",
                                  "species", "mean_lamellae_number", "Binary_coding")
lamellae_actino_nb <- 
  lamellae_actino_nb %>%
  filter(mean_lamellae_number != "NA")

lamellae_actino_chemoreceptor_nb <- 
  left_join(chemoreceptor_table_RayFinnedFishes, lamellae_actino_nb, by='species') %>%
  filter(! is.na(class)) %>%
  filter(class == "Actinopterygii")

lamellae_actino_chemoreceptor_nb_order <- 
  left_join(lamellae_actino_chemoreceptor_nb,sp_tips_order, by="species")


lamellae_actino_chemoreceptor_nb_order <-
  lamellae_actino_chemoreceptor_nb_order %>%
  mutate(total_olfactory_receptor = 
           number.Complete_OR + number.Complete_TAAR + 
           number.Complete_V2R + number.Complete_V1R)






#### Data load -- Lamellae - Rayfinned fishes ---------------------------------

#Import the number of lamellae per species for ray-finned fishes
rayfinned_lamellae <- 
  read.table("Ray_finned_fishes_lamellae.tsv", sep="\t",
             header=TRUE)
rayfinned_lamellae_nb <- 
  rayfinned_lamellae %>%
  filter(! is.na(Lamellae_number_mean))


#Import the  SL of ray-finned fishes and add to the lamellae table

Lengths_df_rayfinned <- 
  read.table("Teleost_length_table.tsv", 
             sep="\t",
             header=TRUE)

Lengths_df_rayfinned <- 
  as.data.frame(Lengths_df_rayfinned %>%
  rowwise() %>%
   mutate(fishtree_name = gsub(" ", "_", Species)))
  
#Import a ray-finned fish phylogeny

rayfinned_phylo <- 
  ape::read.tree("actinopt_12k_treePL.tre")

rayfinned_phylo$tip.label[rayfinned_phylo$tip.label=="Neogobius_fluviatilis"] <- "Neogobius_melanostomus"
rayfinned_phylo$tip.label[rayfinned_phylo$tip.label=="Omobranchus_obliquus"] <- "Omobranchus_elegans"
rayfinned_phylo$tip.label[rayfinned_phylo$tip.label=="Triplophysa_rosa"] <- "Triplophysa_dalaica"
rayfinned_phylo$tip.label[rayfinned_phylo$tip.label=="Lutjanus_ophuysenii"] <- "Hoplopagrus_guentherii"
rayfinned_phylo$tip.label[rayfinned_phylo$tip.label=="Scoloplax_distolothrix"] <- "Scoloplax_empousa"
rayfinned_phylo$tip.label[rayfinned_phylo$tip.label=="Sardinops_sagax"] <- "Sardinops_melanosticta"
rayfinned_phylo$tip.label[rayfinned_phylo$tip.label=="Sarcocheilichthys_soldatovi"] <- "Sarcocheilichthys_variegatus"
rayfinned_phylo$tip.label[rayfinned_phylo$tip.label=="Pagrus_pagrus"] <- "Pagrus_major"
rayfinned_phylo$tip.label[rayfinned_phylo$tip.label=="Halichoeres_leucurus"] <- "Parajulis_poecilepterus"
rayfinned_phylo$tip.label[rayfinned_phylo$tip.label=="Minous_trachycephalus"] <- "Minous_monodactylus"
rayfinned_phylo$tip.label[rayfinned_phylo$tip.label=="Lotella_rhacina"] <- "Lotella_phycis"
rayfinned_phylo$tip.label[rayfinned_phylo$tip.label=="Fistularia_commersonii"] <- "Fistularia_tabacaria"
rayfinned_phylo$tip.label[rayfinned_phylo$tip.label=="Craterocephalus_nouhuysi"] <- "Craterocephalus_sp"
rayfinned_phylo$tip.label[rayfinned_phylo$tip.label=="Pseudotropheus_elongatus"] <- "Pseudotropheus_fuellburnii"
rayfinned_phylo$tip.label[rayfinned_phylo$tip.label=="Benthodesmus_simonyi"] <- "Benthodesmus_tenuis"


#http://www.phytools.org/eqg2015/asr.html

#First prune tree for species present in the table

rayfinned_phylo_splist <- 
  rayfinned_phylo$tip.label

rayfinned_lamellae_nb_inphylo <- as.data.frame(NULL)
for(species in rayfinned_phylo_splist){
  
  current_sp_df <- 
    head(
      rayfinned_lamellae_nb %>%
        filter(., Good_species_name == species | Species_name_in_literature == species | Species_name_in_NCBI == species), 
      1) %>%
    mutate(fishtree_name = species)
  
  rayfinned_lamellae_nb_inphylo <- 
    rbind(rayfinned_lamellae_nb_inphylo, current_sp_df)
  
}

#add the length to the table
rayfinned_lamellae_nb_inphylo <- 
  left_join(rayfinned_lamellae_nb_inphylo, 
            Lengths_df_rayfinned, 
            by="fishtree_name")



#Now filter the tree

rayfinned_lamellae_nb_inphylo_filt <- 
  rayfinned_lamellae_nb_inphylo %>%
  dplyr::select(Lamellae_number_mean, Binary_coding, fishtree_name) %>%
  distinct()

rayfinned_phylo_pruned <- 
  keep.tip(rayfinned_phylo, 
           rayfinned_lamellae_nb_inphylo_filt %>% pull(fishtree_name))




#Now prune only teleost species


teleost_lamellae_nb_inphylo_filt <- 
  dplyr::filter(rayfinned_lamellae_nb_inphylo, !grepl("Acipenser",fishtree_name))
teleost_lamellae_nb_inphylo_filt <- 
  dplyr::filter(teleost_lamellae_nb_inphylo_filt, !grepl("Polypterus",fishtree_name))
teleost_lamellae_nb_inphylo_filt <- 
  dplyr::filter(teleost_lamellae_nb_inphylo_filt, !grepl("calabaricus",fishtree_name))
teleost_lamellae_nb_inphylo_filt <- 
  dplyr::filter(teleost_lamellae_nb_inphylo_filt, !grepl("Amia_calva",fishtree_name))
teleost_lamellae_nb_inphylo_filt <- 
  dplyr::filter(teleost_lamellae_nb_inphylo_filt, !grepl("Lepisosteus",fishtree_name))
teleost_lamellae_nb_inphylo_filt <- 
  dplyr::filter(teleost_lamellae_nb_inphylo_filt, !grepl("Scaphirhynchus",fishtree_name))
teleost_lamellae_nb_inphylo_filt <- 
  dplyr::filter(teleost_lamellae_nb_inphylo_filt, !grepl("Polyodon",fishtree_name))

teleost_phylo_pruned <- 
  keep.tip(rayfinned_phylo, 
           teleost_lamellae_nb_inphylo_filt %>% pull(fishtree_name))




#### Data load - Lamellae -- Cichlids  ---------------------------------

#Import eye table

Eye_size_df <- 
  read.table("Eye_size_cichlids.csv",
             sep=" ",
             header=TRUE)
colnames(Eye_size_df) <- c("ID", "Sp", "Csize", "A_eye", "SL", "body_depth")

Eye_size_df_sp <- 
  as.data.frame(
    Eye_size_df %>%
      group_by(Sp) %>%
      summarise(mean_Csize = mean(Csize),
                mean_Aeye = mean(A_eye),
                relative_Aeye = sqrt(mean_Aeye)/mean_Csize))


#Import lamellae count and info tables 

#Lamellae_count_df <- 
#  read.table("Lamellae_cichlids_count.tsv", sep="\t", header=TRUE)
#colnames(Lamellae_count_df) <- c("ID", "R_lamellae", "L_lamellae")
#

Lamellae_info_df <- 
  read.table("Lamellae_cichlids_info.tsv", sep="\t", header=TRUE)

Lamellae_count_info_df <- 
  read.table("Lamellae_cichlids_count_info.tsv", sep="\t", header=TRUE)


#Remove specimen with no lamellae count 

Lamellae_count_info_df <- 
  Lamellae_count_info_df %>%
  filter(! (is.na(R_lamellae) & is.na(L_lamellae))) 

#Compute the mean number of lamellae per specimen

Lamellae_count_info_df <-
  as.data.frame(Lamellae_count_info_df %>%
  rowwise() %>%
  mutate(mean_lamellae_ind = mean(c(R_lamellae, L_lamellae), na.rm=TRUE)))



#Number of total individuals for which we have at-least one rosette count

nrow(Lamellae_count_info_df)

#Generate a table with species and the number of ind with one or two rosettes counted
Total_ind_per_sp <- 
  as.data.frame(
    Lamellae_count_info_df %>%
      group_by(Sp) %>%
      summarise(count_total = n()))

#Number of species part of the LT radiation
nrow(Total_ind_per_sp %>%
       filter(Sp %in% radiation_tree_wo_Neospl$tip.label))

#Number of outgroup species
nrow(Total_ind_per_sp %>%
       filter(! Sp %in% radiation_tree_wo_Neospl$tip.label))



write.table(Lamellae_count_info_df,
            "Table_detailed_lamellae_pub.tsv", 
            sep="\t", col.names=TRUE, row.names=FALSE,
            quote=FALSE)
#Generate a table with species and the numb of ind with the two rosettes counted

Total_L_and_R <- 
  as.data.frame(
    Lamellae_count_info_df %>%
      filter(! is.na(R_lamellae)) %>%
      filter(! is.na(L_lamellae)) %>%
      group_by(Sp) %>%
      summarise(count = n()))


Table_summary <- left_join(Total_ind_per_sp, Total_L_and_R, by="Sp")
colnames(Table_summary) <- c("Sp", "Nb_ind", "Nb_ind_both_sides")
Table_summary[is.na(Table_summary)] <- 0


Table_summary <- 
  Table_summary %>%
  mutate(Nb_ind_one_side = Nb_ind - Nb_ind_both_sides)



write.table(Table_summary, "Table_summary_lamellae.tsv", 
            sep="\t", col.names=TRUE, row.names=FALSE,
            quote=FALSE)


# Now compute the mean number of lamellae per species

Lamellae_count_per_sp <- 
  as.data.frame(
    Lamellae_count_info_df %>%
      group_by(Sp) %>%
      summarise(mean_lamellae = mean(mean_lamellae_ind, na.rm = TRUE),
                mean_SL = mean(SL, na.rm = TRUE),
                mean_TL = mean(TL, na.rm = TRUE),
                mean_weight = mean(Weight, na.rm = TRUE)) %>%
      ungroup() %>%
      dplyr::select(Sp, mean_lamellae, mean_SL, mean_TL, mean_weight) %>%
      filter(! is.na(mean_lamellae)) %>%
      filter(! is.na(Sp)) %>%
      distinct()
  )

#Add the mean stats for Plestr taken for labkey
Lamellae_count_per_sp[(Lamellae_count_per_sp$Sp == "Plestr"),"mean_SL"] <- 7.806
Lamellae_count_per_sp[(Lamellae_count_per_sp$Sp == "Plestr"),"mean_TL"] <- 9.645
Lamellae_count_per_sp[(Lamellae_count_per_sp$Sp == "Plestr"),"mean_weight"] <- 13.455


write.table(Lamellae_count_per_sp,
            file="~/Olfactory_epithelium_lamellae/Cichlid_rosette_table.tsv",
            sep="\t",
            row.names=FALSE,
            col.names=FALSE,
            quote=FALSE)

#### Cichlids lamellae - Differences between left and right sides ?  ---------------------------------

Lamellae_count_info_df %>%
  filter(! is.na(R_lamellae)) %>%
  filter(! is.na(L_lamellae)) %>%
  mutate(diff_L_R = abs(R_lamellae - L_lamellae)) %>%
  pull(diff_L_R) %>%
  mean()

Lamellae_count_info_df %>%
  filter(! is.na(R_lamellae)) %>%
  filter(! is.na(L_lamellae)) %>%
  mutate(diff_L_R = abs(R_lamellae - L_lamellae)) %>%
  pull(diff_L_R) %>%
  max()


Lamellae_count_info_df %>%
  filter(! is.na(R_lamellae)) %>%
  filter(! is.na(L_lamellae)) %>%
  mutate(diff_L_R = abs(R_lamellae - L_lamellae)) %>%
  pull(diff_L_R) %>%
  min()


#### Cichlids lamellae - Differences between individuals of the same species ?  ---------------------------------

Lamellae_count_info_df <- as.data.frame(Lamellae_count_info_df)


list_sp <- 
  Lamellae_count_info_df %>%
  group_by(Sp) %>%
  summarise(count = n()) %>%
  filter(count >= 2) %>%
  pull(Sp)

vector_mean_diff <- c()
for(curr_sp in list_sp){
  
  current_counts <- 
    Lamellae_count_info_df %>%
    filter(Sp == curr_sp) %>%
    pull(mean_lamellae_ind) 
  
  vector_mean_diff <- c(vector_mean_diff, mean(abs(combn(current_counts,2, FUN = base::diff))))
  
}

mean(vector_mean_diff)


pdf(file = "Raw_R_plots/Lamellae_per_Ind.pdf",width = 8.34,  height = 4.61)

Lamellae_count_info_df %>%
  ggplot(., aes(x=reorder(Sp,mean_lamellae_ind) , y=mean_lamellae_ind, color=Tribe)) +
  geom_boxplot() +
  xlab("Species") +
  ylab("Mean number of lamellae") +
  theme_minimal() +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        axis.text.x=element_blank(),
        plot.subtitle=element_text(size=16),
        legend.position="none") +
  scale_color_manual(values = tribes.colors)

dev.off()




#### Cichlids lamellae - Differences between Sex ?  ---------------------------------

Lamellae_count_info_df_sexknown <- 
  Lamellae_count_info_df %>%
  filter(! is.na(Sex)) %>%
  filter(Sex != "") %>%
  filter(Sex != "M?") 

Lamellae_count_info_df_sexknown_summary <- 
  as.data.frame(
    Lamellae_count_info_df_sexknown %>%
      group_by(Sp, Sex) %>%
      summarise(count = n()))

sp_list <- Lamellae_count_info_df_sexknown_summary$Sp %>% unique()

t_test_stats_vector <- c()
t_test_pvalue_vector <- c()
t_test_species_vector <- c()
for(curr_sp in sp_list){
  nb_M <- 
    Lamellae_count_info_df_sexknown_summary %>%
    filter(Sp == curr_sp) %>%
    filter(Sex == 'M') %>% pull(count)
  
  
  nb_F <- 
    Lamellae_count_info_df_sexknown_summary %>%
    filter(Sp == curr_sp) %>%
    filter(Sex == 'F') %>% pull(count)
  
  if(length(nb_M) == 0){ nb_M = 0 }
  if(length(nb_F) == 0){ nb_F = 0 }
  
  if(nb_M >= 3 & nb_F >= 3){
    
    print(curr_sp)
    
    female_counts <- 
      Lamellae_count_info_df_sexknown %>%
      filter(Sp == curr_sp) %>%
      filter(Sex == 'F') %>% pull(mean_lamellae_ind)
    
    male_counts <- 
      Lamellae_count_info_df_sexknown %>%
      filter(Sp == curr_sp) %>%
      filter(Sex == 'M') %>% pull(mean_lamellae_ind)
    
    t_test_result <- t.test(female_counts, male_counts)
    
    t_test_stat <- t_test_result$statistic
    t_test_pvalue <- t_test_result$p.value
    
    t_test_species_vector <- c(t_test_species_vector, curr_sp)
    t_test_stats_vector <- c(t_test_stats_vector, t_test_stat)
    t_test_pvalue_vector <- c(t_test_pvalue_vector, t_test_pvalue)
  }
  
}


Sex_vs_lamellae_ttest_df <- 
  as.data.frame(cbind(t_test_species_vector, t_test_stats_vector, t_test_pvalue_vector))
colnames(Sex_vs_lamellae_ttest_df) <- c("species", "ttest_stat", "ttest_pvalue")

Sex_colors <- 
  c('F'="#D81B60",
    'M' = "#1E88E5")

#Lamellae_count_info_df_sexknown %>%
#  ggplot(., aes(x=Sp, y=mean_lamellae_ind,color=Sex)) +
#  geom_boxplot() +
#  scale_color_manual(values = Sex_colors) +
#  theme_classic()


pdf(file = "Raw_R_plots/Lamellae_vs_Sex.pdf",width = 8.34,  height = 4.61)

Lamellae_count_info_df_sexknown %>%
  filter(Sp %in% c("Benmel", "Chabri", "Ectdes", "Neochr")) %>%
  ggplot(., aes(x=Sp, y=mean_lamellae_ind,color=Sex)) +
  geom_boxplot() +
  scale_color_manual(values = Sex_colors) +
  theme_classic() +
  xlab("") +
  ylab("Mean number of lamellae") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none")

dev.off()

#### Cichlids lamellae - Scale with SL at the level of intra-species ?  ---------------------------------

head(Lamellae_count_info_df %>%
       arrange(desc(mean_lamellae_ind)),3)

Lamellae_count_info_df_known_SL <- 
  Lamellae_count_info_df %>%
  filter(! is.na(SL)) 

Lamellae_count_info_df_known_SL_summary <- 
  as.data.frame(
    Lamellae_count_info_df_known_SL %>%
      group_by(Sp) %>%
      summarise(count = n()))

sp_list_atleast4 <- 
  Lamellae_count_info_df_known_SL_summary %>%
  filter(count >= 4) %>%
  pull(Sp)


Intraspecies_corr_SL <- as.data.frame(NULL)
for(curr_sp in sp_list_atleast4){
  
  lamellae_count <- 
    Lamellae_count_info_df_known_SL %>%
    filter(Sp == curr_sp) %>%
    pull(mean_lamellae_ind)
  
  SL_measures <- 
    Lamellae_count_info_df_known_SL %>%
    filter(Sp == curr_sp) %>%
    pull(SL)
  
  
  corr_SL_lamellae <- cor.test(lamellae_count, SL_measures, method="pearson")
  
  corr_pvalue <- as.numeric(corr_SL_lamellae$p.value)
  corr_R <- as.numeric(corr_SL_lamellae$estimate)
  
  curr_df <- as.data.frame(cbind(curr_sp, corr_pvalue, corr_R))
  colnames(curr_df) <- c("Sp", "pvalue", "R")
  
  Intraspecies_corr_SL <- 
    rbind(Intraspecies_corr_SL,
          curr_df)
  
}


Intraspecies_corr_SL$pvalue <- as.numeric(Intraspecies_corr_SL$pvalue)
Intraspecies_corr_SL$R <- as.numeric(Intraspecies_corr_SL$R)

Intraspecies_corr_SL

Intraspecies_corr_SL[(is.na(Intraspecies_corr_SL$pvalue)),"pvalue"] <- 1
Intraspecies_corr_SL[(is.na(Intraspecies_corr_SL$R)),"R"] <- 0


Intraspecies_corr_SL %>%
  filter(pvalue < 0.05)

write.table(Intraspecies_corr_SL, "Raw_R_plots/Intraspecies_corr_SL.tsv", 
            sep="\t", col.names=TRUE, row.names=FALSE,
            quote=FALSE)



#### Cichlids lamellae - SL + Sex  ---------------------------------


Lamellae_count_info_df_sexknown_summary <- 
  as.data.frame(
    Lamellae_count_info_df_sexknown %>%
      group_by(Sp, Sex) %>%
      summarise(count = n()))

sp_list <- Lamellae_count_info_df_sexknown_summary$Sp %>% unique()

SL_sex_df <- as.data.frame(NULL)
for(curr_sp in sp_list){
  nb_M <- 
    Lamellae_count_info_df_sexknown_summary %>%
    filter(Sp == curr_sp) %>%
    filter(Sex == 'M') %>% pull(count)
  
  
  nb_F <- 
    Lamellae_count_info_df_sexknown_summary %>%
    filter(Sp == curr_sp) %>%
    filter(Sex == 'F') %>% pull(count)
  
  if(length(nb_M) == 0){ nb_M = 0 }
  if(length(nb_F) == 0){ nb_F = 0 }
  
  if(nb_M >= 3 & nb_F >= 3){
    
    print(curr_sp)
    
    female_counts <- 
      Lamellae_count_info_df_sexknown %>%
      filter(Sp == curr_sp) %>%
      filter(Sex == 'F') %>% pull(mean_lamellae_ind)
    
    male_counts <- 
      Lamellae_count_info_df_sexknown %>%
      filter(Sp == curr_sp) %>%
      filter(Sex == 'M') %>% pull(mean_lamellae_ind)
    
    
    
    lm_sex_sl <- 
      lm(mean_lamellae_ind ~ Sex + SL,
         data= Lamellae_count_info_df_sexknown %>% filter(Sp == curr_sp))
    summary_lm_sex_sl <- summary(lm_sex_sl)
    
    curr_R <- summary_lm_sex_sl$adj.r.squared
    curr_sex_pvalue <- summary_lm_sex_sl$coefficients[11]
    curr_SL_pvalue <- summary_lm_sex_sl$coefficients[12]

    curr_model_pvalue = pf(summary_lm_sex_sl$fstatistic[1], 
                           summary_lm_sex_sl$fstatistic[2], 
                           summary_lm_sex_sl$fstatistic[3], 
                                  lower.tail = FALSE)

    curr_df <- 
      as.data.frame(
        cbind(curr_sp, curr_sex_pvalue, curr_SL_pvalue, curr_model_pvalue, curr_R)
      )
    
    
    SL_sex_df <- rbind(SL_sex_df, curr_df)
    
  }
  
}



write.table(SL_sex_df, "Raw_R_plots/SL_sex_df.tsv", 
            sep="\t", col.names=TRUE, row.names=FALSE,
            quote=FALSE)






#### Cichlids lamellae -  Interesting species level stats ---------------------------------

Lamellae_count_per_sp %>%
  pull(mean_lamellae) %>%
  mean()

Lamellae_count_per_sp %>%
  pull(mean_lamellae) %>%
  min()


Lamellae_count_per_sp %>%
  pull(mean_lamellae) %>%
  max()

Lamellae_count_per_sp %>%
  arrange(desc(mean_lamellae))


#### Cichlids lamellae - Boxplot per tribe  ---------------------------------


Lamellae_count_per_sp_tribe <- 
  left_join(Lamellae_count_per_sp, Tribe_for_tree %>% distinct(), by="Sp")
Lamellae_count_per_sp_tribe$Tribe[Lamellae_count_per_sp_tribe$Tribe=="Oreochromini"] <- "Other"
Lamellae_count_per_sp_tribe$Tribe[Lamellae_count_per_sp_tribe$Tribe=="Tylochromini"] <- "Other"

pdf(file = "Raw_R_plots/BoxPlot_Lamellae_tribe.pdf",width = 6.34,  height = 4.61)

Lamellae_count_per_sp_tribe %>%
  ggplot(., aes(x=reorder(Tribe,mean_lamellae), y=mean_lamellae, fill=Tribe)) +
  geom_boxplot() + 
  scale_fill_manual(values = tribes.colors) +
  scale_color_manual(values = tribes.colors) +
  theme_classic() + 
  xlab("Tribe") +
  ylab("Mean number of lamellae") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.text.x=element_blank(),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position = "none")

dev.off()


as.data.frame(
  Lamellae_count_per_sp_tribe %>%
  group_by(Tribe) %>%
  summarise(mean_L = mean(mean_lamellae)) %>%
  arrange(mean_L)
)


#### Cichlids lamellae - Compare with teleost distributions  ---------------------------------

cichlid_teleost_color <- 
  c(Cichlid="#004D40",
    Teleost = "#FFC107")

teleost_lamellae_nb_inphylo_filt <- 
  teleost_lamellae_nb_inphylo_filt %>%
  mutate(hist_color = "Teleost")

Lamellae_count_per_sp <- 
  Lamellae_count_per_sp %>%
  mutate(hist_color = "Cichlid")


teleost_for_join_df <- 
  teleost_lamellae_nb_inphylo_filt %>%
  dplyr::select(Good_species_name, Lamellae_number_mean, hist_color)
colnames(teleost_for_join_df) <- c("Sp", "mean_lamellae", "hist_color")

Lamellae_cichlid_teleost <- 
  rbind(Lamellae_count_per_sp %>%
          dplyr::select(Sp, mean_lamellae, hist_color), 
        teleost_for_join_df)


pdf(file = "Raw_R_plots/Histo_TeleostCichlid_lamellae.pdf",width = 6.34,  height = 4.61)

Lamellae_cichlid_teleost %>%
  ggplot(., aes(x=mean_lamellae, fill=hist_color)) +
  geom_histogram(bins=60, alpha = 0.5, position = "identity", color="black") +
  theme_classic() +
  scale_fill_manual(values = cichlid_teleost_color) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position = "none") +
  xlab("Number of species") +
  ylab("Mean number of lamellae")

dev.off()



pdf(file = "Raw_R_plots/Histo_Cichlid_lamellae_zoom.pdf",width = 6.34,  height = 4.61)

Lamellae_count_per_sp %>%
  ggplot(., aes(x=mean_lamellae, fill=hist_color)) +
  geom_histogram(bins=20, alpha = 0.5, position = "identity", color="black") +
  theme_classic() +
  scale_fill_manual(values = cichlid_teleost_color) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position = "none") +
  xlab("Number of species") +
  ylab("Mean number of lamellae")

dev.off()



#Now draw a plot of the number of lamellae per tribe

Lamellae_count_per_sp_tribe <- 
  left_join(Lamellae_count_per_sp, tribes_infos, by="Sp")

pdf(file = "Raw_R_plots/Boxplot_Tribe_lamellae.pdf",width = 6.34,  height = 4.61)

Lamellae_count_per_sp_tribe %>%
  filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>%
  ggplot(., aes(x=reorder(Tribe_n, mean_lamellae) , y=mean_lamellae, color=Tribe_n)) +
  geom_boxplot() +
  theme_classic() +
  scale_color_manual(values = tribes.colors) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14), 
        axis.text.x=element_blank()) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position = "none") +
  xlab("Tribe") +
  ylab("Mean number of lamellae")

dev.off()



#### Cichlids lamellae - pGLS between lamellae and SL  ---------------------------------

#first lets test the correlation among cichlids
Tribe_for_tree <- distinct(Tribe_for_tree)
Lamellae_count_per_sp <- 
  left_join(Lamellae_count_per_sp, Tribe_for_tree, by="Sp") 

caper_data_cichlids_SL_lamellae <- 
  comparative.data(phy = b1_tree, 
                   data = Lamellae_count_per_sp,
                   names.col = Sp, vcv = TRUE,
                   na.omit = FALSE, warn.dropped = TRUE)



fit_phylo_lamellae_SL <- pgls(mean_lamellae ~ mean_SL,
                              data = caper_data_cichlids_SL_lamellae, 
                              lambda = "ML")
summary(fit_phylo_lamellae_SL)

sum_fit_phy <- summary(fit_phylo_lamellae_SL)
PGLS_pente =   formatC(sum_fit_phy$coefficients[2], digits = 3)
PGLS_r2 =      formatC( sum_fit_phy$r.squared, digits = 2)
PGLS_pvalue =  formatC(sum_fit_phy$coefficients[8], digits = 3)
if (PGLS_pvalue == "   0"){ PGLS_pvalue = "< 2.2e-16"}
PGLS_cc <- coef(fit_phylo_lamellae_SL)
PGLS_fit_function <- function(x) PGLS_cc[1] + PGLS_cc[2]*x
PGLS_lambda <- sum_fit_phy$param[2]
PGLS_r2


lm_fit <- lm(data = Lamellae_count_per_sp, 
             formula = mean_lamellae ~ mean_SL)

sum_fit <- summary(lm_fit)
GLS_r2 =      formatC( sum_fit$r.squared, digits = 3)
GLS_pente =   formatC(sum_fit$coefficients[2], digits = 3)
GLS_pvalue =  formatC(sum_fit$coefficients[8], digits = 3)
GLS_cc <- coef(lm_fit)
GLS_fit_function <- function(x) GLS_cc[1] + GLS_cc[2]*x

pdf(file = "Raw_R_plots/Lamellae_vs_SL.pdf",width = 8.34,  height = 4.61)

Lamellae_count_per_sp %>% 
  ggplot(aes( x = mean_SL, y= mean_lamellae)) +
  geom_point(aes(color=Tribe)) +
  stat_function(fun = GLS_fit_function, color="black")+
  labs(subtitle = bquote(R^2 ~ "=" ~ .(PGLS_r2) ~ ";" ~ P ~ "=" ~.(PGLS_pvalue) ~ ";" ~ lambda ~ "=" ~ .(PGLS_lambda))) +
  theme_classic() +
  xlab("Mean standart length") +
  ylab("Mean number of lamellae (cm)") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none") + 
  scale_color_manual(values=tribes.colors) 

dev.off()


#Now lets do the same at the scale of teleost 


caper_data_teleost_SL_lamellae <- 
  comparative.data(phy = teleost_phylo_pruned, 
                   data = rayfinned_lamellae_nb_inphylo,
                   names.col = fishtree_name, vcv = TRUE,
                   na.omit = FALSE, warn.dropped = TRUE)

caper_data_teleost_SL_lamellae_nonull <- 
  comparative.data(phy = teleost_phylo_pruned, 
                   data = rayfinned_lamellae_nb_inphylo %>% filter(Lamellae_number_mean > 0),
                   names.col = fishtree_name, vcv = TRUE,
                   na.omit = FALSE, warn.dropped = TRUE)



fit_phylo_lamellae_SL_teleost <- pgls(Lamellae_number_mean ~ Absolute_SL,
                              data = caper_data_teleost_SL_lamellae, 
                              lambda = "ML")
fit_phylo_lamellae_SL_teleost_nonull <- pgls(Lamellae_number_mean ~ Absolute_SL,
                                      data = caper_data_teleost_SL_lamellae_nonull, 
                                      lambda = "ML")


summary(fit_phylo_lamellae_SL_teleost)
summary(fit_phylo_lamellae_SL_teleost_nonull)




sum_fit_phy <- summary(fit_phylo_lamellae_SL_teleost)
PGLS_pente =   formatC(sum_fit_phy$coefficients[2], digits = 3)
PGLS_r2 =      formatC( sum_fit_phy$r.squared, digits = 2)
PGLS_pvalue =  formatC(sum_fit_phy$coefficients[8], digits = 3)
if (PGLS_pvalue == "   0"){ PGLS_pvalue = "< 2.2e-16"}
PGLS_cc <- coef(fit_phylo_lamellae_SL_teleost)
PGLS_fit_function <- function(x) PGLS_cc[1] + PGLS_cc[2]*x
PGLS_lambda <- sum_fit_phy$param[2]
PGLS_r2


lm_fit <- lm(data = rayfinned_lamellae_nb_inphylo, 
             formula = Lamellae_number_mean ~ Absolute_SL)

sum_fit <- summary(lm_fit)
GLS_r2 =      formatC( sum_fit$r.squared, digits = 3)
GLS_pente =   formatC(sum_fit$coefficients[2], digits = 3)
GLS_pvalue =  formatC(sum_fit$coefficients[8], digits = 3)
GLS_cc <- coef(lm_fit)
GLS_fit_function <- function(x) GLS_cc[1] + GLS_cc[2]*x

pdf(file = "Raw_R_plots/Teleost_Lamellae_vs_SL.pdf",width = 8.34,  height = 4.61)

rayfinned_lamellae_nb_inphylo %>% 
  ggplot(aes( x = Absolute_SL, y= Lamellae_number_mean)) +
  geom_point() +
  stat_function(fun = GLS_fit_function, color="black")+
  labs(subtitle = bquote(R^2 ~ "=" ~ .(PGLS_r2) ~ ";" ~ P ~ "=" ~.(PGLS_pvalue) ~ ";" ~ lambda ~ "=" ~ .(PGLS_lambda))) +
  theme_classic() +
  xlab("Mean standart length") +
  ylab("Mean number of lamellae (cm)") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none") 

dev.off()


#### Cichlids lamellae - Estimate ancestral state  ---------------------------------


cichlid_lamellae_matrix <- 
  Lamellae_count_per_sp %>%
  filter(Sp %in% b1_tree$tip.label) %>%
  dplyr::select(Sp, mean_lamellae)

cichlid_lamellae_matrix <- 
  data.frame(cichlid_lamellae_matrix[,2], 
             row.names=cichlid_lamellae_matrix[,1])


colnames(cichlid_lamellae_matrix) <- "lamellae_nb"
cichlid_lamellae_matrix <- 
  as.matrix(cichlid_lamellae_matrix)[,1]


#Estimate ancestral states

b1_tree_lamellae <- 
  keep.tip(radiation_tree_wo_Neospl, 
           Lamellae_count_per_sp %>% filter(Sp %in% b1_tree$tip.label) %>% pull(Sp))


fit_lamellae_cichlids <- 
  fastAnc(b1_tree_lamellae,
          cichlid_lamellae_matrix,
          vars=TRUE,CI=TRUE)


#Compute phylogenetic signal


tip_order <- b1_tree_lamellae$tip.label
lamellae_df_ordered <- 
  enframe(tip_order, name = NULL, value = "Sp") %>% 
  left_join(Lamellae_count_per_sp) %>% 
  column_to_rownames("Sp")

phylosig(tree = b1_tree_lamellae, x = lamellae_df_ordered$mean_lamellae, method = "lambda", test = T)


#### Cichlids lamellae - Visualize the ancestral states   ---------------------------------

#create a node number and node label correspondance table
b1_tree_lamellae <- makeNodeLabel(b1_tree_lamellae, method="number", prefix="Node")
species_test <- "Neomul"
misc_table_cichlid <- as.data.frame(species_test) %>% mutate(value = 1)
colnames(misc_table_cichlid) <- c("node","value")
misc_table_cichlid$node <- as.numeric(misc_table_cichlid$node)

node_label_corresp_cichlid <- 
  left_join(b1_tree_lamellae, misc_table_cichlid, by = 'node') %>%
  dplyr::select(node, label)

node_nb <- 
  node_label_corresp_cichlid %>%
  filter(., grepl("Node", label)) %>%
  pull(node)
species_names <- 
  node_label_corresp_cichlid %>%
  filter(., !grepl("Node", label)) %>%
  pull(label)
species_and_nodes <- c(species_names, node_nb)
node_label_corresp_cichlid <- node_label_corresp_cichlid %>% mutate(tree_names = species_and_nodes)
node_label_corresp_cichlid <- as.data.frame(node_label_corresp_cichlid)



#Visualize the phylogenetic tree of cichlids with the nb of lamellae


values_ancestral_cichlids <- 
  Lamellae_count_per_sp %>%
  filter(Sp %in% b1_tree$tip.label) %>%
  dplyr::select(Sp, mean_lamellae)


test <- as.data.frame(as.table(fit_lamellae_cichlids$ace))
colnames(test) <- c("Sp", "mean_lamellae")
values_ancestral_cichlids <- 
  rbind(values_ancestral_cichlids, test)
colnames(values_ancestral_cichlids) <- c("tree_names", "lamellae_nb")
values_ancestral_cichlids <- 
  left_join(values_ancestral_cichlids, 
            node_label_corresp_cichlid, 
            by="tree_names")
values_ancestral_cichlids <- 
  values_ancestral_cichlids %>% dplyr::select("label", "lamellae_nb")


ggtree(b1_tree_lamellae, layout="circular") %<+% values_ancestral_cichlids +
  aes(color=as.numeric(lamellae_nb)) +
  scale_color_continuous(name='Lamellae number', limits=c(0, 18),
                         oob=scales::squish, low="darkgreen", high="red") 

library(ggnewscale)
offspring.tbl_tree_item <- getFromNamespace("offspring", "tidytree")
assign("offspring.tbl_tree_item", offspring.tbl_tree_item, envir = .GlobalEnv)

pdf(file = "Raw_R_plots/Cichlids_anc_lamellae_nb.pdf",width = 4.34,  height = 6.61)

p1 <- ggtree(b1_tree_lamellae, size=1) %<+% values_ancestral_cichlids +
  aes(color=as.numeric(lamellae_nb)) +
  xlim(0, 21) + 
  scale_color_viridis(option="magma", direction=-1, na.value="white") +
  new_scale_color() + 
  geom_tiplab(size = 2) +
  theme(legend.position="none") 

node_rotate <- getMRCA(b1_tree_lamellae, c("Neobue", "Neomul"))
node_rotate_sec <- getMRCA(b1_tree_lamellae, c("Lepelo", "Neomul"))

p2 <- ggtree::rotate(p1, node_rotate) %>% ggtree::rotate(node_rotate_sec)
print(p2)


p1 <- ggtree(b1_tree_lamellae, size=1) %<+% values_ancestral_cichlids +
  aes(color=as.numeric(lamellae_nb)) +
  xlim(0, 21) + 
  scale_color_viridis(option="magma", direction=-1, na.value="white") +
  new_scale_color() + 
  geom_tiplab(size = 2)

p2 <- ggtree::rotate(p1, node_rotate) %>% ggtree::rotate(node_rotate_sec)
print(p2)

dev.off()


#### Cichlids lamellae - Pulse of diversification and DTT   ---------------------------------

#Compute the lambda of NL

phylosig(b1_tree_lamellae, cichlid_lamellae_matrix,
         test=TRUE, method="lambda")


#Draw phenogram

pdf(file = "Raw_R_plots/Phenogram_lamellae.pdf",width = 8.34,  height = 15.61)

phenogram(b1_tree_lamellae, cichlid_lamellae_matrix)

dev.off()

#### Data load - Brain volume -- Rayfinned fishes  ---------------------------------

#Import ray-finned fishes brain data from Hofmann et al =>  https://doi.org/10.1159/000530243

RayFinnedFishes_BrainVol_df <- 
  read.table("RayFinnedFishes_Brain_vol.tsv", sep="\t", header=TRUE)

#remove a species that is probably wrong (Optic tectum bigger than the total brain volume)
RayFinnedFishes_BrainVol_df <- RayFinnedFishes_BrainVol_df %>% filter(! Species == "Lipophrys_pholis")


#compute the relative volume of brain regions
RayFinnedFishes_BrainVol_df <-  RayFinnedFishes_BrainVol_df %>% dplyr::select(-c(H, SB))
colnames(RayFinnedFishes_BrainVol_df) <- 
  c("Species","Total", "Olfactory_bulb", "Telencephalon", 
    "Optic_tectum", "TL", "TLat", "IL", "Cerebellum", "CC", "FL", "VL") 

RayFinnedFishes_BrainVol_df <- 
  as.data.frame(RayFinnedFishes_BrainVol_df %>%
  dplyr::select(Species, Total, Olfactory_bulb, Telencephalon, Optic_tectum, Cerebellum) %>%
  rowwise() %>%
  mutate(Rest_of_the_brain = (Total - sum(Olfactory_bulb, Telencephalon, Optic_tectum, 
                                          Cerebellum, na.rm=TRUE))))

RayFinnedFishes_BrainRelVol_df <- 
  RayFinnedFishes_BrainVol_df %>%
  mutate(R_OB = Olfactory_bulb/Total,
         R_TEL = Telencephalon/Total,
         R_OT = Optic_tectum/Total,
         R_CB = Cerebellum/Total,
         R_ROB = Rest_of_the_brain/Total,
         R_Total = Total/Total) %>%
  dplyr::select(Species, R_OB, R_TEL, R_OT, R_CB, R_ROB, R_Total)
colnames(RayFinnedFishes_BrainRelVol_df) <- c("Species", "Olfactory_bulb", "Telencephalon",
                                              "Optic_tectum", "Cerebellum", "Rest_of_the_brain", "Total")

RayFinnedFishes_BrainVol_df_long <- 
  as.data.frame(RayFinnedFishes_BrainVol_df %>%
  pivot_longer(!Species, names_to = "Brain_region", values_to = "Volume"))

RayFinnedFishes_BrainRelVol_df_long <- 
  as.data.frame(RayFinnedFishes_BrainRelVol_df %>%
                  pivot_longer(!Species, names_to = "Brain_region", values_to = "Relative_volume"))


#Remove cichlids from the table

RayFinnedFishes_BrainRelVol_df_long_woCichlids <- 
  RayFinnedFishes_BrainRelVol_df_long %>%
  filter(! Species %in% c("Astronotus_ocellatus","Aequidens_pulcher","Heros_efasciatus",
                        "Heros_severus","Pterophyllum_scalare","Thorichthys_meeki",
                        "Apistogramma_agassizii","Mikrogeophagus_ramirezi","Maylandia_zebra",
                        "Julidochromis_marlieri","Neolamprologus_brichardi","Steatocranus_casuarius"))


RayFinnedFishes_BrainVol_df_long_woCichlids <- 
  RayFinnedFishes_BrainVol_df_long %>%
  filter(! Species %in% c("Astronotus_ocellatus","Aequidens_pulcher","Heros_efasciatus",
                          "Heros_severus","Pterophyllum_scalare","Thorichthys_meeki",
                          "Apistogramma_agassizii","Mikrogeophagus_ramirezi","Maylandia_zebra",
                          "Julidochromis_marlieri","Neolamprologus_brichardi","Steatocranus_casuarius"))

RayFinnedFishes_BrainVol_df_long_woCichlids <-
  RayFinnedFishes_BrainVol_df_long_woCichlids %>%
  mutate(Log_vol = log10(Volume))



#### Data load - Brain volume -- Cichlids  ---------------------------------


Policarpo_Cichlids_Brain_df <- 
  read.table("OlfactoryBulb_data.tsv", sep="\t", header=TRUE)
colnames(Policarpo_Cichlids_Brain_df) <- c("ID", "Brain_region", "Vol", "Sp", "Log_vol", "Prop_region")

Policarpo_Cichlids_Brain_df[(Policarpo_Cichlids_Brain_df$Brain_region == "Rest of the brain"),"Brain_region"] <- "Rest_of_the_brain"
Policarpo_Cichlids_Brain_df[(Policarpo_Cichlids_Brain_df$Brain_region == "Olfactory bulb"),"Brain_region"] <- "Olfactory_bulb"
Policarpo_Cichlids_Brain_df[(Policarpo_Cichlids_Brain_df$Brain_region == "Optic tectum"),"Brain_region"] <- "Optic_tectum"
Policarpo_Cichlids_Brain_df <- Policarpo_Cichlids_Brain_df %>% mutate(Study = "Policarpo")


#Add  cichlids data from the paper of Hofmann to our table

RayFinnedFishes_BrainVol_df_long_cichlid <- 
  RayFinnedFishes_BrainVol_df_long %>%
  filter(Species %in% c("Astronotus_ocellatus","Aequidens_pulcher","Heros_efasciatus",
                         "Heros_severus","Pterophyllum_scalare","Thorichthys_meeki",
                         "Apistogramma_agassizii","Mikrogeophagus_ramirezi","Maylandia_zebra",
                         "Julidochromis_marlieri","Neolamprologus_brichardi","Steatocranus_casuarius"))
colnames(RayFinnedFishes_BrainVol_df_long_cichlid) <- c("ID", "Brain_region", "Vol")

RayFinnedFishes_BrainRelVol_df_long_cichlid <- 
  RayFinnedFishes_BrainRelVol_df_long %>%
  filter(Species %in% c("Astronotus_ocellatus","Aequidens_pulcher","Heros_efasciatus",
                        "Heros_severus","Pterophyllum_scalare","Thorichthys_meeki",
                        "Apistogramma_agassizii","Mikrogeophagus_ramirezi","Maylandia_zebra",
                        "Julidochromis_marlieri","Neolamprologus_brichardi","Steatocranus_casuarius"))
colnames(RayFinnedFishes_BrainRelVol_df_long_cichlid) <- c("ID", "Brain_region", "Prop_region")


RayFinnedFishes_BrainVol_df_long_cichlid <- RayFinnedFishes_BrainVol_df_long_cichlid %>%
  mutate(Sp = case_when(
    ID == "Julidochromis_marlieri" ~ "JulmaN",
    ID == "Astronotus_ocellatus" ~ "Astoce",
    ID == "Aequidens_pulcher" ~ "Aeqpul",
    ID == "Heros_efasciatus" ~ "Herefa",
    ID == "Heros_severus" ~ "Hersev",
    ID == "Pterophyllum_scalare" ~ "Petsca",
    ID == "Thorichthys_meeki" ~ "Thomee",
    ID == "Apistogramma_agassizii" ~ "Apiaga",
    ID == "Mikrogeophagus_ramirezi" ~ "Mikram",
    ID == "Maylandia_zebra" ~ "Mayzeb",
    ID == "Neolamprologus_brichardi" ~ "Neobri",
    ID == "Steatocranus_casuarius" ~ "Stecas"
  )) %>%
  rowwise() %>%
  mutate(Log_vol = log10(Vol)) %>%
  mutate(Study = "Hofmann")


RayFinnedFishes_BrainRelVol_df_long_cichlid <- RayFinnedFishes_BrainRelVol_df_long_cichlid %>%
  mutate(Sp = case_when(
    ID == "Julidochromis_marlieri" ~ "JulmaN",
    ID == "Astronotus_ocellatus" ~ "Astoce",
    ID == "Aequidens_pulcher" ~ "Aeqpul",
    ID == "Heros_efasciatus" ~ "Herefa",
    ID == "Heros_severus" ~ "Hersev",
    ID == "Pterophyllum_scalare" ~ "Petsca",
    ID == "Thorichthys_meeki" ~ "Thomee",
    ID == "Apistogramma_agassizii" ~ "Apiaga",
    ID == "Mikrogeophagus_ramirezi" ~ "Mikram",
    ID == "Maylandia_zebra" ~ "Mayzeb",
    ID == "Neolamprologus_brichardi" ~ "Neobri",
    ID == "Steatocranus_casuarius" ~ "Stecas"
  )) %>%
  mutate(Study = "Hofmann")



Hofmann_Cichlids_Brain_df <- 
  left_join(RayFinnedFishes_BrainVol_df_long_cichlid, RayFinnedFishes_BrainRelVol_df_long_cichlid,
            by=c("ID", "Sp", "Brain_region", "Study"))

Hofmann_Cichlids_Brain_df <- 
  Hofmann_Cichlids_Brain_df %>%
  dplyr::select(ID, Brain_region, Vol, Sp, Log_vol, Prop_region, Study)




Policarpo_Hofmann_Cichlids_Brain_df <- 
  rbind(Policarpo_Cichlids_Brain_df,
        Hofmann_Cichlids_Brain_df)


#Add the standart length to the table 

ID_SL <- Lamellae_info_df %>% dplyr::select(ID, SL)

Policarpo_Hofmann_Cichlids_Brain_df <- 
  left_join(Policarpo_Hofmann_Cichlids_Brain_df, ID_SL, by="ID")



Policarpo_Hofmann_Cichlids_Brain_df[(Policarpo_Hofmann_Cichlids_Brain_df$Sp == "JulmaN"),"SL"] <- 6.406 #labkey
Policarpo_Hofmann_Cichlids_Brain_df[(Policarpo_Hofmann_Cichlids_Brain_df$Sp == "Astoce"),"SL"] <- 12.05 #labkey
Policarpo_Hofmann_Cichlids_Brain_df[(Policarpo_Hofmann_Cichlids_Brain_df$Sp == "Aeqpul"),"SL"] <- 12 #fishbase
Policarpo_Hofmann_Cichlids_Brain_df[(Policarpo_Hofmann_Cichlids_Brain_df$Sp == "Herefa"),"SL"] <- 11.6 #labkey
Policarpo_Hofmann_Cichlids_Brain_df[(Policarpo_Hofmann_Cichlids_Brain_df$Sp == "Hersev"),"SL"] <- 10.3 #fishbase
Policarpo_Hofmann_Cichlids_Brain_df[(Policarpo_Hofmann_Cichlids_Brain_df$Sp == "Petsca"),"SL"] <- 6.72 #labkey
Policarpo_Hofmann_Cichlids_Brain_df[(Policarpo_Hofmann_Cichlids_Brain_df$Sp == "Thomee"),"SL"] <- 4.8739 #fishbase
Policarpo_Hofmann_Cichlids_Brain_df[(Policarpo_Hofmann_Cichlids_Brain_df$Sp == "Apiaga"),"SL"] <- 3.1 #fishbase
Policarpo_Hofmann_Cichlids_Brain_df[(Policarpo_Hofmann_Cichlids_Brain_df$Sp == "Mikram"),"SL"] <- 3 #fishbase
Policarpo_Hofmann_Cichlids_Brain_df[(Policarpo_Hofmann_Cichlids_Brain_df$Sp == "Mayzeb"),"SL"] <- 11 #fishbase
Policarpo_Hofmann_Cichlids_Brain_df[(Policarpo_Hofmann_Cichlids_Brain_df$Sp == "Neobri"),"SL"] <- 6.095 #labkey
Policarpo_Hofmann_Cichlids_Brain_df[(Policarpo_Hofmann_Cichlids_Brain_df$Sp == "Stecas"),"SL"] <- 8.39 #fishbase

#Add Tribe info

Policarpo_Hofmann_Cichlids_Brain_df <- 
  left_join(Policarpo_Hofmann_Cichlids_Brain_df, Tribe_for_tree, by="Sp")

Policarpo_Hofmann_Cichlids_Brain_df$Tribe[is.na(Policarpo_Hofmann_Cichlids_Brain_df$Tribe)] <- "Other"


#Make the mean per species

Policarpo_Hofmann_Cichlids_Brain_df_meansp <- 
  as.data.frame(
    Policarpo_Hofmann_Cichlids_Brain_df %>%
      group_by(Sp, Brain_region, Study, Tribe) %>%
      summarise(mean_prop = mean(Prop_region),
                mean_vol = mean(Vol),
                mean_log_vol = mean(Log_vol),
                mean_SL = mean(SL)))



#Make a wide dataframe from the long one

Policarpo_Hofmann_Cichlids_Brain_df_meansp_wide <- 
  as.data.frame(
    Policarpo_Hofmann_Cichlids_Brain_df_meansp %>%
  dplyr::select(Sp, Tribe, Study, Brain_region, mean_prop) %>%
  pivot_wider(names_from = Brain_region, values_from = mean_prop, values_fill = 0)
  )


#write a table for the paper


write.table(Policarpo_Hofmann_Cichlids_Brain_df %>% filter(Study == "Policarpo") %>% dplyr::select(-Study),
            "DataSupp_paper_Brain.tsv", 
            sep="\t",
            col.names=TRUE,
            row.names=FALSE,
            quote=FALSE)

#### Cichlids brain - Compare to other cichlids and to ray-finned fishes---------------------------------

#Iniate a new color palette combining sp and tribe

Tribe_Sp <- Policarpo_Hofmann_Cichlids_Brain_df %>% dplyr::select(Tribe, Sp)
tribes.colors_df <- as.data.frame(tribes.colors)
tribes.colors_df$Tribe <- row.names(tribes.colors_df)
Tribe_Sp_color <- left_join(Tribe_Sp, tribes.colors_df, by="Tribe")

brain.sp_colors_df <- Tribe_Sp_color %>% dplyr::select(Sp, tribes.colors) %>% distinct()
brain.sp_colors <- as.character(brain.sp_colors_df$tribes.colors)
names(brain.sp_colors) <- brain.sp_colors_df$Sp


brain.sp_colors <- c(brain.sp_colors, "Hofmann"="gray")

#First look at the relation between the brain volume and the SL

Policarpo_Hofmann_Cichlids_Brain_df %>%
  filter(Study == "Policarpo") %>%
  dplyr::select(Sp, Tribe) %>%
  distinct()

vol_vs_SL_lm <- lm(Vol ~ SL,
   data=Policarpo_Hofmann_Cichlids_Brain_df %>% 
     filter(Study == "Policarpo") %>% 
     filter(! is.na(SL)) %>% 
     filter(Brain_region == "Total"))

sum_fit <- summary(vol_vs_SL_lm)
GLS_r2 =      formatC( sum_fit$r.squared, digits = 3)
GLS_pente =   formatC(sum_fit$coefficients[2], digits = 3)
GLS_pvalue =  formatC(sum_fit$coefficients[8], digits = 3)
GLS_cc <- coef(vol_vs_SL_lm)
GLS_fit_function <- function(x) GLS_cc[1] + GLS_cc[2]*x


pdf(file = "Raw_R_plots/BrainVolume_vs_SL.pdf",width = 8.34,  height = 4.61)

Policarpo_Hofmann_Cichlids_Brain_df %>%
  filter(Study == "Policarpo") %>%
  filter(Brain_region == "Total") %>%
  ggplot(., aes(x=SL, y=Vol, col=Tribe, shape=Sp)) +
  geom_point(size = 2) +
  stat_function(fun = GLS_fit_function, color="black")+
  theme_classic() +
  scale_color_manual(values = tribes.colors) +
  scale_shape_manual(values = c(
    "Asplep" = 16,
    "Cunlon" = 15,
    "Cphgib" = 16,
    "Pcybri" = 16,
    "Petpol" = 16,
    "Benmel" = 16,
    "Batvit" = 16,
    "Baicen" = 16,
    "Tremar" = 16,
    "Neomul"= 15,
    "Lepken" = 16
  )) +
  theme(legend.position = "none") +
  xlab("Standart length (cm)") +
  ylab("Brain volume (mm3)") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none") + 
  scale_color_manual(values=tribes.colors) 

dev.off()


#Now look at the variation of the relative volumes of the different brain regions

hofman_sp_list <- Policarpo_Hofmann_Cichlids_Brain_df %>% filter(Study != "Policarpo") %>% pull(Sp) %>% unique()

Policarpo_Hofmann_Cichlids_Brain_df$Brain_region <-
  factor(Policarpo_Hofmann_Cichlids_Brain_df$Brain_region,
         levels=c("Olfactory_bulb", "Telencephalon", "Optic_tectum",
                  "Cerebellum", "Rest_of_the_brain", "Total"
         ))


Policarpo_Hofmann_Cichlids_Brain_df$Sp <-
  factor(Policarpo_Hofmann_Cichlids_Brain_df$Sp ,
         levels=c("Tremar", "Batvit", "Pcybri", "Benmel", "Petpol", "Cphgib", "Baicen",
                  "Asplep", "Cunlon","Neomul", "Lepken", hofman_sp_list
                  ))

Policarpo_Hofmann_Cichlids_Brain_df <- 
  Policarpo_Hofmann_Cichlids_Brain_df %>%
  mutate(Sp_or_other = if_else(
    Study == "Policarpo",
    Sp,
    "Hofmann"
  ))


Policarpo_Hofmann_Cichlids_Brain_df$Sp_or_other <-
  factor(Policarpo_Hofmann_Cichlids_Brain_df$Sp_or_other ,
         levels=c("Tremar", "Batvit", "Pcybri", "Benmel", "Petpol", "Cphgib", "Baicen",
                  "Asplep", "Cunlon","Neomul", "Lepken", "Hofmann"
         ))


pdf(file = "Raw_R_plots/Regions_volumes_between_sp.pdf",width = 8.34,  height = 4.61)

Policarpo_Hofmann_Cichlids_Brain_df %>%
  #filter(Study == "Policarpo") %>%
  filter(Brain_region == "Rest_of_the_brain") %>%
  ggplot(., aes(x=Brain_region, y=Prop_region, fill=Sp_or_other)) +
  geom_boxplot() +
  scale_fill_manual(values = brain.sp_colors) +
  theme(legend.position = "none") +
  xlab("") +
  ylab("Relative volume") +
  theme_classic() +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none") 

Policarpo_Hofmann_Cichlids_Brain_df %>%
  #filter(Study == "Policarpo") %>%
  filter(Brain_region == "Olfactory_bulb") %>%
  ggplot(., aes(x=Brain_region, y=Prop_region, fill=Sp_or_other)) +
  geom_boxplot() +
  scale_fill_manual(values = brain.sp_colors) +
  theme(legend.position = "none") +
  xlab("") +
  ylab("Relative volume") +
  theme_classic() +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none") 

Policarpo_Hofmann_Cichlids_Brain_df %>%
  #filter(Study == "Policarpo") %>%
  filter(Brain_region == "Telencephalon") %>%
  ggplot(., aes(x=Brain_region, y=Prop_region, fill=Sp_or_other)) +
  geom_boxplot() +
  scale_fill_manual(values = brain.sp_colors) +
  theme(legend.position = "none") +
  xlab("") +
  ylab("Relative volume") +
  theme_classic() +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none")  

Policarpo_Hofmann_Cichlids_Brain_df %>%
  #filter(Study == "Policarpo") %>%
  filter(Brain_region == "Cerebellum") %>%
  ggplot(., aes(x=Brain_region, y=Prop_region, fill=Sp_or_other)) +
  geom_boxplot() +
  scale_fill_manual(values = brain.sp_colors) +
  theme(legend.position = "none") +
  xlab("") +
  ylab("Relative volume") +
  theme_classic() +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none")  

Policarpo_Hofmann_Cichlids_Brain_df %>%
  #filter(Study == "Policarpo") %>%
  filter(Brain_region == "Optic_tectum") %>%
  ggplot(., aes(x=Brain_region, y=Prop_region, fill=Sp_or_other)) +
  geom_boxplot() +
  scale_fill_manual(values = brain.sp_colors) +
  theme(legend.position = "none") +
  xlab("") +
  ylab("Relative volume") +
  theme_classic() +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none") 



dev.off()



#Now compare brain volumes to other ray-finned fishes (non cichlids only)...

Policarpo_cichlids_meansp_df <- 
  Policarpo_Hofmann_Cichlids_Brain_df_meansp %>%
  filter(Study == "Policarpo") %>%
  dplyr::select(Sp, Brain_region, mean_prop)


colnames(Policarpo_cichlids_meansp_df) <- 
  c("Species", "Brain_region", "Relative_volume")


cichlid_teleost_color <- 
  c(Cichlid="#004D40",
    Teleost = "#FFC107")



RayFinned_vs_Cichlids_brain_df <- 
  rbind(RayFinnedFishes_BrainRelVol_df_long_woCichlids %>% mutate(clade = "Teleost"),
        Policarpo_cichlids_meansp_df %>% mutate(clade = "Cichlid"))

pdf(file = "Raw_R_plots/Cichlids_vs_Teleost_BrainsRegions.pdf",width = 8.34,  height = 4.61)

RayFinned_vs_Cichlids_brain_df %>%
  filter(Brain_region == "Olfactory_bulb") %>%
  ggplot(., aes(x=Relative_volume, fill=clade)) +
  geom_histogram(bins=40, alpha = 0.4, position = "identity", color="black") +
  theme_classic() +
  scale_fill_manual(values = cichlid_teleost_color) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position = "none") +
  xlab("Relative olfactory bulb volume") +
  ylab("Number of species")

  
RayFinned_vs_Cichlids_brain_df %>%
  filter(Brain_region == "Telencephalon") %>%
  ggplot(., aes(x=Relative_volume, fill=clade)) +
  geom_histogram(bins=40, alpha = 0.4, position = "identity", color="black") +
  theme_classic() +
  scale_fill_manual(values = cichlid_teleost_color) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position = "none") +
  xlab("Relative telencephalon volume") +
  ylab("Number of species")


RayFinned_vs_Cichlids_brain_df %>%
  filter(Brain_region == "Optic_tectum") %>%
  ggplot(., aes(x=Relative_volume, fill=clade)) +
  geom_histogram(bins=40, alpha = 0.4, position = "identity", color="black") +
  theme_classic() +
  scale_fill_manual(values = cichlid_teleost_color) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position = "none") +
  xlab("Relative optic tectum volume") +
  ylab("Number of species")


RayFinned_vs_Cichlids_brain_df %>%
  filter(Brain_region == "Cerebellum") %>%
  ggplot(., aes(x=Relative_volume, fill=clade)) +
  geom_histogram(bins=40, alpha = 0.4, position = "identity", color="black") +
  theme_classic() +
  scale_fill_manual(values = cichlid_teleost_color) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position = "none") +
  xlab("Relative cerebellum volume") +
  ylab("Number of species")


RayFinned_vs_Cichlids_brain_df %>%
  filter(Brain_region == "Rest_of_the_brain") %>%
  ggplot(., aes(x=Relative_volume, fill=clade)) +
  geom_histogram(bins=40, alpha = 0.4, position = "identity", color="black") +
  theme_classic() +
  scale_fill_manual(values = cichlid_teleost_color) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position = "none") +
  xlab("Relative rest of the brain volume") +
  ylab("Number of species")

dev.off()




#### Cichlids Brain -  Interesting species level stats ---------------------------------

Policarpo_Hofmann_Cichlids_Brain_df_meansp %>%
  filter(Study == "Policarpo") %>%
  filter(Brain_region == "Olfactory_bulb") %>%
  arrange(mean_prop)



#### Prepare dataframe and caper data for pGLS analysis  ---------------------------------

#Df with number of genes
All_Doc_mean_sp_wide <- 
  All_Doc_mean_sp_wide %>%
  mutate(Total_OLR = Total_OR + Total_TAAR + Total_V1R + Total_V2R) %>%
  filter(Sp != "Neospl")


#Df with brain relative volumes
Policarpo_Hofmann_Cichlids_Brain_df_meansp_wide_prop <- 
  as.data.frame(
    Policarpo_Hofmann_Cichlids_Brain_df_meansp %>%
      dplyr::select(Sp, Tribe, Study, Brain_region, mean_prop) %>%
      pivot_wider(names_from = Brain_region, values_from = mean_prop, values_fill = 0)
  )
colnames(Policarpo_Hofmann_Cichlids_Brain_df_meansp_wide_prop) <- 
  c("Sp", "Tribe", "Study",  "prop_Cerebellum", "prop_Olfactory_bulb", "prop_Optic_tectum",
    "prop_Rest_of_the_brain", "prop_Telencephalon", "prop_Total")

#Df with brain absolute volumes
Policarpo_Hofmann_Cichlids_Brain_df_meansp_wide_vol <- 
  as.data.frame(
    Policarpo_Hofmann_Cichlids_Brain_df_meansp %>%
      dplyr::select(Sp, Tribe, Study, Brain_region, mean_vol) %>%
      pivot_wider(names_from = Brain_region, values_from = mean_vol, values_fill = 0)
  )
colnames(Policarpo_Hofmann_Cichlids_Brain_df_meansp_wide_vol) <- 
  c("Sp", "Tribe", "Study",  "vol_Cerebellum", "vol_Olfactory_bulb", "vol_Optic_tectum",
    "vol_Rest_of_the_brain", "vol_Telencephalon", "vol_Total")


#Join brain dfs
Policarpo_Hofmann_Cichlids_Brain_df_meansp_wide <- 
  left_join(Policarpo_Hofmann_Cichlids_Brain_df_meansp_wide_vol, Policarpo_Hofmann_Cichlids_Brain_df_meansp_wide_prop,
            by=c("Sp", "Tribe", "Study"))


#Merge with number of genes
Cichlids_lamellae_brains_df <- 
  left_join(All_Doc_mean_sp_wide, Policarpo_Hofmann_Cichlids_Brain_df_meansp_wide,
            by=c("Sp", "Tribe"))


#Merge with the SL + number of lamellae
Cichlids_lamellae_brains_df <- 
  left_join(Cichlids_lamellae_brains_df, Lamellae_count_per_sp,
            by=c("Sp", "Tribe"))




caper_data_lamellae_brain <- 
  comparative.data(phy = b1_tree, 
                   data = Cichlids_lamellae_brains_df,
                   names.col = Sp, vcv = TRUE,
                   na.omit = FALSE, warn.dropped = TRUE)


b1_tree_lamellae <- keep.tip(b1_tree, 
                             Cichlids_lamellae_brains_df %>% 
                               filter(Sp %in% b1_tree$tip.label) %>%
                               filter(! Sp %in% c("Oretan", "TelteN", "Tylpol")) %>%
                               pull(Sp))


write.table(Cichlids_lamellae_brains_df,
            "Cichlids_lamellae_brains_df.csv",
            sep=",",
            col.names=TRUE,
            row.names=FALSE,
            quote=FALSE)

#### Cichlids - pGLS Lamellae vs Relative OB  ---------------------------------


fit_phylo_lamellae_OB <- pgls(prop_Olfactory_bulb ~ mean_lamellae,
                              data = caper_data_lamellae_brain, 
                              lambda = "ML")
summary(fit_phylo_lamellae_OB)

sum_fit_phy <- summary(fit_phylo_lamellae_OB)
PGLS_pente =   formatC(sum_fit_phy$coefficients[2], digits = 3)
PGLS_r2 =      formatC( sum_fit_phy$r.squared, digits = 2)
PGLS_pvalue =  formatC(sum_fit_phy$coefficients[8], digits = 3)
if (PGLS_pvalue == "   0"){ PGLS_pvalue = "< 2.2e-16"}
PGLS_cc <- coef(fit_phylo_lamellae_OB)
PGLS_fit_function <- function(x) PGLS_cc[1] + PGLS_cc[2]*x
PGLS_lambda <- sum_fit_phy$param[2]
PGLS_r2


lm_fit <- lm(data = Cichlids_lamellae_brains_df %>% filter(Sp %in% b1_tree_lamellae$tip.label), 
             formula = prop_Olfactory_bulb ~ mean_lamellae)

sum_fit <- summary(lm_fit)
GLS_r2 =      formatC( sum_fit$r.squared, digits = 3)
GLS_pente =   formatC(sum_fit$coefficients[2], digits = 3)
GLS_pvalue =  formatC(sum_fit$coefficients[8], digits = 3)
GLS_cc <- coef(lm_fit)
GLS_fit_function <- function(x) GLS_cc[1] + GLS_cc[2]*x

pdf(file = "Raw_R_plots/OB_vs_lamellae.pdf",width = 8.34,  height = 4.61)

Cichlids_lamellae_brains_df %>% 
  ggplot(aes( x = mean_lamellae, y= prop_Olfactory_bulb)) +
  geom_point(aes(color=Tribe)) +
  stat_function(fun = GLS_fit_function, color="black")+
  labs(subtitle = bquote(R^2 ~ "=" ~ .(PGLS_r2) ~ ";" ~ P ~ "=" ~.(PGLS_pvalue) ~ ";" ~ lambda ~ "=" ~ .(PGLS_lambda))) +
  theme_classic() +
  xlab("Mean number of lamellae") +
  ylab("Mean relative olfactory bulb volume") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none") + 
  scale_color_manual(values=tribes.colors) 

dev.off()



#### Cichlids - pGLS Lamellae vs Relative ROB  ---------------------------------


fit_phylo_lamellae_ROB <- pgls(prop_Rest_of_the_brain ~ mean_lamellae,
                               data = caper_data_lamellae_brain, 
                               lambda = "ML")
summary(fit_phylo_lamellae_ROB)

sum_fit_phy <- summary(fit_phylo_lamellae_ROB)
PGLS_pente =   formatC(sum_fit_phy$coefficients[2], digits = 3)
PGLS_r2 =      formatC( sum_fit_phy$r.squared, digits = 2)
PGLS_pvalue =  formatC(sum_fit_phy$coefficients[8], digits = 3)
if (PGLS_pvalue == "   0"){ PGLS_pvalue = "< 2.2e-16"}
PGLS_cc <- coef(fit_phylo_lamellae_ROB)
PGLS_fit_function <- function(x) PGLS_cc[1] + PGLS_cc[2]*x
PGLS_lambda <- sum_fit_phy$param[2]
PGLS_r2


lm_fit <- lm(data = Cichlids_lamellae_brains_df %>% filter(Sp %in% b1_tree_lamellae$tip.label), 
             formula = prop_Rest_of_the_brain ~ mean_lamellae)

sum_fit <- summary(lm_fit)
GLS_r2 =      formatC( sum_fit$r.squared, digits = 3)
GLS_pente =   formatC(sum_fit$coefficients[2], digits = 3)
GLS_pvalue =  formatC(sum_fit$coefficients[8], digits = 3)
GLS_cc <- coef(lm_fit)
GLS_fit_function <- function(x) GLS_cc[1] + GLS_cc[2]*x

pdf(file = "Raw_R_plots/ROB_vs_lamellae.pdf",width = 8.34,  height = 4.61)

Cichlids_lamellae_brains_df %>% 
  ggplot(aes( x = mean_lamellae, y= prop_Rest_of_the_brain)) +
  geom_point(aes(color=Tribe, size=2)) +
  stat_function(fun = GLS_fit_function, color="black")+
  labs(subtitle = bquote(R^2 ~ "=" ~ .(PGLS_r2) ~ ";" ~ P ~ "=" ~.(PGLS_pvalue) ~ ";" ~ lambda ~ "=" ~ .(PGLS_lambda))) +
  theme_classic() +
  xlab("Mean number of lamellae") +
  ylab("Mean rest of the brain volume") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none") + 
  scale_color_manual(values=tribes.colors) 

dev.off()




#### Cichlids - pGLS Lamellae vs Relative Cerebellum  ---------------------------------


fit_phylo_lamellae_Cerebellum <- pgls(prop_Cerebellum ~ mean_lamellae,
                              data = caper_data_lamellae_brain, 
                              lambda = "ML")
summary(fit_phylo_lamellae_Cerebellum)

sum_fit_phy <- summary(fit_phylo_lamellae_Cerebellum)
PGLS_pente =   formatC(sum_fit_phy$coefficients[2], digits = 3)
PGLS_r2 =      formatC( sum_fit_phy$r.squared, digits = 2)
PGLS_pvalue =  formatC(sum_fit_phy$coefficients[8], digits = 3)
if (PGLS_pvalue == "   0"){ PGLS_pvalue = "< 2.2e-16"}
PGLS_cc <- coef(fit_phylo_lamellae_Cerebellum)
PGLS_fit_function <- function(x) PGLS_cc[1] + PGLS_cc[2]*x
PGLS_lambda <- sum_fit_phy$param[2]
PGLS_r2


lm_fit <- lm(data = Cichlids_lamellae_brains_df %>% filter(Sp %in% b1_tree_lamellae$tip.label), 
             formula = prop_Cerebellum ~ mean_lamellae)

sum_fit <- summary(lm_fit)
GLS_r2 =      formatC( sum_fit$r.squared, digits = 3)
GLS_pente =   formatC(sum_fit$coefficients[2], digits = 3)
GLS_pvalue =  formatC(sum_fit$coefficients[8], digits = 3)
GLS_cc <- coef(lm_fit)
GLS_fit_function <- function(x) GLS_cc[1] + GLS_cc[2]*x

pdf(file = "Raw_R_plots/Cerebellum_vs_lamellae.pdf",width = 8.34,  height = 4.61)

Cichlids_lamellae_brains_df %>% 
  ggplot(aes( x = mean_lamellae, y= prop_Cerebellum)) +
  geom_point(aes(color=Tribe, size=2)) +
  stat_function(fun = GLS_fit_function, color="black")+
  labs(subtitle = bquote(R^2 ~ "=" ~ .(PGLS_r2) ~ ";" ~ P ~ "=" ~.(PGLS_pvalue) ~ ";" ~ lambda ~ "=" ~ .(PGLS_lambda))) +
  theme_classic() +
  xlab("Mean number of lamellae") +
  ylab("Mean relative cerebellum volume") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none") + 
  scale_color_manual(values=tribes.colors) 

dev.off()

#### Cichlids - pGLS Lamellae vs Relative Optic_tectum  ---------------------------------


fit_phylo_lamellae_Optic_tectum <- pgls(prop_Optic_tectum ~ mean_lamellae,
                                        data = caper_data_lamellae_brain, 
                                        lambda = "ML")
summary(fit_phylo_lamellae_Optic_tectum)

sum_fit_phy <- summary(fit_phylo_lamellae_Optic_tectum)
PGLS_pente =   formatC(sum_fit_phy$coefficients[2], digits = 3)
PGLS_r2 =      formatC( sum_fit_phy$r.squared, digits = 2)
PGLS_pvalue =  formatC(sum_fit_phy$coefficients[8], digits = 3)
if (PGLS_pvalue == "   0"){ PGLS_pvalue = "< 2.2e-16"}
PGLS_cc <- coef(fit_phylo_lamellae_Optic_tectum)
PGLS_fit_function <- function(x) PGLS_cc[1] + PGLS_cc[2]*x
PGLS_lambda <- sum_fit_phy$param[2]
PGLS_r2


lm_fit <- lm(data = Cichlids_lamellae_brains_df %>% filter(Sp %in% b1_tree_lamellae$tip.label), 
             formula = prop_Optic_tectum ~ mean_lamellae)

sum_fit <- summary(lm_fit)
GLS_r2 =      formatC( sum_fit$r.squared, digits = 3)
GLS_pente =   formatC(sum_fit$coefficients[2], digits = 3)
GLS_pvalue =  formatC(sum_fit$coefficients[8], digits = 3)
GLS_cc <- coef(lm_fit)
GLS_fit_function <- function(x) GLS_cc[1] + GLS_cc[2]*x

pdf(file = "Raw_R_plots/Optic_tectum_vs_lamellae.pdf",width = 8.34,  height = 4.61)

Cichlids_lamellae_brains_df %>% 
  ggplot(aes( x = mean_lamellae, y= prop_Optic_tectum)) +
  geom_point(aes(color=Tribe, size = 2)) +
  stat_function(fun = GLS_fit_function, color="black")+
  labs(subtitle = bquote(R^2 ~ "=" ~ .(PGLS_r2) ~ ";" ~ P ~ "=" ~.(PGLS_pvalue) ~ ";" ~ lambda ~ "=" ~ .(PGLS_lambda))) +
  theme_classic() +
  xlab("Mean number of lamellae") +
  ylab("Mean relative optic tectum volume") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none") + 
  scale_color_manual(values=tribes.colors) 

dev.off()

#### Cichlids - pGLS Lamellae vs Relative Telencephalon  ---------------------------------


fit_phylo_lamellae_Telencephalon <- pgls(prop_Telencephalon ~ mean_lamellae,
                                         data = caper_data_lamellae_brain, 
                                         lambda = "ML")
summary(fit_phylo_lamellae_Telencephalon)

sum_fit_phy <- summary(fit_phylo_lamellae_Telencephalon)
PGLS_pente =   formatC(sum_fit_phy$coefficients[2], digits = 3)
PGLS_r2 =      formatC( sum_fit_phy$r.squared, digits = 2)
PGLS_pvalue =  formatC(sum_fit_phy$coefficients[8], digits = 3)
if (PGLS_pvalue == "   0"){ PGLS_pvalue = "< 2.2e-16"}
PGLS_cc <- coef(fit_phylo_lamellae_Telencephalon)
PGLS_fit_function <- function(x) PGLS_cc[1] + PGLS_cc[2]*x
PGLS_lambda <- sum_fit_phy$param[2]
PGLS_r2


lm_fit <- lm(data = Cichlids_lamellae_brains_df %>% filter(Sp %in% b1_tree_lamellae$tip.label), 
             formula = prop_Telencephalon ~ mean_lamellae)

sum_fit <- summary(lm_fit)
GLS_r2 =      formatC( sum_fit$r.squared, digits = 3)
GLS_pente =   formatC(sum_fit$coefficients[2], digits = 3)
GLS_pvalue =  formatC(sum_fit$coefficients[8], digits = 3)
GLS_cc <- coef(lm_fit)
GLS_fit_function <- function(x) GLS_cc[1] + GLS_cc[2]*x

pdf(file = "Raw_R_plots/Telencephalon_vs_lamellae.pdf",width = 8.34,  height = 4.61)

Cichlids_lamellae_brains_df %>% 
  ggplot(aes( x = mean_lamellae, y= prop_Telencephalon)) +
  geom_point(aes(color=Tribe, size=2)) +
  stat_function(fun = GLS_fit_function, color="black")+
  labs(subtitle = bquote(R^2 ~ "=" ~ .(PGLS_r2) ~ ";" ~ P ~ "=" ~.(PGLS_pvalue) ~ ";" ~ lambda ~ "=" ~ .(PGLS_lambda))) +
  theme_classic() +
  xlab("Mean number of lamellae") +
  ylab("Mean relative telencephalon volume") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none") + 
  scale_color_manual(values=tribes.colors) 

dev.off()

#### Cichlids - pGLS Lamellae vs Olfactory receptors  ---------------------------------


#Perform the pGLS
fit_phylo_lamellae_OLR <- pgls(Total_OLR ~ mean_lamellae,
                               data = caper_data_lamellae_brain, 
                               lambda = "ML")
fit_phylo_lamellae_OR <- pgls(Total_OR ~ mean_lamellae,
                              data = caper_data_lamellae_brain, 
                              lambda = "ML")
fit_phylo_lamellae_TAAR <- pgls(Total_TAAR ~ mean_lamellae,
                              data = caper_data_lamellae_brain, 
                              lambda = "ML")
fit_phylo_lamellae_V2R <- pgls(Total_V2R ~ mean_lamellae,
                              data = caper_data_lamellae_brain, 
                              lambda = "ML")



#Lets plot OLR vs lamellae
sum_fit_phy <- summary(fit_phylo_lamellae_OLR)
PGLS_pente =   formatC(sum_fit_phy$coefficients[2], digits = 3)
PGLS_r2 =      formatC( sum_fit_phy$r.squared, digits = 2)
PGLS_pvalue =  formatC(sum_fit_phy$coefficients[8], digits = 3)
if (PGLS_pvalue == "   0"){ PGLS_pvalue = "< 2.2e-16"}
PGLS_cc <- coef(fit_phylo_lamellae_OLR)
PGLS_fit_function <- function(x) PGLS_cc[1] + PGLS_cc[2]*x
PGLS_lambda <- sum_fit_phy$param[2]
PGLS_r2


lm_fit <- lm(data = Cichlids_lamellae_brains_df %>% filter(Sp %in% b1_tree_lamellae$tip.label), 
             formula = Total_OLR ~ mean_lamellae)

sum_fit <- summary(lm_fit)
GLS_r2 =      formatC( sum_fit$r.squared, digits = 3)
GLS_pente =   formatC(sum_fit$coefficients[2], digits = 3)
GLS_pvalue =  formatC(sum_fit$coefficients[8], digits = 3)
GLS_cc <- coef(lm_fit)
GLS_fit_function <- function(x) GLS_cc[1] + GLS_cc[2]*x

pdf(file = "Raw_R_plots/OLR_vs_lamellae.pdf",width = 8.34,  height = 4.61)

Cichlids_lamellae_brains_df %>% 
  ggplot(aes( x = mean_lamellae, y= Total_OLR)) +
  geom_point(aes(color=Tribe)) +
  stat_function(fun = GLS_fit_function, color="black")+
  labs(subtitle = bquote(R^2 ~ "=" ~ .(PGLS_r2) ~ ";" ~ P ~ "=" ~.(PGLS_pvalue) ~ ";" ~ lambda ~ "=" ~ .(PGLS_lambda))) +
  theme_classic() +
  xlab("Mean number of lamellae") +
  ylab("Mean number of OLR") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none") + 
  scale_color_manual(values=tribes.colors) 

dev.off()


#Lets plot OR vs lamellae
sum_fit_phy <- summary(fit_phylo_lamellae_OR)
PGLS_pente =   formatC(sum_fit_phy$coefficients[2], digits = 3)
PGLS_r2 =      formatC( sum_fit_phy$r.squared, digits = 2)
PGLS_pvalue =  formatC(sum_fit_phy$coefficients[8], digits = 3)
if (PGLS_pvalue == "   0"){ PGLS_pvalue = "< 2.2e-16"}
PGLS_cc <- coef(fit_phylo_lamellae_OR)
PGLS_fit_function <- function(x) PGLS_cc[1] + PGLS_cc[2]*x
PGLS_lambda <- sum_fit_phy$param[2]
PGLS_r2


lm_fit <- lm(data = Cichlids_lamellae_brains_df %>% filter(Sp %in% b1_tree_lamellae$tip.label), 
             formula = Total_OR ~ mean_lamellae)

sum_fit <- summary(lm_fit)
GLS_r2 =      formatC( sum_fit$r.squared, digits = 3)
GLS_pente =   formatC(sum_fit$coefficients[2], digits = 3)
GLS_pvalue =  formatC(sum_fit$coefficients[8], digits = 3)
GLS_cc <- coef(lm_fit)
GLS_fit_function <- function(x) GLS_cc[1] + GLS_cc[2]*x

pdf(file = "Raw_R_plots/OR_vs_lamellae.pdf",width = 8.34,  height = 4.61)

Cichlids_lamellae_brains_df %>% 
  ggplot(aes( x = mean_lamellae, y= Total_OR)) +
  geom_point(aes(color=Tribe)) +
  stat_function(fun = GLS_fit_function, color="black")+
  labs(subtitle = bquote(R^2 ~ "=" ~ .(PGLS_r2) ~ ";" ~ P ~ "=" ~.(PGLS_pvalue) ~ ";" ~ lambda ~ "=" ~ .(PGLS_lambda))) +
  theme_classic() +
  xlab("Mean number of lamellae") +
  ylab("Mean number of OR") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none") + 
  scale_color_manual(values=tribes.colors) 

dev.off()


#Lets plot V2R vs lamellae
sum_fit_phy <- summary(fit_phylo_lamellae_V2R)
PGLS_pente =   formatC(sum_fit_phy$coefficients[2], digits = 3)
PGLS_r2 =      formatC( sum_fit_phy$r.squared, digits = 2)
PGLS_pvalue =  formatC(sum_fit_phy$coefficients[8], digits = 3)
if (PGLS_pvalue == "   0"){ PGLS_pvalue = "< 2.2e-16"}
PGLS_cc <- coef(fit_phylo_lamellae_V2R)
PGLS_fit_function <- function(x) PGLS_cc[1] + PGLS_cc[2]*x
PGLS_lambda <- sum_fit_phy$param[2]
PGLS_r2


lm_fit <- lm(data = Cichlids_lamellae_brains_df %>% filter(Sp %in% b1_tree_lamellae$tip.label), 
             formula = Total_V2R ~ mean_lamellae)

sum_fit <- summary(lm_fit)
GLS_r2 =      formatC( sum_fit$r.squared, digits = 3)
GLS_pente =   formatC(sum_fit$coefficients[2], digits = 3)
GLS_pvalue =  formatC(sum_fit$coefficients[8], digits = 3)
GLS_cc <- coef(lm_fit)
GLS_fit_function <- function(x) GLS_cc[1] + GLS_cc[2]*x

pdf(file = "Raw_R_plots/V2R_vs_lamellae.pdf",width = 8.34,  height = 4.61)

Cichlids_lamellae_brains_df %>% 
  ggplot(aes( x = mean_lamellae, y= Total_V2R)) +
  geom_point(aes(color=Tribe)) +
  stat_function(fun = GLS_fit_function, color="black")+
  labs(subtitle = bquote(R^2 ~ "=" ~ .(PGLS_r2) ~ ";" ~ P ~ "=" ~.(PGLS_pvalue) ~ ";" ~ lambda ~ "=" ~ .(PGLS_lambda))) +
  theme_classic() +
  xlab("Mean number of lamellae") +
  ylab("Mean number of V2R") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none") + 
  scale_color_manual(values=tribes.colors) 

dev.off()


#Lets plot TAAR vs lamellae
sum_fit_phy <- summary(fit_phylo_lamellae_TAAR)
PGLS_pente =   formatC(sum_fit_phy$coefficients[2], digits = 3)
PGLS_r2 =      formatC( sum_fit_phy$r.squared, digits = 2)
PGLS_pvalue =  formatC(sum_fit_phy$coefficients[8], digits = 3)
if (PGLS_pvalue == "   0"){ PGLS_pvalue = "< 2.2e-16"}
PGLS_cc <- coef(fit_phylo_lamellae_TAAR)
PGLS_fit_function <- function(x) PGLS_cc[1] + PGLS_cc[2]*x
PGLS_lambda <- sum_fit_phy$param[2]
PGLS_r2


lm_fit <- lm(data = Cichlids_lamellae_brains_df %>% filter(Sp %in% b1_tree_lamellae$tip.label), 
             formula = Total_TAAR ~ mean_lamellae)

sum_fit <- summary(lm_fit)
GLS_r2 =      formatC( sum_fit$r.squared, digits = 3)
GLS_pente =   formatC(sum_fit$coefficients[2], digits = 3)
GLS_pvalue =  formatC(sum_fit$coefficients[8], digits = 3)
GLS_cc <- coef(lm_fit)
GLS_fit_function <- function(x) GLS_cc[1] + GLS_cc[2]*x

pdf(file = "Raw_R_plots/TAAR_vs_lamellae.pdf",width = 8.34,  height = 4.61)

Cichlids_lamellae_brains_df %>% 
  ggplot(aes( x = mean_lamellae, y= Total_TAAR)) +
  geom_point(aes(color=Tribe)) +
  stat_function(fun = GLS_fit_function, color="black")+
  labs(subtitle = bquote(R^2 ~ "=" ~ .(PGLS_r2) ~ ";" ~ P ~ "=" ~.(PGLS_pvalue) ~ ";" ~ lambda ~ "=" ~ .(PGLS_lambda))) +
  theme_classic() +
  xlab("Mean number of lamellae") +
  ylab("Mean number of TAAR") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none") + 
  scale_color_manual(values=tribes.colors) 

dev.off()


#### Cichlids - pGLS Olfactory bulb vs Olfactory receptors  ---------------------------------


#Perform the pGLS
fit_phylo_lamellae_OLR <- pgls(Total_OLR ~ prop_Olfactory_bulb,
                               data = caper_data_lamellae_brain, 
                               lambda = "ML")
fit_phylo_lamellae_OR <- pgls(Total_OR ~ prop_Olfactory_bulb,
                              data = caper_data_lamellae_brain, 
                              lambda = "ML")
fit_phylo_lamellae_TAAR <- pgls(Total_TAAR ~ prop_Olfactory_bulb,
                                data = caper_data_lamellae_brain, 
                                lambda = "ML")
fit_phylo_lamellae_V2R <- pgls(Total_V2R ~ prop_Olfactory_bulb,
                               data = caper_data_lamellae_brain, 
                               lambda = "ML")



#Lets plot OLR vs lamellae
sum_fit_phy <- summary(fit_phylo_lamellae_OLR)
PGLS_pente =   formatC(sum_fit_phy$coefficients[2], digits = 3)
PGLS_r2 =      formatC( sum_fit_phy$r.squared, digits = 2)
PGLS_pvalue =  formatC(sum_fit_phy$coefficients[8], digits = 3)
if (PGLS_pvalue == "   0"){ PGLS_pvalue = "< 2.2e-16"}
PGLS_cc <- coef(fit_phylo_lamellae_OLR)
PGLS_fit_function <- function(x) PGLS_cc[1] + PGLS_cc[2]*x
PGLS_lambda <- sum_fit_phy$param[2]
PGLS_r2


lm_fit <- lm(data = Cichlids_lamellae_brains_df %>% filter(Sp %in% b1_tree_lamellae$tip.label), 
             formula = Total_OLR ~ prop_Olfactory_bulb)

sum_fit <- summary(lm_fit)
GLS_r2 =      formatC( sum_fit$r.squared, digits = 3)
GLS_pente =   formatC(sum_fit$coefficients[2], digits = 3)
GLS_pvalue =  formatC(sum_fit$coefficients[8], digits = 3)
GLS_cc <- coef(lm_fit)
GLS_fit_function <- function(x) GLS_cc[1] + GLS_cc[2]*x

pdf(file = "Raw_R_plots/OLR_vs_OBprop.pdf",width = 8.34,  height = 4.61)

Cichlids_lamellae_brains_df %>% 
  ggplot(aes( x = prop_Olfactory_bulb, y= Total_OLR)) +
  geom_point(aes(color=Tribe)) +
  stat_function(fun = GLS_fit_function, color="black")+
  labs(subtitle = bquote(R^2 ~ "=" ~ .(PGLS_r2) ~ ";" ~ P ~ "=" ~.(PGLS_pvalue) ~ ";" ~ lambda ~ "=" ~ .(PGLS_lambda))) +
  theme_classic() +
  xlab("Mean relative olfactory bulb volume") +
  ylab("Mean number of OLR") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none") + 
  scale_color_manual(values=tribes.colors) 

dev.off()


#Lets plot OR vs lamellae
sum_fit_phy <- summary(fit_phylo_lamellae_OR)
PGLS_pente =   formatC(sum_fit_phy$coefficients[2], digits = 3)
PGLS_r2 =      formatC( sum_fit_phy$r.squared, digits = 2)
PGLS_pvalue =  formatC(sum_fit_phy$coefficients[8], digits = 3)
if (PGLS_pvalue == "   0"){ PGLS_pvalue = "< 2.2e-16"}
PGLS_cc <- coef(fit_phylo_lamellae_OR)
PGLS_fit_function <- function(x) PGLS_cc[1] + PGLS_cc[2]*x
PGLS_lambda <- sum_fit_phy$param[2]
PGLS_r2


lm_fit <- lm(data = Cichlids_lamellae_brains_df %>% filter(Sp %in% b1_tree_lamellae$tip.label), 
             formula = Total_OR ~ prop_Olfactory_bulb)

sum_fit <- summary(lm_fit)
GLS_r2 =      formatC( sum_fit$r.squared, digits = 3)
GLS_pente =   formatC(sum_fit$coefficients[2], digits = 3)
GLS_pvalue =  formatC(sum_fit$coefficients[8], digits = 3)
GLS_cc <- coef(lm_fit)
GLS_fit_function <- function(x) GLS_cc[1] + GLS_cc[2]*x

pdf(file = "Raw_R_plots/OR_vs_OBprop.pdf",width = 8.34,  height = 4.61)

Cichlids_lamellae_brains_df %>% 
  ggplot(aes( x = prop_Olfactory_bulb, y= Total_OR)) +
  geom_point(aes(color=Tribe)) +
  stat_function(fun = GLS_fit_function, color="black")+
  labs(subtitle = bquote(R^2 ~ "=" ~ .(PGLS_r2) ~ ";" ~ P ~ "=" ~.(PGLS_pvalue) ~ ";" ~ lambda ~ "=" ~ .(PGLS_lambda))) +
  theme_classic() +
  xlab("Mean relative olfactory bulb volume") +
  ylab("Mean number of OR") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none") + 
  scale_color_manual(values=tribes.colors) 

dev.off()


#Lets plot V2R vs lamellae
sum_fit_phy <- summary(fit_phylo_lamellae_V2R)
PGLS_pente =   formatC(sum_fit_phy$coefficients[2], digits = 3)
PGLS_r2 =      formatC( sum_fit_phy$r.squared, digits = 2)
PGLS_pvalue =  formatC(sum_fit_phy$coefficients[8], digits = 3)
if (PGLS_pvalue == "   0"){ PGLS_pvalue = "< 2.2e-16"}
PGLS_cc <- coef(fit_phylo_lamellae_V2R)
PGLS_fit_function <- function(x) PGLS_cc[1] + PGLS_cc[2]*x
PGLS_lambda <- sum_fit_phy$param[2]
PGLS_r2


lm_fit <- lm(data = Cichlids_lamellae_brains_df %>% filter(Sp %in% b1_tree_lamellae$tip.label), 
             formula = Total_V2R ~ prop_Olfactory_bulb)

sum_fit <- summary(lm_fit)
GLS_r2 =      formatC( sum_fit$r.squared, digits = 3)
GLS_pente =   formatC(sum_fit$coefficients[2], digits = 3)
GLS_pvalue =  formatC(sum_fit$coefficients[8], digits = 3)
GLS_cc <- coef(lm_fit)
GLS_fit_function <- function(x) GLS_cc[1] + GLS_cc[2]*x

pdf(file = "Raw_R_plots/V2R_vs_OBprop.pdf",width = 8.34,  height = 4.61)

Cichlids_lamellae_brains_df %>% 
  ggplot(aes( x = prop_Olfactory_bulb, y= Total_V2R)) +
  geom_point(aes(color=Tribe)) +
  stat_function(fun = GLS_fit_function, color="black")+
  labs(subtitle = bquote(R^2 ~ "=" ~ .(PGLS_r2) ~ ";" ~ P ~ "=" ~.(PGLS_pvalue) ~ ";" ~ lambda ~ "=" ~ .(PGLS_lambda))) +
  theme_classic() +
  xlab("Mean relative olfactory bulb volume") +
  ylab("Mean number of V2R") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none") + 
  scale_color_manual(values=tribes.colors) 

dev.off()


#Lets plot TAAR vs lamellae
sum_fit_phy <- summary(fit_phylo_lamellae_TAAR)
PGLS_pente =   formatC(sum_fit_phy$coefficients[2], digits = 3)
PGLS_r2 =      formatC( sum_fit_phy$r.squared, digits = 2)
PGLS_pvalue =  formatC(sum_fit_phy$coefficients[8], digits = 3)
if (PGLS_pvalue == "   0"){ PGLS_pvalue = "< 2.2e-16"}
PGLS_cc <- coef(fit_phylo_lamellae_TAAR)
PGLS_fit_function <- function(x) PGLS_cc[1] + PGLS_cc[2]*x
PGLS_lambda <- sum_fit_phy$param[2]
PGLS_r2


lm_fit <- lm(data = Cichlids_lamellae_brains_df %>% filter(Sp %in% b1_tree_lamellae$tip.label), 
             formula = Total_TAAR ~ prop_Olfactory_bulb)

sum_fit <- summary(lm_fit)
GLS_r2 =      formatC( sum_fit$r.squared, digits = 3)
GLS_pente =   formatC(sum_fit$coefficients[2], digits = 3)
GLS_pvalue =  formatC(sum_fit$coefficients[8], digits = 3)
GLS_cc <- coef(lm_fit)
GLS_fit_function <- function(x) GLS_cc[1] + GLS_cc[2]*x

pdf(file = "Raw_R_plots/TAAR_vs_OBprop.pdf",width = 8.34,  height = 4.61)

Cichlids_lamellae_brains_df %>% 
  ggplot(aes( x = prop_Olfactory_bulb, y= Total_TAAR)) +
  geom_point(aes(color=Tribe)) +
  stat_function(fun = GLS_fit_function, color="black")+
  labs(subtitle = bquote(R^2 ~ "=" ~ .(PGLS_r2) ~ ";" ~ P ~ "=" ~.(PGLS_pvalue) ~ ";" ~ lambda ~ "=" ~ .(PGLS_lambda))) +
  theme_classic() +
  xlab("Mean relative olfactory bulb volume") +
  ylab("Mean number of TAAR") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none") + 
  scale_color_manual(values=tribes.colors) 

dev.off()


#### Cichlids - pGLS telencephalon vs Olfactory receptors  ---------------------------------


#Perform the pGLS
fit_phylo_Telencephalon_OLR <- pgls(Total_OLR ~ prop_Telencephalon,
                                    data = caper_data_lamellae_brain, 
                                    lambda = "ML", bounds=list(lambda=c(0,1)))
fit_phylo_Telencephalon_OR <- pgls(Total_OR ~ prop_Telencephalon,
                                   data = caper_data_lamellae_brain, 
                                   lambda = "ML")
fit_phylo_Telencephalon_TAAR <- pgls(Total_TAAR ~ prop_Telencephalon,
                                     data = caper_data_lamellae_brain, 
                                     lambda = "ML")
fit_phylo_Telencephalon_V2R <- pgls(Total_V2R ~ prop_Telencephalon,
                                    data = caper_data_lamellae_brain, 
                                    lambda = "ML")



#Lets plot OLR vs Telencephalon
sum_fit_phy <- summary(fit_phylo_Telencephalon_OLR)
PGLS_pente =   formatC(sum_fit_phy$coefficients[2], digits = 3)
PGLS_r2 =      formatC( sum_fit_phy$r.squared, digits = 2)
PGLS_pvalue =  formatC(sum_fit_phy$coefficients[8], digits = 3)
if (PGLS_pvalue == "   0"){ PGLS_pvalue = "< 2.2e-16"}
PGLS_cc <- coef(fit_phylo_Telencephalon_OLR)
PGLS_fit_function <- function(x) PGLS_cc[1] + PGLS_cc[2]*x
PGLS_lambda <- sum_fit_phy$param[2]
PGLS_r2


lm_fit <- lm(data = Cichlids_lamellae_brains_df %>% filter(Sp %in% b1_tree_lamellae$tip.label), 
             formula = Total_OLR ~ prop_Telencephalon)

sum_fit <- summary(lm_fit)
GLS_r2 =      formatC( sum_fit$r.squared, digits = 3)
GLS_pente =   formatC(sum_fit$coefficients[2], digits = 3)
GLS_pvalue =  formatC(sum_fit$coefficients[8], digits = 3)
GLS_cc <- coef(lm_fit)
GLS_fit_function <- function(x) GLS_cc[1] + GLS_cc[2]*x

pdf(file = "Raw_R_plots/OLR_vs_Telencephalonprop.pdf",width = 8.34,  height = 4.61)

Cichlids_lamellae_brains_df %>% 
  ggplot(aes( x = prop_Telencephalon, y= Total_OLR)) +
  geom_point(aes(color=Tribe)) +
  stat_function(fun = GLS_fit_function, color="black")+
  labs(subtitle = bquote(R^2 ~ "=" ~ .(PGLS_r2) ~ ";" ~ P ~ "=" ~.(PGLS_pvalue) ~ ";" ~ lambda ~ "=" ~ .(PGLS_lambda))) +
  theme_classic() +
  xlab("Mean relative telencephalon volume") +
  ylab("Mean number of OLR") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none") + 
  scale_color_manual(values=tribes.colors) 

dev.off()


#Lets plot OR vs Telencephalon
sum_fit_phy <- summary(fit_phylo_Telencephalon_OR)
PGLS_pente =   formatC(sum_fit_phy$coefficients[2], digits = 3)
PGLS_r2 =      formatC( sum_fit_phy$r.squared, digits = 2)
PGLS_pvalue =  formatC(sum_fit_phy$coefficients[8], digits = 3)
if (PGLS_pvalue == "   0"){ PGLS_pvalue = "< 2.2e-16"}
PGLS_cc <- coef(fit_phylo_Telencephalon_OR)
PGLS_fit_function <- function(x) PGLS_cc[1] + PGLS_cc[2]*x
PGLS_lambda <- sum_fit_phy$param[2]
PGLS_r2


lm_fit <- lm(data = Cichlids_lamellae_brains_df %>% filter(Sp %in% b1_tree_lamellae$tip.label), 
             formula = Total_OR ~ prop_Telencephalon)

sum_fit <- summary(lm_fit)
GLS_r2 =      formatC( sum_fit$r.squared, digits = 3)
GLS_pente =   formatC(sum_fit$coefficients[2], digits = 3)
GLS_pvalue =  formatC(sum_fit$coefficients[8], digits = 3)
GLS_cc <- coef(lm_fit)
GLS_fit_function <- function(x) GLS_cc[1] + GLS_cc[2]*x

pdf(file = "Raw_R_plots/OR_vs_Telencephalonprop.pdf",width = 8.34,  height = 4.61)

Cichlids_lamellae_brains_df %>% 
  ggplot(aes( x = prop_Telencephalon, y= Total_OR)) +
  geom_point(aes(color=Tribe)) +
  stat_function(fun = GLS_fit_function, color="black")+
  labs(subtitle = bquote(R^2 ~ "=" ~ .(PGLS_r2) ~ ";" ~ P ~ "=" ~.(PGLS_pvalue) ~ ";" ~ lambda ~ "=" ~ .(PGLS_lambda))) +
  theme_classic() +
  xlab("Mean relative telencephalon volume") +
  ylab("Mean number of OR") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none") + 
  scale_color_manual(values=tribes.colors) 

dev.off()


#Lets plot V2R vs Telencephalon
sum_fit_phy <- summary(fit_phylo_Telencephalon_V2R)
PGLS_pente =   formatC(sum_fit_phy$coefficients[2], digits = 3)
PGLS_r2 =      formatC( sum_fit_phy$r.squared, digits = 2)
PGLS_pvalue =  formatC(sum_fit_phy$coefficients[8], digits = 3)
if (PGLS_pvalue == "   0"){ PGLS_pvalue = "< 2.2e-16"}
PGLS_cc <- coef(fit_phylo_Telencephalon_V2R)
PGLS_fit_function <- function(x) PGLS_cc[1] + PGLS_cc[2]*x
PGLS_lambda <- sum_fit_phy$param[2]
PGLS_r2


lm_fit <- lm(data = Cichlids_lamellae_brains_df %>% filter(Sp %in% b1_tree_lamellae$tip.label), 
             formula = Total_V2R ~ prop_Telencephalon)

sum_fit <- summary(lm_fit)
GLS_r2 =      formatC( sum_fit$r.squared, digits = 3)
GLS_pente =   formatC(sum_fit$coefficients[2], digits = 3)
GLS_pvalue =  formatC(sum_fit$coefficients[8], digits = 3)
GLS_cc <- coef(lm_fit)
GLS_fit_function <- function(x) GLS_cc[1] + GLS_cc[2]*x

pdf(file = "Raw_R_plots/V2R_vs_Telencephalonprop.pdf",width = 8.34,  height = 4.61)

Cichlids_lamellae_brains_df %>% 
  ggplot(aes( x = prop_Telencephalon, y= Total_V2R)) +
  geom_point(aes(color=Tribe)) +
  stat_function(fun = GLS_fit_function, color="black")+
  labs(subtitle = bquote(R^2 ~ "=" ~ .(PGLS_r2) ~ ";" ~ P ~ "=" ~.(PGLS_pvalue) ~ ";" ~ lambda ~ "=" ~ .(PGLS_lambda))) +
  theme_classic() +
  xlab("Mean relative telencephalon volume") +
  ylab("Mean number of V2R") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none") + 
  scale_color_manual(values=tribes.colors) 

dev.off()


#Lets plot TAAR vs Telencephalon
sum_fit_phy <- summary(fit_phylo_Telencephalon_TAAR)
PGLS_pente =   formatC(sum_fit_phy$coefficients[2], digits = 3)
PGLS_r2 =      formatC( sum_fit_phy$r.squared, digits = 2)
PGLS_pvalue =  formatC(sum_fit_phy$coefficients[8], digits = 3)
if (PGLS_pvalue == "   0"){ PGLS_pvalue = "< 2.2e-16"}
PGLS_cc <- coef(fit_phylo_Telencephalon_TAAR)
PGLS_fit_function <- function(x) PGLS_cc[1] + PGLS_cc[2]*x
PGLS_lambda <- sum_fit_phy$param[2]
PGLS_r2


lm_fit <- lm(data = Cichlids_lamellae_brains_df %>% filter(Sp %in% b1_tree_lamellae$tip.label), 
             formula = Total_TAAR ~ prop_Telencephalon)

sum_fit <- summary(lm_fit)
GLS_r2 =      formatC( sum_fit$r.squared, digits = 3)
GLS_pente =   formatC(sum_fit$coefficients[2], digits = 3)
GLS_pvalue =  formatC(sum_fit$coefficients[8], digits = 3)
GLS_cc <- coef(lm_fit)
GLS_fit_function <- function(x) GLS_cc[1] + GLS_cc[2]*x

pdf(file = "Raw_R_plots/TAAR_vs_Telencephalonprop.pdf",width = 8.34,  height = 4.61)

Cichlids_lamellae_brains_df %>% 
  ggplot(aes( x = prop_Telencephalon, y= Total_TAAR)) +
  geom_point(aes(color=Tribe, size=2)) +
  stat_function(fun = GLS_fit_function, color="black")+
  labs(subtitle = bquote(R^2 ~ "=" ~ .(PGLS_r2) ~ ";" ~ P ~ "=" ~.(PGLS_pvalue) ~ ";" ~ lambda ~ "=" ~ .(PGLS_lambda))) +
  theme_classic() +
  xlab("Mean relative telencephalon volume") +
  ylab("Mean number of TAAR") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none") + 
  scale_color_manual(values=tribes.colors) 

dev.off()





#### Cichlids - pGLS Optic_tectum vs Olfactory receptors  ---------------------------------


#Perform the pGLS
fit_phylo_Optic_tectum_OLR <- pgls(Total_OLR ~ prop_Optic_tectum,
                                   data = caper_data_lamellae_brain, 
                                   lambda = "ML", bounds=list(lambda=c(0,1)))
fit_phylo_Optic_tectum_OR <- pgls(Total_OR ~ prop_Optic_tectum,
                                  data = caper_data_lamellae_brain, 
                                  lambda = "ML")
fit_phylo_Optic_tectum_TAAR <- pgls(Total_TAAR ~ prop_Optic_tectum,
                                    data = caper_data_lamellae_brain, 
                                    lambda = "ML")
fit_phylo_Optic_tectum_V2R <- pgls(Total_V2R ~ prop_Optic_tectum,
                                   data = caper_data_lamellae_brain, 
                                   lambda = "ML")



#Lets plot OLR vs Optic_tectum
sum_fit_phy <- summary(fit_phylo_Optic_tectum_OLR)
PGLS_pente =   formatC(sum_fit_phy$coefficients[2], digits = 3)
PGLS_r2 =      formatC( sum_fit_phy$r.squared, digits = 2)
PGLS_pvalue =  formatC(sum_fit_phy$coefficients[8], digits = 3)
if (PGLS_pvalue == "   0"){ PGLS_pvalue = "< 2.2e-16"}
PGLS_cc <- coef(fit_phylo_Optic_tectum_OLR)
PGLS_fit_function <- function(x) PGLS_cc[1] + PGLS_cc[2]*x
PGLS_lambda <- sum_fit_phy$param[2]
PGLS_r2


lm_fit <- lm(data = Cichlids_lamellae_brains_df %>% filter(Sp %in% b1_tree_lamellae$tip.label), 
             formula = Total_OLR ~ prop_Optic_tectum)

sum_fit <- summary(lm_fit)
GLS_r2 =      formatC( sum_fit$r.squared, digits = 3)
GLS_pente =   formatC(sum_fit$coefficients[2], digits = 3)
GLS_pvalue =  formatC(sum_fit$coefficients[8], digits = 3)
GLS_cc <- coef(lm_fit)
GLS_fit_function <- function(x) GLS_cc[1] + GLS_cc[2]*x

pdf(file = "Raw_R_plots/OLR_vs_Optic_tectumprop.pdf",width = 8.34,  height = 4.61)

Cichlids_lamellae_brains_df %>% 
  ggplot(aes( x = prop_Optic_tectum, y= Total_OLR)) +
  geom_point(aes(color=Tribe)) +
  stat_function(fun = GLS_fit_function, color="black")+
  labs(subtitle = bquote(R^2 ~ "=" ~ .(PGLS_r2) ~ ";" ~ P ~ "=" ~.(PGLS_pvalue) ~ ";" ~ lambda ~ "=" ~ .(PGLS_lambda))) +
  theme_classic() +
  xlab("Mean relative optic tectum volume") +
  ylab("Mean number of OLR") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none") + 
  scale_color_manual(values=tribes.colors) 

dev.off()


#Lets plot OR vs Optic_tectum
sum_fit_phy <- summary(fit_phylo_Optic_tectum_OR)
PGLS_pente =   formatC(sum_fit_phy$coefficients[2], digits = 3)
PGLS_r2 =      formatC( sum_fit_phy$r.squared, digits = 2)
PGLS_pvalue =  formatC(sum_fit_phy$coefficients[8], digits = 3)
if (PGLS_pvalue == "   0"){ PGLS_pvalue = "< 2.2e-16"}
PGLS_cc <- coef(fit_phylo_Optic_tectum_OR)
PGLS_fit_function <- function(x) PGLS_cc[1] + PGLS_cc[2]*x
PGLS_lambda <- sum_fit_phy$param[2]
PGLS_r2


lm_fit <- lm(data = Cichlids_lamellae_brains_df %>% filter(Sp %in% b1_tree_lamellae$tip.label), 
             formula = Total_OR ~ prop_Optic_tectum)

sum_fit <- summary(lm_fit)
GLS_r2 =      formatC( sum_fit$r.squared, digits = 3)
GLS_pente =   formatC(sum_fit$coefficients[2], digits = 3)
GLS_pvalue =  formatC(sum_fit$coefficients[8], digits = 3)
GLS_cc <- coef(lm_fit)
GLS_fit_function <- function(x) GLS_cc[1] + GLS_cc[2]*x

pdf(file = "Raw_R_plots/OR_vs_Optic_tectumprop.pdf",width = 8.34,  height = 4.61)

Cichlids_lamellae_brains_df %>% 
  ggplot(aes( x = prop_Optic_tectum, y= Total_OR)) +
  geom_point(aes(color=Tribe)) +
  stat_function(fun = GLS_fit_function, color="black")+
  labs(subtitle = bquote(R^2 ~ "=" ~ .(PGLS_r2) ~ ";" ~ P ~ "=" ~.(PGLS_pvalue) ~ ";" ~ lambda ~ "=" ~ .(PGLS_lambda))) +
  theme_classic() +
  xlab("Mean relative optic tectum volume") +
  ylab("Mean number of OR") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none") + 
  scale_color_manual(values=tribes.colors) 

dev.off()


#Lets plot V2R vs Optic_tectum
sum_fit_phy <- summary(fit_phylo_Optic_tectum_V2R)
PGLS_pente =   formatC(sum_fit_phy$coefficients[2], digits = 3)
PGLS_r2 =      formatC( sum_fit_phy$r.squared, digits = 2)
PGLS_pvalue =  formatC(sum_fit_phy$coefficients[8], digits = 3)
if (PGLS_pvalue == "   0"){ PGLS_pvalue = "< 2.2e-16"}
PGLS_cc <- coef(fit_phylo_Optic_tectum_V2R)
PGLS_fit_function <- function(x) PGLS_cc[1] + PGLS_cc[2]*x
PGLS_lambda <- sum_fit_phy$param[2]
PGLS_r2


lm_fit <- lm(data = Cichlids_lamellae_brains_df %>% filter(Sp %in% b1_tree_lamellae$tip.label), 
             formula = Total_V2R ~ prop_Optic_tectum)

sum_fit <- summary(lm_fit)
GLS_r2 =      formatC( sum_fit$r.squared, digits = 3)
GLS_pente =   formatC(sum_fit$coefficients[2], digits = 3)
GLS_pvalue =  formatC(sum_fit$coefficients[8], digits = 3)
GLS_cc <- coef(lm_fit)
GLS_fit_function <- function(x) GLS_cc[1] + GLS_cc[2]*x

pdf(file = "Raw_R_plots/V2R_vs_Optic_tectumprop.pdf",width = 8.34,  height = 4.61)

Cichlids_lamellae_brains_df %>% 
  ggplot(aes( x = prop_Optic_tectum, y= Total_V2R)) +
  geom_point(aes(color=Tribe)) +
  stat_function(fun = GLS_fit_function, color="black")+
  labs(subtitle = bquote(R^2 ~ "=" ~ .(PGLS_r2) ~ ";" ~ P ~ "=" ~.(PGLS_pvalue) ~ ";" ~ lambda ~ "=" ~ .(PGLS_lambda))) +
  theme_classic() +
  xlab("Mean relative optic tectum volume") +
  ylab("Mean number of V2R") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none") + 
  scale_color_manual(values=tribes.colors) 

dev.off()


#Lets plot TAAR vs Optic_tectum
sum_fit_phy <- summary(fit_phylo_Optic_tectum_TAAR)
PGLS_pente =   formatC(sum_fit_phy$coefficients[2], digits = 3)
PGLS_r2 =      formatC( sum_fit_phy$r.squared, digits = 2)
PGLS_pvalue =  formatC(sum_fit_phy$coefficients[8], digits = 3)
if (PGLS_pvalue == "   0"){ PGLS_pvalue = "< 2.2e-16"}
PGLS_cc <- coef(fit_phylo_Optic_tectum_TAAR)
PGLS_fit_function <- function(x) PGLS_cc[1] + PGLS_cc[2]*x
PGLS_lambda <- sum_fit_phy$param[2]
PGLS_r2


lm_fit <- lm(data = Cichlids_lamellae_brains_df %>% filter(Sp %in% b1_tree_lamellae$tip.label), 
             formula = Total_TAAR ~ prop_Optic_tectum)

sum_fit <- summary(lm_fit)
GLS_r2 =      formatC( sum_fit$r.squared, digits = 3)
GLS_pente =   formatC(sum_fit$coefficients[2], digits = 3)
GLS_pvalue =  formatC(sum_fit$coefficients[8], digits = 3)
GLS_cc <- coef(lm_fit)
GLS_fit_function <- function(x) GLS_cc[1] + GLS_cc[2]*x

pdf(file = "Raw_R_plots/TAAR_vs_Optic_tectumprop.pdf",width = 8.34,  height = 4.61)

Cichlids_lamellae_brains_df %>% 
  ggplot(aes( x = prop_Optic_tectum, y= Total_TAAR)) +
  geom_point(aes(color=Tribe)) +
  stat_function(fun = GLS_fit_function, color="black")+
  labs(subtitle = bquote(R^2 ~ "=" ~ .(PGLS_r2) ~ ";" ~ P ~ "=" ~.(PGLS_pvalue) ~ ";" ~ lambda ~ "=" ~ .(PGLS_lambda))) +
  theme_classic() +
  xlab("Mean relative optic tectum volume") +
  ylab("Mean number of TAAR") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none") + 
  scale_color_manual(values=tribes.colors) 

dev.off()





#### Cichlids - pGLS Rest_of_the_brain vs Olfactory receptors  ---------------------------------


#Perform the pGLS
fit_phylo_Rest_of_the_brain_OLR <- pgls(Total_OLR ~ prop_Rest_of_the_brain,
                                        data = caper_data_lamellae_brain, 
                                        lambda = "ML", bounds=list(lambda=c(0,1)))
fit_phylo_Rest_of_the_brain_OR <- pgls(Total_OR ~ prop_Rest_of_the_brain,
                                       data = caper_data_lamellae_brain, 
                                       lambda = "ML", bounds=list(lambda=c(0,1)))
fit_phylo_Rest_of_the_brain_TAAR <- pgls(Total_TAAR ~ prop_Rest_of_the_brain,
                                         data = caper_data_lamellae_brain, 
                                         lambda = "ML")
fit_phylo_Rest_of_the_brain_V2R <- pgls(Total_V2R ~ prop_Rest_of_the_brain,
                                        data = caper_data_lamellae_brain, 
                                        lambda = "ML")



#Lets plot OLR vs Rest_of_the_brain
sum_fit_phy <- summary(fit_phylo_Rest_of_the_brain_OLR)
PGLS_pente =   formatC(sum_fit_phy$coefficients[2], digits = 3)
PGLS_r2 =      formatC( sum_fit_phy$r.squared, digits = 2)
PGLS_pvalue =  formatC(sum_fit_phy$coefficients[8], digits = 3)
if (PGLS_pvalue == "   0"){ PGLS_pvalue = "< 2.2e-16"}
PGLS_cc <- coef(fit_phylo_Rest_of_the_brain_OLR)
PGLS_fit_function <- function(x) PGLS_cc[1] + PGLS_cc[2]*x
PGLS_lambda <- sum_fit_phy$param[2]
PGLS_r2


lm_fit <- lm(data = Cichlids_lamellae_brains_df %>% filter(Sp %in% b1_tree_lamellae$tip.label), 
             formula = Total_OLR ~ prop_Rest_of_the_brain)

sum_fit <- summary(lm_fit)
GLS_r2 =      formatC( sum_fit$r.squared, digits = 3)
GLS_pente =   formatC(sum_fit$coefficients[2], digits = 3)
GLS_pvalue =  formatC(sum_fit$coefficients[8], digits = 3)
GLS_cc <- coef(lm_fit)
GLS_fit_function <- function(x) GLS_cc[1] + GLS_cc[2]*x

pdf(file = "Raw_R_plots/OLR_vs_Rest_of_the_brainprop.pdf",width = 8.34,  height = 4.61)

Cichlids_lamellae_brains_df %>% 
  ggplot(aes( x = prop_Rest_of_the_brain, y= Total_OLR)) +
  geom_point(aes(color=Tribe)) +
  stat_function(fun = GLS_fit_function, color="black")+
  labs(subtitle = bquote(R^2 ~ "=" ~ .(PGLS_r2) ~ ";" ~ P ~ "=" ~.(PGLS_pvalue) ~ ";" ~ lambda ~ "=" ~ .(PGLS_lambda))) +
  theme_classic() +
  xlab("Mean relative rest of the brain volume") +
  ylab("Mean number of OLR") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none") + 
  scale_color_manual(values=tribes.colors) 

dev.off()


#Lets plot OR vs Rest_of_the_brain
sum_fit_phy <- summary(fit_phylo_Rest_of_the_brain_OR)
PGLS_pente =   formatC(sum_fit_phy$coefficients[2], digits = 3)
PGLS_r2 =      formatC( sum_fit_phy$r.squared, digits = 2)
PGLS_pvalue =  formatC(sum_fit_phy$coefficients[8], digits = 3)
if (PGLS_pvalue == "   0"){ PGLS_pvalue = "< 2.2e-16"}
PGLS_cc <- coef(fit_phylo_Rest_of_the_brain_OR)
PGLS_fit_function <- function(x) PGLS_cc[1] + PGLS_cc[2]*x
PGLS_lambda <- sum_fit_phy$param[2]
PGLS_r2


lm_fit <- lm(data = Cichlids_lamellae_brains_df %>% filter(Sp %in% b1_tree_lamellae$tip.label), 
             formula = Total_OR ~ prop_Rest_of_the_brain)

sum_fit <- summary(lm_fit)
GLS_r2 =      formatC( sum_fit$r.squared, digits = 3)
GLS_pente =   formatC(sum_fit$coefficients[2], digits = 3)
GLS_pvalue =  formatC(sum_fit$coefficients[8], digits = 3)
GLS_cc <- coef(lm_fit)
GLS_fit_function <- function(x) GLS_cc[1] + GLS_cc[2]*x

pdf(file = "Raw_R_plots/OR_vs_Rest_of_the_brainprop.pdf",width = 8.34,  height = 4.61)

Cichlids_lamellae_brains_df %>% 
  ggplot(aes( x = prop_Rest_of_the_brain, y= Total_OR)) +
  geom_point(aes(color=Tribe)) +
  stat_function(fun = GLS_fit_function, color="black")+
  labs(subtitle = bquote(R^2 ~ "=" ~ .(PGLS_r2) ~ ";" ~ P ~ "=" ~.(PGLS_pvalue) ~ ";" ~ lambda ~ "=" ~ .(PGLS_lambda))) +
  theme_classic() +
  xlab("Mean relative rest of the brain volume") +
  ylab("Mean number of OR") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none") + 
  scale_color_manual(values=tribes.colors) 

dev.off()


#Lets plot V2R vs Rest_of_the_brain
sum_fit_phy <- summary(fit_phylo_Rest_of_the_brain_V2R)
PGLS_pente =   formatC(sum_fit_phy$coefficients[2], digits = 3)
PGLS_r2 =      formatC( sum_fit_phy$r.squared, digits = 2)
PGLS_pvalue =  formatC(sum_fit_phy$coefficients[8], digits = 3)
if (PGLS_pvalue == "   0"){ PGLS_pvalue = "< 2.2e-16"}
PGLS_cc <- coef(fit_phylo_Rest_of_the_brain_V2R)
PGLS_fit_function <- function(x) PGLS_cc[1] + PGLS_cc[2]*x
PGLS_lambda <- sum_fit_phy$param[2]
PGLS_r2


lm_fit <- lm(data = Cichlids_lamellae_brains_df %>% filter(Sp %in% b1_tree_lamellae$tip.label), 
             formula = Total_V2R ~ prop_Rest_of_the_brain)

sum_fit <- summary(lm_fit)
GLS_r2 =      formatC( sum_fit$r.squared, digits = 3)
GLS_pente =   formatC(sum_fit$coefficients[2], digits = 3)
GLS_pvalue =  formatC(sum_fit$coefficients[8], digits = 3)
GLS_cc <- coef(lm_fit)
GLS_fit_function <- function(x) GLS_cc[1] + GLS_cc[2]*x

pdf(file = "Raw_R_plots/V2R_vs_Rest_of_the_brainprop.pdf",width = 8.34,  height = 4.61)

Cichlids_lamellae_brains_df %>% 
  ggplot(aes( x = prop_Rest_of_the_brain, y= Total_V2R)) +
  geom_point(aes(color=Tribe)) +
  stat_function(fun = GLS_fit_function, color="black")+
  labs(subtitle = bquote(R^2 ~ "=" ~ .(PGLS_r2) ~ ";" ~ P ~ "=" ~.(PGLS_pvalue) ~ ";" ~ lambda ~ "=" ~ .(PGLS_lambda))) +
  theme_classic() +
  xlab("Mean relative rest of the brain volume") +
  ylab("Mean number of V2R") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none") + 
  scale_color_manual(values=tribes.colors) 

dev.off()


#Lets plot TAAR vs Rest_of_the_brain
sum_fit_phy <- summary(fit_phylo_Rest_of_the_brain_TAAR)
PGLS_pente =   formatC(sum_fit_phy$coefficients[2], digits = 3)
PGLS_r2 =      formatC( sum_fit_phy$r.squared, digits = 2)
PGLS_pvalue =  formatC(sum_fit_phy$coefficients[8], digits = 3)
if (PGLS_pvalue == "   0"){ PGLS_pvalue = "< 2.2e-16"}
PGLS_cc <- coef(fit_phylo_Rest_of_the_brain_TAAR)
PGLS_fit_function <- function(x) PGLS_cc[1] + PGLS_cc[2]*x
PGLS_lambda <- sum_fit_phy$param[2]
PGLS_r2


lm_fit <- lm(data = Cichlids_lamellae_brains_df %>% filter(Sp %in% b1_tree_lamellae$tip.label), 
             formula = Total_TAAR ~ prop_Rest_of_the_brain)

sum_fit <- summary(lm_fit)
GLS_r2 =      formatC( sum_fit$r.squared, digits = 3)
GLS_pente =   formatC(sum_fit$coefficients[2], digits = 3)
GLS_pvalue =  formatC(sum_fit$coefficients[8], digits = 3)
GLS_cc <- coef(lm_fit)
GLS_fit_function <- function(x) GLS_cc[1] + GLS_cc[2]*x

pdf(file = "Raw_R_plots/TAAR_vs_Rest_of_the_brainprop.pdf",width = 8.34,  height = 4.61)

Cichlids_lamellae_brains_df %>% 
  ggplot(aes( x = prop_Rest_of_the_brain, y= Total_TAAR)) +
  geom_point(aes(color=Tribe)) +
  stat_function(fun = GLS_fit_function, color="black")+
  labs(subtitle = bquote(R^2 ~ "=" ~ .(PGLS_r2) ~ ";" ~ P ~ "=" ~.(PGLS_pvalue) ~ ";" ~ lambda ~ "=" ~ .(PGLS_lambda))) +
  theme_classic() +
  xlab("Mean relative rest of the brain volume") +
  ylab("Mean number of TAAR") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none") + 
  scale_color_manual(values=tribes.colors) 

dev.off()





#### Cichlids - pGLS Cerebellum vs Olfactory receptors  ---------------------------------


#Perform the pGLS
fit_phylo_Cerebellum_OLR <- pgls(Total_OLR ~ prop_Cerebellum,
                                 data = caper_data_lamellae_brain, 
                                 lambda = "ML", bounds=list(lambda=c(0,1)))
fit_phylo_Cerebellum_OR <- pgls(Total_OR ~ prop_Cerebellum,
                                data = caper_data_lamellae_brain, 
                                lambda = "ML")
fit_phylo_Cerebellum_TAAR <- pgls(Total_TAAR ~ prop_Cerebellum,
                                  data = caper_data_lamellae_brain, 
                                  lambda = "ML")
fit_phylo_Cerebellum_V2R <- pgls(Total_V2R ~ prop_Cerebellum,
                                 data = caper_data_lamellae_brain, 
                                 lambda = "ML")



#Lets plot OLR vs Cerebellum
sum_fit_phy <- summary(fit_phylo_Cerebellum_OLR)
PGLS_pente =   formatC(sum_fit_phy$coefficients[2], digits = 3)
PGLS_r2 =      formatC( sum_fit_phy$r.squared, digits = 2)
PGLS_pvalue =  formatC(sum_fit_phy$coefficients[8], digits = 3)
if (PGLS_pvalue == "   0"){ PGLS_pvalue = "< 2.2e-16"}
PGLS_cc <- coef(fit_phylo_Cerebellum_OLR)
PGLS_fit_function <- function(x) PGLS_cc[1] + PGLS_cc[2]*x
PGLS_lambda <- sum_fit_phy$param[2]
PGLS_r2


lm_fit <- lm(data = Cichlids_lamellae_brains_df %>% filter(Sp %in% b1_tree_lamellae$tip.label), 
             formula = Total_OLR ~ prop_Cerebellum)

sum_fit <- summary(lm_fit)
GLS_r2 =      formatC( sum_fit$r.squared, digits = 3)
GLS_pente =   formatC(sum_fit$coefficients[2], digits = 3)
GLS_pvalue =  formatC(sum_fit$coefficients[8], digits = 3)
GLS_cc <- coef(lm_fit)
GLS_fit_function <- function(x) GLS_cc[1] + GLS_cc[2]*x

pdf(file = "Raw_R_plots/OLR_vs_Cerebellumprop.pdf",width = 8.34,  height = 4.61)

Cichlids_lamellae_brains_df %>% 
  ggplot(aes( x = prop_Cerebellum, y= Total_OLR)) +
  geom_point(aes(color=Tribe)) +
  stat_function(fun = GLS_fit_function, color="black")+
  labs(subtitle = bquote(R^2 ~ "=" ~ .(PGLS_r2) ~ ";" ~ P ~ "=" ~.(PGLS_pvalue) ~ ";" ~ lambda ~ "=" ~ .(PGLS_lambda))) +
  theme_classic() +
  xlab("Mean relative cerebellum volume") +
  ylab("Mean number of OLR") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none") + 
  scale_color_manual(values=tribes.colors) 

dev.off()


#Lets plot OR vs Cerebellum
sum_fit_phy <- summary(fit_phylo_Cerebellum_OR)
PGLS_pente =   formatC(sum_fit_phy$coefficients[2], digits = 3)
PGLS_r2 =      formatC( sum_fit_phy$r.squared, digits = 2)
PGLS_pvalue =  formatC(sum_fit_phy$coefficients[8], digits = 3)
if (PGLS_pvalue == "   0"){ PGLS_pvalue = "< 2.2e-16"}
PGLS_cc <- coef(fit_phylo_Cerebellum_OR)
PGLS_fit_function <- function(x) PGLS_cc[1] + PGLS_cc[2]*x
PGLS_lambda <- sum_fit_phy$param[2]
PGLS_r2


lm_fit <- lm(data = Cichlids_lamellae_brains_df %>% filter(Sp %in% b1_tree_lamellae$tip.label), 
             formula = Total_OR ~ prop_Cerebellum)

sum_fit <- summary(lm_fit)
GLS_r2 =      formatC( sum_fit$r.squared, digits = 3)
GLS_pente =   formatC(sum_fit$coefficients[2], digits = 3)
GLS_pvalue =  formatC(sum_fit$coefficients[8], digits = 3)
GLS_cc <- coef(lm_fit)
GLS_fit_function <- function(x) GLS_cc[1] + GLS_cc[2]*x

pdf(file = "Raw_R_plots/OR_vs_Cerebellumprop.pdf",width = 8.34,  height = 4.61)

Cichlids_lamellae_brains_df %>% 
  ggplot(aes( x = prop_Cerebellum, y= Total_OR)) +
  geom_point(aes(color=Tribe)) +
  stat_function(fun = GLS_fit_function, color="black")+
  labs(subtitle = bquote(R^2 ~ "=" ~ .(PGLS_r2) ~ ";" ~ P ~ "=" ~.(PGLS_pvalue) ~ ";" ~ lambda ~ "=" ~ .(PGLS_lambda))) +
  theme_classic() +
  xlab("Mean relative cerebellum volume") +
  ylab("Mean number of OR") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none") + 
  scale_color_manual(values=tribes.colors) 

dev.off()


#Lets plot V2R vs Cerebellum
sum_fit_phy <- summary(fit_phylo_Cerebellum_V2R)
PGLS_pente =   formatC(sum_fit_phy$coefficients[2], digits = 3)
PGLS_r2 =      formatC( sum_fit_phy$r.squared, digits = 2)
PGLS_pvalue =  formatC(sum_fit_phy$coefficients[8], digits = 3)
if (PGLS_pvalue == "   0"){ PGLS_pvalue = "< 2.2e-16"}
PGLS_cc <- coef(fit_phylo_Cerebellum_V2R)
PGLS_fit_function <- function(x) PGLS_cc[1] + PGLS_cc[2]*x
PGLS_lambda <- sum_fit_phy$param[2]
PGLS_r2


lm_fit <- lm(data = Cichlids_lamellae_brains_df %>% filter(Sp %in% b1_tree_lamellae$tip.label), 
             formula = Total_V2R ~ prop_Cerebellum)

sum_fit <- summary(lm_fit)
GLS_r2 =      formatC( sum_fit$r.squared, digits = 3)
GLS_pente =   formatC(sum_fit$coefficients[2], digits = 3)
GLS_pvalue =  formatC(sum_fit$coefficients[8], digits = 3)
GLS_cc <- coef(lm_fit)
GLS_fit_function <- function(x) GLS_cc[1] + GLS_cc[2]*x

pdf(file = "Raw_R_plots/V2R_vs_Cerebellumprop.pdf",width = 8.34,  height = 4.61)

Cichlids_lamellae_brains_df %>% 
  ggplot(aes( x = prop_Cerebellum, y= Total_V2R)) +
  geom_point(aes(color=Tribe)) +
  stat_function(fun = GLS_fit_function, color="black")+
  labs(subtitle = bquote(R^2 ~ "=" ~ .(PGLS_r2) ~ ";" ~ P ~ "=" ~.(PGLS_pvalue) ~ ";" ~ lambda ~ "=" ~ .(PGLS_lambda))) +
  theme_classic() +
  xlab("Mean relative cerebellum volume") +
  ylab("Mean number of V2R") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none") + 
  scale_color_manual(values=tribes.colors) 

dev.off()


#Lets plot TAAR vs Cerebellum
sum_fit_phy <- summary(fit_phylo_Cerebellum_TAAR)
PGLS_pente =   formatC(sum_fit_phy$coefficients[2], digits = 3)
PGLS_r2 =      formatC( sum_fit_phy$r.squared, digits = 2)
PGLS_pvalue =  formatC(sum_fit_phy$coefficients[8], digits = 3)
if (PGLS_pvalue == "   0"){ PGLS_pvalue = "< 2.2e-16"}
PGLS_cc <- coef(fit_phylo_Cerebellum_TAAR)
PGLS_fit_function <- function(x) PGLS_cc[1] + PGLS_cc[2]*x
PGLS_lambda <- sum_fit_phy$param[2]
PGLS_r2


lm_fit <- lm(data = Cichlids_lamellae_brains_df %>% filter(Sp %in% b1_tree_lamellae$tip.label), 
             formula = Total_TAAR ~ prop_Cerebellum)

sum_fit <- summary(lm_fit)
GLS_r2 =      formatC( sum_fit$r.squared, digits = 3)
GLS_pente =   formatC(sum_fit$coefficients[2], digits = 3)
GLS_pvalue =  formatC(sum_fit$coefficients[8], digits = 3)
GLS_cc <- coef(lm_fit)
GLS_fit_function <- function(x) GLS_cc[1] + GLS_cc[2]*x

pdf(file = "Raw_R_plots/TAAR_vs_Cerebellumprop.pdf",width = 8.34,  height = 4.61)

Cichlids_lamellae_brains_df %>% 
  ggplot(aes( x = prop_Cerebellum, y= Total_TAAR)) +
  geom_point(aes(color=Tribe)) +
  stat_function(fun = GLS_fit_function, color="black")+
  labs(subtitle = bquote(R^2 ~ "=" ~ .(PGLS_r2) ~ ";" ~ P ~ "=" ~.(PGLS_pvalue) ~ ";" ~ lambda ~ "=" ~ .(PGLS_lambda))) +
  theme_classic() +
  xlab("Mean relative cerebellum volume") +
  ylab("Mean number of TAAR") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none") + 
  scale_color_manual(values=tribes.colors) 

dev.off()




