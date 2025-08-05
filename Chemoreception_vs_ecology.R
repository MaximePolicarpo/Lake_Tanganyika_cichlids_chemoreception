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

#### Data load - Genomics , Lamellae and Brain  ---------------------------------

Cichlids_lamellae_brains_df <- 
  read.table("Cichlids_lamellae_brains_df.csv",sep=",", header=TRUE)

#Remove Neospl
Cichlids_lamellae_brains_df <- 
  Cichlids_lamellae_brains_df %>%
  filter(Sp != "Neospl")

Cichlids_lamellae_brains_df <- 
  Cichlids_lamellae_brains_df %>%
  dplyr::select(- c(mean_SL, mean_TL, mean_weight))


#### Data load -- Subfamilies number of genes ----------

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



All_Doc_mean_sp <- 
  rbind(OR_DoC_df_info_mean_total, TAAR_DoC_df_info_mean_total, 
        V1R_DoC_df_info_mean_total, V2R_DoC_df_info_mean_total,
        T1R_DoC_df_info_mean_total, T2R_DoC_df_info_mean_total)

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



OR_DoC_df_info_mean_wo_total <- 
  OR_DoC_df_info_mean %>%
  filter(Subfamily != "Total_OR")

TAAR_DoC_df_info_mean_wo_total <- 
  TAAR_DoC_df_info_mean %>%
  filter(Subfamily != "Total_TAAR")

V1R_DoC_df_info_mean_wo_total <- 
  V1R_DoC_df_info_mean %>%
  filter(Subfamily != "Total_V1R")

V2R_DoC_df_info_mean_wo_total <- 
  V2R_DoC_df_info_mean %>%
  filter(Subfamily != "Total_V2R")


T1R_DoC_df_info_mean_wo_total <- 
  T1R_DoC_df_info_mean %>%
  filter(Subfamily != "Total_T1R")

T2R_DoC_df_info_mean_wo_total <- 
  T2R_DoC_df_info_mean %>%
  filter(Subfamily != "Total_T2R")


CR_DoC_df_info_mean_wo_total <- 
  rbind(OR_DoC_df_info_mean_wo_total, TAAR_DoC_df_info_mean_wo_total,
        V1R_DoC_df_info_mean_wo_total, V2R_DoC_df_info_mean_wo_total,
        T1R_DoC_df_info_mean_wo_total, T2R_DoC_df_info_mean_wo_total)

CR_DoC_df_info_mean_wo_total[(CR_DoC_df_info_mean_wo_total$Tribe == "Serranochromini"),"Tribe"] <- "Haplochromini"

CR_Doc_sp_wide <- 
  as.data.frame(
    pivot_wider(CR_DoC_df_info_mean_wo_total, 
                names_from = c(Subfamily), values_from = mean_normalized_nb)
  )


#### Data load - Genomics PCA  ---------------------------------

gPCA_df <- 
  read.table("gPC1_2_3.csv",sep=",", header=TRUE)



#### Data load - Species tree  ---------------------------------


#load phylo tree
b1_tree <- read.tree("b1_tree_filt.nwk")
b1_tree_wo_Neospl <- 
  drop.tip(b1_tree, "Neospl")
b1_list_sp <- b1_tree_wo_Neospl$tip.label
b1_tree_wo_Neospl_NodeLabel <- makeNodeLabel(b1_tree_wo_Neospl, method="number", prefix="Node")




#### Data load - Ecology  ---------------------------------


#Import table with ecological variables per species 
ecoCat <- read.table("Species_ecoCat.csv", sep=",")
colnames(ecoCat) <- c("Sp", "Tribe", "Food", "habitat")
ecoCat <- ecoCat %>% filter(! Sp == "Neospl")



#SI table
SI_infos <- read.table("Stable_isotope_table_simplified.tsv",
                       sep="\t", header=TRUE)
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
breed_infos <- read.table("Species_tribe_breeding.tsv",
                          sep="\t", header=TRUE)
colnames(breed_infos) <- c("Sp", "Species", "Tribe", "breeding_type", "breeding_mode")
breed_infos <- breed_infos %>% filter(! Sp == "Neospl")

#### Data load - Behavior  ---------------------------------

Explo_behavior <- read.table("Median_exploratory_behavior.tsv", 
                             sep="\t", header=TRUE)


#### Data load - Standard length  ---------------------------------

Cichlids_length_weight <- 
  read.table("SL_TL_Weight_df.tsv", sep="\t", header=TRUE)
colnames(Cichlids_length_weight) <- c("Sp", "SL", "TL", "Weight")


Cichlids_length_weight_sp <- 
  as.data.frame(
    Cichlids_length_weight %>%
  group_by(Sp) %>%
  summarise(
    mean_SL = mean(SL, na.rm=TRUE),
    mean_TL = mean(TL, na.rm=TRUE),
    mean_weight = mean(Weight, na.rm=TRUE))
  )


#### Combine dataframes and prepare pGLS caper data  ---------------------------------

#Define the tribe colors (same as Ronco et al. 2021)
tribes.colors <- c(Bathybatini = "#242626" , Benthochromini = "#AE262A" , Boulengerochromini = "#59595C" , Cyphotilapiini = "#FDDF13" , Cyprichromini = "#F04D29" , Ectodini = "#9AB9D9" , Eretmodini = "#682E7A" , Lamprologini = "#C588BB" , Limnochromini = "#535CA9" , Perissodini = "#FBAC43" , Trematocarini = "#959170" , Tropheini = "#86C773" , Haplochromini = "#274e13" , outgroup = "gray", Other="gray")
sp.colors <- c(Batmin = "#242626" , Cypfro = "#FDDF13" , Cyplep = "#F04D29" , Cunlon = "#9AB9D9" , Neomul = "#C588BB" , Simdia = "#86C773")


Cichlids_chemoreception_df <- 
  left_join(Cichlids_lamellae_brains_df, ecoCat %>% dplyr::select(-Tribe), by="Sp")
Cichlids_chemoreception_df <- 
  left_join(Cichlids_chemoreception_df, SI_infos_mean, by="Sp")
Cichlids_chemoreception_df <- 
  left_join(Cichlids_chemoreception_df, breed_infos %>% dplyr::select(-c(Tribe, Species)), by="Sp")
Cichlids_chemoreception_df <- 
  left_join(Cichlids_chemoreception_df, Explo_behavior, by="Sp")
Cichlids_chemoreception_df <- 
  left_join(Cichlids_chemoreception_df, gPCA_df, by="Sp")
Cichlids_chemoreception_df <- 
  left_join(Cichlids_chemoreception_df, CR_Doc_sp_wide, by=c("Sp", "Species", "Tribe"))

Cichlids_chemoreception_df <- 
  left_join(Cichlids_chemoreception_df, Cichlids_length_weight_sp, by=c("Sp"))



write.table(Cichlids_chemoreception_df,
            "Cichlids_chemoreception_df.csv",
            quote=FALSE,
            sep=",",
            col.names=TRUE,
            row.names=FALSE)


caper_data_chemoreception <- 
  comparative.data(phy = b1_tree_wo_Neospl_NodeLabel, 
                   data = Cichlids_chemoreception_df,
                   names.col = Sp, vcv = TRUE,
                   na.omit = FALSE, warn.dropped = TRUE)



chemo_responses <- 
  c("Total_OLR","Total_OR", "Total_TAAR", "Total_V2R", "Total_T1R", "prop_Olfactory_bulb", "mean_lamellae",
    "gPC1", "gPC2", "gPC3", 
    colnames(CR_Doc_sp_wide %>% dplyr::select(-c("Sp", "Species", "Tribe", "T2RA", "ORA1", "ORA2",
                                                 "ORA3", "ORA4", "ORA5", "ORA6", "Kappa-1", "Lambda-1",
                                                 "Theta-1", "TAARL-1", "V2RD2-1")))
    )
 

#### Plot a tree with the ecological factors  -------------------


Tribe_table <- Cichlids_chemoreception_df %>% dplyr::select(Sp, Tribe)
SL_table <- Cichlids_chemoreception_df %>% dplyr::select(Sp, mean_SL)
d13C_table <- Cichlids_chemoreception_df %>% dplyr::select(Sp, mean_d13C)
d15N_table <- Cichlids_chemoreception_df %>% dplyr::select(Sp, mean_d15N)
Explo_table <- Cichlids_chemoreception_df %>% dplyr::select(Sp, median_exploration)
Food_table <- Cichlids_chemoreception_df %>% dplyr::select(Sp, Food)
Habitat_table <- Cichlids_chemoreception_df %>% dplyr::select(Sp, habitat)
BM_table <- Cichlids_chemoreception_df %>% dplyr::select(Sp, breeding_mode)
BT_table <- Cichlids_chemoreception_df %>% dplyr::select(Sp, breeding_type)


Food_colors <- 
  c("invert" = "#E69F00",
    "Auf_invert" = "#56B4E9",
    "fish" = "#000000",
    "omni" = "#009E73",
    "Auf_herb" = "#F0E442",
    "herb" = "#0072B2",
    "fish_invert" = "#D55E00",
    "plan" = "#CC79A7",
    "plan_invert" = "#757879",
    "fry_plan" = "#81F9E5",
    "scales" = "#CD00FF")
  
Habitat_colors <- 
  c("deep" = "#E69F00",
    "pelagic" = "#56B4E9",
    "litoral" = "#000000",
    "inter" = "#009E73",
    "Auf_herb" = "#F0E442",
    "shallow" = "#0072B2")

breeding_mode_colors <- 
  c("cave spawning" = "#E69F00",
    "biparental brooding" = "#56B4E9",
    "maternal bower brooding" = "#000000",
    "maternal territorial brooding" = "#009E73",
    "open spawning" = "#F0E442",
    "maternal open water brooding" = "#0072B2",
    "shell spawning" = "#D55E00")



breeding_type_colors <- 
  c("Substrate spawning" = "#D55E00",
    "Mouthbrooding" = "#56B4E9")


pdf(file = "Raw_R_plots/Species_tree_Ecology.pdf",width = 6.34,  height = 4.61)

ggtree(b1_tree_wo_Neospl_NodeLabel, layout="circular", size=0.4) %<+% Tribe_table +
  aes(color=Tribe) +
  geom_tiplab(size=0.6, offset=0.07) +
  scale_color_manual(values = tribes.colors) +
  new_scale_colour() + 
  geom_fruit(
    data=SL_table,
    geom=geom_tile,
    mapping=aes(y=Sp, fill=mean_SL),
    size=0.02,
    width=0.3,
    color="white",
    pwidth=0.1,
    position="auto",
    offset=0.15,
    grid.params=list()) +
  scale_fill_viridis(option="magma", direction=-1, na.value="white") + 
  new_scale_fill()  +
  geom_fruit(
    data=d13C_table,
    geom=geom_tile,
    mapping=aes(y=Sp, fill=mean_d13C),
    size=0.02,
    width=0.3,
    color="white",
    pwidth=0.1,
    position="auto",
    offset=0.05,
    grid.params=list()) +
  scale_fill_viridis(option="magma", direction=-1, na.value="white") + 
  new_scale_fill() +
  geom_fruit(
    data=d15N_table,
    geom=geom_tile,
    mapping=aes(y=Sp, fill=mean_d15N),
    size=0.02,
    width=0.3,
    color="white",
    pwidth=0.1,
    position="auto",
    offset=0.05,
    grid.params=list()) +
  scale_fill_viridis(option="magma", direction=-1, na.value="white") + 
  new_scale_fill()  +
  geom_fruit(
    data=Explo_table,
    geom=geom_tile,
    mapping=aes(y=Sp, fill=median_exploration),
    size=0.02,
    width=0.3,
    color="white",
    pwidth=0.1,
    position="auto",
    offset=0.05,
    grid.params=list()) +
  scale_fill_viridis(option="magma", direction=-1, na.value="white") + 
  new_scale_fill() +
  geom_fruit(
    data=Food_table,
    geom=geom_tile,
    mapping=aes(y=Sp, fill=Food),
    size=0.02,
    width=0.3,
    color="white",
    pwidth=0.1,
    position="auto",
    offset=0.1,
    grid.params=list() 
  ) +
  scale_fill_manual(values = Food_colors) + 
  new_scale_fill() + 
  geom_fruit(
    data=Habitat_table,
    geom=geom_tile,
    mapping=aes(y=Sp, fill=habitat),
    size=0.02,
    width=0.3,
    color="white",
    pwidth=0.1,
    position="auto",
    offset=0.05,
    grid.params=list() 
  ) +
  scale_fill_manual(values = Habitat_colors) + 
  new_scale_fill() + 
  geom_fruit(
    data=BM_table,
    geom=geom_tile,
    mapping=aes(y=Sp, fill=breeding_mode),
    size=0.02,
    width=0.3,
    color="white",
    pwidth=0.1,
    position="auto",
    offset=0.05,
    grid.params=list() 
  ) +
  scale_fill_manual(values = breeding_mode_colors) + 
  new_scale_fill() + 
  geom_fruit(
    data=BT_table,
    geom=geom_tile,
    mapping=aes(y=Sp, fill=breeding_type),
    size=0.02,
    width=0.3,
    color="white",
    pwidth=0.1,
    position="auto",
    offset=0.05,
    grid.params=list() 
  ) +
  scale_fill_manual(values = breeding_type_colors) + 
  new_scale_fill() + 
  theme(legend.position = "none") 


dev.off()


ggtree(b1_tree_wo_Neospl_NodeLabel, layout="circular", size=0.4) %<+% Tribe_table +
  aes(color=Tribe) +
  geom_tiplab(size=0.6, offset=0.07) +
  scale_color_manual(values = tribes.colors) +
  new_scale_colour() + 
  geom_fruit(
    data=SL_table,
    geom=geom_tile,
    mapping=aes(y=Sp, fill=mean_SL),
    size=0.02,
    width=0.3,
    color="white",
    pwidth=0.1,
    position="auto",
    offset=0.1,
    grid.params=list()) +
  scale_fill_viridis(option="magma", direction=-1, na.value="white") + 
  new_scale_fill()  +
  geom_fruit(
    data=d13C_table,
    geom=geom_tile,
    mapping=aes(y=Sp, fill=mean_d13C),
    size=0.02,
    width=0.3,
    color="white",
    pwidth=0.1,
    position="auto",
    offset=0.05,
    grid.params=list()) +
  scale_fill_viridis(option="magma", direction=-1, na.value="white") + 
  new_scale_fill() +
  geom_fruit(
    data=d15N_table,
    geom=geom_tile,
    mapping=aes(y=Sp, fill=mean_d15N),
    size=0.02,
    width=0.3,
    color="white",
    pwidth=0.1,
    position="auto",
    offset=0.05,
    grid.params=list()) +
  scale_fill_viridis(option="magma", direction=-1, na.value="white") + 
  new_scale_fill()  +
  geom_fruit(
    data=Explo_table,
    geom=geom_tile,
    mapping=aes(y=Sp, fill=median_exploration),
    size=0.02,
    width=0.3,
    color="white",
    pwidth=0.1,
    position="auto",
    offset=0.05,
    grid.params=list()) +
  scale_fill_viridis(option="magma", direction=-1, na.value="white") + 
  new_scale_fill() +
  geom_fruit(
    data=Food_table,
    geom=geom_tile,
    mapping=aes(y=Sp, fill=Food),
    size=0.02,
    width=0.3,
    color="white",
    pwidth=0.1,
    position="auto",
    offset=0.1,
    grid.params=list() 
  ) +
  scale_fill_manual(values = Food_colors) + 
  new_scale_fill() + 
  geom_fruit(
    data=Habitat_table,
    geom=geom_tile,
    mapping=aes(y=Sp, fill=habitat),
    size=0.02,
    width=0.3,
    color="white",
    pwidth=0.1,
    position="auto",
    offset=0.05,
    grid.params=list() 
  ) +
  scale_fill_manual(values = Habitat_colors) + 
  new_scale_fill() + 
  geom_fruit(
    data=BM_table,
    geom=geom_tile,
    mapping=aes(y=Sp, fill=breeding_mode),
    size=0.02,
    width=0.3,
    color="white",
    pwidth=0.1,
    position="auto",
    offset=0.05,
    grid.params=list() 
  ) +
  scale_fill_manual(values = breeding_mode_colors) + 
  new_scale_fill() + 
  geom_fruit(
    data=BT_table,
    geom=geom_tile,
    mapping=aes(y=Sp, fill=breeding_type),
    size=0.02,
    width=0.3,
    color="white",
    pwidth=0.1,
    position="auto",
    offset=0.05,
    grid.params=list() 
  ) +
  scale_fill_manual(values = breeding_type_colors) + 
  new_scale_fill() 





ggtree(b1_tree_wo_Neospl_NodeLabel, layout="circular", size=0.4) %<+% Tribe_table +
  aes(color=Tribe) +
  geom_tiplab(size=0.6, offset=0.07) +
  scale_color_manual(values = tribes.colors) +
  new_scale_colour() + 
  new_scale_fill() + 
  geom_fruit(
    data=Habitat_table,
    geom=geom_tile,
    mapping=aes(y=Sp, fill=habitat),
    size=0.02,
    width=0.3,
    color="white",
    pwidth=0.1,
    position="auto",
    offset=0.05,
    grid.params=list() 
  ) +
  scale_fill_manual(values = Habitat_colors) + 
  new_scale_fill() 

#### pGLS with numerical predictors  ---------------------------------


numerical_predictor <- 
  c("mean_SL", "mean_weight", "mean_d15N", "mean_d13C", "median_exploration")

abnorm_termination <- 
  c("Delta-16_mean_d13C", "Zeta-4_mean_d15N", "TAARB4-4_mean_d15N",
    "V2RD7-1_mean_d15N", "Eta-10_mean_SL")

Chemoreception_numerical_pGLS_df <- as.data.frame(NULL)
for(curr_response in chemo_responses){
  for(curr_predictor in numerical_predictor){
    
    current_test_name <- paste(curr_response, curr_predictor, sep="_")
    print(current_test_name)
    
    curr_data_chemoreception <- 
      Cichlids_chemoreception_df %>%
      dplyr::select(Sp, curr_response, curr_predictor)
    
    
    colnames(curr_data_chemoreception) <- c("Sp", "reponse", "predictor")
    
    curr_caper_data <- 
      comparative.data(phy = b1_tree_wo_Neospl_NodeLabel, 
                       data = curr_data_chemoreception,
                       names.col = Sp, vcv = TRUE,
                       na.omit = FALSE, warn.dropped = TRUE)
    
    
    if(current_test_name %in% abnorm_termination){ #reduce the lambda search boundary for failed computations
      curr_pgls_result <-
        pgls(reponse ~ predictor, 
             data = curr_caper_data, 
             lambda = "ML",
             bounds=list(lambda=c(0.7,1))
        )
    } else {
      curr_pgls_result <-
        pgls(reponse ~ predictor, 
             data = curr_caper_data, 
             lambda = "ML",
             bounds=list(lambda=c(0,1))
        )
    }
   

           
           
    sum_fit_phy <- summary(curr_pgls_result)
    PGLS_pente =   formatC(sum_fit_phy$coefficients[2], digits = 3)
    PGLS_r2 =      formatC(sum_fit_phy$adj.r.squared, digits = 2)
    PGLS_pvalue =  formatC(sum_fit_phy$coefficients[8], digits = 3)
    if (PGLS_pvalue == "   0"){ PGLS_pvalue = "< 2.2e-16"}
    lambda <- as.numeric(sum_fit_phy$param[2])
    slope <- formatC(sum_fit_phy$coefficients[2], digits = 3)
    
    print(PGLS_pvalue)
    
    curr_df <- 
      as.data.frame(
        cbind(
          curr_response,
          PGLS_r2,
          PGLS_pvalue,
          lambda,
          slope,
          curr_predictor
        )
      )
    
    colnames(curr_df) <- c("Response", "R2", "pvalue", "lambda", "slope","predictor")
    
    Chemoreception_numerical_pGLS_df <- 
      rbind(Chemoreception_numerical_pGLS_df, curr_df)
  }
  
}

Chemoreception_numerical_pGLS_df %>%
  filter(pvalue < 0.05)

#write.table(Chemoreception_numerical_pGLS_df,
#            "Chemoreception_PGLS_numerical_data.tsv",
#            col.names=TRUE,
#            row.names=FALSE,
#            quote=FALSE,
#            sep=",")
#



#### pGLS with categorical predictors  -------------------


categorical_predictors <- 
  c("Food", "habitat", "breeding_type", "breeding_mode")

abnorm_termination <- 
  c("Beta-2_Food", "Delta-11_breeding_mode", "Epsilon-3_Food",
    "Eta-11_habitat", "V2RD7-1_habitat")


Chemoreception_categorical_pGLS_df <- as.data.frame(NULL)
for(curr_response in chemo_responses){
  for(curr_predictor in categorical_predictors){
    
    
    current_test_name <- paste(curr_response, curr_predictor, sep="_")
    print(current_test_name)
    
    curr_data_chemoreception <- 
      Cichlids_chemoreception_df %>%
      dplyr::select(Sp, curr_response, curr_predictor) %>%
      filter(! is.na(curr_predictor))
    colnames(curr_data_chemoreception) <- c("Sp", "reponse", "predictor")
    
    

    curr_data_caper <- 
      comparative.data(phy = b1_tree_wo_Neospl_NodeLabel, 
                       data = curr_data_chemoreception,
                       names.col = Sp, vcv = TRUE,
                       na.omit = FALSE, warn.dropped = TRUE)
    
    
    
    if(current_test_name %in% abnorm_termination){ #reduce the lambda search boundary for failed computations
      curr_pgls_result <-
        pgls(reponse ~ predictor, 
             data = curr_data_caper, 
             lambda = "ML",
             bounds=list(lambda=c(0.5,1)))
    } else {
      curr_pgls_result <-
        pgls(reponse ~ predictor, 
             data = curr_data_caper, 
             lambda = "ML",
             bounds=list(lambda=c(0,1)))
    }
    
    
    
    
    curr_pgls_summary <- summary(curr_pgls_result)
    
    pGLS_R2 <- as.numeric(curr_pgls_summary$adj.r.squared)
    pGLS_lambda <- as.numeric(curr_pgls_summary$param[2])
    
    pGLS_pvalue <- pf(curr_pgls_summary$fstatistic[1], 
                      curr_pgls_summary$fstatistic[2], 
                      curr_pgls_summary$fstatistic[3], 
                      lower.tail = FALSE)
    
    
    curr_df <- 
      as.data.frame(
        cbind(
          curr_response,
          pGLS_R2,
          pGLS_pvalue,
          pGLS_lambda,
          curr_predictor
        )
      )
    
    colnames(curr_df) <- c("Response", "R2", "pvalue", "lambda","predictor")
    
    
    Chemoreception_categorical_pGLS_df <- 
      rbind(Chemoreception_categorical_pGLS_df,
            curr_df)
    
    
  }
  
}

#write.table(Chemoreception_categorical_pGLS_df,
#            "Chemoreception_PGLS_categorical_data.tsv",
#            col.names=TRUE,
#            row.names=FALSE,
#            quote=FALSE,
#            sep=",")

Chemoreception_categorical_pGLS_df <- 
  read.table("Chemoreception_PGLS_categorical_data.tsv",
             sep=",",
             header=TRUE)



#### pGLS Summary  -------------------

#Import tables

Chemoreception_numerical_pGLS_df <-
  read.table("Chemoreception_PGLS_numerical_data.tsv",
             sep=",",
             header=TRUE)


Chemoreception_categorical_pGLS_df <-
  read.table("Chemoreception_PGLS_categorical_data.tsv",
             sep=",",
             header=TRUE)
Chemoreception_categorical_pGLS_df <- Chemoreception_categorical_pGLS_df %>%  mutate(slope = NA)

Chemoreception_categorical_pGLS_df <- 
  Chemoreception_categorical_pGLS_df %>%
  dplyr::select(Response,R2, pvalue, lambda, slope, predictor)

colnames(Chemoreception_categorical_pGLS_df) <- colnames(Chemoreception_numerical_pGLS_df)


#Merge tables

Chemoreception_pGLS_df <- 
  rbind(Chemoreception_numerical_pGLS_df, Chemoreception_categorical_pGLS_df)


#Summarise the data
as.data.frame(
  Chemoreception_pGLS_df %>%
    group_by(Response) %>%
    summarise(n()))

as.data.frame(
  Chemoreception_pGLS_df %>%
    group_by(predictor) %>%
    summarise(n()))


Chemoreception_pGLS_df <-
  Chemoreception_pGLS_df %>%
  filter(predictor != "mean_weight")
#Add the FDR and Bonferroni p-values 

unique_response <- Chemoreception_pGLS_df$Response %>% unique()
Chemoreception_pGLS_corr_df <- as.data.frame(NULL)
for(curr_response in unique_response){
  
  curr_df <- 
    Chemoreception_pGLS_df %>%
    filter(Response == curr_response)
  
  
  initial_pvalues <- 
    as.numeric(
      curr_df %>% pull(pvalue)
    )
  
  corrected_pvalues_FDR <- p.adjust(initial_pvalues, method = "fdr")
  corrected_pvalues_bonferroni <- p.adjust(initial_pvalues, method = "bonferroni")
  
  curr_df <- 
    curr_df %>%
    mutate(FDR_pvalue = corrected_pvalues_FDR,
           Bonferroni_pvalue = corrected_pvalues_bonferroni)
  
  Chemoreception_pGLS_corr_df <- 
    rbind(Chemoreception_pGLS_corr_df, curr_df)
  
  
}


## Define significance levels with the difference corrections

Chemoreception_pGLS_corr_df <- 
  Chemoreception_pGLS_corr_df %>%
  mutate(FDR_significant = if_else(
    as.numeric(FDR_pvalue) < 0.05,
    "Significant", 
    "N.S"
  ))

Chemoreception_pGLS_corr_df <- 
  Chemoreception_pGLS_corr_df %>%
  mutate(Bonferroni_significant = if_else(
    as.numeric(Bonferroni_pvalue) < 0.05,
    "Significant", 
    "N.S"
  ))

Chemoreception_pGLS_corr_df <- 
  Chemoreception_pGLS_corr_df %>%
  mutate(significant = if_else(
    as.numeric(pvalue) < 0.05,
    "Significant", 
    "N.S"
  ))


## Define the slope ...
Chemoreception_pGLS_corr_df <- 
  Chemoreception_pGLS_corr_df %>%
  mutate(slope_sign = case_when(
    slope <= 0 ~ "minus",
    slope > 0 ~ "plus",
    is.na(slope) ~ "cat",
  ))


## Define a shape for slopes for future graphics
slope_sign_shapes <- 
  c("minus"="\u25BC",
    "plus"="\u25B2",
    "cat"="\u25FE")

#Remove chemoreceptors subfamilies with only one gene per species


Chemoreception_pGLS_corr_df <- 
  Chemoreception_pGLS_corr_df %>%
  filter(! Response %in% 
           c("T1R1A", "ORA1", "ORA2", "ORA3", "ORA4", "ORA5", "ORA6", "Beta-1","Delta-13", 
             "Eta-4", "TAARL-1","TAARB4-5", "V2RD2-1","V2RD11-3", "V2RD11-4", "V2RD11-6", "V2RD5-1"))


#Re-order response for the graphic 

Chemoreception_pGLS_corr_df$Response <-
  factor(Chemoreception_pGLS_corr_df$Response ,
         levels=c("prop_Olfactory_bulb","mean_lamellae",
                  "gPC1","gPC2","gPC3","Total_OLR",
                  
                  "Total_OR",
                  #"Beta-1",
                  "Beta-2","Delta-1","Delta-2",
                  "Delta-3","Delta-4","Delta-5","Delta-6","Delta-7","Delta-8","Delta-9",
                  "Delta-10","Delta-11","Delta-12",
                  #"Delta-13",
                  "Delta-14","Delta-15","Delta-16","Delta-17",
                  "Delta-18","Delta-19","Delta-20",
                  "Epsilon-1","Epsilon-2","Epsilon-3","Epsilon-4",
                  "Eta-1","Eta-2","Eta-3",
                  #"Eta-4",
                  "Eta-6","Eta-7","Eta-8","Eta-9",
                  "Eta-10","Eta-11","Eta-12","Eta-13","Eta-14","Eta-15",
                  "Zeta-1","Zeta-2","Zeta-3","Zeta-4",
                  
                  "Total_TAAR","TAARA2-1","TAARA2-2","TAARB4-1",
                  "TAARB4-2","TAARB4-3","TAARB4-4",
                  #"TAARB4-5",
                  
                  "Total_V2R","V2RD3-1","V2RD4-1",
                  #"V2RD5-1",
                  "V2RD7-1","V2RD8-1","V2RD9-1",
                  "V2RD10-1","V2RD11-1","V2RD11-2",
                  #"V2RD11-3","V2RD11-4",
                  "V2RD11-5",
                  #"V2RD11-6",
                  "V2RD11-7",
                  
                  
                  #"Total_T1R","T1R1A","T1R1B","T1R3"
                  "Total_T1R","T1R1B","T1R3"
         ))


Chemoreception_pGLS_corr_df$predictor <-
  factor(Chemoreception_pGLS_corr_df$predictor ,
         levels=c("mean_SL", 
                  "median_exploration", 
                  "mean_d15N", "Food",
                  "mean_d13C", "habitat",
                  "breeding_type", "breeding_mode"))



Signif_colors <-
  c("N.S" = 0,
    "Significant" = 1)



#Replace negative R2 by R2 of 0

Chemoreception_pGLS_corr_df$R2[Chemoreception_pGLS_corr_df$R2 < 0] <- 0

#Lets plot
pdf(file = "Raw_R_plots/Ecological_pGLS_nocorr.pdf", width = 10.34,  height = 4.61)

ggplot(Chemoreception_pGLS_corr_df, 
       aes(y=predictor, x=Response, fill= R2, alpha=significant, color="black")) + 
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

ggplot(Chemoreception_pGLS_corr_df, 
       aes(y=predictor, x=Response, fill= R2, alpha=significant, color="black")) + 
  geom_tile(color = "black") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme_minimal() +
  scale_fill_viridis_c(option = "plasma") +
  scale_alpha_manual(values = Signif_colors) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("") +
  ylab("") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


ggplot(Chemoreception_pGLS_corr_df, 
       aes(y=predictor, x=Response, fill= R2, alpha=significant, 
           linewidth=FDR_significant, color=FDR_significant)) + 
  geom_tile(color = "black") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_viridis_c(option = "plasma") +
  scale_alpha_manual(values = Signif_colors) +
  scale_discrete_manual("linewidth", values = c("N.S" = 0.05, "Significant" = 0.8)) +
  scale_color_manual(values = c("N.S" = "gray", "Significant" = "black")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("") +
  ylab("") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

Chemoreception_pGLS_corr_df_g <- 
  Chemoreception_pGLS_corr_df %>%
  mutate(FDR_significant_TF = if_else(
    FDR_significant == "N.S",
    "FALSE",
    "TRUE"
  ))

ggplot(Chemoreception_pGLS_corr_df_g, 
       aes(y=predictor, x=Response, fill= R2, alpha=significant)) + 
  geom_tile(color = "black") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme_minimal() +
  geom_text(aes(label = ifelse(FDR_significant_TF, "*", "")), color = "white", size=6) +
  theme(legend.position = "none") +
  scale_fill_viridis_c(option = "plasma") +
  scale_alpha_manual(values = Signif_colors) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("") +
  ylab("") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

dev.off()


#Lets plot
pdf(file = "Raw_R_plots/Ecological_pGLS_FDR.pdf", width = 10.34,  height = 4.61)

ggplot(Chemoreception_pGLS_corr_df, 
       aes(y=predictor, x=Response, fill= R2, alpha=FDR_significant, color="black")) + 
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

ggplot(Chemoreception_pGLS_corr_df, 
       aes(y=predictor, x=Response, fill= R2, alpha=FDR_significant, color="black")) + 
  geom_tile(color = "black") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme_minimal() +
  scale_fill_viridis_c(option = "plasma") +
  scale_alpha_manual(values = Signif_colors) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("") +
  ylab("") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


Chemoreception_pGLS_corr_df_g <- 
  Chemoreception_pGLS_corr_df %>%
  mutate(sign_TF = case_when(
    slope_sign == "plus" ~ "TRUE",
    slope_sign == "minus" ~ "FALSE"
  ))



ggplot(Chemoreception_pGLS_corr_df_g, 
       aes(y=predictor, x=Response, fill= R2, alpha=FDR_significant, color="black")) + 
  geom_tile(color = "black") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_text(aes(label = ifelse(sign_TF, "+", "-")), color = "white") + 
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_viridis_c(option = "plasma") +
  scale_alpha_manual(values = Signif_colors) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("") +
  ylab("") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


dev.off()


Chemoreception_pGLS_corr_df %>%
  filter(Response == "mean_lamellae") %>%
  filter(FDR_significant == "Significant")



Chemoreception_pGLS_corr_df %>%
  filter(Response == "prop_Olfactory_bulb") %>%
  filter(FDR_significant == "Significant")


Chemoreception_pGLS_corr_df %>%
  filter(predictor == "mean_d15N") %>%
  filter(FDR_significant == "Significant")


Chemoreception_pGLS_corr_df %>%
  filter(predictor == "Food") %>%
  filter(FDR_significant == "Significant") %>%
  arrange(R2)

Chemoreception_pGLS_corr_df %>%
  filter(predictor == "habitat") %>%
  filter(FDR_significant == "Significant") %>%
  arrange(R2)


write.table(Chemoreception_pGLS_corr_df %>% filter(! predictor %in% c("breeding_type", "breeding_mode")),
            "Raw_R_plots/Chemoreception_pGLS_corr_df_publi.tsv",
            sep="\t",
            col.names=TRUE,
            row.names=FALSE,
            quote=FALSE)


###For significant categorical variables, lets perform post-hoc test to see
###which categories are differents.

library(nlme)
library(emmeans)

#Habitat x T1R

Habitat_signif_resp <- 
  Chemoreception_pGLS_corr_df %>%
  filter(FDR_significant == "Significant") %>%
  filter(predictor == "habitat") %>%
  pull(Response)

Cichlids_chemoreception_df_ordered <- 
  Cichlids_chemoreception_df %>% filter(! is.na(habitat)) %>%
  filter(Sp %in% b1_tree_wo_Neospl_NodeLabel$tip.label)
b1_tree_wo_Neospl_NodeLabel_habitat <- 
  keep.tip(b1_tree_wo_Neospl_NodeLabel, intersect(b1_tree_wo_Neospl_NodeLabel$tip.label, Cichlids_chemoreception_df_ordered$Sp))
spp <- b1_tree_wo_Neospl_NodeLabel_habitat$tip.label
corBM<-corBrownian(phy=b1_tree_wo_Neospl_NodeLabel_habitat,form=~Sp)

Cichlids_chemoreception_df_ordered <- 
  Cichlids_chemoreception_df_ordered[match(spp, Cichlids_chemoreception_df_ordered$Sp), ]

ancova_habitat.Total_T1R <-
  gls(Total_T1R ~ habitat,
      data=Cichlids_chemoreception_df_ordered,
      correlation=corBM)
emmeans_result <- emmeans(ancova_habitat.Total_T1R, ~ habitat)

pairwise.habitat.Total_T1R <- 
  as.data.frame(contrast(emmeans_result, method = "pairwise", adjust = "tukey"))

pairwise.habitat.Total_T1R <- 
  pairwise.habitat.Total_T1R %>%
  mutate(group1 = gsub(" -.*", "", contrast)) %>%
  mutate(group2 = gsub(".*- ", "", contrast)) %>%
  mutate(predictor = "habitat") %>%
  mutate(response = "Total_T1R")


mean.habitat.Total_T1R <- 
  Cichlids_chemoreception_df_ordered %>%
  group_by(habitat) %>%
  summarise(group_mean = mean(Total_T1R))
colnames(mean.habitat.Total_T1R) <- c("group1", "mean_group1")
pairwise.habitat.Total_T1R <- left_join(pairwise.habitat.Total_T1R, mean.habitat.Total_T1R, by="group1")
colnames(mean.habitat.Total_T1R) <- c("group2", "mean_group2")
pairwise.habitat.Total_T1R <- left_join(pairwise.habitat.Total_T1R, mean.habitat.Total_T1R, by="group2")



#Habitat x Delta-10

Habitat_signif_resp <- 
  Chemoreception_pGLS_corr_df %>%
  filter(FDR_significant == "Significant") %>%
  filter(predictor == "habitat") %>%
  pull(Response)

Cichlids_chemoreception_df_ordered <- 
  Cichlids_chemoreception_df %>% filter(! is.na(habitat)) %>%
  filter(Sp %in% b1_tree_wo_Neospl_NodeLabel$tip.label)
b1_tree_wo_Neospl_NodeLabel_habitat <- 
  keep.tip(b1_tree_wo_Neospl_NodeLabel, intersect(b1_tree_wo_Neospl_NodeLabel$tip.label, Cichlids_chemoreception_df_ordered$Sp))
spp <- b1_tree_wo_Neospl_NodeLabel_habitat$tip.label
corBM<-corBrownian(phy=b1_tree_wo_Neospl_NodeLabel_habitat,form=~Sp)

Cichlids_chemoreception_df_ordered <- 
  Cichlids_chemoreception_df_ordered[match(spp, Cichlids_chemoreception_df_ordered$Sp), ]

colnames(Cichlids_chemoreception_df_ordered)[colnames(Cichlids_chemoreception_df_ordered) == "Delta-10"] <- "Delta_10"

ancova_habitat.Delta_10 <-
  gls(Delta_10 ~ habitat,
      data=Cichlids_chemoreception_df_ordered,
      correlation=corBM)
emmeans_result <- emmeans(ancova_habitat.Delta_10, ~ habitat)

pairwise.habitat.Delta_10 <- 
  as.data.frame(contrast(emmeans_result, method = "pairwise", adjust = "tukey"))

pairwise.habitat.Delta_10 <- 
  pairwise.habitat.Delta_10 %>%
  mutate(group1 = gsub(" -.*", "", contrast)) %>%
  mutate(group2 = gsub(".*- ", "", contrast)) %>%
  mutate(predictor = "habitat") %>%
  mutate(response = "Delta_10")


mean.habitat.Delta_10 <- 
  Cichlids_chemoreception_df_ordered %>%
  group_by(habitat) %>%
  summarise(group_mean = mean(Delta_10))
colnames(mean.habitat.Delta_10) <- c("group1", "mean_group1")
pairwise.habitat.Delta_10 <- left_join(pairwise.habitat.Delta_10, mean.habitat.Delta_10, by="group1")
colnames(mean.habitat.Delta_10) <- c("group2", "mean_group2")
pairwise.habitat.Delta_10 <- left_join(pairwise.habitat.Delta_10, mean.habitat.Delta_10, by="group2")



#Habitat x Eta-10

Habitat_signif_resp <- 
  Chemoreception_pGLS_corr_df %>%
  filter(FDR_significant == "Significant") %>%
  filter(predictor == "habitat") %>%
  pull(Response)

Cichlids_chemoreception_df_ordered <- 
  Cichlids_chemoreception_df %>% filter(! is.na(habitat)) %>%
  filter(Sp %in% b1_tree_wo_Neospl_NodeLabel$tip.label)
b1_tree_wo_Neospl_NodeLabel_habitat <- 
  keep.tip(b1_tree_wo_Neospl_NodeLabel, intersect(b1_tree_wo_Neospl_NodeLabel$tip.label, Cichlids_chemoreception_df_ordered$Sp))
spp <- b1_tree_wo_Neospl_NodeLabel_habitat$tip.label
corBM<-corBrownian(phy=b1_tree_wo_Neospl_NodeLabel_habitat,form=~Sp)

Cichlids_chemoreception_df_ordered <- 
  Cichlids_chemoreception_df_ordered[match(spp, Cichlids_chemoreception_df_ordered$Sp), ]

colnames(Cichlids_chemoreception_df_ordered)[colnames(Cichlids_chemoreception_df_ordered) == "Eta-10"] <- "Eta_10"

ancova_habitat.Eta_10 <-
  gls(Eta_10 ~ habitat,
      data=Cichlids_chemoreception_df_ordered,
      correlation=corBM)
emmeans_result <- emmeans(ancova_habitat.Eta_10, ~ habitat)

pairwise.habitat.Eta_10 <- 
  as.data.frame(contrast(emmeans_result, method = "pairwise", adjust = "tukey"))

pairwise.habitat.Eta_10 <- 
  pairwise.habitat.Eta_10 %>%
  mutate(group1 = gsub(" -.*", "", contrast)) %>%
  mutate(group2 = gsub(".*- ", "", contrast)) %>%
  mutate(predictor = "habitat") %>%
  mutate(response = "Eta_10")


mean.habitat.Eta_10 <- 
  Cichlids_chemoreception_df_ordered %>%
  group_by(habitat) %>%
  summarise(group_mean = mean(Eta_10))
colnames(mean.habitat.Eta_10) <- c("group1", "mean_group1")
pairwise.habitat.Eta_10 <- left_join(pairwise.habitat.Eta_10, mean.habitat.Eta_10, by="group1")
colnames(mean.habitat.Eta_10) <- c("group2", "mean_group2")
pairwise.habitat.Eta_10 <- left_join(pairwise.habitat.Eta_10, mean.habitat.Eta_10, by="group2")


#Habitat x T1R3

Habitat_signif_resp <- 
  Chemoreception_pGLS_corr_df %>%
  filter(FDR_significant == "Significant") %>%
  filter(predictor == "habitat") %>%
  pull(Response)

Cichlids_chemoreception_df_ordered <- 
  Cichlids_chemoreception_df %>% filter(! is.na(habitat)) %>%
  filter(Sp %in% b1_tree_wo_Neospl_NodeLabel$tip.label)
b1_tree_wo_Neospl_NodeLabel_habitat <- 
  keep.tip(b1_tree_wo_Neospl_NodeLabel, intersect(b1_tree_wo_Neospl_NodeLabel$tip.label, Cichlids_chemoreception_df_ordered$Sp))
spp <- b1_tree_wo_Neospl_NodeLabel_habitat$tip.label
corBM<-corBrownian(phy=b1_tree_wo_Neospl_NodeLabel_habitat,form=~Sp)

Cichlids_chemoreception_df_ordered <- 
  Cichlids_chemoreception_df_ordered[match(spp, Cichlids_chemoreception_df_ordered$Sp), ]

ancova_habitat.T1R3 <-
  gls(T1R3 ~ habitat,
      data=Cichlids_chemoreception_df_ordered,
      correlation=corBM)
emmeans_result <- emmeans(ancova_habitat.T1R3, ~ habitat)

pairwise.habitat.T1R3 <- 
  as.data.frame(contrast(emmeans_result, method = "pairwise", adjust = "tukey"))

pairwise.habitat.T1R3 <- 
  pairwise.habitat.T1R3 %>%
  mutate(group1 = gsub(" -.*", "", contrast)) %>%
  mutate(group2 = gsub(".*- ", "", contrast)) %>%
  mutate(predictor = "habitat") %>%
  mutate(response = "T1R3")


mean.habitat.T1R3 <- 
  Cichlids_chemoreception_df_ordered %>%
  group_by(habitat) %>%
  summarise(group_mean = mean(T1R3))
colnames(mean.habitat.T1R3) <- c("group1", "mean_group1")
pairwise.habitat.T1R3 <- left_join(pairwise.habitat.T1R3, mean.habitat.T1R3, by="group1")
colnames(mean.habitat.T1R3) <- c("group2", "mean_group2")
pairwise.habitat.T1R3 <- left_join(pairwise.habitat.T1R3, mean.habitat.T1R3, by="group2")



#Food x mean_lamellae

Chemoreception_pGLS_corr_df %>%
  filter(FDR_significant == "Significant") %>%
  filter(predictor == "Food") %>%
  pull(Response)

Cichlids_chemoreception_df_ordered <- 
  Cichlids_chemoreception_df %>% filter(! is.na(Food)) %>%
  filter(! is.na(mean_lamellae)) %>%
  filter(Sp %in% b1_tree_wo_Neospl_NodeLabel$tip.label)
b1_tree_wo_Neospl_NodeLabel_Food <- 
  keep.tip(b1_tree_wo_Neospl_NodeLabel, intersect(b1_tree_wo_Neospl_NodeLabel$tip.label, Cichlids_chemoreception_df_ordered$Sp))
spp <- b1_tree_wo_Neospl_NodeLabel_Food$tip.label
corBM<-corBrownian(phy=b1_tree_wo_Neospl_NodeLabel_Food,form=~Sp)

Cichlids_chemoreception_df_ordered <- 
  Cichlids_chemoreception_df_ordered[match(spp, Cichlids_chemoreception_df_ordered$Sp), ]

ancova_Food.mean_lamellae <-
  gls(mean_lamellae ~ Food,
      data=Cichlids_chemoreception_df_ordered,
      correlation=corBM)
emmeans_result <- emmeans(ancova_Food.mean_lamellae, ~ Food)

pairwise.Food.mean_lamellae <- 
  as.data.frame(contrast(emmeans_result, method = "pairwise", adjust = "tukey"))

pairwise.Food.mean_lamellae <- 
  pairwise.Food.mean_lamellae %>%
  mutate(group1 = gsub(" -.*", "", contrast)) %>%
  mutate(group2 = gsub(".*- ", "", contrast)) %>%
  mutate(predictor = "Food") %>%
  mutate(response = "mean_lamellae")


mean.Food.mean_lamellae <- 
  Cichlids_chemoreception_df_ordered %>%
  group_by(Food) %>%
  summarise(group_mean = mean(mean_lamellae))
colnames(mean.Food.mean_lamellae) <- c("group1", "mean_group1")
pairwise.Food.mean_lamellae <- left_join(pairwise.Food.mean_lamellae, mean.Food.mean_lamellae, by="group1")
colnames(mean.Food.mean_lamellae) <- c("group2", "mean_group2")
pairwise.Food.mean_lamellae <- left_join(pairwise.Food.mean_lamellae, mean.Food.mean_lamellae, by="group2")


#Food x Delta_6

Chemoreception_pGLS_corr_df %>%
  filter(FDR_significant == "Significant") %>%
  filter(predictor == "Food") %>%
  pull(Response)

Cichlids_chemoreception_df_ordered <- 
  Cichlids_chemoreception_df %>% filter(! is.na(Food)) %>%
  filter(Sp %in% b1_tree_wo_Neospl_NodeLabel$tip.label)
b1_tree_wo_Neospl_NodeLabel_Food <- 
  keep.tip(b1_tree_wo_Neospl_NodeLabel, intersect(b1_tree_wo_Neospl_NodeLabel$tip.label, Cichlids_chemoreception_df_ordered$Sp))
spp <- b1_tree_wo_Neospl_NodeLabel_Food$tip.label
corBM<-corBrownian(phy=b1_tree_wo_Neospl_NodeLabel_Food,form=~Sp)

Cichlids_chemoreception_df_ordered <- 
  Cichlids_chemoreception_df_ordered[match(spp, Cichlids_chemoreception_df_ordered$Sp), ]


colnames(Cichlids_chemoreception_df_ordered)[colnames(Cichlids_chemoreception_df_ordered) == "Delta-6"] <- "Delta_6"


ancova_Food.Delta_6 <-
  gls(Delta_6 ~ Food,
      data=Cichlids_chemoreception_df_ordered,
      correlation=corBM)
emmeans_result <- emmeans(ancova_Food.Delta_6, ~ Food)

pairwise.Food.Delta_6 <- 
  as.data.frame(contrast(emmeans_result, method = "pairwise", adjust = "tukey"))

pairwise.Food.Delta_6 <- 
  pairwise.Food.Delta_6 %>%
  mutate(group1 = gsub(" -.*", "", contrast)) %>%
  mutate(group2 = gsub(".*- ", "", contrast)) %>%
  mutate(predictor = "Food") %>%
  mutate(response = "Delta_6")


mean.Food.Delta_6 <- 
  Cichlids_chemoreception_df_ordered %>%
  group_by(Food) %>%
  summarise(group_mean = mean(Delta_6))
colnames(mean.Food.Delta_6) <- c("group1", "mean_group1")
pairwise.Food.Delta_6 <- left_join(pairwise.Food.Delta_6, mean.Food.Delta_6, by="group1")
colnames(mean.Food.Delta_6) <- c("group2", "mean_group2")
pairwise.Food.Delta_6 <- left_join(pairwise.Food.Delta_6, mean.Food.Delta_6, by="group2")


#Food x Eta_10

Chemoreception_pGLS_corr_df %>%
  filter(FDR_significant == "Significant") %>%
  filter(predictor == "Food") %>%
  pull(Response)

Cichlids_chemoreception_df_ordered <- 
  Cichlids_chemoreception_df %>% filter(! is.na(Food)) %>%
  filter(Sp %in% b1_tree_wo_Neospl_NodeLabel$tip.label)
b1_tree_wo_Neospl_NodeLabel_Food <- 
  keep.tip(b1_tree_wo_Neospl_NodeLabel, intersect(b1_tree_wo_Neospl_NodeLabel$tip.label, Cichlids_chemoreception_df_ordered$Sp))
spp <- b1_tree_wo_Neospl_NodeLabel_Food$tip.label
corBM<-corBrownian(phy=b1_tree_wo_Neospl_NodeLabel_Food,form=~Sp)

Cichlids_chemoreception_df_ordered <- 
  Cichlids_chemoreception_df_ordered[match(spp, Cichlids_chemoreception_df_ordered$Sp), ]


colnames(Cichlids_chemoreception_df_ordered)[colnames(Cichlids_chemoreception_df_ordered) == "Eta-10"] <- "Eta_10"


ancova_Food.Eta_10 <-
  gls(Eta_10 ~ Food,
      data=Cichlids_chemoreception_df_ordered,
      correlation=corBM)
emmeans_result <- emmeans(ancova_Food.Eta_10, ~ Food)

pairwise.Food.Eta_10 <- 
  as.data.frame(contrast(emmeans_result, method = "pairwise", adjust = "tukey"))

pairwise.Food.Eta_10 <- 
  pairwise.Food.Eta_10 %>%
  mutate(group1 = gsub(" -.*", "", contrast)) %>%
  mutate(group2 = gsub(".*- ", "", contrast)) %>%
  mutate(predictor = "Food") %>%
  mutate(response = "Eta_10")


mean.Food.Eta_10 <- 
  Cichlids_chemoreception_df_ordered %>%
  group_by(Food) %>%
  summarise(group_mean = mean(Eta_10))
colnames(mean.Food.Eta_10) <- c("group1", "mean_group1")
pairwise.Food.Eta_10 <- left_join(pairwise.Food.Eta_10, mean.Food.Eta_10, by="group1")
colnames(mean.Food.Eta_10) <- c("group2", "mean_group2")
pairwise.Food.Eta_10 <- left_join(pairwise.Food.Eta_10, mean.Food.Eta_10, by="group2")


#Food x Eta_14

Chemoreception_pGLS_corr_df %>%
  filter(FDR_significant == "Significant") %>%
  filter(predictor == "Food") %>%
  pull(Response)

Cichlids_chemoreception_df_ordered <- 
  Cichlids_chemoreception_df %>% filter(! is.na(Food)) %>%
  filter(Sp %in% b1_tree_wo_Neospl_NodeLabel$tip.label)
b1_tree_wo_Neospl_NodeLabel_Food <- 
  keep.tip(b1_tree_wo_Neospl_NodeLabel, intersect(b1_tree_wo_Neospl_NodeLabel$tip.label, Cichlids_chemoreception_df_ordered$Sp))
spp <- b1_tree_wo_Neospl_NodeLabel_Food$tip.label
corBM<-corBrownian(phy=b1_tree_wo_Neospl_NodeLabel_Food,form=~Sp)

Cichlids_chemoreception_df_ordered <- 
  Cichlids_chemoreception_df_ordered[match(spp, Cichlids_chemoreception_df_ordered$Sp), ]


colnames(Cichlids_chemoreception_df_ordered)[colnames(Cichlids_chemoreception_df_ordered) == "Eta-14"] <- "Eta_14"


ancova_Food.Eta_14 <-
  gls(Eta_14 ~ Food,
      data=Cichlids_chemoreception_df_ordered,
      correlation=corBM)
emmeans_result <- emmeans(ancova_Food.Eta_14, ~ Food)

pairwise.Food.Eta_14 <- 
  as.data.frame(contrast(emmeans_result, method = "pairwise", adjust = "tukey"))

pairwise.Food.Eta_14 <- 
  pairwise.Food.Eta_14 %>%
  mutate(group1 = gsub(" -.*", "", contrast)) %>%
  mutate(group2 = gsub(".*- ", "", contrast)) %>%
  mutate(predictor = "Food") %>%
  mutate(response = "Eta_14")


mean.Food.Eta_14 <- 
  Cichlids_chemoreception_df_ordered %>%
  group_by(Food) %>%
  summarise(group_mean = mean(Eta_14))
colnames(mean.Food.Eta_14) <- c("group1", "mean_group1")
pairwise.Food.Eta_14 <- left_join(pairwise.Food.Eta_14, mean.Food.Eta_14, by="group1")
colnames(mean.Food.Eta_14) <- c("group2", "mean_group2")
pairwise.Food.Eta_14 <- left_join(pairwise.Food.Eta_14, mean.Food.Eta_14, by="group2")



#Food x Eta_8

Chemoreception_pGLS_corr_df %>%
  filter(FDR_significant == "Significant") %>%
  filter(predictor == "Food") %>%
  pull(Response)

Cichlids_chemoreception_df_ordered <- 
  Cichlids_chemoreception_df %>% filter(! is.na(Food)) %>%
  filter(Sp %in% b1_tree_wo_Neospl_NodeLabel$tip.label)
b1_tree_wo_Neospl_NodeLabel_Food <- 
  keep.tip(b1_tree_wo_Neospl_NodeLabel, intersect(b1_tree_wo_Neospl_NodeLabel$tip.label, Cichlids_chemoreception_df_ordered$Sp))
spp <- b1_tree_wo_Neospl_NodeLabel_Food$tip.label
corBM<-corBrownian(phy=b1_tree_wo_Neospl_NodeLabel_Food,form=~Sp)

Cichlids_chemoreception_df_ordered <- 
  Cichlids_chemoreception_df_ordered[match(spp, Cichlids_chemoreception_df_ordered$Sp), ]


colnames(Cichlids_chemoreception_df_ordered)[colnames(Cichlids_chemoreception_df_ordered) == "Eta-8"] <- "Eta_8"


ancova_Food.Eta_8 <-
  gls(Eta_8 ~ Food,
      data=Cichlids_chemoreception_df_ordered,
      correlation=corBM)
emmeans_result <- emmeans(ancova_Food.Eta_8, ~ Food)

pairwise.Food.Eta_8 <- 
  as.data.frame(contrast(emmeans_result, method = "pairwise", adjust = "tukey"))

pairwise.Food.Eta_8 <- 
  pairwise.Food.Eta_8 %>%
  mutate(group1 = gsub(" -.*", "", contrast)) %>%
  mutate(group2 = gsub(".*- ", "", contrast)) %>%
  mutate(predictor = "Food") %>%
  mutate(response = "Eta_8")


mean.Food.Eta_8 <- 
  Cichlids_chemoreception_df_ordered %>%
  group_by(Food) %>%
  summarise(group_mean = mean(Eta_8))
colnames(mean.Food.Eta_8) <- c("group1", "mean_group1")
pairwise.Food.Eta_8 <- left_join(pairwise.Food.Eta_8, mean.Food.Eta_8, by="group1")
colnames(mean.Food.Eta_8) <- c("group2", "mean_group2")
pairwise.Food.Eta_8 <- left_join(pairwise.Food.Eta_8, mean.Food.Eta_8, by="group2")



#Food x Eta_9

Chemoreception_pGLS_corr_df %>%
  filter(FDR_significant == "Significant") %>%
  filter(predictor == "Food") %>%
  pull(Response)

Cichlids_chemoreception_df_ordered <- 
  Cichlids_chemoreception_df %>% filter(! is.na(Food)) %>%
  filter(Sp %in% b1_tree_wo_Neospl_NodeLabel$tip.label)
b1_tree_wo_Neospl_NodeLabel_Food <- 
  keep.tip(b1_tree_wo_Neospl_NodeLabel, intersect(b1_tree_wo_Neospl_NodeLabel$tip.label, Cichlids_chemoreception_df_ordered$Sp))
spp <- b1_tree_wo_Neospl_NodeLabel_Food$tip.label
corBM<-corBrownian(phy=b1_tree_wo_Neospl_NodeLabel_Food,form=~Sp)

Cichlids_chemoreception_df_ordered <- 
  Cichlids_chemoreception_df_ordered[match(spp, Cichlids_chemoreception_df_ordered$Sp), ]


colnames(Cichlids_chemoreception_df_ordered)[colnames(Cichlids_chemoreception_df_ordered) == "Eta-9"] <- "Eta_9"


ancova_Food.Eta_9 <-
  gls(Eta_9 ~ Food,
      data=Cichlids_chemoreception_df_ordered,
      correlation=corBM)
emmeans_result <- emmeans(ancova_Food.Eta_9, ~ Food)

pairwise.Food.Eta_9 <- 
  as.data.frame(contrast(emmeans_result, method = "pairwise", adjust = "tukey"))

pairwise.Food.Eta_9 <- 
  pairwise.Food.Eta_9 %>%
  mutate(group1 = gsub(" -.*", "", contrast)) %>%
  mutate(group2 = gsub(".*- ", "", contrast)) %>%
  mutate(predictor = "Food") %>%
  mutate(response = "Eta_9")


mean.Food.Eta_9 <- 
  Cichlids_chemoreception_df_ordered %>%
  group_by(Food) %>%
  summarise(group_mean = mean(Eta_9))
colnames(mean.Food.Eta_9) <- c("group1", "mean_group1")
pairwise.Food.Eta_9 <- left_join(pairwise.Food.Eta_9, mean.Food.Eta_9, by="group1")
colnames(mean.Food.Eta_9) <- c("group2", "mean_group2")
pairwise.Food.Eta_9 <- left_join(pairwise.Food.Eta_9, mean.Food.Eta_9, by="group2")




#Food x Zeta_2

Chemoreception_pGLS_corr_df %>%
  filter(FDR_significant == "Significant") %>%
  filter(predictor == "Food") %>%
  pull(Response)

Cichlids_chemoreception_df_ordered <- 
  Cichlids_chemoreception_df %>% filter(! is.na(Food)) %>%
  filter(Sp %in% b1_tree_wo_Neospl_NodeLabel$tip.label)
b1_tree_wo_Neospl_NodeLabel_Food <- 
  keep.tip(b1_tree_wo_Neospl_NodeLabel, intersect(b1_tree_wo_Neospl_NodeLabel$tip.label, Cichlids_chemoreception_df_ordered$Sp))
spp <- b1_tree_wo_Neospl_NodeLabel_Food$tip.label
corBM<-corBrownian(phy=b1_tree_wo_Neospl_NodeLabel_Food,form=~Sp)

Cichlids_chemoreception_df_ordered <- 
  Cichlids_chemoreception_df_ordered[match(spp, Cichlids_chemoreception_df_ordered$Sp), ]


colnames(Cichlids_chemoreception_df_ordered)[colnames(Cichlids_chemoreception_df_ordered) == "Zeta-2"] <- "Zeta_2"


ancova_Food.Zeta_2 <-
  gls(Zeta_2 ~ Food,
      data=Cichlids_chemoreception_df_ordered,
      correlation=corBM)
emmeans_result <- emmeans(ancova_Food.Zeta_2, ~ Food)

pairwise.Food.Zeta_2 <- 
  as.data.frame(contrast(emmeans_result, method = "pairwise", adjust = "tukey"))

pairwise.Food.Zeta_2 <- 
  pairwise.Food.Zeta_2 %>%
  mutate(group1 = gsub(" -.*", "", contrast)) %>%
  mutate(group2 = gsub(".*- ", "", contrast)) %>%
  mutate(predictor = "Food") %>%
  mutate(response = "Zeta_2")


mean.Food.Zeta_2 <- 
  Cichlids_chemoreception_df_ordered %>%
  group_by(Food) %>%
  summarise(group_mean = mean(Zeta_2))
colnames(mean.Food.Zeta_2) <- c("group1", "mean_group1")
pairwise.Food.Zeta_2 <- left_join(pairwise.Food.Zeta_2, mean.Food.Zeta_2, by="group1")
colnames(mean.Food.Zeta_2) <- c("group2", "mean_group2")
pairwise.Food.Zeta_2 <- left_join(pairwise.Food.Zeta_2, mean.Food.Zeta_2, by="group2")




#Lets make a summarry of Posthoc tests ! 

pGLS_Food_habitat_posthoc <- 
  rbind(pairwise.habitat.Delta_10,
        pairwise.habitat.Eta_10,
        pairwise.habitat.Total_T1R,
        pairwise.habitat.T1R3,
        
        pairwise.Food.mean_lamellae,
        pairwise.Food.Delta_6,
        pairwise.Food.Eta_10,
        pairwise.Food.Eta_14,
        pairwise.Food.Eta_8,
        pairwise.Food.Eta_9,
        pairwise.Food.Zeta_2)


pGLS_Food_habitat_posthoc <- 
  pGLS_Food_habitat_posthoc %>%
  mutate(Signif = if_else(
    p.value < 0.05,
    "Significant",
    "Non_significant"
  ))
write.table(pGLS_Food_habitat_posthoc, 
            "pGLS_Food_habitat_posthoc.tsv",
            sep="\t",
            col.names=TRUE,
            row.names=FALSE,
            quote=FALSE)

pGLS_Food_habitat_posthoc %>% 
  filter(p.value < 0.05) %>%
  filter(predictor == "habitat")  %>%
  filter(response %in% c("Total_T1R", "T1R3"))
  
pGLS_Food_habitat_posthoc %>%
  filter(p.value < 0.05) %>%
  filter(predictor == "Food")  %>%
  pull(response) %>% unique()

pGLS_Food_habitat_posthoc %>%
  filter(p.value < 0.05) %>%
  filter(predictor == "Food") %>%
  filter(response == "Delta_6")

pGLS_Food_habitat_posthoc %>%
  filter(p.value < 0.05) %>%
  filter(predictor == "Food") %>%
  filter(response == "Eta_10")

pGLS_Food_habitat_posthoc %>%
  filter(p.value < 0.05) %>%
  filter(predictor == "Food") %>%
  filter(response == "Eta_14")


pGLS_Food_habitat_posthoc %>%
  filter(p.value < 0.05) %>%
  filter(predictor == "Food") %>%
  filter(response == "Eta_8")


pGLS_Food_habitat_posthoc %>%
  filter(p.value < 0.05) %>%
  filter(predictor == "Food") %>%
  filter(response == "Eta_9")

Cichlids_chemoreception_df_ordered$Food %>% unique()

pGLS_Food_habitat_posthoc %>%
  filter(p.value < 0.05) %>%
  filter(predictor == "Food") %>%
  filter(response == "mean_lamellae")

#### pGLS All detailed plots  -------------------


significant_correlations <- 
  Chemoreception_pGLS_corr_df %>%
  filter(FDR_significant == "Significant")

tribes.colors <- c(Bathybatini = "#242626" , Benthochromini = "#AE262A" , Boulengerochromini = "#59595C" , Cyphotilapiini = "#FDDF13" , Cyprichromini = "#F04D29" , Ectodini = "#9AB9D9" , Eretmodini = "#682E7A" , Lamprologini = "#C588BB" , Limnochromini = "#535CA9" , Perissodini = "#FBAC43" , Trematocarini = "#959170" , Tropheini = "#86C773" , Haplochromini = "#274e13" , outgroup = "gray", Other="gray")

### First graphs with mean SL

significant_correlation_subset <- 
  significant_correlations %>%
  filter(predictor == "mean_SL") 


plot_list <- list()
all_plot_names <- c()
for (row in 1:nrow(significant_correlation_subset)){
  
  curr_lambda <- significant_correlation_subset[row, "lambda"]
  curr_R2  <- significant_correlation_subset[row, "R2"]
  curr_response  <- significant_correlation_subset[row, "Response"]
  curr_predictor  <- significant_correlation_subset[row, "predictor"]
  curr_pvalue <- significant_correlation_subset[row, "pvalue"]

  curr_chemo_data <- 
    Cichlids_chemoreception_df %>%
    filter(Sp %in% b1_tree_wo_Neospl_NodeLabel$tip.label) %>%
    dplyr::select(Sp, curr_predictor, curr_response, Tribe)
  
  colnames(curr_chemo_data) <- c("Sp", "predictor", "response", "Tribe")
  curr_chemo_data <- 
    curr_chemo_data %>% filter(! is.na(predictor)) %>% filter(! is.na(response))
  

  lm_fit <- lm(data = curr_chemo_data, 
               formula = response ~ predictor)
  sum_fit <- summary(lm_fit)
  GLS_r2 =      formatC( sum_fit$r.squared, digits = 3)
  GLS_pente =   formatC(sum_fit$coefficients[2], digits = 3)
  GLS_pvalue =  formatC(sum_fit$coefficients[8], digits = 3)
  GLS_cc <- coef(lm_fit)
  GLS_fit_function <- function(x) GLS_cc[1] + GLS_cc[2]*x
  
    
  p = curr_chemo_data %>%
    ggplot(., aes(x=predictor, y=response, color=Tribe)) +
    geom_point(size=2) + 
    theme_classic() +
    #stat_function(fun = GLS_fit_function, color="black") + 
    geom_smooth(method = "lm", formula = y ~ x, color="black", se = FALSE, data = curr_chemo_data, aes(group = 1)) + 
    labs(subtitle = bquote(R^2 ~ "=" ~ .(curr_R2) ~ ";" ~ P ~ "=" ~.(curr_pvalue) ~ ";" ~ lambda ~ "=" ~ .(curr_lambda))) +
    xlab(toString(curr_predictor)) +
    ylab(toString(curr_response)) +
    theme(axis.text=element_text(size=16),
          axis.title=element_text(size=18),
          plot.subtitle=element_text(size=16),
          legend.position="none") +
    scale_color_manual(values=tribes.colors) 
  
  curr_list_name <- paste(curr_response, curr_predictor, sep="_")
  
  plot_list[[curr_list_name]] = p
  all_plot_names <- c(all_plot_names, curr_list_name)
}


pdf("Raw_R_plots/Signif_pGLS_SL.pdf",width = 8.34,  height = 4.61)
for (curr_plot in all_plot_names) {
  print(plot_list[[curr_plot]])
}
dev.off()


### Graph with mean d15N

significant_correlation_subset <- 
  significant_correlations %>%
  filter(predictor == "mean_d15N") 


plot_list <- list()
all_plot_names <- c()
for (row in 1:nrow(significant_correlation_subset)){
  
  curr_lambda <- significant_correlation_subset[row, "lambda"]
  curr_R2  <- significant_correlation_subset[row, "R2"]
  curr_response  <- significant_correlation_subset[row, "Response"]
  curr_predictor  <- significant_correlation_subset[row, "predictor"]
  curr_pvalue <- significant_correlation_subset[row, "pvalue"]
  
  curr_chemo_data <- 
    Cichlids_chemoreception_df %>%
    filter(Sp %in% b1_tree_wo_Neospl_NodeLabel$tip.label) %>%
    dplyr::select(Sp, curr_predictor, curr_response, Tribe)
  
  colnames(curr_chemo_data) <- c("Sp", "predictor", "response", "Tribe")
  curr_chemo_data <- 
    curr_chemo_data %>% filter(! is.na(predictor)) %>% filter(! is.na(response))
  
  
  lm_fit <- lm(data = curr_chemo_data, 
               formula = response ~ predictor)
  sum_fit <- summary(lm_fit)
  GLS_r2 =      formatC( sum_fit$r.squared, digits = 3)
  GLS_pente =   formatC(sum_fit$coefficients[2], digits = 3)
  GLS_pvalue =  formatC(sum_fit$coefficients[8], digits = 3)
  GLS_cc <- coef(lm_fit)
  GLS_fit_function <- function(x) GLS_cc[1] + GLS_cc[2]*x
  
  
  p = curr_chemo_data %>%
    ggplot(., aes(x=predictor, y=response, color=Tribe)) +
    geom_point(size=2) + 
    theme_classic() +
    #stat_function(fun = GLS_fit_function, color="black") + 
    geom_smooth(method = "lm", formula = y ~ x, color="black", se = FALSE, data = curr_chemo_data, aes(group = 1)) + 
    labs(subtitle = bquote(R^2 ~ "=" ~ .(curr_R2) ~ ";" ~ P ~ "=" ~.(curr_pvalue) ~ ";" ~ lambda ~ "=" ~ .(curr_lambda))) +
    xlab(toString(curr_predictor)) +
    ylab(toString(curr_response)) +
    theme(axis.text=element_text(size=16),
          axis.title=element_text(size=18),
          plot.subtitle=element_text(size=16),
          legend.position="none") +
    scale_color_manual(values=tribes.colors) 
  
  curr_list_name <- paste(curr_response, curr_predictor, sep="_")
  
  plot_list[[curr_list_name]] = p
  all_plot_names <- c(all_plot_names, curr_list_name)
}


pdf("Raw_R_plots/Signif_pGLS_d15N.pdf",width = 8.34,  height = 4.61)
for (curr_plot in all_plot_names) {
  print(plot_list[[curr_plot]])
}
dev.off()


### Graph with exploratory tendency

significant_correlation_subset <- 
  significant_correlations %>%
  filter(predictor == "median_exploration") 


plot_list <- list()
all_plot_names <- c()
for (row in 1:nrow(significant_correlation_subset)){
  
  curr_lambda <- significant_correlation_subset[row, "lambda"]
  curr_R2  <- significant_correlation_subset[row, "R2"]
  curr_response  <- significant_correlation_subset[row, "Response"]
  curr_predictor  <- significant_correlation_subset[row, "predictor"]
  curr_pvalue <- significant_correlation_subset[row, "pvalue"]
  
  curr_chemo_data <- 
    Cichlids_chemoreception_df %>%
    filter(Sp %in% b1_tree_wo_Neospl_NodeLabel$tip.label) %>%
    dplyr::select(Sp, curr_predictor, curr_response, Tribe)
  
  colnames(curr_chemo_data) <- c("Sp", "predictor", "response", "Tribe")
  curr_chemo_data <- 
    curr_chemo_data %>% filter(! is.na(predictor)) %>% filter(! is.na(response))
  
  
  lm_fit <- lm(data = curr_chemo_data, 
               formula = response ~ predictor)
  sum_fit <- summary(lm_fit)
  GLS_r2 =      formatC( sum_fit$r.squared, digits = 3)
  GLS_pente =   formatC(sum_fit$coefficients[2], digits = 3)
  GLS_pvalue =  formatC(sum_fit$coefficients[8], digits = 3)
  GLS_cc <- coef(lm_fit)
  GLS_fit_function <- function(x) GLS_cc[1] + GLS_cc[2]*x
  
  
  p = curr_chemo_data %>%
    ggplot(., aes(x=predictor, y=response, color=Tribe)) +
    geom_point(size=2) + 
    theme_classic() +
    #stat_function(fun = GLS_fit_function, color="black") + 
    geom_smooth(method = "lm", formula = y ~ x, color="black", se = FALSE, data = curr_chemo_data, aes(group = 1)) + 
    labs(subtitle = bquote(R^2 ~ "=" ~ .(curr_R2) ~ ";" ~ P ~ "=" ~.(curr_pvalue) ~ ";" ~ lambda ~ "=" ~ .(curr_lambda))) +
    xlab(toString(curr_predictor)) +
    ylab(toString(curr_response)) +
    theme(axis.text=element_text(size=16),
          axis.title=element_text(size=18),
          plot.subtitle=element_text(size=16),
          legend.position="none") +
    scale_color_manual(values=tribes.colors) 
  
  curr_list_name <- paste(curr_response, curr_predictor, sep="_")
  
  plot_list[[curr_list_name]] = p
  all_plot_names <- c(all_plot_names, curr_list_name)
}


pdf("Raw_R_plots/Signif_pGLS_Explo.pdf",width = 8.34,  height = 4.61)
for (curr_plot in all_plot_names) {
  print(plot_list[[curr_plot]])
}
dev.off()

### Graph with mean habitat

significant_correlation_subset <- 
  significant_correlations %>%
  filter(predictor == "habitat") 


plot_list <- list()
all_plot_names <- c()
for (row in 1:nrow(significant_correlation_subset)){
  
  curr_lambda <- significant_correlation_subset[row, "lambda"]
  curr_R2  <- significant_correlation_subset[row, "R2"]
  curr_response  <- significant_correlation_subset[row, "Response"]
  curr_predictor  <- significant_correlation_subset[row, "predictor"]
  curr_pvalue <- significant_correlation_subset[row, "pvalue"]
  
  curr_chemo_data <- 
    Cichlids_chemoreception_df %>%
    filter(Sp %in% b1_tree_wo_Neospl_NodeLabel$tip.label) %>%
    dplyr::select(Sp, curr_predictor, curr_response, Tribe)
  
  colnames(curr_chemo_data) <- c("Sp", "predictor", "response", "Tribe")
  curr_chemo_data <- 
    curr_chemo_data %>% filter(! is.na(predictor)) %>% filter(! is.na(response))
  
  curr_chemo_data$predictor <-
    factor(curr_chemo_data$predictor ,
           levels=c("shallow", "litoral","inter", "pelagic", "deep"))
  
  
  p = curr_chemo_data %>%
    ggplot(., aes(x=predictor, y=response)) +
    geom_violin() + 
    geom_jitter(aes(color = Tribe, alpha=0.8), position=position_jitter(0.2)) +
    theme_classic() +
    labs(subtitle = bquote(R^2 ~ "=" ~ .(curr_R2) ~ ";" ~ P ~ "=" ~.(curr_pvalue) ~ ";" ~ lambda ~ "=" ~ .(curr_lambda))) +
    xlab(toString(curr_predictor)) +
    ylab(toString(curr_response)) +
    theme(axis.text=element_text(size=16),
          axis.title=element_text(size=18),
          plot.subtitle=element_text(size=16),
          legend.position="none") +
    scale_color_manual(values=tribes.colors) 
  

  curr_list_name <- paste(curr_response, curr_predictor, sep="_")
  
  plot_list[[curr_list_name]] = p
  all_plot_names <- c(all_plot_names, curr_list_name)
}


pdf("Raw_R_plots/Signif_pGLS_habitat.pdf",width = 8.34,  height = 4.61)
for (curr_plot in all_plot_names) {
  print(plot_list[[curr_plot]])
}
dev.off()





### Graph with mean Food

significant_correlation_subset <- 
  significant_correlations %>%
  filter(predictor == "Food") 


plot_list <- list()
all_plot_names <- c()
for (row in 1:nrow(significant_correlation_subset)){
  
  curr_lambda <- significant_correlation_subset[row, "lambda"]
  curr_R2  <- significant_correlation_subset[row, "R2"]
  curr_response  <- significant_correlation_subset[row, "Response"]
  curr_predictor  <- significant_correlation_subset[row, "predictor"]
  curr_pvalue <- significant_correlation_subset[row, "pvalue"]
  
  curr_chemo_data <- 
    Cichlids_chemoreception_df %>%
    filter(Sp %in% b1_tree_wo_Neospl_NodeLabel$tip.label) %>%
    dplyr::select(Sp, curr_predictor, curr_response, Tribe)
  
  colnames(curr_chemo_data) <- c("Sp", "predictor", "response", "Tribe")
  curr_chemo_data <- 
    curr_chemo_data %>% filter(! is.na(predictor)) %>% filter(! is.na(response))
  
  
  curr_chemo_data$predictor <-
    factor(curr_chemo_data$predictor ,
           levels=c("herb", "Auf_herb", "plan","plan_invert", "invert", "Auf_invert",
                    "fry_plan", "scales", "fish_invert", "fish","omni"))
  
  
  p = curr_chemo_data %>%
    ggplot(., aes(x=predictor, y=response)) +
    geom_violin() + 
    geom_jitter(aes(color = Tribe, alpha=0.8), position=position_jitter(0.2)) +
    theme_classic() +
    labs(subtitle = bquote(R^2 ~ "=" ~ .(curr_R2) ~ ";" ~ P ~ "=" ~.(curr_pvalue) ~ ";" ~ lambda ~ "=" ~ .(curr_lambda))) +
    xlab(toString(curr_predictor)) +
    ylab(toString(curr_response)) +
    theme(axis.text=element_text(size=16),
          axis.title=element_text(size=18),
          plot.subtitle=element_text(size=16),
          legend.position="none") +
    scale_color_manual(values=tribes.colors) 
  
  
  curr_list_name <- paste(curr_response, curr_predictor, sep="_")
  
  plot_list[[curr_list_name]] = p
  all_plot_names <- c(all_plot_names, curr_list_name)
}


pdf("Raw_R_plots/Signif_pGLS_Food.pdf",width = 8.34,  height = 4.61)
for (curr_plot in all_plot_names) {
  print(plot_list[[curr_plot]])
}
dev.off()




### Graph with breeding mode

significant_correlation_subset <- 
  significant_correlations %>%
  filter(predictor == "breeding_mode") 


plot_list <- list()
all_plot_names <- c()
for (row in 1:nrow(significant_correlation_subset)){
  
  curr_lambda <- significant_correlation_subset[row, "lambda"]
  curr_R2  <- significant_correlation_subset[row, "R2"]
  curr_response  <- significant_correlation_subset[row, "Response"]
  curr_predictor  <- significant_correlation_subset[row, "predictor"]
  curr_pvalue <- significant_correlation_subset[row, "pvalue"]
  
  curr_chemo_data <- 
    Cichlids_chemoreception_df %>%
    filter(Sp %in% b1_tree_wo_Neospl_NodeLabel$tip.label) %>%
    dplyr::select(Sp, curr_predictor, curr_response, Tribe)
  
  colnames(curr_chemo_data) <- c("Sp", "predictor", "response", "Tribe")
  curr_chemo_data <- 
    curr_chemo_data %>% filter(! is.na(predictor)) %>% filter(! is.na(response))
  
  
  curr_chemo_data$predictor <-
    factor(curr_chemo_data$predictor ,
           levels=c("maternal bower brooding", "maternal territorial brooding", "maternal open water brooding",
                    "biparental brooding","cave spawning", "shell spawning", "open spawning"))
  
  
  p = curr_chemo_data %>%
    ggplot(., aes(x=predictor, y=response)) +
    geom_violin() + 
    geom_jitter(aes(color = Tribe, alpha=0.8), position=position_jitter(0.2)) +
    theme_classic() +
    labs(subtitle = bquote(R^2 ~ "=" ~ .(curr_R2) ~ ";" ~ P ~ "=" ~.(curr_pvalue) ~ ";" ~ lambda ~ "=" ~ .(curr_lambda))) +
    xlab(toString(curr_predictor)) +
    ylab(toString(curr_response)) +
    theme(axis.text=element_text(size=16),
          axis.title=element_text(size=18),
          plot.subtitle=element_text(size=16),
          legend.position="none") +
    scale_color_manual(values=tribes.colors) 
  
  
  curr_list_name <- paste(curr_response, curr_predictor, sep="_")
  
  plot_list[[curr_list_name]] = p
  all_plot_names <- c(all_plot_names, curr_list_name)
}


pdf("Raw_R_plots/Signif_pGLS_BM.pdf",width = 8.34,  height = 4.61)
for (curr_plot in all_plot_names) {
  print(plot_list[[curr_plot]])
}
dev.off()



### Graph with breeding type

significant_correlation_subset <- 
  significant_correlations %>%
  filter(predictor == "breeding_type") 


plot_list <- list()
all_plot_names <- c()
for (row in 1:nrow(significant_correlation_subset)){
  
  curr_lambda <- significant_correlation_subset[row, "lambda"]
  curr_R2  <- significant_correlation_subset[row, "R2"]
  curr_response  <- significant_correlation_subset[row, "Response"]
  curr_predictor  <- significant_correlation_subset[row, "predictor"]
  curr_pvalue <- significant_correlation_subset[row, "pvalue"]
  
  curr_chemo_data <- 
    Cichlids_chemoreception_df %>%
    filter(Sp %in% b1_tree_wo_Neospl_NodeLabel$tip.label) %>%
    dplyr::select(Sp, curr_predictor, curr_response, Tribe)
  
  colnames(curr_chemo_data) <- c("Sp", "predictor", "response", "Tribe")
  curr_chemo_data <- 
    curr_chemo_data %>% filter(! is.na(predictor)) %>% filter(! is.na(response))
  
  
  curr_chemo_data$predictor <-
    factor(curr_chemo_data$predictor ,
           levels=c("Mouthbrooding", "Substrate spawning"))
  
  
  p = curr_chemo_data %>%
    ggplot(., aes(x=predictor, y=response)) +
    geom_violin() + 
    geom_jitter(aes(color = Tribe, alpha=0.8), position=position_jitter(0.2)) +
    theme_classic() +
    labs(subtitle = bquote(R^2 ~ "=" ~ .(curr_R2) ~ ";" ~ P ~ "=" ~.(curr_pvalue) ~ ";" ~ lambda ~ "=" ~ .(curr_lambda))) +
    xlab(toString(curr_predictor)) +
    ylab(toString(curr_response)) +
    theme(axis.text=element_text(size=16),
          axis.title=element_text(size=18),
          plot.subtitle=element_text(size=16),
          legend.position="none") +
    scale_color_manual(values=tribes.colors) 
  
  
  curr_list_name <- paste(curr_response, curr_predictor, sep="_")
  
  plot_list[[curr_list_name]] = p
  all_plot_names <- c(all_plot_names, curr_list_name)
}


pdf("Raw_R_plots/Signif_pGLS_BT.pdf",width = 8.34,  height = 4.61)
for (curr_plot in all_plot_names) {
  print(plot_list[[curr_plot]])
}
dev.off()



#### OLD - Multiple pGLS with the number of lamellae  -------------------



caper_data_chemoreception_noNA <- 
  comparative.data(
    phy = b1_tree_wo_Neospl_NodeLabel,
    data = Cichlids_chemoreception_df %>% filter(! is.na(Food)) %>% 
      filter(! is.na(breeding_mode)) %>% filter(! is.na(mean_SL)),
                   names.col = Sp, vcv = TRUE,
                   na.omit = FALSE, warn.dropped = TRUE)





Lamellae_vs_SL <-
  pgls(mean_lamellae ~ mean_SL, 
       data = caper_data_chemoreception_noNA , 
       lambda = "ML",
       bounds=list(lambda=c(0,1))
  )


Lamellae_vs_Food_SL <-
  pgls(mean_lamellae ~ Food + mean_SL, 
       data = caper_data_chemoreception_noNA , 
       lambda = "ML",
       bounds=list(lambda=c(0,1))
  )

Lamellae_vs_BM_SL <-
  pgls(mean_lamellae ~ breeding_mode + mean_SL, 
       data = caper_data_chemoreception_noNA , 
       lambda = "ML",
       bounds=list(lambda=c(0,1))
  )

Lamellae_vs_Food_BM_SL <-
  pgls(mean_lamellae ~ Food + breeding_mode + mean_SL, 
       data = caper_data_chemoreception_noNA , 
       lambda = "ML",
       bounds=list(lambda=c(0,1))
  )
       
       
       
Lamellae_AIC_df <- 
  AIC(Lamellae_vs_SL, Lamellae_vs_Food_SL, Lamellae_vs_BM_SL, Lamellae_vs_Food_BM_SL)





caper_data_chemoreception_noNA <- 
  comparative.data(
    phy = b1_tree_wo_Neospl_NodeLabel,
    data = Cichlids_chemoreception_df %>% filter(! is.na(Food)) %>% 
      filter(! is.na(breeding_mode)) %>% filter(! is.na(mean_SL))  %>% filter(! is.na(median_exploration)),
    names.col = Sp, vcv = TRUE,
    na.omit = FALSE, warn.dropped = TRUE)


Lamellae_vs_Exploratory <-
  pgls(mean_lamellae ~ median_exploration, 
       data = caper_data_chemoreception_noNA , 
       lambda = "ML",
       bounds=list(lambda=c(0,1))
  )
summary(Lamellae_vs_Exploratory)


Lamellae_vs_SL <-
  pgls(mean_lamellae ~ mean_SL, 
       data = caper_data_chemoreception_noNA , 
       lambda = "ML",
       bounds=list(lambda=c(0,1))
  )
summary(Lamellae_vs_Exploratory)

Lamellae_vs_Exploratory_SL <-
  pgls(mean_lamellae ~ median_exploration + mean_SL, 
       data = caper_data_chemoreception_noNA , 
       lambda = "ML",
       bounds=list(lambda=c(0,1))
  )
summary(Lamellae_vs_Exploratory_SL)



Lamellae_AIC_df <- 
  AIC(Lamellae_vs_SL, Lamellae_vs_Exploratory_SL)


#### NEW - Multiple pGLS with the number of lamellae  -------------------


caper_data_chemoreception_noNA <- 
  comparative.data(
    phy = b1_tree_wo_Neospl_NodeLabel,
    data = Cichlids_chemoreception_df %>% filter(! is.na(Food)) %>% 
      filter(! is.na(mean_SL))  %>% filter(! is.na(median_exploration)),
    names.col = Sp, vcv = TRUE,
    na.omit = FALSE, warn.dropped = TRUE)


Lamellae_vs_Exploratory <-
  pgls(mean_lamellae ~ median_exploration, 
       data = caper_data_chemoreception_noNA , 
       lambda = "ML",
       bounds=list(lambda=c(0,1))
  )
Lamellae_vs_SL <-
  pgls(mean_lamellae ~ mean_SL, 
       data = caper_data_chemoreception_noNA , 
       lambda = "ML",
       bounds=list(lambda=c(0,1))
  )
Lamellae_vs_Food <-
  pgls(mean_lamellae ~ Food, 
       data = caper_data_chemoreception_noNA , 
       lambda = "ML",
       bounds=list(lambda=c(0,1))
  )
Lamellae_vs_All <-
  pgls(mean_lamellae ~ Food + median_exploration + mean_SL, 
       data = caper_data_chemoreception_noNA , 
       lambda = "ML",
       bounds=list(lambda=c(0,1))
  )

Lamellae_AIC_df <- 
  AIC(Lamellae_vs_SL, Lamellae_vs_Exploratory, Lamellae_vs_Food, Lamellae_vs_All)



#### pGLS between brain regions and ecological factors  -------------------


brain_responses <- c("prop_Olfactory_bulb","prop_Telencephalon","prop_Optic_tectum", "prop_Cerebellum",
                     "prop_Rest_of_the_brain")
  
numerical_predictor <- 
  c("mean_SL", "mean_weight", "mean_d15N", "mean_d13C", "median_exploration")


abnorm_termination <- c()

Brains_numerical_pGLS_df <- as.data.frame(NULL)
for(curr_response in brain_responses){
  for(curr_predictor in numerical_predictor){
    
    current_test_name <- paste(curr_response, curr_predictor, sep="_")
    print(current_test_name)
    
    curr_data_chemoreception <- 
      Cichlids_chemoreception_df %>%
      dplyr::select(Sp, curr_response, curr_predictor)
    
    
    colnames(curr_data_chemoreception) <- c("Sp", "reponse", "predictor")
    
    curr_caper_data <- 
      comparative.data(phy = b1_tree_wo_Neospl_NodeLabel, 
                       data = curr_data_chemoreception,
                       names.col = Sp, vcv = TRUE,
                       na.omit = FALSE, warn.dropped = TRUE)
    
    
    if(current_test_name %in% abnorm_termination){ #reduce the lambda search boundary for failed computations
      curr_pgls_result <-
        pgls(reponse ~ predictor, 
             data = curr_caper_data, 
             lambda = "ML",
             bounds=list(lambda=c(0.7,1))
        )
    } else {
      curr_pgls_result <-
        pgls(reponse ~ predictor, 
             data = curr_caper_data, 
             lambda = "ML",
             bounds=list(lambda=c(0,1))
        )
    }
    
    
    
    
    sum_fit_phy <- summary(curr_pgls_result)
    PGLS_pente =   formatC(sum_fit_phy$coefficients[2], digits = 3)
    PGLS_r2 =      formatC(sum_fit_phy$adj.r.squared, digits = 2)
    PGLS_pvalue =  formatC(sum_fit_phy$coefficients[8], digits = 3)
    if (PGLS_pvalue == "   0"){ PGLS_pvalue = "< 2.2e-16"}
    lambda <- as.numeric(sum_fit_phy$param[2])
    slope <- formatC(sum_fit_phy$coefficients[2], digits = 3)
    
    print(PGLS_pvalue)
    
    curr_df <- 
      as.data.frame(
        cbind(
          curr_response,
          PGLS_r2,
          PGLS_pvalue,
          lambda,
          slope,
          curr_predictor
        )
      )
    
    colnames(curr_df) <- c("Response", "R2", "pvalue", "lambda", "slope","predictor")
    
    Brains_numerical_pGLS_df <- 
      rbind(Brains_numerical_pGLS_df, curr_df)
  }
  
}


#Now with categorical predictors
categorical_predictors <- 
  c("Food", "habitat", "breeding_type", "breeding_mode")

abnorm_termination <- 
  c()


Brains_categorical_pGLS_df <- as.data.frame(NULL)
for(curr_response in brain_responses){
  for(curr_predictor in categorical_predictors){
    
    
    current_test_name <- paste(curr_response, curr_predictor, sep="_")
    print(current_test_name)
    
    curr_data_chemoreception <- 
      Cichlids_chemoreception_df %>%
      dplyr::select(Sp, curr_response, curr_predictor) %>%
      filter(! is.na(curr_predictor))
    colnames(curr_data_chemoreception) <- c("Sp", "reponse", "predictor")
    
    
    
    curr_data_caper <- 
      comparative.data(phy = b1_tree_wo_Neospl_NodeLabel, 
                       data = curr_data_chemoreception,
                       names.col = Sp, vcv = TRUE,
                       na.omit = FALSE, warn.dropped = TRUE)
    
    
    
    if(current_test_name %in% abnorm_termination){ #reduce the lambda search boundary for failed computations
      curr_pgls_result <-
        pgls(reponse ~ predictor, 
             data = curr_data_caper, 
             lambda = "ML",
             bounds=list(lambda=c(0.5,1)))
    } else {
      curr_pgls_result <-
        pgls(reponse ~ predictor, 
             data = curr_data_caper, 
             lambda = "ML",
             bounds=list(lambda=c(0,1)))
    }
    
    
    
    
    curr_pgls_summary <- summary(curr_pgls_result)
    
    pGLS_R2 <- as.numeric(curr_pgls_summary$adj.r.squared)
    pGLS_lambda <- as.numeric(curr_pgls_summary$param[2])
    
    pGLS_pvalue <- pf(curr_pgls_summary$fstatistic[1], 
                      curr_pgls_summary$fstatistic[2], 
                      curr_pgls_summary$fstatistic[3], 
                      lower.tail = FALSE)
    
    
    curr_df <- 
      as.data.frame(
        cbind(
          curr_response,
          pGLS_R2,
          pGLS_pvalue,
          pGLS_lambda,
          curr_predictor
        )
      )
    
    colnames(curr_df) <- c("Response", "R2", "pvalue", "lambda","predictor")
    
    
    Brains_categorical_pGLS_df <- 
      rbind(Brains_categorical_pGLS_df,
            curr_df)
    
    
  }
  
}





Brains_categorical_pGLS_df <- Brains_categorical_pGLS_df %>%  mutate(slope = NA)

Brains_categorical_pGLS_df <- 
  Brains_categorical_pGLS_df %>%
  dplyr::select(Response,R2, pvalue, lambda, slope, predictor)

colnames(Brains_categorical_pGLS_df) <- colnames(Brains_numerical_pGLS_df)


#Merge tables

Brains_pGLS_df <- 
  rbind(Brains_numerical_pGLS_df, Brains_categorical_pGLS_df)


#Summarise the data
as.data.frame(
  Brains_pGLS_df %>%
    group_by(Response) %>%
    summarise(n()))

as.data.frame(
  Brains_pGLS_df %>%
    group_by(predictor) %>%
    summarise(n()))


Brains_pGLS_df <-
  Brains_pGLS_df %>%
  filter(predictor != "mean_weight")

Brains_pGLS_df$pvalue <- as.numeric(Brains_pGLS_df$pvalue)
Brains_pGLS_df$R2 <- as.numeric(Brains_pGLS_df$R2)

#Add the FDR and Bonferroni p-values 

unique_response <- Brains_pGLS_df$Response %>% unique()
Brains_pGLS_corr_df <- as.data.frame(NULL)
for(curr_response in unique_response){
  
  curr_df <- 
    Brains_pGLS_df %>%
    filter(Response == curr_response)
  
  
  initial_pvalues <- 
    as.numeric(
      curr_df %>% pull(pvalue)
    )
  
  corrected_pvalues_FDR <- p.adjust(initial_pvalues, method = "fdr")
  corrected_pvalues_bonferroni <- p.adjust(initial_pvalues, method = "bonferroni")
  
  curr_df <- 
    curr_df %>%
    mutate(FDR_pvalue = corrected_pvalues_FDR,
           Bonferroni_pvalue = corrected_pvalues_bonferroni)
  
  Brains_pGLS_corr_df <- 
    rbind(Brains_pGLS_corr_df, curr_df)
  
  
}


## Define significance levels with the difference corrections

Brains_pGLS_corr_df <- 
  Brains_pGLS_corr_df %>%
  mutate(FDR_significant = if_else(
    as.numeric(FDR_pvalue) < 0.05,
    "Significant", 
    "N.S"
  ))

Brains_pGLS_corr_df <- 
  Brains_pGLS_corr_df %>%
  mutate(Bonferroni_significant = if_else(
    as.numeric(Bonferroni_pvalue) < 0.05,
    "Significant", 
    "N.S"
  ))

Brains_pGLS_corr_df <- 
  Brains_pGLS_corr_df %>%
  mutate(significant = if_else(
    as.numeric(pvalue) < 0.05,
    "Significant", 
    "N.S"
  ))


## Define the slope ...
Brains_pGLS_corr_df <- 
  Brains_pGLS_corr_df %>%
  mutate(slope_sign = case_when(
    slope <= 0 ~ "minus",
    slope > 0 ~ "plus",
    is.na(slope) ~ "cat",
  ))


## Define a shape for slopes for future graphics
slope_sign_shapes <- 
  c("minus"="\u25BC",
    "plus"="\u25B2",
    "cat"="\u25FE")

#Lets make a summary graphic 

Brains_pGLS_corr_df$Response <-
  factor(Brains_pGLS_corr_df$Response ,
         levels=c("prop_Olfactory_bulb", "prop_Telencephalon", "prop_Optic_tectum", 
                  "prop_Cerebellum", "prop_Rest_of_the_brain"))


Brains_pGLS_corr_df$predictor <-
  factor(Brains_pGLS_corr_df$predictor ,
         levels=c("mean_SL", 
                  "median_exploration", 
                  "mean_d15N", "Food",
                  "mean_d13C", "habitat",
                  "breeding_type", "breeding_mode"))



Signif_colors <-
  c("N.S" = 0,
    "Significant" = 1)



#Replace negative R2 by R2 of 0


Brains_pGLS_corr_df$R2[Brains_pGLS_corr_df$R2 < 0] <- 0


#Lets plot
pdf(file = "Raw_R_plots/Brains_pGLS_FDR.pdf", width = 10.34,  height = 4.61)

ggplot(Brains_pGLS_corr_df, 
       aes(y=predictor, x=Response, fill= R2, alpha=FDR_significant, color="black")) + 
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

ggplot(Brains_pGLS_corr_df, 
       aes(y=predictor, x=Response, fill= R2, alpha=FDR_significant, color="black")) + 
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




write.table(Brains_pGLS_corr_df,
            "Raw_R_plots/Brains_pGLS_FDR.tsv",
            sep="\t",
            col.names=TRUE,
            row.names=FALSE,
            quote=FALSE)


#lets plots the details

signif_Brains_pGLS_corr_df <- 
  Brains_pGLS_corr_df %>% filter(FDR_pvalue < 0.05)


pdf(file = "Raw_R_plots/Brains_pGLS_details.pdf", width = 8.34,  height = 4.61)

Cichlids_chemoreception_df %>%
  filter(! is.na(prop_Optic_tectum)) %>%
  ggplot(., aes(x=mean_SL, y=prop_Optic_tectum, color=Tribe)) +
  geom_point(size=2) + 
  theme_classic() +
  #stat_function(fun = GLS_fit_function, color="black") + 
  geom_smooth(method = "lm", formula = y ~ x, color="black", se = FALSE, data = Cichlids_chemoreception_df, aes(group = 1)) + 
  labs(subtitle = bquote(R^2 ~ "=" ~ .(0.5000000) ~ ";" ~ P ~ "=" ~.(0.004330000) ~ ";" ~ lambda ~ "=" ~ .(0))) +
  xlab("mean SL") +
  ylab("prop_Optic_tectum") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none") +
  scale_color_manual(values=tribes.colors)


Cichlids_chemoreception_df %>%
  ggplot(., aes(x=breeding_type, y=prop_Telencephalon)) +
  geom_violin() + 
  geom_jitter(aes(color = Tribe, alpha=0.8), position=position_jitter(0.2)) +
  theme_classic() +
  labs(subtitle = bquote(R^2 ~ "=" ~ .(0.5065741) ~ ";" ~ P ~ "=" ~.(0.003823239) ~ ";" ~ lambda ~ "=" ~ .(0))) +
  xlab("breeding_type") +
  ylab("prop_Telencephalon") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none") +
  scale_color_manual(values=tribes.colors) 


dev.off()

