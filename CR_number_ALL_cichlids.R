#### Libraries  ---------------------------------

set.seed(2712)

options(bitmapType = "cairo")


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
#library("MoreTreeTools")
library("ggstar")
library(MASS)
#library(PCAtest)
library("viridis")
library(ggnewscale)

#### Data load  ---------------------------------

rm(list = ls())

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
  c(Total_OR="#E69F00",
    Total_TAAR = "#2995D2", 
    Total_V1R = "#FF7474",
    Total_V2R = "#22A481",
    Total_T1R = "gray65", 
    Total_T2R = "#F0E442")




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


#Initale subfam colors

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

#### Plot the number of OR gene per species / per specimen -- DOC  ---------------------------------

#Orderer by the mean number of OR genes computed

OR_DoC_df_info_total <- 
  OR_DoC_df_info %>%
  filter(Subfamily == "Total_OR") 


OR_DoC_df_info_total %>%
  filter(Sp == "Neospl")

OR_DoC_df_info_total %>%
  filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>%
  group_by(Sp) %>%
  mutate(count = n()) %>%
  dplyr::select(Sp, count) %>%
  filter(count < 2)


ID_for_Rer <- 
  OR_DoC_df_info_total %>%
  filter(Sp != "Neospl") %>%
  filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>%
  dplyr::select(ID, Sp) %>%
  group_by(Sp) %>%
  slice(1) %>%
  ungroup()

pdf(file = "DoC_OR_per_ind_wo_legend.pdf",width = 6.34,  height = 4.61)

OR_DoC_df_info %>%
  filter(Subfamily == "Total_OR") %>%
  filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>%
  filter(Sp != "Neospl") %>%
  ggplot(aes(x= reorder(Sp, normalized_nb), y = normalized_nb, color=(Tribe))) +
  geom_line() +
  geom_point() +
  scale_color_manual(values=tribes.colors) +
  theme_classic() +
  ylab("# of OR genes - DoC") +
  xlab("Species") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.text.x=element_blank(),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position = "none")

dev.off()



pdf(file = "DoC_OR_per_ind_with_legend.pdf",width = 6.34,  height = 4.61)

OR_DoC_df_info %>%
  filter(Subfamily == "Total_OR") %>%
  filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>%
  filter(Sp != "Neospl") %>%
  ggplot(aes(x= reorder(Sp, normalized_nb), y = normalized_nb, color=(Tribe))) +
  geom_line() +
  geom_point() +
  scale_color_manual(values=tribes.colors) +
  theme_classic() +
  ylab("# of OR genes - DoC") +
  xlab("Species") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.text.x=element_blank(),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) 

dev.off()

#Ordered according to the phylogeny

OR_DoC_df_info %>%
  filter(Subfamily == "Total_OR") %>%
  filter(Sp == "Cteben") 


OR_DoC_df_info_fortree <- 
  OR_DoC_df_info %>%
  filter(Subfamily == "Total_OR") %>%
  dplyr::select(Sp, normalized_nb) %>%
  mutate(Species = Sp)


p <- 
  ggtree(radiation_tree_wo_Neospl, size=0.6) %<+% Tribe_for_tree +
  aes(color=Tribe) +
  geom_tiplab(size=0.8) +
  scale_color_manual(values = tribes.colors) +
  theme(legend.position = "none")

facet_plot(p, panel='OR', data=OR_DoC_df_info_fortree, geom=geom_point, 
           aes(x=normalized_nb,  group=Species), size=0.5) +
  theme_tree2()


#### Plot the number of TAAR gene per species / per specimen -- DOC  ---------------------------------

#Orderer by the mean number of TAAR genes computed

TAAR_DoC_df_info_total <- 
  TAAR_DoC_df_info %>%
  filter(Subfamily == "Total_TAAR") 

pdf(file = "DoC_TAAR_per_ind_wo_legend.pdf",width = 6.34,  height = 4.61)

TAAR_DoC_df_info %>%
  filter(Subfamily == "Total_TAAR") %>%
  filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>%
  filter(Sp != "Neospl") %>%
  ggplot(aes(x= reorder(Sp, normalized_nb), y = normalized_nb, color=(Tribe))) +
  geom_line() +
  geom_point() +
  scale_color_manual(values=tribes.colors) +
  theme_classic() +
  ylab("# of TAAR genes - DoC") +
  xlab("Species") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.text.x=element_blank(),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position = "none")

dev.off()

#Ordered according to the phylogeny

TAAR_DoC_df_info_fortree <- 
  TAAR_DoC_df_info %>%
  filter(Subfamily == "Total_TAAR") %>%
  dplyr::select(Sp, normalized_nb) %>%
  mutate(Species = Sp)


p <- 
  ggtree(b1_tree_wo_Neospl, size=0.6) %<+% Tribe_for_tree +
  aes(color=Tribe) +
  geom_tiplab(size=0.8) +
  scale_color_manual(values = tribes.colors) +
  theme(legend.position = "none")

facet_plot(p, panel='TAAR', data=TAAR_DoC_df_info_fortree, geom=geom_line, 
           aes(x=normalized_nb,  group=Species), size=0.5) +
  theme_tree2()

TAAR_DoC_df_info %>%
  filter(Subfamily == "Total_TAAR") %>%
  filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>% 
  pull(normalized_nb) %>% min()

TAAR_DoC_df_info %>%
  filter(Subfamily == "Total_TAAR") %>%
  filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>% 
  pull(normalized_nb) %>% max()

#### Plot the number of V2R gene per species / per specimen -- DOC  ---------------------------------

#Orderer by the mean number of V2R genes computed

V2R_DoC_df_info_total <- 
  V2R_DoC_df_info %>%
  filter(Subfamily == "Total_V2R") 


pdf(file = "DoC_V2R_per_ind_wo_legend.pdf",width = 6.34,  height = 4.61)

V2R_DoC_df_info %>%
  filter(Subfamily == "Total_V2R") %>%
  filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>%
  filter(Sp != "Neospl") %>%
  ggplot(aes(x= reorder(Sp, normalized_nb), y = normalized_nb, color=(Tribe))) +
  geom_line() +
  geom_point() +
  scale_color_manual(values=tribes.colors) +
  theme_classic() +
  ylab("# of V2R genes - DoC") +
  xlab("Species") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.text.x=element_blank(),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position = "none")

dev.off()


#Ordered according to the phylogeny

V2R_DoC_df_info_fortree <- 
  V2R_DoC_df_info %>%
  filter(Subfamily == "Total_V2R") %>%
  dplyr::select(Sp, normalized_nb) %>%
  mutate(Species = Sp)


p <- 
  ggtree(b1_tree_wo_Neospl, size=0.6) %<+% Tribe_for_tree +
  aes(color=Tribe) +
  geom_tiplab(size=0.8) +
  scale_color_manual(values = tribes.colors) +
  theme(legend.position = "none")

facet_plot(p, panel='V2R', data=V2R_DoC_df_info_fortree, geom=geom_line, 
           aes(x=normalized_nb,  group=Species), size=0.5) +
  theme_tree2()

#### Plot the number of V1R gene per species / per specimen -- DOC  ---------------------------------

#Orderer by the mean number of V1R genes computed

V1R_DoC_df_info_total <- 
  V1R_DoC_df_info %>%
  filter(Subfamily == "Total_V1R") 



pdf(file = "DoC_V1R_per_ind_wo_legend.pdf",width = 6.34,  height = 4.61)

V1R_DoC_df_info %>%
  filter(Subfamily == "Total_V1R") %>%
  filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>%
  filter(Sp != "Neospl") %>%
  ggplot(aes(x= reorder(Sp, normalized_nb), y = normalized_nb, color=(Tribe))) +
  geom_line() +
  geom_point() +
  scale_color_manual(values=tribes.colors) +
  theme_classic() +
  ylab("# of V1R genes - DoC") +
  xlab("Species") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.text.x=element_blank(),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position = "none") +
  ylim(4, 7)

dev.off()




#Ordered according to the phylogeny

V1R_DoC_df_info_fortree <- 
  V1R_DoC_df_info %>%
  filter(Subfamily == "Total_V1R") %>%
  dplyr::select(Sp, normalized_nb) %>%
  mutate(Species = Sp)


p <- 
  ggtree(b1_tree_wo_Neospl, size=0.6) %<+% Tribe_for_tree +
  aes(color=Tribe) +
  geom_tiplab(size=0.8) +
  scale_color_manual(values = tribes.colors) +
  theme(legend.position = "none")

facet_plot(p, panel='V1R', data=V1R_DoC_df_info_fortree, geom=geom_line, 
           aes(x=normalized_nb,  group=Species), size=0.5) +
  theme_tree2()

#### Plot the number of T1R gene per species / per specimen -- DOC  ---------------------------------

#Orderer by the mean number of T1R genes computed

T1R_DoC_df_info_total <- 
  T1R_DoC_df_info %>%
  filter(Subfamily == "Total_T1R") 

T1R_DoC_df_info_total %>%
  arrange(desc(normalized_nb))


pdf(file = "DoC_T1R_per_ind_wo_legend.pdf",width = 6.34,  height = 4.61)

T1R_DoC_df_info %>%
  filter(Subfamily == "Total_T1R") %>%
  filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>%
  filter(Sp != "Neospl") %>%
  ggplot(aes(x= reorder(Sp, normalized_nb), y = normalized_nb, color=(Tribe))) +
  geom_line() +
  geom_point() +
  scale_color_manual(values=tribes.colors) +
  theme_classic() +
  ylab("# of T1R genes - DoC") +
  xlab("Species") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.text.x=element_blank(),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position = "none")

dev.off()





#Ordered according to the phylogeny

T1R_DoC_df_info_fortree <- 
  T1R_DoC_df_info %>%
  filter(Subfamily == "Total_T1R") %>%
  dplyr::select(Sp, normalized_nb) %>%
  mutate(Species = Sp)


p <- 
  ggtree(b1_tree_wo_Neospl, size=0.6) %<+% Tribe_for_tree +
  aes(color=Tribe) +
  geom_tiplab(size=0.8) +
  scale_color_manual(values = tribes.colors) +
  theme(legend.position = "none")

facet_plot(p, panel='T1R', data=T1R_DoC_df_info_fortree, geom=geom_line, 
           aes(x=normalized_nb,  group=Species), size=0.5) +
  theme_tree2()


#### Plot the number of T2R gene per species / per specimen -- DOC  ---------------------------------

#Orderer by the mean number of T2R genes computed

T2R_DoC_df_info_total <- 
  T2R_DoC_df_info %>%
  filter(Subfamily == "Total_T2R") 


pdf(file = "DoC_T2R_per_ind_wo_legend.pdf",width = 6.34,  height = 4.61)

T2R_DoC_df_info %>%
  filter(Subfamily == "Total_T2R") %>%
  filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>%
  filter(Sp != "Neospl") %>%
  ggplot(aes(x= reorder(Sp, normalized_nb), y = normalized_nb, color=(Tribe))) +
  geom_line() +
  geom_point() +
  scale_color_manual(values=tribes.colors) +
  theme_classic() +
  ylab("# of T2R genes - DoC") +
  xlab("Species") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.text.x=element_blank(),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position = "none") +
  ylim(0, 1.5)

dev.off()



#Ordered according to the phylogeny

T2R_DoC_df_info_fortree <- 
  T2R_DoC_df_info %>%
  filter(Subfamily == "Total_T2R") %>%
  dplyr::select(Sp, normalized_nb) %>%
  mutate(Species = Sp)







p <- 
  ggtree(b1_tree_wo_Neospl, size=0.6) %<+% Tribe_for_tree +
  aes(color=Tribe) +
  geom_tiplab(size=0.8) +
  scale_color_manual(values = tribes.colors) +
  theme(legend.position = "none")

facet_plot(p, panel='T2R', data=T2R_DoC_df_info_fortree, geom=geom_line, 
           aes(x=normalized_nb,  group=Species), size=0.5) +
  theme_tree2()


#### Boxplot -  OR gene per Tribe -- DOC  ---------------------------------

OR_DoC_df_info_mean_total[(OR_DoC_df_info_mean_total$Tribe == "Serranochromini"),"Tribe"] <- "Haplochromini"


t_pacbio_olfactory_receptor %>%
  filter(Species == "Bathybates_minor") %>%
  mutate(ALL_OR = OR_Complete + OR_Pseudogene + OR_Truncated) %>%
  pull(ALL_OR)

pdf(file = "BoxPlot_Tribe_OR.pdf",width = 6.34,  height = 4.61)

OR_DoC_df_info_mean_total %>%
  filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>%
  ggplot(., aes(x=reorder(Tribe,mean_normalized_nb), y=mean_normalized_nb, fill=Tribe)) +
  geom_boxplot() + 
  scale_fill_manual(values = tribes.colors) +
  scale_color_manual(values = tribes.colors) +
  theme_classic() + 
  xlab("Tribe") +
  ylab("# of OR genes") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.text.x=element_blank(),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position = "none")


dev.off()


pdf(file = "BoxPlot_Tribe_OR.pdf",width = 6.34,  height = 4.61)

OR_DoC_df_info_mean_total %>%
  filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>%
  ggplot(., aes(x=reorder(Tribe,mean_normalized_nb), y=mean_normalized_nb, fill=Tribe)) +
  geom_boxplot() + 
  scale_fill_manual(values = tribes.colors) +
  scale_color_manual(values = tribes.colors) +
  theme_classic() + 
  xlab("Tribe") +
  ylab("# of OR genes") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.text.x=element_blank(),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position = "none")

dev.off()






#### Boxplot -  TAAR gene per Tribe -- DOC  ---------------------------------


TAAR_DoC_df_info_mean_total[(TAAR_DoC_df_info_mean_total$Tribe == "Serranochromini"),"Tribe"] <- "Haplochromini"


pdf(file = "BoxPlot_Tribe_TAAR.pdf",width = 6.34,  height = 4.61)

TAAR_DoC_df_info_mean_total %>%
  filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>%
  ggplot(., aes(x=reorder(Tribe,mean_normalized_nb), y=mean_normalized_nb, fill=Tribe)) +
  geom_boxplot() + 
  scale_fill_manual(values = tribes.colors) +
  scale_color_manual(values = tribes.colors) +
  theme_classic() + 
  xlab("Tribe") +
  ylab("# of TAAR genes") +
  geom_star(aes(x=
                  TAAR_DoC_df_info_mean_total %>% 
                  filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>%
                  group_by(Tribe) %>%
                  summarise(mean_TAAR = mean(mean_normalized_nb)) %>%
                  arrange(mean_TAAR) %>%
                  mutate(order = c(1,2,3,4,5,6,7,8,9,10,11,12,13)) %>%
                  filter(Tribe == "Bathybatini") %>%
                  pull(order), 
                y=
                  t_pacbio_olfactory_receptor %>%
                  filter(Species == "Bathybates_minor") %>%
                  mutate(ALL_TAAR = TAAR_Complete + TAAR_Pseudogene + TAAR_Truncated) %>%
                  pull(ALL_TAAR)),
            color="black",
            fill="red",
            size=4) +
  geom_star(aes(x=
                  TAAR_DoC_df_info_mean_total %>% 
                  filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>%
                  group_by(Tribe) %>%
                  summarise(mean_TAAR = mean(mean_normalized_nb)) %>%
                  arrange(mean_TAAR) %>%
                  mutate(order = c(1,2,3,4,5,6,7,8,9,10,11,12,13)) %>%
                  filter(Tribe == "Cyphotilapiini") %>%
                  pull(order), 
                y=
                  t_pacbio_olfactory_receptor %>%
                  filter(Species == "Cyphotilapia_frontosa") %>%
                  mutate(ALL_TAAR = TAAR_Complete + TAAR_Pseudogene + TAAR_Truncated) %>%
                  pull(ALL_TAAR)),
            color="black",
            fill="red",
            size=4) +
  geom_star(aes(x=
                  TAAR_DoC_df_info_mean_total %>% 
                  filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>%
                  group_by(Tribe) %>%
                  summarise(mean_TAAR = mean(mean_normalized_nb)) %>%
                  arrange(mean_TAAR) %>%
                  mutate(order = c(1,2,3,4,5,6,7,8,9,10,11,12,13)) %>%
                  filter(Tribe == "Ectodini") %>%
                  pull(order), 
                y=
                  t_pacbio_olfactory_receptor %>%
                  filter(Species == "Cunningtonia_longiventralis") %>%
                  mutate(ALL_TAAR = TAAR_Complete + TAAR_Pseudogene + TAAR_Truncated) %>%
                  pull(ALL_TAAR)),
            color="black",
            fill="red",
            size=4) +
  geom_star(aes(x=
                  TAAR_DoC_df_info_mean_total %>% 
                  filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>%
                  group_by(Tribe) %>%
                  summarise(mean_TAAR = mean(mean_normalized_nb)) %>%
                  arrange(mean_TAAR) %>%
                  mutate(order = c(1,2,3,4,5,6,7,8,9,10,11,12,13)) %>%
                  filter(Tribe == "Cyprichromini") %>%
                  pull(order), 
                y=
                  t_pacbio_olfactory_receptor %>%
                  filter(Species == "Cyprichromis_leptosoma") %>%
                  mutate(ALL_TAAR = TAAR_Complete + TAAR_Pseudogene + TAAR_Truncated) %>%
                  pull(ALL_TAAR)),
            color="black",
            fill="red",
            size=4) +
  geom_star(aes(x=
                  TAAR_DoC_df_info_mean_total %>% 
                  filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>%
                  group_by(Tribe) %>%
                  summarise(mean_TAAR = mean(mean_normalized_nb)) %>%
                  arrange(mean_TAAR) %>%
                  mutate(order = c(1,2,3,4,5,6,7,8,9,10,11,12,13)) %>%
                  filter(Tribe == "Tropheini") %>%
                  pull(order), 
                y=
                  t_pacbio_olfactory_receptor %>%
                  filter(Species == "Simochromis_diagramma") %>%
                  mutate(ALL_TAAR = TAAR_Complete + TAAR_Pseudogene + TAAR_Truncated) %>%
                  pull(ALL_TAAR)),
            color="black",
            fill="red",
            size=4) +
  geom_star(aes(x=
                  TAAR_DoC_df_info_mean_total %>% 
                  filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>%
                  group_by(Tribe) %>%
                  summarise(mean_TAAR = mean(mean_normalized_nb)) %>%
                  arrange(mean_TAAR) %>%
                  mutate(order = c(1,2,3,4,5,6,7,8,9,10,11,12,13)) %>%
                  filter(Tribe == "Lamprologini") %>%
                  pull(order), 
                y=
                  t_pacbio_olfactory_receptor %>%
                  filter(Species == "Neolamprologus_multifasciatus") %>%
                  mutate(ALL_TAAR = TAAR_Complete + TAAR_Pseudogene + TAAR_Truncated) %>%
                  pull(ALL_TAAR)),
            color="black",
            fill="red",
            size=4) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.text.x=element_blank(),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position = "none")

dev.off()


pdf(file = "BoxPlot_Tribe_TAAR.pdf",width = 6.34,  height = 4.61)

TAAR_DoC_df_info_mean_total %>%
  filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>%
  ggplot(., aes(x=reorder(Tribe,mean_normalized_nb), y=mean_normalized_nb, fill=Tribe)) +
  geom_boxplot() + 
  scale_fill_manual(values = tribes.colors) +
  scale_color_manual(values = tribes.colors) +
  theme_classic() + 
  xlab("Tribe") +
  ylab("# of TAAR genes") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.text.x=element_blank(),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position = "none")

dev.off()



#### Boxplot -  V1R gene per Tribe -- DOC  ---------------------------------


V1R_DoC_df_info_mean_total[(V1R_DoC_df_info_mean_total$Tribe == "Serranochromini"),"Tribe"] <- "Haplochromini"


pdf(file = "BoxPlot_Tribe_V1R.pdf",width = 6.34,  height = 4.61)

V1R_DoC_df_info_mean_total %>%
  filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>%
  ggplot(., aes(x=reorder(Tribe,mean_normalized_nb), y=mean_normalized_nb, fill=Tribe)) +
  geom_boxplot() + 
  scale_fill_manual(values = tribes.colors) +
  scale_color_manual(values = tribes.colors) +
  theme_classic() + 
  xlab("Tribe") +
  ylab("# of V1R genes") +
  geom_star(aes(x=
                  V1R_DoC_df_info_mean_total %>% 
                  filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>%
                  group_by(Tribe) %>%
                  summarise(mean_V1R = mean(mean_normalized_nb)) %>%
                  arrange(mean_V1R) %>%
                  mutate(order = c(1,2,3,4,5,6,7,8,9,10,11,12,13)) %>%
                  filter(Tribe == "Bathybatini") %>%
                  pull(order), 
                y=
                  t_pacbio_olfactory_receptor %>%
                  filter(Species == "Bathybates_minor") %>%
                  mutate(ALL_V1R = V1R_Complete + V1R_Pseudogene + V1R_Truncated) %>%
                  pull(ALL_V1R)),
            color="black",
            fill="red",
            size=4) +
  geom_star(aes(x=
                  V1R_DoC_df_info_mean_total %>% 
                  filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>%
                  group_by(Tribe) %>%
                  summarise(mean_V1R = mean(mean_normalized_nb)) %>%
                  arrange(mean_V1R) %>%
                  mutate(order = c(1,2,3,4,5,6,7,8,9,10,11,12,13)) %>%
                  filter(Tribe == "Cyphotilapiini") %>%
                  pull(order), 
                y=
                  t_pacbio_olfactory_receptor %>%
                  filter(Species == "Cyphotilapia_frontosa") %>%
                  mutate(ALL_V1R = V1R_Complete + V1R_Pseudogene + V1R_Truncated) %>%
                  pull(ALL_V1R)),
            color="black",
            fill="red",
            size=4) +
  geom_star(aes(x=
                  V1R_DoC_df_info_mean_total %>% 
                  filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>%
                  group_by(Tribe) %>%
                  summarise(mean_V1R = mean(mean_normalized_nb)) %>%
                  arrange(mean_V1R) %>%
                  mutate(order = c(1,2,3,4,5,6,7,8,9,10,11,12,13)) %>%
                  filter(Tribe == "Ectodini") %>%
                  pull(order), 
                y=
                  t_pacbio_olfactory_receptor %>%
                  filter(Species == "Cunningtonia_longiventralis") %>%
                  mutate(ALL_V1R = V1R_Complete + V1R_Pseudogene + V1R_Truncated) %>%
                  pull(ALL_V1R)),
            color="black",
            fill="red",
            size=4) +
  geom_star(aes(x=
                  V1R_DoC_df_info_mean_total %>% 
                  filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>%
                  group_by(Tribe) %>%
                  summarise(mean_V1R = mean(mean_normalized_nb)) %>%
                  arrange(mean_V1R) %>%
                  mutate(order = c(1,2,3,4,5,6,7,8,9,10,11,12,13)) %>%
                  filter(Tribe == "Cyprichromini") %>%
                  pull(order), 
                y=
                  t_pacbio_olfactory_receptor %>%
                  filter(Species == "Cyprichromis_leptosoma") %>%
                  mutate(ALL_V1R = V1R_Complete + V1R_Pseudogene + V1R_Truncated) %>%
                  pull(ALL_V1R)),
            color="black",
            fill="red",
            size=4) +
  geom_star(aes(x=
                  V1R_DoC_df_info_mean_total %>% 
                  filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>%
                  group_by(Tribe) %>%
                  summarise(mean_V1R = mean(mean_normalized_nb)) %>%
                  arrange(mean_V1R) %>%
                  mutate(order = c(1,2,3,4,5,6,7,8,9,10,11,12,13)) %>%
                  filter(Tribe == "Tropheini") %>%
                  pull(order), 
                y=
                  t_pacbio_olfactory_receptor %>%
                  filter(Species == "Simochromis_diagramma") %>%
                  mutate(ALL_V1R = V1R_Complete + V1R_Pseudogene + V1R_Truncated) %>%
                  pull(ALL_V1R)),
            color="black",
            fill="red",
            size=4) +
  geom_star(aes(x=
                  V1R_DoC_df_info_mean_total %>% 
                  filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>%
                  group_by(Tribe) %>%
                  summarise(mean_V1R = mean(mean_normalized_nb)) %>%
                  arrange(mean_V1R) %>%
                  mutate(order = c(1,2,3,4,5,6,7,8,9,10,11,12,13)) %>%
                  filter(Tribe == "Lamprologini") %>%
                  pull(order), 
                y=
                  t_pacbio_olfactory_receptor %>%
                  filter(Species == "Neolamprologus_multifasciatus") %>%
                  mutate(ALL_V1R = V1R_Complete + V1R_Pseudogene + V1R_Truncated) %>%
                  pull(ALL_V1R)),
            color="black",
            fill="red",
            size=4) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.text.x=element_blank(),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position = "none")

dev.off()


pdf(file = "BoxPlot_Tribe_V1R.pdf",width = 6.34,  height = 4.61)

V1R_DoC_df_info_mean_total %>%
  filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>%
  ggplot(., aes(x=reorder(Tribe,mean_normalized_nb), y=mean_normalized_nb, fill=Tribe)) +
  geom_boxplot() + 
  scale_fill_manual(values = tribes.colors) +
  scale_color_manual(values = tribes.colors) +
  theme_classic() + 
  xlab("Tribe") +
  ylab("# of V1R genes") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.text.x=element_blank(),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position = "none") +
  ylim(c(5, 6.5))

dev.off()



#### Boxplot -  V2R gene per Tribe -- DOC  ---------------------------------


V2R_DoC_df_info_mean_total[(V2R_DoC_df_info_mean_total$Tribe == "Serranochromini"),"Tribe"] <- "Haplochromini"


pdf(file = "BoxPlot_Tribe_V2R.pdf",width = 6.34,  height = 4.61)

V2R_DoC_df_info_mean_total %>%
  filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>%
  ggplot(., aes(x=reorder(Tribe,mean_normalized_nb), y=mean_normalized_nb, fill=Tribe)) +
  geom_boxplot() + 
  scale_fill_manual(values = tribes.colors) +
  scale_color_manual(values = tribes.colors) +
  theme_classic() + 
  xlab("Tribe") +
  ylab("# of V2R genes") +
  geom_star(aes(x=
                  V2R_DoC_df_info_mean_total %>% 
                  filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>%
                  group_by(Tribe) %>%
                  summarise(mean_V2R = mean(mean_normalized_nb)) %>%
                  arrange(mean_V2R) %>%
                  mutate(order = c(1,2,3,4,5,6,7,8,9,10,11,12,13)) %>%
                  filter(Tribe == "Bathybatini") %>%
                  pull(order), 
                y=
                  t_pacbio_olfactory_receptor %>%
                  filter(Species == "Bathybates_minor") %>%
                  mutate(ALL_V2R = V2R_Complete + V2R_Pseudogene + V2R_Truncated) %>%
                  pull(ALL_V2R)),
            color="black",
            fill="red",
            size=4) +
  geom_star(aes(x=
                  V2R_DoC_df_info_mean_total %>% 
                  filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>%
                  group_by(Tribe) %>%
                  summarise(mean_V2R = mean(mean_normalized_nb)) %>%
                  arrange(mean_V2R) %>%
                  mutate(order = c(1,2,3,4,5,6,7,8,9,10,11,12,13)) %>%
                  filter(Tribe == "Cyphotilapiini") %>%
                  pull(order), 
                y=
                  t_pacbio_olfactory_receptor %>%
                  filter(Species == "Cyphotilapia_frontosa") %>%
                  mutate(ALL_V2R = V2R_Complete + V2R_Pseudogene + V2R_Truncated) %>%
                  pull(ALL_V2R)),
            color="black",
            fill="red",
            size=4) +
  geom_star(aes(x=
                  V2R_DoC_df_info_mean_total %>% 
                  filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>%
                  group_by(Tribe) %>%
                  summarise(mean_V2R = mean(mean_normalized_nb)) %>%
                  arrange(mean_V2R) %>%
                  mutate(order = c(1,2,3,4,5,6,7,8,9,10,11,12,13)) %>%
                  filter(Tribe == "Ectodini") %>%
                  pull(order), 
                y=
                  t_pacbio_olfactory_receptor %>%
                  filter(Species == "Cunningtonia_longiventralis") %>%
                  mutate(ALL_V2R = V2R_Complete + V2R_Pseudogene + V2R_Truncated) %>%
                  pull(ALL_V2R)),
            color="black",
            fill="red",
            size=4) +
  geom_star(aes(x=
                  V2R_DoC_df_info_mean_total %>% 
                  filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>%
                  group_by(Tribe) %>%
                  summarise(mean_V2R = mean(mean_normalized_nb)) %>%
                  arrange(mean_V2R) %>%
                  mutate(order = c(1,2,3,4,5,6,7,8,9,10,11,12,13)) %>%
                  filter(Tribe == "Cyprichromini") %>%
                  pull(order), 
                y=
                  t_pacbio_olfactory_receptor %>%
                  filter(Species == "Cyprichromis_leptosoma") %>%
                  mutate(ALL_V2R = V2R_Complete + V2R_Pseudogene + V2R_Truncated) %>%
                  pull(ALL_V2R)),
            color="black",
            fill="red",
            size=4) +
  geom_star(aes(x=
                  V2R_DoC_df_info_mean_total %>% 
                  filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>%
                  group_by(Tribe) %>%
                  summarise(mean_V2R = mean(mean_normalized_nb)) %>%
                  arrange(mean_V2R) %>%
                  mutate(order = c(1,2,3,4,5,6,7,8,9,10,11,12,13)) %>%
                  filter(Tribe == "Tropheini") %>%
                  pull(order), 
                y=
                  t_pacbio_olfactory_receptor %>%
                  filter(Species == "Simochromis_diagramma") %>%
                  mutate(ALL_V2R = V2R_Complete + V2R_Pseudogene + V2R_Truncated) %>%
                  pull(ALL_V2R)),
            color="black",
            fill="red",
            size=4) +
  geom_star(aes(x=
                  V2R_DoC_df_info_mean_total %>% 
                  filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>%
                  group_by(Tribe) %>%
                  summarise(mean_V2R = mean(mean_normalized_nb)) %>%
                  arrange(mean_V2R) %>%
                  mutate(order = c(1,2,3,4,5,6,7,8,9,10,11,12,13)) %>%
                  filter(Tribe == "Lamprologini") %>%
                  pull(order), 
                y=
                  t_pacbio_olfactory_receptor %>%
                  filter(Species == "Neolamprologus_multifasciatus") %>%
                  mutate(ALL_V2R = V2R_Complete + V2R_Pseudogene + V2R_Truncated) %>%
                  pull(ALL_V2R)),
            color="black",
            fill="red",
            size=4) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.text.x=element_blank(),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position = "none")

dev.off()


pdf(file = "BoxPlot_Tribe_V2R.pdf",width = 6.34,  height = 4.61)

V2R_DoC_df_info_mean_total %>%
  filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>%
  ggplot(., aes(x=reorder(Tribe,mean_normalized_nb), y=mean_normalized_nb, fill=Tribe)) +
  geom_boxplot() + 
  scale_fill_manual(values = tribes.colors) +
  scale_color_manual(values = tribes.colors) +
  theme_classic() + 
  xlab("Tribe") +
  ylab("# of V2R genes") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.text.x=element_blank(),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position = "none")

dev.off()


#### Boxplot -  T1R gene per Tribe -- DOC  ---------------------------------


T1R_DoC_df_info_mean_total[(T1R_DoC_df_info_mean_total$Tribe == "Serranochromini"),"Tribe"] <- "Haplochromini"


pdf(file = "BoxPlot_Tribe_T1R.pdf",width = 6.34,  height = 4.61)

T1R_DoC_df_info_mean_total %>%
  filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>%
  ggplot(., aes(x=reorder(Tribe,mean_normalized_nb), y=mean_normalized_nb, fill=Tribe)) +
  geom_boxplot() + 
  scale_fill_manual(values = tribes.colors) +
  scale_color_manual(values = tribes.colors) +
  theme_classic() + 
  xlab("Tribe") +
  ylab("# of T1R genes") +
  geom_star(aes(x=
                  T1R_DoC_df_info_mean_total %>% 
                  filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>%
                  group_by(Tribe) %>%
                  summarise(mean_T1R = mean(mean_normalized_nb)) %>%
                  arrange(mean_T1R) %>%
                  mutate(order = c(1,2,3,4,5,6,7,8,9,10,11,12,13)) %>%
                  filter(Tribe == "Bathybatini") %>%
                  pull(order), 
                y=
                  t_pacbio_olfactory_receptor %>%
                  filter(Species == "Bathybates_minor") %>%
                  mutate(ALL_T1R = T1R_Complete + T1R_Pseudogene + T1R_Truncated) %>%
                  pull(ALL_T1R)),
            color="black",
            fill="red",
            size=4) +
  geom_star(aes(x=
                  T1R_DoC_df_info_mean_total %>% 
                  filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>%
                  group_by(Tribe) %>%
                  summarise(mean_T1R = mean(mean_normalized_nb)) %>%
                  arrange(mean_T1R) %>%
                  mutate(order = c(1,2,3,4,5,6,7,8,9,10,11,12,13)) %>%
                  filter(Tribe == "Cyphotilapiini") %>%
                  pull(order), 
                y=
                  t_pacbio_olfactory_receptor %>%
                  filter(Species == "Cyphotilapia_frontosa") %>%
                  mutate(ALL_T1R = T1R_Complete + T1R_Pseudogene + T1R_Truncated) %>%
                  pull(ALL_T1R)),
            color="black",
            fill="red",
            size=4) +
  geom_star(aes(x=
                  T1R_DoC_df_info_mean_total %>% 
                  filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>%
                  group_by(Tribe) %>%
                  summarise(mean_T1R = mean(mean_normalized_nb)) %>%
                  arrange(mean_T1R) %>%
                  mutate(order = c(1,2,3,4,5,6,7,8,9,10,11,12,13)) %>%
                  filter(Tribe == "Ectodini") %>%
                  pull(order), 
                y=
                  t_pacbio_olfactory_receptor %>%
                  filter(Species == "Cunningtonia_longiventralis") %>%
                  mutate(ALL_T1R = T1R_Complete + T1R_Pseudogene + T1R_Truncated) %>%
                  pull(ALL_T1R)),
            color="black",
            fill="red",
            size=4) +
  geom_star(aes(x=
                  T1R_DoC_df_info_mean_total %>% 
                  filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>%
                  group_by(Tribe) %>%
                  summarise(mean_T1R = mean(mean_normalized_nb)) %>%
                  arrange(mean_T1R) %>%
                  mutate(order = c(1,2,3,4,5,6,7,8,9,10,11,12,13)) %>%
                  filter(Tribe == "Cyprichromini") %>%
                  pull(order), 
                y=
                  t_pacbio_olfactory_receptor %>%
                  filter(Species == "Cyprichromis_leptosoma") %>%
                  mutate(ALL_T1R = T1R_Complete + T1R_Pseudogene + T1R_Truncated) %>%
                  pull(ALL_T1R)),
            color="black",
            fill="red",
            size=4) +
  geom_star(aes(x=
                  T1R_DoC_df_info_mean_total %>% 
                  filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>%
                  group_by(Tribe) %>%
                  summarise(mean_T1R = mean(mean_normalized_nb)) %>%
                  arrange(mean_T1R) %>%
                  mutate(order = c(1,2,3,4,5,6,7,8,9,10,11,12,13)) %>%
                  filter(Tribe == "Tropheini") %>%
                  pull(order), 
                y=
                  t_pacbio_olfactory_receptor %>%
                  filter(Species == "Simochromis_diagramma") %>%
                  mutate(ALL_T1R = T1R_Complete + T1R_Pseudogene + T1R_Truncated) %>%
                  pull(ALL_T1R)),
            color="black",
            fill="red",
            size=4) +
  geom_star(aes(x=
                  T1R_DoC_df_info_mean_total %>% 
                  filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>%
                  group_by(Tribe) %>%
                  summarise(mean_T1R = mean(mean_normalized_nb)) %>%
                  arrange(mean_T1R) %>%
                  mutate(order = c(1,2,3,4,5,6,7,8,9,10,11,12,13)) %>%
                  filter(Tribe == "Lamprologini") %>%
                  pull(order), 
                y=
                  t_pacbio_olfactory_receptor %>%
                  filter(Species == "Neolamprologus_multifasciatus") %>%
                  mutate(ALL_T1R = T1R_Complete + T1R_Pseudogene + T1R_Truncated) %>%
                  pull(ALL_T1R)),
            color="black",
            fill="red",
            size=4) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.text.x=element_blank(),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position = "none")

dev.off()


pdf(file = "BoxPlot_Tribe_T1R.pdf",width = 6.34,  height = 4.61)

T1R_DoC_df_info_mean_total %>%
  filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>%
  ggplot(., aes(x=reorder(Tribe,mean_normalized_nb), y=mean_normalized_nb, fill=Tribe)) +
  geom_boxplot() + 
  scale_fill_manual(values = tribes.colors) +
  scale_color_manual(values = tribes.colors) +
  theme_classic() + 
  xlab("Tribe") +
  ylab("# of T1R genes") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.text.x=element_blank(),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position = "none")

dev.off()

#### Boxplot -  T2R gene per Tribe -- DOC  ---------------------------------


T2R_DoC_df_info_mean_total[(T2R_DoC_df_info_mean_total$Tribe == "Serranochromini"),"Tribe"] <- "Haplochromini"


pdf(file = "BoxPlot_Tribe_T2R.pdf",width = 6.34,  height = 4.61)

T2R_DoC_df_info_mean_total %>%
  filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>%
  ggplot(., aes(x=reorder(Tribe,mean_normalized_nb), y=mean_normalized_nb, fill=Tribe)) +
  geom_boxplot() + 
  scale_fill_manual(values = tribes.colors) +
  scale_color_manual(values = tribes.colors) +
  theme_classic() + 
  xlab("Tribe") +
  ylab("# of T2R genes") +
  geom_star(aes(x=
                  T2R_DoC_df_info_mean_total %>% 
                  filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>%
                  group_by(Tribe) %>%
                  summarise(mean_T2R = mean(mean_normalized_nb)) %>%
                  arrange(mean_T2R) %>%
                  mutate(order = c(1,2,3,4,5,6,7,8,9,10,11,12,13)) %>%
                  filter(Tribe == "Bathybatini") %>%
                  pull(order), 
                y=
                  t_pacbio_olfactory_receptor %>%
                  filter(Species == "Bathybates_minor") %>%
                  mutate(ALL_T2R = T2R_Complete + T2R_Pseudogene + T2R_Truncated) %>%
                  pull(ALL_T2R)),
            color="black",
            fill="red",
            size=4) +
  geom_star(aes(x=
                  T2R_DoC_df_info_mean_total %>% 
                  filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>%
                  group_by(Tribe) %>%
                  summarise(mean_T2R = mean(mean_normalized_nb)) %>%
                  arrange(mean_T2R) %>%
                  mutate(order = c(1,2,3,4,5,6,7,8,9,10,11,12,13)) %>%
                  filter(Tribe == "Cyphotilapiini") %>%
                  pull(order), 
                y=
                  t_pacbio_olfactory_receptor %>%
                  filter(Species == "Cyphotilapia_frontosa") %>%
                  mutate(ALL_T2R = T2R_Complete + T2R_Pseudogene + T2R_Truncated) %>%
                  pull(ALL_T2R)),
            color="black",
            fill="red",
            size=4) +
  geom_star(aes(x=
                  T2R_DoC_df_info_mean_total %>% 
                  filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>%
                  group_by(Tribe) %>%
                  summarise(mean_T2R = mean(mean_normalized_nb)) %>%
                  arrange(mean_T2R) %>%
                  mutate(order = c(1,2,3,4,5,6,7,8,9,10,11,12,13)) %>%
                  filter(Tribe == "Ectodini") %>%
                  pull(order), 
                y=
                  t_pacbio_olfactory_receptor %>%
                  filter(Species == "Cunningtonia_longiventralis") %>%
                  mutate(ALL_T2R = T2R_Complete + T2R_Pseudogene + T2R_Truncated) %>%
                  pull(ALL_T2R)),
            color="black",
            fill="red",
            size=4) +
  geom_star(aes(x=
                  T2R_DoC_df_info_mean_total %>% 
                  filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>%
                  group_by(Tribe) %>%
                  summarise(mean_T2R = mean(mean_normalized_nb)) %>%
                  arrange(mean_T2R) %>%
                  mutate(order = c(1,2,3,4,5,6,7,8,9,10,11,12,13)) %>%
                  filter(Tribe == "Cyprichromini") %>%
                  pull(order), 
                y=
                  t_pacbio_olfactory_receptor %>%
                  filter(Species == "Cyprichromis_leptosoma") %>%
                  mutate(ALL_T2R = T2R_Complete + T2R_Pseudogene + T2R_Truncated) %>%
                  pull(ALL_T2R)),
            color="black",
            fill="red",
            size=4) +
  geom_star(aes(x=
                  T2R_DoC_df_info_mean_total %>% 
                  filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>%
                  group_by(Tribe) %>%
                  summarise(mean_T2R = mean(mean_normalized_nb)) %>%
                  arrange(mean_T2R) %>%
                  mutate(order = c(1,2,3,4,5,6,7,8,9,10,11,12,13)) %>%
                  filter(Tribe == "Tropheini") %>%
                  pull(order), 
                y=
                  t_pacbio_olfactory_receptor %>%
                  filter(Species == "Simochromis_diagramma") %>%
                  mutate(ALL_T2R = T2R_Complete + T2R_Pseudogene + T2R_Truncated) %>%
                  pull(ALL_T2R)),
            color="black",
            fill="red",
            size=4) +
  geom_star(aes(x=
                  T2R_DoC_df_info_mean_total %>% 
                  filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>%
                  group_by(Tribe) %>%
                  summarise(mean_T2R = mean(mean_normalized_nb)) %>%
                  arrange(mean_T2R) %>%
                  mutate(order = c(1,2,3,4,5,6,7,8,9,10,11,12,13)) %>%
                  filter(Tribe == "Lamprologini") %>%
                  pull(order), 
                y=
                  t_pacbio_olfactory_receptor %>%
                  filter(Species == "Neolamprologus_multifasciatus") %>%
                  mutate(ALL_T2R = T2R_Complete + T2R_Pseudogene + T2R_Truncated) %>%
                  pull(ALL_T2R)),
            color="black",
            fill="red",
            size=4) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.text.x=element_blank(),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position = "none")

dev.off()


pdf(file = "BoxPlot_Tribe_T2R.pdf",width = 6.34,  height = 4.61)

T2R_DoC_df_info_mean_total %>%
  filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>%
  ggplot(., aes(x=reorder(Tribe,mean_normalized_nb), y=mean_normalized_nb, fill=Tribe)) +
  geom_boxplot() + 
  scale_fill_manual(values = tribes.colors) +
  scale_color_manual(values = tribes.colors) +
  theme_classic() + 
  xlab("Tribe") +
  ylab("# of T2R genes") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.text.x=element_blank(),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position = "none") +
  ylim(c(0, 1.5))

dev.off()
#### Boxplot -  OLR gene per Tribe -- DOC  ---------------------------------


All_Doc_mean_sp <- 
  rbind(OR_DoC_df_info_mean_total, TAAR_DoC_df_info_mean_total, 
        V1R_DoC_df_info_mean_total, V2R_DoC_df_info_mean_total)

All_Doc_mean_sp_TOTAL <- 
  as.data.frame(All_Doc_mean_sp %>%
  group_by(Sp, Tribe) %>%
  summarise(total_OLR = sum(mean_normalized_nb)))

All_Doc_mean_sp_TOTAL %>%
  filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>%
  arrange(total_OLR)

All_Doc_mean_sp_TOTAL %>%
  filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>%
  ggplot(., aes(x=reorder(Tribe,total_OLR), y=total_OLR, fill=Tribe)) +
  geom_boxplot() + 
  scale_fill_manual(values = tribes.colors) +
  scale_color_manual(values = tribes.colors) +
  theme_classic() + 
  xlab("Tribe") +
  ylab("# of OLR genes") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.text.x=element_blank(),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position = "none")


as.data.frame(
  All_Doc_mean_sp_TOTAL %>%
  filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>%
  group_by(Tribe) %>%
  summarise(mean_OLR = mean(total_OLR)) %>%
    arrange(mean_OLR)
)



#### Compute the mean diff between individuals of the same species + interspecies  ---------------------------------

#Compute the mean diff intra species + inter species

list_sp <- 
  b1_tree_wo_Neospl$tip.label

OR_DoC_df_info %>% 
  filter(Subfamily == "Total_OR") %>%
  filter(Sp == "Neomul")


vector_intrasp_OR <- c()
vector_intrasp_TAAR <- c()
vector_intrasp_V2R <- c()
vector_intrasp_T1R <- c()

vector_intersp_OR <- c()
vector_intersp_TAAR <- c()
vector_intersp_V2R <- c()
vector_intersp_T1R <- c()



mean(vector_intrasp_OR)
mean(vector_intrasp_TAAR)
mean(vector_intrasp_V2R)
mean(vector_intrasp_T1R)



mean(vector_intersp_OR)
mean(vector_intersp_TAAR)
mean(vector_intersp_V2R)
mean(vector_intersp_T1R)


OR_df_intra <- 
  as.data.frame(vector_intrasp_OR) %>%
  mutate(comp = "intra-species")
OR_df_inter <- 
  as.data.frame(vector_intersp_OR) %>%
  mutate(comp = "inter-species")
colnames(OR_df_intra) <- c("diff_nb", "comp")
colnames(OR_df_inter) <- c("diff_nb", "comp")

rbind(OR_df_intra, OR_df_inter) %>%
  ggplot(., aes(x=comp, y=diff_nb)) +
  geom_boxplot() +
  theme_classic()


#### Compute the phylogenetic signal of each gene family  ---------------------------------

Radiation_OR <- 
  OR_DoC_df_info_mean_total %>%
  filter(Sp %in% radiation_tree_wo_Neospl$tip.label)
tip_order <- radiation_tree_wo_Neospl$tip.label
Radiation_OR_ordered <- 
  enframe(tip_order, name = NULL, value = "Sp") %>% 
  left_join(Radiation_OR) %>% 
  column_to_rownames("Sp")

phylosig(tree = radiation_tree_wo_Neospl, x = Radiation_OR_ordered$mean_normalized_nb, method = "lambda", test = T)
phylosig(tree = radiation_tree_wo_Neospl, x = Radiation_OR_ordered$mean_normalized_nb, method = "K", test = T)

observed_phylosig_OR <- 
  phylosig(tree = radiation_tree_wo_Neospl, x = Radiation_OR_ordered$mean_normalized_nb, method = "lambda", test = T)$lambda

#### 

Radiation_TAAR <- 
  TAAR_DoC_df_info_mean_total %>%
  filter(Sp %in% radiation_tree_wo_Neospl$tip.label)
tip_order <- radiation_tree_wo_Neospl$tip.label
Radiation_TAAR_ordered <- 
  enframe(tip_order, name = NULL, value = "Sp") %>% 
  left_join(Radiation_TAAR) %>% 
  column_to_rownames("Sp")

phylosig(tree = radiation_tree_wo_Neospl, x = Radiation_TAAR_ordered$mean_normalized_nb, method = "lambda", test = T)


####


Radiation_V2R <- 
  V2R_DoC_df_info_mean_total %>%
  filter(Sp %in% radiation_tree_wo_Neospl$tip.label)
tip_order <- radiation_tree_wo_Neospl$tip.label
Radiation_V2R_ordered <- 
  enframe(tip_order, name = NULL, value = "Sp") %>% 
  left_join(Radiation_V2R) %>% 
  column_to_rownames("Sp")

phylosig(tree = radiation_tree_wo_Neospl, x = Radiation_V2R_ordered$mean_normalized_nb, method = "lambda", test = T)


####


Radiation_T1R <- 
  T1R_DoC_df_info_mean_total %>%
  filter(Sp %in% radiation_tree_wo_Neospl$tip.label)
tip_order <- radiation_tree_wo_Neospl$tip.label
Radiation_T1R_ordered <- 
  enframe(tip_order, name = NULL, value = "Sp") %>% 
  left_join(Radiation_T1R) %>% 
  column_to_rownames("Sp")

phylosig(tree = radiation_tree_wo_Neospl, x = Radiation_T1R_ordered$mean_normalized_nb, method = "lambda", test = T)



#### T2R


Radiation_T2R <- 
  T2R_DoC_df_info_mean_total %>%
  filter(Sp %in% radiation_tree_wo_Neospl$tip.label)
tip_order <- radiation_tree_wo_Neospl$tip.label
Radiation_T2R_ordered <- 
  enframe(tip_order, name = NULL, value = "Sp") %>% 
  left_join(Radiation_T2R) %>% 
  column_to_rownames("Sp")

phylosig(tree = radiation_tree_wo_Neospl, x = Radiation_T2R_ordered$mean_normalized_nb, method = "lambda", test = T)



#### V1R


Radiation_V1R <- 
  V1R_DoC_df_info_mean_total %>%
  filter(Sp %in% radiation_tree_wo_Neospl$tip.label)
tip_order <- radiation_tree_wo_Neospl$tip.label
Radiation_V1R_ordered <- 
  enframe(tip_order, name = NULL, value = "Sp") %>% 
  left_join(Radiation_V1R) %>% 
  column_to_rownames("Sp")

phylosig(tree = radiation_tree_wo_Neospl, x = Radiation_V1R_ordered$mean_normalized_nb, method = "lambda", test = T)


#### Import teleost chemoreceptor data ---------------------------------

species_taxa <- read.table("species_list.tsv", header=TRUE, sep="\t")
busco_results <- read.table("BUSCO_results_table.tsv", header=TRUE, sep="\t")
species_taxa <- left_join(species_taxa, busco_results, by="genome_name")
species_taxa <-
  species_taxa %>%
  mutate(busco_perc_complete = Complete/(Complete + Fragmented + Missing))
sp_tips_class <-
  species_taxa %>%
  dplyr::select(species, class)
sp_tips_order <-
  species_taxa %>%
  dplyr::select(species, order_ncbi)
chemoreceptor_tidy_df <- read.table("Table_count_chemoreceptors.csv", sep=",", header=FALSE)
colnames(chemoreceptor_tidy_df) <- c("species", "gene_clade", "gene_type", "gene_family","number")


chemoreceptor_tidy_df_total_all <- 
  chemoreceptor_tidy_df %>% 
  filter(gene_clade == "Total") %>% 
  filter(gene_type %in% c("Complete", "Pseudogene", "Truncated"))


chemoreceptor_tidy_df_total_all <- 
  as.data.frame(
    chemoreceptor_tidy_df_total_all %>%
      group_by(species, gene_clade, gene_family) %>%
      summarise(number_T = sum(number)))

chemoreceptor_tidy_df_total_all$gene_family <- 
  factor(chemoreceptor_tidy_df_total_all$gene_family, levels=c("T2R",'T1R','V2R',"V1R",'TAAR','OR'))



teleost_species <- 
  species_taxa %>% 
  filter(class == "Actinopterygii") %>% 
  filter(busco_perc_complete >= 0.8) %>%
  pull(species) %>% 
  unique()

teleost_chemoreceptors <- 
  chemoreceptor_tidy_df_total_all %>%
  filter(species %in% teleost_species)


# Remove cichlids 

teleost_to_remove <- 
  species_taxa %>% 
  filter(class == "Actinopterygii") %>%
  filter(order_ncbi == "Cichliformes") %>%
  pull(species) %>% 
  unique()


teleost_chemoreceptors <- 
  teleost_chemoreceptors %>%
  filter(! species %in% teleost_to_remove)


# Remove non teleost species

non_teleost_clades <- 
  c("Acipenseriformes", 
    "Polypteriformes", 
    "Amiiformes", 
    "Lepisosteiformes", 
    "Semionotiformes")


non_teleost_species <- 
  species_taxa %>% 
  filter(class == "Actinopterygii") %>%
  filter(order_ncbi %in% non_teleost_clades) %>%
  pull(species) %>% 
  unique()

teleost_chemoreceptors <- 
  teleost_chemoreceptors %>%
  filter(! species %in% non_teleost_species)


teleost_chemoreceptors$species %>% unique() #451 species 

order_ncbi_sp <- 
  species_taxa %>%
  dplyr::select(species, order_ncbi)


cichlid_teleost_color <- 
  c(Cichlid="#004D40",
    Teleost = "#FFC107")


#### Histogram teleost vs cichlid -- OR  ---------------------------------

OR_teleost <- 
  teleost_chemoreceptors %>%
  filter(gene_family == "OR") %>%
  dplyr::select(species, number_T) %>%
  mutate(hist_color = "Teleost")
OR_teleost <- 
  left_join(OR_teleost, order_ncbi_sp, by="species") 
colnames(OR_teleost) <- c("Species", "Number_OR", "Color_hist", "order")
OR_cichlids <- 
  OR_DoC_df_info_mean_total %>%
  dplyr::select(Sp, mean_normalized_nb) %>%
  mutate(hist_color = "Cichlid") %>%
  mutate(order = "Cichliformes")

colnames(OR_cichlids) <- c("Species", "Number_OR", "Color_hist", "order")

OR_cichlid_teleost <- 
  rbind(OR_teleost, OR_cichlids)


pdf(file = "Histo_TeleostCichlid_OR.pdf",width = 6.34,  height = 4.61)

OR_cichlid_teleost %>%
  ggplot(., aes(x=Number_OR, fill=Color_hist)) +
  geom_histogram(bins=40, alpha = 0.5, position = "identity", color="black") +
  theme_classic() +
  scale_fill_manual(values = cichlid_teleost_color) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position = "none") +
  xlab("Number of OR genes") +
  ylab("Number of species")

dev.off()


# lets compute mean per teleost, per cichlids, or per family ...

OR_cichlid_teleost %>%
  group_by(order) %>%
  summarise(number_species = n(),
            mean_OR = mean(Number_OR)) %>%
  arrange(desc(mean_OR))



OR_cichlid_teleost %>%
  group_by(Color_hist) %>%
  summarise(number_species = n(),
            mean_OR = mean(Number_OR)) %>%
  arrange(desc(mean_OR))



#### Histogram teleost vs cichlid -- TAAR  ---------------------------------

TAAR_teleost <- 
  teleost_chemoreceptors %>%
  filter(gene_family == "TAAR") %>%
  dplyr::select(species, number_T) %>%
  mutate(hist_color = "Teleost")
TAAR_teleost <- 
  left_join(TAAR_teleost, order_ncbi_sp, by="species") 
colnames(TAAR_teleost) <- c("Species", "Number_TAAR", "Color_hist", "order")
TAAR_cichlids <- 
  TAAR_DoC_df_info_mean_total %>%
  dplyr::select(Sp, mean_normalized_nb) %>%
  mutate(hist_color = "Cichlid") %>%
  mutate(order = "Cichliformes")

colnames(TAAR_cichlids) <- c("Species", "Number_TAAR", "Color_hist", "order")

TAAR_cichlid_teleost <- 
  rbind(TAAR_teleost, TAAR_cichlids)


pdf(file = "Histo_TeleostCichlid_TAAR.pdf",width = 6.34,  height = 4.61)

TAAR_cichlid_teleost %>%
  ggplot(., aes(x=Number_TAAR, fill=Color_hist)) +
  geom_histogram(bins=40, alpha = 0.5, position = "identity", color="black") +
  theme_classic() +
  scale_fill_manual(values = cichlid_teleost_color) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position = "none") +
  xlab("Number of TAAR genes") +
  ylab("Number of species")

dev.off()

# lets compute mean per teleost, per cichlids, or per family ...

TAAR_cichlid_teleost %>%
  group_by(order) %>%
  summarise(number_species = n(),
            mean_TAAR = mean(Number_TAAR)) %>%
  arrange(desc(mean_TAAR))



TAAR_cichlid_teleost %>%
  group_by(Color_hist) %>%
  summarise(number_species = n(),
            mean_TAAR = mean(Number_TAAR),
            min_TAAR = min(Number_TAAR),
            max_TAAR = max(Number_TAAR)) %>%
  arrange(desc(mean_TAAR))

#### Histogram teleost vs cichlid -- V2R  ---------------------------------

V2R_teleost <- 
  teleost_chemoreceptors %>%
  filter(gene_family == "V2R") %>%
  dplyr::select(species, number_T) %>%
  mutate(hist_color = "Teleost")
V2R_teleost <- 
  left_join(V2R_teleost, order_ncbi_sp, by="species") 
colnames(V2R_teleost) <- c("Species", "Number_V2R", "Color_hist", "order")
V2R_cichlids <- 
  V2R_DoC_df_info_mean_total %>%
  dplyr::select(Sp, mean_normalized_nb) %>%
  mutate(hist_color = "Cichlid") %>%
  mutate(order = "Cichliformes")

colnames(V2R_cichlids) <- c("Species", "Number_V2R", "Color_hist", "order")

V2R_cichlid_teleost <- 
  rbind(V2R_teleost, V2R_cichlids)


pdf(file = "Histo_TeleostCichlid_V2R.pdf",width = 6.34,  height = 4.61)

V2R_cichlid_teleost %>%
  ggplot(., aes(x=Number_V2R, fill=Color_hist)) +
  geom_histogram(bins=40, alpha = 0.5, position = "identity", color="black") +
  theme_classic() +
  scale_fill_manual(values = cichlid_teleost_color) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position = "none") +
  xlab("Number of V2R genes") +
  ylab("Number of species")

dev.off()

# lets compute mean per teleost, per cichlids, or per family ...

V2R_cichlid_teleost %>%
  group_by(order) %>%
  summarise(number_species = n(),
            mean_V2R = mean(Number_V2R)) %>%
  arrange(desc(mean_V2R))



V2R_cichlid_teleost %>%
  group_by(Color_hist) %>%
  summarise(number_species = n(),
            mean_V2R = mean(Number_V2R)) %>%
  arrange(desc(mean_V2R))


V2R_cichlid_teleost %>%
  filter(Color_hist == "Cichlid") %>%
  pull(Species) %>% unique()


#### Histogram teleost vs cichlid -- T1R  ---------------------------------

T1R_teleost <- 
  teleost_chemoreceptors %>%
  filter(gene_family == "T1R") %>%
  dplyr::select(species, number_T) %>%
  mutate(hist_color = "Teleost")
T1R_teleost <- 
  left_join(T1R_teleost, order_ncbi_sp, by="species") 
colnames(T1R_teleost) <- c("Species", "Number_T1R", "Color_hist", "order")
T1R_cichlids <- 
  T1R_DoC_df_info_mean_total %>%
  dplyr::select(Sp, mean_normalized_nb) %>%
  mutate(hist_color = "Cichlid") %>%
  mutate(order = "Cichliformes")

colnames(T1R_cichlids) <- c("Species", "Number_T1R", "Color_hist", "order")

T1R_cichlid_teleost <- 
  rbind(T1R_teleost, T1R_cichlids)


pdf(file = "Histo_TeleostCichlid_T1R.pdf",width = 6.34,  height = 4.61)

T1R_cichlid_teleost %>%
  ggplot(., aes(x=Number_T1R, fill=Color_hist)) +
  geom_histogram(bins=40, alpha = 0.5, position = "identity", color="black") +
  theme_classic() +
  scale_fill_manual(values = cichlid_teleost_color) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position = "none") +
  xlab("Number of T1R genes") +
  ylab("Number of species")

dev.off()

# lets compute mean per teleost, per cichlids, or per family ...

T1R_cichlid_teleost %>%
  group_by(order) %>%
  summarise(number_species = n(),
            mean_T1R = mean(Number_T1R)) %>%
  arrange(desc(mean_T1R))



T1R_cichlid_teleost %>%
  group_by(Color_hist) %>%
  summarise(number_species = n(),
            mean_T1R = mean(Number_T1R)) %>%
  arrange(desc(mean_T1R))

#### Histogram teleost vs cichlid -- V1R  ---------------------------------

V1R_teleost <- 
  teleost_chemoreceptors %>%
  filter(gene_family == "V1R") %>%
  dplyr::select(species, number_T) %>%
  mutate(hist_color = "Teleost")
V1R_teleost <- 
  left_join(V1R_teleost, order_ncbi_sp, by="species") 
colnames(V1R_teleost) <- c("Species", "Number_V1R", "Color_hist", "order")
V1R_cichlids <- 
  V1R_DoC_df_info_mean_total %>%
  dplyr::select(Sp) %>%
  mutate(number_V1R = 6) %>%
  mutate(hist_color = "Cichlid") %>%
  mutate(order = "Cichliformes")

colnames(V1R_cichlids) <- c("Species", "Number_V1R", "Color_hist", "order")

V1R_cichlid_teleost <- 
  rbind(V1R_teleost, V1R_cichlids)


pdf(file = "Histo_TeleostCichlid_V1R.pdf",width = 6.34,  height = 4.61)

V1R_cichlid_teleost %>%
  ggplot(., aes(x=Number_V1R, fill=Color_hist)) +
  geom_histogram(bins=40, alpha = 0.5, position = "identity", color="black") +
  theme_classic() +
  scale_fill_manual(values = cichlid_teleost_color) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position = "none") +
  xlab("Number of V1R genes") +
  ylab("Number of species")

dev.off()

# lets compute mean per teleost, per cichlids, or per family ...

V1R_cichlid_teleost %>%
  group_by(order) %>%
  summarise(number_species = n(),
            mean_V1R = mean(Number_V1R)) %>%
  arrange(desc(mean_V1R))



V1R_cichlid_teleost %>%
  group_by(Color_hist) %>%
  summarise(number_species = n(),
            mean_V1R = mean(Number_V1R)) %>%
  arrange(desc(mean_V1R))








#### Histogram teleost vs cichlid -- T2R  ---------------------------------

T2R_teleost <- 
  teleost_chemoreceptors %>%
  filter(gene_family == "T2R") %>%
  dplyr::select(species, number_T) %>%
  mutate(hist_color = "Teleost")
T2R_teleost <- 
  left_join(T2R_teleost, order_ncbi_sp, by="species") 
colnames(T2R_teleost) <- c("Species", "Number_T2R", "Color_hist", "order")
T2R_cichlids <- 
  T2R_DoC_df_info_mean_total %>%
  dplyr::select(Sp, mean_normalized_nb) %>%
  mutate(hist_color = "Cichlid") %>%
  mutate(order = "Cichliformes")

colnames(T2R_cichlids) <- c("Species", "Number_T2R", "Color_hist", "order")

T2R_cichlid_teleost <- 
  rbind(T2R_teleost, T2R_cichlids)


pdf(file = "Histo_TeleostCichlid_T2R.pdf",width = 6.34,  height = 4.61)

T2R_cichlid_teleost %>%
  ggplot(., aes(x=Number_T2R, fill=Color_hist)) +
  geom_histogram(bins=40, alpha = 0.5, position = "identity", color="black") +
  theme_classic() +
  scale_fill_manual(values = cichlid_teleost_color) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position = "none") +
  xlab("Number of T2R genes") +
  ylab("Number of species")

dev.off()

# lets compute mean per teleost, per cichlids, or per family ...

T2R_cichlid_teleost %>%
  group_by(order) %>%
  summarise(number_species = n(),
            mean_T2R = mean(Number_T2R)) %>%
  arrange(desc(mean_T2R))



T2R_cichlid_teleost %>%
  group_by(Color_hist) %>%
  summarise(number_species = n(),
            mean_T2R = mean(Number_T2R)) %>%
  arrange(desc(mean_T2R))

#### Histogram teleost vs cichlid -- OLR  ---------------------------------

OLR_teleost <- 
  teleost_chemoreceptors %>%
  filter(gene_family %in% c("OR", "TAAR", "V2R", "V1R")) %>%
  group_by(species) %>%
  summarise(total_number_T = sum(number_T)) %>%
  mutate(hist_color = "Teleost")
OLR_teleost <- 
  left_join(OLR_teleost, order_ncbi_sp, by="species") 
colnames(OLR_teleost) <- c("Species", "Number_OLR", "Color_hist", "order")

OLR_cichlids <- 
  OLR_Doc_mean_sp %>%
  group_by(Sp) %>%
  summarise(total_number_T = sum(mean_normalized_nb)) %>%
  dplyr::select(Sp, total_number_T) %>%
  mutate(hist_color = "Cichlid") %>%
  mutate(order = "Cichliformes")

colnames(OLR_cichlids) <- c("Species", "Number_OLR", "Color_hist", "order")

OLR_cichlid_teleost <- 
  rbind(OLR_teleost, OLR_cichlids)


pdf(file = "Histo_TeleostCichlid_OLR.pdf",width = 6.34,  height = 4.61)

OLR_cichlid_teleost %>%
  ggplot(., aes(x=Number_OLR, fill=Color_hist)) +
  geom_histogram(bins=40, alpha = 0.5, position = "identity", color="black") +
  theme_classic() +
  scale_fill_manual(values = cichlid_teleost_color) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position = "none") +
  xlab("Number of olfactory receptor genes") +
  ylab("Number of species")

dev.off()

# lets compute mean per teleost, per cichlids, or per family ...

OLR_cichlid_teleost %>%
  group_by(order) %>%
  summarise(number_species = n(),
            mean_OLR = mean(Number_OLR)) %>%
  arrange(desc(mean_OLR))



OLR_cichlid_teleost %>%
  group_by(Color_hist) %>%
  summarise(number_species = n(),
            mean_OLR = mean(Number_OLR)) %>%
  arrange(desc(mean_OLR))




#### Histogram teleost vs cichlid -- TR  ---------------------------------

TR_teleost <- 
  teleost_chemoreceptors %>%
  filter(gene_family %in% c("T1R", "T2R")) %>%
  group_by(species) %>%
  summarise(total_number_T = sum(number_T)) %>%
  mutate(hist_color = "Teleost")
TR_teleost <- 
  left_join(TR_teleost, order_ncbi_sp, by="species") 
colnames(TR_teleost) <- c("Species", "Number_TR", "Color_hist", "order")

TR_cichlids <- 
  TR_Doc_mean_sp %>%
  group_by(Sp) %>%
  summarise(total_number_T = sum(mean_normalized_nb)) %>%
  dplyr::select(Sp, total_number_T) %>%
  mutate(hist_color = "Cichlid") %>%
  mutate(order = "Cichliformes")

colnames(TR_cichlids) <- c("Species", "Number_TR", "Color_hist", "order")

TR_cichlid_teleost <- 
  rbind(TR_teleost, TR_cichlids)


pdf(file = "Histo_TeleostCichlid_TR.pdf",width = 6.34,  height = 4.61)

TR_cichlid_teleost %>%
  ggplot(., aes(x=Number_TR, fill=Color_hist)) +
  geom_histogram(bins=40, alpha = 0.5, position = "identity", color="black") +
  theme_classic() +
  scale_fill_manual(values = cichlid_teleost_color) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position = "none") +
  xlab("Number of taste receptor genes") +
  ylab("Number of species")

dev.off()

# lets compute mean per teleost, per cichlids, or per family ...

TR_cichlid_teleost %>%
  group_by(order) %>%
  summarise(number_species = n(),
            mean_TR = mean(Number_TR)) %>%
  arrange(desc(mean_TR))



TR_cichlid_teleost %>%
  group_by(Color_hist) %>%
  summarise(number_species = n(),
            mean_TR = mean(Number_TR)) %>%
  arrange(desc(mean_TR))



#### Map the species tree with the number of OLR/TR/CR -- DOC  ---------------------------------

#Prepare files
All_Doc_mean_sp <- 
  rbind(OR_DoC_df_info_mean_total, TAAR_DoC_df_info_mean_total, 
        V1R_DoC_df_info_mean_total, V2R_DoC_df_info_mean_total,
        T1R_DoC_df_info_mean_total, T2R_DoC_df_info_mean_total)

Tribe_for_tree <- 
  OLR_Doc_mean_sp %>%
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




pdf(file = "Species_tree_OLR.pdf",width = 6.34,  height = 4.61)

ggtree(radiation_tree_wo_Neospl, layout="circular", size=0.6) %<+% Tribe_for_tree +
  aes(color=Tribe) +
  geom_tiplab(size=0.8) +
  scale_color_manual(values = tribes.colors) +
  geom_fruit(
    data=OLR_Doc_mean_sp,
    geom=geom_bar,
    mapping=aes(y=Sp, x=mean_normalized_nb, fill=Subfamily),
    pwidth=0.5,
    stat="identity",
    orientation="y", 
    position="auto",
    offset=0.1,  #distance depuis les tips
    axis.params=list(
      axis="x", 
      text.angle=-45, 
      hjust=0  
    ),
    grid.params=list() 
  ) +
  scale_fill_manual(values = chemoreceptor_family_colors) + 
  theme(legend.position = "none")

dev.off()



pdf(file = "Species_tree_TR.pdf",width = 6.34,  height = 4.61)

ggtree(radiation_tree_wo_Neospl, layout="circular", size=0.6) %<+% Tribe_for_tree +
  aes(color=Tribe) +
  geom_tiplab(size=0.8) +
  scale_color_manual(values = tribes.colors) +
  geom_fruit(
    data=TR_Doc_mean_sp,
    geom=geom_bar,
    mapping=aes(y=Sp, x=mean_normalized_nb, fill=Subfamily),
    pwidth=0.5,
    stat="identity",
    orientation="y", 
    position="auto",
    offset=0.1,  #distance depuis les tips
    axis.params=list(
      axis="x", 
      text.angle=-45, 
      hjust=0  
    ),
    grid.params=list() 
  ) +
  scale_fill_manual(values = chemoreceptor_family_colors) + 
  theme(legend.position = "none")

dev.off()




OLR_Doc_mean_sp$Subfamily <-
  factor(OLR_Doc_mean_sp$Subfamily,
         levels=rev(c("Total_OR", "Total_TAAR", "Total_V1R", 
                  "Total_V2R"
         )))


TR_Doc_mean_sp$Subfamily <-
  factor(TR_Doc_mean_sp$Subfamily,
         levels=rev(c("Total_T1R", "Total_T2R"
         )))



chemoreceptor_family_colors <- 
  c(Total_OR="#E69F00",
    Total_TAAR = "#2995D2", 
    Total_V1R = "#D55E00",
    Total_V2R = "#22A481",
    Total_T1R = "gray65", 
    Total_T2R = "#CAC367")

pdf(file = "Species_tree_OLR_TR.pdf",width = 6.34,  height = 4.61)

ggtree(radiation_tree_wo_Neospl, layout="circular", size=0.4) %<+% Tribe_for_tree +
  aes(color=Tribe) +
  geom_tiplab(size=0.6, offset=0.07) +
  scale_color_manual(values = tribes.colors) +
  geom_fruit(
    data=OLR_Doc_mean_sp,
    geom=geom_bar,
    mapping=aes(y=Sp, x=mean_normalized_nb, fill=Subfamily),
    pwidth=0.5,
    stat="identity",
    orientation="y", 
    position="auto",
    offset=0.16,  #distance depuis les tips
    axis.params=list(
      axis="x", 
      text.angle=-45, 
      hjust=0  
    ),
    grid.params=list() 
  ) +
  geom_fruit(
    data=TR_Doc_mean_sp,
    geom=geom_bar,
    mapping=aes(y=Sp, x=mean_normalized_nb, fill=Subfamily),
    pwidth=0.5,
    stat="identity",
    orientation="y", 
    position="auto",
    offset=0.1,  #distance depuis les tips
    axis.params=list(
      axis="x", 
      text.angle=-45, 
      hjust=0  
    ),
    grid.params=list() 
  ) +
  scale_fill_manual(values = chemoreceptor_family_colors) + 
  theme(legend.position = "none")

dev.off()



## Now lets try to add RNA-seq and morphology data



Sp_RNA_list <- c("Altfas","Asplep","Auldew","Batvit","Benmel","Cphgib","Ctehor",
                 "Erecya","Gnaper","Gnapfe","Gralem","Julorn","Lepatt","Lepelo","Neocau",
                 "Neomul","Neopul","Neosav","Neotet","Ophven","Pcybri","Permic","Petpol",
                 "Simdia","Telvit","Tremar","Xenspi")

species_RNA_df <- 
  OLR_Doc_mean_sp %>%
  dplyr::select(Sp) %>%
  mutate(RNA_data = 
           if_else(Sp %in% Sp_RNA_list,
                   "RNA",
                   "noRNA")
  )

RNA_color <- 
  c(RNA="darkgreen",
    noRNA = "white")


Sp_PacBio_list <- c("Cunlon", "Batmin", "Cphfro", "Cyplep", "Neomul", "Simdia")


species_PacBio_df <- 
  OLR_Doc_mean_sp %>%
  dplyr::select(Sp) %>%
  mutate(Genomic_data = 
           if_else(Sp %in% Sp_PacBio_list,
                   "PacBio",
                   "noPacBio")
  )


PacBio_color <- 
  c(PacBio="darkgreen",
    noPacBio = "white")


species_lamellae_df <-
  left_join(OLR_Doc_mean_sp, Lamellae_per_sp, by="Sp") %>%
  dplyr::select(Sp, Lamellae_number) %>%
  distinct()




pdf(file = "Species_tree_OLR_TR_info.pdf",width = 6.34,  height = 4.61)

ggtree(radiation_tree_wo_Neospl, layout="circular", size=0.4) %<+% Tribe_for_tree +
  aes(color=Tribe) +
  geom_tiplab(size=0.6, offset=0.07) +
  scale_color_manual(values = tribes.colors) +
  new_scale_colour() + 
  geom_fruit(
    data=species_PacBio_df,
    geom=geom_tile,
    mapping=aes(y=Sp, fill=Genomic_data),
    orientation="x",
    size=0.02,
    width=0.3,
    color="white",
    pwidth=0.1,
    position="auto",
    offset=0.16,
    grid.params=list() 
  ) +
  scale_fill_manual(values = PacBio_color) + 
  new_scale_fill() + 
  geom_fruit(
    data=species_RNA_df,
    geom=geom_tile,
    mapping=aes(y=Sp, fill=RNA_data),
    orientation="x",
    size=0.02,
    width=0.3,
    color="white",
    pwidth=0.1,
    position="auto",
    offset=0.1,
    grid.params=list() 
  ) +
  scale_fill_manual(values = RNA_color) + 
  new_scale_fill() + 
  #geom_fruit(
  #  data=species_lamellae_df,
  #  geom=geom_tile,
  #  mapping=aes(y=Sp, fill=Lamellae_number),
  #  orientation="x",
  #  size=0.02,
  #  width=0.3,
  #  color="white",
  #  pwidth=0.1,
  #  position="auto",
  #  offset=0.1,
  #  grid.params=list() 
  #) +
  #scale_fill_viridis(option="magma", direction=-1, na.value="white") + 
  new_scale_fill() + 
  geom_fruit(
    data=OLR_Doc_mean_sp,
    geom=geom_bar,
    mapping=aes(y=Sp, x=mean_normalized_nb, fill=Subfamily),
    pwidth=0.5,
    stat="identity",
    orientation="y", 
    position="auto",
    offset=0.1,  #distance depuis les tips
    axis.params=list(
      axis="x", 
      text.angle=-45, 
      hjust=0  
    ),
    grid.params=list() 
  ) +
  geom_fruit(
    data=TR_Doc_mean_sp,
    geom=geom_bar,
    mapping=aes(y=Sp, x=mean_normalized_nb, fill=Subfamily),
    pwidth=0.5,
    stat="identity",
    orientation="y", 
    position="auto",
    offset=0.1,  #distance depuis les tips
    axis.params=list(
      axis="x", 
      text.angle=-45, 
      hjust=0  
    ),
    grid.params=list() 
  ) +
  scale_fill_manual(values = chemoreceptor_family_colors) +
  theme(legend.position = "none") 

dev.off()





pdf(file = "Species_tree_OLR_TR_legend.pdf",width = 6.34,  height = 4.61)

ggtree(radiation_tree_wo_Neospl, layout="circular", size=0.4) %<+% Tribe_for_tree +
  geom_tiplab(size=0.6, offset=0.07) +
  new_scale_colour() + 
  geom_fruit(
    data=species_PacBio_df,
    geom=geom_tile,
    mapping=aes(y=Sp, fill=Genomic_data),
    orientation="x",
    size=0.02,
    width=0.3,
    color="white",
    pwidth=0.1,
    position="auto",
    offset=0.16,
    grid.params=list() 
  ) +
  scale_fill_manual(values = RNA_color) + 
  new_scale_fill() + 
  geom_fruit(
    data=species_lamellae_df,
    geom=geom_tile,
    mapping=aes(y=Sp, fill=Lamellae_number),
    orientation="x",
    size=0.02,
    width=0.3,
    color="white",
    pwidth=0.1,
    position="auto",
    offset=0.1,
    grid.params=list() 
  ) +
  scale_fill_viridis(option="magma", direction=-1, na.value="white") 

dev.off()


pdf(file = "Species_tree_OLR_TR_legend_tribes.pdf",width = 6.34,  height = 4.61)

ggtree(radiation_tree_wo_Neospl, layout="circular", size=0.4) %<+% Tribe_for_tree +
  aes(color=Tribe) +
  geom_tiplab(size=0.6, offset=0.07) +
  scale_color_manual(values = tribes.colors)

dev.off()


#### Map the species tree with the number of taste receptor -- DOC  ---------------------------------

b1_tree <- read.tree("b1_tree_filt.nwk")
b1_list_sp <- b1_tree$tip.label
b1_tree_wo_Neospl <- 
  drop.tip(b1_tree, "Neospl") 

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

pdf(file = "Species_tree_OLR.pdf",width = 6.34,  height = 4.61)

ggtree(b1_tree_wo_Neospl, layout="circular", size=0.6) %<+% Tribe_for_tree +
  aes(color=Tribe) +
  geom_tiplab(size=0.8) +
  scale_color_manual(values = tribes.colors) +
  geom_fruit(
    data=OLR_Doc_mean_sp,
    geom=geom_bar,
    mapping=aes(y=Sp, x=mean_normalized_nb, fill=Subfamily),
    pwidth=0.5,
    stat="identity",
    orientation="y", 
    position="auto",
    offset=0.1,  #distance depuis les tips
    axis.params=list(
      axis="x", 
      text.angle=-45, 
      hjust=0  
    ),
    grid.params=list() 
  ) +
  scale_fill_manual(values = chemoreceptor_family_colors) + 
  theme(legend.position = "none")

dev.off()


pdf(file = "Species_tree_TR.pdf",width = 6.34,  height = 4.61)

ggtree(b1_tree_wo_Neospl, layout="circular", size=0.6) %<+% Tribe_for_tree +
  aes(color=Tribe) +
  geom_tiplab(size=0.8) +
  scale_color_manual(values = tribes.colors) +
  geom_fruit(
    data=TR_Doc_mean_sp,
    geom=geom_bar,
    mapping=aes(y=Sp, x=mean_normalized_nb, fill=Subfamily),
    pwidth=0.5,
    stat="identity",
    orientation="y", 
    position="auto",
    offset=0.1,  #distance depuis les tips
    axis.params=list(
      axis="x", 
      text.angle=-45, 
      hjust=0  
    ),
    grid.params=list() 
  ) +
  scale_fill_manual(values = chemoreceptor_family_colors) + 
  theme(legend.position = "none")

dev.off()






#### Map the species tree with the number of OR genes -- DOC  ---------------------------------



OR_DoC_df_info_mean_wototal <- OR_DoC_df_info_mean %>% filter(Subfamily != "Total_OR")

library(RColorBrewer)
n <- length(OR_DoC_df_info_mean_wototal %>% pull(Subfamily) %>% unique())
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector <- col_vector[1:47]
vector_subfam <- as.vector(OR_DoC_df_info_mean_wototal %>% pull(Subfamily) %>% unique())
names(col_vector) <- vector_subfam

OR_DoC_df_info_mean_wototal$Subfamily <- 
  factor(OR_DoC_df_info_mean_wototal$Subfamily, 
         levels=rev(c("Lambda-1","Kappa-1","Theta-1","Eta-1","Eta-2","Eta-3","Eta-4","Eta-6","Eta-7","Eta-8","Eta-9","Eta-10","Eta-11","Eta-12","Eta-13","Eta-14","Eta-15","Beta-1","Beta-2","Epsilon-1","Epsilon-2","Epsilon-3","Epsilon-4","Zeta-1","Zeta-2","Zeta-3","Zeta-4","Delta-1","Delta-2","Delta-3","Delta-4","Delta-5","Delta-6","Delta-7","Delta-8","Delta-9","Delta-10","Delta-11","Delta-12","Delta-13","Delta-14","Delta-15","Delta-16","Delta-17","Delta-18","Delta-19","Delta-20")))


pdf(file = "Species_tree_OR.pdf",width = 6.34,  height = 4.61)

ggtree(radiation_tree_wo_Neospl, layout="circular", size=0.4) %<+% Tribe_for_tree +
  aes(color=Tribe) +
  geom_tiplab(size=0.6, offset=0.07) +
  scale_color_manual(values = tribes.colors) +
  geom_fruit(
    data=OR_DoC_df_info_mean_wototal,
    geom=geom_bar,
    mapping=aes(y=Sp, x=mean_normalized_nb, fill=Subfamily),
    pwidth=0.5,
    stat="identity",
    orientation="y", 
    position="auto",
    offset=0.16,  #distance depuis les tips
    axis.params=list(
      axis="x", 
      text.angle=-45, 
      hjust=0  
    ),
    grid.params=list() 
  ) +
  scale_fill_manual(values = col_vector) + 
  theme(legend.position = "none")

dev.off()

#### Map the species tree with the number of TAAR genes -- DOC  ---------------------------------

TAAR_DoC_df_info_mean_wototal <- TAAR_DoC_df_info_mean %>% filter(Subfamily != "Total_TAAR")
TAAR_DoC_df_info_mean_wototal$Subfamily <- 
  factor(TAAR_DoC_df_info_mean_wototal$Subfamily, levels=rev(c("TAARL-1",'TAARA2-1',"TAARA2-2",'TAARB4-1',"TAARB4-2","TAARB4-3",'TAARB4-4','TAARB4-5')))

pdf(file = "Species_tree_TAAR.pdf",width = 6.34,  height = 4.61)

ggtree(radiation_tree_wo_Neospl, layout="circular", size=0.4) %<+% Tribe_for_tree +
  aes(color=Tribe) +
  geom_tiplab(size=0.6, offset=0.07) +
  scale_color_manual(values = tribes.colors) +
  geom_fruit(
    data=TAAR_DoC_df_info_mean_wototal,
    geom=geom_bar,
    mapping=aes(y=Sp, x=mean_normalized_nb, fill=Subfamily),
    pwidth=0.5,
    stat="identity",
    orientation="y", 
    position="auto",
    offset=0.16,  #distance depuis les tips
    axis.params=list(
      axis="x", 
      text.angle=-45, 
      hjust=0  
    ),
    grid.params=list() 
  ) +
  theme(legend.position = "none") +
  scale_fill_manual(values = TAAR_subfam_colors)

dev.off()


#### Map the species tree with the number of V2R genes -- DOC  ---------------------------------

V2R_DoC_df_info_mean_wototal <- V2R_DoC_df_info_mean %>% filter(Subfamily != "Total_V2R")
V2R_DoC_df_info_mean_wototal$Subfamily <- 
  factor(V2R_DoC_df_info_mean_wototal$Subfamily, levels=rev(c("V2RD2-1", "V2RD3-1", "V2RD4-1", "V2RD5-1", "V2RD7-1", "V2RD8-1", "V2RD9-1", "V2RD10-1", "V2RD11-1", "V2RD11-2", "V2RD11-3", "V2RD11-4", "V2RD11-5", "V2RD11-6", "V2RD11-7")))

pdf(file = "Species_tree_V2R.pdf",width = 6.34,  height = 4.61)

ggtree(radiation_tree_wo_Neospl, layout="circular", size=0.4) %<+% Tribe_for_tree +
  aes(color=Tribe) +
  geom_tiplab(size=0.6, offset=0.07) +
  scale_color_manual(values = tribes.colors) +
  geom_fruit(
    data=V2R_DoC_df_info_mean_wototal,
    geom=geom_bar,
    mapping=aes(y=Sp, x=mean_normalized_nb, fill=Subfamily),
    pwidth=0.5,
    stat="identity",
    orientation="y", 
    position="auto",
    offset=0.16,  #distance depuis les tips
    axis.params=list(
      axis="x", 
      text.angle=-45, 
      hjust=0  
    ),
    grid.params=list() 
  ) +
  theme(legend.position = "none") +
  scale_fill_manual(values = V2R_subfam_colors)

dev.off()



#### Map the species tree with the number of V1R genes -- DOC  ---------------------------------

V1R_DoC_df_info_mean_wototal <- V1R_DoC_df_info_mean %>% filter(Subfamily != "Total_V1R")
V1R_DoC_df_info_mean_wototal$Subfamily <- 
  factor(V1R_DoC_df_info_mean_wototal$Subfamily, levels=rev(c("ORA1", "ORA2", "ORA3", "ORA4", "ORA5", "ORA6")))

pdf(file = "Species_tree_V1R.pdf",width = 6.34,  height = 4.61)

ggtree(radiation_tree_wo_Neospl, layout="circular", size=0.4) %<+% Tribe_for_tree +
  aes(color=Tribe) +
  geom_tiplab(size=0.6, offset=0.07) +
  scale_color_manual(values = tribes.colors) +
  geom_fruit(
    data=V1R_DoC_df_info_mean_wototal,
    geom=geom_bar,
    mapping=aes(y=Sp, x=mean_normalized_nb, fill=Subfamily),
    pwidth=0.5,
    stat="identity",
    orientation="y", 
    position="auto",
    offset=0.16,  #distance depuis les tips
    axis.params=list(
      axis="x", 
      text.angle=-45, 
      hjust=0  
    ),
    grid.params=list() 
  ) +
  theme(legend.position = "none") +
  scale_fill_manual(values = V1R_subfam_colors)

dev.off()




#### Map the species tree with the number of T1R genes -- DOC  ---------------------------------

T1R_DoC_df_info_mean_wototal <- T1R_DoC_df_info_mean %>% filter(Subfamily != "Total_T1R")
T1R_DoC_df_info_mean_wototal$Subfamily <- 
  factor(T1R_DoC_df_info_mean_wototal$Subfamily, levels=rev(c("T1R1A", "T1R1B", "T1R3")))

pdf(file = "Species_tree_T1R.pdf",width = 6.34,  height = 4.61)

ggtree(radiation_tree_wo_Neospl, layout="circular", size=0.4) %<+% Tribe_for_tree +
  aes(color=Tribe) +
  geom_tiplab(size=0.6, offset=0.07) +
  scale_color_manual(values = tribes.colors) +
  geom_fruit(
    data=T1R_DoC_df_info_mean_wototal,
    geom=geom_bar,
    mapping=aes(y=Sp, x=mean_normalized_nb, fill=Subfamily),
    pwidth=0.5,
    stat="identity",
    orientation="y", 
    position="auto",
    offset=0.16,  #distance depuis les tips
    axis.params=list(
      axis="x", 
      text.angle=-45, 
      hjust=0  
    ),
    grid.params=list() 
  ) +
  theme(legend.position = "none") +
  scale_fill_manual(values = T1R_subfam_colors)

dev.off()

#### Map the species tree with the number of T2R genes -- DOC  ---------------------------------

T2R_DoC_df_info_mean_wototal <- T2R_DoC_df_info_mean %>% filter(Subfamily != "Total_T2R")
T2R_DoC_df_info_mean_wototal$Subfamily <- 
  factor(T2R_DoC_df_info_mean_wototal$Subfamily, levels=rev(c("T2RA")))

pdf(file = "Species_tree_T2R.pdf",width = 6.34,  height = 4.61)

ggtree(radiation_tree_wo_Neospl, layout="circular", size=0.4) %<+% Tribe_for_tree +
  aes(color=Tribe) +
  geom_tiplab(size=0.6, offset=0.07) +
  scale_color_manual(values = tribes.colors) +
  geom_fruit(
    data=T2R_DoC_df_info_mean_wototal,
    geom=geom_bar,
    mapping=aes(y=Sp, x=mean_normalized_nb, fill=Subfamily),
    pwidth=0.5,
    stat="identity",
    orientation="y", 
    position="auto",
    offset=0.16,  #distance depuis les tips
    axis.params=list(
      axis="x", 
      text.angle=-45, 
      hjust=0  
    ),
    grid.params=list() 
  ) +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("T2RA" = "#004D40"))

dev.off()


#### Compute some statistics - OLR  ---------------------------------



Stats_OLR <- 
  All_Doc_mean_sp_wide %>%
  filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>%
  mutate(Total_OLR = Total_OR + Total_TAAR + Total_V1R + Total_V2R) %>%
  summarise(
    mean_OLR = mean(Total_OLR),
    min_OLR = min(Total_OLR),
    max_OLR = max(Total_OLR)
  )
min_OLR <- Stats_OLR %>% pull(min_OLR)
max_OLR <- Stats_OLR %>% pull(max_OLR)
mean_OLR <- Stats_OLR %>% pull(mean_OLR)

All_Doc_mean_sp_wide %>%
  filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>%
  mutate(Total_OLR = Total_OR + Total_TAAR + Total_V1R + Total_V2R) %>%
  filter(Total_OLR %in% c(min_OLR, max_OLR)) %>%
  dplyr::select(Sp, Species, Total_OLR)




#### Compute some statistics - TR  ---------------------------------

#All_Doc_mean_sp_wide <- 
#  as.data.frame(
#    pivot_wider(All_Doc_mean_sp, 
#                names_from = c(Subfamily), values_from = mean_normalized_nb)
#  )



Stats_TR <- 
  All_Doc_mean_sp_wide %>%
  filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>%
  mutate(Total_TR = Total_T1R + Total_T2R) %>%
  summarise(
    mean_TR = mean(Total_TR),
    min_TR = min(Total_TR),
    max_TR = max(Total_TR)
  )
min_TR <- Stats_TR %>% pull(min_TR)
max_TR <- Stats_TR %>% pull(max_TR)
mean_TR <- Stats_TR %>% pull(mean_TR)

All_Doc_mean_sp_wide %>%
  filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>%
  mutate(Total_TR = Total_T1R + Total_T2R) %>%
  filter(Total_TR %in% c(min_TR, max_TR)) %>%
  dplyr::select(Sp, Species, Total_TR)



#### Compute some statistics - OR  ---------------------------------

Stats_OR <- 
  OR_DoC_df_info_mean_total %>%
  filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>%
  summarise(
    mean_OR = mean(mean_normalized_nb),
    min_OR = min(mean_normalized_nb),
    max_OR = max(mean_normalized_nb)
  )
min_OR <- Stats_OR %>% pull(min_OR)
max_OR <- Stats_OR %>% pull(max_OR)

OR_DoC_df_info_mean_total %>%
  filter(mean_normalized_nb %in% c(min_OR, max_OR)) %>%
  dplyr::select(Sp, Species, mean_normalized_nb)



#Which are the subfamilies in single copy ?

as.data.frame(
  OR_DoC_df_info_mean %>%
    group_by(Subfamily) %>%
    filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>%
    summarise(
      mean = mean(mean_normalized_nb),
      min = min(mean_normalized_nb),
      max = max(mean_normalized_nb)
    )) %>%
  mutate(diff = max - min) %>%
  arrange(diff)



  
OR_teleost %>%
  summarise(
    mean_OR = mean(Number_OR),
    min_OR = min(Number_OR),
    max_OR = max(Number_OR)
  )

#### Compute some statistics - TAAR  ---------------------------------

Stats_TAAR <- 
  TAAR_DoC_df_info_mean_total %>%
  filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>%
  summarise(
    mean_TAAR = mean(mean_normalized_nb),
    min_TAAR = min(mean_normalized_nb),
    max_TAAR = max(mean_normalized_nb)
  )
min_TAAR <- Stats_TAAR %>% pull(min_TAAR)
max_TAAR <- Stats_TAAR %>% pull(max_TAAR)

TAAR_DoC_df_info_mean_total %>%
  filter(mean_normalized_nb %in% c(min_TAAR, max_TAAR)) %>%
  dplyr::select(Sp, Species, mean_normalized_nb)



#Which are the subfamilies in single copy ?

as.data.frame(
  TAAR_DoC_df_info_mean %>%
    group_by(Subfamily) %>%
    filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>%
    summarise(
      mean = mean(mean_normalized_nb),
      min = min(mean_normalized_nb),
      max = max(mean_normalized_nb)
    )) %>%
  mutate(diff = max - min) %>%
  arrange(diff)


TAAR_teleost %>%
  summarise(
    mean_TAAR = mean(Number_TAAR),
    min_TAAR = min(Number_TAAR),
    max_TAAR = max(Number_TAAR)
  )

#### Compute some statistics - V2R  ---------------------------------

Stats_V2R <- 
  V2R_DoC_df_info_mean_total %>%
  filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>%
  summarise(
    mean_V2R = mean(mean_normalized_nb),
    min_V2R = min(mean_normalized_nb),
    max_V2R = max(mean_normalized_nb)
  )
min_V2R <- Stats_V2R %>% pull(min_V2R)
max_V2R <- Stats_V2R %>% pull(max_V2R)

V2R_DoC_df_info_mean_total %>%
  filter(mean_normalized_nb %in% c(min_V2R, max_V2R)) %>%
  dplyr::select(Sp, Species, mean_normalized_nb)


#Which are the subfamilies in single copy ?

as.data.frame(
  V2R_DoC_df_info_mean %>%
    group_by(Subfamily) %>%
    filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>%
    summarise(
      mean = mean(mean_normalized_nb),
      min = min(mean_normalized_nb),
      max = max(mean_normalized_nb)
    ))  %>%
  mutate(diff = max - min) %>%
  arrange(diff)


V2R_teleost %>%
  summarise(
    mean_V2R = mean(Number_V2R),
    min_V2R = min(Number_V2R),
    max_V2R = max(Number_V2R)
  )

#### Compute some statistics - T1R  ---------------------------------

Stats_T1R <- 
  T1R_DoC_df_info_mean_total %>%
  filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>%
  summarise(
    mean_T1R = mean(mean_normalized_nb),
    min_T1R = min(mean_normalized_nb),
    max_T1R = max(mean_normalized_nb)
  )
min_T1R <- Stats_T1R %>% pull(min_T1R)
max_T1R <- Stats_T1R %>% pull(max_T1R)

T1R_DoC_df_info_mean_total %>%
  filter(mean_normalized_nb %in% c(min_T1R, max_T1R)) %>%
  dplyr::select(Sp, Species, mean_normalized_nb)

#### Correlations between gene number in families  --------


All_Doc_mean_sp_wide <- 
  All_Doc_mean_sp_wide %>%
  mutate(Total_OLR = Total_OR + Total_TAAR  + Total_V2R + Total_V1R) %>%
  mutate(Total_TR = Total_T1R + Total_T2R)

#The caper data has 264 species because we removed Neospl
caper_data_cichlids <- 
  comparative.data(phy = radiation_tree_wo_Neospl, 
                   data = All_Doc_mean_sp_wide,
                   names.col = Sp, vcv = TRUE,
                   na.omit = FALSE, warn.dropped = TRUE)



## Test the phylogenetic signal of each chemoreceptor family
tip_order <- radiation_tree_wo_Neospl$tip.label

Ordered_data_receptor <- 
  enframe(tip_order, name = NULL, value = "Sp") %>% 
  left_join(All_Doc_mean_sp_wide) %>% 
  column_to_rownames("Sp")



#### OR vs V2R

fit_phylo_OR_V2R <- pgls(Total_V2R ~ Total_OR,
                         data = caper_data_cichlids, lambda = "ML",
                         bounds=list(lambda=c(0,1)))

sum_fit_phy <- summary(fit_phylo_OR_V2R)
pente =   formatC(sum_fit_phy$coefficients[2], digits = 3)
r2 =      formatC( sum_fit_phy$adj.r.squared, digits = 2)
pvalue =  formatC(sum_fit_phy$coefficients[8], digits = 3)
PGLS_cc <- coef(fit_phylo_OR_V2R)
PGLS_lambda <- sum_fit_phy$param[2]
fit_function <- function(x) PGLS_cc[1] + PGLS_cc[2]*x
r2
if (pvalue == "   0"){ pvalue = "< 2.2e-16"}
lm_fit <- lm(data = All_Doc_mean_sp_wide %>% filter(Sp %in% radiation_tree_wo_Neospl$tip.label), formula = Total_V2R~ Total_OR)
GLS_cc <- coef(lm_fit)
GLS_fit_function <- function(x) GLS_cc[1] + GLS_cc[2]*x

pdf(file = "OR_vs_V2R.pdf",width = 7.34,  height = 4.61)

All_Doc_mean_sp_wide %>%
  filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>%
  filter(Sp != "Neospl") %>%
  ggplot(., aes( y = Total_V2R, x= Total_OR))+
  geom_point(aes(color=Tribe))+
  stat_function(fun = GLS_fit_function, color="black")+
  labs(subtitle = bquote(R^2 ~ "=" ~ .(r2) ~ ";" ~ P ~ .(pvalue) ~ ";" ~ n ~ "=" ~ "264" ~ ";" ~ lambda ~ "=" ~ .(PGLS_lambda))) +
  theme_classic() +
  ylab("Number of V2R genes") +
  xlab("Number of OR genes") +
  #geom_abline(slope=1, color="black", linetype="dashed") + 
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none") +
  scale_size_identity() +
  scale_color_manual(values=tribes.colors)

dev.off()


#### OR vs TAAR

fit_phylo_OR_TAAR <- pgls(Total_TAAR ~ Total_OR,
                          data = caper_data_cichlids, lambda = "ML")

sum_fit_phy <- summary(fit_phylo_OR_TAAR)
pente =   formatC(sum_fit_phy$coefficients[2], digits = 3)
r2 =      formatC( sum_fit_phy$adj.r.squared, digits = 2)
pvalue =  formatC(sum_fit_phy$coefficients[8], digits = 3)
PGLS_cc <- coef(fit_phylo_OR_TAAR)
PGLS_lambda <- sum_fit_phy$param[2]
fit_function <- function(x) PGLS_cc[1] + PGLS_cc[2]*x
r2
if (pvalue == "   0"){ pvalue = "< 2.2e-16"}
lm_fit <- lm(data = All_Doc_mean_sp_wide %>% filter(Sp %in% radiation_tree_wo_Neospl$tip.label), formula = Total_TAAR~ Total_OR)
GLS_cc <- coef(lm_fit)
GLS_fit_function <- function(x) GLS_cc[1] + GLS_cc[2]*x

pdf(file = "OR_vs_TAAR.pdf",width = 7.34,  height = 4.61)

All_Doc_mean_sp_wide %>%
  filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>%
  filter(Sp != "Neospl") %>%
  ggplot(., aes( y = Total_TAAR, x= Total_OR))+
  geom_point(aes(color=Tribe))+
  stat_function(fun = GLS_fit_function, color="black")+
  labs(subtitle = bquote(R^2 ~ "=" ~ .(r2) ~ ";" ~ P ~ "=" ~ .(pvalue) ~ ";" ~ n ~ "=" ~ "264" ~ ";" ~ lambda ~ "=" ~ .(PGLS_lambda))) +
  theme_classic() +
  ylab("Number of TAAR genes") +
  xlab("Number of OR genes") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none") +
  #geom_abline(slope=1, color="black", linetype="dashed") + 
  scale_size_identity() +
  scale_color_manual(values=tribes.colors)

dev.off()


#### V2R vs TAAR


fit_phylo_V2R_TAAR <- pgls(Total_TAAR ~ Total_V2R,
                           data = caper_data_cichlids, lambda = "ML")

sum_fit_phy <- summary(fit_phylo_V2R_TAAR)
pente =   formatC(sum_fit_phy$coefficients[2], digits = 3)
r2 =      formatC( sum_fit_phy$r.squared, digits = 2)
pvalue =  formatC(sum_fit_phy$coefficients[8], digits = 3)
PGLS_cc <- coef(fit_phylo_V2R_TAAR)
PGLS_lambda <- sum_fit_phy$param[2]
fit_function <- function(x) PGLS_cc[1] + PGLS_cc[2]*x
r2
if (pvalue == "   0"){ pvalue = "< 2.2e-16"}
lm_fit <- lm(data = All_Doc_mean_sp_wide %>% filter(Sp %in% radiation_tree_wo_Neospl$tip.label), formula = Total_TAAR~ Total_V2R)
GLS_cc <- coef(lm_fit)
GLS_fit_function <- function(x) GLS_cc[1] + GLS_cc[2]*x

pdf(file = "V2R_vs_TAAR.pdf",width = 7.34,  height = 4.61)

All_Doc_mean_sp_wide %>%
  filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>%
  filter(Sp != "Neospl") %>%
  ggplot(., aes( y = Total_TAAR, x= Total_V2R))+
  geom_point(aes(color=Tribe))+
  stat_function(fun = GLS_fit_function, color="black")+
  labs(subtitle = bquote(R^2 ~ "=" ~ .(r2) ~ ";" ~ P ~ "=" ~ .(pvalue) ~ ";" ~ n ~ "=" ~ "264" ~ ";" ~ lambda ~ "=" ~ .(PGLS_lambda))) +
  theme_classic() +
  ylab("Number of TAAR genes") +
  xlab("Number of V2R genes") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none") +
  #geom_abline(slope=1, color="black", linetype="dashed") + 
  scale_size_identity() +
  scale_color_manual(values=tribes.colors)

dev.off()



#### Taste receptor vs Olfactory receptor


fit_phylo_TR_OLR <- pgls(Total_TR ~ Total_OLR,
                         data = caper_data_cichlids, lambda = "ML")

sum_fit_phy <- summary(fit_phylo_TR_OLR)
pente =   formatC(sum_fit_phy$coefficients[2], digits = 3)
r2 =      formatC( sum_fit_phy$r.squared, digits = 2)
pvalue =  formatC(sum_fit_phy$coefficients[8], digits = 3)
PGLS_cc <- coef(fit_phylo_TR_OLR)
PGLS_lambda <- sum_fit_phy$param[2]
fit_function <- function(x) PGLS_cc[1] + PGLS_cc[2]*x
r2
if (pvalue == "   0"){ pvalue = "< 2.2e-16"}
lm_fit <- lm(data = All_Doc_mean_sp_wide %>% filter(Sp %in% radiation_tree_wo_Neospl$tip.label), formula = Total_TR ~ Total_OLR)
GLS_cc <- coef(lm_fit)
GLS_fit_function <- function(x) GLS_cc[1] + GLS_cc[2]*x
GLS_fit_function <- function(x) PGLS_cc[1] + PGLS_cc[2]*x

pdf(file = "OLR_vs_TR.pdf",width = 7.34,  height = 4.61)

All_Doc_mean_sp_wide %>%
  filter(Sp %in% radiation_tree_wo_Neospl$tip.label) %>%
  filter(Sp != "Neospl") %>%
  ggplot(., aes( y = Total_TR, x= Total_OLR))+
  geom_point(aes(color=Tribe))+
  stat_function(fun = GLS_fit_function, color="black")+
  labs(subtitle = bquote(R^2 ~ "=" ~ .(r2) ~ ";" ~ P ~ "=" ~ .(pvalue) ~ ";" ~ n ~ "=" ~ "264" ~ ";" ~ lambda ~ "=" ~ .(PGLS_lambda))) +
  theme_classic() +
  ylab("Number of taste receptor genes") +
  xlab("Number of olfactory receptor genes") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none") +
  #geom_abline(slope=1, color="black", linetype="dashed") + 
  scale_size_identity() +
  scale_color_manual(values=tribes.colors)

dev.off()


#### PCA OR nb in  subfamilies   ---------------------------------


OR_DoC_df_info_mean_wo_total <- 
  OR_DoC_df_info_mean %>%
  filter(Subfamily != "Total_OR")

OR_Doc_sp_wide <- 
  as.data.frame(
    pivot_wider(OR_DoC_df_info_mean_wo_total, 
                names_from = c(Subfamily), values_from = mean_normalized_nb)
  )


species_in_tree <- b1_tree$tip.label
OR_Doc_sp_wide <- 
  OR_Doc_sp_wide %>% 
  filter(Sp != "Neospl") %>% 
  filter(Sp %in% radiation_tree_wo_Neospl$tip.label)


OR_Doc_sp_wide_woinfo_radiation <- 
  OR_Doc_sp_wide_woinfo %>%
  filter

OR_Doc_sp_wide_woinfo <- 
  OR_Doc_sp_wide %>%
  dplyr::select(- c(Species, Tribe))

rownames(OR_Doc_sp_wide_woinfo) <- OR_Doc_sp_wide_woinfo[,1]
OR_Doc_sp_wide_woinfo <- OR_Doc_sp_wide_woinfo[,-1]

tribe <- as.factor(OR_Doc_sp_wide$Tribe[1:nrow(OR_Doc_sp_wide)])
tribes.colors <- c(Bathybatini = "#242626" , Benthochromini = "#AE262A" , Boulengerochromini = "#59595C" , Cyphotilapiini = "#FDDF13" , Cyprichromini = "#F04D29" , Ectodini = "#9AB9D9" , Eretmodini = "#682E7A" , Lamprologini = "#C588BB" , Limnochromini = "#535CA9" , Perissodini = "#FBAC43" , Trematocarini = "#959170" , Tropheini = "#86C773" , Haplochromini = "#274e13" , outgroup = "gray", Other="gray")


PCA_OR <- prcomp(OR_Doc_sp_wide_woinfo, scale = TRUE)


autoplot(PCA_OR, data = OR_Doc_sp_wide, colour = 'Tribe') +
  scale_color_manual(values = tribes.colors) +
  theme_classic()


PC1 <- PCA_OR$x[,1]
PC2 <- PCA_OR$x[,2]

ggplot(OR_Doc_sp_wide,
       aes(x = PC1,
           y = PC2,
           color=Tribe)) +
  geom_point() +
  stat_ellipse(level = 0.95) +
  scale_color_manual(values = tribes.colors) +
  theme_classic()



OR_Doc_sp_wide_ecocat <- left_join(OR_Doc_sp_wide, ecoCat, by="Sp")
OR_Doc_sp_wide_ecocat_breeding <- 
  left_join(OR_Doc_sp_wide_ecocat, breed_infos, by="Sp")

ggplot(OR_Doc_sp_wide,
       aes(x = PC1,
           y = PC2,
           color=Tribe)) +
  geom_point() +
  stat_ellipse(level = 0.95) +
  scale_color_manual(values = tribes.colors) +
  theme_classic()

ggplot(OR_Doc_sp_wide_ecocat,
       aes(x = PC1,
           y = PC2,
           color=Food)) +
  geom_point() +
  stat_ellipse(level = 0.95) +
  theme_classic()


ggplot(OR_Doc_sp_wide_ecocat,
       aes(x = PC1,
           y = PC2,
           color=habitat)) +
  geom_point() +
  stat_ellipse(level = 0.95) +
  theme_classic()



ggplot(OR_Doc_sp_wide_ecocat_breeding,
       aes(x = PC1,
           y = PC2,
           color=breeding_type)) +
  geom_point() +
  stat_ellipse(level = 0.95) +
  theme_classic()

ggplot(OR_Doc_sp_wide_ecocat_breeding,
       aes(x = PC1,
           y = PC2,
           color=breeding_mode)) +
  geom_point() +
  stat_ellipse(level = 0.95) +
  theme_classic()



#### PCA TAAR nb in  subfamilies   ---------------------------------


TAAR_DoC_df_info_mean_wo_total <- 
  TAAR_DoC_df_info_mean %>%
  filter(Subfamily != "Total_TAAR")

TAAR_Doc_sp_wide <- 
  as.data.frame(
    pivot_wider(TAAR_DoC_df_info_mean_wo_total, 
                names_from = c(Subfamily), values_from = mean_normalized_nb)
  )


species_in_tree <- b1_tree$tip.label
TAAR_Doc_sp_wide <- 
  TAAR_Doc_sp_wide %>% 
  filter(Sp != "Neospl") %>% 
  filter(Sp %in% radiation_tree_wo_Neospl$tip.label)



TAAR_Doc_sp_wide_woinfo <- 
  TAAR_Doc_sp_wide %>%
  dplyr::select(- c(Species, Tribe))

rownames(TAAR_Doc_sp_wide_woinfo) <- TAAR_Doc_sp_wide_woinfo[,1]
TAAR_Doc_sp_wide_woinfo <- TAAR_Doc_sp_wide_woinfo[,-1]

tribe <- as.factor(TAAR_Doc_sp_wide$Tribe[1:nrow(TAAR_Doc_sp_wide)])
tribes.colors <- c(Bathybatini = "#242626" , Benthochromini = "#AE262A" , Boulengerochromini = "#59595C" , Cyphotilapiini = "#FDDF13" , Cyprichromini = "#F04D29" , Ectodini = "#9AB9D9" , Eretmodini = "#682E7A" , Lamprologini = "#C588BB" , Limnochromini = "#535CA9" , Perissodini = "#FBAC43" , Trematocarini = "#959170" , Tropheini = "#86C773" , Haplochromini = "#274e13" , outgroup = "gray", Other="gray")


PCA_TAAR <- prcomp(TAAR_Doc_sp_wide_woinfo, scale = TRUE)

autoplot(PCA_TAAR, data = TAAR_Doc_sp_wide, colour = 'Tribe') +
  scale_color_manual(values = tribes.colors) +
  theme_classic()


PC1 <- PCA_TAAR$x[,1]
PC2 <- PCA_TAAR$x[,2]

ggplot(TAAR_Doc_sp_wide,
       aes(x = PC1,
           y = PC2,
           color=Tribe)) +
  geom_point() +
  stat_ellipse(level = 0.95) +
  scale_color_manual(values = tribes.colors) +
  theme_classic()



TAAR_Doc_sp_wide_ecocat <- left_join(TAAR_Doc_sp_wide, ecoCat, by="Sp")
TAAR_Doc_sp_wide_ecocat_breeding <- 
  left_join(TAAR_Doc_sp_wide_ecocat, breed_infos, by="Sp")

ggplot(TAAR_Doc_sp_wide,
       aes(x = PC1,
           y = PC2,
           color=Tribe)) +
  geom_point() +
  stat_ellipse(level = 0.95) +
  scale_color_manual(values = tribes.colors) +
  theme_classic()

ggplot(TAAR_Doc_sp_wide_ecocat,
       aes(x = PC1,
           y = PC2,
           color=Food)) +
  geom_point() +
  stat_ellipse(level = 0.95) +
  theme_classic()


ggplot(TAAR_Doc_sp_wide_ecocat,
       aes(x = PC1,
           y = PC2,
           color=habitat)) +
  geom_point() +
  stat_ellipse(level = 0.95) +
  theme_classic()



ggplot(TAAR_Doc_sp_wide_ecocat_breeding,
       aes(x = PC1,
           y = PC2,
           color=breeding_type)) +
  geom_point() +
  stat_ellipse(level = 0.95) +
  theme_classic()

ggplot(TAAR_Doc_sp_wide_ecocat_breeding,
       aes(x = PC1,
           y = PC2,
           color=breeding_mode)) +
  geom_point() +
  stat_ellipse(level = 0.95) +
  theme_classic()


#### PCA V2R nb in  subfamilies   ---------------------------------


V2R_DoC_df_info_mean_wo_total <- 
  V2R_DoC_df_info_mean %>%
  filter(Subfamily != "Total_V2R")

V2R_Doc_sp_wide <- 
  as.data.frame(
    pivot_wider(V2R_DoC_df_info_mean_wo_total, 
                names_from = c(Subfamily), values_from = mean_normalized_nb)
  )


species_in_tree <- b1_tree$tip.label
V2R_Doc_sp_wide <- 
  V2R_Doc_sp_wide %>% 
  filter(Sp != "Neospl") %>% 
  filter(Sp %in% radiation_tree_wo_Neospl$tip.label)

V2R_Doc_sp_wide_woinfo <- 
  V2R_Doc_sp_wide %>%
  dplyr::select(- c(Species, Tribe))

rownames(V2R_Doc_sp_wide_woinfo) <- V2R_Doc_sp_wide_woinfo[,1]
V2R_Doc_sp_wide_woinfo <- V2R_Doc_sp_wide_woinfo[,-1]

tribe <- as.factor(V2R_Doc_sp_wide$Tribe[1:nrow(V2R_Doc_sp_wide)])
tribes.colors <- c(Bathybatini = "#242626" , Benthochromini = "#AE262A" , Boulengerochromini = "#59595C" , Cyphotilapiini = "#FDDF13" , Cyprichromini = "#F04D29" , Ectodini = "#9AB9D9" , Eretmodini = "#682E7A" , Lamprologini = "#C588BB" , Limnochromini = "#535CA9" , Perissodini = "#FBAC43" , Trematocarini = "#959170" , Tropheini = "#86C773" , Haplochromini = "#274e13" , outgroup = "gray", Other="gray")


PCA_V2R <- prcomp(V2R_Doc_sp_wide_woinfo, scale = TRUE)

autoplot(PCA_V2R, data = V2R_Doc_sp_wide, colour = 'Tribe') +
  scale_color_manual(values = tribes.colors) +
  theme_classic()


PC1 <- PCA_V2R$x[,1]
PC2 <- PCA_V2R$x[,2]

ggplot(V2R_Doc_sp_wide,
       aes(x = PC1,
           y = PC2,
           color=Tribe,
           label=Sp)) +
  geom_point() +
  stat_ellipse(level = 0.95) +
  scale_color_manual(values = tribes.colors) +
  theme_classic() +
  
  
  #label="Sp"
  #geom_text(hjust=0, vjust=0)
  
  
  V2R_Doc_sp_wide_ecocat <- left_join(V2R_Doc_sp_wide, ecoCat, by="Sp")
V2R_Doc_sp_wide_ecocat_breeding <- 
  left_join(V2R_Doc_sp_wide_ecocat, breed_infos, by="Sp")

ggplot(V2R_Doc_sp_wide,
       aes(x = PC1,
           y = PC2,
           color=Tribe)) +
  geom_point() +
  stat_ellipse(level = 0.95) +
  scale_color_manual(values = tribes.colors) +
  theme_classic()

ggplot(V2R_Doc_sp_wide_ecocat,
       aes(x = PC1,
           y = PC2,
           color=Food)) +
  geom_point() +
  stat_ellipse(level = 0.95) +
  theme_classic()


ggplot(V2R_Doc_sp_wide_ecocat,
       aes(x = PC1,
           y = PC2,
           color=habitat)) +
  geom_point() +
  stat_ellipse(level = 0.95) +
  theme_classic()



ggplot(V2R_Doc_sp_wide_ecocat_breeding,
       aes(x = PC1,
           y = PC2,
           color=breeding_type)) +
  geom_point() +
  stat_ellipse(level = 0.95) +
  theme_classic()

ggplot(V2R_Doc_sp_wide_ecocat_breeding,
       aes(x = PC1,
           y = PC2,
           color=breeding_mode)) +
  geom_point() +
  stat_ellipse(level = 0.95) +
  theme_classic()



#### PCA T1R nb in  subfamilies   ---------------------------------


T1R_DoC_df_info_mean_wo_total <- 
  T1R_DoC_df_info_mean %>%
  filter(Subfamily != "Total_T1R")

T1R_Doc_sp_wide <- 
  as.data.frame(
    pivot_wider(T1R_DoC_df_info_mean_wo_total, 
                names_from = c(Subfamily), values_from = mean_normalized_nb)
  )


species_in_tree <- b1_tree$tip.label
T1R_Doc_sp_wide <- 
  T1R_Doc_sp_wide %>% 
  filter(Sp != "Neospl") %>% 
  filter(Sp %in% radiation_tree_wo_Neospl$tip.label)

T1R_Doc_sp_wide_woinfo <- 
  T1R_Doc_sp_wide %>%
  dplyr::select(- c(Species, Tribe))

rownames(T1R_Doc_sp_wide_woinfo) <- T1R_Doc_sp_wide_woinfo[,1]
T1R_Doc_sp_wide_woinfo <- T1R_Doc_sp_wide_woinfo[,-1]

tribe <- as.factor(T1R_Doc_sp_wide$Tribe[1:nrow(T1R_Doc_sp_wide)])
tribes.colors <- c(Bathybatini = "#242626" , Benthochromini = "#AE262A" , Boulengerochromini = "#59595C" , Cyphotilapiini = "#FDDF13" , Cyprichromini = "#F04D29" , Ectodini = "#9AB9D9" , Eretmodini = "#682E7A" , Lamprologini = "#C588BB" , Limnochromini = "#535CA9" , Perissodini = "#FBAC43" , Trematocarini = "#959170" , Tropheini = "#86C773" , Haplochromini = "#274e13" , outgroup = "gray", Other="gray")


PCA_T1R <- prcomp(T1R_Doc_sp_wide_woinfo, scale = TRUE)

autoplot(PCA_T1R, data = T1R_Doc_sp_wide, colour = 'Tribe') +
  scale_color_manual(values = tribes.colors) +
  theme_classic()


#### PCA all olfactory receptors nb in  subfamilies   ---------------------------------

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

OLR_DoC_df_info_mean_wo_total <- 
  rbind(OR_DoC_df_info_mean_wo_total, TAAR_DoC_df_info_mean_wo_total,
        V1R_DoC_df_info_mean_wo_total, V2R_DoC_df_info_mean_wo_total)

OLR_Doc_sp_wide <- 
  as.data.frame(
    pivot_wider(OLR_DoC_df_info_mean_wo_total, 
                names_from = c(Subfamily), values_from = mean_normalized_nb)
  )

species_in_tree <- b1_tree$tip.label
OLR_Doc_sp_wide <- 
  OLR_Doc_sp_wide %>% 
  filter(Sp != "Neospl") %>%
  filter(Sp %in% radiation_tree_wo_Neospl$tip.label)
OLR_Doc_sp_wide_woinfo <- 
  OLR_Doc_sp_wide %>%
  dplyr::select(- c(Species, Tribe))

rownames(OLR_Doc_sp_wide_woinfo) <- OLR_Doc_sp_wide_woinfo[,1]
OLR_Doc_sp_wide_woinfo <- OLR_Doc_sp_wide_woinfo[,-1]

tribe <- as.factor(OLR_Doc_sp_wide$Tribe[1:nrow(OLR_Doc_sp_wide)])
tribes.colors <- c(Bathybatini = "#242626" , Benthochromini = "#AE262A" , Boulengerochromini = "#59595C" , Cyphotilapiini = "#FDDF13" , Cyprichromini = "#F04D29" , Ectodini = "#9AB9D9" , Eretmodini = "#682E7A" , Lamprologini = "#C588BB" , Limnochromini = "#535CA9" , Perissodini = "#FBAC43" , Trematocarini = "#959170" , Tropheini = "#86C773" , Haplochromini = "#274e13" , outgroup = "gray", Other="gray")


PCA_OLR <- prcomp(OLR_Doc_sp_wide_woinfo, scale = TRUE)

autoplot(PCA_OLR, data = OLR_Doc_sp_wide, colour = 'Tribe') +
  scale_color_manual(values = tribes.colors) +
  theme_classic()



PC1 <- PCA_OLR$x[,1]
PC2 <- PCA_OLR$x[,2]





pdf(file = "PCA_olfactory_receptors.pdf",width = 8.34,  height = 4.61)


ggplot(OLR_Doc_sp_wide,
       aes(x = PC1,
           y = PC2,
           color=Tribe)) +
  geom_point() +
  stat_ellipse(level = 0.95) +
  scale_color_manual(values = tribes.colors) +
  theme_classic() +
  xlab("PC1 (15.49%)") +
  ylab("PC2 (10.43%)") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none") 


dev.off()





#### PCA all chemoreceptors nb in  subfamilies   ---------------------------------

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

species_in_tree <- b1_tree$tip.label
CR_Doc_sp_wide <- CR_Doc_sp_wide %>% filter(Sp != "Neospl") %>% filter(Sp %in% radiation_tree_wo_Neospl$tip.label)
CR_Doc_sp_wide_woinfo <- 
  CR_Doc_sp_wide %>%
  dplyr::select(- c(Species, Tribe))

rownames(CR_Doc_sp_wide_woinfo) <- CR_Doc_sp_wide_woinfo[,1]
CR_Doc_sp_wide_woinfo <- CR_Doc_sp_wide_woinfo[,-1]

tribe <- as.factor(CR_Doc_sp_wide$Tribe[1:nrow(CR_Doc_sp_wide)])
tribes.colors <- c(Bathybatini = "#242626" , Benthochromini = "#AE262A" , Boulengerochromini = "#59595C" , Cyphotilapiini = "#FDDF13" , Cyprichromini = "#F04D29" , Ectodini = "#9AB9D9" , Eretmodini = "#682E7A" , Lamprologini = "#C588BB" , Limnochromini = "#535CA9" , Perissodini = "#FBAC43" , Trematocarini = "#959170" , Tropheini = "#86C773" , Haplochromini = "#274e13" , outgroup = "gray", Other="gray")


PCA_CR <- prcomp(CR_Doc_sp_wide_woinfo, scale = TRUE)

autoplot(PCA_CR, data = CR_Doc_sp_wide, colour = 'Tribe') +
  scale_color_manual(values = tribes.colors) +
  theme_classic()



PC1 <- PCA_CR$x[,1]
PC2 <- PCA_CR$x[,2]

ggplot(CR_Doc_sp_wide,
       aes(x = PC1,
           y = PC2,
           color=Tribe,
           label=Sp)) +
  geom_point() +
  stat_ellipse(level = 0.95) +
  scale_color_manual(values = tribes.colors) +
  theme_classic() +
  geom_text(hjust=0, vjust=0) +
  xlab("PC1 (15.4%)") +
  ylab("PC2 (10.21%)")


#plotly::ggplotly(PLOTNAME)


pdf(file = "PCA_chemoreceptors.pdf",width = 8.34,  height = 4.61)


ggplot(CR_Doc_sp_wide,
       aes(x = PC1,
           y = PC2,
           color=Tribe)) +
  geom_point() +
  stat_ellipse(level = 0.95) +
  scale_color_manual(values = tribes.colors) +
  theme_classic() +
  xlab("PC1 (15.4%)") +
  ylab("PC2 (10.21%)") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none") 


dev.off()


CR_Doc_sp_wide_ecocat <- left_join(CR_Doc_sp_wide, ecoCat, by="Sp")
CR_Doc_sp_wide_ecocat_breeding <- 
  left_join(CR_Doc_sp_wide_ecocat, breed_infos, by="Sp")

ggplot(CR_Doc_sp_wide,
       aes(x = PC1,
           y = PC2,
           color=Tribe)) +
  geom_point() +
  stat_ellipse(level = 0.95) +
  scale_color_manual(values = tribes.colors) +
  theme_classic()

ggplot(CR_Doc_sp_wide_ecocat,
       aes(x = PC1,
           y = PC2,
           color=Food)) +
  geom_point() +
  theme_classic()


ggplot(CR_Doc_sp_wide_ecocat,
       aes(x = PC1,
           y = PC2,
           color=habitat)) +
  geom_point() +
  stat_ellipse(level = 0.95) +
  theme_classic()



ggplot(CR_Doc_sp_wide_ecocat_breeding,
       aes(x = PC1,
           y = PC2,
           color=breeding_type)) +
  geom_point() +
  stat_ellipse(level = 0.95) +
  theme_classic()

ggplot(CR_Doc_sp_wide_ecocat_breeding,
       aes(x = PC1,
           y = PC2,
           color=breeding_mode)) +
  geom_point() +
  stat_ellipse(level = 0.95) +
  theme_classic()


ggplot(CR_Doc_sp_wide_ecocat_breeding,
       aes(x = PC1,
           y = PC2,
           color=Food)) +
  geom_point() +
  theme_classic()

# Phylogenetic signal of pc1 qnd pc2

PC1_df <- as.data.frame(PCA_CR$x[,1])
PC2_df <- as.data.frame(PCA_CR$x[,2])
PC3_df <- as.data.frame(PCA_CR$x[,3])

PC1_df$Sp <- row.names(PC1_df)                     
PC2_df$Sp <- row.names(PC2_df)         
PC3_df$Sp <- row.names(PC3_df)         

colnames(PC1_df) <- c("PC1", "Sp")
colnames(PC2_df) <- c("PC2", "Sp")
colnames(PC3_df) <- c("PC3", "Sp")

PC1_2 <- left_join(PC1_df, PC2_df, by="Sp")
PC1_2_3 <- left_join(PC1_2, PC3_df, by="Sp")

CR_Doc_sp_wide <- 
  left_join(CR_Doc_sp_wide, 
            PC1_2_3, by="Sp")

CR_Doc_sp_wide_rad <- 
  CR_Doc_sp_wide %>%
  filter(Sp %in% radiation_tree_wo_Neospl$tip.label)
tip_order <- radiation_tree_wo_Neospl$tip.label
CR_Doc_sp_wide_rad_ordered <- 
  enframe(tip_order, name = NULL, value = "Sp") %>% 
  left_join(CR_Doc_sp_wide_rad) %>% 
  column_to_rownames("Sp")

phylosig(tree = radiation_tree_wo_Neospl, x = CR_Doc_sp_wide_rad_ordered$PC1, method = "lambda", test = T)

phylosig(tree = radiation_tree_wo_Neospl, x = CR_Doc_sp_wide_rad_ordered$PC2, method = "lambda", test = T)


## test to put arrows of variance inside the plot 

p1 <- ggplot(CR_Doc_sp_wide,
             aes(x = PC1,
                 y = PC2,
                 color=Tribe)) +
  geom_point() +
  scale_color_manual(values = tribes.colors) +
  theme_classic() +
  xlab("PC1 (15.4%)") +
  ylab("PC2 (10.21%)") +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16),
        legend.position="none") 
  

PCA_CR_loadings <- PCA_CR$rotation
PCA_CR_loadings <- as.data.frame(PCA_CR_loadings)

loadings_scale_factor <- 40
PCA_CR_loadings_scaled <- PCA_CR_loadings * loadings_scale_factor
PCA_CR_loadings_scaled$varname <- rownames(PCA_CR_loadings_scaled)


pdf(file = "PCA_chemoreceptors_loading.pdf",width = 16.34,  height = 10.61)


ggplot(CR_Doc_sp_wide, aes(x = PC1, y = PC2, color = tribe)) +
  geom_point() +
  scale_color_manual(values = tribes.colors) +
  geom_segment(data = PCA_CR_loadings_scaled, aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.2, "cm")), color = "red") +
  geom_text(data = PCA_CR_loadings_scaled, aes(x = PC1, y = PC2, label = varname),
            color = "red", size = 5, vjust = 1.5) +
  theme_classic() +
  xlab(paste0("PC1 (", round(PCA_CR$sdev[1]^2 / sum(PCA_CR$sdev^2) * 100, 2), "%)")) +
  ylab(paste0("PC2 (", round(PCA_CR$sdev[2]^2 / sum(PCA_CR$sdev^2) * 100, 2), "%)")) +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        plot.subtitle = element_text(size = 16),
        legend.position = "none")


dev.off()

# Test the significance of the PCA

test_PCA_signif <- 
  PCAtest(CR_Doc_sp_wide_woinfo, varcorr=FALSE, counter=FALSE, plot=FALSE)

PCAtest(CR_Doc_sp_wide_woinfo, varcorr=FALSE, counter=FALSE, plot=TRUE)


PCA1_signif_columns <- c(2,4,5,8,9,10,11,12,13,14,15,16,18,19,20,21,23,25,27,28,29,30,32,33,34,35,38,39,40,46,47,48,49,50,51,53,62,63,64,67,72,75,76,77,78,79)
PCA2_signif_columns <- c(3,5,6,8,9,10,11,14,15,16,18,20,22,24,25,26,27,28,29,30,31,32,33,35,37,38,39,40,45,46,47,48,49,51,59,63,64,65,67,76,77,78,79)




colnames(CR_Doc_sp_wide_woinfo[,PCA1_signif_columns])
colnames(CR_Doc_sp_wide_woinfo[,PCA2_signif_columns])




#### pGLS PC1 and PC2 ?   ---------------------------------

PC1_df <- as.data.frame(PCA_CR$x[,1])
PC2_df <- as.data.frame(PCA_CR$x[,2])
PC3_df <- as.data.frame(PCA_CR$x[,3])

PC1_df$Sp <- row.names(PC1_df)                     
PC2_df$Sp <- row.names(PC2_df)         
PC3_df$Sp <- row.names(PC3_df)         

colnames(PC1_df) <- c("PC1", "Sp")
colnames(PC2_df) <- c("PC2", "Sp")
colnames(PC3_df) <- c("PC3", "Sp")

PC1_2 <- left_join(PC1_df, PC2_df, by="Sp")
PC1_2_3 <- left_join(PC1_2, PC3_df, by="Sp")
colnames(PC1_2_3) <- c("gPC1", "Sp", "gPC2", "gPC3")


All_Doc_mean_sp_wide_PCA <- 
  left_join(All_Doc_mean_sp_wide, 
            PC1_2_3, by="Sp")


All_Doc_mean_sp_wide_PCA <- 
  All_Doc_mean_sp_wide_PCA %>% 
  mutate(Total_OLR = Total_TAAR+Total_V2R+Total_OR+Total_V1R)



#### Table TE Kimura  ---------------------------------

kimura_dist_TE_df <-
  read.table("Kimura_dist_summary.tsv",
             sep="\t",
             header=TRUE)


kimura_dist_TE_df <- 
  as.data.frame(kimura_dist_TE_df %>%
  rowwise() %>%
  mutate(Significant = case_when(
    (P.value < 0.05) & (Mean_Kimura_distance_Cluster > Mean_Kimura_distance_GenomeWide) ~ "Cluster > GenomeWide",
    (P.value < 0.05) & (Mean_Kimura_distance_Cluster < Mean_Kimura_distance_GenomeWide) ~ "Cluster < GenomeWide",
    P.value > 0.05 ~ "Cluster = GenomeWide"
  )))

kimura_dist_TE_df %>%
  filter(Gene.family %in% c("OR", "TAAR")) %>%
  group_by(Significant) %>%
  summarise(count = n())
