##### Libraries  ---------------------------------

rm(list=ls())

set.seed(2712)


library("caper")
library("ape")
library(dplyr)
library(ggplot2)
library(tidyverse)
library(data.table)
library(phytools)
library(purrr)
library("RERconverge")
library(ggtree)
library(ggtreeExtra)
split_tibble <- function(tibble, column = 'col') {
  tibble %>% split(., .[,column]) %>% lapply(., function(x) x[,setdiff(names(x),column)])
}




#### My functions + Colors palettes  ---------------------------------

args = commandArgs(trailingOnly=TRUE)


PGLS_pvalue <- function(pgls_rslt) {
  sum_cor <- summary(pgls_rslt)
  pvalue =  formatC(sum_cor$coefficients[8], digits = 3)
  if (pvalue == "   0"){ pvalue = 2.2e-16}
  return(pvalue)
}

PGLS_R2 <- function(pgls_rslt) {
  sum_cor <- summary(pgls_rslt)
  r2 = formatC(sum_cor$r.squared, digits = 2)
  return(r2)
}

PGLS_lambda <- function(pgls_rslt) {
  sum_cor <- summary(pgls_rslt)
  PGLS_lambda = sum_cor$param[2]
  return(PGLS_lambda)
}


GLS_function <- function(gls_rslt) {
  GLS_cc <- coef(gls_rslt)
  GLS_f <- function(x) GLS_cc[1] + GLS_cc[2]*x
  return(GLS_f)
}


##### Data load -  GO terms   -------------------

GO_annots.BP <- read.gmt("BP.Onil.gmt")
GO_annots.CC <- read.gmt("CC.Onil.gmt")
GO_annots.MF <- read.gmt("MF.Onil.gmt")
annots.Reactome <- read.gmt("REACTOME.Onil.gmt")
annots.KEGG.Dr <- read.gmt("KEGG.Onil.gmt")

GO_annots_list <- list(GO_annots.BP, GO_annots.CC, GO_annots.MF, annots.Reactome, annots.KEGG.Dr)
names(GO_annots_list) <- 
  c("biological_process", "cellular_component", "molecular_function", "REACTOME", "KEGG")




##### Data load - Species trees   ---------------------------------

species_tree <- 
  read.tree("b1_tree_wo_Neospl.nwk")
cichlid_species <- species_tree$tip.label

##### Data load -  Number of chemoreceptors    ---------------------------------


OR_DoC_df_info <- read.table("OR_DoC_df_info.csv", header=TRUE,  sep=",")
OR_DoC_df_info <- OR_DoC_df_info %>% filter(! ID %in% c("LBE2", "LBE9"))
OR_DoC_df_info_mean <- read.table("OR_DoC_df_info_mean.csv", header=TRUE, sep=",")
OR_DoC_df_info_mean_total <- OR_DoC_df_info_mean %>% filter(Subfamily == "Total_OR") 
OR_number_df <- OR_DoC_df_info_mean_total %>% filter(Sp %in% cichlid_species) %>% dplyr::select(Sp, mean_normalized_nb)
colnames(OR_number_df) <- c("Sp", "OR")


TAAR_DoC_df_info <- read.table("TAAR_DoC_df_info.csv", header=TRUE,  sep=",")
TAAR_DoC_df_info <- TAAR_DoC_df_info %>% filter(! ID %in% c("LBE2", "LBE9"))
TAAR_DoC_df_info_mean <- read.table("TAAR_DoC_df_info_mean.csv", header=TRUE, sep=",")
TAAR_DoC_df_info_mean_total <- TAAR_DoC_df_info_mean %>% filter(Subfamily == "Total_TAAR") 
TAAR_number_df <- TAAR_DoC_df_info_mean_total %>% filter(Sp %in% cichlid_species) %>% dplyr::select(Sp, mean_normalized_nb)
colnames(TAAR_number_df) <- c("Sp", "TAAR")

V2R_DoC_df_info <- read.table("V2R_DoC_df_info.csv", header=TRUE,  sep=",")
V2R_DoC_df_info <- V2R_DoC_df_info %>% filter(! ID %in% c("LBE2", "LBE9"))
V2R_DoC_df_info_mean <- read.table("V2R_DoC_df_info_mean.csv", header=TRUE, sep=",")
V2R_DoC_df_info_mean_total <- V2R_DoC_df_info_mean %>% filter(Subfamily == "Total_V2R") 
V2R_number_df <- V2R_DoC_df_info_mean_total %>% dplyr::filter(Sp %in% cichlid_species) %>% dplyr::select(Sp, mean_normalized_nb)
colnames(V2R_number_df) <- c("Sp", "V2R")



V1R_DoC_df_info <- read.table("V1R_DoC_df_info.csv", header=TRUE,  sep=",")
V1R_DoC_df_info <- V1R_DoC_df_info %>% filter(! ID %in% c("LBE2", "LBE9"))
V1R_DoC_df_info_mean <- read.table("V1R_DoC_df_info_mean.csv", header=TRUE, sep=",")
V1R_DoC_df_info_mean_total <- V1R_DoC_df_info_mean %>% filter(Subfamily == "Total_V1R") 
V1R_number_df <- V1R_DoC_df_info_mean_total %>% dplyr::filter(Sp %in% cichlid_species) %>% dplyr::select(Sp, mean_normalized_nb)
colnames(V1R_number_df) <- c("Sp", "V1R")




OLR_number_df <- left_join(OR_number_df, TAAR_number_df, by="Sp")
OLR_number_df <- left_join(OLR_number_df, V2R_number_df, by="Sp")
OLR_number_df <- left_join(OLR_number_df, V1R_number_df, by="Sp")

OLR_number_df <- OLR_number_df %>% mutate(OLR = OR + TAAR + V2R + V1R) %>% dplyr::select(Sp, OLR)


##### Data load -  Number of lamellae    ---------------------------------

Lamellae_number_df <- 
  read.table("Cichlid_rosette_table.tsv",
            sep="\t",
            header=FALSE)
colnames(Lamellae_number_df) <- c("Sp", "mean_lamellae", "mean_SL", "mean_TL", "mean_weight")
Lamellae_number_df <- Lamellae_number_df %>% dplyr::select(Sp, mean_lamellae) %>% filter(Sp %in% cichlid_species)


##### Extract the lamellae residual against size - Lm    ---------------------------------

Lamellae_df <- 
  read.table("Cichlid_rosette_table.tsv",
             sep="\t",
             header=FALSE)
colnames(Lamellae_df) <- c("Sp", "mean_lamellae", "mean_SL", "mean_TL", "mean_weight")
Lamellae_df <- Lamellae_df %>% filter(Sp %in% species_tree$tip.label)

fit_lamellae_SL <- lm(data = Lamellae_df, 
             formula = mean_lamellae ~ mean_SL)


summary(fit_lamellae_SL)
sum_fit_phy <- summary(fit_lamellae_SL)

lm_residuals <- as.data.frame(sum_fit_phy$residuals)
colnames(lm_residuals) <- "lm_residuals"
lm_residuals$Sp <- Lamellae_df$Sp

lm_residuals_df <- lm_residuals %>% dplyr::select(Sp, lm_residuals)


##### Extract the lamellae residual against size - PGLS    ---------------------------------

Lamellae_df <- 
  read.table("Cichlid_rosette_table.tsv",
             sep="\t",
             header=FALSE)
colnames(Lamellae_df) <- c("Sp", "mean_lamellae", "mean_SL", "mean_TL", "mean_weight")
Lamellae_df <- Lamellae_df %>% filter(Sp %in% species_tree$tip.label)
caper_data_cichlids_SL_lamellae <- 
  comparative.data(phy = species_tree, 
                   data = Lamellae_df,
                   names.col = Sp, vcv = TRUE,
                   na.omit = FALSE, warn.dropped = TRUE)



fit_phylo_lamellae_SL <- pgls(mean_lamellae ~ mean_SL,
                              data = caper_data_cichlids_SL_lamellae, 
                              lambda = "ML")
summary(fit_phylo_lamellae_SL)
sum_fit_phy <- summary(fit_phylo_lamellae_SL)

pgls_residuals <- as.data.frame(sum_fit_phy$residuals)
colnames(pgls_residuals) <- "pGLS_residuals"
pgls_residuals$Sp <- Lamellae_df$Sp

pgls_residuals_df <- pgls_residuals %>% dplyr::select(Sp, pGLS_residuals)



##### Read trees and form a MasterTree   ---------------------------------

GW_trees.all <- 
  RERconverge::readTrees(
    "AllTrees.Phangorn.nwk",
    reestimateBranches = FALSE,
    minSpecs=80)



##### Launch RERconverge   ---------------------------------
 
olfactory_receptors.RERw.all <-
  getAllResiduals(GW_trees.all,
                  transform = "sqrt",
                  useSpecies = intersect(OLR_number_df %>% pull(Sp) , GW_trees.all$masterTree$tip.label),
                  weighted = T,
                  scale = T,
                  min.sp = 80)


lamellae.RERw.all <-
  getAllResiduals(GW_trees.all,
                  transform = "sqrt",
                  useSpecies = intersect(Lamellae_number_df %>% pull(Sp) , GW_trees.all$masterTree$tip.label),
                  weighted = T,
                  scale = T,
                  min.sp = 50)




## Convert the trait vector to a RER matrix path

OR.trait_vector <- setNames(OR_number_df[[2]], OR_number_df[[1]])
OR.charpaths <- char2Paths(OR.trait_vector, GW_trees.all)

TAAR.trait_vector <- setNames(TAAR_number_df[[2]], TAAR_number_df[[1]])
TAAR.charpaths <- char2Paths(TAAR.trait_vector, GW_trees.all)

V2R.trait_vector <- setNames(V2R_number_df[[2]], V2R_number_df[[1]])
V2R.charpaths <- char2Paths(V2R.trait_vector, GW_trees.all)

OLR.trait_vector <- setNames(OLR_number_df[[2]], OLR_number_df[[1]])
OLR.charpaths <- char2Paths(OLR.trait_vector, GW_trees.all)

lamellae.trait_vector <- setNames(Lamellae_number_df[[2]], Lamellae_number_df[[1]])
lamellae.charpaths <- char2Paths(lamellae.trait_vector, GW_trees.all)

pGLS_residual_lamellae.trait_vector <- setNames(pgls_residuals_df[[2]], pgls_residuals_df[[1]])
pGLS_residual_lamellae.charpaths <- char2Paths(pGLS_residual_lamellae.trait_vector, GW_trees.all)


lm_residual_lamellae.trait_vector <- setNames(lm_residuals_df[[2]], lm_residuals_df[[1]])
lm_residual_lamellae.charpaths <- char2Paths(lm_residual_lamellae.trait_vector, GW_trees.all)


## Find correlations between genes rates of evolution and the rate of change of the phenotype

OR.all.res <- 
  correlateWithContinuousPhenotype(
    olfactory_receptors.RERw.all, OR.charpaths,
    min.sp = 80, winsorizeRER = 5, winsorizetrait = 5)
OR.all.res$gene_name <- row.names(OR.all.res)
OR.all.res <- OR.all.res %>% mutate(trait = "OR")
OR.all.res$stat <- -log10(OR.all.res$P) * sign(OR.all.res$Rho)


TAAR.all.res <- 
  correlateWithContinuousPhenotype(
    olfactory_receptors.RERw.all, TAAR.charpaths,
    min.sp = 80, winsorizeRER = 5, winsorizetrait = 5)
TAAR.all.res$gene_name <- row.names(TAAR.all.res)
TAAR.all.res <- TAAR.all.res %>% mutate(trait = "TAAR")
TAAR.all.res$stat <- -log10(TAAR.all.res$P) * sign(TAAR.all.res$Rho)

V2R.all.res <- 
  correlateWithContinuousPhenotype(
    olfactory_receptors.RERw.all, V2R.charpaths,
    min.sp = 80, winsorizeRER = 5, winsorizetrait = 5)
V2R.all.res$gene_name <- row.names(V2R.all.res)
V2R.all.res <- V2R.all.res %>% mutate(trait = "V2R")
V2R.all.res$stat <- -log10(V2R.all.res$P) * sign(V2R.all.res$Rho)


OLR.all.res <- 
  correlateWithContinuousPhenotype(
    olfactory_receptors.RERw.all, OLR.charpaths,
    min.sp = 80, winsorizeRER = 5, winsorizetrait = 5)
OLR.all.res$gene_name <- row.names(OLR.all.res)
OLR.all.res <- OLR.all.res %>% mutate(trait = "OLR")
OLR.all.res$stat <- -log10(OLR.all.res$P) * sign(OLR.all.res$Rho)

lamellae.all.res <- 
  correlateWithContinuousPhenotype(
    lamellae.RERw.all, lamellae.charpaths,
    min.sp = 50, winsorizeRER = 0, winsorizetrait = 0)
lamellae.all.res$gene_name <- row.names(lamellae.all.res)
lamellae.all.res <- lamellae.all.res %>% mutate(trait = "Lamellae")
lamellae.all.res$stat <- -log10(lamellae.all.res$P) * sign(lamellae.all.res$Rho)


pGLS_residual_lamellae.all.res <- 
  correlateWithContinuousPhenotype(
    lamellae.RERw.all, pGLS_residual_lamellae.charpaths,
    min.sp = 50, winsorizeRER = 0, winsorizetrait = 0)
pGLS_residual_lamellae.all.res$gene_name <- row.names(pGLS_residual_lamellae.all.res)
pGLS_residual_lamellae.all.res <- pGLS_residual_lamellae.all.res %>% mutate(trait = "Lamellae")
pGLS_residual_lamellae.all.res$stat <- -log10(pGLS_residual_lamellae.all.res$P) * sign(pGLS_residual_lamellae.all.res$Rho)


lm_residual_lamellae.all.res <- 
  correlateWithContinuousPhenotype(
    lamellae.RERw.all, lm_residual_lamellae.charpaths,
    min.sp = 50, winsorizeRER = 0, winsorizetrait = 0)
lm_residual_lamellae.all.res$gene_name <- row.names(lm_residual_lamellae.all.res)
lm_residual_lamellae.all.res <- lm_residual_lamellae.all.res %>% mutate(trait = "Lamellae")
lm_residual_lamellae.all.res$stat <- -log10(lm_residual_lamellae.all.res$P) * sign(lm_residual_lamellae.all.res$Rho)

print("RERconverge in numeric mode ran successfully")

#Save results to dataframe

print("Launching numerical permulation analysis")

#Generate a rooted master tree
rooted_master_tree <- ape::root(GW_trees.all$masterTree, "Boumic", resolve.root= TRUE)

#Perform enrichment analysis on the raw results

OR.my_stat <- sign(OR.all.res$Rho) * (-log10(OR.all.res$P))
names(OR.my_stat) <- rownames(OR.all.res)
genenames <- names(OR.my_stat)
OR.my_stat <- OR.my_stat[!is.na(OR.my_stat)]
OR.enrichment <- fastwilcoxGMTall(OR.my_stat, GO_annots_list, outputGeneVals=T, num.g=10)

TAAR.my_stat <- sign(TAAR.all.res$Rho) * (-log10(TAAR.all.res$P))
names(TAAR.my_stat) <- rownames(TAAR.all.res)
genenames <- names(TAAR.my_stat)
TAAR.my_stat <- TAAR.my_stat[!is.na(TAAR.my_stat)]
TAAR.enrichment <- fastwilcoxGMTall(TAAR.my_stat, GO_annots_list, outputGeneVals=T, num.g=10)

V2R.my_stat <- sign(V2R.all.res$Rho) * (-log10(V2R.all.res$P))
names(V2R.my_stat) <- rownames(V2R.all.res)
genenames <- names(V2R.my_stat)
V2R.my_stat <- V2R.my_stat[!is.na(V2R.my_stat)]
V2R.enrichment <- fastwilcoxGMTall(V2R.my_stat, GO_annots_list, outputGeneVals=T, num.g=10)

OLR.my_stat <- sign(OLR.all.res$Rho) * (-log10(OLR.all.res$P))
names(OLR.my_stat) <- rownames(OLR.all.res)
genenames <- names(OLR.my_stat)
OLR.my_stat <- OLR.my_stat[!is.na(OLR.my_stat)]
OLR.enrichment <- fastwilcoxGMTall(OLR.my_stat, GO_annots_list, outputGeneVals=T, num.g=10)

lamellae.my_stat <- sign(lamellae.all.res$Rho) * (-log10(lamellae.all.res$P))
names(lamellae.my_stat) <- rownames(lamellae.all.res)
genenames <- names(lamellae.my_stat)
lamellae.my_stat <- lamellae.my_stat[!is.na(lamellae.my_stat)]
lamellae.enrichment <- fastwilcoxGMTall(lamellae.my_stat, GO_annots_list, outputGeneVals=T, num.g=10)

pGLS_residual_lamellae.my_stat <- sign(pGLS_residual_lamellae.all.res$Rho) * (-log10(pGLS_residual_lamellae.all.res$P))
names(pGLS_residual_lamellae.my_stat) <- rownames(pGLS_residual_lamellae.all.res)
genenames <- names(pGLS_residual_lamellae.my_stat)
pGLS_residual_lamellae.my_stat <- pGLS_residual_lamellae.my_stat[!is.na(pGLS_residual_lamellae.my_stat)]
pGLS_residual_lamellae.enrichment <- fastwilcoxGMTall(pGLS_residual_lamellae.my_stat, GO_annots_list, outputGeneVals=T, num.g=10)


lm_residual_lamellae.my_stat <- sign(lm_residual_lamellae.all.res$Rho) * (-log10(lm_residual_lamellae.all.res$P))
names(lm_residual_lamellae.my_stat) <- rownames(lm_residual_lamellae.all.res)
genenames <- names(lm_residual_lamellae.my_stat)
lm_residual_lamellae.my_stat <- lm_residual_lamellae.my_stat[!is.na(lm_residual_lamellae.my_stat)]
lm_residual_lamellae.enrichment <- fastwilcoxGMTall(lm_residual_lamellae.my_stat, GO_annots_list, outputGeneVals=T, num.g=10)



#Lets launch 1,000 permulations of our numerical phenotype  + permulation enrichment analysis
my_perms_continuous_OR <- 
  getPermsContinuous(
    numperms = 1000,
    traitvec = OR.trait_vector,
    RERmat = olfactory_receptors.RERw.all,
    annotlist= GO_annots_list,
    trees = GW_trees.all,
    mastertree = rooted_master_tree,
    type = "simperm",
    winT = 5,
    winR = 5,
    method = "p",
    calculateenrich = FALSE)

my_perms_continuous_TAAR <- 
  getPermsContinuous(
    numperms = 1000,
    traitvec = TAAR.trait_vector,
    RERmat = olfactory_receptors.RERw.all,
    annotlist= GO_annots_list,
    trees = GW_trees.all,
    mastertree = rooted_master_tree,
    type = "simperm",
    winT = 5,
    winR = 5,
    method = "p",
    calculateenrich = FALSE)

my_perms_continuous_V2R <- 
  getPermsContinuous(
    numperms = 1000,
    traitvec = V2R.trait_vector,
    RERmat = olfactory_receptors.RERw.all,
    annotlist= GO_annots_list,
    trees = GW_trees.all,
    mastertree = rooted_master_tree,
    type = "simperm",
    winT = 5,
    winR = 5,
    method = "p",
    calculateenrich = FALSE)


my_perms_continuous_OLR <- 
  getPermsContinuous(
    numperms = 1000,
    traitvec = OLR.trait_vector,
    RERmat = olfactory_receptors.RERw.all,
    annotlist= GO_annots_list,
    trees = GW_trees.all,
    mastertree = rooted_master_tree,
    type = "simperm",
    winT = 5,
    winR = 5,
    method = "p",
    calculateenrich = FALSE)


my_perms_continuous_lamellae <- 
  getPermsContinuous(
    numperms = 1000,
    traitvec = lamellae.trait_vector,
    RERmat = lamellae.RERw.all,
    annotlist= GO_annots_list,
    trees = GW_trees.all,
    mastertree = rooted_master_tree,
    type = "simperm",
    winT = 5,
    winR = 5,
    method = "p",
    calculateenrich = FALSE)


my_perms_continuous_pGLS_residual_lamellae <- 
  getPermsContinuous(
    numperms = 1000,
    traitvec = pGLS_residual_lamellae.trait_vector,
    RERmat = lamellae.RERw.all,
    annotlist= GO_annots_list,
    trees = GW_trees.all,
    mastertree = rooted_master_tree,
    type = "simperm",
    winT = 5,
    winR = 5,
    method = "p",
    calculateenrich = FALSE)


my_perms_continuous_lm_residual_lamellae <- 
  getPermsContinuous(
    numperms = 1000,
    traitvec = lm_residual_lamellae.trait_vector,
    RERmat = lamellae.RERw.all,
    annotlist= GO_annots_list,
    trees = GW_trees.all,
    mastertree = rooted_master_tree,
    type = "simperm",
    winT = 5,
    winR = 5,
    method = "p",
    calculateenrich = FALSE)



OR.permswithenrich <- getEnrichPerms(my_perms_continuous_OR, OR.enrichment, GO_annots_list)
TAAR.permswithenrich <- getEnrichPerms(my_perms_continuous_TAAR, TAAR.enrichment, GO_annots_list)
V2R.permswithenrich <- getEnrichPerms(my_perms_continuous_V2R, V2R.enrichment, GO_annots_list)
OLR.permswithenrich <- getEnrichPerms(my_perms_continuous_OLR, OLR.enrichment, GO_annots_list)
lamellae.permswithenrich <- getEnrichPerms(my_perms_continuous_lamellae, lamellae.enrichment, GO_annots_list)
pGLS_residual_lamellae.permswithenrich <- getEnrichPerms(my_perms_continuous_pGLS_residual_lamellae, pGLS_residual_lamellae.enrichment, GO_annots_list)
lm_residual_lamellae.permswithenrich <- getEnrichPerms(my_perms_continuous_lm_residual_lamellae, lm_residual_lamellae.enrichment, GO_annots_list)

#Compute permulations p-values
OR.corpermpvals <- permpvalcor(OR.all.res, my_perms_continuous_OR) 
TAAR.corpermpvals <- permpvalcor(TAAR.all.res, my_perms_continuous_TAAR) 
V2R.corpermpvals <- permpvalcor(V2R.all.res, my_perms_continuous_V2R) 
OLR.corpermpvals <- permpvalcor(OLR.all.res, my_perms_continuous_OLR) 
lamellae.corpermpvals <- permpvalcor(lamellae.all.res, my_perms_continuous_lamellae) 
pGLS_residual_lamellae.corpermpvals <- permpvalcor(pGLS_residual_lamellae.all.res, my_perms_continuous_pGLS_residual_lamellae) 
lm_residual_lamellae.corpermpvals <- permpvalcor(lm_residual_lamellae.all.res, my_perms_continuous_lm_residual_lamellae) 


#Compute permulations enrichment p-values
OR.enrichpermpvals <- permpvalenrich(OR.enrichment, OR.permswithenrich)
TAAR.enrichpermpvals <- permpvalenrich(TAAR.enrichment, TAAR.permswithenrich)
V2R.enrichpermpvals <- permpvalenrich(V2R.enrichment, V2R.permswithenrich)
OLR.enrichpermpvals <- permpvalenrich(OLR.enrichment, OLR.permswithenrich)
lamellae.enrichpermpvals <- permpvalenrich(lamellae.enrichment, lamellae.permswithenrich)
pGLS_residual_lamellae.enrichpermpvals <- permpvalenrich(pGLS_residual_lamellae.enrichment, pGLS_residual_lamellae.permswithenrich)
lm_residual_lamellae.enrichpermpvals <- permpvalenrich(lm_residual_lamellae.enrichment, lm_residual_lamellae.permswithenrich)


#Add the enrichment permulations p-values to the real enrichment list

count=1
while(count<=length(OR.enrichment)){
  OR.enrichment[[count]]$permpval <- OR.enrichpermpvals[[count]][match(rownames(OR.enrichment[[count]]),names(OR.enrichpermpvals[[count]]))]
  OR.enrichment[[count]]$permpvaladj <- p.adjust(OR.enrichment[[count]]$permpval, method="BH")
  count=count+1 
}

count=1
while(count<=length(TAAR.enrichment)){
  TAAR.enrichment[[count]]$permpval <- TAAR.enrichpermpvals[[count]][match(rownames(TAAR.enrichment[[count]]),names(TAAR.enrichpermpvals[[count]]))]
  TAAR.enrichment[[count]]$permpvaladj <- p.adjust(TAAR.enrichment[[count]]$permpval, method="BH")
  count=count+1 
}

count=1
while(count<=length(V2R.enrichment)){
  V2R.enrichment[[count]]$permpval <- V2R.enrichpermpvals[[count]][match(rownames(V2R.enrichment[[count]]),names(V2R.enrichpermpvals[[count]]))]
  V2R.enrichment[[count]]$permpvaladj <- p.adjust(V2R.enrichment[[count]]$permpval, method="BH")
  count=count+1 
}

count=1
while(count<=length(OLR.enrichment)){
  OLR.enrichment[[count]]$permpval <- OLR.enrichpermpvals[[count]][match(rownames(OLR.enrichment[[count]]),names(OLR.enrichpermpvals[[count]]))]
  OLR.enrichment[[count]]$permpvaladj <- p.adjust(OLR.enrichment[[count]]$permpval, method="BH")
  count=count+1 
}

count=1
while(count<=length(lamellae.enrichment)){
  lamellae.enrichment[[count]]$permpval <- lamellae.enrichpermpvals[[count]][match(rownames(lamellae.enrichment[[count]]),names(lamellae.enrichpermpvals[[count]]))]
  lamellae.enrichment[[count]]$permpvaladj <- p.adjust(lamellae.enrichment[[count]]$permpval, method="BH")
  count=count+1 
}


count=1
while(count<=length(pGLS_residual_lamellae.enrichment)){
  pGLS_residual_lamellae.enrichment[[count]]$permpval <- pGLS_residual_lamellae.enrichpermpvals[[count]][match(rownames(pGLS_residual_lamellae.enrichment[[count]]),names(pGLS_residual_lamellae.enrichpermpvals[[count]]))]
  pGLS_residual_lamellae.enrichment[[count]]$permpvaladj <- p.adjust(pGLS_residual_lamellae.enrichment[[count]]$permpval, method="BH")
  count=count+1 
}

count=1
while(count<=length(lm_residual_lamellae.enrichment)){
  lm_residual_lamellae.enrichment[[count]]$permpval <- lm_residual_lamellae.enrichpermpvals[[count]][match(rownames(lm_residual_lamellae.enrichment[[count]]),names(lm_residual_lamellae.enrichpermpvals[[count]]))]
  lm_residual_lamellae.enrichment[[count]]$permpvaladj <- p.adjust(lm_residual_lamellae.enrichment[[count]]$permpval, method="BH")
  count=count+1 
}




#Add the permulations p-values and statistics (just -log10(p-values) * sign of Rho) to the results 

OR.all.res.perm <- 
  left_join(
    rownames_to_column(OR.all.res), rownames_to_column(OR.corpermpvals), by="rowname") %>% #for permulations_results_df if not enrichment analysis done
  column_to_rownames(var = "rowname")

TAAR.all.res.perm <- 
  left_join(
    rownames_to_column(TAAR.all.res), rownames_to_column(TAAR.corpermpvals), by="rowname") %>%
  column_to_rownames(var = "rowname")


V2R.all.res.perm <- 
  left_join(
    rownames_to_column(V2R.all.res), rownames_to_column(V2R.corpermpvals), by="rowname") %>%
  column_to_rownames(var = "rowname")


OLR.all.res.perm <- 
  left_join(
    rownames_to_column(OLR.all.res), rownames_to_column(OLR.corpermpvals), by="rowname") %>%
  column_to_rownames(var = "rowname")


lamellae.all.res.perm <- 
  left_join(
    rownames_to_column(lamellae.all.res), rownames_to_column(lamellae.corpermpvals), by="rowname") %>%
  column_to_rownames(var = "rowname")

pGLS_residual_lamellae.all.res.perm <- 
  left_join(
    rownames_to_column(pGLS_residual_lamellae.all.res), rownames_to_column(pGLS_residual_lamellae.corpermpvals), by="rowname") %>%
  column_to_rownames(var = "rowname")


lm_residual_lamellae.all.res.perm <- 
  left_join(
    rownames_to_column(lm_residual_lamellae.all.res), rownames_to_column(lm_residual_lamellae.corpermpvals), by="rowname") %>%
  column_to_rownames(var = "rowname")


#Compute the adjusted p-value :) 
OR.all.res.perm$permpvaladj <- p.adjust(OR.all.res.perm$permpval, method="BH")
TAAR.all.res.perm$permpvaladj <- p.adjust(TAAR.all.res.perm$permpval, method="BH")
V2R.all.res.perm$permpvaladj <- p.adjust(V2R.all.res.perm$permpval, method="BH")
OLR.all.res.perm$permpvaladj <- p.adjust(OLR.all.res.perm$permpval, method="BH")
lamellae.all.res.perm$permpvaladj <- p.adjust(lamellae.all.res.perm$permpval, method="BH")
pGLS_residual_lamellae.all.res.perm$permpvaladj <- p.adjust(pGLS_residual_lamellae.all.res.perm$permpval, method="BH")
lm_residual_lamellae.all.res.perm$permpvaladj <- p.adjust(lm_residual_lamellae.all.res.perm$permpval, method="BH")



