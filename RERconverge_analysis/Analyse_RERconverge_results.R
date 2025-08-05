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




##### Data load -  gene names   -------------------

gene_names_df <- 
  read.table("Names.Onil.tsv",
             sep="\t",
             header=FALSE)

colnames(gene_names_df) <- c("gene_name", "dr_gene_name","gene_symbol", "entrez_gene", 
                             "human_gene_symbol")

gene_names_df <- 
  gene_names_df %>%
  distinct(gene_name, .keep_all = TRUE)


#Keep only one row per gene



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


##### Extract the lamellae residual against size - PGLS    
##### Load RERconverge associations results    ---------------------------------

#enrichment with negative stat : decelerated 
#enrichment with positive stat : accelerated

#OR
RER_OR <- 
  read.table("RERconverge_Results/OR.all.res.perm",
             sep=",",
             header=TRUE)
#TAAR
RER_TAAR <- 
  read.table("RERconverge_Results/TAAR.all.res.perm",
             sep=",",
             header=TRUE)
#V2R
RER_V2R <- 
  read.table("RERconverge_Results/V2R.all.res.perm",
             sep=",",
             header=TRUE)
#OLR
RER_OLR <- 
  read.table("RERconverge_Results/OLR.all.res.perm",
             sep=",",
             header=TRUE)

#Lamellae
RER_lamellae <- 
  read.table("RERconverge_Results/lamellae.all.res.perm",
             sep=",",
             header=TRUE)

#Res lamellae
RER_res_lamellae <- 
  read.table("RERconverge_Results/lm_residual_lamellae.all.res.perm",
             sep=",",
             header=TRUE)


#Process tables to remove NA 
RER_OR <- RER_OR %>% filter(! is.na(Rho))
RER_TAAR <- RER_TAAR %>% filter(! is.na(Rho))
RER_V2R <- RER_V2R %>% filter(! is.na(Rho))
RER_OLR <- RER_OLR %>% filter(! is.na(Rho))
RER_lamellae <- RER_lamellae %>% filter(! is.na(Rho))
RER_res_lamellae <- RER_res_lamellae %>% filter(! is.na(Rho))


#Combine tables with gene names

RER_OR <- left_join(RER_OR, gene_names_df, by="gene_name")
RER_TAAR <- left_join(RER_TAAR, gene_names_df, by="gene_name")
RER_V2R <- left_join(RER_V2R, gene_names_df, by="gene_name")
RER_OLR <- left_join(RER_OLR, gene_names_df, by="gene_name")
RER_lamellae <- left_join(RER_lamellae, gene_names_df, by="gene_name")
RER_res_lamellae <- left_join(RER_res_lamellae, gene_names_df, by="gene_name")

#Add the significance
RER_OR <- 
  RER_OR %>%
  mutate(category = case_when(
    Rho > 0 & p.adj < 0.05 ~ "Accelerated",
    Rho < 0 & p.adj < 0.05 ~ "Decelerated",
    p.adj >= 0.05 ~ "Non_signif",
  )) %>%
  mutate(category_perm = case_when(
    Rho > 0 & permpval < 0.05 ~ "Accelerated",
    Rho < 0 & permpval < 0.05 ~ "Decelerated",
    permpval >= 0.05 ~ "Non_signif",
  )) %>%
  mutate(both_category = case_when(
    category_perm == "Accelerated" & category == "Accelerated" ~ "Accelerated",
    category_perm == "Accelerated" & category == "Non_signif" ~ "Non_signif",
    category_perm == "Non_signif" & category == "Accelerated" ~ "Non_signif",
    category_perm == "Decelerated" & category == "Decelerated" ~ "Decelerated",
    category_perm == "Non_signif" & category == "Decelerated" ~ "Non_signif",
    category_perm == "Decelerated" & category == "Non_signif" ~ "Non_signif",
    category_perm == "Non_signif" & category == "Non_signif" ~ "Non_signif"
  ))



RER_TAAR <- 
  RER_TAAR %>%
  mutate(category = case_when(
    Rho > 0 & p.adj < 0.05 ~ "Accelerated",
    Rho < 0 & p.adj < 0.05 ~ "Decelerated",
    p.adj >= 0.05 ~ "Non_signif",
  )) %>%
  mutate(category_perm = case_when(
    Rho > 0 & permpval < 0.05 ~ "Accelerated",
    Rho < 0 & permpval < 0.05 ~ "Decelerated",
    permpval >= 0.05 ~ "Non_signif",
  )) %>%
  mutate(both_category = case_when(
    category_perm == "Accelerated" & category == "Accelerated" ~ "Accelerated",
    category_perm == "Accelerated" & category == "Non_signif" ~ "Non_signif",
    category_perm == "Non_signif" & category == "Accelerated" ~ "Non_signif",
    category_perm == "Decelerated" & category == "Decelerated" ~ "Decelerated",
    category_perm == "Non_signif" & category == "Decelerated" ~ "Non_signif",
    category_perm == "Decelerated" & category == "Non_signif" ~ "Non_signif",
    category_perm == "Non_signif" & category == "Non_signif" ~ "Non_signif"
  ))



RER_V2R <- 
  RER_V2R %>%
  mutate(category = case_when(
    Rho > 0 & p.adj < 0.05 ~ "Accelerated",
    Rho < 0 & p.adj < 0.05 ~ "Decelerated",
    p.adj >= 0.05 ~ "Non_signif",
  )) %>%
  mutate(category_perm = case_when(
    Rho > 0 & permpval < 0.05 ~ "Accelerated",
    Rho < 0 & permpval < 0.05 ~ "Decelerated",
    permpval >= 0.05 ~ "Non_signif",
  )) %>%
  mutate(both_category = case_when(
    category_perm == "Accelerated" & category == "Accelerated" ~ "Accelerated",
    category_perm == "Accelerated" & category == "Non_signif" ~ "Non_signif",
    category_perm == "Non_signif" & category == "Accelerated" ~ "Non_signif",
    category_perm == "Decelerated" & category == "Decelerated" ~ "Decelerated",
    category_perm == "Non_signif" & category == "Decelerated" ~ "Non_signif",
    category_perm == "Decelerated" & category == "Non_signif" ~ "Non_signif",
    category_perm == "Non_signif" & category == "Non_signif" ~ "Non_signif"
  ))



RER_OLR <- 
  RER_OLR %>%
  mutate(category = case_when(
    Rho > 0 & p.adj < 0.05 ~ "Accelerated",
    Rho < 0 & p.adj < 0.05 ~ "Decelerated",
    p.adj >= 0.05 ~ "Non_signif",
  )) %>%
  mutate(category_perm = case_when(
    Rho > 0 & permpval < 0.05 ~ "Accelerated",
    Rho < 0 & permpval < 0.05 ~ "Decelerated",
    permpval >= 0.05 ~ "Non_signif",
  )) %>%
  mutate(both_category = case_when(
    category_perm == "Accelerated" & category == "Accelerated" ~ "Accelerated",
    category_perm == "Accelerated" & category == "Non_signif" ~ "Non_signif",
    category_perm == "Non_signif" & category == "Accelerated" ~ "Non_signif",
    category_perm == "Decelerated" & category == "Decelerated" ~ "Decelerated",
    category_perm == "Non_signif" & category == "Decelerated" ~ "Non_signif",
    category_perm == "Decelerated" & category == "Non_signif" ~ "Non_signif",
    category_perm == "Non_signif" & category == "Non_signif" ~ "Non_signif"
  ))



RER_lamellae <- 
  RER_lamellae %>%
  mutate(category = case_when(
    Rho > 0 & p.adj < 0.05 ~ "Accelerated",
    Rho < 0 & p.adj < 0.05 ~ "Decelerated",
    p.adj >= 0.05 ~ "Non_signif",
  )) %>%
  mutate(category_perm = case_when(
    Rho > 0 & permpval < 0.05 ~ "Accelerated",
    Rho < 0 & permpval < 0.05 ~ "Decelerated",
    permpval >= 0.05 ~ "Non_signif",
  )) %>%
  mutate(both_category = case_when(
    category_perm == "Accelerated" & category == "Accelerated" ~ "Accelerated",
    category_perm == "Accelerated" & category == "Non_signif" ~ "Non_signif",
    category_perm == "Non_signif" & category == "Accelerated" ~ "Non_signif",
    category_perm == "Decelerated" & category == "Decelerated" ~ "Decelerated",
    category_perm == "Non_signif" & category == "Decelerated" ~ "Non_signif",
    category_perm == "Decelerated" & category == "Non_signif" ~ "Non_signif",
    category_perm == "Non_signif" & category == "Non_signif" ~ "Non_signif"
  ))




RER_res_lamellae <- 
  RER_res_lamellae %>%
  mutate(category = case_when(
    Rho > 0 & p.adj < 0.05 ~ "Accelerated",
    Rho < 0 & p.adj < 0.05 ~ "Decelerated",
    p.adj >= 0.05 ~ "Non_signif",
  )) %>%
  mutate(category_perm = case_when(
    Rho > 0 & permpval < 0.05 ~ "Accelerated",
    Rho < 0 & permpval < 0.05 ~ "Decelerated",
    permpval >= 0.05 ~ "Non_signif",
  )) %>%
  mutate(both_category = case_when(
    category_perm == "Accelerated" & category == "Accelerated" ~ "Accelerated",
    category_perm == "Accelerated" & category == "Non_signif" ~ "Non_signif",
    category_perm == "Non_signif" & category == "Accelerated" ~ "Non_signif",
    category_perm == "Decelerated" & category == "Decelerated" ~ "Decelerated",
    category_perm == "Non_signif" & category == "Decelerated" ~ "Non_signif",
    category_perm == "Decelerated" & category == "Non_signif" ~ "Non_signif",
    category_perm == "Non_signif" & category == "Non_signif" ~ "Non_signif"
  ))



##### Load RERconverge enrichment results    ---------------------------------

#OR
enrich_OR <- 
  readRDS("RERconverge_Results/OR.enrichment.RDS")

#TAAR
enrich_TAAR <- 
  readRDS("RERconverge_Results/TAAR.enrichment.RDS")

#V2R
enrich_V2R <- 
  readRDS("RERconverge_Results/V2R.enrichment.RDS")

#OLR
enrich_OLR <- 
  readRDS("RERconverge_Results/OLR.enrichment.RDS")

#Lamellae
enrich_lamellae <- 
  readRDS("RERconverge_Results/lamellae.enrichment.RDS")

#Residual Lamellae
enrich_res_lamellae <- 
  readRDS("RERconverge_Results/lm_residual_lamellae.enrichment.RDS")



#Convert enrichment list to dataframes

nb_enrich_group <- length(enrich_OR)
enrich_OR_df <- as.data.frame(NULL)
for(count in 1:nb_enrich_group){
  curr_df <- as.data.frame(enrich_OR[[count]])
  enrich_OR_df <- rbind(enrich_OR_df, curr_df)
}


nb_enrich_group <- length(enrich_TAAR)
enrich_TAAR_df <- as.data.frame(NULL)
for(count in 1:nb_enrich_group){
  curr_df <- as.data.frame(enrich_TAAR[[count]])
  enrich_TAAR_df <- rbind(enrich_TAAR_df, curr_df)
}


nb_enrich_group <- length(enrich_V2R)
enrich_V2R_df <- as.data.frame(NULL)
for(count in 1:nb_enrich_group){
  curr_df <- as.data.frame(enrich_V2R[[count]])
  enrich_V2R_df <- rbind(enrich_V2R_df, curr_df)
}


nb_enrich_group <- length(enrich_OLR)
enrich_OLR_df <- as.data.frame(NULL)
for(count in 1:nb_enrich_group){
  curr_df <- as.data.frame(enrich_OLR[[count]])
  enrich_OLR_df <- rbind(enrich_OLR_df, curr_df)
}


nb_enrich_group <- length(enrich_lamellae)
enrich_lamellae_df <- as.data.frame(NULL)
for(count in 1:nb_enrich_group){
  curr_df <- as.data.frame(enrich_lamellae[[count]])
  enrich_lamellae_df <- rbind(enrich_lamellae_df, curr_df)
}

nb_enrich_group <- length(enrich_res_lamellae)
enrich_res_lamellae_df <- as.data.frame(NULL)
for(count in 1:nb_enrich_group){
  curr_df <- as.data.frame(enrich_res_lamellae[[count]])
  enrich_res_lamellae_df <- rbind(enrich_res_lamellae_df, curr_df)
}


#Add databases names

rownames_enrich_OR_df <- row.names(enrich_OR_df)
rownames_enrich_OR_df <- gsub("_.*", "", rownames_enrich_OR_df)
enrich_OR_df <- 
  enrich_OR_df %>%
  mutate(database = rownames_enrich_OR_df)

rownames_enrich_TAAR_df <- row.names(enrich_TAAR_df)
rownames_enrich_TAAR_df <- gsub("_.*", "", rownames_enrich_TAAR_df)
enrich_TAAR_df <- 
  enrich_TAAR_df %>%
  mutate(database = rownames_enrich_TAAR_df)


rownames_enrich_V2R_df <- row.names(enrich_V2R_df)
rownames_enrich_V2R_df <- gsub("_.*", "", rownames_enrich_V2R_df)
enrich_V2R_df <- 
  enrich_V2R_df %>%
  mutate(database = rownames_enrich_V2R_df)


rownames_enrich_OLR_df <- row.names(enrich_OLR_df)
rownames_enrich_OLR_df <- gsub("_.*", "", rownames_enrich_OLR_df)
enrich_OLR_df <- 
  enrich_OLR_df %>%
  mutate(database = rownames_enrich_OLR_df)


rownames_enrich_lamellae_df <- row.names(enrich_lamellae_df)
rownames_enrich_lamellae_df <- gsub("_.*", "", rownames_enrich_lamellae_df)
enrich_lamellae_df <- 
  enrich_lamellae_df %>%
  mutate(database = rownames_enrich_lamellae_df)


rownames_enrich_res_lamellae_df <- row.names(enrich_res_lamellae_df)
rownames_enrich_res_lamellae_df <- gsub("_.*", "", rownames_enrich_res_lamellae_df)
enrich_res_lamellae_df <- 
  enrich_res_lamellae_df %>%
  mutate(database = rownames_enrich_res_lamellae_df)




##### Analyse RERconverge associations results -- Lamellae ---------------------------------

#Extract significant associations

lamellae_accelerated_df %>% arrange(p.adj)
lamellae_accelerated_df <- RER_lamellae %>% filter(category == "Accelerated")
lamellae_decelerated_df <- RER_lamellae %>% filter(category == "Decelerated")
lamellae_accelerated_genes <- RER_lamellae %>% filter(category == "Accelerated") %>% pull(gene_name)
lamellae_decelerated_genes <- RER_lamellae %>% filter(category == "Decelerated") %>% pull(gene_name)


RER_lamellae %>%
  ggplot(., aes(x=Rho, y=-log10(p.adj), color=category)) + 
  geom_point() +
  scale_color_manual(values = c("Non_signif" = "black", "Accelerated" = "#E69F00", 
                                "Decelerated" = "#56B4E9")) +
  geom_hline(yintercept = -log10(0.05), color="black", linetype="dashed") +
  geom_vline(xintercept = 0, color="black", linetype="dashed") +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position = "none")

RER_lamellae %>% filter(both_category == "Accelerated")
RER_lamellae %>% filter(both_category == "Decelerated")


RER_lamellae %>%
  filter(category == "Accelerated") %>%
  arrange(desc(Rho))

RER_lamellae %>%
  filter(category == "Decelerated") %>%
  arrange((Rho))


##### Analyse RERconverge enrichment results -- Lamellae ---------------------------------

enrich_lamellae_df %>% 
  filter(p.adj < 0.05) %>%
  group_by(database) %>%
  summarise(count = n())


enrich_lamellae_df %>% 
  filter(permpvaladj < 0.05) %>%
  group_by(database) %>%
  summarise(count = n())


enrich_lamellae_df %>% filter(stat < 0) %>% filter(p.adj < 0.05)
enrich_lamellae_df %>% filter(stat > 0) %>% filter(p.adj < 0.05)

enrich_lamellae_df %>% filter(stat < 0) %>% filter(p.adj < 0.05) %>% filter(permpval < 0.05) %>% arrange(p.adj)
enrich_lamellae_df %>% filter(stat > 0) %>% filter(p.adj < 0.05) %>% filter(permpval < 0.05) %>% arrange(p.adj)



enrich_lamellae_df %>% 
  filter(stat < 0) %>% 
  filter(p.adj < 0.05) %>%
  ggplot(., aes(x=p.adj, y=stat)) +
  geom_point()
  

enrich_lamellae_df$term_name <- rownames(enrich_lamellae_df)

enrich_lamellae_df %>% 
  filter(stat < 0) %>% 
  filter(p.adj < 0.05) %>%
  ggplot(., aes(x=term_name, y=p.adj, fill=p.adj)) +
  geom_raster()



enrich_test_decelerated <- 
  enrich_lamellae_df %>% 
  filter(stat < 0) %>% 
  filter(p.adj < 0.05) %>% 
  filter(permpval < 0.05) %>% 
  arrange(p.adj) %>%
  dplyr::select(p.adj, permpval)
enrich_test_decelerated$term <- row.names(enrich_test_decelerated)


enrich_test_decelerated_long <- 
  enrich_test_decelerated %>%
  pivot_longer(!term, names_to = "pvalue_cat", values_to = "value") %>%
  mutate(term_reduced = gsub("GOCC_", "", term))



enrich_test_decelerated_long %>%
  arrange(desc(term), value) %>% 
  ggplot(., 
       aes(y=term_reduced, x=pvalue_cat, fill= value)) + 
  geom_tile(color = "black") +
  geom_text(aes(label=round(value, 6))) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_gradient(low = "cadetblue3", high="dodgerblue4") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("") +
  ylab("") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


enrich_lamellae_df %>% 
  filter(p.adj < 0.05) %>% 
  filter(permpval < 0.05) %>%
  dplyr::select(stat, p.adj, permpval) %>%
  arrange(stat) %>%
  mutate(round_pvalue = round(p.adj, 3))

##### Analyse RERconverge associations results -- Res - Lamellae ---------------------------------

#Extract significant associations

lamellae_accelerated_df <- RER_res_lamellae %>% filter(category == "Accelerated")
lamellae_decelerated_df <- RER_res_lamellae %>% filter(category == "Decelerated")
lamellae_accelerated_genes <- RER_res_lamellae %>% filter(category == "Accelerated") %>% pull(gene_name)
lamellae_decelerated_genes <- RER_res_lamellae %>% filter(category == "Decelerated") %>% pull(gene_name)


RER_res_lamellae %>%
  ggplot(., aes(x=Rho, y=-log10(p.adj), color=category)) + 
  geom_point() +
  scale_color_manual(values = c("Non_signif" = "black", "Accelerated" = "#E69F00", 
                                "Decelerated" = "#56B4E9")) +
  geom_hline(yintercept = -log10(0.05), color="black", linetype="dashed") +
  geom_vline(xintercept = 0, color="black", linetype="dashed") +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        plot.subtitle=element_text(size=16)) +
  theme(legend.position = "none")


RER_res_lamellae %>% filter(both_category == "Accelerated")
RER_res_lamellae %>% filter(both_category == "Decelerated")


RER_res_lamellae %>%
  filter(category == "Accelerated") %>%
  arrange(desc(Rho))

RER_res_lamellae %>%
  filter(category == "Decelerated") %>%
  arrange((Rho))


##### Analyse RERconverge enrichment results --Res -  Lamellae ---------------------------------

enrich_res_lamellae_df %>% 
  filter(p.adj < 0.05) %>%
  group_by(database) %>%
  summarise(count = n())


enrich_res_lamellae_df %>% 
  filter(permpvaladj < 0.05) %>%
  group_by(database) %>%
  summarise(count = n())


enrich_res_lamellae_df %>% filter(stat < 0) %>% filter(p.adj < 0.05)
enrich_res_lamellae_df %>% filter(stat > 0) %>% filter(p.adj < 0.05)

enrich_res_lamellae_df %>% filter(stat < 0) %>% filter(p.adj < 0.05) %>% filter(permpval < 0.05)
enrich_res_lamellae_df %>% filter(stat > 0) %>% filter(p.adj < 0.05) %>% filter(permpval < 0.05)


enrich_lamellae_df %>% 
  filter(stat < 0) %>% 
  filter(p.adj < 0.05) %>%
  ggplot(., aes(x=p.adj, y=stat)) +
  geom_point()


enrich_lamellae_df$term_name <- rownames(enrich_lamellae_df)

enrich_lamellae_df %>% 
  filter(stat < 0) %>% 
  filter(p.adj < 0.05) %>%
  ggplot(., aes(x=term_name, y=p.adj, fill=p.adj)) +
  geom_raster()


