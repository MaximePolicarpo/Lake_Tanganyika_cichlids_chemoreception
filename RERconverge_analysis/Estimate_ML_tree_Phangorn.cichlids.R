rm(list=ls())

set.seed(2712)

library("ape")
library(dplyr)
library(tidyverse)
library(data.table)
library(phytools)
library(purrr)
library("RERconverge")

args = commandArgs(trailingOnly=TRUE)

my_aln_file <- args[1]
species_tree <- "b1_tree_wo_Neospl.nwk"


MLrecon_phangorn <- 
  estimatePhangornTree(
    my_aln_file, 
    species_tree,
    submodel = "LG",
    type = "AA",
    format = "fasta",
    k = 4)

MLrecon_phangorn_tree <- MLrecon_phangorn$tree.opt

write.tree(MLrecon_phangorn_tree, args[2])