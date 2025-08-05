#### We will try to find genes for which selective pressure is associated to the number of olfactory receptors
#### and to the number of lamellae in the epithelium. 


###==================================================================================================================================
###==================================================================================================================================
###==================================================================================================================================
###=========================================== Prepare files ====================================================================
###==================================================================================================================================
###==================================================================================================================================
###==================================================================================================================================


#First prepare exons of O. nil

grep "	CDS	"GCF_001858045.2_O_niloticus_UMD_NMBU_genomic.gff.LI > GCF_001858045.2_O_niloticus_UMD_NMBU_genomic.gff.exons

cut -f1 GCF_001858045.2_O_niloticus_UMD_NMBU_genomic.gff.exons > exons_scaffold_name.txt
cut -f4,5 GCF_001858045.2_O_niloticus_UMD_NMBU_genomic.gff.exons | tr "\t" "-" > exons_locations.txt
paste -d ":" exons_scaffold_name.txt exons_locations.txt > exons_scaffolds_locations.txt
cut -f9 GCF_001858045.2_O_niloticus_UMD_NMBU_genomic.gff.exons | sed 's/.*Parent=rna-//g' | sed 's/.*Parent=id-//g' | sed 's/;.*//g' > exons_XM_id.txt
cut -f9 GCF_001858045.2_O_niloticus_UMD_NMBU_genomic.gff.exons | sed 's/.*ID=cds-//g' | sed 's/;.*//g' > exons_XP_id.txt
cut -f7 GCF_001858045.2_O_niloticus_UMD_NMBU_genomic.gff.exons > strand_exons

paste -d "," exons_scaffolds_locations.txt exons_XM_id.txt exons_XP_id.txt > exons_positions_infos.csv
paste -d "," exons_positions_infos.csv strand_exons > temp ; mv temp exons_positions_infos.csv
awk -F',' 'BEGIN {OFS=","} {print $0, NR}' exons_positions_infos.csv > temp ; mv temp exons_positions_infos.csv



cut -f2,3 -d "," exons_positions_infos.csv | sed 's/,/-/g' | sort | uniq > uniq_Onil_genes.XP-XM.txt
cut -f2 -d "," exons_positions_infos.csv | sort | uniq > uniq_Onil_genes.XM.txt


cut -f2,4,5 -d "," exons_positions_infos.csv | sed 's/,/_/g' > new_names.txt
paste -d "\t" exons_scaffolds_locations.txt new_names.txt  > names_mapping.tsv
sed -i 's/:/-/g' names_mapping.tsv



#Prepare one ID per species


cut -f1 -d "," Table_ID_RER.csv > List_RER_ID.txt



###==================================================================================================================================
###==================================================================================================================================
###==================================================================================================================================
###===========================================Map reads and extract gene sequences from every individuals ===========================
###==================================================================================================================================
###==================================================================================================================================
###==================================================================================================================================


rm -r LogFiles ; mkdir LogFiles

#Align and call variant for each species

for curr_ID in `cat List_RER_ID.txt` ; do 
	sbatch --qos=6hours -c 40 --mem=50G -e LogFiles/error.$curr_ID.out -o LogFiles/slurm.$curr_ID.out Align_and_call_variants.sh $curr_ID
done

mkdir BCF_files

mv *.bcf BCF_files/ ; mv *.bcf.csi BCF_files/

#Create nice gene files

mkdir Genes_position_files
for gene in `cat uniq_Onil_genes.XM.txt` ; do
	grep ",$gene," exons_positions_infos.csv | cut -f1 -d "," > Genes_position_files/$gene.pos
	sed 's/:/-/g' Genes_position_files/$gene.pos > Genes_position_files/$gene.tiret.pos
	tac Genes_position_files/$gene.tiret.pos > Genes_position_files/$gene.tiret.pos.rev
done 


#Generate a consensus gene sequence for each species and for each gene


split -l 10000 --numeric-suffixes uniq_Onil_genes.XM.txt splitted_uniq_Onil_genes
split -l 40 --numeric-suffixes List_RER_ID.txt splitted_List_RER_ID.

rm -r Gene_LogFiles ; mkdir Gene_LogFiles
rm -r Consensus_Gene_Sequences/ ; mkdir Consensus_Gene_Sequences


#for curr_gene in `cat splitted_uniq_Onil_genes02` ; do 
#	sbatch --qos=6hours -c 2 --mem=5G -e Gene_LogFiles/error.$curr_gene.out -o Gene_LogFiles/slurm.$curr_gene.out Generate_consensus_exons.sh $curr_gene
#done

sbatch submit_array_generate_exons.sh 


rm -f Consensus_Gene_Sequences/*.rev
rm -f Consensus_Gene_Sequences/*.temp.fa.rev
rm -f Consensus_Gene_Sequences/*.temp.fa
rm -f Consensus_Gene_Sequences/*.exons.reformat.fa*
rm -f Consensus_Gene_Sequences/*.exons.fa
rm -f Consensus_Gene_Sequences/core*
rm -f Consensus_Gene_Sequences/sed*

grep -c ">" Consensus_Gene_Sequences/*.consensus.final.fa | grep -v ":264" | sed 's/Consensus_Gene_Sequences.//g' | sed 's/.consensus.final.*//g'  > genes_to_relaunch.txt

sbatch submit_array_generate_exons.relaunch.sh

for gene in `cat genes_to_relaunch.txt` ; do x=`grep -c ">" Consensus_Gene_Sequences/$gene.consensus.final.fa` ; echo "$gene,$x" ; done | grep -v ",264" | cut -f1 -d ","  > genes_to_relaunch.2.txt
cp  genes_to_relaunch.2.txt genes_to_relaunch.txt


sbatch submit_array_generate_exons.relaunch.sh


for gene in `cat genes_to_relaunch.txt` ; do x=`grep -c ">" Consensus_Gene_Sequences/$gene.consensus.final.fa` ; echo "$gene,$x" ; done | awk -F, '$2 > 200' | cut -f1 -d ","  > genes_to_relaunch.2.txt


sbatch submit_array_generate_exons.relaunch.sh


for gene in `cat genes_to_relaunch.2.txt` ; do grep -c ">" Consensus_Gene_Sequences/$gene.consensus.final.fa ; done



#====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
################## Make alignment and ML trees with phangorn #########################################################################c
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================



cp -r COFFRE_Consensus_Gene_Sequences/ Consensus_Gene_Sequences/

ls -l Consensus_Gene_Sequences/ | grep ".consensus.final.fa"  | sed 's/.* //g' | sed 's/.consensus.final.fa//g' > gene_to_align.txt



sbatch submit_array_align_trim.sh









======== submit_array_align_trim.sh ============================================

#!/bin/bash

#SBATCH --job-name=Align_trim 
#SBATCH --qos=6hours
#SBATCH -c 4
#SBATCH --mem=10G
#SBATCH --error=Gene_LogFiles/error.%A_%a.out
#SBATCH --output=Gene_LogFiles/slurm.%A_%a.out
#SBATCH --array=1-29903%2000

#Remove existing log files
rm -f Gene_LogFiles/error.*.out
rm -f Gene_LogFiles/slurm.*.out

# Extract the gene corresponding to the array task ID
gene=$(sed -n "${SLURM_ARRAY_TASK_ID}p" gene_to_align.txt)

# Call the script to process each gene
./RER_align_trim_phangorn.sh $gene




=============================================================================




rm genes_to_relaunch.txt
for gene in `cat gene_to_align.txt` ; do if [ ! -f Consensus_Gene_Sequences/$gene.consensus.final.prot.aln.trimmed.nwk ] ; then echo "$gene" >> genes_to_relaunch.txt ; fi ; done

rm really_relaunch.txt
for gene in `cat genes_to_relaunch.txt` ; do 
	nb_seq=`grep -c ">" Consensus_Gene_Sequences/$gene.consensus.final.prot`
	if [ $nb_seq -ge 3 ] ; then
		echo "$gene" >> really_relaunch.txt
	fi
done

for gene in `cat really_relaunch.txt` ; do rm Consensus_Gene_Sequences/$gene.* ; done
for gene in `cat really_relaunch.txt` ; do cp COFFRE_Consensus_Gene_Sequences/$gene.consensus.final.fa Consensus_Gene_Sequences/$gene.consensus.final.fa ; done

sbatch submit_array_align_trim.relaunch.sh



for gene in `cat really_relaunch.txt` ; do wc -l Consensus_Gene_Sequences/$gene.consensus.final.prot.aln.trimmed.nwk  ; done



rm AllTrees.Phangorn.nwk
for gene in `cat gene_to_align.txt` ; do
	newick_tree=$(cat "Consensus_Gene_Sequences/$gene.consensus.final.prot.aln.trimmed.nwk")
	echo -e "$gene\t$newick_tree" >> AllTrees.Phangorn.nwk
done
awk -F'\t' '$2 != ""' AllTrees.Phangorn.nwk > temp ; mv temp AllTrees.Phangorn.nwk



======== Estimate_ML_tree_Phangorn.cichlids.sh

#!/bin/bash


#SBATCH --job-name=Phangorn  # Job name


module purge
module load R/4.2.1-foss-2022a


aln_file=$1

Rscript Estimate_ML_tree_Phangorn.cichlids.R $aln_file $aln_file.nwk

echo "Nwk file generated"



======== Estimate_ML_tree_Phangorn.cichlids.R



##### Libraries  ---------------------------------

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



======================================================







#====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
################## Assign genes to GOterms / Gene names / KEGG pathways  #############################################################
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================


#First, find correspondances between Oreochromis niloticus genes and Zebrafish genes

for gene in `cat uniq_Onil_genes.XM.txt` ; do
	gene_ID=`echo "$gene" | sed 's/\..*//g'`
	danio_ID=`grep "$gene_ID" /scicore/home/salzburg/polica0000/Horizontal_transfer_project/OrthoFinder_ALL_annotations/Proteomes_BUSCO80/OrthoFinder/Results_Apr18/Phylogenetic_Hierarchical_Orthogroups/N5.tsv | cut -f68 | cut -f1 -d "," | sed 's/Danio_rerio---id-//g' | sed 's/Danio_rerio---rna-//g' | sed 's/^ *//g'`

	echo "$gene,$danio_ID" >> Onil_Zebrafish_ID.csv
done


#Find Goterms and pathways of Zebrafish genes

module load R/4.3.0-foss-2021a

>R
>library(msigdbr)
>library(dplyr)
>library(data.table)
>as.data.frame(msigdbr_collections())
>reactome_gene_sets <- msigdbr(species = "zebrafish", category = "C2", subcategory = "CP:REACTOME")
>BP_gene_sets <- msigdbr(species = "zebrafish", category = "C5", subcategory = "GO:BP")
>CC_gene_sets <- msigdbr(species = "zebrafish", category = "C5", subcategory = "GO:CC")
>MF_gene_sets <- msigdbr(species = "zebrafish", category = "C5", subcategory = "GO:MF")
>KEGG_gene_sets <- msigdbr(species = "zebrafish", category = "C2", subcategory = "CP:KEGG")
>reactome_gene_sets <- as.data.frame(reactome_gene_sets)
>BP_gene_sets <- as.data.frame(BP_gene_sets)
>CC_gene_sets <- as.data.frame(CC_gene_sets)
>MF_gene_sets <- as.data.frame(MF_gene_sets)
>KEGG_gene_sets <- as.data.frame(KEGG_gene_sets)

all_genes_sets <- msigdbr(species = "zebrafish")
all_genes_sets <- as.data.frame(all_genes_sets)
all_genes_sets <- all_genes_sets %>% dplyr::select(gene_symbol, entrez_gene, human_gene_symbol) %>% distinct()


>reactome_gene_sets <- reactome_gene_sets %>% dplyr::select(-c(ortholog_sources,num_ortholog_sources,taxon_id,gs_url,gs_geoid,gs_pmid,gs_id, gs_cat, gs_subcat))
>BP_gene_sets <- BP_gene_sets %>% dplyr::select(-c(ortholog_sources,num_ortholog_sources,taxon_id,gs_url,gs_geoid,gs_pmid,gs_id, gs_cat, gs_subcat))
>CC_gene_sets <- CC_gene_sets %>% dplyr::select(-c(ortholog_sources,num_ortholog_sources,taxon_id,gs_url,gs_geoid,gs_pmid,gs_id, gs_cat, gs_subcat))
>MF_gene_sets <- MF_gene_sets %>% dplyr::select(-c(ortholog_sources,num_ortholog_sources,taxon_id,gs_url,gs_geoid,gs_pmid,gs_id, gs_cat, gs_subcat))
>KEGG_gene_sets <- KEGG_gene_sets %>% dplyr::select(-c(ortholog_sources,num_ortholog_sources,taxon_id,gs_url,gs_geoid,gs_pmid,gs_id, gs_cat, gs_subcat))
>
>fwrite(reactome_gene_sets, "REACTOME.Danio_rerio.tsv", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
>fwrite(BP_gene_sets, "GO_BP.Danio_rerio.tsv", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
>fwrite(CC_gene_sets, "GO_CC.Danio_rerio.tsv", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
>fwrite(MF_gene_sets, "GO_MF.Danio_rerio.tsv", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
>fwrite(KEGG_gene_sets, "KEGG.Danio_rerio.tsv", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
>fwrite(all_genes_sets, "Danio_rerio.names.tsv", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)


cut -f3 REACTOME.Danio_rerio.tsv | tail -n+2 | sort | uniq > list_entrez.txt
cut -f3 GO_BP.Danio_rerio.tsv | tail -n+2 | sort | uniq >> list_entrez.txt
cut -f3 GO_CC.Danio_rerio.tsv | tail -n+2 | sort | uniq >> list_entrez.txt
cut -f3 GO_MF.Danio_rerio.tsv | tail -n+2 | sort | uniq >> list_entrez.txt
cut -f3 KEGG.Danio_rerio.tsv | tail -n+2 | sort | uniq >> list_entrez.txt
sort list_entrez.txt | uniq > temp ; mv temp list_entrez.txt


rm Danio_rerio.EntrezMapping.csv
for curr_ID in `cat list_entrez.txt` ; do 

	gene_name=`grep "GeneID:$curr_ID,"  /scicore/home/salzburg/polica0000/Horizontal_transfer_project/Genomic_data/Danio_rerio/GCF_000002035.6_GRCz11_genomic.gff.LI | grep -m1 "	mRNA	" | sed 's/;Parent.*//g' | sed 's/.*=//g'`

	echo "$curr_ID,$gene_name" >> Danio_rerio.EntrezMapping.csv

done

sed -i "s/'//g" REACTOME.Danio_rerio.tsv
sed -i "s/'//g" GO_BP.Danio_rerio.tsv
sed -i "s/'//g" GO_CC.Danio_rerio.tsv
sed -i "s/'//g" GO_MF.Danio_rerio.tsv
sed -i "s/'//g" KEGG.Danio_rerio.tsv



>R
>library(dplyr)
>library(tidyr)

entrez_mapping_df <- read.table("Danio_rerio.EntrezMapping.csv", sep=",", header=FALSE)
colnames(entrez_mapping_df) <- c("entrez_gene", "dr_gene_name")
entrez_mapping_df$dr_gene_name <- gsub("\\.*", "", entrez_mapping_df$dr_gene_name)
entrez_mapping_df$dr_gene_name <- gsub("rna-", "", entrez_mapping_df$dr_gene_name)

onil_drerio_mapping_df <- read.table("Onil_Zebrafish_ID.csv", sep=",", header=FALSE)
colnames(onil_drerio_mapping_df) <- c("Onil_gene_name", "dr_gene_name")
onil_drerio_mapping_df$dr_gene_name <- gsub("\\.*", "", onil_drerio_mapping_df$dr_gene_name)
onil_drerio_mapping_df <- onil_drerio_mapping_df %>% filter(! is.na(dr_gene_name)) %>% filter(dr_gene_name != "")



reactome_df <- read.table("REACTOME.Danio_rerio.tsv", sep="\t", header=TRUE)
reactome_df.entrez  <- left_join(entrez_mapping_df, reactome_df, by="entrez_gene")
reactome_df.entrez.onil  <- left_join(reactome_df.entrez, onil_drerio_mapping_df, by="dr_gene_name")
reactome_df.entrez.onil.filt <- reactome_df.entrez.onil %>% dplyr::select(Onil_gene_name, gs_description, gs_name) %>% distinct()

BP_df <- read.table("GO_BP.Danio_rerio.tsv", sep="\t", header=TRUE)
BP_df.entrez  <- left_join(entrez_mapping_df, BP_df, by="entrez_gene")
BP_df.entrez.onil  <- left_join(BP_df.entrez, onil_drerio_mapping_df, by="dr_gene_name")
BP_df.entrez.onil.filt <- BP_df.entrez.onil %>% dplyr::select(Onil_gene_name, gs_description, gs_name) %>% distinct()

CC_df <- read.table("GO_CC.Danio_rerio.tsv", sep="\t", header=TRUE)
CC_df.entrez  <- left_join(entrez_mapping_df, CC_df, by="entrez_gene")
CC_df.entrez.onil  <- left_join(CC_df.entrez, onil_drerio_mapping_df, by="dr_gene_name")
CC_df.entrez.onil.filt <- CC_df.entrez.onil %>% dplyr::select(Onil_gene_name, gs_description, gs_name) %>% distinct()

MF_df <- read.table("GO_MF.Danio_rerio.tsv", sep="\t", header=TRUE)
MF_df.entrez  <- left_join(entrez_mapping_df, MF_df, by="entrez_gene")
MF_df.entrez.onil  <- left_join(MF_df.entrez, onil_drerio_mapping_df, by="dr_gene_name")
MF_df.entrez.onil.filt <- MF_df.entrez.onil %>% dplyr::select(Onil_gene_name, gs_description, gs_name) %>% distinct()

KEGG_df <- read.table("KEGG.Danio_rerio.tsv", sep="\t", header=TRUE)
KEGG_df.entrez  <- left_join(entrez_mapping_df, KEGG_df, by="entrez_gene")
KEGG_df.entrez.onil  <- left_join(KEGG_df.entrez, onil_drerio_mapping_df, by="dr_gene_name")
KEGG_df.entrez.onil.filt <- KEGG_df.entrez.onil %>% dplyr::select(Onil_gene_name, gs_description, gs_name) %>% distinct()


reactome_df.entrez.onil.filt <- 
reactome_df.entrez.onil.filt %>% 
dplyr::filter(! is.na(Onil_gene_name)) %>% 
dplyr::filter(! is.na(gs_name)) %>% 
dplyr::filter(! is.na(gs_description)) %>% 
dplyr::select(gs_name, gs_description, Onil_gene_name) %>%
distinct()


BP_df.entrez.onil.filt <- 
BP_df.entrez.onil.filt %>% 
dplyr::filter(! is.na(Onil_gene_name)) %>% 
dplyr::filter(! is.na(gs_name)) %>% 
dplyr::filter(! is.na(gs_description)) %>% 
dplyr::select(gs_name, gs_description, Onil_gene_name) %>%
distinct()


CC_df.entrez.onil.filt <- 
CC_df.entrez.onil.filt %>% 
dplyr::filter(! is.na(Onil_gene_name)) %>% 
dplyr::filter(! is.na(gs_name)) %>% 
dplyr::filter(! is.na(gs_description)) %>% 
dplyr::select(gs_name, gs_description, Onil_gene_name) %>%
distinct()


MF_df.entrez.onil.filt <- 
MF_df.entrez.onil.filt %>% 
dplyr::filter(! is.na(Onil_gene_name)) %>% 
dplyr::filter(! is.na(gs_name)) %>% 
dplyr::filter(! is.na(gs_description)) %>% 
dplyr::select(gs_name, gs_description, Onil_gene_name) %>%
distinct()


KEGG_df.entrez.onil.filt <- 
KEGG_df.entrez.onil.filt %>% 
dplyr::filter(! is.na(Onil_gene_name)) %>% 
dplyr::filter(! is.na(gs_name)) %>% 
dplyr::filter(! is.na(gs_description)) %>% 
dplyr::select(gs_name, gs_description, Onil_gene_name) %>%
distinct()

reactome.GMT <- as.data.frame(reactome_df.entrez.onil.filt %>% group_by(gs_name, gs_description) %>%  summarize(gene_name = paste(Onil_gene_name, collapse = "\t")))
BP.GMT <- as.data.frame(BP_df.entrez.onil.filt %>% group_by(gs_name, gs_description) %>%  summarize(gene_name = paste(Onil_gene_name, collapse = "\t")))
CC.GMT <- as.data.frame(CC_df.entrez.onil.filt %>% group_by(gs_name, gs_description) %>%  summarize(gene_name = paste(Onil_gene_name, collapse = "\t")))
MF.GMT <- as.data.frame(MF_df.entrez.onil.filt %>% group_by(gs_name, gs_description) %>%  summarize(gene_name = paste(Onil_gene_name, collapse = "\t")))
KEGG.GMT <- as.data.frame(KEGG_df.entrez.onil.filt %>% group_by(gs_name, gs_description) %>%  summarize(gene_name = paste(Onil_gene_name, collapse = "\t")))

write.table(reactome.GMT, "REACTOME.Onil.gmt", col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
write.table(BP.GMT, "BP.Onil.gmt", col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
write.table(CC.GMT, "CC.Onil.gmt", col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
write.table(MF.GMT, "MF.Onil.gmt", col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
write.table(KEGG.GMT, "KEGG.Onil.gmt", col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")


write.table(reactome_df.entrez.onil.filt, "REACTOME.Onil.tsv", col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
write.table(BP_df.entrez.onil.filt, "BP.Onil.tsv", col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
write.table(CC_df.entrez.onil.filt, "CC.Onil.tsv", col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
write.table(MF_df.entrez.onil.filt, "MF.Onil.tsv", col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
write.table(KEGG_df.entrez.onil.filt, "KEGG.Onil.tsv", col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")


danio_names_df <- read.table("Danio_rerio.names.tsv", sep="\t", header=TRUE)
danio_names_df  <- left_join(danio_names_df, entrez_mapping_df, by="entrez_gene")
Onil_names_df  <- left_join(onil_drerio_mapping_df, danio_names_df, by="dr_gene_name")
write.table(Onil_names_df, "Names.Onil.tsv", col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")



#====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
################## Launch RERconverge  #############################################################
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================
##====================================================================================================================================



rm -rf RERconverge_Results ; mkdir RERconverge_Results


sbatch --qos=1week -c 12 --mem=30G launch_RER.sh
sbatch --qos=1day -c 12 --mem=30G launch_RER.OR.sh
sbatch --qos=1day -c 12 --mem=30G launch_RER.TAAR.sh
sbatch --qos=1day -c 12 --mem=30G launch_RER.V2R.sh
sbatch --qos=1day -c 12 --mem=30G launch_RER.OLR.sh
sbatch --qos=1day -c 12 --mem=30G launch_RER.lamellae.sh
sbatch --qos=1day -c 12 --mem=30G launch_RER.pGLSlamellae.sh
sbatch --qos=1day -c 12 --mem=30G launch_RER.Lmlamellae.sh





======= launch_RER.sh ==========================================

#!/bin/bash


#SBATCH --job-name=RERconverge  # Job name


module load R/4.4.1-foss-2023b

Rscript RER_cichlids.R

===============================================================
