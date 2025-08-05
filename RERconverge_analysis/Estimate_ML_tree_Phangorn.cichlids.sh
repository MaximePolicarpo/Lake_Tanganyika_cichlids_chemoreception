#!/bin/bash


#SBATCH --job-name=Phangorn  # Job name


module purge
module load R/4.2.1-foss-2022a


aln_file=$1

Rscript Estimate_ML_tree_Phangorn.cichlids.R $aln_file $aln_file.nwk

echo "Nwk file generated"