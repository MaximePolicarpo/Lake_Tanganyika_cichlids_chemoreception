#!/bin/bash


#SBATCH --job-name=Consensus_Exons  # Job name



eval "$(conda shell.bash hook)"
conda activate miniprot

LMOD_DISABLE_SAME_NAME_AUTOSWAP=no

module purge  ; module load bwa-mem2


ID=$1
ref_genome=GCF_001858045.2_O_niloticus_UMD_NMBU_genomic.fna



R1_A_file=`ls -l /scicore/home/salzburg/GROUP/genomedata_cichlidX/Sample_${ID}_A/ | grep ".*R1.*fastq.gz" | sed 's/.* //g'`
R1_B_file=`ls -l /scicore/home/salzburg/GROUP/genomedata_cichlidX/Sample_${ID}_B/ | grep ".*R1.*fastq.gz" | sed 's/.* //g'`
R2_A_file=`ls -l /scicore/home/salzburg/GROUP/genomedata_cichlidX/Sample_${ID}_A/ | grep ".*R2.*fastq.gz" | sed 's/.* //g'`
R2_B_file=`ls -l /scicore/home/salzburg/GROUP/genomedata_cichlidX/Sample_${ID}_B/ | grep ".*R2.*fastq.gz" | sed 's/.* //g'`

cat /scicore/home/salzburg/GROUP/genomedata_cichlidX/Sample_${ID}_A/$R1_A_file /scicore/home/salzburg/GROUP/genomedata_cichlidX/Sample_${ID}_B/$R1_B_file > $ID.R1.fastq.gz
cat /scicore/home/salzburg/GROUP/genomedata_cichlidX/Sample_${ID}_A/$R2_A_file /scicore/home/salzburg/GROUP/genomedata_cichlidX/Sample_${ID}_B/$R2_B_file > $ID.R2.fastq.gz


bwa-mem2 mem -t 40 $ref_genome $ID.R1.fastq.gz $ID.R2.fastq.gz > $ID.Onil.sam


rm $ID.R1.fastq.gz ; rm $ID.R2.fastq.gz



module purge ; module load SAMtools

#Convert SAM to BAM
samtools sort $ID.Onil.sam -o $ID.sorted.bam ; samtools index $ID.sorted.bam ; rm $ID.Onil.sam
