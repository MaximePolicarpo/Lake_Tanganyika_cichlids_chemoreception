#!/bin/bash


#SBATCH --job-name=Consensus_Exons  # Job name


eval "$(conda shell.bash hook)"
conda activate miniprot

LMOD_DISABLE_SAME_NAME_AUTOSWAP=no

module purge ; module load SAMtools 


gene=$1

#Reformat the gene file and rename ID to Sp
fasta_formatter -i Consensus_Gene_Sequences/$gene.consensus.final.fa -o Consensus_Gene_Sequences/$gene.consensus.final.reformat.fa -w 60 
mv Consensus_Gene_Sequences/$gene.consensus.final.reformat.fa Consensus_Gene_Sequences/$gene.consensus.final.fa

transeq Consensus_Gene_Sequences/$gene.consensus.final.fa Consensus_Gene_Sequences/$gene.consensus.final.prot ; sed -i 's/_1$//g' Consensus_Gene_Sequences/$gene.consensus.final.prot
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' Consensus_Gene_Sequences/$gene.consensus.final.prot | sed 's/\*$//g' | awk -F '\t'  '!($2 ~ /\*/)' | tr "\t" "\n" > Consensus_Gene_Sequences/$gene.consensus.final.prot.nostop
grep ">" Consensus_Gene_Sequences/$gene.consensus.final.prot.nostop | sed 's/>//g' > Consensus_Gene_Sequences/$gene.nostop.id
xargs samtools faidx Consensus_Gene_Sequences/$gene.consensus.final.fa < Consensus_Gene_Sequences/$gene.nostop.id > Consensus_Gene_Sequences/$gene.consensus.final.fa.nostop
rm Consensus_Gene_Sequences/$gene.nostop.id

mv Consensus_Gene_Sequences/$gene.consensus.final.fa.nostop Consensus_Gene_Sequences/$gene.consensus.final.fa 
mv Consensus_Gene_Sequences/$gene.consensus.final.prot.nostop Consensus_Gene_Sequences/$gene.consensus.final.prot
rm -f Consensus_Gene_Sequences/$gene.consensus.final.fa.fai



awk 'NR==FNR{a[$1]=$2;next} NF==2{$2=a[$2]; print ">" $2;next} 1' FS=',' Table_ID_RER.csv FS='>' Consensus_Gene_Sequences/$gene.consensus.final.fa > Consensus_Gene_Sequences/$gene.consensus.final.renamed.fa
mv Consensus_Gene_Sequences/$gene.consensus.final.renamed.fa Consensus_Gene_Sequences/$gene.consensus.final.fa


#awk -F',' '{print "s/"$1"/"$2"/g"}' Table_ID_RER.csv | xargs -I {} sed -i '{}' Consensus_Gene_Sequences/$gene.consensus.final.fa

transeq Consensus_Gene_Sequences/$gene.consensus.final.fa Consensus_Gene_Sequences/$gene.consensus.final.prot ; sed -i 's/_1$//g' Consensus_Gene_Sequences/$gene.consensus.final.prot


#for line in `cat Table_ID_RER.csv` ; do
#	ID=`echo "$line" | cut -f1 -d ","`
#	Sp=`echo "$line" | cut -f2 -d ","`
#
#	sed -i "s/>$ID/>$Sp/g" Consensus_Gene_Sequences/$gene.consensus.final.prot
#	sed -i "s/>$ID/>$Sp/g" Consensus_Gene_Sequences/$gene.consensus.final.fa
#done



#Perform alignment and trimming

grep ">" Consensus_Gene_Sequences/$gene.consensus.final.prot | sed 's/>//g' | sort > Consensus_Gene_Sequences/$gene.allid.txt
grep ">" Consensus_Gene_Sequences/$gene.consensus.final.prot | shuf -n25 | sed 's/>//g' | sort > Consensus_Gene_Sequences/$gene.backboneid.txt
comm -23 Consensus_Gene_Sequences/$gene.allid.txt Consensus_Gene_Sequences/$gene.backboneid.txt  > Consensus_Gene_Sequences/$gene.toaddid.txt

xargs samtools faidx Consensus_Gene_Sequences/$gene.consensus.final.prot < Consensus_Gene_Sequences/$gene.backboneid.txt > Consensus_Gene_Sequences/$gene.backboneid.prot
xargs samtools faidx Consensus_Gene_Sequences/$gene.consensus.final.prot < Consensus_Gene_Sequences/$gene.toaddid.txt > Consensus_Gene_Sequences/$gene.toaddid.prot 

rm -f Consensus_Gene_Sequences/$gene.backboneid.txt ; rm -f Consensus_Gene_Sequences/$gene.allid.txt ; rm -f Consensus_Gene_Sequences/$gene.toaddid.txt


/scicore/home/salzburg/polica0000/Non_visual_opsins_Project/Opsins_analysis_2024/muscle5.1.linux_intel64 -align Consensus_Gene_Sequences/$gene.backboneid.prot -output Consensus_Gene_Sequences/$gene.backboneid.prot.aln
rm -f Consensus_Gene_Sequences/$gene.backboneid.prot

mafft --add Consensus_Gene_Sequences/$gene.toaddid.prot --keeplength Consensus_Gene_Sequences/$gene.backboneid.prot.aln > Consensus_Gene_Sequences/$gene.consensus.final.prot.aln
rm -f Consensus_Gene_Sequences/$gene.backboneid.prot.aln ; rm -f Consensus_Gene_Sequences/$gene.toaddid.prot

trimal -in Consensus_Gene_Sequences/$gene.consensus.final.prot.aln -automated1 -out Consensus_Gene_Sequences/$gene.consensus.final.prot.aln.trimmed


rm Consensus_Gene_Sequences/$gene.consensus.final.prot.fai

#Estilmate ML tree with phangorn


module purge ; module load R/4.2.1-foss-2022a

./Estimate_ML_tree_Phangorn.cichlids.sh Consensus_Gene_Sequences/$gene.consensus.final.prot.aln.trimmed 


