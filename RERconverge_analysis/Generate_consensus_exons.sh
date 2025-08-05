#!/bin/bash


#SBATCH --job-name=Consensus_Exons  # Job name



eval "$(conda shell.bash hook)"
conda activate miniprot

LMOD_DISABLE_SAME_NAME_AUTOSWAP=no


gene=$1 #gene=LOC100695241
ref_genome=/scicore/home/salzburg/polica0000/Cichlids_Genomes/Oreochromis_niloticus/GCF_001858045.2_O_niloticus_UMD_NMBU_genomic.fna


##Generate a consensus of exon sequences for each gene and each species, and merge exons

module purge ; module load SAMtools ; module load BCFtools


#rm Consensus_Gene_Sequences/$gene.consensus.final.fa


#extract the exon positions of the gene
grep ",$gene," exons_positions_infos.csv | cut -f1 -d "," > $gene.pos
sed 's/:/-/g' $gene.pos > $gene.tiret.pos
tac $gene.tiret.pos > $gene.tiret.pos.rev

#extract the gene strand
strand=`grep "$gene" GCF_001858045.2_O_niloticus_UMD_NMBU_genomic.gff.exons | cut -f7 | head -1`


rm -f Gene_LogFiles/slurm.$gene.out
rm -f Gene_LogFiles/error.$gene.out

for ID in `cat splitted_List_RER_ID.00` ; do 

	if grep -q ">$ID" Consensus_Gene_Sequences/$gene.consensus.final.fa ; then

		echo "$ID already in file"

	else 

		#Extract current exons consensus in the current individual
		xargs samtools faidx $ref_genome < $gene.pos | bcftools consensus $ID.norm.flt-indels.bcf >> Consensus_Gene_Sequences/$gene.$ID.exons.fa
		
		#Reformat the fasta file for samtools
		sed -i 's/:/-/g'  Consensus_Gene_Sequences/$gene.$ID.exons.fa
		fasta_formatter -i Consensus_Gene_Sequences/$gene.$ID.exons.fa -o Consensus_Gene_Sequences/$gene.$ID.exons.reformat.fa -w 60 ; rm Consensus_Gene_Sequences/$gene.$ID.exons.fa
	
	
		if [ $strand == "+" ] ; then
			xargs samtools faidx Consensus_Gene_Sequences/$gene.$ID.exons.reformat.fa < $gene.tiret.pos | grep -v ">" | sed "1i\>$ID"  >> Consensus_Gene_Sequences/$gene.consensus.final.fa
		else
	
			for curr_exon in `cat $gene.tiret.pos.rev` ; do 
	
				samtools faidx Consensus_Gene_Sequences/$gene.$ID.exons.reformat.fa $curr_exon | grep -v ">" > Consensus_Gene_Sequences/$gene.$ID.temp.fa
				revseq Consensus_Gene_Sequences/$gene.$ID.temp.fa Consensus_Gene_Sequences/$gene.$ID.temp.fa.rev
				grep -v ">" Consensus_Gene_Sequences/$gene.$ID.temp.fa.rev >> Consensus_Gene_Sequences/$gene.$ID.rev
	
			done
	
			sed "1i\>$ID" Consensus_Gene_Sequences/$gene.$ID.rev >> Consensus_Gene_Sequences/$gene.consensus.final.fa
	
			rm Consensus_Gene_Sequences/$gene.$ID.rev
			rm Consensus_Gene_Sequences/$gene.$ID.temp.fa.rev
			rm Consensus_Gene_Sequences/$gene.$ID.temp.fa
		fi
	
	
		rm Consensus_Gene_Sequences/$gene.$ID.exons.reformat.fa*


	fi

done



rm $gene.tiret.pos.rev
rm $gene.pos
rm $gene.tiret.pos



##Generate a consensus of exon sequences for each gene and each species, and merge exons
#module purge ; module load SAMtools ; module load BCFtools
#
#
#rm Consensus_Gene_Sequences/$ID.ConsensusGenes.fa
#xargs samtools faidx $ref_genome < exons_scaffolds_locations.txt | bcftools consensus $ID.norm.flt-indels.bcf >> Consensus_Gene_Sequences/$ID.ConsensusGenes.fa
#sed -i 's/:/-/g'  Consensus_Gene_Sequences/$ID.ConsensusGenes.fa
#perl rename_fasta.pl names_mapping.tsv Consensus_Gene_Sequences/$ID.ConsensusGenes.fa > Consensus_Gene_Sequences/$ID.ConsensusGenes.renamed.fa
#fasta_formatter -i Consensus_Gene_Sequences/$ID.ConsensusGenes.renamed.fa -o Consensus_Gene_Sequences/$ID.ConsensusGenes.renamed.reformat.fa -w 60
#
#
#for gene in `cat uniq_Onil_genes.XM.txt` ; do 
#
#	if grep -q "$gene" Consensus_Gene_Sequences/$ID.ConsensusGenes.renamed.fa ; then
#		grep ">$gene" Consensus_Gene_Sequences/$ID.ConsensusGenes.renamed.fa | sed 's/>//g' > curr_gene.$ID.id
#		strand=`grep "$gene" GCF_001858045.2_O_niloticus_UMD_NMBU_genomic.gff.exons | cut -f7 | head -1`
#		
#		if [ $strand == "-" ] ; then
#			tac curr_gene.$ID.id > curr_gene.$ID.id.rev
#		
#			for curr_exon in `cat curr_gene.$ID.id.rev` ; do 
#		
#				samtools faidx Consensus_Gene_Sequences/$ID.ConsensusGenes.renamed.reformat.fa $curr_exon | grep -v ">" > curr_gene.$ID.fa
#				revseq curr_gene.$ID.fa curr_gene.$ID.fa.rev 
#				grep -v ">" curr_gene.$ID.fa.rev >> curr_gene.$ID.finalrev
#				
#			done
#		
#			sed -i "1i\>$gene" curr_gene.$ID.finalrev
#			cat curr_gene.$ID.finalrev >> Consensus_Gene_Sequences/$ID.ConsensusGenes.final.fa 
#			rm curr_gene.$ID.id.rev ; rm curr_gene.$ID.fa ; rm curr_gene.$ID.fa.rev  ; rm curr_gene.$ID.finalrev 
#	
#		else 
#		
#			xargs samtools faidx Consensus_Gene_Sequences/$ID.ConsensusGenes.renamed.reformat.fa < curr_gene.$ID.id | grep -v ">" | sed "1i\>$gene"  >> Consensus_Gene_Sequences/$ID.ConsensusGenes.final.fa 
#		fi
#	
#		rm curr_gene.$ID.id
#	fi
#done
#
#
#fasta_formatter -i Consensus_Gene_Sequences/$ID.ConsensusGenes.final.fa  -o Consensus_Gene_Sequences/$ID.ConsensusGenes.final.reformat.fa  -w 60 
#mv Consensus_Gene_Sequences/$ID.ConsensusGenes.final.reformat.fa Consensus_Gene_Sequences/$ID.ConsensusGenes.final.fa
#transeq Consensus_Gene_Sequences/$ID.ConsensusGenes.final.fa Consensus_Gene_Sequences/$ID.ConsensusGenes.final.prot ; sed -i 's/_1$//g' Consensus_Gene_Sequences/$ID.ConsensusGenes.final.prot
#awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' Consensus_Gene_Sequences/$ID.ConsensusGenes.final.prot | sed 's/\*$//g' | awk -F '\t'  '!($2 ~ /\*/)' | tr "\t" "\n" > Consensus_Gene_Sequences/$ID.ConsensusGenes.final.nostop.prot
#grep ">" Consensus_Gene_Sequences/$ID.ConsensusGenes.final.nostop.prot | sed 's/>//g' > $ID.nostop.id
#xargs samtools faidx Consensus_Gene_Sequences/$ID.ConsensusGenes.final.fa < $ID.nostop.id > Consensus_Gene_Sequences/$ID.ConsensusGenes.final.nostop.fa
#
#rm $ID.nostop.id
#rm $ID.norm.flt-indels.bcf
#rm $ID.norm.flt-indels.bcf.csi
#