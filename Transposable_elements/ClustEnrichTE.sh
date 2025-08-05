#!/bin/bash


#SBATCH --job-name=EnrClust  # Job name


LMOD_DISABLE_SAME_NAME_AUTOSWAP=no
module purge
module load R


RepeatMaskerOutput=$1
GenomeFAI=$2
ClusterBed=$3
RepeatAlign=$4
ExonTSV=$5

### First lets create a TE bed file
grep -v "Simple_repeat" $RepeatMaskerOutput | grep -v "snRNA" | grep -v "Satellite" | grep -v "Low_complexity" | grep -v "tRNA" | grep -v "rRNA" |  sed 's/  */	/g' | sed 's/^	*//g' | cut -f5,6,7,11 | sed 's/\//_/g' | tail -n+4  > transposons.bed

#Lets create a genome BED file
cut -f1 $GenomeFAI | sed 's/$/	0/g' > scaffold_start
cut -f2 $GenomeFAI > stops
paste -d "\t" scaffold_start stops > Genome.bed ; rm stops ; rm scaffold_start

#Lets create a kimura distance table
grep "^[0-9]" $RepeatAlign | sed 's/ /,/g'  | cut -f5,6,7 -d "," > list_transposons_loc.csv
grep -A1 "Matrix" $RepeatAlign | sed 's/--//g' | sed '/^$/d' | grep -v "Matrix" | sed 's/Transitions.*//g' | sed 's/.*= //g' > list_kimura_distances.csv
paste -d "," list_transposons_loc.csv list_kimura_distances.csv > transposons_kimura_dist.csv

Rscript ClustEnrichTE.R $ClusterBed $ExonTSV

echo "Cluster Enrichment analysis done"

