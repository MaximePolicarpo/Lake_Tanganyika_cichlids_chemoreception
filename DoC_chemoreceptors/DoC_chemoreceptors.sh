###############################################################################################################################################
###############################################################################################################################################
###############################################################################################################################################
####################################################Compute the depth of coverage for every inviduals #########################################
###############################################################################################################################################
###############################################################################################################################################

module load GATK/4.2.6.1-foss-2018b-Java-1.8
module load SAMtools/1.15-GCC-10.3.0
module load bwa-mem2

mkdir Onil_ALL_bam_coverages_stats

ID=$1

IFS=$'\n'

R1_A_file=`ls -l Sample_${ID}_A/ | grep ".*R1.*fastq.gz" | sed 's/.* //g'`
R1_B_file=`ls -l Sample_${ID}_B/ | grep ".*R1.*fastq.gz" | sed 's/.* //g'`
R2_A_file=`ls -l Sample_${ID}_A/ | grep ".*R2.*fastq.gz" | sed 's/.* //g'`
R2_B_file=`ls -l Sample_${ID}_B/ | grep ".*R2.*fastq.gz" | sed 's/.* //g'`

cat Sample_${ID}_A/$R1_A_file Sample_${ID}_B/$R1_B_file > $ID.R1.fastq.gz
cat Sample_${ID}_A/$R2_A_file Sample_${ID}_B/$R2_B_file > $ID.R2.fastq.gz


#Align reads with bwa-mem
bwa-mem2 mem -t 40 GCF_001858045.2_O_niloticus_UMD_NMBU_genomic.fna $ID.R1.fastq.gz $ID.R2.fastq.gz > $ID.Onil.sam
rm $ID.R1.fastq.gz ; rm $ID.R2.fastq.gz
samtools view -Sb $ID.Onil.sam > $ID.Onil.bam
rm $ID.Onil.sam ; rm $ID.OnilMasked.sam
samtools sort -o $ID.Onil.sorted.bam $ID.Onil.bam 
samtools index $ID.Onil.sorted.bam
rm $ID.Onil.bam ; rm $ID.OnilMasked.bam


#Remove duplicate reads
gatk --java-options "-Xmx5g" MarkDuplicates -I $ID.Onil.sorted.bam -O $ID.Onil.sorted.markdup.bam --REMOVE_DUPLICATES TRUE -M Onil_Bam_duplicates_removed/$ID.metrics_dup --TMP_DIR Temporary_FILES/ 

#Get the total number of reads in the BAM file
total_reads=`samtools view -c $ID.Onil.sorted.bam`
#Get the total number of mapped reads in the BAM file
mapped_reads=`samtools view -c -F 260 $ID.Onil.sorted.bam`
#Percentage of mapped reads
calc(){ awk "BEGIN { print "$*" }"; }
perc_mapped=`calc $mapped_reads/$total_reads`


echo "$total_reads,$mapped_reads,$perc_mapped" > Onil_ALL_bam_coverages_stats/$ID.nbreads.txt


#Get the coverage on OR genes
samtools depth -a $ID.Onil.sorted.markdup.bam -q 0 -b Onil.OR.all.pos  > Onil_ALL_bam_coverages_stats/$ID.OR.unmasked.raw.depth

rm Onil_ALL_bam_coverages_stats/$ID.OR.unmasked.mean.depth
for line in `cat Oreochromis_nilotius_OR_exons_coordinates.tsv` ; do
	scaffold=`echo "$line" | cut -f6`
	start=`echo "$line" | cut -f7`
	end=`echo "$line" | cut -f8`
	mean_coverage=`awk -v myvar_scaff="$scaffold" '{ if ($1 == myvar_scaff) { print } }' Onil_ALL_bam_coverages_stats/$ID.OR.unmasked.raw.depth | awk -v myvar_start="$start" '{ if ($2 >= myvar_start) { print } }' | awk -v myvar_end="$end" '{ if ($2 <= myvar_end) { print } }'  | awk '{D+=int($3)} END{print D/(NR*1.0);}'`
	echo "$line	$ID	$mean_coverage" >> Onil_ALL_bam_coverages_stats/$ID.OR.unmasked.mean.depth
done



#Get the coverage on TAAR genes
samtools depth -a $ID.Onil.sorted.markdup.bam -q 0 -b Onil.TAAR.all.pos > Onil_ALL_bam_coverages_stats/$ID.TAAR.unmasked.raw.depth

rm Onil_ALL_bam_coverages_stats/$ID.TAAR.unmasked.mean.depth
for line in `cat Oreochromis_nilotius_TAAR_exons_coordinates.tsv` ; do
	scaffold=`echo "$line" | cut -f6`
	start=`echo "$line" | cut -f7`
	end=`echo "$line" | cut -f8`
	mean_coverage=`awk -v myvar_scaff="$scaffold" '{ if ($1 == myvar_scaff) { print } }' Onil_ALL_bam_coverages_stats/$ID.TAAR.unmasked.raw.depth | awk -v myvar_start="$start" '{ if ($2 >= myvar_start) { print } }' | awk -v myvar_end="$end" '{ if ($2 <= myvar_end) { print } }'  | awk '{D+=int($3)} END{print D/(NR*1.0);}'`
	echo "$line	$ID	$mean_coverage" >> Onil_ALL_bam_coverages_stats/$ID.TAAR.unmasked.mean.depth
done


#Get the coverage on V2R genes
samtools depth -a $ID.Onil.sorted.markdup.bam -q 0 -b Onil.V2R.all.pos  > Onil_ALL_bam_coverages_stats/$ID.V2R.unmasked.raw.depth

rm Onil_ALL_bam_coverages_stats/$ID.V2R.unmasked.mean.depth
for line in `cat Oreochromis_nilotius_V2R_exons_coordinates.tsv` ; do
	scaffold=`echo "$line" | cut -f6`
	start=`echo "$line" | cut -f7`
	end=`echo "$line" | cut -f8`
	mean_coverage=`awk -v myvar_scaff="$scaffold" '{ if ($1 == myvar_scaff) { print } }' Onil_ALL_bam_coverages_stats/$ID.V2R.unmasked.raw.depth | awk -v myvar_start="$start" '{ if ($2 >= myvar_start) { print } }' | awk -v myvar_end="$end" '{ if ($2 <= myvar_end) { print } }'  | awk '{D+=int($3)} END{print D/(NR*1.0);}'`
	echo "$line	$ID	$mean_coverage" >> Onil_ALL_bam_coverages_stats/$ID.V2R.unmasked.mean.depth
done


#Get the coverage on V1R genes
samtools depth -a $ID.Onil.sorted.markdup.bam -q 0 -b Onil.V1R.all.pos  > Onil_ALL_bam_coverages_stats/$ID.V1R.unmasked.raw.depth

rm Onil_ALL_bam_coverages_stats/$ID.V1R.unmasked.mean.depth
for line in `cat Oreochromis_nilotius_V1R_exons_coordinates.tsv` ; do
	scaffold=`echo "$line" | cut -f6`
	start=`echo "$line" | cut -f7`
	end=`echo "$line" | cut -f8`
	mean_coverage=`awk -v myvar_scaff="$scaffold" '{ if ($1 == myvar_scaff) { print } }' Onil_ALL_bam_coverages_stats/$ID.V1R.unmasked.raw.depth | awk -v myvar_start="$start" '{ if ($2 >= myvar_start) { print } }' | awk -v myvar_end="$end" '{ if ($2 <= myvar_end) { print } }'  | awk '{D+=int($3)} END{print D/(NR*1.0);}'`
	echo "$line	$ID	$mean_coverage" >> Onil_ALL_bam_coverages_stats/$ID.V1R.unmasked.mean.depth
done

#Get the coverage on T1R genes
samtools depth -a $ID.Onil.sorted.markdup.bam -q 0 -b Onil.T1R.all.pos > Onil_ALL_bam_coverages_stats/$ID.T1R.unmasked.raw.depth

rm Onil_ALL_bam_coverages_stats/$ID.T1R.unmasked.mean.depth
for line in `cat Oreochromis_nilotius_T1R_exons_coordinates.tsv` ; do
	scaffold=`echo "$line" | cut -f6`
	start=`echo "$line" | cut -f7`
	end=`echo "$line" | cut -f8`
	mean_coverage=`awk -v myvar_scaff="$scaffold" '{ if ($1 == myvar_scaff) { print } }' Onil_ALL_bam_coverages_stats/$ID.T1R.unmasked.raw.depth | awk -v myvar_start="$start" '{ if ($2 >= myvar_start) { print } }' | awk -v myvar_end="$end" '{ if ($2 <= myvar_end) { print } }'  | awk '{D+=int($3)} END{print D/(NR*1.0);}'`
	echo "$line	$ID	$mean_coverage" >> Onil_ALL_bam_coverages_stats/$ID.T1R.unmasked.mean.depth
done

#Get the coverage on T2R genes
samtools depth -a $ID.Onil.sorted.markdup.bam -q 0 -b Onil.T2R.all.pos > Onil_ALL_bam_coverages_stats/$ID.T2R.unmasked.raw.depth

rm Onil_ALL_bam_coverages_stats/$ID.T2R.unmasked.mean.depth
for line in `cat Oreochromis_nilotius_T2R_exons_coordinates.tsv` ; do
	scaffold=`echo "$line" | cut -f6`
	start=`echo "$line" | cut -f7`
	end=`echo "$line" | cut -f8`
	mean_coverage=`awk -v myvar_scaff="$scaffold" '{ if ($1 == myvar_scaff) { print } }' Onil_ALL_bam_coverages_stats/$ID.T2R.unmasked.raw.depth | awk -v myvar_start="$start" '{ if ($2 >= myvar_start) { print } }' | awk -v myvar_end="$end" '{ if ($2 <= myvar_end) { print } }'  | awk '{D+=int($3)} END{print D/(NR*1.0);}'`
	echo "$line	$ID	$mean_coverage" >> Onil_ALL_bam_coverages_stats/$ID.T2R.unmasked.mean.depth
done


#Get the coverage on all BUSCO genes
samtools depth -a $ID.Onil.sorted.markdup.bam -q 0 -b Oreochromis_BUSCO_genes_cds_locations_c123.tsv | awk '{D+=int($3)} END{print D/(NR*1.0);}' > Onil_ALL_bam_coverages_stats/$ID.AllBUSCO.unmasked.depth

#Get the coverage on conserved BUSCO genes
samtools depth -a $ID.Onil.sorted.markdup.bam -q 0 -b Oreochromis_Conserved_BUSCO_genes_cds_locations_c123.tsv | awk '{D+=int($3)} END{print D/(NR*1.0);}' > Onil_ALL_bam_coverages_stats/$ID.ConservedBUSCO.unmasked.depth

#Get the GW coverage
samtools depth -a $ID.Onil.sorted.markdup.bam -q 0 | awk '{D+=int($3)} END{print D/(NR*1.0);}' > Onil_ALL_bam_coverages_stats/$ID.GW.unmasked.depth


rm $ID.Onil.sorted.bam ; rm $ID.Onil.sorted.markdup.bam ; rm $ID.Onil.sorted.bam.bai
rm $ID.OnilMasked.sorted.bam ; rm $ID.OnilMasked.sorted.markdup.bam ; rm $ID.OnilMasked.sorted.bam.bai
