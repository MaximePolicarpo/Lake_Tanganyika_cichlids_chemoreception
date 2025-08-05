###==================================================================================================================================
###==================================================================================================================================
###==================================================================================================================================
########################################## Lets start with reads cleaning ! #########################################################
###==================================================================================================================================
###==================================================================================================================================
###==================================================================================================================================


mkdir FastQC_Results ; mkdir FastP_results

for folder in `cat RNA_seq_folders.txt` ; do sbatch --qos=6hours -c 4 --mem=10G Fastp_RNA_reads.sh $folder ; done



============== Fastp_RNA_reads.sh  ======================================================================
#!/bin/bash


#SBATCH --job-name=Fastp_RNA_reads   # Job name


eval "$(conda shell.bash hook)"
conda activate olfactory


LMOD_DISABLE_SAME_NAME_AUTOSWAP=no

module purge
module load fastp 
module load FastQC


folder=$1 


R1_reads=`ls -l $folder | grep "fastq.gz" | sed 's/.* //g' | grep "_R1_"`
R2_reads=`ls -l $folder | grep "fastq.gz" | sed 's/.* //g' | grep "_R2_"`

R1_reads_cleaned=`echo "$R1_reads" | sed 's/.fastq.gz/.fastpcleaned.fastq.gz/g'`
R2_reads_cleaned=`echo "$R2_reads" | sed 's/.fastq.gz/.fastpcleaned.fastq.gz/g'`

report_file=`echo "$R1_reads" | sed 's/.fastq.gz//g'`


cp $folder/$R1_reads ./
cp $folder/$R2_reads ./

fastp -i $R1_reads -I $R2_reads -o $R1_reads_cleaned -O $R2_reads_cleaned -h $report_file.html -j $report_file.json

fastqc -t 4 --noextract -o FastQC_Results $R1_reads_cleaned
fastqc -t 4 --noextract -o FastQC_Results $R2_reads_cleaned
fastqc -t 4 --noextract -o FastQC_Results $R1_reads
fastqc -t 4 --noextract -o FastQC_Results $R2_reads

rm $R1_reads ; rm $R2_reads
mv $report_file.html FastP_results/
mv $report_file.json FastP_results/

==============================================================================================================================



###==================================================================================================================================
###==================================================================================================================================
###==================================================================================================================================
########################################## Reads mapping against O. niloticus ! #####################################################
###==================================================================================================================================
###==================================================================================================================================
###==================================================================================================================================


#first create a index of the genome

module load STAR #https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf
sbatch --qos=1day -c 8 --mem=50G --job-name=STAR_genome_generate --wrap="module load STAR ; STAR --runThreadN 8 --runMode genomeGenerate --genomeDir /scicore/home/salzburg/polica0000/Cichlids_Genomes/Oreochromis_niloticus/ --genomeFastaFiles GCF_001858045.2_O_niloticus_UMD_NMBU_genomic.fna --sjdbGTFfile GCF_001858045.2_O_niloticus_UMD_NMBU_genomic.gtf --genomeSAindexNbases 13"


mkdir STAR_vs_GCF_001858045.2
sbatch --qos=1day -c 30 --mem=200G STAR_align.sh


========================================= STAR_align.sh ===========================================================================================================================

#!/bin/bash


#SBATCH --job-name=STAR_alignment   # Job name


eval "$(conda shell.bash hook)"
conda activate olfactory

module load STAR


for folder in `cat RNA_seq_folders.txt` ; do

	R1_reads=`ls -l $folder | grep "fastq.gz" | sed 's/.* //g' | grep "_R1_"`
	R2_reads=`ls -l $folder | grep "fastq.gz" | sed 's/.* //g' | grep "_R2_"`

	R1_reads_cleaned=`echo "$R1_reads" | sed 's/.fastq.gz/.fastpcleaned.fastq.gz/g'`
	R2_reads_cleaned=`echo "$R2_reads" | sed 's/.fastq.gz/.fastpcleaned.fastq.gz/g'`


	sample_name=`echo "$R1_reads_cleaned" | cut -f4 -d "_"`

	STAR --runMode alignReads --runThreadN 30 --genomeDir /scicore/home/salzburg/polica0000/Cichlids_Genomes/Oreochromis_niloticus/ --readFilesIn $R1_reads_cleaned $R2_reads_cleaned --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --readFilesCommand zcat --outFileNamePrefix $sample_name --outSAMunmapped Within --outSAMmultNmax 1

done


mv *.sortedByCoord.out.bam STAR_vs_GCF_001858045.2/
mv *.out.tab STAR_vs_GCF_001858045.2/
mv *Log.progress.out STAR_vs_GCF_001858045.2/
mv *Log.out STAR_vs_GCF_001858045.2/
mv *.final.out STAR_vs_GCF_001858045.2/


======================================================================================================================================================================================================================================================


#Create an index for bam files


for file in *.bam ; do sbatch --qos=6hours -c 8 --mem=8G samtools_index.sh $file ; done



====================== samtools_index.sh ==========================================================================================================================================================


#!/bin/bash


#SBATCH --job-name=samtools_index   # Job name


eval "$(conda shell.bash hook)"
conda activate olfactory

module load SAMtools

file=$1

samtools index $file


======================================================================================================================================================================================================================================================


###==================================================================================================================================
###==================================================================================================================================
###==================================================================================================================================
############ Count the reads on each annotated genes + manually annotated olfactory receptors of O. niloticus #######################
###==================================================================================================================================
###==================================================================================================================================
###==================================================================================================================================

mkdir HTSeq_Count_Results

for i in *.bam ; do sbatch --qos=6hours -c 8 --mem=10G HTSeq_count.sh $i ; done


============================ HTSeq_count.sh ================================================================================================================

#!/bin/bash


#SBATCH --job-name=HTSeq_count   # Job name


eval "$(conda shell.bash hook)"
conda activate HTSEQ_env

bam_file=$1
sample_name=`echo "$bam_file" | sed 's/OEAligned.*//g'`

htseq-count -f bam -r pos -m union -s reverse -t exon $bam_file Olfactory_receptors.gff3 --idattr=ID > HTSeq_Count_Results/$sample_name.Olfactory_receptors.count

htseq-count -f bam -r pos -m union -s reverse -t exon $bam_file GCF_001858045.2_O_niloticus_UMD_NMBU_genomic.gff.LI.corrected --idattr=ID > HTSeq_Count_Results/$sample_name.ALL_genes.count


================================================================================================================================================================================================================================


for i in *.count ; do head -n -5 $i > $i.parsed ; done
 
for i in *.parsed ; do sample=`echo "$i" | sed 's/\..*//g'` ; sed -i "s/^/$sample,/g" $i ; done

cat *.Olfactory_receptors.count.parsed > Read_count_all_samples.olfactory_receptors.csv
cat *.ALL_genes.count.parsed > Read_count_all_samples.all_genes.csv

sed -i 's/OR-exon/OR,OR-exon/g' Read_count_all_samples.olfactory_receptors.csv
sed -i 's/TAAR-exon/TAAR,TAAR-exon/g' Read_count_all_samples.olfactory_receptors.csv
sed -i 's/V2R-exon/V2R,V2R-exon/g' Read_count_all_samples.olfactory_receptors.csv
sed -i 's/V1R-exon/V1R,V1R-exon/g' Read_count_all_samples.olfactory_receptors.csv

sed -i 's/-[0-9]*	/,/g' Read_count_all_samples.olfactory_receptors.csv
sed -i 's/-[0-9]*	/,/g' Read_count_all_samples.all_genes.csv


# Add the name of the OGG for each gene

cut -f3 -d "," Read_count_all_samples.olfactory_receptors.csv  | sed 's/V2R-exon-//g' | sed 's/V1R-exon-//g'  | sed 's/OR-exon-//g' | sed 's/TAAR-exon-//g' > col3_onil_gene_name.txt



sed -i 's/OR-exon-//g' Read_count_all_samples.olfactory_receptors.csv
sed -i 's/TAAR-exon-//g' Read_count_all_samples.olfactory_receptors.csv
sed -i 's/V1R-exon-//g' Read_count_all_samples.olfactory_receptors.csv
sed -i 's/V2R-exon-//g' Read_count_all_samples.olfactory_receptors.csv


rm Read_count_all_samples.olfactory_receptors.OGG.csv
IFS=$'\n'
for line in `cat Read_count_all_samples.olfactory_receptors.csv` ; do
	gene=`echo "$line" | cut -f3 -d ","`
	gene_type=`grep "$gene" /scicore/home/salzburg/polica0000/Phylogeny_analysis_Olfactory_receptor/ALL_genes_classification/ALL_Chemoreceptors_classification.csv | cut -f1 -d ","`
	subfamily=`grep "$gene" /scicore/home/salzburg/polica0000/Phylogeny_analysis_Olfactory_receptor/ALL_genes_classification/ALL_Chemoreceptors_classification.csv | cut -f3 -d ","`
	clade=`grep "$gene" /scicore/home/salzburg/polica0000/Phylogeny_analysis_Olfactory_receptor/ALL_genes_classification/ALL_Chemoreceptors_classification.csv | cut -f4 -d ","`
	OGG=`grep "$gene" /scicore/home/salzburg/polica0000/Phylogeny_analysis_Olfactory_receptor/ALL_genes_classification/ALL_Chemoreceptors_classification.csv | cut -f5 -d ","`

	echo "$line,$gene_type,$clade,$subfamily,$OGG" >> Read_count_all_samples.olfactory_receptors.OGG.csv

done


grep "Oreochromis" ALL_TAAR.fa.fai | cut -f1,2 | sed 's/Oreochromis_niloticus---//g' > gene_length.tsv
grep "Oreochromis" ALL_OR.fa.fai | cut -f1,2 | sed 's/Oreochromis_niloticus---//g' >> gene_length.tsv
grep "Oreochromis" Cichlids_V1R_Complete.fa.fai | cut -f1,2 | sed 's/Oreochromis_niloticus---//g' >> gene_length.tsv
grep "Oreochromis" ALL_V2R.fa.fai | cut -f1,2 | sed 's/Oreochromis_niloticus---//g' >> gene_length.tsv


