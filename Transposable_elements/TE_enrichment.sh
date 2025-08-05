#####################################################################################################################################################################
#####################################################################################################################################################################
#####################################################################################################################################################################
#####################################################################################################################################################################
############################################################# V2R ##############################################################################################################
#####################################################################################################################################################################
#####################################################################################################################################################################
#####################################################################################################################################################################




#import V2R genes bed file and generate clusters manually
cd Bathybates_minor/ ; rm * ; grep ">Bathy" ALL_V2R.fa | sed 's/>//g' | sed 's/Bathybates_minor---//g' | sed 's/---.*//g' | sed 's/-/,/g' | sort -k1,1 -k2n -t "," > V2R_genes.bed ; cd ../
cd Cunningtonia_longiventralis/ ; rm * ; grep ">Cunningtonia" ALL_V2R.fa | sed 's/>//g' | sed 's/Cunningtonia_longiventralis---//g' | sed 's/---.*//g' | sed 's/-/,/g' | sort -k1,1 -k2n -t "," > V2R_genes.bed ; cd ../
cd Cyphotilapia_frontosa/ ; rm * ; grep ">Cyphotilapia" ALL_V2R.fa | sed 's/>//g' | sed 's/Cyphotilapia_frontosa---//g' | sed 's/---.*//g' | sed 's/-/,/g' | sort -k1,1 -k2n -t "," > V2R_genes.bed ; cd ../
cd Cyprichromis_leptosoma/ ; rm * ; grep ">Cyprichromis" ALL_V2R.fa | sed 's/>//g' | sed 's/Cyprichromis_leptosoma---//g' | sed 's/---.*//g' | sed 's/-/,/g' | sort -k1,1 -k2n -t "," > V2R_genes.bed ; cd ../
cd Neolamprologus_multifasciatus/ ; rm * ; grep ">Neolamprologus" ALL_V2R.fa | sed 's/>//g' | sed 's/Neolamprologus_multifasciatus---//g' | sed 's/---.*//g' | sed 's/-/,/g' | sort -k1,1 -k2n -t "," > V2R_genes.bed ; cd ../
cd Simochromis_diagramma/ ; rm * ; grep ">Simochromis" ALL_V2R.fa | sed 's/>//g' | sed 's/Simochromis_diagramma---//g' | sed 's/---.*//g' | sed 's/-/,/g' | sort -k1,1 -k2n -t "," > V2R_genes.bed ; cd ../
cd Oreochromis_niloticus/ ; rm * ; grep ">Oreochromis_niloticus" ALL_V2R.fa | sed 's/>//g' | sed 's/Oreochromis_niloticus---//g' | sed 's/---.*//g' | sed 's/-/,/g' | sort -k1,1 -k2n -t "," > V2R_genes.bed ; cd ../


cd Bathybates_minor/ ; sbatch --job-name=V2R_Cluster --qos=1day -c 8 --mem=15G ClustEnrichTE.sh fBatMin1.PB.asm1.purge1.scaff1.polish2.decontaminated.fa.out /scicore/home/salzburg/polica0000/Cichlids_Genomes/PacBio_Genomes_Cichlids/Bathybates_minor/fBatMin1.PB.asm1.purge1.scaff1.polish2.decontaminated.fa.fai V2R_genes.bed fBatMin1.PB.asm1.purge1.scaff1.polish2.decontaminated.fa.align cds_coordinates_V2R_genes.Bathybates_minor.c123.csv ; cd ../
cd Cunningtonia_longiventralis/ ; sbatch --job-name=V2R_Cluster --qos=1day -c 8 --mem=15G ClustEnrichTE.sh fCunLon1.PB.asm1.purge1.fa.out /scicore/home/salzburg/polica0000/Cichlids_Genomes/PacBio_Genomes_Cichlids/Cunningtonia_longiventralis/fCunLon1.PB.asm1.purge1.fa.fai V2R_genes.bed fCunLon1.PB.asm1.purge1.fa.align cds_coordinates_V2R_genes.Cunningtonia_longiventralis.c123.csv ; cd ../
cd Cyphotilapia_frontosa/ ; sbatch --job-name=V2R_Cluster --qos=1day -c 8 --mem=15G ClustEnrichTE.sh fCypFro1.PB.asm1.purge1.scaff1.polish2.decontaminated.fa.out /scicore/home/salzburg/polica0000/Cichlids_Genomes/PacBio_Genomes_Cichlids/Cyphotilapia_frontosa/fCypFro1.PB.asm1.purge1.scaff1.polish2.decontaminated.fa.fai V2R_genes.bed fCypFro1.PB.asm1.purge1.scaff1.polish2.decontaminated.fa.align cds_coordinates_V2R_genes.Cyphotilapia_frontosa.c123.csv ; cd ../
cd Cyprichromis_leptosoma/ ; sbatch --job-name=V2R_Cluster --qos=1day -c 8 --mem=15G ClustEnrichTE.sh fCypLep1.PB.asm1.purge1.scaff1.polish2.decontaminated.fa.out /scicore/home/salzburg/polica0000/Cichlids_Genomes/PacBio_Genomes_Cichlids/Cyprichromis_leptosoma/fCypLep1.PB.asm1.purge1.scaff1.polish2.decontaminated.fa.fai V2R_genes.bed fCypLep1.PB.asm1.purge1.scaff1.polish2.decontaminated.fa.align cds_coordinates_V2R_genes.Cyprichromis_leptosoma.c123.csv ; cd ../
cd Neolamprologus_multifasciatus/ ; sbatch --job-name=V2R_Cluster --qos=1week -c 8 --mem=50G ClustEnrichTE.sh fNeoMul12_final.fa.out /scicore/home/salzburg/polica0000/Cichlids_Genomes/PacBio_Genomes_Cichlids/Neolamprologus_multifasciatus/fNeoMul12_final.fa.fai V2R_genes.bed fNeoMul12_final.fa.align cds_coordinates_V2R_genes.Neolamprologus_multifasciatus.c123.csv ; cd ../
cd Simochromis_diagramma/ ; sbatch --job-name=V2R_Cluster --qos=1day -c 8 --mem=15G ClustEnrichTE.sh GCF_900408965.1_fSimDia1.1_genomic.fna.out /scicore/home/salzburg/polica0000/Cichlids_Genomes/PacBio_Genomes_Cichlids/Simochromis_diagramma/GCF_900408965.1_fSimDia1.1_genomic.fna.fai V2R_genes.bed GCF_900408965.1_fSimDia1.1_genomic.fna.align cds_coordinates_V2R_genes.Simochromis_diagramma.c123.csv ; cd ../																																														 
cd Oreochromis_niloticus/ ; sbatch --job-name=V2R_Cluster --qos=1week -c 8 --mem=100G ClustEnrichTE.sh GCF_001858045.2_O_niloticus_UMD_NMBU_genomic.fna.out /scicore/home/salzburg/polica0000/Cichlids_Genomes/Oreochromis_niloticus/GCF_001858045.2_O_niloticus_UMD_NMBU_genomic.fna.fai V2R_genes.bed GCF_001858045.2_O_niloticus_UMD_NMBU_genomic.fna.align cds_coordinates_V2R_genes.Oreochromis_niloticus.c123.csv ; cd ../




#####################################################################################################################################################################
#####################################################################################################################################################################
#####################################################################################################################################################################
#####################################################################################################################################################################
############################################################# TAAR ##############################################################################################################
#####################################################################################################################################################################
#####################################################################################################################################################################
#####################################################################################################################################################################


mkdir Bathybates_minor
mkdir Cunningtonia_longiventralis
mkdir Cyphotilapia_frontosa
mkdir Cyprichromis_leptosoma
mkdir Neolamprologus_multifasciatus
mkdir Simochromis_diagramma
mkdir Oreochromis_niloticus

#import TAAR genes bed file and generate clusters manually
cd Bathybates_minor/ ; rm * ; grep ">Bathy" ALL_TAAR.fa | sed 's/>//g' | sed 's/Bathybates_minor---//g' | sed 's/---.*//g' | sed 's/-/,/g' | sort -k1,1 -k2n -t "," > TAAR_genes.bed ; cd ../
cd Cunningtonia_longiventralis/ ; rm * ; grep ">Cunningtonia" ALL_TAAR.fa | sed 's/>//g' | sed 's/Cunningtonia_longiventralis---//g' | sed 's/---.*//g' | sed 's/-/,/g' | sort -k1,1 -k2n -t "," > TAAR_genes.bed ; cd ../
cd Cyphotilapia_frontosa/ ; rm * ; grep ">Cyphotilapia" ALL_TAAR.fa | sed 's/>//g' | sed 's/Cyphotilapia_frontosa---//g' | sed 's/---.*//g' | sed 's/-/,/g' | sort -k1,1 -k2n -t "," > TAAR_genes.bed ; cd ../
cd Cyprichromis_leptosoma/ ; rm * ; grep ">Cyprichromis" ALL_TAAR.fa | sed 's/>//g' | sed 's/Cyprichromis_leptosoma---//g' | sed 's/---.*//g' | sed 's/-/,/g' | sort -k1,1 -k2n -t "," > TAAR_genes.bed ; cd ../
cd Neolamprologus_multifasciatus/ ; rm * ; grep ">Neolamprologus" ALL_TAAR.fa | sed 's/>//g' | sed 's/Neolamprologus_multifasciatus---//g' | sed 's/---.*//g' | sed 's/-/,/g' | sort -k1,1 -k2n -t "," > TAAR_genes.bed ; cd ../
cd Simochromis_diagramma/ ; rm * ; grep ">Simochromis" ALL_TAAR.fa | sed 's/>//g' | sed 's/Simochromis_diagramma---//g' | sed 's/---.*//g' | sed 's/-/,/g' | sort -k1,1 -k2n -t "," > TAAR_genes.bed ; cd ../
cd Oreochromis_niloticus/ ; rm * ; grep ">Oreochromis_niloticus" ALL_TAAR.fa | sed 's/>//g' | sed 's/Oreochromis_niloticus---//g' | sed 's/---.*//g' | sed 's/-/,/g' | sort -k1,1 -k2n -t "," > TAAR_genes.bed ; cd ../

cd Bathybates_minor/ ; sbatch --job-name=TAAR_Cluster --qos=1day -c 8 --mem-per-cpu=15g ClustEnrichTE.sh fBatMin1.PB.asm1.purge1.scaff1.polish2.decontaminated.fa.out /scicore/home/salzburg/polica0000/Cichlids_Genomes/PacBio_Genomes_Cichlids/Bathybates_minor/fBatMin1.PB.asm1.purge1.scaff1.polish2.decontaminated.fa.fai TAAR_genes.bed fBatMin1.PB.asm1.purge1.scaff1.polish2.decontaminated.fa.align cds_coordinates_TAAR_genes.Bathybates_minor.c123.csv ; cd ../
cd Cunningtonia_longiventralis/ ; sbatch --job-name=TAAR_Cluster --qos=1day -c 8 --mem-per-cpu=15g ClustEnrichTE.sh fCunLon1.PB.asm1.purge1.fa.out /scicore/home/salzburg/polica0000/Cichlids_Genomes/PacBio_Genomes_Cichlids/Cunningtonia_longiventralis/fCunLon1.PB.asm1.purge1.fa.fai TAAR_genes.bed fCunLon1.PB.asm1.purge1.fa.align cds_coordinates_TAAR_genes.Cunningtonia_longiventralis.c123.csv ; cd ../
cd Cyphotilapia_frontosa/ ; sbatch --job-name=TAAR_Cluster --qos=1day -c 8 --mem-per-cpu=15g ClustEnrichTE.sh fCypFro1.PB.asm1.purge1.scaff1.polish2.decontaminated.fa.out /scicore/home/salzburg/polica0000/Cichlids_Genomes/PacBio_Genomes_Cichlids/Cyphotilapia_frontosa/fCypFro1.PB.asm1.purge1.scaff1.polish2.decontaminated.fa.fai TAAR_genes.bed fCypFro1.PB.asm1.purge1.scaff1.polish2.decontaminated.fa.align cds_coordinates_TAAR_genes.Cyphotilapia_frontosa.c123.csv ; cd ../
cd Cyprichromis_leptosoma/ ; sbatch --job-name=TAAR_Cluster --qos=1day -c 8 --mem-per-cpu=15g ClustEnrichTE.sh fCypLep1.PB.asm1.purge1.scaff1.polish2.decontaminated.fa.out /scicore/home/salzburg/polica0000/Cichlids_Genomes/PacBio_Genomes_Cichlids/Cyprichromis_leptosoma/fCypLep1.PB.asm1.purge1.scaff1.polish2.decontaminated.fa.fai TAAR_genes.bed fCypLep1.PB.asm1.purge1.scaff1.polish2.decontaminated.fa.align cds_coordinates_TAAR_genes.Cyprichromis_leptosoma.c123.csv ; cd ../
cd Neolamprologus_multifasciatus/ ; sbatch --job-name=TAAR_Cluster --qos=1day -c 8 --mem-per-cpu=15g ClustEnrichTE.sh fNeoMul12_final.fa.out /scicore/home/salzburg/polica0000/Cichlids_Genomes/PacBio_Genomes_Cichlids/Neolamprologus_multifasciatus/fNeoMul12_final.fa.fai TAAR_genes.bed fNeoMul12_final.fa.align cds_coordinates_TAAR_genes.Neolamprologus_multifasciatus.c123.csv ; cd ../
cd Simochromis_diagramma/ ; sbatch --job-name=TAAR_Cluster --qos=1day -c 8 --mem-per-cpu=15g ClustEnrichTE.sh GCF_900408965.1_fSimDia1.1_genomic.fna.out /scicore/home/salzburg/polica0000/Cichlids_Genomes/PacBio_Genomes_Cichlids/Simochromis_diagramma/GCF_900408965.1_fSimDia1.1_genomic.fna.fai TAAR_genes.bed GCF_900408965.1_fSimDia1.1_genomic.fna.align cds_coordinates_TAAR_genes.Simochromis_diagramma.c123.csv ; cd ../
cd Oreochromis_niloticus/ ; sbatch --job-name=TAAR_Cluster --qos=1day -c 8 --mem-per-cpu=15g ClustEnrichTE.sh GCF_001858045.2_O_niloticus_UMD_NMBU_genomic.fna.out /scicore/home/salzburg/polica0000/Cichlids_Genomes/Oreochromis_niloticus/GCF_001858045.2_O_niloticus_UMD_NMBU_genomic.fna.fai TAAR_genes.bed GCF_001858045.2_O_niloticus_UMD_NMBU_genomic.fna.align cds_coordinates_TAAR_genes.Oreochromis_niloticus.c123.csv ; cd ../

#####################################################################################################################################################################
#####################################################################################################################################################################
#####################################################################################################################################################################
#####################################################################################################################################################################
############################################################# OR ##############################################################################################################
#####################################################################################################################################################################
#####################################################################################################################################################################
#####################################################################################################################################################################


rm -r Bathybates_minor
rm -r Cunningtonia_longiventralis
rm -r Cyphotilapia_frontosa
rm -r Cyprichromis_leptosoma
rm -r Neolamprologus_multifasciatus
rm -r Simochromis_diagramma
rm -r Oreochromis_niloticus

#import OR genes bed file and generate clusters manually
cd Bathybates_minor/ ; rm * ; grep ">Bathy" ALL_OR.fa | sed 's/>//g' | sed 's/Bathybates_minor---//g' | sed 's/---.*//g' | sed 's/-/,/g' | sort -k1,1 -k2n -t "," > OR_genes.bed ; cd ../
cd Cunningtonia_longiventralis/ ; rm * ; grep ">Cunningtonia" ALL_OR.fa | sed 's/>//g' | sed 's/Cunningtonia_longiventralis---//g' | sed 's/---.*//g' | sed 's/-/,/g' | sort -k1,1 -k2n -t "," > OR_genes.bed ; cd ../
cd Cyphotilapia_frontosa/ ; rm * ; grep ">Cyphotilapia" ALL_OR.fa | sed 's/>//g' | sed 's/Cyphotilapia_frontosa---//g' | sed 's/---.*//g' | sed 's/-/,/g' | sort -k1,1 -k2n -t "," > OR_genes.bed ; cd ../
cd Cyprichromis_leptosoma/ ; rm * ; grep ">Cyprichromis" ALL_OR.fa | sed 's/>//g' | sed 's/Cyprichromis_leptosoma---//g' | sed 's/---.*//g' | sed 's/-/,/g' | sort -k1,1 -k2n -t "," > OR_genes.bed ; cd ../
cd Neolamprologus_multifasciatus/ ; rm * ; grep ">Neolamprologus" ALL_OR.fa | sed 's/>//g' | sed 's/Neolamprologus_multifasciatus---//g' | sed 's/---.*//g' | sed 's/-/,/g' | sort -k1,1 -k2n -t "," > OR_genes.bed ; cd ../
cd Simochromis_diagramma/ ; rm * ; grep ">Simochromis" ALL_OR.fa | sed 's/>//g' | sed 's/Simochromis_diagramma---//g' | sed 's/---.*//g' | sed 's/-/,/g' | sort -k1,1 -k2n -t "," > OR_genes.bed ; cd ../
cd Oreochromis_niloticus/ ; rm * ; grep ">Oreochromis_niloticus" ALL_OR.fa | sed 's/>//g' | sed 's/Oreochromis_niloticus---//g' | sed 's/---.*//g' | sed 's/-/,/g' | sort -k1,1 -k2n -t "," > OR_genes.bed ; cd ../


cd Bathybates_minor/ ; sbatch --job-name=OR_Cluster --qos=1day -c 8 --mem-per-cpu=15g ClustEnrichTE.sh fBatMin1.PB.asm1.purge1.scaff1.polish2.decontaminated.fa.out /scicore/home/salzburg/polica0000/Cichlids_Genomes/PacBio_Genomes_Cichlids/Bathybates_minor/fBatMin1.PB.asm1.purge1.scaff1.polish2.decontaminated.fa.fai OR_genes.bed fBatMin1.PB.asm1.purge1.scaff1.polish2.decontaminated.fa.align cds_coordinates_OR_genes.Bathybates_minor.c123.csv ; cd ../
cd Cunningtonia_longiventralis/ ; sbatch --job-name=OR_Cluster --qos=1day -c 8 --mem-per-cpu=15g ClustEnrichTE.sh fCunLon1.PB.asm1.purge1.fa.out /scicore/home/salzburg/polica0000/Cichlids_Genomes/PacBio_Genomes_Cichlids/Cunningtonia_longiventralis/fCunLon1.PB.asm1.purge1.fa.fai OR_genes.bed fCunLon1.PB.asm1.purge1.fa.align cds_coordinates_OR_genes.Cunningtonia_longiventralis.c123.csv ; cd ../
cd Cyphotilapia_frontosa/ ; sbatch --job-name=OR_Cluster --qos=1day -c 8 --mem-per-cpu=15g ClustEnrichTE.sh fCypFro1.PB.asm1.purge1.scaff1.polish2.decontaminated.fa.out /scicore/home/salzburg/polica0000/Cichlids_Genomes/PacBio_Genomes_Cichlids/Cyphotilapia_frontosa/fCypFro1.PB.asm1.purge1.scaff1.polish2.decontaminated.fa.fai OR_genes.bed fCypFro1.PB.asm1.purge1.scaff1.polish2.decontaminated.fa.align cds_coordinates_OR_genes.Cyphotilapia_frontosa.c123.csv ; cd ../
cd Cyprichromis_leptosoma/ ; sbatch --job-name=OR_Cluster --qos=1day -c 8 --mem-per-cpu=15g ClustEnrichTE.sh fCypLep1.PB.asm1.purge1.scaff1.polish2.decontaminated.fa.out /scicore/home/salzburg/polica0000/Cichlids_Genomes/PacBio_Genomes_Cichlids/Cyprichromis_leptosoma/fCypLep1.PB.asm1.purge1.scaff1.polish2.decontaminated.fa.fai OR_genes.bed fCypLep1.PB.asm1.purge1.scaff1.polish2.decontaminated.fa.align cds_coordinates_OR_genes.Cyprichromis_leptosoma.c123.csv ; cd ../
cd Neolamprologus_multifasciatus/ ; sbatch --job-name=OR_Cluster --qos=1week -c 8 --mem-per-cpu=15g ClustEnrichTE.sh fNeoMul12_final.fa.out /scicore/home/salzburg/polica0000/Cichlids_Genomes/PacBio_Genomes_Cichlids/Neolamprologus_multifasciatus/fNeoMul12_final.fa.fai OR_genes.bed fNeoMul12_final.fa.align cds_coordinates_OR_genes.Neolamprologus_multifasciatus.c123.csv ; cd ../
cd Simochromis_diagramma/ ; sbatch --job-name=OR_Cluster --qos=1day -c 8 --mem-per-cpu=15g ClustEnrichTE.sh GCF_900408965.1_fSimDia1.1_genomic.fna.out /scicore/home/salzburg/polica0000/Cichlids_Genomes/PacBio_Genomes_Cichlids/Simochromis_diagramma/GCF_900408965.1_fSimDia1.1_genomic.fna.fai OR_genes.bed GCF_900408965.1_fSimDia1.1_genomic.fna.align cds_coordinates_OR_genes.Simochromis_diagramma.c123.csv ; cd ../																																														 
cd Oreochromis_niloticus/ ; sbatch --job-name=OR_Cluster --qos=1week -c 8 --mem-per-cpu=25g ClustEnrichTE.sh GCF_001858045.2_O_niloticus_UMD_NMBU_genomic.fna.out /scicore/home/salzburg/polica0000/Cichlids_Genomes/Oreochromis_niloticus/GCF_001858045.2_O_niloticus_UMD_NMBU_genomic.fna.fai OR_genes.bed GCF_001858045.2_O_niloticus_UMD_NMBU_genomic.fna.align cds_coordinates_OR_genes.Oreochromis_niloticus.c123.csv ; cd ../




#####################################################################################################################################################################
#####################################################################################################################################################################
#####################################################################################################################################################################
#####################################################################################################################################################################
############################################################# V1R ##############################################################################################################
#####################################################################################################################################################################
#####################################################################################################################################################################
#####################################################################################################################################################################


mkdir Bathybates_minor
mkdir Cunningtonia_longiventralis
mkdir Cyphotilapia_frontosa
mkdir Cyprichromis_leptosoma
mkdir Neolamprologus_multifasciatus
mkdir Simochromis_diagramma
mkdir Oreochromis_niloticus

#import V1R genes bed file and generate clusters manually
cd Bathybates_minor/ ; rm * ; grep ">Bathy" ALL_V1R.prot | sed 's/>//g' | sed 's/Bathybates_minor---//g' | sed 's/---.*//g' | sed 's/-/,/g' | sort -k1,1 -k2n -t "," > V1R_genes.bed ; cd ../
cd Cunningtonia_longiventralis/ ; rm * ; grep ">Cunningtonia" ALL_V1R.prot | sed 's/>//g' | sed 's/Cunningtonia_longiventralis---//g' | sed 's/---.*//g' | sed 's/-/,/g' | sort -k1,1 -k2n -t "," > V1R_genes.bed ; cd ../
cd Cyphotilapia_frontosa/ ; rm * ; grep ">Cyphotilapia" ALL_V1R.prot | sed 's/>//g' | sed 's/Cyphotilapia_frontosa---//g' | sed 's/---.*//g' | sed 's/-/,/g' | sort -k1,1 -k2n -t "," > V1R_genes.bed ; cd ../
cd Cyprichromis_leptosoma/ ; rm * ; grep ">Cyprichromis" ALL_V1R.prot | sed 's/>//g' | sed 's/Cyprichromis_leptosoma---//g' | sed 's/---.*//g' | sed 's/-/,/g' | sort -k1,1 -k2n -t "," > V1R_genes.bed ; cd ../
cd Neolamprologus_multifasciatus/ ; rm * ; grep ">Neolamprologus" ALL_V1R.prot | sed 's/>//g' | sed 's/Neolamprologus_multifasciatus---//g' | sed 's/---.*//g' | sed 's/-/,/g' | sort -k1,1 -k2n -t "," > V1R_genes.bed ; cd ../
cd Simochromis_diagramma/ ; rm * ; grep ">Simochromis" ALL_V1R.prot | sed 's/>//g' | sed 's/Simochromis_diagramma---//g' | sed 's/---.*//g' | sed 's/-/,/g' | sort -k1,1 -k2n -t "," > V1R_genes.bed ; cd ../
cd Oreochromis_niloticus/ ; rm * ; grep ">Oreochromis_niloticus" ALL_V1R.prot | sed 's/>//g' | sed 's/Oreochromis_niloticus---//g' | sed 's/---.*//g' | sed 's/-/,/g' | sort -k1,1 -k2n -t "," > V1R_genes.bed ; cd ../

cd Bathybates_minor/ ; sbatch --job-name=V1R_Cluster --qos=1day -c 8 --mem-per-cpu=15g ClustEnrichTE.sh fBatMin1.PB.asm1.purge1.scaff1.polish2.decontaminated.fa.out /scicore/home/salzburg/polica0000/Cichlids_Genomes/PacBio_Genomes_Cichlids/Bathybates_minor/fBatMin1.PB.asm1.purge1.scaff1.polish2.decontaminated.fa.fai V1R_genes.bed fBatMin1.PB.asm1.purge1.scaff1.polish2.decontaminated.fa.align cds_coordinates_V1R_genes.Bathybates_minor.c123.csv ; cd ../
cd Cunningtonia_longiventralis/ ; sbatch --job-name=V1R_Cluster --qos=1day -c 8 --mem-per-cpu=15g ClustEnrichTE.sh fCunLon1.PB.asm1.purge1.fa.out /scicore/home/salzburg/polica0000/Cichlids_Genomes/PacBio_Genomes_Cichlids/Cunningtonia_longiventralis/fCunLon1.PB.asm1.purge1.fa.fai V1R_genes.bed fCunLon1.PB.asm1.purge1.fa.align cds_coordinates_V1R_genes.Cunningtonia_longiventralis.c123.csv ; cd ../
cd Cyphotilapia_frontosa/ ; sbatch --job-name=V1R_Cluster --qos=1day -c 8 --mem-per-cpu=15g ClustEnrichTE.sh fCypFro1.PB.asm1.purge1.scaff1.polish2.decontaminated.fa.out /scicore/home/salzburg/polica0000/Cichlids_Genomes/PacBio_Genomes_Cichlids/Cyphotilapia_frontosa/fCypFro1.PB.asm1.purge1.scaff1.polish2.decontaminated.fa.fai V1R_genes.bed fCypFro1.PB.asm1.purge1.scaff1.polish2.decontaminated.fa.align cds_coordinates_V1R_genes.Cyphotilapia_frontosa.c123.csv ; cd ../
cd Cyprichromis_leptosoma/ ; sbatch --job-name=V1R_Cluster --qos=1day -c 8 --mem-per-cpu=15g ClustEnrichTE.sh fCypLep1.PB.asm1.purge1.scaff1.polish2.decontaminated.fa.out /scicore/home/salzburg/polica0000/Cichlids_Genomes/PacBio_Genomes_Cichlids/Cyprichromis_leptosoma/fCypLep1.PB.asm1.purge1.scaff1.polish2.decontaminated.fa.fai V1R_genes.bed fCypLep1.PB.asm1.purge1.scaff1.polish2.decontaminated.fa.align cds_coordinates_V1R_genes.Cyprichromis_leptosoma.c123.csv ; cd ../
cd Neolamprologus_multifasciatus/ ; sbatch --job-name=V1R_Cluster --qos=1day -c 8 --mem-per-cpu=15g ClustEnrichTE.sh fNeoMul12_final.fa.out /scicore/home/salzburg/polica0000/Cichlids_Genomes/PacBio_Genomes_Cichlids/Neolamprologus_multifasciatus/fNeoMul12_final.fa.fai V1R_genes.bed fNeoMul12_final.fa.align cds_coordinates_V1R_genes.Neolamprologus_multifasciatus.c123.csv ; cd ../
cd Simochromis_diagramma/ ; sbatch --job-name=V1R_Cluster --qos=1day -c 8 --mem-per-cpu=15g ClustEnrichTE.sh GCF_900408965.1_fSimDia1.1_genomic.fna.out /scicore/home/salzburg/polica0000/Cichlids_Genomes/PacBio_Genomes_Cichlids/Simochromis_diagramma/GCF_900408965.1_fSimDia1.1_genomic.fna.fai V1R_genes.bed GCF_900408965.1_fSimDia1.1_genomic.fna.align cds_coordinates_V1R_genes.Simochromis_diagramma.c123.csv ; cd ../																																														 
cd Oreochromis_niloticus/ ; sbatch --job-name=V1R_Cluster --qos=1day -c 8 --mem-per-cpu=25g ClustEnrichTE.sh GCF_001858045.2_O_niloticus_UMD_NMBU_genomic.fna.out /scicore/home/salzburg/polica0000/Cichlids_Genomes/Oreochromis_niloticus/GCF_001858045.2_O_niloticus_UMD_NMBU_genomic.fna.fai V1R_genes.bed GCF_001858045.2_O_niloticus_UMD_NMBU_genomic.fna.align cds_coordinates_V1R_genes.Oreochromis_niloticus.c123.csv ; cd ../




#####################################################################################################################################################################
#####################################################################################################################################################################
#####################################################################################################################################################################
#####################################################################################################################################################################
############################################################# T1R ##############################################################################################################
#####################################################################################################################################################################
#####################################################################################################################################################################
#####################################################################################################################################################################



mkdir Bathybates_minor
mkdir Cunningtonia_longiventralis
mkdir Cyphotilapia_frontosa
mkdir Cyprichromis_leptosoma
mkdir Neolamprologus_multifasciatus
mkdir Simochromis_diagramma
mkdir Oreochromis_niloticus

#import OR genes bed file and generate clusters manually
cd Bathybates_minor/ ; rm * ; grep ">Bathy" /scicore/home/salzburg/polica0000/Phylogeny_analysis_Olfactory_receptor/T1R_final/Orthogroup_cluserting/ALL_T1R.fa | sed 's/>//g' | sed 's/Bathybates_minor---//g' | sed 's/---.*//g' | sed 's/-/,/g' | sort -k1,1 -k2n -t "," > T1R_genes.bed ; cd ../
cd Cunningtonia_longiventralis/ ; rm * ; grep ">Cunningtonia" /scicore/home/salzburg/polica0000/Phylogeny_analysis_Olfactory_receptor/T1R_final/Orthogroup_cluserting/ALL_T1R.fa | sed 's/>//g' | sed 's/Cunningtonia_longiventralis---//g' | sed 's/---.*//g' | sed 's/-/,/g' | sort -k1,1 -k2n -t "," > T1R_genes.bed ; cd ../
cd Cyphotilapia_frontosa/ ; rm * ; grep ">Cyphotilapia" /scicore/home/salzburg/polica0000/Phylogeny_analysis_Olfactory_receptor/T1R_final/Orthogroup_cluserting/ALL_T1R.fa | sed 's/>//g' | sed 's/Cyphotilapia_frontosa---//g' | sed 's/---.*//g' | sed 's/-/,/g' | sort -k1,1 -k2n -t "," > T1R_genes.bed ; cd ../
cd Cyprichromis_leptosoma/ ; rm * ; grep ">Cyprichromis" /scicore/home/salzburg/polica0000/Phylogeny_analysis_Olfactory_receptor/T1R_final/Orthogroup_cluserting/ALL_T1R.fa | sed 's/>//g' | sed 's/Cyprichromis_leptosoma---//g' | sed 's/---.*//g' | sed 's/-/,/g' | sort -k1,1 -k2n -t "," > T1R_genes.bed ; cd ../
cd Neolamprologus_multifasciatus/ ; rm * ; grep ">Neolamprologus" /scicore/home/salzburg/polica0000/Phylogeny_analysis_Olfactory_receptor/T1R_final/Orthogroup_cluserting/ALL_T1R.fa | sed 's/>//g' | sed 's/Neolamprologus_multifasciatus---//g' | sed 's/---.*//g' | sed 's/-/,/g' | sort -k1,1 -k2n -t "," > T1R_genes.bed ; cd ../
cd Simochromis_diagramma/ ; rm * ; grep ">Simochromis" /scicore/home/salzburg/polica0000/Phylogeny_analysis_Olfactory_receptor/T1R_final/Orthogroup_cluserting/ALL_T1R.fa | sed 's/>//g' | sed 's/Simochromis_diagramma---//g' | sed 's/---.*//g' | sed 's/-/,/g' | sort -k1,1 -k2n -t "," > T1R_genes.bed ; cd ../
cd Oreochromis_niloticus/ ; rm * ; grep ">Oreochromis_niloticus" /scicore/home/salzburg/polica0000/Phylogeny_analysis_Olfactory_receptor/T1R_final/Orthogroup_cluserting/ALL_T1R.fa | sed 's/>//g' | sed 's/Oreochromis_niloticus---//g' | sed 's/---.*//g' | sed 's/-/,/g' | sort -k1,1 -k2n -t "," > T1R_genes.bed ; cd ../


cd Bathybates_minor/ ; sbatch --job-name=T1R_Cluster --qos=1day -c 8 --mem-per-cpu=15g ClustEnrichTE.sh fBatMin1.PB.asm1.purge1.scaff1.polish2.decontaminated.fa.out /scicore/home/salzburg/polica0000/Cichlids_Genomes/PacBio_Genomes_Cichlids/Bathybates_minor/fBatMin1.PB.asm1.purge1.scaff1.polish2.decontaminated.fa.fai T1R_genes.bed fBatMin1.PB.asm1.purge1.scaff1.polish2.decontaminated.fa.align cds_coordinates_T1R_genes.Bathybates_minor.c123.csv ; cd ../
cd Cunningtonia_longiventralis/ ; sbatch --job-name=T1R_Cluster --qos=1day -c 8 --mem-per-cpu=15g ClustEnrichTE.sh fCunLon1.PB.asm1.purge1.fa.out /scicore/home/salzburg/polica0000/Cichlids_Genomes/PacBio_Genomes_Cichlids/Cunningtonia_longiventralis/fCunLon1.PB.asm1.purge1.fa.fai T1R_genes.bed fCunLon1.PB.asm1.purge1.fa.align cds_coordinates_T1R_genes.Cunningtonia_longiventralis.c123.csv ; cd ../
cd Cyphotilapia_frontosa/ ; sbatch --job-name=T1R_Cluster --qos=1day -c 8 --mem-per-cpu=15g ClustEnrichTE.sh fCypFro1.PB.asm1.purge1.scaff1.polish2.decontaminated.fa.out /scicore/home/salzburg/polica0000/Cichlids_Genomes/PacBio_Genomes_Cichlids/Cyphotilapia_frontosa/fCypFro1.PB.asm1.purge1.scaff1.polish2.decontaminated.fa.fai T1R_genes.bed fCypFro1.PB.asm1.purge1.scaff1.polish2.decontaminated.fa.align cds_coordinates_T1R_genes.Cyphotilapia_frontosa.c123.csv ; cd ../
cd Cyprichromis_leptosoma/ ; sbatch --job-name=T1R_Cluster --qos=1day -c 8 --mem-per-cpu=15g ClustEnrichTE.sh fCypLep1.PB.asm1.purge1.scaff1.polish2.decontaminated.fa.out /scicore/home/salzburg/polica0000/Cichlids_Genomes/PacBio_Genomes_Cichlids/Cyprichromis_leptosoma/fCypLep1.PB.asm1.purge1.scaff1.polish2.decontaminated.fa.fai T1R_genes.bed fCypLep1.PB.asm1.purge1.scaff1.polish2.decontaminated.fa.align cds_coordinates_T1R_genes.Cyprichromis_leptosoma.c123.csv ; cd ../
cd Neolamprologus_multifasciatus/ ; sbatch --job-name=T1R_Cluster --qos=1day -c 8 --mem-per-cpu=15g ClustEnrichTE.sh fNeoMul12_final.fa.out /scicore/home/salzburg/polica0000/Cichlids_Genomes/PacBio_Genomes_Cichlids/Neolamprologus_multifasciatus/fNeoMul12_final.fa.fai T1R_genes.bed fNeoMul12_final.fa.align cds_coordinates_T1R_genes.Neolamprologus_multifasciatus.c123.csv ; cd ../
cd Simochromis_diagramma/ ; sbatch --job-name=T1R_Cluster --qos=1day -c 8 --mem-per-cpu=15g ClustEnrichTE.sh GCF_900408965.1_fSimDia1.1_genomic.fna.out /scicore/home/salzburg/polica0000/Cichlids_Genomes/PacBio_Genomes_Cichlids/Simochromis_diagramma/GCF_900408965.1_fSimDia1.1_genomic.fna.fai T1R_genes.bed GCF_900408965.1_fSimDia1.1_genomic.fna.align cds_coordinates_T1R_genes.Simochromis_diagramma.c123.csv ; cd ../																																														 
cd Oreochromis_niloticus/ ; sbatch --job-name=T1R_Cluster --qos=1day -c 8 --mem-per-cpu=25g ClustEnrichTE.sh GCF_001858045.2_O_niloticus_UMD_NMBU_genomic.fna.out /scicore/home/salzburg/polica0000/Cichlids_Genomes/Oreochromis_niloticus/GCF_001858045.2_O_niloticus_UMD_NMBU_genomic.fna.fai T1R_genes.bed GCF_001858045.2_O_niloticus_UMD_NMBU_genomic.fna.align cds_coordinates_T1R_genes.Oreochromis_niloticus.c123.csv ; cd ../




