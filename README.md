# Lake_Tanganyika_cichlids_chemoreception

## I - Chemoreceptor number estimation with the depth-of-coverage method 

The depth-of-covergae method to estimate the number of chemoreceptors per species can be reproduced using the following scripts:

- [DoC_chemoreceptors.sh](DoC_chemoreceptors/DoC_chemoreceptors.sh) : Reads alignment and cleaning, coverage computation on chemoreceptor genes exons as well as on BUSCO genes exons. Reads files and the Oreochromis niloticus reference genome (GCF_001858045.2) can be retrieved on GenBank. 
- [DoC_chemoreceptors.R](DoC_chemoreceptors/DoC_chemoreceptors.R) : Process coverage values to estimate the number of genes

All files necessary to run these scripts are available on Zenodo (), in the directory "DoC_chemoreceptors"

## II -Transposable elements

- [TE_enrichment.sh](Transposable_elements/TE_enrichment.sh) : Generate bed files with chemoreceptors coordinates and launch the script ClustEnrichTE.sh
- [ClustEnrichTE.sh](Transposable_elements/ClustEnrichTE.sh) : Generate bed files of transposable elements coordinates, generate a kimura distance table and launch ClustEnrichTE.R
- [ClustEnrichTE.R](Transposable_elements/ClustEnrichTE.R) : Analysis of transposable elements enrichment in chemoreceptor genes clusters, and analysis of their age

All files necessary to run these scripts, including RepeatMasker outputs, are available on Zenodo (), in the directory "Transposable_elements"


## III - Lamellae number and brain volume analysis

- [TCichlids_lamellae.R](Cichlids_lamellae.R) : Generate bed files with chemoreceptors coordinates and launch the script ClustEnrichTE.sh

Pictures of olfactory epithelium, and all files necessary to run this script are available on Zenodo (), in the directory "Olfactory_epithelium" and in the directory "Dataset_and_various_data_Cichlids"


## IV - Chemoreception analysis

- [CR_number_ALL_cichlids.R](CR_number_ALL_cichlids.R) : Various plots and analysis of the number of chemoreceptors in cichlids. 
- [Chemoreception_vs_ecology.R](Chemoreception_vs_ecology.R) : Correlations between chemoreception (chemoreceptor numbers and morphology measures) and ecological factors

All files necessary to run these scripts are available on Zenodo (), in the directory "Dataset_and_various_data_Cichlids"

## V - RERconverge analysis

- [Launch_RER_cichlids.sh](RERconverge_analysis/Launch_RER_cichlids.sh) : Prepare files to extract gene sequences from every cichlid species and run RERconverge
- [Align_and_call_variants.sh](RERconverge_analysis/Align_and_call_variants.sh) : Align cichlid reads to the O. niloticus reference genome 
- [Generate_consensus_exons.sh](RERconverge_analysis/Generate_consensus_exons.sh) : Extract exon sequences for every gene and every cichlid species, based on reads mapping
- [Align_trim_phangorn.sh](RERconverge_analysis/Align_trim_phangorn.sh) : Align genes and generate ML trees with phangorn ( with the accessory scripts [Estimate_ML_tree_Phangorn.cichlids.sh] (RERconverge_analysis/Estimate_ML_tree_Phangorn.cichlids.sh) and [Estimate_ML_tree_Phangorn.cichlids.R] (RERconverge_analysis/Estimate_ML_tree_Phangorn.cichlids.R) )

- [RER_cichlids.R](RERconverge_analysis/RER_cichlids.R) : Launch RERconverge
- [Analyse_RERconverge_results.R](RERconverge_analysis/Analyse_RERconverge_results.R) : Analyse the results of RERconverge


All files necessary to run these scripts are available on Zenodo (), in the directory "RERconverge" and in the directory "Dataset_and_various_data_Cichlids"


## VI - Olfactory epithelium RNA-seq analysis

- [RNAseq_data_processing.sh](RNAseq_data_processing.sh) : Reads cleaning, alignment and counts. 
- [Cichlids_RNAseq.R](Cichlids_RNAseq.R) : Analysis of olfactory epithelium RNA-seq data. 

All files necessary to run these scripts are available on Zenodo (), in the directory "Dataset_and_various_data_Cichlids" and in the directory "OE_RNAseq". Reads files are available on GenBank. 