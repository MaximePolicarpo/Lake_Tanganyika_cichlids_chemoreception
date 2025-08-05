#### Set environement and load packages ---------------------------------

set.seed(2712)
library("dplyr")
library("ggpubr")
library("stringr")

args = commandArgs(trailingOnly=TRUE)

#### First define clusters based on the gene bed file ---------------------------------

library("plyranges")
library("GenomicRanges")

gene_coords <- read.table(args[1], sep=",", header=FALSE)
colnames(gene_coords) <- c("seqnames", "begin", "stop")


#extend each gene 30kb upstream and downstream
gene_coords <- gene_coords %>% mutate(start = begin - 30000) %>% mutate(end = stop + 30000)
#Put the table as a grange object
genes_irange <- gene_coords %>% as_granges()
#Reduce the table to merge overlapping results
genes_disjoin <- reduce(genes_irange,with.revmap=TRUE)
#transform the table in a data.frame
clusters_regions <- as.data.frame(genes_disjoin)

#Define the clusters and the genes in each cluster
genes_clustered_df <- data.frame(NULL)
cluster_nb=0
singleton_nb=0
for(i in 1:nrow(clusters_regions)){
  
  curr_line=clusters_regions[i,]
  curr_scaff=curr_line$seqnames
  curr_start=curr_line$start
  curr_stop=curr_line$end
  
  curr_genes <- 
    gene_coords %>% 
    filter(seqnames == curr_scaff) %>%
    filter(begin >= curr_start) %>%
    filter(stop <= curr_stop) %>%
    dplyr::select(seqnames, begin, stop) 
  
  if( nrow(curr_genes) >= 2) {
    cluster_nb = cluster_nb+1
    curr_genes <- curr_genes %>% mutate(cluster_name = paste("cluster", cluster_nb, sep=""))
  }
  
  if( nrow(curr_genes) < 2) {
      singleton_nb = singleton_nb+1
      curr_genes <- curr_genes %>% mutate(cluster_name = paste("singleton", singleton_nb, sep=""))
  }
    
  
  genes_clustered_df <- rbind(genes_clustered_df, curr_genes)
  
}


gene_coords <- genes_clustered_df
colnames(gene_coords) <- c("scaffold", "start", "stop", "cluster_name")


#detach the plyrange package that creates conflicts ..
detach("package:plyranges", unload=TRUE)



write.table(gene_coords, 
            file="gene_coords_clustered.tsv", 
            quote=FALSE, sep=",", 
            row.names = FALSE, 
            col.names = FALSE)


#### Graphics - Density of TE and Olfactory genes on scaffolds  ---------------------------------


#Choose a window and step size for the graph
windows_size <- 200000
step_size <- 10000
options(scipen = 999)

#Use bedops to generate the sliding windows 
system(paste("module purge ; module load BEDOPS ; bedops --chop", windows_size, "--stagger", step_size, "Genome.bed > sliding_windows.bed"))

#Import a summary table of the genome (each line : scaffold , 0, scaffold length)
Genome_coord <- read.table("Genome.bed", sep="\t", header=FALSE)
colnames(Genome_coord) <- c("scaffold", "start", "end")

#Import a table with the location of TE
transposons_coords <- read.table("transposons.bed", sep="\t", header=FALSE)
colnames(transposons_coords) <- c("scaffold", "start", "stop", "transposon_name")

#Remove TE that overlap with our genes
#args[2] <- "/scicore/home/salzburg/polica0000/Phylogeny_analysis_Olfactory_receptor/Chemoreceptors_exons_positions/V2R/cds_coordinates_V2R_genes.Oreochromis_niloticus.c123.csv"
exons_coord <- read.table(args[2], header=FALSE, sep="\t")

#exons_coord <- read.table("/scicore/home/salzburg/polica0000/Phylogeny_analysis_Olfactory_receptor/Chemoreceptors_exons_positions/V2R/cds_coordinates_V2R_genes.Bathybates_minor.c123.csv",
#                          header=FALSE,
#                          sep="\t")

colnames(exons_coord) <- c("scaffold", "start", "end")

overlapped_transposons <- as.data.frame(NULL)
for (row in 1:nrow(exons_coord)) {
  curr_scaffold <- exons_coord[row, "scaffold"]
  curr_start_exon <- exons_coord[row, "start"]
  curr_end_exon  <- exons_coord[row, "end"]
  
  curr_overlapped_transposons <- 
    transposons_coords %>%
    filter(scaffold == curr_scaffold & start >= curr_start_exon & stop <= curr_end_exon |
           scaffold == curr_scaffold & start <= curr_start_exon & stop >= curr_start_exon | 
           scaffold == curr_scaffold & start <= curr_end_exon & stop >= curr_end_exon | 
           scaffold == curr_scaffold & start <= curr_start_exon & stop >= curr_end_exon)
  
  overlapped_transposons <- 
    rbind(overlapped_transposons, curr_overlapped_transposons)
  
  
}

transposons_coords <- 
  anti_join(transposons_coords, overlapped_transposons,
            by = c("scaffold", "start", "stop"))



#Import sliding windows table generated with bedops
sliding_windows <- read.table("sliding_windows.bed", sep="\t", header=FALSE)
colnames(sliding_windows) <- c("scaffold", "start", "stop")

#Import the kimura distance values of each TE element
kimura_dists <- read.table("transposons_kimura_dist.csv", sep=",", header=FALSE)
colnames(kimura_dists) <- c("scaffold", "start", "stop", "kimura_dist")

#Merge the TE kimura distance and TE locations
kimura_dists <- left_join(transposons_coords, kimura_dists, by=c("scaffold", "start", "stop"))

#Remove TE for which RepeatMasker does not give the kimura distance
kimura_dists <- kimura_dists %>% filter(! is.na(kimura_dist))

#Compute TE density over each cluster
cluster_list <- gene_coords %>% pull(cluster_name)
cluster_list <- table(cluster_list)
cluster_list <- names(cluster_list)[cluster_list >= 2]

#Define the cluster coordinate and size and add 20,000 bp to its size
cluster_coords <- 
  as.data.frame(gene_coords %>% 
  filter(cluster_name %in% cluster_list) %>%
  group_by(scaffold, cluster_name) %>%
  summarise(start_coord = min(start), stop_coord=max(stop))) %>%
  mutate(cluster_size = stop_coord - start_coord + 20000)
  
#Lets compute the TE density (in the cluster and 10,000 bp around it)

colnames(cluster_coords) <- c("scaffold", "cluster_name", "start_cluster", "stop_cluster", "cluster_size")
colnames(transposons_coords) <- c("scaffold", "start_TE", "stop_TE", "transposon_name")

Observed_te_densities_df <- 
  left_join(cluster_coords, transposons_coords, by="scaffold") %>%
  filter(as.numeric(start_TE) >= (as.numeric(start_cluster) - 10000)) %>%
  filter(as.numeric(stop_TE) <= (as.numeric(stop_cluster) + 10000)) %>%
  mutate(TE_length = as.numeric(stop_TE) - as.numeric(start_TE))  %>%
  group_by(scaffold, cluster_name, start_cluster, stop_cluster, cluster_size) %>%
  summarise(TE_length = sum(TE_length)) %>%
  mutate(TE_density = TE_length/cluster_size)
  

#Draw the density graph per clusters 

colnames(gene_coords) <- c("scaffold", "start_gene", "stop_gene", "cluster_name")
colnames(sliding_windows) <- c("scaffold", "start_window", "stop_window")


scaffold_list <- Observed_te_densities_df %>% pull(scaffold) %>% unique()
transposon_nb_per_wd_all <- as.data.frame(NULL)
gene_nb_per_wd_all <- as.data.frame(NULL)
for(curr_scaffold in scaffold_list){
  sliding_windows_curr <- sliding_windows %>% filter(scaffold==curr_scaffold) #import sliding windows on the specified scaffold
  transposons_coords_curr <- transposons_coords %>% filter(scaffold==curr_scaffold) #extract TE present on the scaff
  sliding_windows_transposon <- left_join(sliding_windows_curr, transposons_coords_curr, by = "scaffold")
  sliding_windows_genes <- left_join(sliding_windows_curr, gene_coords, by = "scaffold")
  
  transposon_density_per_wd <-  #extract the number of TE and their length over windows
    as.data.frame(sliding_windows_transposon %>%
                    filter(as.numeric(start_TE) >= as.numeric(start_window)) %>%
                    filter(as.numeric(stop_TE) <= as.numeric(stop_window)) %>%
                    mutate(TE_length = stop_TE - start_TE) %>% 
                    group_by(scaffold, start_window, stop_window) %>%
                    dplyr::summarise(transposon_nb = n(), 
                                     transposon_length = sum(TE_length)
                                     ) %>%
                    ungroup())
  
  gene_density_per_wd <- #extract the number of gene and their length over windows
    as.data.frame(sliding_windows_genes %>%
                    filter(as.numeric(start_gene) >= as.numeric(start_window)) %>%
                    filter(as.numeric(stop_gene) <= as.numeric(stop_window)) %>%
                    mutate(gene_length = stop_gene - start_gene) %>% 
                    group_by(scaffold, start_window, stop_window) %>%
                    dplyr::summarise(gene_nb = n(), 
                                     gene_length = sum(gene_length)
                                     ) %>%
                    ungroup())
  
  
  transposon_nb_per_wd_all <- rbind(transposon_nb_per_wd_all, transposon_density_per_wd)
  gene_nb_per_wd_all <- rbind(gene_nb_per_wd_all, gene_density_per_wd)
  
  
}

#extract windows of interest on the results
sliding_windows_subset <- 
  sliding_windows %>%
  filter(scaffold %in% scaffold_list)
colnames(sliding_windows_subset) <- c("scaffold", "start_window", "stop_window")

transposon_nb_per_wd_all <- 
  full_join(transposon_nb_per_wd_all, 
            sliding_windows_subset, by=c("scaffold", "start_window", "stop_window"))
transposon_nb_per_wd_all[is.na(transposon_nb_per_wd_all)] <- 0

gene_nb_per_wd_all <- 
  full_join(gene_nb_per_wd_all, 
            sliding_windows_subset, by=c("scaffold", "start_window", "stop_window"))
gene_nb_per_wd_all[is.na(gene_nb_per_wd_all)] <- 0


#Lets make summary tables
Summary_table_gene <- gene_nb_per_wd_all %>% mutate(type = "gene")
colnames(Summary_table_gene) <- c("scaffold", "start_window", "stop_window", "number", "length","type")
Summary_table_TE <- transposon_nb_per_wd_all %>% mutate(type = "TE")
colnames(Summary_table_TE) <- c("scaffold", "start_window", "stop_window", "number", "length","type")


#combine TE and gene tables
Summary_TE_gene <- rbind(Summary_table_TE, Summary_table_gene)

#compute the density (= length of gene or TE divided by window size in bp)
Summary_TE_gene <- 
  Summary_TE_gene %>% 
  mutate(bp_density = length/windows_size) 


#Graphic of density in bp length

plot_list = list()
for (curr_scaff in scaffold_list) {
  p = ggplot((Summary_TE_gene %>% 
                filter(scaffold == curr_scaff)), aes(x=start_window, y=bp_density, color=type))+
    geom_line() +
    theme_bw() +
    ylab("Density (%bp)") +
    xlab(curr_scaff)
  plot_list[[curr_scaff]] = p
}

pdf("./Density_TE_genes_along_scaffolds.pdf",width = 6.34,  height = 4.61)
for (curr_scaff in scaffold_list) {
  print(plot_list[[curr_scaff]])
}
dev.off()

#Classify TEs based on their name

transposon_table_class <- 
  transposons_coords %>%
  mutate(TE_class = case_when(
    str_detect(transposon_name, "DNA") ~ "DNA",
    str_detect(transposon_name, "LINE") ~ "LINE",
    str_detect(transposon_name, "SINE") ~ "SINE",
    str_detect(transposon_name, "LTR") ~ "LTR",
    str_detect(transposon_name, "Helitron") ~ "DNA",
    str_detect(transposon_name, "Retroposon") ~ "Retroposon",
    str_detect(transposon_name, "Unspecified") ~ "Unspecified"
  )) 

TE_class_df <- 
  left_join(cluster_coords, transposon_table_class, by="scaffold") %>%
  filter(as.numeric(start_TE) >= (as.numeric(start_cluster) -10000)) %>%
  filter(as.numeric(stop_TE) <= (as.numeric(stop_cluster) + 10000)) %>%
  mutate(TE_length = as.numeric(stop_TE) - as.numeric(start_TE))  %>%
  group_by(cluster_name,TE_class) %>%
  dplyr::summarise(count = n(), 
                   sum_length = sum(TE_length)
                   ) %>%
  ungroup()



## number of element per class in the cluster(s)

pdf("./TE_class_number_per_cluster.pdf",width = 6.34,  height = 4.61)

TE_class_df %>%
  ggplot(., aes(y=count, x=cluster_name, fill=TE_class)) +
  geom_col(position="dodge") +
  theme_bw() +
  xlab("TE copy number")

dev.off()




#### Are TE "young" ?  ---------------------------------

#Classify TE in the kimura distance table
kimura_dists <- 
  kimura_dists %>%
  mutate(TE_class = case_when(
    str_detect(transposon_name, "DNA") ~ "DNA",
    str_detect(transposon_name, "LINE") ~ "LINE",
    str_detect(transposon_name, "SINE") ~ "SINE",
    str_detect(transposon_name, "LTR") ~ "LTR",
    str_detect(transposon_name, "Helitron") ~ "DNA",
    str_detect(transposon_name, "Retroposon") ~ "Retroposon",
    str_detect(transposon_name, "Unspecified") ~ "Unspecified"
  ))

colnames(kimura_dists) <- c("scaffold", "start_TE", "stop_TE", "transposon_name", "kimura_dist", "TE_class")

#Extract TE in cluster (-/+ 10,000 bp) vs TE outside cluster

kimura_dists_in_cluster <- 
  left_join(cluster_coords, kimura_dists, by="scaffold") %>%
  filter(as.numeric(start_TE) >= (as.numeric(start_cluster) - 10000)) %>%
  filter(as.numeric(stop_TE) <= (as.numeric(stop_cluster) + 10000))  %>%
  mutate(Group = "Cluster") 

#Extract TE genome wide and remove TE in our cluster(s) .. 

kimura_dists_gw <- 
  kimura_dists %>% 
  mutate(Group = "Genome-wide")
kimura_dists_gw <- 
  anti_join(kimura_dists_gw, kimura_dists_in_cluster, by = c("scaffold", "start_TE", "stop_TE"))

kimura_dists_in_cluster <- 
  kimura_dists_in_cluster %>% dplyr::select(colnames(kimura_dists_gw))


#Only retain TE class for which there is at-least 3 element in clusters

TE_class_both <- 
  kimura_dists_in_cluster %>% 
  group_by(TE_class) %>% 
  summarise(count = n()) %>% 
  filter(count > 2) %>% 
  pull(TE_class) %>% unique()



#Represent the kimura distance distrib on clusters vs outside clusters and
#perform wilcox.test


options(scipen=0)

pdf("./Comparison_TE_ages.pdf",width = 15.34,  height = 4.61)

rbind(kimura_dists_in_cluster,kimura_dists_gw) %>%
  filter(TE_class %in% TE_class_both) %>%
  ggplot(., aes(x=TE_class, y=kimura_dist , fill=Group)) +
  geom_boxplot() +
  theme_classic() +
  xlab("TE class") +
  ylab("Kimura distance") +
  stat_compare_means(method = "wilcox.test")

dev.off()


options(scipen = 999)



#### Are the clusters enriched in TEs ? ---------------------------------


### 1 - Compute the TE density over each cluster

#First extract cluster names
cluster_list <- gene_coords %>% pull(cluster_name) 
cluster_list <- table(cluster_list)
cluster_list <- names(cluster_list)[cluster_list >= 2]

#Now extract cluster and add 20,000 bp to its size
cluster_coords <- 
  as.data.frame(gene_coords %>% 
                  filter(cluster_name %in% cluster_list) %>%
                  group_by(scaffold, cluster_name) %>%
                  summarise(start_cluster = min(start_gene), stop_cluster=max(stop_gene))) %>%
  mutate(cluster_size = stop_cluster - start_cluster + 20000)

#Extract the TE density on the region
Observed_te_densities_df <- 
  left_join(cluster_coords, transposons_coords, by="scaffold") %>%
  filter(as.numeric(start_TE) >= (as.numeric(start_cluster) - 10000)) %>%
  filter(as.numeric(stop_TE) <= (as.numeric(stop_cluster) + 10000)) %>%
  mutate(TE_length = as.numeric(stop_TE) - as.numeric(start_TE))  %>%
  group_by(scaffold, cluster_name, start_cluster, stop_cluster, cluster_size) %>%
  summarise(TE_length = sum(TE_length)) %>%
  mutate(TE_density = TE_length/cluster_size)

#Compute the TE density over the whole genome

genome_wide_bp <- Genome_coord %>% pull(end) %>% sum()

transposon_length_gw <- 
  transposons_coords %>% 
  mutate(length = as.numeric(stop_TE) - as.numeric(start_TE)) %>% 
  pull(length) %>% sum()


density_GW <- transposon_length_gw/genome_wide_bp

#Now compute the density over every possible windows in the genome

options(scipen = 999)
cluster_list <- gene_coords %>% pull(cluster_name)
cluster_list <- table(cluster_list)
cluster_list <- names(cluster_list)[cluster_list >= 2]


probabilities_values <- c()
Summary_results_wd <- data.frame(NULL)
All_results_wd <- data.frame(NULL)
for(cluster in cluster_list){
  
  #Create sliding windows of the size of the cluster and with a step size of 10,000 bp
  
  step_size <- 10000
  cluster_size <- cluster_coords %>% filter(cluster_name == cluster) %>% pull(cluster_size)
  system(paste("module purge ; module load BEDOPS ; bedops --chop", cluster_size, "--stagger", step_size, "Genome.bed > cluster_sliding_windows.bed"))
  sliding_windows <- read.table("cluster_sliding_windows.bed", sep="\t", header=FALSE)
  colnames(sliding_windows) <- c("scaffold", "start_window", "stop_window")
  
  #remove windows that are below the cluster size
  sliding_windows <- 
    sliding_windows %>% 
    mutate(length_window = stop_window-start_window) %>% 
    filter(length_window >= cluster_size)

  #remove every windows overlapping with the V2R region
  scaffold_cluster <- 
    cluster_coords %>% 
    filter(cluster_name == cluster) %>% 
    pull(scaffold)
  
  begin_cluster <- 
    cluster_coords %>% 
    filter(cluster_name == cluster) %>% 
    pull(start_cluster)
  end_cluster <- 
    cluster_coords %>% 
    filter(cluster_name == cluster) %>% 
    pull(stop_cluster)
  
  sliding_windows <- 
    sliding_windows %>% 
    filter(!(scaffold == scaffold_cluster & start_window >= (begin_cluster-10005) & start_window <= (end_cluster+10005)))
  sliding_windows <- 
    sliding_windows %>% 
    filter(!(scaffold == scaffold_cluster & stop_window >= (begin_cluster-10005) & stop_window <= (end_cluster+10005)))
  sliding_windows <- 
    sliding_windows %>% filter(!(scaffold == scaffold_cluster & start_window <= (end_cluster+10005) & stop_window >= (end_cluster+10005)))

  #Compute the TE density in every windows (do that scaff by scaff to reduce memory consumption)

  scaffold_list <- 
    sliding_windows %>%
    pull(scaffold) %>% 
    unique()

  TE_density_GW <- data.frame(NULL)
  for(curr_scaffold in scaffold_list){
    sliding_windows_curr <- sliding_windows %>% filter(scaffold==curr_scaffold)
    transposons_coords_curr <- transposons_coords %>% filter(scaffold==curr_scaffold)
    sliding_windows_transposon <- left_join(sliding_windows_curr, transposons_coords_curr, by = "scaffold")
  
    #compute the number of TE + the TE total length over each window
    transposon_density_per_wd <- 
      as.data.frame(sliding_windows_transposon %>%
                    filter(as.numeric(start_TE) >= as.numeric(start_window)) %>%
                    filter(as.numeric(stop_TE) <= as.numeric(stop_window)) %>%
                    mutate(transposon_length = stop_TE - start_TE) %>% 
                    group_by(scaffold, start_window, stop_window) %>%
                    dplyr::summarise(transposon_nb = n(), 
                                     TE_length = sum(transposon_length)
                                     ) %>%
                    ungroup())
  
    TE_density_GW <- rbind(TE_density_GW, transposon_density_per_wd)
  
  
  }

  #supress useless variables that are heavy in memory
  rm(transposon_density_per_wd)
  rm(sliding_windows_transposon)
  rm(transposons_coords_curr)
  rm(sliding_windows_curr)

  TE_density_GW <- full_join(TE_density_GW, sliding_windows, by=c("scaffold", "start_window", "stop_window"))
  #If one window is not present in the TE_density_GW table, it means it does not contain TE, so length and nb = 0
  TE_density_GW[is.na(TE_density_GW)] <- 0
  #Compute TE density (TE length / window size)
  TE_density_GW <- TE_density_GW %>% mutate(TE_density = TE_length/length_window)



  ##Compute the precentage of windows with a TE density  
  #above the observed one and make a graph

  current_observed_density <- 
    Observed_te_densities_df %>% 
    filter(cluster_name == cluster) %>% 
    pull(TE_density)
  
  TE_densities_values <- 
    TE_density_GW %>% 
    pull(TE_density)
  
  #the probability = number of window with an equal or superior TE density that the
  #observed one, divided by the total number of windows
  probability_density <- 
    length(TE_densities_values[TE_densities_values >= current_observed_density]) / length(TE_densities_values)
  probabilities_values <- c(probabilities_values,probability_density)


  ##Make a Summary of GW, genes and ALL windows densities

  TE_densities_values_df <- 
    as.data.frame(TE_densities_values) %>% 
    mutate(Region = "All genome windows") %>% 
    mutate(cluster_name = cluster)
  colnames(TE_densities_values_df) <- c("TE_density", "Region", "Cluster")

  density_GW_df <- 
    as.data.frame(density_GW) %>% 
    mutate(Region = "Genome") %>% 
    mutate(cluster_name = cluster)
  colnames(density_GW_df) <- c("TE_density", "Region", "Cluster")

  observed_density <- 
    as.data.frame(current_observed_density) %>%
    mutate(Region = "observed cluster") %>%
    mutate(cluster_name = cluster)
  colnames(observed_density) <- c("TE_density", "Region", "Cluster")


  all_computed_densities <- rbind(TE_densities_values_df, density_GW_df, observed_density)
  
  All_results_wd <- rbind(All_results_wd, all_computed_densities)
  
  all_computed_densities <- 
    all_computed_densities %>% 
    group_by(Region, Cluster) %>%
    summarise(mean_value = mean(TE_density), sd_value = sd(TE_density))
  
  Summary_results_wd <- rbind(Summary_results_wd, all_computed_densities)
  
} 
  
#Make a nice table for the probability of each cluster
probabilities_values_df <- 
  as.data.frame(probabilities_values) %>% 
  mutate(cluster_name = cluster_list)

#Write the table and write the results tables
write.table(probabilities_values_df, file="probabilities_clusters.tsv", quote=FALSE, sep="\t", row.names = FALSE, col.names = FALSE)
write.table(Summary_results_wd, file="Summary_results_wd.tsv", quote=FALSE, sep="\t", row.names = FALSE, col.names = FALSE)
write.table(All_results_wd, file="All_results_wd.tsv", quote=FALSE, sep="\t", row.names = FALSE, col.names = FALSE)



#Draw a boxplot with the TE density over windows + a point for the cluster TE density
pdf("./Density_comparison_all_genome_windows.pdf",width = 6.34,  height = 4.61)

All_results_wd %>%
  filter(Region == "All genome windows") %>%
  ggplot(., aes(y=TE_density, x=Cluster)) +
  geom_boxplot(color="gray") +
  theme_classic() +
  ylab("TE density (%bp)") +
  geom_point(data = (All_results_wd %>% filter(Region == "observed cluster")), mapping= aes(y=TE_density, x=Cluster), color="red", size=3)
  #geom_point(data = (All_results_wd %>% filter(Region == "Genome")), mapping= aes(y=TE_density, x=Cluster), color="blue", size=3) +
  #geom_hline(yintercept = density_GW, color="blue", linetype="dashed")
  
dev.off()



#For each cluster, draw the distribution of genome wide windows TE densities,
#and draw a line where the observed value of the gene cluster is

plot_list = list()
for (curr_cluster in cluster_list) {
  p = All_results_wd %>%
    filter(Region == "All genome windows") %>%
    filter(Cluster == curr_cluster) %>%
    ggplot(., aes(x=TE_density)) +
    geom_histogram() +
    theme_classic() +
    xlab("TE density (%bp)") +
    ylab("Cluster size region number") +
    geom_vline(xintercept= density_GW, color="black") +
    geom_vline(xintercept= All_results_wd %>% filter(Region == "observed cluster") %>% filter(Cluster == curr_cluster) %>% pull(TE_density), color="red") +
    geom_vline(xintercept=quantile((All_results_wd %>% filter(Region == "All genome windows") %>% filter(Cluster == curr_cluster) %>% pull(TE_density)),probs=c(0.025)), color="gray", linetype="dashed") + 
    geom_vline(xintercept=quantile((All_results_wd %>% filter(Region == "All genome windows") %>% filter(Cluster == curr_cluster) %>% pull(TE_density)),probs=c(0.975)), color="gray", linetype="dashed") +
    annotate("text", x = (All_results_wd %>% filter(Region == "observed cluster") %>% filter(Cluster == curr_cluster) %>% pull(TE_density) - 0.02), y = Inf, label = as.character(round(1-probabilities_values_df %>% filter(cluster_name == curr_cluster) %>% pull(probabilities_values), 3)), color="red", angle=90, vjust = "inward", hjust = "inward") +
    coord_cartesian(clip = "off")
  plot_list[[curr_cluster]] = p
}

pdf("./Histogram_TEdensity_all_regions.pdf",width = 6.34,  height = 4.61)
for (curr_cluster in cluster_list) {
  print(plot_list[[curr_cluster]])
}
dev.off()


### 2 - If there is more than 1 cluster, then we will also use a second method
### 2 - This method makes a probability at the level of the gene family for the TE density
### 2 - instead of by cluster independently.

if(length(cluster_list) >= 2){
  
  possible_scaffold_per_cluster <- as.data.frame(NULL)
  
  for(cluster in cluster_list){

    #Extract the TE density computed on the cluster
    Observed_mean_family_density <- 
     Observed_te_densities_df %>%
     pull(TE_density) %>%
     mean()

    #Extract the cluster size
    cluster_size <- 
      cluster_coords %>% 
      filter(cluster_name == cluster) %>% 
      pull(cluster_size)
    
    #for each cluster, define on which scaffold it could be placed 
    #Remove scaffold with a length < cluster size
    
    Possible_scaffolds_df <- 
      as.data.frame(
        Genome_coord %>%
          filter(end > (cluster_size+1)) %>%
          pull(scaffold)) %>%
      mutate(Cluster = cluster) %>%
      mutate(Cluster_size = cluster_size)
    colnames(Possible_scaffolds_df) <- c("scaffold", "Cluster", "Cluster_size") 
    
    possible_scaffold_per_cluster <- rbind(possible_scaffold_per_cluster, Possible_scaffolds_df)
    
  }
  
  #Initiate some variables
  all_simulated_densities <- c()
  all_simulations_densities <- data.frame(NULL)
  simulation_nb = 1
  
  #Lets perform 10,000 simulations per cluster
  while(nrow(all_simulations_densities) < (length(cluster_list) * 10000)){
    
    #First, extract one random scaffold for each cluster
    random_positions <- as.data.frame(NULL)
    for(cluster in cluster_list){
      
      random_picked_scaff <- 
        sample_n((possible_scaffold_per_cluster %>% 
                    filter(Cluster == cluster)), 1,
                 replace = TRUE)
      
      random_positions <- rbind(random_positions, random_picked_scaff)
    }
    
    
    random_picked_scaffolds <- 
      left_join(random_positions, 
                Genome_coord, 
                by="scaffold")
    
    #For each selected scaffold, define a random windows (size = cluster size)
    Simulated_regions <-
      random_picked_scaffolds %>%
      mutate(min_cluster_start = end - Cluster_size - 1) %>%
      rowwise() %>%
      mutate(random_start = runif(1, min = 1 , max = min_cluster_start)) %>%
      ungroup() %>%
      mutate(random_end = random_start+Cluster_size)
    
    #Verify that there is no overlap in the selected regions 
    
    library("plyranges")
    library("GenomicRanges")
    
    verification_non_overlap_regions <- 
      as.data.frame(Simulated_regions %>% 
                      dplyr::select(scaffold, random_start, random_end))
    colnames(verification_non_overlap_regions) <- c("seqnames", "start", "end")
    verification_irange <- verification_non_overlap_regions %>% as_granges()
    verification_disjoin <- reduce(verification_irange,with.revmap=TRUE)
    number_nonoverlap <- nrow(as.data.frame(verification_disjoin))
  
    detach("package:plyranges", unload=TRUE)
  
    #Keep on if there is no overlap
    if(number_nonoverlap == length(cluster_list)){
    
      simulation_nb = simulation_nb+1
      
      #Compute the TE density over each random regions
      simulated_densities_per_region <- 
        left_join(Simulated_regions, transposons_coords, by="scaffold") %>%
        filter(start_TE >= random_start) %>%
        filter(stop_TE <= random_end) %>%
        mutate(transposon_length = stop_TE - start_TE) %>% 
        group_by(scaffold, random_start, random_end, Cluster_size) %>%
        dplyr::summarise(transposon_nb = n(), TE_length = sum(transposon_length)) %>%
        ungroup() %>%
        mutate(TE_density = TE_length / Cluster_size) 
    
      simulated_densities_per_region <- 
        full_join(Simulated_regions, 
                  simulated_densities_per_region, 
                  by=c("scaffold", "random_start", "random_end"))
      
      #If a region is not present, it means there is no TE, so TEdensity = 0 and TEnb=0
      simulated_densities_per_region[is.na(simulated_densities_per_region)] <- 0
    
      #Extract the TE densities values and put them in a dataframe with the simulation nb
      simulated_densities <- simulated_densities_per_region %>% pull(TE_density)
      all_simulated_densities <- c(all_simulated_densities, simulated_densities)
      all_simulated_densities <- 
        as.data.frame(all_simulated_densities) %>%
        mutate(Simulation_nb = simulation_nb)
    
      all_simulations_densities <- 
        rbind(all_simulations_densities, all_simulated_densities)
      
      #re-initialize the density vector
      all_simulated_densities <- c()
    }
  }
  
  #Make the mean of each simulation for all clusters
  mean_densities_df <- 
    all_simulations_densities %>%
    group_by(Simulation_nb) %>%
    summarise(mean_TE_density = mean(all_simulated_densities))
  
  #Extract the mean density values
  mean_densities <- 
    mean_densities_df %>% 
    pull(mean_TE_density)
  
  #Lets compute a probability = number of simulation where the mean density is 
  #higher or equal to the observed one / total number of simulation
  
  probability_density_sim <- 
    length(mean_densities[mean_densities > Observed_mean_family_density]) / length(mean_densities)
  
  probability_family <- 
    as.data.frame(probability_density_sim) %>% mutate(cluster_name = "Family")
  
  write.table(probability_family, file="probability_family.tsv", quote=FALSE, sep="\t", row.names = FALSE, col.names = FALSE)
  
  #Now lets make a graph combined with single cluster graph
    
  mean_densities_df <- 
    mean_densities_df %>%
    dplyr::select(mean_TE_density) %>%
    mutate(Region = "All genome windows") %>%
    mutate(Cluster = "Family")
  colnames(mean_densities_df) <- colnames(All_results_wd)
  
  All_results_wd_2 <- rbind(All_results_wd, mean_densities_df)
    
    
  observed_density <- as.data.frame(Observed_mean_family_density) %>% 
    mutate(Region = "All clusters") %>% 
    mutate(cluster_name = "Family")
  
  colnames(observed_density) <- c("TE_density", "Region", "Cluster")
  
  All_results_wd_2 <- rbind(All_results_wd_2, observed_density)


  pdf("./Density_comparison_Family_summary.pdf",width = 6.34,  height = 4.61)
  
  
  All_results_wd_2 %>%
    filter(Region == "All genome windows") %>%
    ggplot(., aes(y=TE_density, x=Cluster)) +
    geom_boxplot(color="gray") +
    theme_classic() +
    ylab("TE density (%bp)") +
    geom_point(data = (All_results_wd_2 %>% filter(Region == "observed cluster")), mapping= aes(y=TE_density, x=Cluster), color="red", size=3)+
    #geom_point(data = (All_results_wd_2 %>% filter(Region == "Genome")), mapping= aes(y=TE_density, x=Cluster), color="blue", size=3) +
    geom_point(data = (All_results_wd_2 %>% filter(Region == "All clusters")), mapping= aes(y=TE_density, x=Cluster), color="red", size=3) #+
    #geom_hline(yintercept = density_GW, color="blue", linetype="dashed")
  
    
  dev.off()
  
  write.table(All_results_wd_2, file="All_results_wd_2.tsv", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

  #Now lets draw the distribution graph .. 
  
  density_mean_sim <- mean(mean_densities)
  

  pdf("./Histogram_simulations_results_family.pdf",width = 6.34,  height = 4.61)
  
  All_results_wd_2 %>%
    filter(Region == "All genome windows") %>%
    filter(Cluster == "Family") %>%
    ggplot(., aes(x=TE_density)) +
    geom_histogram() +
    theme_classic() +
    xlab("TE density (%bp)") +
    ylab("Number of simulations") +
    geom_vline(xintercept= Observed_mean_family_density, color="red") +
    geom_vline(xintercept=quantile(All_results_wd_2 %>% filter(Region == "All genome windows") %>% filter(Cluster == "Family") %>% pull(TE_density),probs=c(0.025)), color="gray") +
    geom_vline(xintercept=quantile(All_results_wd_2 %>% filter(Region == "All genome windows") %>% filter(Cluster == "Family") %>% pull(TE_density),probs=c(0.975)), color="gray") +
    annotate("text", x = (Observed_mean_family_density - 0.02), y = Inf, label = as.character(round(1-probability_density_sim, 3)), color="red", angle=90, vjust = "inward", hjust = "inward") +
    geom_vline(xintercept= density_mean_sim, color="black") +
    coord_cartesian(clip = "off")
  
  
  dev.off()
  
  
  
}


pdf("./Density_comparison_Family_summary.pdf",width = 6.34,  height = 4.61)


All_results_wd_2 %>%
  filter(Region == "All genome windows") %>%
  ggplot(., aes(y=TE_density, x=Cluster)) +
  geom_boxplot(color="gray") +
  theme_classic() +
  ylab("TE density (%bp)") +
  geom_point(data = (All_results_wd_2 %>% filter(Region == "observed cluster")), mapping= aes(y=TE_density, x=Cluster), color="red", size=3)+
  #geom_point(data = (All_results_wd_2 %>% filter(Region == "Genome")), mapping= aes(y=TE_density, x=Cluster), color="blue", size=3) +
  geom_point(data = (All_results_wd_2 %>% filter(Region == "All clusters")), mapping= aes(y=TE_density, x=Cluster), color="red", size=3) #+
#geom_hline(yintercept = density_GW, color="blue", linetype="dashed")


dev.off()


pdf("./Histogram_simulations_results_family.pdf",width = 6.34,  height = 4.61)

All_results_wd_2 %>%
  filter(Region == "All genome windows") %>%
  filter(Cluster == "Family") %>%
  ggplot(., aes(x=TE_density)) +
  geom_histogram() +
  theme_classic() +
  xlab("TE density (%bp)") +
  ylab("Number of simulations") +
  geom_vline(xintercept= Observed_mean_family_density, color="red") +
  geom_vline(xintercept=quantile(All_results_wd_2 %>% filter(Region == "All genome windows") %>% filter(Cluster == "Family") %>% pull(TE_density),probs=c(0.025)), color="gray") +
  geom_vline(xintercept=quantile(All_results_wd_2 %>% filter(Region == "All genome windows") %>% filter(Cluster == "Family") %>% pull(TE_density),probs=c(0.975)), color="gray") +
  annotate("text", x = (Observed_mean_family_density - 0.02), y = Inf, label = as.character(round(1-probability_density_sim, 3)), color="red", angle=90, vjust = "inward", hjust = "inward") +
  geom_vline(xintercept= density_mean_sim, color="black") +
  coord_cartesian(clip = "off")


dev.off()


#### Where are localized TE relative to gene positions ? ---------------------------------


# First compute the intragenic TE density (TE in introns ...) 

gene_coords_name <- 
  gene_coords %>% 
  mutate(gene_name = paste(scaffold, start_gene, stop_gene, sep="-"))
gene_TE_joint <- 
  left_join(gene_coords_name, transposons_coords, by="scaffold")

Gene_TE_dens_df <-
  gene_TE_joint %>%
  mutate(TE_class = case_when(
  str_detect(transposon_name, "DNA") ~ "DNA",
  str_detect(transposon_name, "LINE") ~ "LINE",
  str_detect(transposon_name, "SINE") ~ "SINE",
  str_detect(transposon_name, "LTR") ~ "LTR",
  str_detect(transposon_name, "Helitron") ~ "DNA",
  str_detect(transposon_name, "Retroposon") ~ "Retroposon",
  str_detect(transposon_name, "Unspecified") ~ "Unspecified")) %>%
  filter((as.numeric(start_TE) >= as.numeric(start_gene) & as.numeric(stop_TE) <= as.numeric(stop_gene)) | 
           (as.numeric(stop_TE) >= as.numeric(start_gene) & as.numeric(stop_TE) <= as.numeric(stop_gene)) |
           (as.numeric(start_TE) >= as.numeric(start_gene) & as.numeric(start_TE) <= as.numeric(stop_gene)) |
           (as.numeric(start_TE) <= as.numeric(start_gene) & as.numeric(stop_TE) >= as.numeric(stop_gene))) %>%
  mutate(TE_new_start = if_else(
    start_TE < start_gene,
    start_gene,
    start_TE)) %>%
  mutate(TE_new_stop = if_else(
    stop_TE > stop_gene,
    stop_gene,
    stop_TE)) %>%
  mutate(TE_length = TE_new_stop - TE_new_start) %>%
  mutate(TE_density = TE_length/(stop_gene-start_gene)) %>%
  group_by(scaffold, start_gene, stop_gene, gene_name, TE_class) %>%
  dplyr::summarise(transposon_nb = n(), 
                   transposon_length = sum(TE_length), 
                   TE_density = sum(TE_density)
                   ) %>%
  ungroup()



gene_coords_name_renamed <- gene_coords_name
colnames(gene_coords_name_renamed) <- c("scaffold", "start_gene", "stop_gene", "cluster_name", "gene_name")
complete_gene_TE_joint <- 
  left_join(gene_coords_name_renamed, 
            Gene_TE_dens_df, by=c("gene_name"))


complete_gene_TE_joint$transposon_length[is.na(complete_gene_TE_joint$transposon_length)] <- 0
complete_gene_TE_joint$TE_density[is.na(complete_gene_TE_joint$TE_density)] <- 0

TE_class <- unique(complete_gene_TE_joint$TE_class)
gene_name <- unique(complete_gene_TE_joint$gene_name)
complete_gene_TE_joint_all <- 
  expand.grid(TE_class=TE_class, gene_name=gene_name, stringsAsFactors=FALSE) %>%
  left_join(complete_gene_TE_joint, by=c("gene_name", "TE_class")) %>% 
  dplyr::select(TE_class, gene_name, transposon_length, TE_density)  %>%
  filter(! is.na(TE_class))
complete_gene_TE_joint_all[is.na(complete_gene_TE_joint_all)] <- 0



Class_colors <- c(DNA = "#FFC107" , LINE = "#576FBB" , SINE = "#E351E3" , LTR = "#004D40" , Helitron = "#737D62" , Retroposon = "#A27058" , Unspecified = "gray", Total = "black")


pdf("./Intragenic_TE_density.pdf",width = 6.34,  height = 4.61)


complete_gene_TE_joint_all %>% 
  group_by(TE_class) %>%
  summarise(mean_density = mean(TE_density)) %>%
  bind_rows(summarise_all(., ~if(is.numeric(.)) sum(.) else "Total")) %>%
  ggplot(., aes(x=1, y=mean_density, color=TE_class)) +
  geom_point() +
  theme_classic() +
  ylab("TE density") +
  xlab("Intragenic")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_color_manual(values=Class_colors)

dev.off()

write.table(complete_gene_TE_joint_all, file="TE_table_intragenic.tsv", quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)

#Now lets do the same kind of thing but to see in upstream --  5kb 

gene_coords_named <- gene_coords %>% mutate(gene_name = paste(scaffold, start_gene, stop_gene, sep="-"))

TE_class_df <- 
  transposons_coords %>%
  mutate(TE_class = case_when(
    str_detect(transposon_name, "DNA") ~ "DNA",
    str_detect(transposon_name, "LINE") ~ "LINE",
    str_detect(transposon_name, "SINE") ~ "SINE",
    str_detect(transposon_name, "LTR") ~ "LTR",
    str_detect(transposon_name, "Helitron") ~ "DNA",
    str_detect(transposon_name, "Retroposon") ~ "Retroposon",
    str_detect(transposon_name, "Unspecified") ~ "Unspecified"
  )) %>%
  dplyr::select(TE_class) %>% distinct()

TE_class_df <- TE_class_df %>% filter(! is.na(TE_class))

downstream_df <- data.frame(NULL)
for(pos in 1:5001){
  
  #Define the start - XX kb position
  gene_coords_ext <- 
    gene_coords_named %>%
    mutate(ext_start = start_gene-pos) %>%
    dplyr::select(scaffold, ext_start)
  
  gene_TE_joint <- 
    left_join(gene_coords_ext, 
              transposons_coords, 
              by="scaffold")
  
  current_pos_df <- 
    as.data.frame(
    gene_TE_joint %>%
    mutate(TE_class = case_when(
      str_detect(transposon_name, "DNA") ~ "DNA",
      str_detect(transposon_name, "LINE") ~ "LINE",
      str_detect(transposon_name, "SINE") ~ "SINE",
      str_detect(transposon_name, "LTR") ~ "LTR",
      str_detect(transposon_name, "Helitron") ~ "DNA",
      str_detect(transposon_name, "Retroposon") ~ "Retroposon",
      str_detect(transposon_name, "Unspecified") ~ "Unspecified")) %>%
    filter(as.numeric(start_TE) <= as.numeric(ext_start)) %>%
    filter(as.numeric(stop_TE) >= as.numeric(ext_start)) %>%
    group_by(TE_class) %>%
    summarise(count = n()) %>%
    mutate(TE_density = count/nrow(gene_coords_ext)) %>%
    ungroup())


  current_pos_df <- left_join(TE_class_df, current_pos_df, by="TE_class") 
  current_pos_df[is.na(current_pos_df)] <- 0

  current_pos_df <- current_pos_df %>% mutate(position = pos)

  downstream_df <- rbind(downstream_df,current_pos_df)

}

Present_TE_classes <- 
  downstream_df %>% filter(TE_density != 0) %>% 
  group_by(TE_class) %>% 
  summarise(count = n())%>%
  pull(TE_class) %>%
  unique()



pdf("./Upstream_TE_density.pdf",width = 6.34,  height = 4.61)

downstream_df %>%
  filter(TE_class %in% Present_TE_classes) %>%
  ggplot(., aes(x=position, y=TE_density, color=TE_class)) +
  geom_line() +
  geom_point() +
  theme_classic() +
  geom_point(data = downstream_df %>% filter(TE_class %in% Present_TE_classes) %>% group_by(position) %>% summarise(Total = sum(TE_density)), mapping= aes(x=position, y=Total), color="black") +
  geom_line(data = downstream_df %>% filter(TE_class %in% Present_TE_classes) %>% group_by(position) %>% summarise(Total = sum(TE_density)), mapping= aes(x=position, y=Total), color="black") +
  scale_x_reverse() +
  xlab("BP Upstream") +
  ylab("TE density") +
  scale_color_manual(values=Class_colors)



dev.off()

write.table(downstream_df, file="TE_table_upstream.tsv", quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)


# Now lets try to see in downstream 5kb 


TE_class_df <- 
  transposons_coords %>%
  mutate(TE_class = case_when(
    str_detect(transposon_name, "DNA") ~ "DNA",
    str_detect(transposon_name, "LINE") ~ "LINE",
    str_detect(transposon_name, "SINE") ~ "SINE",
    str_detect(transposon_name, "LTR") ~ "LTR",
    str_detect(transposon_name, "Helitron") ~ "DNA",
    str_detect(transposon_name, "Retroposon") ~ "Retroposon",
    str_detect(transposon_name, "Unspecified") ~ "Unspecified"
  )) %>%
  dplyr::select(TE_class) %>% distinct()

TE_class_df <- TE_class_df %>% filter(! is.na(TE_class))

downstream_df <- data.frame(NULL)
for(pos in 1:5001){
  gene_coords_ext <- 
    gene_coords_named %>%
    mutate(ext_stop = stop_gene+pos) %>%
    dplyr::select(scaffold, ext_stop)
  
  gene_TE_joint <- left_join(gene_coords_ext, transposons_coords, by="scaffold")
  
  current_pos_df <- 
    as.data.frame(
      gene_TE_joint %>%
        mutate(TE_class = case_when(
          str_detect(transposon_name, "DNA") ~ "DNA",
          str_detect(transposon_name, "LINE") ~ "LINE",
          str_detect(transposon_name, "SINE") ~ "SINE",
          str_detect(transposon_name, "LTR") ~ "LTR",
          str_detect(transposon_name, "Helitron") ~ "DNA",
          str_detect(transposon_name, "Retroposon") ~ "Retroposon",
          str_detect(transposon_name, "Unspecified") ~ "Unspecified")) %>%
        filter(as.numeric(start_TE) <= as.numeric(ext_stop)) %>%
        filter(as.numeric(stop_TE) >= as.numeric(ext_stop)) %>%
        group_by(TE_class) %>%
        summarise(count = n()) %>%
        mutate(TE_density = count/nrow(gene_coords_ext)) %>%
        ungroup())


  current_pos_df <- left_join(TE_class_df, current_pos_df, by="TE_class") 
  current_pos_df[is.na(current_pos_df)] <- 0

  current_pos_df <- current_pos_df %>% mutate(position = pos)

  downstream_df <- rbind(downstream_df,current_pos_df)

}

Present_TE_classes <- 
  downstream_df %>% filter(TE_density != 0) %>% 
  group_by(TE_class) %>% 
  summarise(count = n())%>%
  pull(TE_class) %>%
  unique()


pdf("./Downstream_TE_density.pdf",width = 6.34,  height = 4.61)


downstream_df %>%
  filter(TE_class %in% Present_TE_classes) %>%
  ggplot(., aes(x=position, y=TE_density, color=TE_class)) +
  geom_line() +
  geom_point() +
  theme_classic() +
  xlab("BP Downstream") +
  ylab("TE density") +
  geom_point(data = downstream_df %>% filter(TE_class %in% Present_TE_classes) %>% group_by(position) %>% summarise(Total = sum(TE_density)), mapping= aes(x=position, y=Total), color="black") +
  geom_line(data = downstream_df %>% filter(TE_class %in% Present_TE_classes) %>% group_by(position) %>% summarise(Total = sum(TE_density)), mapping= aes(x=position, y=Total), color="black") +
  scale_color_manual(values=Class_colors)


dev.off()


write.table(downstream_df, file="TE_table_downstream.tsv", quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)

