install.packages("dplyr")
options(repos = c(CRAN = "https://cloud.r-project.org"))
install.packages('biclust')
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("bioDist")

install.packages('dendextend')
library(dendextend)
install.packages('ape')
library(ape)

install.packages(c("htmlwidgets", "png", "base64enc"))

BiocManager::install('ComplexHeatmap', force = TRUE)
BiocManager::install('InteractiveComplexHeatmap', force = TRUE)
BiocManager::install('cola', force = TRUE)


install.packages('e1071')
library('e1071')
install.packages("Cairo")
library(Cairo)
library(cola)
library(ComplexHeatmap)
library(InteractiveComplexHeatmap)
install.packages('/Users/gautham/microbiome_analysis/d3heatmap_0.6.0.tar.gz', repos = NULL, type = 'source')
library(d3heatmap)
install.packages('dynamicTreeCut')
library(dynamicTreeCut)
install.packages("withr")
library(biclust)
library(dplyr)
install.packages("pheatmap")
install.packages("gridExtra")
library(pheatmap)
library(gridExtra)
install.packages("infotheo")
library(infotheo)
install.packages('proxy')
library(proxy)
install.packages("entropy")
library(entropy)

install.packages("Information")
library(Information)

install.packages('dendroextras')
library(dendroextras)

df <- test_array_output
selected_columns <- df %>%
  select(contains("GCF"))

tcid_acc_cols <- df[c("TCID", "Acc")]
tc_list <- list()

for (i in 1:nrow(tcid_acc_cols)){
  
  tcid <- tcid_acc_cols[i, 1]
  acc <- tcid_acc_cols[i, 2]
  
  tc_acc <- paste(tcid, acc, sep = "-")
  tc_list[[i]] <- tc_acc
}


genome_matrix_df <- selected_columns %>% 
  mutate_all(~as.numeric(ifelse(. == '+', 1, ifelse(. == '-', 0, .))))

fam_df <- family_matrix
selected_genomes <- fam_df %>%
  select(contains("GCF"))

family_matrix_df <- selected_genomes %>% 
  mutate_all(~as.numeric(ifelse(. == '+', 1, ifelse(. == '-', 0, .))))

family_matrix_df <- family_matrix[,-1]
rownames(family_matrix_df) <- family_matrix[,1]

fam_percent_df <- family_percent[,-1]

rownames(fam_percent_df) <- family_percent[,1]

rownames(genome_matrix_df) <- tc_list
genome_matrix <- as.matrix(genome_matrix_df)
fam_matrix <- as.matrix(family_matrix_df)
percent_matrix <- as.matrix(fam_percent_df)

res <- biclust(x=genome_matrix, method=BCBimax(), minr=4, minc=4, number=17)

# this is using manhattan distance measurement as default for both rows and columns using kmeans clustering
#original_heatmap <- pheatmap(genome_matrix, cluster_rows = TRUE, cluster_cols = TRUE, scale='none', main = 'Original heatmap', kmeans_k = 17, clustering_distance_rows = 'manhattan', clustering_distance_cols = 'manhattan')
# temp <- proxy::dist(genome_matrix, method = 'Manhattan')
# hc <- hclust(temp, method = 'complete')
# plot(hc)
# dendrogram <- hclust(dist(genome_matrix, method = "manhattan"))
# plot(dendrogram)
# 
# row_clusters <- cutree(dendrogram, k=17)
# col_clusters <- cutree(dendrogram, k=17)
# 
# new_data <- genome_matrix[row_clusters, col_clusters]
# heatmap <- pheatmap(new_data, clustering_distance_rows = "manhattan",
#          clustering_distance_cols = "manhattan",
#          cluster_rows = TRUE, cluster_cols = TRUE,
#          main = "Hierarchical Clustering Heatmap (Manhattan Distance)")

original_heatmap <- pheatmap(genome_matrix, cluster_rows = TRUE, cluster_cols = TRUE, scale='none', main = 'Original heatmap', clustering_method = "complete", clustering_distance_rows = 'manhattan', clustering_distance_cols = 'manhattan')

row_dendogran <- original_heatmap$tree_row
row_groups <- cutree(row_dendogran, h = 1)
plot(row_dendogran)
rect.hclust(row_dendogran, k=17, border = 1:17)
protein_clusters = cutree(row_dendogran, k=17)

col_dendogram <- original_heatmap$tree_col
col_groups <- cutree(col_dendogram, h = 1)
plot(col_dendogram)
rect.hclust(col_dendogram, k=6, border = 1:17)
genome_clusters = cutree(col_dendogram, k=17)

# row_kmeans_clusters <- kmeans(genome_matrix, centers = 17)$cluster
# col_kmeans_clusters <- kmeans(t(genome_matrix), centers = 17)$cluster

# Create dataframes for row and column clusters
# row_clusters_df <- data.frame(Row = rownames(genome_matrix), RowCluster = row_kmeans_clusters)
# col_clusters_df <- data.frame(Column = colnames(genome_matrix), ColCluster = col_kmeans_clusters)

mi_matrix <- matrix(0, ncol(genome_matrix), ncol(genome_matrix))
# for (i in 1:(ncol(genome_matrix) - 1)) {
#   for (j in (i + 1):ncol(genome_matrix)) {
#     mi_value <- mutinformation(genome_matrix[, i], genome_matrix[, j])
#     mi_matrix[i, j] <- mi_value
#     mi_matrix[j, i] <- mi_value
#   }
# }
transpose_percent_matrix <- as.data.frame(t(fam_percent_df))
transpose_genome_matrix <- as.data.frame(t(genome_matrix_df))
transpose_fam_matrix <- as.data.frame(t(family_matrix_df))
                       
mi_matrix <- mutinformation(transpose_genome_matrix, method='emp') # big runtime
mi_matrix_genomes <- mutinformation(genome_matrix_df, method='emp')


mi_matrix_fams <- mutinformation(transpose_fam_matrix, method='emp') # big runtime
mi_matrix_fam_genomes <- mutinformation(family_matrix_df, method='emp')

mi_matrix_perc_genomes <- mutinformation(fam_percent_df, method='emp')
mi_matrix_percents <- mutinformation(transpose_percent_matrix, method='emp') # big runtime

# Convert mutual information to a distance measure (e.g., using 1 - MI)
distance_matrix_protein <- 1 - mi_matrix
distance_matrix_genome <- 1 - mi_matrix_genomes

distance_matrix_fam <- 1 - mi_matrix_fams
distance_matrix_fam_genome <- 1 - mi_matrix_fam_genomes

dist_object_row <- as.dist(distance_matrix_protein)
dist_object_famrow <- as.dist(distance_matrix_fam)

dist_object_col <- as.dist(distance_matrix_genome)
dist_object_famcol <- as.dist(distance_matrix_fam_genome)

custom_heatmap <- pheatmap(
  genome_matrix,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  scale = 'none',
  main = 'Using mutual information',
  clustering_method = "complete",
  clustering_distance_rows = dist_object_row,
  clustering_distance_cols = dist_object_col
)


percent_heatmap <- pheatmap(
  percent_matrix,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  scale = 'none',
  main = 'Using mutual information',
  clustering_method = "complete",
  clustering_distance_rows = 'euclidean',
  clustering_distance_cols = 'euclidean'
)

hamming_mat <- hamming.distance(fam_matrix)
hamming_dist_fams <- as.dist(hamming_mat)

hamming_mat_genomes <- hamming.distance(as.matrix(transpose_fam_matrix))
#hamming_mat_genomes <- proxy::dist(as.matrix(transpose_fam_matrix), method = "hamming")

hamming_dist_genomes <- as.dist(hamming_mat_genomes)

fam_heatmap <- Heatmap(fam_matrix,
                       name = "test",
                       clustering_distance_rows = 'manhattan',
                       clustering_distance_columns = 'manhattan')

mut_info_fam <- Heatmap(fam_matrix,
                         name = "test",
                          clustering_method_rows = 'ward.D2',
                          clustering_method_columns = 'ward.D2',
                         clustering_distance_rows = dist_object_famrow,
                         clustering_distance_columns = dist_object_famcol)

percent_heatmap <- Heatmap(percent_matrix,
                       name = "test",
                       clustering_distance_rows = 'euclidean',
                       clustering_distance_columns = 'euclidean')

inter_heatmap <- Heatmap(genome_matrix,
                   name = "test",
                   clustering_distance_rows = dist_object_row,
                   clustering_distance_columns = dist_object_col)

hamming_heatmap <- Heatmap(fam_matrix,
                           name = "test",
                           clustering_method_rows = 'ward.D2',
                          clustering_method_columns = 'ward.D2',
                           clustering_distance_rows = hamming_dist_fams,
                           clustering_distance_columns = hamming_dist_genomes)


pheatmap_hamming <- pheatmap(fam_matrix,
                            name = "test",
                            clustering_distance_rows = hamming_dist_fams,
                            clustering_distance_columns = hamming_dist_genomes)

# Draw the Heatmap
mut_info_fam <- draw(mut_info_fam)
htShiny(mut_info_fam)

mut_infotree <- as.phylo(row_dend(mut_info_fam))
write.tree(mut_infotree, file='mutinfo_tree.tree')

hamming_colclust <- as.phylo(column_dend(mut_info_fam))
write.tree(hamming_colclust, file='hamming_genomes.tree')

mut_colclust <- as.phylo(column_dend(mut_info_fam))
write.tree(mut_colclust, file='mut_genomes.tree')

hamming_heatmap <- draw(hamming_heatmap)
htShiny(hamming_heatmap)

hamming_tree <- as.phylo(row_dend(hamming_heatmap))
write.tree(hamming_tree, file='hamming_rowtree.tree')



inter_heatmap <- draw(inter_heatmap)
htShiny(inter_heatmap)

fam_heatmap <- draw(fam_heatmap)
htShiny(fam_heatmap)

row_Ded <- row_dend(fam_heatmap)
man_tree <- as.phylo(row_Ded)
plot(man_tree)

write.tree(man_tree, file='man_tree.tree')

percent_heatmap <- draw(percent_heatmap)
htShiny(percent_heatmap)
# # test <- Heatmap(genome_matrix)
# # test <- draw(test)
# 
#htShiny(test)

genome_matrix.clust <- cbind(genome_matrix, cluster = cutree(custom_heatmap$tree_row, k = 17))

for (i in 1:17) {
  subset_data <- genome_matrix.clust[genome_matrix.clust[,"cluster"] == i, -ncol(genome_matrix.clust)]
  pheatmap(
    subset_data,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    scale = 'none',
    main = paste("Cluster", i),
    clustering_method = "complete",
    clustering_distance_rows = dist_object_row,
    clustering_distance_cols = dist_object_col
  )
  
}

my_tree <- as.phylo(custom_heatmap$tree_row)

write.tree(phy=my_tree, file="file_name.tree")  
cd <- cutreeDynamic(dendro = custom_heatmap$tree_row, minClusterSize = 20)

row_dendrogram <- custom_heatmap$tree_row
dend <- as.dendrogram(row_dendrogram)

for (i in seq_along(cd)) {
  # Color the branches of the dendrogram based on the current cluster assignment
  dend_colored <- color_branches(dend, k = cd[[i]])
  
  # Plot the dendrogram with colored clusters
  plot(dend_colored, main = paste("Cluster", i))
  break
}



row_dendogran <- custom_heatmap$tree_row
row_groups <- cutree(row_dendogran, h = 0.9)
plot(row_dendogran)
rect.hclust(row_dendogran, k=6, border = 1:17)
protein_clusters = cutree(row_dendogran, k=17)

col_dendogram <- custom_heatmap$tree_col
col_groups <- cutree(col_dendogram, h = 0.9)
plot(col_dendogram)
rect.hclust(col_dendogram, k=6, border = 1:17)
genome_clusters = cutree(col_dendogram, k=6)

# cluster_heatmaps <- list()
# 
# for (i in 1:max(protein_clusters)) {
#   for (j in 1:max(genome_clusters)) {
#     # Identify the rows and columns belonging to the current cluster
#     rows_in_cluster <- which(protein_clusters == i)
#     cols_in_cluster <- which(genome_clusters == j)
#     
#     # Subset the data for the current cluster
#     subset_data <- genome_matrix[rows_in_cluster, cols_in_cluster]
#     
#     # Create a heatmap for the current cluster
#     heatmap <- pheatmap(
#       subset_data,
#       scale = 'none',
#       main = paste("Row Cluster", i, " - ", "Column Cluster", j),
#       clustering_method = "complete"
#     )
#     
#     row_dendogran <- heatmap$tree_row
#     row_groups <- cutree(row_dendogran, h = 1)
#     plot(row_dendogran)
#     rect.hclust(row_dendogran, k=6, border = 1:17)
#     temp_p = cutree(row_dendogran, k=6)
#     
#     col_dendogram <- heatmap$tree_col
#     col_groups <- cutree(col_dendogram, h = 4)
#     plot(col_dendogram)
#     rect.hclust(col_dendogram, k=2, border = 1:17)
#     temp_g = cutree(col_dendogram, k=6)
#     # Store the heatmap in the list
#     cluster_heatmaps[[paste("RowCluster_", i, "_ColumnCluster_", j)]] <- heatmap
#   }
# }

cluster_heatmaps <- list()

# Iterate through the clusters and create individual heatmaps
for (i in 1:max(protein_clusters)) {
  for (j in 1:max(genome_clusters)) {
    # Subset the data for the current cluster
    subset_data <- genome_matrix[row_groups == i & col_groups == j, ]
    
    if (is.null(subset_data) || length(subset_data) == 0) {
      cat("subset_data is empty or contains no data.\n")
    } 
    else {
      heatmap <- pheatmap(
        subset_data,
        cluster_rows = FALSE,  # Avoid re-clustering within the cluster
        cluster_cols = FALSE,  # Avoid re-clustering within the cluster
        scale = 'none',
        main = paste("Row Cluster", i, " - ", "Column Cluster", j),
        clustering_method = "none"  # No clustering within the cluster
      )
      
      # Store the heatmap in the list
      cluster_heatmaps[[paste("RowCluster_", i, "_ColumnCluster_", j)]] <- heatmap
      
    }
    # Create a heatmap for the current cluster
      
    }
}

# test <- biclust(x=genome_matrix, method=BCBimax(), minr=4, minc=4, number=17, dist = list(row = dist_object_row, col = dist_object_col))
