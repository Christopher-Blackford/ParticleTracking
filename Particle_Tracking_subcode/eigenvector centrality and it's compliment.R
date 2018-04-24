#Graph
con_matrix <- matrix(data = sample(1:25, size = 25), nrow = 5, ncol = 5, byrow = FALSE, dimnames = NULL)

right_eigenvalue <- eigen(con_matrix)

left_eigenvalue <- eigen(t(con_matrix))

left_eigenvalue


con_matrix <- eigen_centrality(con_matrix, directed=TRUE, weights=E(con_matrix)$Freq)$vector
df_Eigen_Centrality <- as.data.frame(df_Eigen_Centrality)
df_Eigen_Centrality$Poly_ID <- row.names(df_Eigen_Centrality)
row.names(df_Eigen_Centrality) <- 1:nrow(df_Eigen_Centrality)
df_Eigen_Centrality <- dplyr::rename(df_Eigen_Centrality, Eigenvector_centrality = df_Eigen_Centrality)



con_matrix <- eigen_centrality(con_matrix, directed=TRUE, weights=E(con_matrix)$Freq)
con_matrix
