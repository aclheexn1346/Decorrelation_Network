# 
source("libraries.R")
source("helperFunc.R")
source("single_cell_funcs.R")
library(dplyr)
targetgene <- readRDS("single_cell_data/sig_genes_log_val.rds")
sig_genes_log_val <- readRDS("single_cell_data/sig_genes_log_val_full.rds")
idx_goodgene <- apply(sig_genes_log_val, 1, function(x) sd(x)/abs(mean(x)) > 0.25) %>% which()
goodgene <- sig_genes_log_val[idx_goodgene,]
goodgene %>% dim()
set.seed(12)
othergenes <- sample_n(as.data.frame(goodgene[1:7000,]), 2000) %>% as.matrix() %>% t()
# use 2000 genes to cluster the 1018 cells


# default clustering ------------------------------------------------------
sc_block_idx_full <- readRDS("single_cell_data/single_cell_block_idx_full.rds")
for(i in 1:length(sc_block_idx_full)){
  names(sc_block_idx_full[[i]]) <- colnames(targetgene)[sc_block_idx_full[[i]]]
}
# Xp <- readRDS(file = "data/single_cell_data/sig_genes_log_val.rds")
Xp <- t(targetgene)
Xp %>% dim()
# randomly sample 20 cells from each cell type and merge the indices

full_log_vals = Xp

# Everything below here was used in my code but you probably won't need it.
# targetgene will be the 51 genes and othergenes are the background genes
# # clustering --------------------------------------------------------------
d <- dist(scale(othergenes), method = "euclidean")
hc1 <- hclust(d, method = "complete")
plot(hc1, cex = 0.6, cex.axis = 1.4, cex.lab = 1.4, cex.main = 1.4)



sub_grp <- cutree(hc1, h=70)
sub_grp %>% table()
length(which(table(sub_grp) > 15))
sub_grp_subset <- sub_grp[!(sub_grp %in%  which(table(sub_grp) ==1))]

#plot(hc1)
clusters = unique(sub_grp_subset)
clusters = clusters[-c(142,143,145,146,147, 153,158)] # Has HFF which we want to remove.

#cluster_size = 15 to 30 sampling from full data
tot_sample = c()
block_sizes = c()
for(i in clusters){
  cluster_sub_grp = sub_grp_subset[which(sub_grp_subset == i)]
  if(length(cluster_sub_grp) > 15){
    if(length(cluster_sub_grp) > 30){
      cluster_size = 30
    }
    else{
      cluster_size = length(cluster_sub_grp)
    }
    cluster_sample = sample(cluster_sub_grp, cluster_size)
    block_sizes = c(block_sizes, cluster_size)
    tot_sample = c(tot_sample, cluster_sample)
  }
}

# discretizing the data
discretize_data = function(df){
  new_df = df
  breakpoints = c()
  for(i in 1:ncol(df)){
    breakpoint = arules::discretize(df[,i], method = "cluster", breaks = 2)
    new_df[,i] = as.numeric(breakpoint)
    new_df[,i][new_df[,i] == 1] = 0
    new_df[,i][new_df[,i] == 2] = 1
  }
  return(list(new_df, breakpoints))
}

disc_full_log_vals = discretize_data(full_log_vals)
# full log vals is the 51 target genes and then 1018 is the number of cells we have


obtain_data_from_sample = function(tot_sample, orig_data_values){
  extracted_cell_names = names(tot_sample)
  extracted_samples = orig_data_values[extracted_cell_names,]
  df = data.frame(extracted_samples)
  df$cell_name = extracted_cell_names
  df$cluster_num = tot_sample
  df = df[c(52, 53, 1:51)]
  return(df)
}

df = obtain_data_from_sample(tot_sample, disc_full_log_vals[[1]])

save(df, file = "sub_sample_discrete.RData")

# sample 100 of 2000 othergenes to get covariance estimate

# estimating covariance from background genes
background_gene_cov_estimation = function(bgr_gene_df, n = 100, tot_sample, block_sizes){
  remove_columns = c() # Need to have a decent number of factors
  for(j in 1:ncol(bgr_gene_df)){
    if(length(levels(as.factor(bgr_gene_df[,j]))) < 10){
      remove_columns = c(remove_columns, j)
    }
  }
  sample_from_bgr = setdiff(1:ncol(bgr_gene_df), remove_columns)
  sampled = sample(sample_from_bgr, size = n)
  extracted_cell_names = names(tot_sample)
  bgr_extracted_samples = bgr_gene_df[extracted_cell_names,sampled]
  bgr_data = discretize_data(bgr_extracted_samples)[[1]]
  init_epsilon = matrix(0, nrow = nrow(bgr_data), ncol = ncol(bgr_data))
  remove_column_indices = c() # remove any columns with less than two factors since every data is the same
  for(i in 1:ncol(bgr_data)){
    if(length(levels(as.factor(bgr_data[,i]))) < 2){
      remove_column_indices = c(remove_column_indices, i)
    }
  }
  if(length(remove_column_indices) > 0){
    bgr_data = bgr_data[,-remove_column_indices]
  }
  bgr_data = as.matrix(bgr_data)
  rownames(bgr_data) = 1:dim(bgr_data)[1]
  colnames(bgr_data) = 1:dim(bgr_data)[2]

  init_beta = init_data_dag(data = bgr_data)
  init_beta[which(init_beta != 0)] = 0.0001

  init_epsilon = matrix(0, nrow = nrow(bgr_data), ncol = ncol(bgr_data))

  # Getting estimated covariance and the uncorrelated data for the concensus cpdag estimate
  beta_est_ident_lr <- beta_est_loop_double_cov_est(data = as.matrix(bgr_data), init_beta = init_beta, init_epsilon = init_epsilon, block_sizes = block_sizes, loops = 21, trueB = init_beta)
  est_Sigma = beta_est_ident_lr$Sigma

  return(est_Sigma)
}

bgr_est_Sigma = background_gene_cov_estimation(bgr_gene_df = othergenes, n = 100, tot_sample = tot_sample, block_sizes = block_sizes)

# Data preparation
data = df[,c(-1,-2)]
data = as.matrix(data)
rownames(data) = 1:dim(data)[1]
colnames(data) = 1:dim(data)[2]
init_epsilon = matrix(0, nrow = nrow(data), ncol = ncol(data))

# Cross validation likelihood by cluster

CV_likelihood_all = function(data, block_sizes, bgr_Estimated_Sigma){

  likelihood_folds_init_2 = rep(0,length(block_sizes))
  likelihood_folds_concensus_2 = rep(0, length(block_sizes))
  likelihood_folds_concensus_nocov_2 = rep(0, length(block_sizes))
  likelihood_folds_complete_nocov = rep(0, length(block_sizes))
  likelihood_folds_complete_cov = rep(0, length(block_sizes))
  
  for(cluster_number in 1:length(block_sizes)){
    print(cluster_number)
    init_index = (sum(block_sizes[0:(cluster_number-1)]))+1
    end_index = sum(block_sizes[1:cluster_number])
    test_data = data[init_index:end_index,]
    train_data = data[-(init_index:end_index),]
    remove_column_indices = c() # remove any columns with less than two factors since every data is the same
    for(i in 1:ncol(train_data)){
      if(length(levels(as.factor(train_data[,i]))) < 2){
        remove_column_indices = c(remove_column_indices, i)
      }
    }
    if(length(remove_column_indices) > 0){
      train_data = train_data[,-remove_column_indices]
      test_data = test_data[,-remove_column_indices]
    }
    train_data = as.matrix(train_data)
    rownames(train_data) = 1:dim(train_data)[1]
    colnames(train_data) = 1:dim(train_data)[2]

    init_beta = init_data_dag(data = train_data)
    init_beta[which(init_beta != 0)] = 0.0001

    init_epsilon = matrix(0, nrow = nrow(train_data), ncol = ncol(train_data))

    # Getting estimated covariance and the uncorrelated data for the concensus cpdag estimate
    beta_est_ident_lr <- beta_est_loop4(data = as.matrix(train_data), init_beta = init_beta, init_epsilon = init_epsilon, block_sizes = block_sizes[-cluster_number], loops = 21, trueB = init_beta)

    # getting the final cpdag from the initial data
    init_dag = init_data_dag(train_data)
    init_g1 <- as(t(init_dag), "graphNEL")
    init_cpdag <- dag2cpdag(init_g1)
    init_adj_mat = (as(init_cpdag, "matrix"))

    # Getting the final cpdag from the concensus data
    concensus_data = beta_est_ident_lr$hold_data
    concensus_mat = matrix(0, nrow = ncol(train_data), ncol = ncol(train_data))
    sampled_datasets = sample(2:21, size = 10)
    for(i in sampled_datasets){
      #print(i)
      parents = neighborhood_selection_mmhc_get_dag_decor(concensus_data[[i]])
      parents = dag2cpdag(parents)
      concensus_mat = concensus_mat + parents
    }
    concensus_mat[which(concensus_mat <= 4)] = 0
    concensus_mat[which(concensus_mat >= 5)] = 1
    concensus_g1 <- as(t(concensus_mat), "graphNEL")
    concensus_cpdag = dag2cpdag(concensus_g1)
    concensus_adj_mat = (as(concensus_cpdag, "matrix"))

    # Getting the final beta using the final cpdag structure obtained
    # First is using the final cpdag from the concensus and using the estimated covariance to make one last beta estimate
    concensus_adj_mat[concensus_adj_mat != 0] = 0.0001
    conc_final_beta = beta_est_loop_input_covariance(data = train_data, adj_mat = concensus_adj_mat, cov_est = beta_est_ident_lr$Sigma, init_epsilon = init_epsilon, block_sizes = block_sizes[-cluster_number])[[1]]

    # Using the initial cpdag and assumption of identity matrix but obtaining the final beta using same method as above
    init_final_beta = beta_est_loop_input_covariance(data = train_data, adj_mat = init_beta, cov_est = diag(nrow(train_data)), init_epsilon = init_epsilon, block_sizes = block_sizes[-cluster_number])[[1]]
    # Under complete graph
    complete_beta = matrix(1, nrow = nrow(init_beta), ncol = ncol(init_beta)) # put adjacency matrix as 1 for everything
    complete_final_beta = beta_est_loop_input_covariance(data = train_data, adj_mat = complete_beta, cov_est = diag(nrow(train_data)), init_epsilon = init_epsilon, block_sizes = block_sizes[-cluster_number])[[1]]

    # Finding the log likelihood on the test data for cross_validation
    test_est_Sig = Sig_Estimate_DAG_test(X = test_data, beta = conc_final_beta, block_sizes = c(nrow(test_data)))
    conc_ll2 = t_normal_likelihood(testdata = test_data, Sigma_estimate = bgr_Estimated_Sigma[init_index:end_index, init_index:end_index], est_beta = conc_final_beta)
    conc_nocov_ll2 = t_normal_likelihood(testdata = test_data, Sigma_estimate = diag(nrow(test_data)), est_beta = conc_final_beta)
    init_ll2 = t_normal_likelihood(testdata = test_data, Sigma_estimate = diag(nrow(test_data)), est_beta = init_final_beta)
    complete_ll2 = t_normal_likelihood(testdata = test_data, Sigma_estimate = diag(nrow(test_data)), est_beta = complete_final_beta)
    complete_estsig_ll2 = t_normal_likelihood(testdata = test_data, Sigma_estimate = bgr_Estimated_Sigma[init_index:end_index, init_index:end_index], est_beta = complete_final_beta)
    
    
    likelihood_folds_concensus_2[cluster_number] = conc_ll2[1]
    likelihood_folds_concensus_nocov_2[cluster_number] = conc_nocov_ll2[1]
    likelihood_folds_init_2[cluster_number] = init_ll2[1]
    likelihood_folds_complete_nocov[cluster_number] = complete_ll2[1]
    likelihood_folds_complete_cov[cluster_number] = complete_estsig_ll2[1]
    

  }
  return(list(
              "init_mvt" = likelihood_folds_init_2,
              "consensus_mvt" = likelihood_folds_concensus_2,
              "consensus_nocov_mvt" = likelihood_folds_concensus_nocov_2,
              "complete_nocov_mvt" = likelihood_folds_complete_nocov,
              "complete_cov_mvt" = likelihood_folds_complete_cov,
              "Final_consensus_CPDAG" = concensus_adj_mat))
}

CV_res = CV_likelihood_all(data = data, block_sizes = block_sizes, bgr_Estimated_Sigma = bgr_est_Sigma)

#mvt_data = data.frame(complete_cov = CV_res$complete_cov_mvt, consensus = CV_res$consensus_mvt)
#boxplot(mvt_data, xlab = "Method", ylab = "Log-Likelihood", main = "Mvt Normal log-likelihood for each method on each fold")
#est_cpdag = CV_res$Final_consensus_CPDAG
#colnames(est_cpdag) = rownames(est_cpdag) = colnames(Xp)
#est_cpdag[which(est_cpdag != 0)] = 1
#dat = struc_diff_1_mixed$uncor_Z_list[[9]][[9]][[10]]

#row_centered <- t(apply(dat, 1, function(row) row - mean(row)))
#row_normalized <- t(apply(row_centered, 1, function(row) row / sd(row)))
# Compute the row-covariance matrix
#row_cov_matrix <- (row_normalized %*% t(row_normalized)) / (ncol(dat) - 1)

#true_sig = struc_diff_1_mixed$true_sigma[[9]]


#data3 = c(data3, true_sig[true_sig != 0])
#data4 = c(data4, row_cov_matrix[true_sig != 0])
# breaks <- pretty(c(data1, data2), n = 20)  # Define common breaks
# hist3 <- hist(data3, breaks = breaks, plot = FALSE)
# hist4 <- hist(data4, breaks = breaks, plot = FALSE)
# 
# 
# hist1 <- hist(data1, breaks = breaks, plot = FALSE)
# hist2 <- hist(data2, breaks = breaks, plot = FALSE)
# 
# # Find the maximum frequency for scaling
# max_freq <- max(c(hist1$counts, hist2$counts))
# 
# # Plot the first histogram
# plot(hist3, col = rgb(1, 0, 0, 0.5), ylim = c(0, max_freq), xlim = range(breaks), main = "Combined Histogram", xlab = "Value", ylab = "Frequency")
# plot(hist4, col = rgb(0, 0, 1, 0.5), add = TRUE)
# legend("topright", legend = c("Row Sample Covariance of De-correlated Z", "True Sigma Covariance"), fill = c(rgb(1, 0, 0, 0.5), rgb(0, 0, 1, 0.5)))
# 
# 
# plot(hist1, col = rgb(1, 0, 0, 0.5), ylim = c(0, max_freq), xlim = range(breaks), main = "Correlation before and after De-correlation", xlab = "Correlation", ylab = "Frequency")
# 
# # Add the second histogram
# plot(hist2, col = rgb(0, 0, 1, 0.5), add = TRUE)
# 
# # Add a legend
# legend("topright", legend = c("True Sigma Covariance","Row Sample Correlation of De-correlated Z"), fill = c(rgb(1, 0, 0, 0.5), rgb(0, 0, 1, 0.5)))
