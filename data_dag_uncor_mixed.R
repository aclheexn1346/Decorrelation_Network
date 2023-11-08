#Using MMHC to find dag structure first to see if there's improvement
# and use other method to find difference in structure improvement
source("helperFunc.R")
source("libraries.R")
decor_struct_mmhc = function(n, p, cov_struc, reps, vers){
  #sig_2 = block_diag_sep(n = n, block_size = 2, struc_matrix = matrix(c(1,rho,rho,1), nrow = 2, ncol = 2, byrow = T))
  #block_obj = block_diag_sep_var(n = n, min_block_size = 5, max_block_size = 10, cov_struc = cov_struc)
  #sig_2 = block_obj[[1]]
  #block_sizes = block_obj[[2]]
  omg.sq = rep(1, p)
  # f1_score_pre = numeric(reps)
  # f1_score_post = numeric(reps)
  sig_2_list = list()
  hold_uncor_datasets_list = list()
  hold_uncor_datasets_list_lr = list()
  est_betas_list = list()
  
  hold_true_uncor_datasets_list = list()
  data_list = list()
  Beta_list = list()
  est_sigma_list = list()
  update_beta_list_lr = list()
  update_beta_list_bayes = list()
  
  between_beta_list_lr = list()
  between_beta_list_bayes = list()
  
  error_TrueB_list_lr = list()
  error_TrueB_list_bayes = list()
  
  mse_list_lr = list()
  mse_list_bayes = list()
  
  for(i in 1:reps){
    print(i)
    block_obj = block_diag_sep_var(n = n, min_block_size = 5, max_block_size = 10, cov_struc = cov_struc)
    sig_2 = block_obj[[1]]
    block_sizes = block_obj[[2]]
    sig_2_list[[i]] = sig_2
    hold_uncor_datasets = list()
    hold_uncor_datasets = list()
    hold_true_uncor_datasets = list()
    vers = vers+1
    B = gen.B(p = p, seed = vers)
    sim_result_nplarge = sim_X_LUM(vers, p, n, omg.sq, sig = sig_2, b = B$b)
    
    data = sim_result_nplarge$X[1:n,]
    data_prop = numeric(ncol(data))
    for(j in 1:ncol(data)){
      data_prop[j] = length(which(data[,j] == 1))/nrow(data)
    }
    while(min(data_prop) == 0 | max(data_prop == 1)){
      print("data prop is all 1")
      vers = vers+1
      sim_result_nplarge = sim_X_LUM(vers, p, n, omg.sq, sig = sig_2, b = B$b)
      data = sim_result_nplarge$X
      data_prop = numeric(ncol(data))
      for(k in 1:ncol(data)){
        data_prop[k] = length(which(data[,k] == 0))/nrow(data)
      }
    }
    true_epsilon = sim_result_nplarge$eps_mat
    parents = neighborhood_selection_mmhc_get_dag(data = data)
    # Get F1 score with respect to parents
    # f1_score_discrete = get_F1_dag_parents(parents, B$b)
    init_beta = init_data_dag(data = data)
    true_beta = B$b
    #true_beta[which(true_beta != 0)] = 0.001
    init_epsilon = matrix(0, nrow = n, ncol = p)
    beta_est_ident_b0 = NULL
    while(is.null(beta_est_ident_b0)){
      try({
        #beta_est_p6_10_blockdraw <- beta_est_loop(data = data, init_beta = init_beta, init_epsilon = init_epsilon, block_size = 5, loops = 10, trueB = B$b)
        beta_est_ident_b0 <- beta_est_loop_bayesian_draws_ridge(data = data, init_beta = init_beta, init_epsilon = init_epsilon, block_sizes = block_obj[[2]], loops = 21, trueB = B$b)
        beta_est_ident_lr <- beta_est_loop4(data = data, init_beta = init_beta, init_epsilon = init_epsilon, block_sizes = block_obj[[2]], loops = 21, trueB = B$b)
        
        
        #beta_est_ident_b0_true <- beta_est_loop_decor_z(data = data, init_beta = true_beta, init_epsilon = init_epsilon, block_sizes = block_obj[[2]], loops = 21, trueB = B$b)
        #beta_est_ident_lasso <- beta_est_loop5(data = data, init_beta = init_beta, init_epsilon = init_epsilon, block_sizes = block_obj[[2]], loops = 21, trueB = B$b)
        
        #beta_est_p6_10_blockdraw_old <- beta_est_loop_old(data = data, init_beta = init_beta, init_epsilon = init_epsilon, block_size = 2, loops = 100, trueB = B$b)
        # beta_est_blockdraw_truesig = beta_est_loop_truesig(data = data, init_beta = init_beta, init_epsilon = init_epsilon, block_size = 2, loops = 100, trueB = B$b, truesig = sig_2)
      })
      if(is.null(beta_est_ident_b0)){
        vers = vers+1
        sim_result_nplarge = sim_X_LUM(vers, p, n, omg.sq, sig = sig_2, b = B$b)
        data = sim_result_nplarge$X
        true_epsilon = sim_result_nplarge$eps_mat
      }
    }
    # large_b = which(abs(beta_est_p6_10_blockdraw$beta) > 3)
    # indices = unique(round_any(large_b, 100))/100
    # data_proportions = c(data_proportions, data_prop[indices])
    #uncor_data = beta_est_ident_lasso$hold_data
    # hold_uncor_datasets[[i]] = uncor_data
    # true_uncor_data = beta_est_ident_b0_true$Decor_data
    # hold_true_uncor_datasets[[i]] = true_uncor_data
    hold_uncor_bayes = beta_est_ident_b0$hold_data
    hold_uncor_lr = beta_est_ident_lr$hold_data
    
    est_betas_lr = beta_est_ident_lr$est_betas
    est_betas_bayes = beta_est_ident_b0$est_betas
    
    diff_betas_lr = beta_est_ident_lr$between
    diff_betas_bayes = beta_est_ident_b0$between
    
    error_trueB_lr = beta_est_ident_lr$Error_trueB
    error_trueB_bayes = beta_est_ident_b0$Error_trueB
    
    mse_lr = beta_est_ident_lr$Error
    mse_bayes = beta_est_ident_b0$Error
    
    update_beta_list_lr[[i]] = est_betas_lr
    update_beta_list_bayes[[i]] = est_betas_bayes
    
    between_beta_list_lr[[i]] = diff_betas_lr
    between_beta_list_bayes[[i]] = diff_betas_bayes
    
    error_TrueB_list_lr[[i]] = error_trueB_lr
    error_TrueB_list_bayes[[i]] = error_trueB_bayes
    
    mse_list_lr[[i]] = mse_lr
    mse_list_bayes[[i]] = mse_bayes
    

    hold_uncor_datasets_list[[i]] = hold_uncor_bayes
    hold_uncor_datasets_list_lr[[i]] = hold_uncor_lr
    #hold_true_uncor_datasets_list[[i]] = hold_true_uncor_datasets
    data_list[[i]] = data
    Beta_list[[i]] = B$b
    est_sigma_list[[i]] = beta_est_ident_b0$Sigma
  }
  return(list("Data" = data_list, 
              "Beta" = Beta_list,
              "est_betas" = est_betas_list,
              "est_Sigmas" = est_sigma_list,
              "true_sigma" = sig_2_list,
              "uncor_Z_list" = hold_uncor_datasets_list,
              "uncor_Z_list_lr" = hold_uncor_datasets_list_lr,
              "update_beta_list_lr" = update_beta_list_lr,
              "update_beta_list_bayes" = update_beta_list_bayes,
              "between_beta_list_lr" = between_beta_list_lr,
              "between_beta_list_bayes" = between_beta_list_bayes,
              "error_TrueB_list_lr" = error_TrueB_list_lr,
              "error_TrueB_list_bayes" = error_TrueB_list_bayes,
              "mse_list_lr" = mse_list_lr,
              "mse_list_bayes" = mse_list_bayes))
}



suppressWarnings({
  struc_diff_1_mixed = decor_struct_mmhc(n = 100, p = 100, cov_struc = c(toeplitz_struc, const_struc), reps = 1, vers = 0)
  save(struc_diff_1_mixed, file = "data_dag_uncor_mixed_1_bayes_ridge_lr.RData")
  struc_diff_2_mixed = decor_struct_mmhc(n = 100, p = 500, cov_struc = c(toeplitz_struc, const_struc), reps = 10, vers = 0)
  save(struc_diff_2_mixed, file = "data_dag_uncor_mixed_2_bayes_ridge_lr.RData")
  struc_diff_3_mixed = decor_struct_mmhc(n = 100, p = 1000, cov_struc = c(toeplitz_struc, const_struc), reps = 5, vers = 12345)
  save(struc_diff_3_mixed, file = "data_dag_uncor_mixed_3_1_bayes_ridge_lr.RData")
  #struc_diff_4_mixed = decor_struct_mmhc(n = 300, p = 100, cov_struc = c(toeplitz_struc, const_struc), reps = 10, vers = 0)
  #save(struc_diff_4_mixed, file = "data_dag_uncor_mixed_4.RData")
  #struc_diff_5_mixed = decor_struct_mmhc(n = 300, p = 500, cov_struc = c(toeplitz_struc, const_struc), reps = 10, vers = 0)
  #save(struc_diff_5_mixed, file = "data_dag_uncor_mixed_5.RData")
  #struc_diff_6_mixed = decor_struct_mmhc(n = 300, p = 1000, cov_struc = c(toeplitz_struc, const_struc), reps = 10, vers = 0)
  #save(struc_diff_6_mixed, file = "data_dag_uncor_mixed_6.RData")
  struc_diff_7_mixed = decor_struct_mmhc(n = 500, p = 100, cov_struc = c(toeplitz_struc, const_struc), reps = 10, vers = 0)
  save(struc_diff_7_mixed, file = "data_dag_uncor_mixed_7_bayes_ridge_lr.RData")
})

par(mfrow = c(1,2))

non_zero_betas_lr = c()
for(i in 1:21){
  non_zero_betas_lr = c(non_zero_betas_lr, )
}

par(mfrow = c(1,2))
plot(1:21, struc_diff_9_mixed$between_beta_list_lr[[1]], xlab = "Iterations", ylab = "RMSE Between Betas", type = "l", main = "Linear Reg")
plot(1:100, struc_diff_9_mixed$between_beta_list_bayes[[1]], ylab = "", xlab = "Iterations", type = "l", main = "Bayes")


