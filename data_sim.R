#Using MMHC to find dag structure first and de-correlate data in order to
#test for causal structure improvement
source("helperFunc.R")
source("libraries.R")
decor_struct_mmhc = function(n, p, cov_struc, reps, vers){
  sig_2_list = list()
  omg.sq = rep(1, p)
  hold_uncor_datasets_list = list()
  hold_true_uncor_datasets_list = list()
  data_list = list()
  Beta_list = list()
  est_sigma_list = list()
  # reps is number of simulations
  for(i in 1:reps){
    # create block structure
    block_obj = block_diag_sep_var(n = n, min_block_size = 5, max_block_size = 10, cov_struc = cov_struc)
    sig_2 = block_obj[[1]]
    block_sizes = block_obj[[2]]
    sig_2_list[[i]] = sig_2
    hold_uncor_datasets = list()
    hold_true_uncor_datasets = list()
    vers = vers+1
    # Generate beta coefficient matrix between variables
    B = gen.B(p = p, seed = vers)
    # simulate data based on beta matrix
    sim_result_nplarge = sim_X_LUM(vers, p, n, omg.sq, sig = sig_2, b = B$b)
    
    ##########. In case the simulation results in all 0s or all 1s column ######
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
    #########################################################################
    true_epsilon = sim_result_nplarge$eps_mat
    # Initial DAG
    parents = neighborhood_selection_mmhc_get_dag(data = data)
    # Get initial beta
    init_beta = init_beta_0_mb(data = data, mb = parents)
    true_beta = B$b
    true_beta[which(true_beta != 0)] = 0.001
    init_epsilon = matrix(0, nrow = n, ncol = p)
    beta_est_ident_b0 = NULL
    while(is.null(beta_est_ident_b0)){
      try({
        # Perform Z imputation and de-correlation. Output is de-correlated Zs
        beta_est_ident_b0 <- beta_est_loop4(data = data, init_beta = init_beta, init_epsilon = init_epsilon, block_sizes = block_obj[[2]], loops = 21, trueB = B$b)
      })
      if(is.null(beta_est_ident_b0)){
        vers = vers+1
        sim_result_nplarge = sim_X_LUM_multi(vers, p, n, omg.sq, sig = sig_2, b = B$b)
        data = sim_result_nplarge$X
        true_epsilon = sim_result_nplarge$eps_mat
      }
    }
    uncor_data = beta_est_ident_b0$hold_data
    hold_uncor_datasets[[i]] = uncor_data
    hold_uncor_datasets_list[[i]] = hold_uncor_datasets
    #hold_true_uncor_datasets_list[[i]] = hold_true_uncor_datasets
    data_list[[i]] = data
    Beta_list[[i]] = B$b
    est_sigma_list[[i]] = beta_est_ident_b0$Sigma
  }
  return(list("Data" = data_list, 
              "Beta" = Beta_list,
              "est_Sigmas" = est_sigma_list,
              "true_sigma" = sig_2_list,
              "uncor_Z_list" = hold_uncor_datasets_list))
}

suppressWarnings({
  struc_diff_1_mixed = decor_struct_mmhc(n = 100, p = 100, cov_struc = c(toeplitz_struc, const_struc), reps = 10, vers = 0)
})