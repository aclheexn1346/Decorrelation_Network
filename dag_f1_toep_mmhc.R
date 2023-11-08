source("helperFunc.R")
source("libraries.R")

load("data_dag_uncor_toeplitz_1.RData")
load("data_dag_uncor_toeplitz_2.RData")
load("data_dag_uncor_toeplitz_3.RData")

load("data_dag_uncor_toeplitz_7.RData")


mmhc_pre_post_f1 = function(struc_diff_object){
  # Convert true DAG to CPDAG
  init_f1_list = list()
  current_f1_list = list()
  concensus_mat_list = list()
  true_beta_f1_list = list()
  for(j in 1:length(struc_diff_object$Beta)){
    true_beta = struc_diff_object$Beta[[j]]
    true_beta[which(true_beta != 0)] = 1
    colnames(true_beta) <- rownames(true_beta) <- 1:nrow(true_beta)
    g1 <- as(t(true_beta), "graphNEL")
    true_cpdag <- dag2cpdag(g1)
    ############################
    # Getting initial CPDAG estimate from data
    ###################################
    disc_data = struc_diff_object$Data[[j]]
    init_dag = init_data_dag(disc_data)
    init_g1 <- as(t(init_dag), "graphNEL")
    init_cpdag <- dag2cpdag(init_g1)
    init_f1 = get_F1_dag_parents(init_cpdag@edgeL, true_cpdag@edgeL)
    init_f1_list[[j]] = init_f1
    # Getting post CPDAG from current which uses last uncor_z
    #############################
    current_data= struc_diff_object$uncor_Z_list[[j]][[j]][[21]]
    current_dag = neighborhood_selection_mmhc_get_dag_decor(data = current_data)
    current_g1 <- as(t(current_dag), "graphNEL")
    current_cpdag <- dag2cpdag(current_g1)
    current_f1 = get_F1_dag_parents(current_cpdag@edgeL, true_cpdag@edgeL)
    current_f1_list[[j]] = current_f1
    # Getting post CPDAG from concensus estimate
    concensus_data = struc_diff_object$uncor_Z_list[[j]][[j]]
    concensus_mat = matrix(0, nrow = ncol(disc_data), ncol = ncol(disc_data))
    sampled_datasets = sample(2:21, size = 10)
    for(i in sampled_datasets){
      #print(i)
      parents = neighborhood_selection_mmhc_get_dag_decor(concensus_data[[i]])
      parents = dag.to.cpdag(parents)
      concensus_mat = concensus_mat + parents
    }
    concensus_mat_list[[j]] = concensus_mat
    #############################################
    # Post CPDAG from consensus that used true Dag structure
    true_beta_data = struc_diff_object$true_uncor_Z[[j]][[j]]
    true_beta_dag = neighborhood_selection_mmhc_get_dag_decor(data = true_beta_data)
    true_beta_g1 <- as(t(true_beta_dag), "graphNEL")
    true_beta_cpdag <- dag2cpdag(true_beta_g1)
    true_beta_f1 = get_F1_dag_parents(true_beta_cpdag@edgeL, true_cpdag@edgeL)
    true_beta_f1_list[[j]] = true_beta_f1
  }
  # concensus_mat[which(concensus_mat <= 1)] = 0
  # concensus_mat[which(concensus_mat >= 2)] = 1
  # concensus_g1 <- as(t(concensus_mat), "graphNEL")
  # concensus_cpdag <- dag2cpdag(concensus_g1)
  return(list("init" = init_f1_list, "current_f1" = current_f1_list, "concensus_f1" = concensus_mat_list, "true_beta_f1" = true_beta_f1_list))
}




f1_toep_1 = mmhc_pre_post_f1(struc_diff_1_toep)
save(f1_toep_1, file = "dag_f1_scores_toep_1_mmhc.RData")
f1_toep_2 = mmhc_pre_post_f1(struc_diff_2_toep)
save(f1_toep_2, file = "dag_f1_scores_toep_2_mmhc.RData")

f1_toep_7 = mmhc_pre_post_f1(struc_diff_7_toep)
save(f1_toep_7, file = "dag_f1_scores_toep_7_mmhc.RData")


############################



