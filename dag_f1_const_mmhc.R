source("helperFunc.R")
source("libraries.R")

load("data_dag_uncor_toeplitz_1.RData")
load("data_dag_uncor_toeplitz_2.RData")
load("data_dag_uncor_toeplitz_3.RData")

load("data_dag_uncor_toeplitz_7.RData")

MMHC_pre_post_f1 = function(struc_diff_object){
  # Convert true DAG to CPDAG
  init_f1_list = list()
  init_f1_list = list()
  Consensus_mat_list = list()
  true_beta_f1_list = list()
  avg_mat_list = list()
  for(j in 1:length(struc_diff_object$Beta)){
    true_beta = struc_diff_object$Beta[[j]]
    true_beta[which(true_beta != 0)] = 1
    colnames(true_beta) <- rownames(true_beta) <- 1:nrow(true_beta)
    g1 <- as(t(true_beta), "graphNEL")
    true_cpdag <- dag2cpdag(g1)
    ############################
    # Getting initial CPDAG estimate from data
    ###################################
    # disc_data = struc_diff_object$Data[[j]]
    # init_dag = init_data_dag(disc_data)
    # init_g1 <- as(t(init_dag), "graphNEL")
    # init_cpdag <- dag2cpdag(init_g1)
    # init_f1 = get_F1_dag_parents(init_cpdag@edgeL, true_cpdag@edgeL)
    # init_f1_list[[j]] = init_f1
    # # Getting post CPDAG from init which uses last uncor_z
    # #############################
    # init_data= struc_diff_object$uncor_Z_list[[j]][[j]][[21]]
    # init_dag = neighborhood_selection_mmhc_get_dag_decor(data = init_data)
    # init_g1 <- as(t(init_dag), "graphNEL")
    # init_cpdag <- dag2cpdag(init_g1)
    # init_f1 = get_F1_dag_parents(init_cpdag@edgeL, true_cpdag@edgeL)
    # init_f1_list[[j]] = init_f1
    # # Getting post CPDAG from Consensus estimate
    # Consensus_data = struc_diff_object$uncor_Z_list[[j]][[j]]
    # Consensus_mat = matrix(0, nrow = ncol(disc_data), ncol = ncol(disc_data))
    # sampled_datasets = sample(2:21, size = 10)
    # for(i in sampled_datasets){
    #   #print(i)
    #   parents = neighborhood_selection_mmhc_get_dag_decor(Consensus_data[[i]])
    #   parents = dag.to.cpdag(parents)
    #   concensus_mat = concensus_mat + parents
    # }
    # concensus_mat_list[[j]] = concensus_mat
    # #############################################
    # # Post CPDAG from consensus that used true Dag structure
    # true_beta_data = struc_diff_object$true_uncor_Z[[j]][[j]]
    # true_beta_dag = neighborhood_selection_mmhc_get_dag_decor(data = true_beta_data)
    # true_beta_g1 <- as(t(true_beta_dag), "graphNEL")
    # true_beta_cpdag <- dag2cpdag(true_beta_g1)
    # true_beta_f1 = get_F1_dag_parents(true_beta_cpdag@edgeL, true_cpdag@edgeL)
    # true_beta_f1_list[[j]] = true_beta_f1
    #####################
    # Average out decorrelated data
    avg_mat = matrix(data = 0, nrow = nrow(struc_diff_object$Data[[j]]), ncol = ncol(struc_diff_object$Data[[j]]))
    for(k in 2:length(struc_diff_object$uncor_Z_list[[j]][[j]])){
      avg_mat = avg_mat + struc_diff_object$uncor_Z_list[[j]][[j]][[k]]
    }
    avg_mat = avg_mat/20
    avg_dag = neighborhood_selection_mmhc_get_dag_decor(data = avg_mat)
    true_beta_g1 <- as(t(avg_dag), "graphNEL")
    true_beta_cpdag <- dag2cpdag(true_beta_g1)
    true_beta_f1 = get_F1_dag_parents(true_beta_cpdag@edgeL, true_cpdag@edgeL)
    #true_beta_f1_list[[j]] = true_beta_f1
    avg_mat_list[[j]] = true_beta_f1
  }
  # concensus_mat[which(concensus_mat <= 1)] = 0
  # concensus_mat[which(concensus_mat >= 2)] = 1
  # concensus_g1 <- as(t(concensus_mat), "graphNEL")
  # concensus_cpdag <- dag2cpdag(concensus_g1)
  return(list("f1_avg_mat" = avg_mat_list))
}




f1_toep_1 = MMHC_pre_post_f1(struc_diff_1_toep)
save(f1_toep_1, file = "dag_f1_scores_toep_1_avg.RData")
f1_toep_2 = MMHC_pre_post_f1(struc_diff_2_toep)
save(f1_toep_2, file = "dag_f1_scores_toep_2_avg.RData")
f1_toep_3 = MMHC_pre_post_f1(struc_diff_3_toep)
save(f1_toep_3, file = "dag_f1_scores_toep_3_avg.RData")

f1_toep_7 = MMHC_pre_post_f1(struc_diff_7_toep)
save(f1_toep_7, file = "dag_f1_scores_toep_7_avg.RData")
############################

get_f1_consensus = function(struc_diff_object, dag_f1_object){
  true_beta = struc_diff_object$Beta[[1]]
  true_beta[which(true_beta != 0)] = 1
  colnames(true_beta) <- rownames(true_beta) <- 1:nrow(true_beta)
  g1 <- as(t(true_beta), "graphNEL")
  true_cpdag <- dag2cpdag(g1)
  consensus_vals = c()
  for(i in c(4)){
    concensus_mat = dag_f1_object$concensus_f1[[1]]
    concensus_mat[which(concensus_mat <= i)] = 0
    concensus_mat[which(concensus_mat >= i+1)] = 1
    concensus_g1 <- as(t(concensus_mat), "graphNEL")
    concensus_cpdag <- dag2cpdag(concensus_g1)
    concensus_f1 = get_F1_dag_parents(concensus_cpdag@edgeL, true_cpdag@edgeL)
    consensus_vals = c(consensus_vals, concensus_f1$f1_score)
  }
  return(consensus_vals)
}


source("helperFunc.R")
source("libraries.R")
avg = 0
for(i in 1:length(mat_list)){
  concensus_mat = mat_list[[i]]
  true_beta = B_mat[[i]]
  true_beta[which(true_beta != 0)] = 1
  colnames(true_beta) <- rownames(true_beta) <- 1:nrow(true_beta)
  g1 <- as(t(true_beta), "graphNEL")
  true_cpdag <- dag2cpdag(g1)
  concensus_mat[which(concensus_mat <= 4)] = 0
  concensus_mat[which(concensus_mat >= 5)] = 1
  concensus_g1 <- as(t(concensus_mat), "graphNEL")
  concensus_cpdag <- dag2cpdag(concensus_g1)
  concensus_f1 = get_F1_dag_parents(concensus_cpdag@edgeL, true_cpdag@edgeL)
  print(concensus_f1)
  avg = avg + concensus_f1$f1_score
  
}
avg/10


for(i in 1:22){
  print(length(which(beta_est_ident_lasso$hold_beta[[i]] != 0)))
}

avg = 0
mixed_9_concensus_f1 = c()
for(j in 1:10){
  print(j)
  concensus_mat = mat_list[[j]]
  true_beta = B_mat
  true_beta[which(true_beta != 0)] = 1
  colnames(true_beta) <- rownames(true_beta) <- 1:nrow(true_beta)
  g1 <- as(t(true_beta), "graphNEL")
  true_cpdag <- dag2cpdag(g1)
  concensus_mat[which(concensus_mat <= 4)] = 0
  concensus_mat[which(concensus_mat >= 5)] = 1
  concensus_g1 <- as(t(concensus_mat), "graphNEL")
  concensus_cpdag <- dag2cpdag(concensus_g1)
  concensus_f1 = get_F1_dag_parents(concensus_cpdag@edgeL, true_cpdag@edgeL)
  mixed_9_concensus_f1 = c(mixed_9_concensus_f1, concensus_f1$f1_score)
}

mat_list = f1_mixed_1$concensus_f1
mixed_andes_1_concensus_f1 = mixed_9_concensus_f1
save(mixed_andes_1_concensus_f1, file = "andes_1_consensus_pc.RData")
save(f1_mixed_munin1_1_avg, file = "munin1_1_avg_pc.RData")
save(f1_mixed_munin1_1_init, file = "munin1_1_init_pc.RData")

mean(f1_mixed_pigs_1_init)
mean(f1_mixed_pigs_1_avg)
mean(mixed_pigs_1_concensus_f1)

mean(f1_mixed_pigs_2_init)
mean(f1_mixed_pigs_2_avg)
mean(mixed_pigs_2_concensus_f1)

mat_list = list()
mat_list[[1]] = f1_mixed_9$concensus_f1[[1]]
mat_list[[2]] = f1_mixed_9$concensus_f1[[2]]
mat_list[[3]] = concensus_mat_list[[3]]
mat_list[[4]] = concensus_mat_list[[4]]
mat_list[[5]] = concensus_mat_list[[5]]
mat_list[[6]] = concensus_mat_list[[6]]
mat_list[[7]] = concensus_mat_list[[7]]
mat_list[[8]] = concensus_mat_list[[8]]
mat_list[[9]] = concensus_mat_list[[9]]
mat_list[[10]] = concensus_mat_list[[10]]


f1_mixed_munin1_1_init = c()
for(i in 1:10){
  f1_mixed_munin1_1_init = c(f1_mixed_munin1_1_init, f1_mixed_1$init[[i]]$f1_score)
}

f1_mixed_munin1_1_avg = c()
for(i in 1:10){
  f1_mixed_munin1_1_avg = c(f1_mixed_munin1_1_avg, f1_mixed_1$avg_f1[[i]]$f1_score)
}


# one iteration of algorithm then clustering algorithm

f1_mixed_9_1_true_beta = true_beta_f1_list[[1]]$f1_score
f1_mixed_9_2_true_beta = true_beta_f1_list[[2]]$f1_score
f1_mixed_9_3_true_beta = true_beta_f1_list[[3]]$f1_score
f1_mixed_9_4_true_beta = true_beta_f1_list[[4]]$f1_score
f1_mixed_9_5_true_beta = true_beta_f1_list[[5]]$f1_score
f1_mixed_9_6_true_beta = true_beta_f1_list[[6]]$f1_score
f1_mixed_9_7_true_beta = true_beta_f1_list[[7]]$f1_score
f1_mixed_9_8_true_beta = true_beta_f1_list[[8]]$f1_score
f1_mixed_9_9_true_beta = true_beta_f1_list[[9]]$f1_score
f1_mixed_9_10_true_beta = true_beta_f1_list[[10]]$f1_score


f1_mixed_9_true_beta = c(f1_mixed_9_1_true_beta,
                    f1_mixed_9_2_true_beta,
                    f1_mixed_9_3_true_beta,
                    f1_mixed_9_4_true_beta,
                    f1_mixed_9_5_true_beta,
                    f1_mixed_9_6_true_beta,
                    f1_mixed_9_7_true_beta,
                    f1_mixed_9_8_true_beta,
                    f1_mixed_9_9_true_beta,
                    f1_mixed_9_10_true_beta)
save(f1_mixed_9_true_beta, file = "mixed_9_true_beta_f1_pchc.RData")





save(f1_cons_8_pc_true_beta, file = "cons_8_true_beta_f1_pc.RData")
f1_mixed_1_init_mmhc = c()
for(i in 1:length(f1_mixed_1$init)){
  f1_mixed_1_init_mmhc[i] = f1_mixed_1$init[[i]]$f1_score 
}
f1_mixed_1_avg_mmhc = c()
for(i in 1:length(f1_mixed_1$init)){
  f1_mixed_1_avg_mmhc[i] = f1_mixed_1$avg_f1[[i]]$f1_score 
}


f1_mixed_3_concensus_mmhc = c()


concensus_mat = f1_scores$concensus_f1

true_beta = struc_diff_3_mixed$Beta[[10]]
true_beta[which(true_beta != 0)] = 1
colnames(true_beta) <- rownames(true_beta) <- 1:nrow(true_beta)
g1 <- as(t(true_beta), "graphNEL")
true_cpdag <- dag2cpdag(g1)

concensus_mat[which(concensus_mat <= 4)] = 0
concensus_mat[which(concensus_mat >= 5)] = 1
concensus_g1 <- as(t(concensus_mat), "graphNEL")
concensus_cpdag <- dag2cpdag(concensus_g1)
concensus_f1 = get_F1_dag_parents(concensus_cpdag@edgeL, true_cpdag@edgeL)
f1_mixed_3_concensus_mmhc = c(f1_mixed_3_concensus_mmhc, concensus_f1$f1_score)

  #avg = avg + concensus_f1$f1_score
f1_mixed_3_init_mmhc = c()
f1_mixed_3_avg_mmhc = c()
f1_mixed_3_avg_mmhc[10] = f1_scores$avg_f1$f1_score # need to this one when it's done running
f1_mixed_3_init_mmhc[10] = f1_scores$init$f1_score

f1_final_scores_3 = list(f1_mixed_3_init_mmhc,
                         f1_mixed_3_avg_mmhc,
                         f1_mixed_3_concensus_mmhc)
save(f1_final_scores_3, file = "final_f1_scores_3.RData")


struc_diff_9_mixed_1 = struc_diff_9_mixed
struc_diff_9_mixed_1$uncor_Z_list[[2]] = struc_diff_9_mixed$uncor_Z_list[[1]]
struc_diff_9_mixed_1$uncor_Z_list_lr[[2]] = struc_diff_9_mixed$uncor_Z_list_lr[[1]]

struc_diff_9_mixed_1$Data[[2]] = struc_diff_9_mixed$Data[[1]]
struc_diff_9_mixed_1$Beta[[2]] = struc_diff_9_mixed$Beta[[1]]


struc_diff_9_mixed_1$uncor_Z_list[[3]] = struc_diff_9_mixed$uncor_Z_list[[1]]
struc_diff_9_mixed_1$uncor_Z_list_lr[[3]] = struc_diff_9_mixed$uncor_Z_list_lr[[1]]

struc_diff_9_mixed_1$Data[[3]] = struc_diff_9_mixed$Data[[1]]
struc_diff_9_mixed_1$Beta[[3]] = struc_diff_9_mixed$Beta[[1]]

struc_diff_9_mixed_1$uncor_Z_list[[4]] = struc_diff_9_mixed$uncor_Z_list[[1]]
struc_diff_9_mixed_1$uncor_Z_list_lr[[4]] = struc_diff_9_mixed$uncor_Z_list_lr[[1]]

struc_diff_9_mixed_1$Data[[4]] = struc_diff_9_mixed$Data[[1]]
struc_diff_9_mixed_1$Beta[[4]] = struc_diff_9_mixed$Beta[[1]]

struc_diff_9_mixed_1$uncor_Z_list[[5]] = struc_diff_9_mixed$uncor_Z_list[[1]]
struc_diff_9_mixed_1$uncor_Z_list_lr[[5]] = struc_diff_9_mixed$uncor_Z_list_lr[[1]]

struc_diff_9_mixed_1$Data[[5]] = struc_diff_9_mixed$Data[[1]]
struc_diff_9_mixed_1$Beta[[5]] = struc_diff_9_mixed$Beta[[1]]

struc_diff_9_mixed_1$uncor_Z_list[[6]] = struc_diff_9_mixed$uncor_Z_list[[1]]
struc_diff_9_mixed_1$uncor_Z_list_lr[[6]] = struc_diff_9_mixed$uncor_Z_list_lr[[1]]

struc_diff_9_mixed_1$Data[[6]] = struc_diff_9_mixed$Data[[1]]
struc_diff_9_mixed_1$Beta[[6]] = struc_diff_9_mixed$Beta[[1]]

struc_diff_9_mixed_1$uncor_Z_list[[7]] = struc_diff_9_mixed$uncor_Z_list[[1]]
struc_diff_9_mixed_1$uncor_Z_list_lr[[7]] = struc_diff_9_mixed$uncor_Z_list_lr[[1]]

struc_diff_9_mixed_1$Data[[7]] = struc_diff_9_mixed$Data[[1]]
struc_diff_9_mixed_1$Beta[[7]] = struc_diff_9_mixed$Beta[[1]]

struc_diff_9_mixed_1$uncor_Z_list[[8]] = struc_diff_9_mixed$uncor_Z_list[[1]]
struc_diff_9_mixed_1$uncor_Z_list_lr[[8]] = struc_diff_9_mixed$uncor_Z_list_lr[[1]]

struc_diff_9_mixed_1$Data[[8]] = struc_diff_9_mixed$Data[[1]]
struc_diff_9_mixed_1$Beta[[8]] = struc_diff_9_mixed$Beta[[1]]

struc_diff_9_mixed_1$uncor_Z_list[[9]] = struc_diff_9_mixed$uncor_Z_list[[1]]
struc_diff_9_mixed_1$uncor_Z_list_lr[[9]] = struc_diff_9_mixed$uncor_Z_list_lr[[1]]

struc_diff_9_mixed_1$Data[[9]] = struc_diff_9_mixed$Data[[1]]
struc_diff_9_mixed_1$Beta[[9]] = struc_diff_9_mixed$Beta[[1]]

struc_diff_9_mixed_1$uncor_Z_list[[10]] = struc_diff_9_mixed$uncor_Z_list[[1]]
struc_diff_9_mixed_1$uncor_Z_list_lr[[10]] = struc_diff_9_mixed$uncor_Z_list_lr[[1]]

struc_diff_9_mixed_1$Data[[10]] = struc_diff_9_mixed$Data[[1]]
struc_diff_9_mixed_1$Beta[[10]] = struc_diff_9_mixed$Beta[[1]]

struc_diff_9_mixed = struc_diff_9_mixed_1

save(struc_diff_9_mixed, file = "struc_diff_9_mixed.RData")

f1_mixed_3_init = c()
for(i in 1:length(f1_mixed_7$init)){
  f1_mixed_9_init[i] = f1_mixed_9$init[[i]]$f1_score
}

struc_diff_3_mixed = f1_mixed_3

struc_diff_3_mixed[[1]][[1]] = f1_mixed_3[[1]][[1]]
struc_diff_3_mixed[[2]][[1]] = f1_mixed_3[[2]][[1]]
struc_diff_3_mixed[[3]][[1]] = f1_mixed_3[[3]][[1]]

f1_9_mixed_init = c()
for(i in 1:10){
  f1_9_mixed_init = c(f1_9_mixed_init, struc_diff_9_mixed[[1]][[i]]$f1_score)
}

f1_9_mixed_lr = c()
for(i in 1:10){
  f1_9_mixed_lr = c(f1_9_mixed_lr, struc_diff_9_mixed[[2]][[i]]$f1_score)
}

f1_9_mixed_bayes = c()
for(i in 1:10){
  f1_9_mixed_bayes = c(f1_9_mixed_bayes, struc_diff_9_mixed[[3]][[i]]$f1_score)
}

f1_9_mixed_bayes_avg = c()
for(i in 6:10){
  f1_9_mixed_bayes_avg = c(f1_9_mixed_bayes_avg, f1_mixed_9$bayes[[i]]$f1_score)
}


f1_data_1_bayes_ridge = data.frame("f1_score" = c(f1_1_mixed_init,
                                           f1_1_mixed_bayes,
                                           f1_1_mixed_bayes_avg,
                                      f1_1_mixed_lr),
                       "Data" = c(rep("Baseline", 10),
                                  rep("Bayes", 10),
                                  rep("Avg Bayes", 10),
                                  rep("Linear", 10)))

f1_data_1_bayes_ridge$Data = factor(f1_data_1_bayes_ridge$Data, levels = c("Baseline","Bayes","Avg Bayes", "Linear"))

f1_data_3_bayes_ridge = data.frame("f1_score" = c(f1_3_mixed_init,
                                      f1_3_mixed_bayes,
                                      f1_3_mixed_bayes_avg,
                                      f1_3_mixed_lr),
                       "Data" = c(rep("Baseline", 10),
                                  rep("Bayes", 10),
                                  rep("Avg Bayes", 10),
                                  rep("Linear", 10)))
f1_data_3_bayes_ridge$Data = factor(f1_data_3_bayes_ridge$Data, levels = c("Baseline", "Bayes","Avg Bayes", "Linear"))

f1_data_9_bayes_ridge = data.frame("f1_score" = c(f1_9_mixed_init,
                                      f1_9_mixed_bayes,
                                      f1_9_mixed_bayes_avg,
                                      f1_9_mixed_lr),
                       "Data" = c(rep("Baseline", 10),
                                  rep("Bayes", 10),
                                  rep("Avg Bayes", 10),
                                  rep("Linear", 10)))
f1_data_9_bayes_ridge$Data = factor(f1_data_9_bayes_ridge$Data, levels = c("Baseline", "Bayes","Avg Bayes", "Linear"))

f1_data_7_bayes_ridge = data.frame("f1_score" = c(f1_7_mixed_init,
                                      f1_7_mixed_bayes,
                                      f1_7_mixed_bayes_avg,
                                      f1_7_mixed_lr),
                       "Data" = c(rep("Baseline", 10),
                                  rep("Bayes", 10),
                                  rep("Avg Bayes", 10),
                                  rep("Linear", 10)))
f1_data_7_bayes_ridge$Data = factor(f1_data_7_bayes_ridge$Data, levels = c("Baseline", "Bayes","Avg Bayes", "Linear"))



p1 = ggplot(f1_data_1, aes(x = Data, y = f1_score, fill = Data)) +
  geom_boxplot() +
  xlab('Data') + ylab('F1-Score') + theme_bw()+ ggtitle("Toeplitz MMHC, n = 100, p = 100")

p2 = ggplot(f1_data_2, aes(x = Data, y = f1_score, fill = Data)) +
  geom_boxplot() +
  xlab('Data') + ylab('F1-Score') + theme_bw()+
  ggtitle("toeplitz, n = 100, p = 500")

p3 = ggplot(f1_data_3, aes(x = Data, y = f1_score, fill = Data)) +
  geom_boxplot() +
  xlab('Data') + ylab('F1-Score') + theme_bw()+
  theme(text = element_text(size=16),
        axis.text.x = element_text(angle=90, hjust=.5, vjust = 0.5),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))+
  scale_y_continuous(breaks = seq(0, 0.4, by=0.1), limits=c(0,.4))+
  ggtitle("toeplitz, n = 100, p = 1000")

p4 = ggplot(f1_data_7, aes(x = Data, y = f1_score, fill = Data)) +
  geom_boxplot() +
  xlab('Data') + ylab('F1-Score') + theme_bw()+
  theme(text = element_text(size=16),
        axis.text.x = element_text(angle=90, hjust=.5, vjust = 0.5),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))+
  scale_y_continuous(breaks = seq(0.4, 0.9, by=0.1), limits=c(0.45,.9))+
  ggtitle("toeplitz, n = 500, p = 100")

# p <- ggplot(f1_data_7, aes(Data, f1_score,fill=Data))
# p + geom_boxplot() + labs(title = "CMP")

p5 = ggplot(f1_data_8, aes(x = Data, y = f1_score, fill = Data)) +
  geom_boxplot() +
  xlab('Data') + ylab('F1-Score') + theme_bw()+
  theme(text = element_text(size=16),
        axis.text.x = element_text(angle=90, hjust=.5, vjust = 0.5),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))+
  scale_y_continuous(breaks = seq(0.4, 0.9, by=0.1), limits=c(0.45,.9))+
  ggtitle("toeplitz, n = 500, p = 500")

p6 = ggplot(f1_data_9, aes(x = Data, y = f1_score, fill = Data)) +
  geom_boxplot() +
  xlab('Data') + ylab('F1-Score') + theme_bw()+
  theme(text = element_text(size=16),
        axis.text.x = element_text(angle=90, hjust=.5, vjust = 0.5),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))+
  scale_y_continuous(breaks = seq(0.4, 0.9, by=0.1), limits=c(0.45,.75))+
  ggtitle("toeplitz, n = 500, p = 1000")

f1_data_1_bayes$type = rep(c("Bayes"), 30)
#f1_data_2_mmhc$type = rep(c("mmhc"), 25)
f1_data_3_bayes$type = rep(c("Bayes"), 30)
f1_data_7_bayes$type = rep(c("Bayes"), 30)
#f1_data_8_mmhc$type = rep(c("mmhc"), 30)
f1_data_9_bayes$type = rep(c("Bayes"), 30)


f1_data_1_pc$type = rep(c("PC"), 30)



f1_data_1_tot = rbind(f1_data_1_bayes, f1_data_1_linreg)
#f1_data_2_tot = rbind(f1_data_2_mmhc, f1_data_2_pc)
f1_data_3_tot = rbind(f1_data_3_bayes, f1_data_3_linreg)
f1_data_7_tot = rbind(f1_data_7_bayes, f1_data_7_linreg)
#f1_data_8_tot = rbind(f1_data_8_mmhc, f1_data_8_pc)
f1_data_9_tot = rbind(f1_data_9_bayes, f1_data_9_linreg)


f1_data$type = factor(f1_data$type, levels = c("(100,100)","(100,500)","(100,1000)","(500,100)","(500,500)","(500,1000)"))
f1_data = f1_data[-which(f1_data$f1_score <= 0.05),]
library(ggplot2)
p1 <- ggplot(f1_data_1_bayes_ridge, aes(Data, f1_score,fill=Data)) + geom_boxplot() + labs(title = "n = 100, p = 100") +  ylab("F1-Score") + theme_bw() + theme(legend.position = "none",axis.title.x=element_blank(), axis.text = element_text(size = 18), plot.title = element_text(size=18, face = "bold"),
                                                                                                                                                        
                                                                                                                                                        axis.title.y = element_text(size=18, face = "bold"), axis.text.x = element_blank())+ theme(legend.title=element_blank()) + theme(legend.text=element_text(size=18))
p3 <- ggplot(f1_data_3_bayes_ridge, aes(Data, f1_score,fill=Data)) + geom_boxplot() + labs(title = "n = 100, p = 1000") +  ylab("F1-Score") + theme_bw()+ theme(legend.position = "none",axis.title.x=element_blank(), axis.text = element_text(size = 18), plot.title = element_text(size=18, face = "bold"),
                                                                                                                                                      
                                                                                                                                                        axis.title.y = element_text(size=18, face = "bold"), axis.text.x = element_blank())+ theme(legend.title=element_blank())+ theme(legend.text=element_text(size=18))
p7 <- ggplot(f1_data_7_bayes_ridge, aes(Data, f1_score,fill=Data)) + geom_boxplot() + labs(title = "n = 500, p = 100") +  ylab("F1-Score") + theme_bw()+ theme(legend.position = "none",axis.title.x=element_blank(), axis.text = element_text(size = 18), plot.title = element_text(size=18, face = "bold"),
                                                                                                                                                       axis.title.y = element_text(size=18, face = "bold"), axis.text.x = element_blank())+ theme(legend.title=element_blank())+ theme(legend.text=element_text(size=18))
p9 <- ggplot(f1_data_9_bayes_ridge, aes(Data, f1_score,fill=Data)) + geom_boxplot() + labs(title = "n = 500, p = 1000") +  ylab("F1-Score") + theme_bw() + theme(axis.title.x=element_blank(),axis.text = element_text(size = 18), plot.title = element_text(size=18, face = "bold"),
                                                                                                                                                         
                                                                                                                                                         axis.title.y = element_text(size=18, face = "bold"), axis.text.x = element_blank()) + theme(legend.title=element_blank())+ theme(legend.text=element_text(size=18))



library(patchwork)
combined = p1 + p3 + p7 + p9 + plot_layout(ncol = 4) & theme(legend.position = "bottom")
combined + plot_layout(guides = "collect")


library(bnlearn)

hailfinder = bnlearn::hailfinder
modelstring = paste0("[N07muVerMo][SubjVertMo][QGVertMotion][SatContMoist][RaoContMoist]",
                     "[VISCloudCov][IRCloudCover][AMInstabMt][WndHodograph][MorningBound][LoLevMoistAd][Date]",
                     "[MorningCIN][LIfr12ZDENSd][AMDewptCalPl][LatestCIN][LLIW]",
                     "[CombVerMo|N07muVerMo:SubjVertMo:QGVertMotion][CombMoisture|SatContMoist:RaoContMoist]",
                     "[CombClouds|VISCloudCov:IRCloudCover][Scenario|Date][CurPropConv|LatestCIN:LLIW]",
                     "[AreaMesoALS|CombVerMo][ScenRelAMCIN|Scenario][ScenRelAMIns|Scenario][ScenRel34|Scenario]",
                     "[ScnRelPlFcst|Scenario][Dewpoints|Scenario][LowLLapse|Scenario][MeanRH|Scenario]",
                     "[MidLLapse|Scenario][MvmtFeatures|Scenario][RHRatio|Scenario][SfcWndShfDis|Scenario]",
                     "[SynForcng|Scenario][TempDis|Scenario][WindAloft|Scenario][WindFieldMt|Scenario]",
                     "[WindFieldPln|Scenario][AreaMoDryAir|AreaMesoALS:CombMoisture]",
                     "[AMCINInScen|ScenRelAMCIN:MorningCIN][AMInsWliScen|ScenRelAMIns:LIfr12ZDENSd:AMDewptCalPl]",
                     "[CldShadeOth|AreaMesoALS:AreaMoDryAir:CombClouds][InsInMt|CldShadeOth:AMInstabMt]",
                     "[OutflowFrMt|InsInMt:WndHodograph][CldShadeConv|InsInMt:WndHodograph][MountainFcst|InsInMt]",
                     "[Boundaries|WndHodograph:OutflowFrMt:MorningBound][N34StarFcst|ScenRel34:PlainsFcst]",
                     "[CompPlFcst|AreaMesoALS:CldShadeOth:Boundaries:CldShadeConv][CapChange|CompPlFcst]",
                     "[InsChange|CompPlFcst:LoLevMoistAd][CapInScen|CapChange:AMCINInScen]",
                     "[InsSclInScen|InsChange:AMInsWliScen][R5Fcst|MountainFcst:N34StarFcst]",
                     "[PlainsFcst|CapInScen:InsSclInScen:CurPropConv:ScnRelPlFcst]")
dag = model2network(modelstring)


var_names = names(bn)
B_mat = matrix(0, nrow = length(var_names), ncol = length(var_names))
colnames(B_mat) = rownames(B_mat) = var_names
for(i in 1:length(bn)){
  B_mat[bn[[i]]$parents, names(bn)[[i]]] = 1
}

save(B_mat, file = "bnlearn_pigs_B.RData")


diabetes_1_pc = data.frame("f1_score" = c(f1_mixed_diabetes_1_init,
                                          f1_mixed_diabetes_1_avg,
                                          mixed_diabetes_1_concensus_f1),
                          "Data" = c(rep("Baseline", 10),
                                     rep("Average", 10),
                                     rep("Consensus", 10)))
diabetes_1_pc$Data = factor(diabetes_1_pc$Data, levels = c("Baseline", "Average", "Consensus"))


diabetes_2_pc = data.frame("f1_score" = c(f1_mixed_diabetes_2_init,
                                          f1_mixed_diabetes_2_avg,
                                          mixed_diabetes_2_concensus_f1),
                          "Data" = c(rep("Baseline", 10),
                                     rep("Average", 10),
                                     rep("Consensus", 10)))
diabetes_2_pc$Data = factor(diabetes_2_pc$Data, levels = c("Baseline", "Average", "Consensus"))

diabetes_1_mmhc = data.frame("f1_score" = c(f1_mixed_diabetes_1_init,
                                         f1_mixed_diabetes_1_avg,
                                         mixed_diabetes_1_concensus_f1),
                          "Data" = c(rep("Baseline", 10),
                                     rep("Average", 10),
                                     rep("Consensus", 10)))
diabetes_1_mmhc$Data = factor(diabetes_1_mmhc$Data, levels = c("Baseline", "Average", "Consensus"))

diabetes_2_mmhc = data.frame("f1_score" = c(f1_mixed_diabetes_2_init,
                                            f1_mixed_diabetes_2_avg,
                                            mixed_disease_2_concensus_f1),
                          "Data" = c(rep("Baseline", 10),
                                     rep("Average", 10),
                                     rep("Consensus", 10)))
diabetes_2_mmhc$Data = factor(diabetes_2_mmhc$Data, levels = c("Baseline", "Average", "Consensus"))

####################################################################################################

munin1_1_pc = data.frame("f1_score" = c(f1_mixed_munin1_1_init,
                                          f1_mixed_munin1_1_avg,
                                          mixed_munin1_1_concensus_f1),
                           "Data" = c(rep("Baseline", 10),
                                      rep("Average", 10),
                                      rep("Consensus", 10)))
munin1_1_pc$Data = factor(munin1_1_pc$Data, levels = c("Baseline", "Average", "Consensus"))


munin1_2_pc = data.frame("f1_score" = c(f1_mixed_munin1_2_init,
                                          f1_mixed_munin1_2_avg,
                                          mixed_munin1_2_concensus_f1),
                           "Data" = c(rep("Baseline", 10),
                                      rep("Average", 10),
                                      rep("Consensus", 10)))
munin1_2_pc$Data = factor(munin1_2_pc$Data, levels = c("Baseline", "Average", "Consensus"))

munin1_1_mmhc = data.frame("f1_score" = c(f1_mixed_munin1_1_init,
                                            f1_mixed_munin1_1_avg,
                                            mixed_munin1_1_concensus_f1),
                             "Data" = c(rep("Baseline", 10),
                                        rep("Average", 10),
                                        rep("Consensus", 10)))
munin1_1_mmhc$Data = factor(munin1_1_mmhc$Data, levels = c("Baseline", "Average", "Consensus"))

munin1_2_mmhc = data.frame("f1_score" = c(f1_mixed_munin1_2_init,
                                            f1_mixed_munin1_2_avg,
                                            mixed_munin1_2_concensus_f1),
                             "Data" = c(rep("Baseline", 10),
                                        rep("Average", 10),
                                        rep("Consensus", 10)))
munin1_2_mmhc$Data = factor(munin1_2_mmhc$Data, levels = c("Baseline", "Average", "Consensus"))

####################################################################################################

hailfinder_1_pc = data.frame("f1_score" = c(f1_mixed_hailfinder_1_init,
                                        f1_mixed_hailfinder_1_avg,
                                        mixed_hailfinder_1_concensus_f1),
                         "Data" = c(rep("Baseline", 10),
                                    rep("Average", 10),
                                    rep("Consensus", 10)))
hailfinder_1_pc$Data = factor(hailfinder_1_pc$Data, levels = c("Baseline", "Average", "Consensus"))


hailfinder_2_pc = data.frame("f1_score" = c(f1_mixed_hailfinder_2_init,
                                        f1_mixed_hailfinder_2_avg,
                                        mixed_hailfinder_2_concensus_f1),
                         "Data" = c(rep("Baseline", 10),
                                    rep("Average", 10),
                                    rep("Consensus", 10)))
hailfinder_2_pc$Data = factor(hailfinder_2_pc$Data, levels = c("Baseline", "Average", "Consensus"))

hailfinder_1_mmhc = data.frame("f1_score" = c(f1_mixed_hailfinder_1_init,
                                          f1_mixed_hailfinder_1_avg,
                                          mixed_hailfinder_1_concensus_f1),
                           "Data" = c(rep("Baseline", 10),
                                      rep("Average", 10),
                                      rep("Consensus", 10)))
hailfinder_1_mmhc$Data = factor(hailfinder_1_mmhc$Data, levels = c("Baseline", "Average", "Consensus"))

hailfinder_2_mmhc = data.frame("f1_score" = c(f1_mixed_hailfinder_2_init,
                                          f1_mixed_hailfinder_2_avg,
                                          mixed_hailfinder_2_concensus_f1),
                           "Data" = c(rep("Baseline", 10),
                                      rep("Average", 10),
                                      rep("Consensus", 10)))
hailfinder_2_mmhc$Data = factor(hailfinder_2_mmhc$Data, levels = c("Baseline", "Average", "Consensus"))

####################################################################################################

hepar2_1_pc = data.frame("f1_score" = c(f1_mixed_hepar2_1_init,
                                        f1_mixed_hepar2_1_avg,
                                        mixed_hepar2_1_concensus_f1),
                         "Data" = c(rep("Baseline", 10),
                                    rep("Average", 10),
                                    rep("Consensus", 10)))
hepar2_1_pc$Data = factor(hepar2_1_pc$Data, levels = c("Baseline", "Average", "Consensus"))


hepar2_2_pc = data.frame("f1_score" = c(f1_mixed_hepar2_2_init,
                                        f1_mixed_hepar2_2_avg,
                                        mixed_hepar2_2_concensus_f1),
                         "Data" = c(rep("Baseline", 10),
                                    rep("Average", 10),
                                    rep("Consensus", 10)))
hepar2_2_pc$Data = factor(hepar2_2_pc$Data, levels = c("Baseline", "Average", "Consensus"))

hepar2_1_mmhc = data.frame("f1_score" = c(f1_mixed_hepar2_1_init,
                                          f1_mixed_hepar2_1_avg,
                                          mixed_hepar2_1_concensus_f1[1:10]),
                           "Data" = c(rep("Baseline", 10),
                                      rep("Average", 10),
                                      rep("Consensus", 10)))
hepar2_1_mmhc$Data = factor(hepar2_1_mmhc$Data, levels = c("Baseline", "Average", "Consensus"))

hepar2_2_mmhc = data.frame("f1_score" = c(f1_mixed_hepar2_2_init,
                                          f1_mixed_hepar2_2_avg,
                                          mixed_hepar2_2_concensus_f1),
                           "Data" = c(rep("Baseline", 10),
                                      rep("Average", 10),
                                      rep("Consensus", 10)))
hepar2_2_mmhc$Data = factor(hepar2_2_mmhc$Data, levels = c("Baseline", "Average", "Consensus"))

####################################################################################################

link_1_pc = data.frame("f1_score" = c(f1_mixed_link_1_init,
                                        f1_mixed_link_1_avg,
                                        mixed_link_1_concensus_f1),
                         "Data" = c(rep("Baseline", 10),
                                    rep("Average", 10),
                                    rep("Consensus", 10)))
link_1_pc$Data = factor(link_1_pc$Data, levels = c("Baseline", "Average", "Consensus"))


link_2_pc = data.frame("f1_score" = c(f1_mixed_link_2_init,
                                        f1_mixed_link_2_avg,
                                        mixed_link_2_concensus_f1),
                         "Data" = c(rep("Baseline", 11),
                                    rep("Average", 11),
                                    rep("Consensus", 11)))
link_2_pc$Data = factor(link_2_pc$Data, levels = c("Baseline", "Average", "Consensus"))

link_1_mmhc = data.frame("f1_score" = c(f1_mixed_link_1_init,
                                          f1_mixed_link_1_avg,
                                          mixed_link_1_concensus_f1),
                           "Data" = c(rep("Baseline", 10),
                                      rep("Average", 10),
                                      rep("Consensus", 10)))
link_1_mmhc$Data = factor(link_1_mmhc$Data, levels = c("Baseline", "Average", "Consensus"))

link_2_mmhc = data.frame("f1_score" = c(f1_mixed_link_2_init,
                                          f1_mixed_link_2_avg,
                                          mixed_link_2_concensus_f1),
                           "Data" = c(rep("Baseline", 11),
                                      rep("Average", 11),
                                      rep("Consensus", 11)))
link_2_mmhc$Data = factor(link_2_mmhc$Data, levels = c("Baseline", "Average", "Consensus"))

####################################################################################################

pigs_1_pc = data.frame("f1_score" = c(f1_mixed_pigs_1_init,
                                        f1_mixed_pigs_1_avg,
                                        mixed_pigs_1_concensus_f1),
                         "Data" = c(rep("Baseline", 10),
                                    rep("Average", 10),
                                    rep("Consensus", 10)))
pigs_1_pc$Data = factor(pigs_1_pc$Data, levels = c("Baseline", "Average", "Consensus"))


pigs_2_pc = data.frame("f1_score" = c(f1_mixed_pigs_2_init,
                                        f1_mixed_pigs_2_avg,
                                        mixed_pigs_2_concensus_f1),
                         "Data" = c(rep("Baseline", 10),
                                    rep("Average", 10),
                                    rep("Consensus", 10)))
pigs_2_pc$Data = factor(pigs_2_pc$Data, levels = c("Baseline", "Average", "Consensus"))

pigs_1_mmhc = data.frame("f1_score" = c(f1_mixed_pigs_1_init,
                                          f1_mixed_pigs_1_avg,
                                          mixed_pigs_1_concensus_f1),
                           "Data" = c(rep("Baseline", 10),
                                      rep("Average", 10),
                                      rep("Consensus", 10)))
pigs_1_mmhc$Data = factor(pigs_1_mmhc$Data, levels = c("Baseline", "Average", "Consensus"))

pigs_2_mmhc = data.frame("f1_score" = c(f1_mixed_pigs_2_init,
                                          f1_mixed_pigs_2_avg,
                                          mixed_pigs_2_concensus_f1),
                           "Data" = c(rep("Baseline", 10),
                                      rep("Average", 10),
                                      rep("Consensus", 10)))
pigs_2_mmhc$Data = factor(pigs_2_mmhc$Data, levels = c("Baseline", "Average", "Consensus"))

####################################################################################################

win95pts_1_pc = data.frame("f1_score" = c(f1_mixed_win95pts_1_init,
                                        f1_mixed_win95pts_1_avg,
                                        mixed_win95pts_1_concensus_f1),
                         "Data" = c(rep("Baseline", 10),
                                    rep("Average", 10),
                                    rep("Consensus", 10)))
win95pts_1_pc$Data = factor(win95pts_1_pc$Data, levels = c("Baseline", "Average", "Consensus"))


win95pts_2_pc = data.frame("f1_score" = c(f1_mixed_win95pts_2_init,
                                        f1_mixed_win95pts_2_avg,
                                        mixed_win95pts_2_concensus_f1),
                         "Data" = c(rep("Baseline", 10),
                                    rep("Average", 10),
                                    rep("Consensus", 10)))
win95pts_2_pc$Data = factor(win95pts_2_pc$Data, levels = c("Baseline", "Average", "Consensus"))

win95pts_1_mmhc = data.frame("f1_score" = c(f1_mixed_win95pts_1_init,
                                          f1_mixed_win95pts_1_avg,
                                          mixed_win95pts_1_concensus_f1),
                           "Data" = c(rep("Baseline", 10),
                                      rep("Average", 10),
                                      rep("Consensus", 10)))
win95pts_1_mmhc$Data = factor(win95pts_1_mmhc$Data, levels = c("Baseline", "Average", "Consensus"))

win95pts_2_mmhc = data.frame("f1_score" = c(f1_mixed_win95pts_2_init,
                                          f1_mixed_win95pts_2_avg,
                                          mixed_win95pts_2_concensus_f1),
                           "Data" = c(rep("Baseline", 10),
                                      rep("Average", 10),
                                      rep("Consensus", 10)))
win95pts_2_mmhc$Data = factor(win95pts_2_mmhc$Data, levels = c("Baseline", "Average", "Consensus"))

####################################################################################################

andes_1_pc = data.frame("f1_score" = c(f1_mixed_andes_1_init,
                                          f1_mixed_andes_1_avg,
                                          mixed_andes_1_concensus_f1),
                           "Data" = c(rep("Baseline", 10),
                                      rep("Average", 10),
                                      rep("Consensus", 10)))
andes_1_pc$Data = factor(andes_1_pc$Data, levels = c("Baseline", "Average", "Consensus"))


andes_2_pc = data.frame("f1_score" = c(f1_mixed_andes_2_init,
                                          f1_mixed_andes_2_avg,
                                          mixed_andes_2_Consensus_f1),
                           "Data" = c(rep("Baseline", 10),
                                      rep("Average", 10),
                                      rep("Consensus", 10)))
andes_2_pc$Data = factor(andes_2_pc$Data, levels = c("Baseline", "Average", "Consensus"))

andes_1_mmhc = data.frame("f1_score" = c(f1_mixed_andes_1_init,
                                            f1_mixed_andes_1_avg,
                                            mixed_andes_1_concensus_f1),
                             "Data" = c(rep("Baseline", 10),
                                        rep("Average", 10),
                                        rep("Consensus", 10)))
andes_1_mmhc$Data = factor(andes_1_mmhc$Data, levels = c("Baseline", "Average", "Consensus"))

andes_2_mmhc = data.frame("f1_score" = c(f1_mixed_andes_2_init,
                                            f1_mixed_andes_2_avg,
                                            mixed_andes_2_Consensus_f1),
                             "Data" = c(rep("Baseline", 10),
                                        rep("Average", 10),
                                        rep("Consensus", 10)))
andes_2_mmhc$Data = factor(andes_2_mmhc$Data, levels = c("Baseline", "Average", "Consensus"))



diabetes_1_pc$type = rep(c("PC"), 30)
diabetes_1_mmhc$type = rep(c("MMHC"),30)
diabetes_1_tot = rbind(diabetes_1_mmhc, diabetes_1_pc)

diabetes_2_pc$type = rep(c("PC"), 30)
diabetes_2_mmhc$type = rep(c("MMHC"),30)
diabetes_2_tot = rbind(diabetes_2_mmhc, diabetes_2_pc)

###########################################################
munin1_1_pc$type = rep(c("PC"), 30)
munin1_1_mmhc$type = rep(c("MMHC"),30)
munin1_1_tot = rbind(munin1_1_mmhc, munin1_1_pc)

munin1_2_pc$type = rep(c("PC"), 30)
munin1_2_mmhc$type = rep(c("MMHC"),30)
munin1_2_tot = rbind(munin1_2_mmhc, munin1_2_pc)

###########################################################
hailfinder_1_pc$type = rep(c("PC"), 30)
hailfinder_1_mmhc$type = rep(c("MMHC"),30)
hailfinder_1_tot = rbind(hailfinder_1_mmhc, hailfinder_1_pc)

hailfinder_2_pc$type = rep(c("PC"), 30)
hailfinder_2_mmhc$type = rep(c("MMHC"),30)
hailfinder_2_tot = rbind(hailfinder_2_mmhc, hailfinder_2_pc)

###########################################################
hepar2_1_pc$type = rep(c("PC"), 30)
hepar2_1_mmhc$type = rep(c("MMHC"),30)
hepar2_1_tot = rbind(hepar2_1_mmhc, hepar2_1_pc)

hepar2_2_pc$type = rep(c("PC"), 30)
hepar2_2_mmhc$type = rep(c("MMHC"),30)
hepar2_2_tot = rbind(hepar2_2_mmhc, hepar2_2_pc)

###########################################################


link_1_pc$type = rep(c("PC"), 30)
link_1_mmhc$type = rep(c("MMHC"),30)
link_1_tot = rbind(link_1_mmhc, link_1_pc)

link_2_pc$type = rep(c("PC"), 33)
link_2_mmhc$type = rep(c("MMHC"),33)
link_2_tot = rbind(link_2_mmhc, link_2_pc)

###########################################################

pigs_1_pc$type = rep(c("PC"), 30)
pigs_1_mmhc$type = rep(c("MMHC"),30)
pigs_1_tot = rbind(pigs_1_mmhc, pigs_1_pc)

pigs_2_pc$type = rep(c("PC"), 30)
pigs_2_mmhc$type = rep(c("MMHC"),30)
pigs_2_tot = rbind(pigs_2_mmhc, pigs_2_pc)

###########################################################


win95pts_1_pc$type = rep(c("PC"), 30)
win95pts_1_mmhc$type = rep(c("MMHC"),30)
win95pts_1_tot = rbind(win95pts_1_mmhc, win95pts_1_pc)

win95pts_2_pc$type = rep(c("PC"), 30)
win95pts_2_mmhc$type = rep(c("MMHC"),30)
win95pts_2_tot = rbind(win95pts_2_mmhc, win95pts_2_pc)

###########################################################


andes_1_pc$type = rep(c("PC"), 30)
andes_1_mmhc$type = rep(c("MMHC"),30)
andes_1_tot = rbind(andes_1_mmhc, andes_1_pc)

andes_2_pc$type = rep(c("PC"), 30)
andes_2_mmhc$type = rep(c("MMHC"),30)
andes_2_tot = rbind(andes_2_mmhc, andes_2_pc)

#munin1_1_tot$Data = fct_recode(munin1_1_tot$Data, Consensus = "Concensus")
###########################################################
library(ggplot2)
p1 <- ggplot(andes_1_tot, aes(type, f1_score,fill=Data)) + geom_boxplot() + labs(title = "Andes (p = 223)") +  ylab("F1-Score, n = 100") + theme_bw() + theme(legend.position = "none",axis.title.x=element_blank(), axis.text = element_text(size = 18), plot.title = element_text(size=18, face = "bold"),axis.title.y = element_text(size=18, face = "bold"))+ theme(legend.title=element_blank()) + theme(legend.text=element_text(size=18))
p2 <- ggplot(andes_2_tot, aes(type, f1_score,fill=Data)) + geom_boxplot() + labs(title = "Andes (p = 223)") +  ylab("F1-Score, n = 500") + theme_bw()+ theme(legend.position = "none",axis.title.x=element_blank(), axis.text = element_text(size = 18), plot.title = element_text(size=18, face = "bold"),axis.title.y = element_text(size=18, face = "bold"))+ theme(legend.title=element_blank())+ theme(legend.text=element_text(size=18))
p3 <- ggplot(pigs_1_tot, aes(type, f1_score,fill=Data)) + geom_boxplot() + labs(title = "Pigs (p = 441)") +  ylab("F1-Score, n = 100") + theme_bw()+ theme(legend.position = "none",axis.title.x=element_blank(), axis.text = element_text(size = 18), plot.title = element_text(size=18, face = "bold"),axis.title.y = element_text(size=18, face = "bold"))+ theme(legend.title=element_blank())+ theme(legend.text=element_text(size=18))
p4 <- ggplot(pigs_2_tot, aes(type, f1_score,fill=Data)) + geom_boxplot() + labs(title = "Pigs (p = 441)") +  ylab("F1-Score, n = 500") + theme_bw()+ theme(legend.position = "none",axis.title.x=element_blank(), axis.text = element_text(size = 18), plot.title = element_text(size=18, face = "bold"),axis.title.y = element_text(size=18, face = "bold"))+ theme(legend.title=element_blank())+ theme(legend.text=element_text(size=18))
p5 <- ggplot(diabetes_1_tot, aes(type, f1_score,fill=Data)) + geom_boxplot() + labs(title = "Diabetes (p = 413)") +  ylab("F1-Score, n = 100") + theme_bw()+ theme(legend.position = "none",axis.title.x=element_blank(), axis.text = element_text(size = 18), plot.title = element_text(size=18, face = "bold"),axis.title.y = element_text(size=18, face = "bold"))+ theme(legend.title=element_blank())+ theme(legend.text=element_text(size=18))
p6 <- ggplot(diabetes_2_tot, aes(type, f1_score,fill=Data)) + geom_boxplot() + labs(title = "Diabetes (p = 413)") +  ylab("F1-Score, n = 500") + theme_bw()+ theme(legend.position = "none",axis.title.x=element_blank(), axis.text = element_text(size = 18), plot.title = element_text(size=18, face = "bold"),axis.title.y = element_text(size=18, face = "bold"))+ theme(legend.title=element_blank())+ theme(legend.text=element_text(size=18))
p7 <- ggplot(link_1_tot, aes(type, f1_score,fill=Data)) + geom_boxplot() + labs(title = "Link (p = 724)") +  ylab("F1-Score, n = 100") + theme_bw()+ theme(legend.position = "none",axis.title.x=element_blank(), axis.text = element_text(size = 18), plot.title = element_text(size=18, face = "bold"),axis.title.y = element_text(size=18, face = "bold"))+ theme(legend.title=element_blank())+ theme(legend.text=element_text(size=18))
p8 <- ggplot(link_2_tot, aes(type, f1_score,fill=Data)) + geom_boxplot() + labs(title = "Link (p = 724)") +  ylab("F1-Score, n = 500") + theme_bw()+ theme(legend.position = "none",axis.title.x=element_blank(), axis.text = element_text(size = 18), plot.title = element_text(size=18, face = "bold"),axis.title.y = element_text(size=18, face = "bold"))+ theme(legend.title=element_blank())+ theme(legend.text=element_text(size=18))


p9 <- ggplot(hepar2_1_tot, aes(type, f1_score,fill=Data)) + geom_boxplot() + labs(title = "Hepar2 (p = 70)") +  ylab("F1-Score, n = 100") + theme_bw() + theme(legend.position = "none",axis.title.x=element_blank(), axis.text = element_text(size = 18), plot.title = element_text(size=18, face = "bold"),axis.title.y = element_text(size=18, face = "bold"))+ theme(legend.title=element_blank()) + theme(legend.text=element_text(size=18))
p10 <- ggplot(hepar2_2_tot, aes(type, f1_score,fill=Data)) + geom_boxplot() + labs(title = "Hepar2 (p = 70)") +  ylab("F1-Score, n = 500") + theme_bw()+ theme(legend.position = "none",axis.title.x=element_blank(), axis.text = element_text(size = 18), plot.title = element_text(size=18, face = "bold"),axis.title.y = element_text(size=18, face = "bold"))+ theme(legend.title=element_blank())+ theme(legend.text=element_text(size=18))
p11 <- ggplot(win95pts_1_tot, aes(type, f1_score,fill=Data)) + geom_boxplot() + labs(title = "Win95pts (p = 76)") +  ylab("F1-Score, n = 100") + theme_bw()+ theme(legend.position = "none",axis.title.x=element_blank(), axis.text = element_text(size = 18), plot.title = element_text(size=18, face = "bold"),axis.title.y = element_text(size=18, face = "bold"))+ theme(legend.title=element_blank())+ theme(legend.text=element_text(size=18))
p12 <- ggplot(win95pts_2_tot, aes(type, f1_score,fill=Data)) + geom_boxplot() + labs(title = "Win95pts (p = 76)") +  ylab("F1-Score, n = 500") + theme_bw()+ theme(legend.position = "none",axis.title.x=element_blank(), axis.text = element_text(size = 18), plot.title = element_text(size=18, face = "bold"),axis.title.y = element_text(size=18, face = "bold"))+ theme(legend.title=element_blank())+ theme(legend.text=element_text(size=18))
p13 <- ggplot(hailfinder_1_tot, aes(type, f1_score,fill=Data)) + geom_boxplot() + labs(title = "Hailfinder (p = 56)") +  ylab("F1-Score, n = 100") + theme_bw()+ theme(legend.position = "none",axis.title.x=element_blank(), axis.text = element_text(size = 18), plot.title = element_text(size=18, face = "bold"),axis.title.y = element_text(size=18, face = "bold"))+ theme(legend.title=element_blank())+ theme(legend.text=element_text(size=18))
p14 <- ggplot(hailfinder_2_tot, aes(type, f1_score,fill=Data)) + geom_boxplot() + labs(title = "Hailfinder (p = 56)") +  ylab("F1-Score, n = 500") + theme_bw()+ theme(legend.position = "none",axis.title.x=element_blank(), axis.text = element_text(size = 18), plot.title = element_text(size=18, face = "bold"),axis.title.y = element_text(size=18, face = "bold"))+ theme(legend.title=element_blank())+ theme(legend.text=element_text(size=18))
p15 <- ggplot(munin1_1_tot, aes(type, f1_score,fill=Data)) + geom_boxplot() + labs(title = "Munin1 (p = 186)") +  ylab("F1-Score, n = 100") + theme_bw()+ theme(legend.position = "none",axis.title.x=element_blank(), axis.text = element_text(size = 18), plot.title = element_text(size=18, face = "bold"),axis.title.y = element_text(size=18, face = "bold"))+ theme(legend.title=element_blank())+ theme(legend.text=element_text(size=18))
p16 <- ggplot(munin1_2_tot, aes(type, f1_score,fill=Data)) + geom_boxplot() + labs(title = "Munin1 (p = 186)") +  ylab("F1-Score, n = 500") + theme_bw()+ theme(legend.position = "none",axis.title.x=element_blank(), axis.text = element_text(size = 18), plot.title = element_text(size=18, face = "bold"),axis.title.y = element_text(size=18, face = "bold"))+ theme(legend.title=element_blank())+ theme(legend.text=element_text(size=18))


combined = p1 + p3 + p5 +  p7 + p2 + p4 + p6 + p8 +  plot_layout(ncol = 4) & theme(legend.position = "bottom")
combined + plot_layout(guides = "collect")


combined2 = p9 +  p11 + p13 + p15 + p10 + p12 + p14 + p16 +  plot_layout(ncol = 4) & theme(legend.position = "bottom")
combined2 + plot_layout(guides = "collect")


rmse_btw_beta = function(beta_list){
  init_beta = matrix(0, nrow = nrow(beta_list[[1]]), ncol = ncol(beta_list[[1]]))
  between_beta = c()
  for(i in 1:10){
    if(i == 1){
      diff = beta_list[[i]] - init_beta
      diff = diff[diff != 0]
      rmse = sqrt(mean(diff^2))
      between_beta = c(between_beta, rmse)
    }
    else{
      diff = c(beta_list[[i]] - beta_list[[i-1]])
      diff = diff[diff != 0]
      rmse = sqrt(mean(diff^2))
      between_beta = c(between_beta, rmse)
    }
  }
  return(between_beta)
}

rmse_1 = rmse_btw_beta(struc_diff_1_mixed$Beta)



library(patchwork)
p3 <- ggplot(f1_data_1_tot, aes(type, f1_score,fill=Data)) + geom_boxplot() + labs(title = "n = 100, p = 100") +  ylab("F1-Score") + theme_bw()+ theme(legend.position = "none",axis.title.x=element_blank(), axis.text = element_text(size = 18), plot.title = element_text(size=18, face = "bold"),axis.title.y = element_text(size=18, face = "bold"))+ theme(legend.title=element_blank())+ theme(legend.text=element_text(size=18))
p4 <- ggplot(f1_data_9_tot, aes(type, f1_score,fill=Data)) + geom_boxplot() + labs(title = "n = 500, p = 1000") +  ylab("F1-Score") + theme_bw()+ theme(legend.position = "none",axis.title.x=element_blank(), axis.text = element_text(size = 18), plot.title = element_text(size=18, face = "bold"),axis.title.y = element_text(size=18, face = "bold"))+ theme(legend.title=element_blank())+ theme(legend.text=element_text(size=18))
p1 <- ggplot(munin1_1_tot, aes(type, f1_score,fill=Data)) + geom_boxplot() + labs(title = "n = 100, Munin(p = 186)") +  ylab("F1-Score") + theme_bw() + theme(legend.position = "none",axis.title.x=element_blank(), axis.text = element_text(size = 18), plot.title = element_text(size=18, face = "bold"),axis.title.y = element_text(size=18, face = "bold"))+ theme(legend.title=element_blank()) + theme(legend.text=element_text(size=18))
p2 <- ggplot(andes_2_tot, aes(type, f1_score,fill=Data)) + geom_boxplot() + labs(title = "n = 500, Andes(p = 223)") +  ylab("F1-Score") + theme_bw()+ theme(legend.position = "none",axis.title.x=element_blank(), axis.text = element_text(size = 18), plot.title = element_text(size=18, face = "bold"),axis.title.y = element_text(size=18, face = "bold"))+ theme(legend.title=element_blank())+ theme(legend.text=element_text(size=18))
combined = p3 + p4 + p1 + p2 + plot_layout(ncol = 2) & theme(legend.position = "bottom")
combined + plot_layout(guides = "collect")


p7 <- ggplot(f1_data_7_tot, aes(type, f1_score,fill=Data)) + geom_boxplot() + labs(title = "n = 500, p = 100") +  ylab("F1-Score") + theme_bw()+ theme(legend.position = "none",axis.title.x=element_blank(), axis.text = element_text(size = 18), plot.title = element_text(size=18, face = "bold"),axis.title.y = element_text(size=18, face = "bold"))+ theme(legend.title=element_blank())+ theme(legend.text=element_text(size=18))
p8 <- ggplot(f1_data_8_tot, aes(type, f1_score,fill=Data)) + geom_boxplot() + labs(title = "Toeplitz, n = 500, p = 500") + ylab("F1-Score") + theme_bw()
p9 <- ggplot(f1_data_9_tot, aes(type, f1_score,fill=Data)) + 
  geom_boxplot() + 
  labs(title = "n = 500, p = 1000") +  
  ylab("F1-Score") + 
  theme_bw() + 
  theme(axis.title.x=element_blank(),axis.text = element_text(size = 18), plot.title = element_text(size=18, face = "bold"),axis.title.y = element_text(size=18, face = "bold")) + theme(legend.title=element_blank())+ theme(legend.text=element_text(size=18))

library(patchwork)
combined = p1 + p3 + p5 + p7 + p2 + p4 + p6 + p8 + plot_layout(ncol = 4) & theme(legend.position = "bottom")
combined + plot_layout(guides = "collect")


# Hailfinder is a Bayesian network designed to forecast severe summer hail 
# in northeastern Colorado
# Hepar2 - Diagnosing various liver disorders
# Win95pts - An expert system for printer troubleshooting in Windows 95
# Andes - Probabilistic model for a tutoring system with physics
# Diabetes - Decision-based network using blood glucose profile for insulin therapy
# Link - Human Pedigree
# Munin - Electromyography: muscle response to a nerve stimulant to a muscle
# Pigs - 



