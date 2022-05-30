rm(list=ls());gc()
library(MCMCglmm)
library(ggplot2)
library(dplyr)
library(data.table)
library(matrixStats)
library(boot)
library(Rmisc)
library(nlme)
library(parallel)

load("Output_files/RData/VCV_A6140.RData")

run_parallel_MCMC <- function(list_param){
  i=list_param[[1]]
  temp_final=list_param[[2]]
  nb_trait=6; vect_P_traits <- c("T12", "T13", "T21", "T23", "T31", "T32")
  
  library(MCMCglmm)
  library(dae)
  
  phen.var = diag(nb_trait) * diag(var(subset(temp_final, select = vect_P_traits)))
  prior_mod <- list(G = list(G1 = list(V = phen.var/3, n = nb_trait), G2 = list(V = phen.var/3, n = nb_trait)),
                    R = list(V = phen.var/3, n = nb_trait))
  
  model_MCMC <- MCMCglmm(cbind(c(T12, T13, T21, T23, T31, T32)) ~ (temperature+rel_humidity+logD)^3 + is_2012 + trait - 1, random = ~us(trait):pop_label + us(trait):date_str,
                         rcov = ~us(trait):units,
                         family = rep("gaussian", nb_trait), data = temp_final, prior = prior_mod, verbose = FALSE,nitt=2100000, burnin=100000,thin=2000)
  
  post_dist = posterior.mode(model_MCMC$VCV)
  
  VCV_mat_temp=list(Population = i, N_measurement = nrow(temp_final), G1_mat = matrix(post_dist[1:nb_trait^2],
                                                                                      nb_trait, nb_trait), G2_mat = matrix(post_dist[(nb_trait^2 + 1):(2 * nb_trait^2)], nb_trait, nb_trait),
                    R_mat = matrix(post_dist[(2 * nb_trait^2 + 1):(3 * nb_trait^2)], nb_trait, nb_trait), VCV_Mat = model_MCMC$VCV)
  return(VCV_mat_temp)
  
}

uniq_A6140 <- as.character(unique(final_A6140$pop_label)) # 188

param_list=list()
for(i in 1:100) param_list[[i]] <- list(i="A6140",temp_final = subset(final_A6140,   pop_label %in% sample(uniq_A6140,60)))
clust <- makeCluster(25)
List_output <-parLapply(clust, param_list , run_parallel_MCMC)
stopCluster(clust)
rm(param_list);rm(run_parallel_MCMC);gc()
save(list=ls(),file=paste0("Output_files/RData/TRUE_G_Analysis_Cemee_Pop_WI_A6140_subset_60.RData"))

df_G1 <- NULL
for(i in 1:length(List_output)){
  df_G1 <- rbind(df_G1,c(List_output[[i]]$G1_mat))
}

rm(List_output);gc()
save(list=ls(),file=paste0("Output_files/RData/TRUE_G_Analysis_Cemee_Pop_WI_A6140_subset_60_G1_matrices.RData"))



