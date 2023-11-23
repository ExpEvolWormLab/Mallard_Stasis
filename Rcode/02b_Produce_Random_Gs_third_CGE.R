
# This code produce randomized G matrices of the A6140 matrix for a subset of lines that have been phenotyped as the same time as the 
# CA100 lines (called 'third' common garde here). This is to directly compare populations when phenotyped in the same common garden experiment.

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
load('Output_files/RData/VCV_A6140_third.RData')

# ### We will then use only 600.000 iterations as in the previous code
run_parallel_MCMC <- function(list_param){
  i=list_param[[1]]
  temp_final=list_param[[2]]
  nb_trait=6; vect_P_traits <- c("T12", "T13", "T21", "T23", "T31", "T32")

  library(MCMCglmm)
  library(dae)

  temp_final$pop_label <- sample(temp_final$pop_label)
  temp_final$date_str <- sample(temp_final$date_str)

  phen.var = diag(nb_trait) * diag(var(subset(temp_final, select = vect_P_traits)))
  prior_mod <- list(G = list(G1 = list(V = phen.var/3, n = nb_trait), G2 = list(V = phen.var/3, n = nb_trait)),
                    R = list(V = phen.var/3, n = nb_trait))

  model_MCMC <- MCMCglmm(cbind(c(T12, T13, T21, T23, T31, T32)) ~ (temperature+rel_humidity+logD)^3 + is_2012 + trait - 1, random = ~us(trait):pop_label + us(trait):date_str,
                         rcov = ~us(trait):units,
                         family = rep("gaussian", nb_trait), data = temp_final, prior = prior_mod, verbose = FALSE,nitt=700000, burnin=100000,thin=2000)

  post_dist = posterior.mode(model_MCMC$VCV)

  VCV_mat_temp=list(Population = i, N_measurement = nrow(temp_final), G1_mat = matrix(post_dist[1:nb_trait^2],
                                                                                      nb_trait, nb_trait), G2_mat = matrix(post_dist[(nb_trait^2 + 1):(2 * nb_trait^2)], nb_trait, nb_trait),
                    R_mat = matrix(post_dist[(2 * nb_trait^2 + 1):(3 * nb_trait^2)], nb_trait, nb_trait), VCV_Mat = model_MCMC$VCV)
  return(VCV_mat_temp)

}

dim(final_A6140_third)
clust <- makeCluster(25)
# # A6140
# 
param_list=list()
for(i in 1:1000) param_list[[i]] <- list(i=k,temp_final = final_A6140_third)
List_output <-parLapply(clust, param_list , run_parallel_MCMC)
save(list=ls(),file=paste0("Random_G_Analysis_Cemee_Pop_WI_A6140_THIRD_CGE.RData"))

df_G1 <- NULL
for(i in 1:length(List_output)){
  df_G1 <- rbind(df_G1,c(List_output[[i]]$G1_mat))
}

rm(List_output);gc()
save(list=ls(),file=paste0("Random_G_Analysis_Cemee_Pop_WI_A6140_THIRD_CGE_G1_matrices.RData"))

stopCluster(clust)
