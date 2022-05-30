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

### Script to estimate the convergence of the posterior mean

pdf(file='plots/MCMC_posterior_mean_convergence_A6140.pdf')
all_means=NULL
for(i in 2:1000){

  all_means <- rbind(all_means,as.numeric(colMeans(VCV_mat_A6140[[1]]$VCV_Mat[1:i,])))
}
for(i in 1:108){
  all_means[,i] <- (all_means[,i]-all_means[999,i])/all_means[999,i]
}

plot(all_means[,1],type="l",ylim=c(-.5,.5),las=1,bty="n")
k=0
for(i in c(2:6,8:12,15:18,22:24,29:30,36)){
  k=k+1
  lines(all_means[,i],col=rainbow(22)[k])
}
abline(h=0.05,lty=2)
abline(h=-0.05,lty=2)
abline(v=300)
dev.off()

# ### We will then use only 600.000 iterations.
# 
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
 clust <- makeCluster(25)
# # A6140
# 
param_list=list()
for(i in 1:1000) param_list[[i]] <- list(i=k,temp_final = final_A6140)
List_output <-parLapply(clust, param_list , run_parallel_MCMC)
save(list=ls(),file=paste0("Random_G_Analysis_Cemee_Pop_WI_A6140_.RData"))

df_G1 <- NULL
for(i in 1:length(List_output)){
  df_G1 <- rbind(df_G1,c(List_output[[i]]$G1_mat))
}

rm(List_output);gc()
save(list=ls(),file=paste0("Random_G_Analysis_Cemee_Pop_WI_A6140_G1_matrices.RData"))

## Subset of 50 A6140 RILs
load("Output_files/RData/VCV_A6140.RData")
uniq_A6140 <- as.character(unique(final_A6140$pop_label)) # 188

param_list=list()
for(i in 1:1000) param_list[[i]] <- list(i="A6140",temp_final = subset(final_A6140,   pop_label %in% sample(uniq_A6140,50)))
List_output <-parLapply(clust, param_list , run_parallel_MCMC)
rm(param_list)
save(list=ls(),file=paste0("Random_G_Analysis_Cemee_Pop_WI_A6140_subset.RData"))

df_G1 <- NULL
for(i in 1:length(List_output)){
  df_G1 <- rbind(df_G1,c(List_output[[i]]$G1_mat))
}

rm(List_output);gc()
save(list=ls(),file=paste0("Random_G_Analysis_Cemee_Pop_WI_A6140_subset_G1_matrices.RData"))

## Update the function to remove is_2012
 
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
   
   model_MCMC <- MCMCglmm(cbind(c(T12, T13, T21, T23, T31, T32)) ~ (temperature+rel_humidity+logD)^3  + trait - 1, random = ~us(trait):pop_label + us(trait):date_str,
                          rcov = ~us(trait):units,
                          family = rep("gaussian", nb_trait), data = temp_final, prior = prior_mod, verbose = FALSE,nitt=700000, burnin=100000,thin=2000)
   
   post_dist = posterior.mode(model_MCMC$VCV)
   
   VCV_mat_temp=list(Population = i, N_measurement = nrow(temp_final), G1_mat = matrix(post_dist[1:nb_trait^2],
                                                                                       nb_trait, nb_trait), G2_mat = matrix(post_dist[(nb_trait^2 + 1):(2 * nb_trait^2)], nb_trait, nb_trait),
                     R_mat = matrix(post_dist[(2 * nb_trait^2 + 1):(3 * nb_trait^2)], nb_trait, nb_trait), VCV_Mat = model_MCMC$VCV)
   return(VCV_mat_temp)
   
 }
 
###### Evolved populations - CA 50s
ls2 <- ls()[! ls()%in%c("run_parallel_MCMC","clust")]
rm(list=c(ls2,"ls2"))

load("Output_files/RData/VCV_CA50.RData")
for(k in c("CA150","CA250","CA350")){
  
  param_list=list()
  for(i in 1:1000) param_list[[i]] <- list(i=k,temp_final = subset(final_CA50,population==k))
  List_output <-parLapply(clust, param_list , run_parallel_MCMC)
  save(list=ls(),file=paste0("Random_G_Analysis_Cemee_Pop_WI_",k,"_.RData"))
  
  df_G1 <- NULL
  for(i in 1:length(List_output)){
    df_G1 <- rbind(df_G1,c(List_output[[i]]$G1_mat))
  }
  
  rm(List_output);gc()
  save(list=ls(),file=paste0("Random_G_Analysis_Cemee_Pop_WI_",k,"_G1_matrices.RData"))
  
}

###### Evolved populations - CA 100s
ls2 <- ls()[! ls()%in%c("run_parallel_MCMC","clust")]
rm(list=c(ls2,"ls2"))

load("Output_files/RData/VCV_CA100.RData")
for(k in c("CA1100","CA2100","CA3100")){
  
  param_list=list()
  for(i in 1:1000) param_list[[i]] <- list(i=k,temp_final = subset(final_CA100,population==k))
  List_output <-parLapply(clust, param_list , run_parallel_MCMC)
  save(list=ls(),file=paste0("Random_G_Analysis_Cemee_Pop_WI_",k,"_.RData"))
  
  df_G1 <- NULL
  for(i in 1:length(List_output)){
    df_G1 <- rbind(df_G1,c(List_output[[i]]$G1_mat))
  }
  
  rm(List_output);gc()
  save(list=ls(),file=paste0("Random_G_Analysis_Ccd .emee_Pop_WI_",k,"_G1_matrices.RData"))
  
}


stopCluster(clust)
