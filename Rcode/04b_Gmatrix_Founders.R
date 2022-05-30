rm(list=ls());gc()
library(MCMCglmm)
library(psych)
library(ggplot2)
library(dplyr)
library(gplots)
library(data.table)
library(matrixStats)
library(boot)
library(corrplot)
library(Rmisc)
library(nlme)

# This code compute G matrices for different sets of lines/populations including some phenotyping data that are not used in the manuscript.

load("Output_files/RData/Analysis_Cemee_Pop_WI_TR_estimates.RData")

## We need to estimate the G-matrix of the WI
WI_for_G_estimation <- subset(all_data_raw,population2=="WI")

for(j in c('temperature',"rel_humidity","logD")){
WI_for_G_estimation[,j] <- (WI_for_G_estimation[,j]-mean(WI_for_G_estimation[,j]))/sd(WI_for_G_estimation[,j])
}

	phen.var = diag(nb_trait) * diag(var(subset(WI_for_G_estimation, select = vect_P_traits)))
	prior_mod <- list(G = list(G1 = list(V = phen.var/3, n = nb_trait), G2 = list(V = phen.var/3, n = nb_trait)), 
		R = list(V = phen.var/3, n = nb_trait))

	model_MCMC_WI <- MCMCglmm(cbind(c(T12, T13, T21, T23, T31, T32)) ~ (logD +rel_humidity+ temperature)+ trait - 1, random = ~us(trait):pop_label + us(trait):date_str,
		 rcov = ~us(trait):units, 
		 family = rep("gaussian", nb_trait), data = WI_for_G_estimation, prior = prior_mod, verbose = FALSE,nitt=2100000, burnin=100000,thin=2000)

post_dist = posterior.mode(model_MCMC_WI$VCV)
WI_G_mat=matrix(post_dist[1:nb_trait^2],nb_trait, nb_trait)

	WI_VCV_mat=list(Population = "WI", G1_mat = matrix(post_dist[1:nb_trait^2], 
		nb_trait, nb_trait), G2_mat = matrix(post_dist[(nb_trait^2 + 1):(2 * nb_trait^2)], nb_trait, nb_trait), 
		R_mat = matrix(post_dist[(2 * nb_trait^2 + 1):(3 * nb_trait^2)], nb_trait, nb_trait), VCV_Mat = model_MCMC_WI$VCV)


save(list=c("WI_VCV_mat"),file="Output_files/RData/WI_G-matrix-all.RData")