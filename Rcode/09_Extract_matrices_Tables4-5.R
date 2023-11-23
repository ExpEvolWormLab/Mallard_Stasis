
# This code extract values for Table 4 and Table 5
rm(list = ls())
gc()
library(MCMCglmm)
library(psych)
library(ggplot2)
library(dplyr)
library(gplots)
library(data.table)
library(matrixStats)
library(boot)
library(Rmisc)
library(dae)
library(nlme)
library(parallel)
library(RColorBrewer)

load('Output_files/RData/VCV_A6140.RData')

### Output the matrices
ev_A6=NULL
for(i in 1:1000){
ev_A6=rbind(ev_A6,eigen(matrix(VCV_mat_A6140[[1]]$VCV_Mat[i,1:36],6,6)/2)$values)
}
temp=round(t(cbind(posterior.mode(as.mcmc(ev_A6)),HPDinterval(as.mcmc(ev_A6)))),digits=3)
temp2 = rbind(temp,round(temp[1,]/sum(temp[1,]),digits=3))
temp3 = as.data.frame(rbind(temp2,round(eigen(VCV_mat_A6140[[1]]$G1_mat/2)$vectors,digits=3)))

names(temp3)=paste0("g",c("max",2:6))
rownames(temp3)=c("Eigenvalues","HPD lower","HPD upper","Proportion",vect_P_traits)

# Values of Table 4
write.csv(temp3,file="Output_files/G_mat_tables/EV_A6140.csv")


##### And gamma

rm(list=ls())
load("Output_files/RData/Analysis_Cemee_Pop_WI_with_gamma_with_Gas_only_Lisbon.RData")

ev_gamma=NULL
for(i in 1:5000){
gamma_i=model_MCMC$Sol[i,vect_retained_effects]
temp_vect= gamma_i[c(7:21,1:6)]
gamma_i <- matrix(c(temp_vect[16]*2,temp_vect[1:5],
                  temp_vect[1],temp_vect[17]*2,temp_vect[6:9],
                  temp_vect[c(2,6)],temp_vect[18]*2,temp_vect[10:12],
                  temp_vect[c(3,7,10)],temp_vect[19]*2,temp_vect[13:14],
                  temp_vect[c(4,8,11,13)],temp_vect[20]*2,temp_vect[15],
                  temp_vect[c(5,9,12,14,15)],temp_vect[21]*2),6,6)

ev_gamma=rbind(ev_gamma,eigen(gamma_i)$values)
}

temp=round(t(cbind(posterior.mode(as.mcmc(ev_gamma)),HPDinterval(as.mcmc(ev_gamma)))),digits=3)
temp2 = rbind(temp,round(abs(temp[1,])/sum(abs(temp[1,])),digits=3))
temp3 = as.data.frame(rbind(temp2,round(eigen(gamma)$vectors,digits=3)))

names(temp3)=paste0("y",c(1:6))
rownames(temp3)=c("Eigenvalues","HPD lower","HPD upper","Proportion",vect_P_traits)


write.csv(temp3,file="Output_files/G_mat_tables/EV_Gamma.csv")

