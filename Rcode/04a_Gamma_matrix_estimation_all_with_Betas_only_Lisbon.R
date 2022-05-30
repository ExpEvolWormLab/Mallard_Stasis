rm(list=ls())
gc()
#load libraries
library(coda)
library(matrixStats)
library(data.table)
library(lme4)
library(MCMCglmm)
library(parallel)

load('Output_files/RData/VCV_A6140.RData')

final_merged$is_2012 <- (final_merged$year=="2012")
final_merged <- subset(final_merged,location_label=='Lisbon')

Lines_means=NULL
for(i in 1:nb_trait){
  
  temp_mod1 <- lmer(final_merged[, vect_P_traits[i]]~ (temperature+rel_humidity+logD)^3 + is_2012 + pop_label-1+(1|date_str),data= final_merged)
  Lines_means <- cbind(Lines_means,as.numeric(summary(temp_mod1)$coef[,1]))
}
row.names(Lines_means)=tstrsplit(names(summary(temp_mod1)$coef[,1]),"pop_label")[[2]]
Lines_means=Lines_means[!is.na(row.names(Lines_means)),]
dim(Lines_means) # 351

mean_phen_values <- data.frame(line=rownames(Lines_means),
                               T12= Lines_means[,1],T13= Lines_means[,2],T21= Lines_means[,3],
                               T23= Lines_means[,4],T31= Lines_means[,5],T32= Lines_means[,6])

#Load the fertility table from Noble 2017
#load("fertility.rda")
load("data/fertility.rda")
fertility <- subset(ecoefs,env=="NGM")
final_fertility <- merge(fertility, mean_phen_values)
dim(final_fertility) # 230

final_fertility$w <- ( exp(final_fertility$fertility)/mean(exp(final_fertility $fertility)))


# We need to center the phenotypic measurements
for(i in vect_P_traits){
  final_fertility[,i] <- final_fertility[,i] - mean(final_fertility[,i])
}

Y=final_fertility$w

vect_retained_effects <- c(2:22)
Y_permuted <- list()
for(i in 1:100) Y_permuted[[i]] <- sample(Y)

final_fertility$population=tstrsplit(final_fertility$line,"L")[[1]]
write.table(final_fertility,file='Output_files/txt/final_fertility_with_GAs_Means_only_Lisbon.txt',sep='\t',row.names=FALSE)

run_random_MCMC <- function(Y_permuted){
  library(MCMCglmm)
  vect_retained_effects <- c(2:22)
  final_fertility <-  read.table('Output_files/txt/final_fertility_with_GAs_Means_only_Lisbon.txt',sep='\t',h=TRUE)
  final_fertility$Y_temp= Y_permuted
  model_MCMC_Rdm <- MCMCglmm(Y_temp ~1+(T12+T13+T21+T23+T31+T32)^2 - (T12+T13+T21+T23+T31+T32) + 
                               I(T12^2)+ I(T13^2)+I(T21^2)+ I(T23^2)+I(T31^2)+ I(T32^2)
                                 , data = final_fertility ,verbose = FALSE,nitt=150000, burnin=100000,thin=10)
  if(length(model_MCMC_Rdm)==20){
    return(posterior.mode(model_MCMC_Rdm$Sol)[vect_retained_effects])
  }else{
    return(NULL)
  }
}

clust <- makeCluster(20)
post_mod_gamma_Rdm <-parLapply(clust, Y_permuted , run_random_MCMC)
stopCluster(clust)

##############################################################################
################## Model with the non-permuted data ##########################
model_MCMC <- MCMCglmm(Y ~1+(T12+T13+T21+T23+T31+T32)^2-(T12+T13+T21+T23+T31+T32)+ I(T12^2)+ I(T13^2)+I(T21^2)+ I(T23^2)+I(T31^2)+ I(T32^2)
                           , data = final_fertility ,verbose = TRUE,nitt= 150000, burnin= 100000 ,thin=10)

##############################################################################
##############################################################################

post_mod_gamma_Rdm_unlist <- unlist(post_mod_gamma_Rdm)
post_mod_gamma_Rdm <- (t(matrix(post_mod_gamma_Rdm_unlist,nrow=21)))
names(posterior.mode(model_MCMC$Sol)[vect_retained_effects])

post_mod_gamma=posterior.mode(model_MCMC$Sol)[vect_retained_effects]
temp_vect= post_mod_gamma[c(7:21,1:6)]
gamma <- matrix(c(temp_vect[16]*2,temp_vect[1:5],
                  temp_vect[1],temp_vect[17]*2,temp_vect[6:9],
                  temp_vect[c(2,6)],temp_vect[18]*2,temp_vect[10:12],
                  temp_vect[c(3,7,10)],temp_vect[19]*2,temp_vect[13:14],
                  temp_vect[c(4,8,11,13)],temp_vect[20]*2,temp_vect[15],
                  temp_vect[c(5,9,12,14,15)],temp_vect[21]*2),6,6)

eigen(gamma)
proj_M = eigen(gamma)$vectors

extract_gamma_EV <- function(temp_vect,temp_proj_M=proj_M){
  temp_vect= temp_vect[c(7:21,1:6)]
  gamma_rdm <- matrix(c(temp_vect[16]*2,temp_vect[1:5],
                        temp_vect[1],temp_vect[17]*2,temp_vect[6:9],
                        temp_vect[c(2,6)],temp_vect[18]*2,temp_vect[10:12],
                        temp_vect[c(3,7,10)],temp_vect[19]*2,temp_vect[13:14],
                        temp_vect[c(4,8,11,13)],temp_vect[20]*2,temp_vect[15],
                        temp_vect[c(5,9,12,14,15)],temp_vect[21]*2),6,6)
  
  # Rotation
  rot_gamma_rdm <- I(t(temp_proj_M)%*% gamma_rdm %*%temp_proj_M)
  # Selection gradient along the trait
  return(diag(rot_gamma_rdm))
}


rdm_EV_along_y=t(apply(post_mod_gamma_Rdm,1, extract_gamma_EV))

## 95% CI of the true gamma
temp_post_mod_gamma=NULL
for(i in sample(1:nrow(model_MCMC$Sol),1000)){
  temp_post_mod_gamma=rbind(temp_post_mod_gamma,model_MCMC$Sol[i,vect_retained_effects])
}

temp_post_mod_gamma <- t(apply(temp_post_mod_gamma,1, extract_gamma_EV))
apply(temp_post_mod_gamma,2,function(x){HPDinterval(as.mcmc(x))})

pdf("plots/Gamma_Means_only_Lisbon/Gamma_EV.pdf")
plot(1:6,eigen(gamma)$values,xlab=expression(paste("Canonical axis of selection (",y[i],")")),ylab=expression(paste("Strength of selection (eigenvalues ",lambda[i],")")),cex.axis=1.2,las=1,cex.lab=1.1,space=0,width=1,xlim=c(0.5,6.5),ylim=c(-13,13),bty="n",pch=8,type="n")
abline(h=0,lty=2)


all_rdm_EV_along_y_90 = t(apply(rdm_EV_along_y,2,function(x){
  HPDinterval(as.mcmc(x),prob=.95)}))

all_rdm_EV_along_y_80 = t(apply(rdm_EV_along_y,2,function(x){
  HPDinterval(as.mcmc(x),prob=.8)}))

arrows((1:6), all_rdm_EV_along_y_90[,1],(1:6), all_rdm_EV_along_y_90[,2],code=3,angle=90,length=.05,col="grey")
arrows((1:6), all_rdm_EV_along_y_80[,1],(1:6), all_rdm_EV_along_y_80[,2],code=3,angle=90,length=0,col="red",lwd=2)
points(c(1:6),colMedians(rdm_EV_along_y),pch=16,col="black")

points(1:6,eigen(gamma)$values,pch=8)
legend(0.3,14,c("Null dist. (median, 80%\n and 95% CI)","True est."),lwd=1,pch=c(16,8),bty="n",cex=.85)
dev.off()

# 95 % CI of the eige of the true gamma
all_gamma_EV = NULL
for(i in 1:nrow(model_MCMC$Sol)){
  
  post_mod_gamma_temp=model_MCMC$Sol[i,vect_retained_effects]
  temp_vect = post_mod_gamma_temp[c(7:21,1:6)]
  gamma_temp <- matrix(c(temp_vect[16]*2,temp_vect[1:5],
                    temp_vect[1],temp_vect[17]*2,temp_vect[6:9],
                    temp_vect[c(2,6)],temp_vect[18]*2,temp_vect[10:12],
                    temp_vect[c(3,7,10)],temp_vect[19]*2,temp_vect[13:14],
                    temp_vect[c(4,8,11,13)],temp_vect[20]*2,temp_vect[15],
                    temp_vect[c(5,9,12,14,15)],temp_vect[21]*2),6,6)
  
  all_gamma_EV=rbind(all_gamma_EV,eigen(gamma_temp)$values)
}

HPDinterval(as.mcmc((all_gamma_EV)),prob=.95)

#         lower       upper
#var1   0.238982860 12.29538220
#var2   0.008448791  1.09398733
#var3  -0.140069543  0.23420940
#var4  -1.031206542  0.03603031
#var5  -3.455695912 -0.46009604
#var6 -18.075600927 -3.89971954




####### NOW THE GAMMA ESTIMATES



### Intervals ###
post_dist_all <-HPDinterval(model_MCMC$Sol,prob=.95)[vect_retained_effects,]
post_dist_all= post_dist_all[c(7:21,1:6),]

post_dist_all2 <-HPDinterval(model_MCMC$Sol,prob=.8)[vect_retained_effects,]
post_dist_all2= post_dist_all2[c(7:21,1:6),]

mode_post_Gamma <- posterior.mode(model_MCMC$Sol)[vect_retained_effects]
mode_post_Gamma = mode_post_Gamma[c(7:21,1:6)]

##########################
pdf("plots/Gamma_Means_only_Lisbon/Gamma_values_no_Betas.pdf")
plot(mode_post_Gamma,c(21:1),yaxt="n",bty="n",ylab="",xlim=c(-6,7),xlab="")

abline(v=0)
axis(side=2,at=21:1,labels=c("SF*SB","SF*FS","SF*FB","SF*BS","SF*BF",
                             "SB*FS","SB*FB","SB*BS","SB*BF","FS*FB",
                             "FS*BS","FS*BF","FB*BS","FB*BF","BS*BF",
                             "SF*SF","SB*SB","FS*FS","FB*FB","BS*BS","BF*BF"),las=1)

arrows(post_dist_all[,1],c(21:1),post_dist_all[,2],c(21:1),code=3,length=.05,angle=90)
arrows(post_dist_all2[,1],c(21:1),post_dist_all2[,2],c(21:1),code=3,length=0,angle=90,lwd=2,col="red")

points(mode_post_Gamma,c(21:1),pch=21,bg="black")
dev.off()



########################
### Gammas + betas   ###
########################

model_MCMC_with_betas <- MCMCglmm(Y ~1+(T12+T13+T21+T23+T31+T32)^2+ I(T12^2)+ I(T13^2)+I(T21^2)+ I(T23^2)+I(T31^2)+ I(T32^2)
                                     , data = final_fertility ,verbose = TRUE,nitt= 150000, burnin= 100000 ,thin=10)

vect_retained_effects_with_betas <- c(2:28)

## Random estimates
run_random_MCMC_with_betas <- function(Y_permuted){
  library(MCMCglmm)
  vect_retained_effects_with_betas <- c(2:28)
  final_fertility <-  read.table('Output_files/txt/final_fertility_with_GAs_Means_only_Lisbon.txt',sep='\t',h=TRUE)
  final_fertility$Y_temp= Y_permuted
  model_MCMC_Rdm <- MCMCglmm(Y_temp ~ 1+(T12+T13+T21+T23+T31+T32)^2+ I(T12^2)+ I(T13^2)+I(T21^2)+ I(T23^2)+I(T31^2)+ I(T32^2)
                                , data = final_fertility ,verbose = FALSE,nitt=150000, burnin=100000,thin=10)
  if(length(model_MCMC_Rdm)==20){
    return(posterior.mode(model_MCMC_Rdm$Sol)[vect_retained_effects_with_betas])
  }else{
    return(NULL)
  }
}

clust <- makeCluster(20)
post_mod_gamma_Rdm_with_betas <-parLapply(clust, Y_permuted , run_random_MCMC_with_betas )
stopCluster(clust)

post_mod_gamma_Rdm_with_betas_unlist <- unlist(post_mod_gamma_Rdm_with_betas)
post_mod_gamma_Rdm_with_betas <- (t(matrix(post_mod_gamma_Rdm_with_betas_unlist,nrow=27)))

### Still needs modifications here
names(posterior.mode(model_MCMC_with_betas$Sol))

post_mod_beta=posterior.mode(model_MCMC_with_betas $Sol)[2:7]
post_mod_gamma_with_betas=posterior.mode(model_MCMC_with_betas $Sol)[vect_retained_effects_with_betas][7:27]
post_mod_gamma_with_betas= post_mod_gamma_with_betas[c(7:21,1:6)]

temp_vect= post_mod_gamma_with_betas[c(7:21,1:6)]

gamma_with_betas <- matrix(c(temp_vect[16]*2,temp_vect[1:5],
                  temp_vect[1],temp_vect[17]*2,temp_vect[6:9],
                  temp_vect[c(2,6)],temp_vect[18]*2,temp_vect[10:12],
                  temp_vect[c(3,7,10)],temp_vect[19]*2,temp_vect[13:14],
                  temp_vect[c(4,8,11,13)],temp_vect[20]*2,temp_vect[15],
                  temp_vect[c(5,9,12,14,15)],temp_vect[21]*2),6,6)


post_dist_all_gamma <-HPDinterval(model_MCMC_with_betas$Sol,prob=.9)[vect_retained_effects_with_betas,][7:27,]
post_dist_all_gamma= post_dist_all_gamma[c(7:21,1:6),]

post_dist_all2_gamma <-HPDinterval(model_MCMC_with_betas$Sol,prob=.8)[vect_retained_effects_with_betas,][7:27,]
post_dist_all2_gamma= post_dist_all2_gamma[c(7:21,1:6),]

post_dist_all_beta <-HPDinterval(model_MCMC_with_betas$Sol,prob=.9)[2:7,]
post_dist_all2_beta <-HPDinterval(model_MCMC_with_betas$Sol,prob=.8)[2:7,]

post_dist_all <- rbind(post_dist_all_beta, post_dist_all_gamma)
post_dist_all2 <- rbind(post_dist_all2_beta, post_dist_all2_gamma)

pdf("plots/Gamma_Means_only_Lisbon/Gamma_values_with_Betas.pdf")
plot(c(post_mod_beta,post_mod_gamma_with_betas),c(27:1),yaxt="n",bty="n",ylab="",xlim=c(-6,7),xlab="")
abline(v=0)
axis(side=2,at=27:1,labels=c("SF","SB","FS","FB","BS","BF",
                             "SF*SB","SF*FS","SF*FB","SF*BS","SF*BF",
                             "SB*FS","SB*FB","SB*BS","SB*BF","FS*FB",
                             "FS*BS","FS*BF","FB*BS","FB*BF","BS*BF",
                             "SF*SF","SB*SB","FS*FS","FB*FB","BS*BS","BF*BF"),las=1)

arrows(post_dist_all[,1],c(27:1),post_dist_all[,2],c(27:1),code=3,length=.05,angle=90)
arrows(post_dist_all2[,1],c(27:1),post_dist_all2[,2],c(27:1),code=3,length=0,angle=90,lwd=2,col="red")
points(c(post_mod_beta,post_mod_gamma_with_betas),c(27:1),pch=21,bg="black")
dev.off()

### And the EV with Betas
eigen(gamma_with_betas)
proj_M_with_betas = eigen(gamma_with_betas)$vectors

rdm_EV_along_y_with_betas=t(apply(post_mod_gamma_Rdm_with_betas,1, extract_gamma_EV,proj_M_with_betas))

pdf("plots/Gamma_Means_only_Lisbon/Gamma_EV_with_betas.pdf")
plot(1:6,eigen(gamma_with_betas)$values,xlab=expression(paste("Canonical axis of selection (",y[i],")")),ylab=expression(paste("Strength of selection (eigenvalues ",lambda[i],")")),cex.axis=1.2,las=1,cex.lab=1.1,space=0,width=1,xlim=c(0.5,6.5),ylim=c(-13,13),bty="n",pch=8,type="n")
abline(h=0,lty=2)

all_rdm_EV_along_y_90 = t(apply(rdm_EV_along_y_with_betas,2,function(x){
  HPDinterval(as.mcmc(x),prob=.95)}))

all_rdm_EV_along_y_80 = t(apply(rdm_EV_along_y_with_betas,2,function(x){
  HPDinterval(as.mcmc(x),prob=.8)}))

arrows((1:6), all_rdm_EV_along_y_90[,1],(1:6), all_rdm_EV_along_y_90[,2],code=3,angle=90,length=.05,col="grey")
arrows((1:6), all_rdm_EV_along_y_80[,1],(1:6), all_rdm_EV_along_y_80[,2],code=3,angle=90,length=0,col="red",lwd=2)
points(c(1:6),colMedians(rdm_EV_along_y_with_betas),pch=16,col="black")

points(1:6,eigen(gamma_with_betas)$values,pch=8)
legend(0.3,14,c("Null dist. (median, 80%\n and 95% CI)","True est."),lwd=1,pch=c(16,8),bty="n",cex=.85)
dev.off()



save(list=ls(),file="Output_files/RData/Analysis_Cemee_Pop_WI_with_gamma_with_Gas_only_Lisbon.RData")

