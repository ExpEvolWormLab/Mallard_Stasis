# Similar thte 05a code - but now restricts the A6140 population to the lines phenotyped in the third CGE (same as CA100 populations)
# Produce Figure S9 B,C

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

load('Output_files/RData/VCV_A6140_third.RData')
load('Output_files/RData/VCV_CA100.RData')

angle_eigenV <- function(x, y) {
	dot.prod <- x %*% y
	norm.x <- norm(x, type = "2")
	norm.y <- norm(y, type = "2")
	theta <- acos(dot.prod/(norm.x * norm.y))
	as.numeric(theta)
}

VCV_mat <- list()
VCV_mat[[1]] <- VCV_mat_A6140_third
for(i in 1:3) VCV_mat[[i+1]] <- VCV_mat_CA100[[i]]

vect_Pops=c("A6140","CA1100","CA2100","CA3100")

####
MCMCtot <- nrow(VCV_mat[[1]]$VCV_Mat)
MCMCsamp <- 1000 
n <- 6 #number of traits
m <- 4 #number of matrices to compare
r <- 3 #number of random effects specified in the model.
traitnames <- vect_P_traits #trait names
Gnames <- vect_Pops

MCMCarray <- array(, c(MCMCsamp, (n^2) * r, m)) #empty array
MCMCarray[, , 1] <- as.matrix(VCV_mat[[1]]$VCV_Mat[sample(1: MCMCtot, MCMCsamp),])
MCMCarray[, , 2] <- as.matrix(VCV_mat[[2]]$VCV_Mat[sample(1: MCMCtot, MCMCsamp),])
MCMCarray[, , 3] <- as.matrix(VCV_mat[[3]]$VCV_Mat[sample(1: MCMCtot, MCMCsamp),])
MCMCarray[, , 4] <- as.matrix(VCV_mat[[4]]$VCV_Mat[sample(1: MCMCtot, MCMCsamp),])

Garray <- array(, c(n, n, m, MCMCsamp))
dimnames(Garray) <- list(traitnames, traitnames, Gnames)
Parray <- array(, c(n, n, m, MCMCsamp))
dimnames(Parray) <- list(traitnames, traitnames, Gnames)

Earray1 <- array(, c(n, n, m, MCMCsamp))
dimnames(Earray1) <- list(traitnames, traitnames, Gnames)
Earray2 <- array(, c(n, n, m, MCMCsamp))
dimnames(Earray2) <- list(traitnames, traitnames, Gnames)

for (i in 1:m) {
	for (j in 1:MCMCsamp) {
		G <- matrix(MCMCarray[j, 1:(n^2), i], ncol = n)
		R1 <- matrix(MCMCarray[j, ((n^2) + 1):((n^2) * 2), i], ncol = n)
		R2 <- matrix(MCMCarray[j, (((n^2) * 2) + 1):((n^2) * 3), i], ncol = n)
		Garray[, , i, j] <- G
		Earray1[, , i, j] <- R1
		Earray2[, , i, j] <- R2	
		Parray[, , i, j] <- G + R1 + R2
	}
}

source('Rcode/functions_tensor.R', chdir = TRUE)


HHGarray <- array(, c(n, n, m, MCMCsamp))
for (k in 1:MCMCsamp) {
	for (j in 1:m) {
		P <- inv.rootP(Parray[, , j, k])
		HHGarray[, , j, k] <- P %*% Garray[, , j, k] %*% P
	}
}

df_for_tensor= rbind(final_A6140_third, final_CA100)

ped_all = rbind(
#first all the lines
data.frame(id = as.character(unique(df_for_tensor$pop_label)), dam = NA, sire = NA,stringsAsFactors=FALSE),
#then all the phenotyped lines
data.frame(id = 1:nrow(df_for_tensor), dam = as.character(df_for_tensor$pop_label), sire = as.character(df_for_tensor$pop_label),stringsAsFactors=FALSE)
)
for(i in 1:3) ped_all[,i]=as.factor(ped_all[,i])

population_for_ped <- data.frame(pop_label= c(as.character(unique(df_for_tensor$pop_label)),as.character(df_for_tensor$pop_label)),stringsAsFactors=FALSE)
population_for_ped$population=NA
for(i in 1:nrow(population_for_ped)) population_for_ped$population[i] = as.character(subset(unique(df_for_tensor[,c("population","pop_label")]),pop_label==population_for_ped$pop_label[i])$population)

rand.Garray <- array(, c(n, n, m, MCMCsamp))
rand.Garray_corrected <- array(, c(n, n, m, MCMCsamp))

dimnames(rand.Garray) <- list(traitnames, traitnames, Gnames)

df_for_tensor$population=as.factor(as.character(df_for_tensor$population))
rm(i)
rm(VCV_mat_A6140,VCV_mat_CA100,VCV_mat);gc()

# Here we save a file that could be used on a server to compute the randomized
# eigentensors that are computationaly demanding.

save(list=c("ped_all","population_for_ped","n","m","Gnames","Garray","df_for_tensor","Earray2","traitnames","nb_trait"),file='Output_files/RData/File_for_parallel_processing_noCA50_A6140_THIRD.RData')

run_parallel_MCMC <- function(i){
  
	library(MCMCglmm)
	library(dae)
  library(data.table)
	load('Output_files/RData/File_for_parallel_processing_noCA50_A6140_THIRD.RData')
	rand.Garray_corrected_parallel <- array(, c(n, n, m))
	nb_trait=6
	
	##### Generate phenotypic values from the pedigree
	A6140.bv <- rbv(subset(ped_all, population_for_ped$population == Gnames[1]), Garray[, , 1, i]/2)
	
	CA1100.bv <- rbv(subset(ped_all, population_for_ped$population == Gnames[2]), Garray[, , 2, i]/2)
	CA2100.bv <- rbv(subset(ped_all, population_for_ped$population == Gnames[3]), Garray[, , 3, i]/2)
	CA3100.bv <- rbv(subset(ped_all, population_for_ped$population == Gnames[4]), Garray[, , 4, i]/2)

  a.pop <- cumsum(as.numeric(tapply(ped_all$id, population_for_ped$population, length)))
	pop.bv <- rbind(A6140.bv,  CA1100.bv, CA2100.bv, CA3100.bv)

	rand.pop.bv <- pop.bv[sample(dim(pop.bv)[1], replace = F), ]
	
	## Here we have to compute the random Garray using the morissey technique and save them in a list
	
	ped_lines_all <- subset(ped_all,!is.na(dam))
	rand.model_MCMC=list()
	for(k in 1:4){

	k_pop=Gnames[k]
	ped_lines_current <- subset(ped_lines_all, df_for_tensor$population==k_pop)
	
	
	if(k ==1 ){
	  sire.bvs <- rand.pop.bv[substring(ped_all$id,1,5)==k_pop & is.na(ped_all$dam),]
	  # Add the proper sire label, that will differ from the row name (because of the shuffling)
	  sire.bvs <- cbind(data.frame(sire=ped_all$id[substring(ped_all$id,1,5)==k_pop & is.na(ped_all$dam)]),sire.bvs)
	}else{
	sire.bvs <- rand.pop.bv[substring(ped_all$id,1,6)==k_pop & is.na(ped_all$dam),]
	# Add the proper sire label, that will differ from the row name (because of the shuffling)
	sire.bvs <- cbind(data.frame(sire=ped_all$id[substring(ped_all$id,1,6)==k_pop & is.na(ped_all$dam)]),sire.bvs)	
	}
	
	ped_lines_current <- merge(ped_lines_current,sire.bvs)[,c(1,4:9)]
	
	#Vectors of E variance	
	z <- t(apply(ped_lines_current[,2:7],1,function(x){rmvnorm(x+rep(0,6),Earray2[,,k,i])}))
	ped_lines_current[,2:7] <- z
	names(ped_lines_current) <- c("pop_label",traitnames)
	phen.var = diag(nb_trait) * diag(var(z))
	prior_mod <- list(G = list(G1 = list(V = phen.var/3, n = nb_trait)), R = list(V = phen.var/3, n = nb_trait))


	rand.model_MCMC.temp <- MCMCglmm(cbind(c(T12, T13, T21, T23, T31, T32)) ~  trait - 1, random = ~us(trait):pop_label ,
		 rcov = ~us(trait):units, 
		 family = rep("gaussian", nb_trait), data = ped_lines_current, prior = prior_mod, verbose = FALSE,nitt=15000, burnin=5000,thin=10)
	
		rand.Garray_corrected_parallel[,,k] <- matrix(posterior.mode(rand.model_MCMC.temp$VCV)[1:36], ncol = n)
	}
	
return(rand.Garray_corrected_parallel)	
}


clust <- makeCluster(20)
param_list=list()
for(i in 1:100) param_list[[i]] <- i
List_output <-parLapply(clust, param_list , run_parallel_MCMC)
stopCluster(clust)


save(list=ls(),file="Output_files/RData/Tensor_processed_noCA50_A6140_THIRD.Rdata")

### End of the parallel thread
rm(list=ls())
gc()

source('Rcode/functions_tensor.R', chdir = TRUE)
load("Output_files/RData/Tensor_processed_noCA50_A6140_THIRD.Rdata")

for(i in 1: MCMCsamp){
	for(k in 1:m){

		if(i%%100!=0) rand.Garray_corrected[,,k,i] <- matrix(List_output[[i%%100]][,,k], ncol = n)
		if(i%%100==0) rand.Garray_corrected[,,k,i] <- matrix(List_output[[100]][,,k], ncol = n)
	}
}
dimnames(rand.Garray_corrected) <- list(traitnames, traitnames, Gnames)
MCMC.covtensor <- covtensor(Garray)

nnonzero <- min(n * (n + 1)/2, m - 1)
MCMC.covtensor.rand <- covtensor(rand.Garray_corrected)

HPD.eT.val <- cbind(HPDinterval(as.mcmc(MCMC.covtensor$MCMC.S.val[, 1:nnonzero]), prob = 0.95),
HPDinterval(as.mcmc(MCMC.covtensor.rand$MCMC.S.val[, 1:nnonzero]), prob = 0.95))
round(HPD.eT.val, 3)

HPD.eT.val_80 <- cbind(HPDinterval(as.mcmc(MCMC.covtensor$MCMC.S.val[, 1:nnonzero]), prob = 0.80),
HPDinterval(as.mcmc(MCMC.covtensor.rand$MCMC.S.val[, 1:nnonzero]), prob = 0.80))
round(HPD.eT.val_80, 3)

HPD.tensor.coord <- array(,c(m,2,nnonzero))
dimnames(HPD.tensor.coord) <- list(Gnames,c("lower","upper"), paste("E",1:nnonzero,sep=" "))
for (i in 1:m){
  for (j in 1:nnonzero){
    HPD.tensor.coord[i,,j] <- HPDinterval(as.mcmc(MCMC.covtensor$MCMC.G.coord[i,j,]),prob=0.95)[1:2]
  }
}


HPD.tensor.coord_80 <- array(,c(m,2,nnonzero))
dimnames(HPD.tensor.coord_80) <- list(Gnames,c("lower","upper"), paste("E",1:nnonzero,sep=" "))
for (i in 1:m){
  for (j in 1:nnonzero){
    HPD.tensor.coord_80[i,,j] <- HPDinterval(as.mcmc(MCMC.covtensor$MCMC.G.coord[i,j,]),prob=0.80)[1:2]
  }
}

#Figure A2
pdf("plots/FigureS9B.pdf")

par(mfrow=c(1,2))
ylim_vect=rbind(c(-1,0),c(-.4,.6))
for (k in 1:2){
plot(1:m,MCMC.covtensor$av.G.coord[,k,],ylab="",xlab="",pch=16,xaxt="n",frame.plot=F,xlim=c(0.5,m+.5),ylim=ylim_vect[k,],main = "",type="n")
  abline(h=0,lty=2)
axis(1,at=1:m,labels=Gnames,las=2)
arrows(1:m,HPD.tensor.coord[,2,k],1:m,HPD.tensor.coord[,1,k],length=0.05,angle=90,code=3)

arrows(1:m,HPD.tensor.coord_80[,2,k],1:m, HPD.tensor.coord_80[,1,k],length=0.05,angle=90,code=0,col="red",lwd=3)

points(1:m,MCMC.covtensor$av.G.coord[,k,],pch=21,bg="grey")
mtext(dimnames(MCMC.covtensor$av.G.coord)[[2]][k],side=3,at=0,font=2)

}

dev.off()

## First Eigentensor's weight
unique(MCMC.covtensor$tensor.summary[,1])[1]/sum(unique(MCMC.covtensor$tensor.summary[,1]))
# 67%
## Second Eigentensor's weight
unique(MCMC.covtensor$tensor.summary[,1])[2]/sum(unique(MCMC.covtensor$tensor.summary[,1]))
#15%

round(MCMC.covtensor$tensor.summary[1:(n*2),2:dim(MCMC.covtensor$tensor.summary)[2]], 3)

# How much variation is explained ?
abs(MCMC.covtensor$tensor.summary[1,2])/sum(abs(MCMC.covtensor$tensor.summary[1:n,2]))
#e11: 90%
abs(MCMC.covtensor$tensor.summary[2,2])/sum(abs(MCMC.covtensor$tensor.summary[1:n,2]))
#e12: 7%
abs(MCMC.covtensor$tensor.summary[3,2])/sum(abs(MCMC.covtensor$tensor.summary[1:n,2]))
#e13: 2%

abs(MCMC.covtensor$tensor.summary[7,2])/sum(abs(MCMC.covtensor$tensor.summary[(n+1):(2*n),2]))
#e21: 58%
abs(MCMC.covtensor$tensor.summary[8,2])/sum(abs(MCMC.covtensor$tensor.summary[(n+1):(2*n),2]))
#e22:37%
abs(MCMC.covtensor$tensor.summary[9,2])/sum(abs(MCMC.covtensor$tensor.summary[(n+1):(2*n),2]))
#e23: 2%

e11 <- c(as.numeric(MCMC.covtensor$tensor.summary[1,3:dim(MCMC.covtensor$tensor.summary)[2]]))
e12 <- c(as.numeric(MCMC.covtensor$tensor.summary[2,3:dim(MCMC.covtensor$tensor.summary)[2]]))
e13 <- c(as.numeric(MCMC.covtensor$tensor.summary[3,3:dim(MCMC.covtensor$tensor.summary)[2]]))

e11.proj <- apply(Garray, 3:4, proj, b = e11)
e12.proj <- apply(Garray, 3:4, proj, b = e12)
e13.proj <- apply(Garray, 3:4, proj, b = e13)
HPD.e11 <- HPDinterval(t(as.mcmc(e11.proj)),prob = 0.95)
HPD.e12 <- HPDinterval(t(as.mcmc(e12.proj)),prob = 0.95)
HPD.e13 <- HPDinterval(t(as.mcmc(e13.proj)),prob = 0.95)

e21 <- c(as.numeric(MCMC.covtensor$tensor.summary[n+1,3:dim(MCMC.covtensor$tensor.summary)[2]]))
e22 <- c(as.numeric(MCMC.covtensor$tensor.summary[n+2,3:dim(MCMC.covtensor$tensor.summary)[2]]))
e23 <- c(as.numeric(MCMC.covtensor$tensor.summary[n+3,3:dim(MCMC.covtensor$tensor.summary)[2]]))

e21.proj <- apply(Garray, 3:4, proj, b = e21)
e22.proj <- apply(Garray, 3:4, proj, b = e22)
e23.proj <- apply(Garray, 3:4, proj, b = e23)
HPD.e21 <- HPDinterval(t(as.mcmc(e21.proj)),prob = 0.95)
HPD.e22 <- HPDinterval(t(as.mcmc(e22.proj)),prob = 0.95)
HPD.e23 <- HPDinterval(t(as.mcmc(e23.proj)),prob = 0.95)


HPD.e11_80 <- HPDinterval(t(as.mcmc(e11.proj)),prob = 0.80)
HPD.e12_80 <- HPDinterval(t(as.mcmc(e12.proj)),prob = 0.80)
HPD.e13_80 <- HPDinterval(t(as.mcmc(e13.proj)),prob = 0.80)

HPD.e21_80 <- HPDinterval(t(as.mcmc(e21.proj)),prob = 0.80)
HPD.e22_80 <- HPDinterval(t(as.mcmc(e22.proj)),prob = 0.80)
HPD.e23_80 <- HPDinterval(t(as.mcmc(e23.proj)),prob = 0.80)


pdf("plots/FigureS9C.pdf")


par(mar=c(5,4,4,2))
plot(1:m,rowMeans(e11.proj),ylab="lambda",xlab="",pch=16,cex = 1,xaxt="n",frame.plot=F,xlim=c(0,7),ylim=c(0,(max(HPD.e11))),las=2,type="n")
axis(1,at=1:m,labels=Gnames,las=2)
arrows(1:m,HPD.e11[,2],1:m,HPD.e11[,1],length=0.1,angle=90,code=3)
arrows(1:m,HPD.e11_80[,2],1:m,HPD.e11_80[,1],length=0.1,angle=90,code=0,lwd=3,col="red")
points(1:m,rowMeans(e11.proj),pch=21,bg="grey")
mtext("e11",side=3,at=0,font=2)


dev.off()

# Now I need to compute confidence intervals for the correlations

env_Gamma <- new.env()
load("Output_files/RData/Analysis_Cemee_Pop_WI_with_gamma_with_Gas_only_Lisbon.RData",envir=env_Gamma)
# Sample in the posterior of the gamma est. MCMC model
vect_sampling <- 1:1000
cor_dist <- array(dim=c(4,length(vect_sampling),6))

for(i in 1:length(vect_sampling)){
# Pick a Gamma matrix from the posterior distribution
temp_vect =(env_Gamma$model_MCMC$Sol[vect_sampling[i],2:22])[c(7:21,1:6)]

rdm.gamma <- matrix(c(temp_vect[16]*2,temp_vect[1:5],
temp_vect[1],temp_vect[17]*2,temp_vect[6:9],
temp_vect[c(2,6)],temp_vect[18]*2,temp_vect[10:12],
temp_vect[c(3,7,10)],temp_vect[19]*2,temp_vect[13:14],
temp_vect[c(4,8,11,13)],temp_vect[20]*2,temp_vect[15],
temp_vect[c(5,9,12,14,15)],temp_vect[21]*2),6,6)

cor_dist[1,i,] <- cor(e11,eigen(rdm.gamma)$vectors)
cor_dist[2,i,] <- cor(e12,eigen(rdm.gamma)$vectors)
cor_dist[3,i,] <- cor(e21,eigen(rdm.gamma)$vectors)
cor_dist[4,i,] <- cor(e22,eigen(rdm.gamma)$vectors)

}
for(k in 1:4){
for(i in 1:ncol(cor_dist[k,,])){
	cor_dist[k,,i] <- sort(abs(cor_dist[k,,i]))
	}
}

library(MetBrewer)
colorblind.friendly("Hiroshige")
vect_col <- met.brewer (n = 6, name = 'Hiroshige')

pdf('plots/EigenT_correlations_with_gamma_THIRD.pdf',w=5)
par(mfrow=c(2,1))
for(k in c(1,2)){
plot(density(cor_dist[k,,6]),xlim=c(0,1),las=1,bty="n",col=vect_col[6],main="",xlab="Correlation")
for(i in 1:5) lines(density(cor_dist[k,,i]),col=vect_col[i])
legend(.2,c(15,6,8,6)[k],expression(y[1], y[2], y[3], y[4], y[5], y[6]),col=vect_col,lwd=1,bty='n',ncol=2,cex=.8)
}
dev.off()

save(list=ls(),file="Output_files/RData/Tensor_processed_noCA50_A6140_THIRD.Rdata")