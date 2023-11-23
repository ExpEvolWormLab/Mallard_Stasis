# This code computes the Tensor analysis comparing the G matrix of the CA[1-3]50 populations 
# Produce Figure 3A,B,C and Figure S13

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
load('Output_files/RData/VCV_CA50.RData')
load('Output_files/RData/VCV_CA100.RData')

angle_eigenV <- function(x, y) {
	dot.prod <- x %*% y
	norm.x <- norm(x, type = "2")
	norm.y <- norm(y, type = "2")
	theta <- acos(dot.prod/(norm.x * norm.y))
	as.numeric(theta)
}

vect_Pops=c("CA150","CA250","CA350")
VCV_mat <- list()
VCV_mat[[1]] <- VCV_mat_A6140[[1]]
for(i in 1:3) VCV_mat[[i+1]] <-  VCV_mat_CA50[[i]]
for(i in 1:3) VCV_mat[[i+4]] <- VCV_mat_CA100[[i]]

####
MCMCtot <- nrow(VCV_mat[[1]]$VCV_Mat)
MCMCsamp <- 1000 
n <- 6 #number of traits
m <- 3 #number of matrices to compare
r <- 3 #number of random effects specified in the model.
traitnames <- vect_P_traits #trait names
Gnames <- vect_Pops

MCMCarray <- array(, c(MCMCsamp, (n^2) * r, m)) #empty array


for(k in 1:m){
MCMCarray[, , k] <- as.matrix(VCV_mat[[(k+1)]]$VCV_Mat[sample(1: MCMCtot, MCMCsamp),])
}

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

## Create a pedigree

df_for_tensor= final_CA50

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

# Here we save a file that could be used on a server to compute the randomized
# eigentensors that are computationaly demanding.

save(list=c("ped_all","population_for_ped","n","m","Gnames","Garray","df_for_tensor","Earray2","traitnames","nb_trait"),file='Output_files/RData/File_for_parallel_processing_onlyCA50.RData')

run_parallel_MCMC <- function(i){

	library(MCMCglmm)
	library(dae)

  load('Output_files/RData/File_for_parallel_processing_onlyCA50.RData')
  
	
	rand.Garray_corrected_parallel <- array(, c(n, n, m))
	
	CA150.bv <- rbv(subset(ped_all, population_for_ped$population == Gnames[1]), Garray[, , 1, i]/2)
	CA250.bv <- rbv(subset(ped_all, population_for_ped$population == Gnames[2]), Garray[, , 2, i]/2)
	CA350.bv <- rbv(subset(ped_all, population_for_ped$population == Gnames[3]), Garray[, , 3, i]/2)
	
	a.pop <- cumsum(as.numeric(tapply(ped_all$id, population_for_ped$population, length)))
	pop.bv <- rbind( CA150.bv, CA250.bv, CA350.bv)

	rand.pop.bv <- pop.bv[sample(dim(pop.bv)[1], replace = F), ]
	
	## Here we have to compute the random Garray using the morissey technique and save them in a list
	
	ped_lines_all <- subset(ped_all,!is.na(dam))
	rand.model_MCMC=list()
	for(k in 1:3){

	k_pop=Gnames[k]
	ped_lines_current <- subset(ped_lines_all, df_for_tensor$population==k_pop)
	if(k %in% 1:3) sire.bvs <- rand.pop.bv[substring(ped_all$id,1,5)==k_pop & is.na(ped_all$dam),]
	if(k %in% 4:6) sire.bvs <- rand.pop.bv[substring(ped_all$id,1,6)==k_pop & is.na(ped_all$dam),]
	# Add the proper sire label, that will differ from the row name (because of the shuffling)
	if(k %in% 1:3) sire.bvs <- cbind(data.frame(sire=ped_all$id[substring(ped_all$id,1,5)==k_pop & is.na(ped_all$dam)]),sire.bvs)
	if(k %in% 4:6) sire.bvs <- cbind(data.frame(sire=ped_all$id[substring(ped_all$id,1,6)==k_pop & is.na(ped_all$dam)]),sire.bvs)	
	
	ped_lines_current <- merge(ped_lines_current,sire.bvs)[,c(1,4:9)]
	#Vectors of E variance	
	z <- t(apply(ped_lines_current[,2:7],1,function(x){rmvnorm(x+rep(0,6),Earray2[,,k,i])}))
	ped_lines_current[,2:7] <- z
	names(ped_lines_current) <- c("pop_label",traitnames)
	phen.var = diag(nb_trait) * diag(var(z))
	prior_mod <- list(G = list(G1 = list(V = phen.var/3, n = nb_trait)), R = list(V = phen.var/3, n = nb_trait))


	rand.model_MCMC.temp <- MCMCglmm(cbind(c(T12, T13, T21, T23, T31, T32)) ~  trait - 1, random = ~us(trait):pop_label ,
		 rcov = ~us(trait):units, 
		 family = rep("gaussian", nb_trait), data = ped_lines_current, prior = prior_mod, verbose = FALSE,nitt=150000, burnin=50000,thin=100)
	
		rand.Garray_corrected_parallel[,,k] <- matrix(posterior.mode(rand.model_MCMC.temp$VCV)[1:36], ncol = n)
	}
	
return(rand.Garray_corrected_parallel)	
}

clust <- makeCluster(22)
param_list=list()
for(i in 1:1000) param_list[[i]] <- i

List_output <-parLapply(clust, param_list, run_parallel_MCMC)

stopCluster(clust)
save(list=ls(),file="Output_files/RData/Tensor_processed_onlyCA50.Rdata")

### End of the parallel thread


rm(list=ls());gc()
#### Back on local computer
source('Rcode/functions_tensor.R', chdir = TRUE)
load('Output_files/RData/Tensor_processed_onlyCA50.Rdata')

for(i in 1:MCMCsamp){
	for(k in 1:m){
		 rand.Garray_corrected[,,k,i] <- matrix(List_output[[i]][,,k], ncol = n)
	}
}
dimnames(rand.Garray_corrected) <- list(traitnames, traitnames, Gnames)
MCMC.covtensor <- covtensor(Garray)
MCMC.covtensor$tensor.summary
nnonzero <- min(n * (n + 1)/2, m - 1)
MCMC.covtensor.rand <- covtensor(rand.Garray_corrected)

HPD.eT.val <- cbind(HPDinterval(as.mcmc(MCMC.covtensor$MCMC.S.val[, 1:nnonzero]), prob = 0.95), 
HPDinterval(as.mcmc(MCMC.covtensor.rand$MCMC.S.val[, 1:nnonzero]), prob = 0.95))
HPD.eT.val_80 <- cbind(HPDinterval(as.mcmc(MCMC.covtensor$MCMC.S.val[, 1:nnonzero]), prob = 0.8), 
                       HPDinterval(as.mcmc(MCMC.covtensor.rand$MCMC.S.val[, 1:nnonzero]), prob = 0.8))


round(HPD.eT.val, 3)
#
# Figure A1
pdf("plots/FigureS10A.pdf")
par(mfrow=c(1,1))

plot((1:nnonzero)-0.2,unique(MCMC.covtensor$tensor.summary[1:(n*nnonzero),1]),xlab="",ylab=expression(alpha),pch=16,cex=1,xaxt="n",frame.plot=F,xlim=c(0.5,3.5),ylim=c(0,max(HPD.eT.val)),type="n")
axis(1,at=1:nnonzero,labels=c(paste("E",rep(1:nnonzero),sep="")))

arrows((1:nnonzero)-0.2, HPD.eT.val[,2],(1:nnonzero)-0.2,HPD.eT.val[,1],length=0.1,angle=90,code=3)
arrows((1:nnonzero)-0.2, HPD.eT.val_80[,2],(1:nnonzero)-0.2, HPD.eT.val_80[,1],length=0.1,angle=90,code=0,lwd=2,col="red")

arrows((1:nnonzero)+0.2, HPD.eT.val[,4],(1:nnonzero)+0.2,HPD.eT.val[,3],length=0.1,angle=90,lty=5,code=3)
arrows((1:nnonzero)+0.2, HPD.eT.val_80[,4],(1:nnonzero)+0.2, HPD.eT.val_80[,3],length=0.1,angle=90,code=0,lwd=2,col="orange")


points((1:nnonzero)+0.2, unique(MCMC.covtensor.rand$tensor.summary[1:(n*nnonzero),1]),pch=1,cex=1)
points((1:nnonzero)-0.2,unique(MCMC.covtensor$tensor.summary[1:(n*nnonzero),1]),pch=16)

legend(2.5,.013,legend=c("observed","randomised"),lty=c(1,5),pch=c(16,1),cex=1,bty="n")

dev.off()


HPD.tensor.coord <- array(,c(m,2,nnonzero))
dimnames(HPD.tensor.coord) <- list(Gnames,c("lower","upper"), paste("E",1:2,sep=" "))
for (i in 1:m){
  for (j in 1:nnonzero){
    HPD.tensor.coord[i,,j] <- HPDinterval(as.mcmc(MCMC.covtensor$MCMC.G.coord[i,j,]),prob=0.95)[1:2]
  }
}

#Figure A2
pdf("plots/FigureS10B.pdf")

k=1
plot(1:m,MCMC.covtensor$av.G.coord[,k,],ylab="",xlab="",pch=16,xaxt="n",frame.plot=F,xlim=c(0.5,m+.5),ylim=c(-.3,.5),main = "")
abline(h=0,lty=2)
axis(1,at=1:m,labels=Gnames,las=2)
arrows(1:m,MCMC.covtensor$av.G.coord[,k,],1:m,HPD.tensor.coord[,1,k],length=0.1,angle=90)
arrows(1:m,MCMC.covtensor$av.G.coord[,k,],1:m,HPD.tensor.coord[,2,k],length=0.1,angle=90)
mtext(dimnames(MCMC.covtensor$av.G.coord)[[2]][k],side=3,at=0,font=2)

dev.off()

round(MCMC.covtensor$tensor.summary[1:(n*2),2:dim(MCMC.covtensor$tensor.summary)[2]], 3)

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

## First Eigentensor's weight
unique(MCMC.covtensor$tensor.summary[,1])[1]/sum(unique(MCMC.covtensor$tensor.summary[,1]))
#53%

#Eigenvectors' weight
abs(MCMC.covtensor$tensor.summary[1,2])/sum(abs(MCMC.covtensor$tensor.summary[1:n,2]))
abs(MCMC.covtensor$tensor.summary[2,2])/sum(abs(MCMC.covtensor$tensor.summary[1:n,2]))

pdf("plots/FigureS10C.pdf")

par(mfrow=c(1,2))
#par(mar=c(5,4,4,2))
plot(1:m,rowMeans(e11.proj),ylab="lambda",xlab="",pch=16,cex = 1,xaxt="n",frame.plot=F,xlim=c(0,7),ylim=c(0,(max(HPD.e11))),las=2)
axis(1,at=1:m,labels=Gnames,las=2)
arrows(1:m,rowMeans(e11.proj),1:m,HPD.e11[,1],length=0.1,angle=90)
arrows(1:m,rowMeans(e11.proj),1:m,HPD.e11[,2],length=0.1,angle=90)
mtext("e11",side=3,at=0,font=2)

plot(1:m,rowMeans(e12.proj),ylab="lambda",xlab="",pch=16,cex = 1,xaxt="n",frame.plot=F,xlim=c(0,7),ylim=c(0,(max(HPD.e12))),las=2)
axis(1,at=1:m,labels=Gnames,las=2)
arrows(1:m,rowMeans(e12.proj),1:m,HPD.e12[,1],length=0.1,angle=90)
arrows(1:m,rowMeans(e12.proj),1:m,HPD.e12[,2],length=0.1,angle=90)
mtext("e12",side=3,at=0,font=2)

dev.off()



save(list=ls(),file='~/Projets_Recherches/Celegans/G_matrix_manuscript/2021_01/New_Tensor/Analysis_Cemee_Pop_WI_tensor_onlyCA50.RData')

