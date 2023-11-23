# This code computes angles between the gmax vector of A6140 and the first three eigenvectors of the CA[1-3]100 populations then produce Figure 2D

rm(list=ls());gc()
library(MCMCglmm)
load('Output_files/RData/VCV_A6140.RData')
load('Output_files/RData/VCV_CA50.RData')
load('Output_files/RData/VCV_CA100.RData')

VCV_mat <- list()
VCV_mat[[1]] <- VCV_mat_A6140[[1]]
for(i in 1:3) VCV_mat[[i+1]] <-  VCV_mat_CA100[[i]]
for(i in 1:3) VCV_mat[[i+4]] <- VCV_mat_CA100[[i]]

#Load random matrices
temp_rand=new.env()
VCV_mat_rand=list()

p=0
for(k in vect_populations){
  p=p+1
  load(paste0("Output_files/RData/Random_G_Analysis_Cemee_Pop_WI_",k,"_G1_matrices.RData"),envir=temp_rand)
  VCV_mat_rand[[p]]=temp_rand$df_G1
}
rm(temp_rand);gc()

#########
angle_theta <- function(x, y) {
  dot.prod <- x %*% y
  norm.x <- norm(x, type = "2")
  norm.y <- norm(y, type = "2")
  theta <- 180/pi * as.numeric(acos(dot.prod/(norm.x * norm.y)))
  as.numeric(theta)
}

gmax <- eigen(VCV_mat_A6140[[1]]$G1_mat)$vectors[,1]

vect_rand_piE_CA100 <- NULL
vect_rand_piE_CA2100 <- NULL
vect_rand_piE_CA3100 <- NULL

vect_rand_piE_CA100_NULL <- NULL
vect_rand_piE_CA2100_NULL <- NULL
vect_rand_piE_CA3100_NULL <- NULL

vect_rand_thetaE_CA100 <- NULL
vect_rand_thetaE_CA100_2 <- NULL
vect_rand_thetaE_CA100_3 <- NULL

vect_rand_thetaE_CA2100 <- NULL
vect_rand_thetaE_CA2100_2 <- NULL
vect_rand_thetaE_CA2100_3 <- NULL

vect_rand_thetaE_CA3100 <- NULL
vect_rand_thetaE_CA3100_2 <- NULL
vect_rand_thetaE_CA3100_3 <- NULL

n_iter=1000
spld_idx <- sample(1:nrow(VCV_mat_A6140[[1]]$VCV_Mat),n_iter)
k=0
for(i in spld_idx){
  k=k+1
  
  # PI
  temp_CA100 <- matrix(VCV_mat_CA100[[1]]$VCV_Mat[i,1:36],6,6)
  temp_CA100_NULL <- matrix(VCV_mat_rand[[5]][k,],6,6)
  
  temp_CA2100 <- matrix(VCV_mat_CA100[[2]]$VCV_Mat[i,1:36],6,6)
  temp_CA2100_NULL <- matrix(VCV_mat_rand[[6]][k,],6,6)
  
  temp_CA3100 <- matrix(VCV_mat_CA100[[3]]$VCV_Mat[i,1:36],6,6)
  temp_CA3100_NULL <- matrix(VCV_mat_rand[[7]][k,],6,6)
  
  temp_gmax_proj <- (t(gmax)%*%(temp_CA100/2)%*%gmax)/sum(gmax^2)
  temp_gmax_proj_NULL <- (t(gmax)%*%(temp_CA100_NULL/2)%*%gmax)/sum(gmax^2)
  vect_rand_piE_CA100 <- c(vect_rand_piE_CA100 , temp_gmax_proj/eigen(temp_CA100/2)$values[1])	
  vect_rand_piE_CA100_NULL <- c(vect_rand_piE_CA100_NULL , temp_gmax_proj_NULL/eigen(temp_CA100_NULL/2)$values[1])	
  
  temp_gmax_proj <- (t(gmax)%*%(temp_CA2100/2)%*%gmax)/sum(gmax^2)
  temp_gmax_proj_NULL <- (t(gmax)%*%(temp_CA2100_NULL/2)%*%gmax)/sum(gmax^2)
  vect_rand_piE_CA2100 <- c(vect_rand_piE_CA2100 , temp_gmax_proj/eigen(temp_CA2100/2)$values[1])	
  vect_rand_piE_CA2100_NULL <- c(vect_rand_piE_CA2100_NULL , temp_gmax_proj_NULL/eigen(temp_CA2100_NULL/2)$values[1])
  
  temp_gmax_proj <- (t(gmax)%*%(temp_CA3100/2)%*%gmax)/sum(gmax^2)
  temp_gmax_proj_NULL <- (t(gmax)%*%(temp_CA3100_NULL/2)%*%gmax)/sum(gmax^2)
  vect_rand_piE_CA3100 <- c(vect_rand_piE_CA3100 , temp_gmax_proj/eigen(temp_CA3100/2)$values[1])	
  vect_rand_piE_CA3100_NULL <- c(vect_rand_piE_CA3100_NULL , temp_gmax_proj_NULL/eigen(temp_CA3100_NULL/2)$values[1])
  
  ## Angles
  vect_rand_thetaE_CA100 <- c(vect_rand_thetaE_CA100,angle_theta(gmax,eigen(temp_CA100/2)$vector[,1]))
  vect_rand_thetaE_CA100_2 <- c(vect_rand_thetaE_CA100_2,angle_theta(gmax,eigen(temp_CA100/2)$vector[,2]))
  vect_rand_thetaE_CA100_3 <- c(vect_rand_thetaE_CA100_3,angle_theta(gmax,eigen(temp_CA100/2)$vector[,3]))
  
  ## Angles
  vect_rand_thetaE_CA2100 <- c(vect_rand_thetaE_CA2100,angle_theta(gmax,eigen(temp_CA2100/2)$vector[,1]))
  vect_rand_thetaE_CA2100_2 <- c(vect_rand_thetaE_CA2100_2,angle_theta(gmax,eigen(temp_CA2100/2)$vector[,2]))
  vect_rand_thetaE_CA2100_3 <- c(vect_rand_thetaE_CA2100_3,angle_theta(gmax,eigen(temp_CA2100/2)$vector[,3]))
  
  ## Angles
  vect_rand_thetaE_CA3100 <- c(vect_rand_thetaE_CA3100,angle_theta(gmax,eigen(temp_CA3100/2)$vector[,1]))
  vect_rand_thetaE_CA3100_2 <- c(vect_rand_thetaE_CA3100_2,angle_theta(gmax,eigen(temp_CA3100/2)$vector[,2]))
  vect_rand_thetaE_CA3100_3 <- c(vect_rand_thetaE_CA3100_3,angle_theta(gmax,eigen(temp_CA3100/2)$vector[,3]))
}


Angle_rand=NULL
for(i in 1:1000){
  Angle_rand=c(Angle_rand,angle_theta(runif(7,min=(-1),max=1), runif(7,min=(-1),max=1)))
}

## All angles range [0-90]
Angle_rand[Angle_rand>90]= 180-Angle_rand[Angle_rand>90]
vect_rand_thetaE_CA100[vect_rand_thetaE_CA100>90] = 180 - vect_rand_thetaE_CA100[vect_rand_thetaE_CA100>90]
vect_rand_thetaE_CA100_2[vect_rand_thetaE_CA100_2>90] = 180 - vect_rand_thetaE_CA100_2[vect_rand_thetaE_CA100_2>90]
vect_rand_thetaE_CA100_3[vect_rand_thetaE_CA100_3>90] = 180 - vect_rand_thetaE_CA100_3[vect_rand_thetaE_CA100_3>90]

vect_rand_thetaE_CA2100[vect_rand_thetaE_CA2100>90] = 180 - vect_rand_thetaE_CA2100[vect_rand_thetaE_CA2100>90]
vect_rand_thetaE_CA2100_2[vect_rand_thetaE_CA2100_2>90] = 180 - vect_rand_thetaE_CA2100_2[vect_rand_thetaE_CA2100_2>90]
vect_rand_thetaE_CA2100_3[vect_rand_thetaE_CA2100_3>90] = 180 - vect_rand_thetaE_CA2100_3[vect_rand_thetaE_CA2100_3>90]

vect_rand_thetaE_CA3100[vect_rand_thetaE_CA3100>90] = 180 - vect_rand_thetaE_CA3100[vect_rand_thetaE_CA3100>90]
vect_rand_thetaE_CA3100_2[vect_rand_thetaE_CA3100_2>90] = 180 - vect_rand_thetaE_CA3100_2[vect_rand_thetaE_CA3100_2>90]
vect_rand_thetaE_CA3100_3[vect_rand_thetaE_CA3100_3>90] = 180 - vect_rand_thetaE_CA3100_3[vect_rand_thetaE_CA3100_3>90]

vect_col=c("brown1","brown2","brown3")
pdf(file='plots/Figure2D.pdf',w=5,h=7)

plot(1,1,type="n",ylim=c(0,90),xlab="",xlim=c(1.2,2.5),bty="n",las=1,yaxt="n",ylab=expression(paste("Angle with A6140 ",g[max])),xaxt="n")
axis(2,at=c(0,45,90))
mtext(side=3,"B",at=0,cex=2)

temp_95 <- HPDinterval(as.mcmc(vect_rand_thetaE_CA100))
temp_80 <- HPDinterval(as.mcmc(vect_rand_thetaE_CA100),prob=0.83)

arrows(1.5,temp_95[1,1],1.5,temp_95[1,2],code=3,angle=90,length=0.05)
arrows(1.5,temp_80[1,1],1.5,temp_80[1,2],code=3,angle=90,length=0,col=vect_col[1],lwd=2)
points(1.5,mean(vect_rand_thetaE_CA100),pch=16)

temp_95 <- HPDinterval(as.mcmc(vect_rand_thetaE_CA100_2))
temp_80 <- HPDinterval(as.mcmc(vect_rand_thetaE_CA100_2),prob=0.83)

arrows(1.75,temp_95[1,1],1.75,temp_95[1,2],code=3,angle=90,length=0.05)
arrows(1.75,temp_80[1,1],1.75,temp_80[1,2],code=3,angle=90,length=0,col=vect_col[1],lwd=2)
points(1.75,mean(vect_rand_thetaE_CA100_2),pch=16)

temp_95 <- HPDinterval(as.mcmc(vect_rand_thetaE_CA100_3))
temp_80 <- HPDinterval(as.mcmc(vect_rand_thetaE_CA100_3),prob=0.83)

arrows(2,temp_95[1,1],2,temp_95[1,2],code=3,angle=90,length=0.05)
arrows(2,temp_80[1,1],2,temp_80[1,2],code=3,angle=90,length=0,col=vect_col[1],lwd=2)
points(2,mean(vect_rand_thetaE_CA100_3),pch=16)

##### CA 2

temp_95 <- HPDinterval(as.mcmc(vect_rand_thetaE_CA2100))
temp_80 <- HPDinterval(as.mcmc(vect_rand_thetaE_CA2100),prob=0.83)

arrows(1.55,temp_95[1,1],1.55,temp_95[1,2],code=3,angle=90,length=0.05)
arrows(1.55,temp_80[1,1],1.55,temp_80[1,2],code=3,angle=90,length=0,col=vect_col[2],lwd=2)
points(1.55,mean(vect_rand_thetaE_CA2100),pch=16)

temp_95 <- HPDinterval(as.mcmc(vect_rand_thetaE_CA2100_2))
temp_80 <- HPDinterval(as.mcmc(vect_rand_thetaE_CA2100_2),prob=0.83)

arrows(1.8,temp_95[1,1],1.8,temp_95[1,2],code=3,angle=90,length=0.05)
arrows(1.8,temp_80[1,1],1.8,temp_80[1,2],code=3,angle=90,length=0,col=vect_col[2],lwd=2)
points(1.8,mean(vect_rand_thetaE_CA2100_2),pch=16)

temp_95 <- HPDinterval(as.mcmc(vect_rand_thetaE_CA2100_3))
temp_80 <- HPDinterval(as.mcmc(vect_rand_thetaE_CA2100_3),prob=0.83)

arrows(2.05,temp_95[1,1],2.05,temp_95[1,2],code=3,angle=90,length=0.05)
arrows(2.05,temp_80[1,1],2.05,temp_80[1,2],code=3,angle=90,length=0,col=vect_col[2],lwd=2)
points(2.05,mean(vect_rand_thetaE_CA2100_3),pch=16)

#### CA3

temp_95 <- HPDinterval(as.mcmc(vect_rand_thetaE_CA100))
temp_80 <- HPDinterval(as.mcmc(vect_rand_thetaE_CA100),prob=0.83)

arrows(1.6,temp_95[1,1],1.6,temp_95[1,2],code=3,angle=90,length=0.05)
arrows(1.6,temp_80[1,1],1.6,temp_80[1,2],code=3,angle=90,length=0,col=vect_col[3],lwd=2)
points(1.6,mean(vect_rand_thetaE_CA3100),pch=16)

temp_95 <- HPDinterval(as.mcmc(vect_rand_thetaE_CA3100_2))
temp_80 <- HPDinterval(as.mcmc(vect_rand_thetaE_CA3100_2),prob=0.83)

arrows(1.85,temp_95[1,1],1.85,temp_95[1,2],code=3,angle=90,length=0.05)
arrows(1.85,temp_80[1,1],1.85,temp_80[1,2],code=3,angle=90,length=0,col=vect_col[3],lwd=2)
points(1.85,mean(vect_rand_thetaE_CA3100_2),pch=16)

temp_95 <- HPDinterval(as.mcmc(vect_rand_thetaE_CA100_3))
temp_80 <- HPDinterval(as.mcmc(vect_rand_thetaE_CA100_3),prob=0.83)

arrows(2.1,temp_95[1,1],2.1,temp_95[1,2],code=3,angle=90,length=0.05)
arrows(2.1,temp_80[1,1],2.1,temp_80[1,2],code=3,angle=90,length=0,col=vect_col[3],lwd=2)
points(2.1,mean(vect_rand_thetaE_CA3100_3),pch=16)

### RANDOM


temp_95 <- HPDinterval(as.mcmc(Angle_rand))
temp_80 <- HPDinterval(as.mcmc(Angle_rand),prob=0.83)

arrows(1.3,temp_95[1,1],1.3,temp_95[1,2],code=3,angle=90,length=0.05)
arrows(1.3,temp_80[1,1],1.3,temp_80[1,2],code=3,angle=90,length=0,col="orange",lwd=2)

axis(side=1,at=c(2,1.75,1.5,1.3),c(expression(g[3]),expression(g[2]),expression(g[max]),"Null"),las=2)
legend(1.2,20,c("CA1100","CA2100","CA3100"),lwd=3,col=vect_col[1:3])
dev.off()