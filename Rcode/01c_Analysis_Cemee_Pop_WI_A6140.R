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

#Load populations data
final_merged =read.table("data/Final_merged_data_NGM.txt",h=TRUE,sep="\t")


#We remove the last block which is very specific and will be analyzed separately
final_merged=subset(final_merged,!data_group_name=='B400')

vect_WI <- c("AB1","CB4507", "CB4852","CB4855" , "CB4856" ,"CB4858", "JU319","JU345" , "JU400" , "MY1","MY16","N2anc","PB306" ,"PX174","PX179" , "RC301")
vect_populations <- c("A6140", "CA150", "CA250", "CA350", "CA1100", "CA2100", "CA3100")

final_merged$population = as.factor((tstrsplit(final_merged$pop_label, split = "L", fixed = TRUE)[[1]]))

## We will select only A6140 populations to fit the G separately

table(final_merged$population)
final_merged $year <- substring(final_merged $date_str,1,4)
final_A6140 <-  subset(final_merged,population=="A6140" & location_label=='Lisbon')
final_A6140 $population <- as.factor(as.character(final_A6140 $population))
table(final_A6140$population)

final_A6140$is_2012 <- (final_A6140$year=="2012")
table(final_A6140$is_2012)
table(final_A6140[,c("year","population","location_label")])

for(j in c('temperature',"rel_humidity","logD")){
final_A6140[,j] <- (final_A6140[,j]-mean(final_A6140[,j]))/sd(final_A6140[,j])
}

vect_P_traits <- c("T12","T13","T21","T23","T31","T32")

VCV_mat_A6140 = NULL
nb_trait = length(vect_P_traits)
k=0

plot.acfs <- function(x) {
	n <- dim(x)[2]
	par(mfrow = c(ceiling(n/2), 2), mar = c(3, 2, 3, 0.5))
	for (i in 1:n) {
		acf(x[, i], lag.max = 100, main = colnames(x)[i], ci.type = "ma", ylim = c(-0.15, 0.15))
		grid()
	}
}

	i="A6140"
	temp_final = final_A6140

	phen.var = diag(nb_trait) * diag(var(subset(temp_final, select = vect_P_traits)))
	prior_mod <- list(G = list(G1 = list(V = phen.var/3, n = nb_trait), G2 = list(V = phen.var/3, n = nb_trait)), 
		R = list(V = phen.var/3, n = nb_trait))

	model_MCMC <- MCMCglmm(cbind(c(T12, T13, T21, T23, T31, T32)) ~ (temperature+rel_humidity+logD)^3 + is_2012 + trait - 1, random = ~us(trait):pop_label + us(trait):date_str,
		 rcov = ~us(trait):units, 
		 family = rep("gaussian", nb_trait), data = temp_final, prior = prior_mod, verbose = TRUE,nitt=2100000, burnin=100000,thin=2000)

	pdf(file = paste0("Output_files/auto_corr_plots_MCMCglmm/Model_pdf_MCMC_", 
		i, ".pdf"), height = 20)
	par(mfrow = c(10, 2), mar = c(2, 2, 1, 0))
	plot(model_MCMC$Sol, auto.layout = F)
	dev.off()

	pdf(file = paste0("Output_files/auto_corr_plots_MCMCglmm/Model_pdf_MCMC_autocorr_", 
		i, ".pdf"), height = 10)
	plot.acfs(model_MCMC$Sol)
	dev.off()

	post_dist = posterior.mode(model_MCMC$VCV)
	k=k+1
	VCV_mat_A6140[[k]]=list(Population = i, N_measurement = nrow(temp_final), G1_mat = matrix(post_dist[1:nb_trait^2], 
		nb_trait, nb_trait), G2_mat = matrix(post_dist[(nb_trait^2 + 1):(2 * nb_trait^2)], nb_trait, nb_trait), 
		R_mat = matrix(post_dist[(2 * nb_trait^2 + 1):(3 * nb_trait^2)], nb_trait, nb_trait), VCV_Mat = model_MCMC$VCV)


save(list=ls(),file='Output_files/RData/VCV_A6140.RData')


rm(list=ls())
gc()
load('Output_files/RData/VCV_A6140.RData')
load('Output_files/RData/VCV_CA50.RData')
load('Output_files/RData/VCV_CA100.RData')

### Then we should plot them

v_col = c("chartreuse","cadetblue1", "cornflowerblue", "slateblue2", "violetred1","darkorange","firebrick")

pdf(file='plots/G_mat_ALL.pdf',h=8,w=6)

par(mar=c(5,7,4,2))
vect_Var <- c(2:6,9:12,16:18,23,24,30,1,8,15,22,29,36)
vProb <- .95

plot(c(VCV_mat_A6140[[1]]$G1_mat/2)[vect_Var],c(24:10,6:1),yaxt="n",bty="n",xlim=c(-.08,.18),xlab="Genetic (co-)variances",xaxt="n",type='n',ylab="",cex.lab=1.2)

mtext(side=2,"Phenotypic traits    \n Diagonal                                  Off-diagonal            ",padj=-2,cex=1.2)
lines(c(0,0),c(24.5,8.5))
lines(c(0,0),c(.5,5.5),col="red")

axis(side=1,pos=0)

axis(side=2,at=c(24:10,6:1),labels=c("SF*SB","SF*FS","SF*FB","SF*BS","SF*BF",
"SB*FS","SB*FB","SB*BS","SB*BF","FS*FB",
"FS*BS","FS*BF","FB*BS","FB*BF","BS*BF",
"SF","SB","FS","FB","BS","BF"),las=1)

i=1
temp_95 <- HPDinterval(VCV_mat_A6140[[i]]$VCV_Mat[,1:36]/2,prob=.95)
temp_80 <- HPDinterval(VCV_mat_A6140[[i]]$VCV_Mat[,1:36]/2,prob=.8)

arrows(temp_95[vect_Var,1],c(24:10,6:1)+(.25*(i-1)),temp_95[vect_Var,2],c(24:10,6:1)+(.25*(i-1)),code=3,length=.02,angle=90)
arrows(temp_80[vect_Var,1],c(24:10,6:1)+(.25*(i-1)),
temp_80[vect_Var,2],c(24:10,6:1)+(.25*(i-1)),code=3,length=0,angle=90,lwd=2,col="chartreuse")

points(c(VCV_mat_A6140[[i]]$G1_mat/2)[vect_Var],c(24:10,6:1)+(.25*(i-1)),pch=21,bg="black",cex=.6)


for(i in 1:3){

temp_95 <- HPDinterval(VCV_mat_CA50[[i]]$VCV_Mat[,1:36]/2,prob=.95)
temp_80 <- HPDinterval(VCV_mat_CA50[[i]]$VCV_Mat[,1:36]/2,prob=.8)

arrows(temp_95[vect_Var,1],c(24:10,6:1)+.2+(.05*(i-1)),temp_95[vect_Var,2],c(24:10,6:1)+.2+(.05*(i-1)),code=3,length=.02,angle=90)
arrows(temp_80[vect_Var,1],c(24:10,6:1)+.2+(.05*(i-1)),
temp_80[vect_Var,2],c(24:10,6:1)+.2+(.05*(i-1)),code=3,length=0,angle=90,lwd=2,col=v_col[i+1])

temp_95 <- HPDinterval(VCV_mat_CA100[[i]]$VCV_Mat[,1:36]/2,prob=.95)
temp_80 <- HPDinterval(VCV_mat_CA100[[i]]$VCV_Mat[,1:36]/2,prob=.8)

arrows(temp_95[vect_Var,1],c(24:10,6:1)+.5+(.05*(i-1)),temp_95[vect_Var,2],c(24:10,6:1)+.5+(.05*(i-1)),code=3,length=.02,angle=90)
arrows(temp_80[vect_Var,1],c(24:10,6:1)+.5+(.05*(i-1)),
temp_80[vect_Var,2],c(24:10,6:1)+.5+(.05*(i-1)),code=3,length=0,angle=90,lwd=2,col=v_col[i+4])

points(c(VCV_mat_CA50[[i]]$G1_mat/2)[vect_Var],c(24:10,6:1)+.2+(.05*(i-1)),pch=21,bg="black",cex=.6)
points(c(VCV_mat_CA100[[i]]$G1_mat/2)[vect_Var],c(24:10,6:1)+.5+(.05*(i-1)),pch=21,bg="black",cex=.6)


}

legend(.06,23,c("A6140","CA150","CA250","CA350","CA1100","CA2100","CA3100"),ncol=2,v_col,bty="n")

dev.off()



#### Then traces

VCV_mat <- list()
VCV_mat[[1]] <- VCV_mat_A6140[[1]]
for(i in 1:3) VCV_mat[[i+1]] <-  VCV_mat_CA50[[i]]
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

load("Output_files/RData/Random_G_Analysis_Cemee_Pop_WI_A6140_subset_G1_matrices.RData",envir=temp_rand)
A6140_mat_subset_rand <- temp_rand$df_G1

## And finally the true A6140 subset
load("Output_files/RData/TRUE_G_Analysis_Cemee_Pop_WI_A6140_subset_60_G1_matrices.RData",envir=temp_rand)
A6140_mat_subset <- temp_rand$df_G1
rm(temp_rand);gc()

all_traces <- array(,c(nrow(VCV_mat[[1]]$VCV_Mat),7))
all_traces_rand <- array(,c(nrow(VCV_mat[[1]]$VCV_Mat),7))
traces_A6140_subset <- NULL
traces_A6140_subset_rand <- NULL
for(i in 1:nrow(VCV_mat[[1]]$VCV_Mat)){
		for(k in 1:length(VCV_mat)){

		all_traces[i,k] <-  sum(diag(matrix(VCV_mat[[k]]$VCV_Mat[i,1:36],6,6)/2))
		all_traces_rand[i,k] <-  sum(diag(matrix(VCV_mat_rand[[k]][i,1:36],6,6)/2))
		}
  if(i<=100){
  traces_A6140_subset <- c(traces_A6140_subset,sum(diag((matrix(A6140_mat_subset[i,],6,6)/2))))
  traces_A6140_subset_rand <- c(traces_A6140_subset_rand,sum(diag((matrix(A6140_mat_subset_rand[i,],6,6)/2))))
  }
}

pdf(file="plots/Traces.pdf",w=4)

plot(1:7,colMeans(all_traces),ylim=c(0,.6),bty="n",las=1,xaxt="n",ylab="G-matrix trace",xlab="",xlim=c(0.8,8))
axis(side=1,at=1:7,c("A6140","CA150","CA250","CA350","CA1100","CA2100","CA3100"),las=2)
int_95 <- apply(all_traces,2,function(x){
	HPDinterval(as.mcmc(x))	
})

int_80 <- apply(all_traces,2,function(x){
	HPDinterval(as.mcmc(x),prob=.8)	
})

arrows(1:7,int_95[1,],1:7,int_95[2,],code=3,length=.05,angle=90)
arrows(1:7, int_80[1,],1:7, int_80[2,],code=3,length=0,angle=90,lwd=2,col="firebrick3")
points(1:7,colMeans(all_traces),pch=16)

int_95_rand <- apply(all_traces_rand,2,function(x){
  HPDinterval(as.mcmc(x))	
})

int_80_rand <- apply(all_traces_rand,2,function(x){
  HPDinterval(as.mcmc(x),prob=.8)	
})

arrows(1:7+.25,int_95_rand[1,],1:7+.25,int_95_rand[2,],code=3,length=.05,angle=90)
arrows(1:7+.25, int_80_rand[1,],1:7+.25, int_80_rand[2,],code=3,length=0,angle=90,lwd=2,col="orange")
points(1:7+.25,colMeans(all_traces_rand),pch=16)


legend(2,.6,c("Observed","Null"),lwd=2,col=c("firebrick3","orange"),bty="n",cex=1.2)

dev.off()
############################################################################################
## Additional plot to compare the A6140 matrix with all line / with a subset of 60 lines
############################################################################################

pdf(file="plots/Traces_A6140_subset.pdf",w=4)

plot(1,colMeans(all_traces)[1],ylim=c(0,.75),bty="n",las=1,xaxt="n",ylab="G-matrix trace",xlab="",xlim=c(0.8,3))
axis(side=1,at=1:2,c("A6140 \n all Lines","A6140 \n 60 random lines"),las=1,padj=1)

arrows(1,int_95[1,1],1,int_95[2,1],code=3,length=.05,angle=90)
arrows(1, int_80[1,1],1, int_80[2,1],code=3,length=0,angle=90,lwd=2,col="firebrick3")
points(1,colMeans(all_traces)[1],pch=16)

arrows(1+.25,int_95_rand[1,1],1+.25,int_95_rand[2,1],code=3,length=.05,angle=90)
arrows(1+.25, int_80_rand[1,1],1+.25, int_80_rand[2,1],code=3,length=0,angle=90,lwd=2,col="orange")
points(1+.25,colMeans(all_traces_rand)[1],pch=16)

A6_subset_80 <-  HPDinterval(as.mcmc(traces_A6140_subset),prob=.8)
A6_subset_95 <-  HPDinterval(as.mcmc(traces_A6140_subset),prob=.95)

arrows(2,A6_subset_95[1],2,A6_subset_95[2],code=3,length=.05,angle=90)
arrows(2, A6_subset_80[1],2, A6_subset_80[2],code=3,length=0,angle=90,lwd=2,col="magenta")
points(2,mean(traces_A6140_subset),pch=16)

#A6_subset_80_rand <-  HPDinterval(as.mcmc(traces_A6140_subset_rand),prob=.8)
#A6_subset_95_rand <-  HPDinterval(as.mcmc(traces_A6140_subset_rand),prob=.95)

#arrows(2.25,A6_subset_95_rand[1],2.25,A6_subset_95_rand[2],code=3,length=.05,angle=90)
#arrows(2.25, A6_subset_80_rand[1],2.25, A6_subset_80_rand[2],code=3,length=0,angle=90,lwd=2,col="lightblue")
#points(2.25,mean(traces_A6140_subset_rand),pch=16)

dev.off()


A6_eigens <- eigen(VCV_mat[[1]]$G1_mat/2)$vectors

all_coVar_eigen <- array(,c(nrow(VCV_mat[[1]]$VCV_Mat),7,6))
all_Var_eigen <- array(,c(nrow(VCV_mat[[1]]$VCV_Mat),7,6))

all_rand_Var_eigen <- array(,c(nrow(VCV_mat[[1]]$VCV_Mat),7,6))
all_rand_coVar_eigen <- array(,c(nrow(VCV_mat[[1]]$VCV_Mat),7,6))

A6140_subset_Var_eigen <- array(,c(nrow(A6140_mat_subset),6))
A6140_subset_coVar_eigen <- array(,c(nrow(A6140_mat_subset),6))

A6140_subset_rand_Var_eigen <- array(,c(nrow(A6140_mat_subset),6))
A6140_subset_rand_coVar_eigen <- array(,c(nrow(A6140_mat_subset),6))


for(i in 1:nrow(VCV_mat[[1]]$VCV_Mat)){
		for(k in 1:length(VCV_mat)){

	all_Var_eigen[i,k,] <-  diag(t(A6_eigens)%*%(matrix(VCV_mat[[k]]$VCV_Mat[i,1:36],6,6)/2)%*% A6_eigens)
	all_coVar_eigen[i,k,] <-  
		colSums(abs(t(A6_eigens)%*%(matrix(VCV_mat[[k]]$VCV_Mat[i,1:36],6,6)/2)%*% A6_eigens))

	all_rand_Var_eigen[i,k,] <-  diag(t(A6_eigens)%*%(matrix(VCV_mat_rand[[k]][i,1:36],6,6)/2)%*% A6_eigens)
	all_rand_coVar_eigen[i,k,] <-  
	  colSums(abs(t(A6_eigens)%*%(matrix(VCV_mat_rand[[k]][i,1:36],6,6)/2)%*% A6_eigens))
	
		}
  
if(i<=100){
  A6140_subset_Var_eigen[i,] <- diag(t(A6_eigens)%*%(matrix(A6140_mat_subset[i,1:36],6,6)/2)%*% A6_eigens)
  A6140_subset_coVar_eigen <- colSums(abs(t(A6_eigens)%*%(matrix(A6140_mat_subset[i,1:36],6,6)/2)%*% A6_eigens))
  
  A6140_subset_rand_Var_eigen[i,] <- diag(t(A6_eigens)%*%(matrix(A6140_mat_subset_rand[i,1:36],6,6)/2)%*% A6_eigens)
  A6140_subset_rand_coVar_eigen <- colSums(abs(t(A6_eigens)%*%(matrix(A6140_mat_subset_rand[i,1:36],6,6)/2)%*% A6_eigens))
}
  }

pdf(file="plots/Variance_decomposition_along_EV.pdf")
par(mfrow=c(2,3))
for(j in 1:6){
#plot(1:7,colMeans(all_Var_eigen[,,j]),ylim=c(-.1,.5))

int_95 <- apply(all_Var_eigen[,,j],2,function(x){
	HPDinterval(as.mcmc(x))	
})

int_80 <- apply(all_Var_eigen[,,j],2,function(x){
  HPDinterval(as.mcmc(x),prob=.8)	
})

int_95_rand <- apply(all_rand_Var_eigen[,,j],2,function(x){
  HPDinterval(as.mcmc(x))	
})

int_80_rand <- apply(all_rand_Var_eigen[,,j],2,function(x){
  HPDinterval(as.mcmc(x),prob=.8)	
})

plot(1:7,colMeans(all_Var_eigen[,,j]),ylim=c(min(int_95[1,]),max(int_95[2,])),bty="n",las=1,xaxt="n",ylab="Variance",xlab="",xlim=c(.8,8))
if(j==1) mtext(side=3,expression(g[max]),adj=0,cex=.8)
if(j>1) mtext(side=3,bquote(g[.(j)]),adj=0,cex=.8)


axis(side=1,at=1:7,c("A6140","CA150","CA250","CA350","CA1100","CA2100","CA3100"),las=2)

arrows(1:7,int_95[1,],1:7,int_95[2,],code=3,length=.05,angle=90)
arrows(1:7, int_80[1,],1:7, int_80[2,],code=3,length=0,angle=90,lwd=2,col="firebrick3")
points(1:7,colMeans(all_Var_eigen[,,j]),pch=16)

arrows(1:7+.25,int_95_rand[1,],1:7+.25,int_95_rand[2,],code=3,length=.05,angle=90)
arrows(1:7+.25, int_80_rand[1,],1:7+.25, int_80_rand[2,],code=3,length=0,angle=90,lwd=2,col="orange")
points(1:7+.25,colMeans(all_rand_Var_eigen[,,j]),pch=16)


}

dev.off()

################################################################
#   Additional plot - Variance decomposition each matrix separately
################################################################

temp_Var_eigen <- array(,c(nrow(VCV_mat[[1]]$VCV_Mat),6))
temp_rand_Var_eigen <- array(,c(nrow(VCV_mat[[1]]$VCV_Mat),6))

pdf(file="plots/Variance_decomposition_along_EV_Separately.pdf")
par(mfrow=c(3,3))
for(k in 1:7){

  temp_VCV <- VCV_mat[[k]]
  temp_rand_VCV <- VCV_mat_rand[[k]]
  
  temp_eigen_decomposition <- eigen(temp_VCV$G1_mat/2)
  for(i in 1:nrow(temp_VCV$VCV_Mat)){
  temp_Var_eigen[i,] <-  diag(t(temp_eigen_decomposition$vectors)%*%(matrix(temp_VCV$VCV_Mat[i,1:36],6,6)/2)%*% temp_eigen_decomposition$vectors)
  temp_rand_Var_eigen[i,] <-  diag(t(temp_eigen_decomposition$vectors)%*%(matrix(temp_rand_VCV[i,1:36],6,6)/2)%*% temp_eigen_decomposition$vectors)
  }

  int_95 <- HPDinterval(as.mcmc(temp_Var_eigen))	
  int_80 <- HPDinterval(as.mcmc(temp_Var_eigen),prob=.8)	
  
  int_95_rand <- HPDinterval(as.mcmc(temp_rand_Var_eigen))	
  int_80_rand <- HPDinterval(as.mcmc(temp_rand_Var_eigen),prob=.8)	
  


  plot(1:6,temp_eigen_decomposition$values,ylim=c(min(int_95[,1]),max(int_95[,2])),bty="n",las=1,xaxt="n",ylab="Variance along eigenV",xlab="",xlim=c(.8,8),type="n")  
  axis(side=1,at=1:6,c("gmax","g2","g3","g4","g5","g6"),las=2)
  #arrows(1:6,int_95[,1],1:6,int_95[,2],code=3,length=.05,angle=90)
  #arrows(1:6, int_80[,1],1:6, int_80[,2],code=3,length=0,angle=90,lwd=2,col="firebrick3")
  points(1:6,colMeans(temp_Var_eigen),pch=16)

  arrows(1:6+.25,int_95_rand[,1],1:6+.25,int_95_rand[,2],code=3,length=.05,angle=90)
  arrows(1:6+.25, int_80_rand[,1],1:6+.25, int_80_rand[,2],code=3,length=0,angle=90,lwd=2,col="orange")
  points(1:6+.25,colMeans(temp_rand_Var_eigen),pch=16)
  
}
}
dev.off()





################################################################
#   Additional plot - Variance decomposition A6140 subset
################################################################


pdf(file="plots/Variance_decomposition_along_EV_A6140_subset.pdf")
par(mfrow=c(2,3))
for(j in 1:6){

  int_95 <- apply(all_Var_eigen[,,j],2,function(x){
    HPDinterval(as.mcmc(x))	
  })
  
  int_80 <- apply(all_Var_eigen[,,j],2,function(x){
    HPDinterval(as.mcmc(x),prob=.8)	
  })
  
  int_95_rand <- apply(all_rand_Var_eigen[,,j],2,function(x){
    HPDinterval(as.mcmc(x))	
  })
  
  int_80_rand <- apply(all_rand_Var_eigen[,,j],2,function(x){
    HPDinterval(as.mcmc(x),prob=.8)	
  })
  plot(1,mean(all_Var_eigen[,1,j]),ylim=c(0,max(A6140_subset_Var_eigen[,j])),bty="n",las=1,xaxt="n",ylab="Variance along eigenV",xlab="",xlim=c(.8,2.3))
  axis(side=1,at=1:2,c("A6140 \n all Lines","A6140 \n 60 random lines"),las=1,padj=1)
  
  arrows(1,int_95[1,1],1,int_95[2,1],code=3,length=.05,angle=90)
  arrows(1, int_80[1,1],1, int_80[2,1],code=3,length=0,angle=90,lwd=2,col="firebrick3")
  points(1,mean(all_Var_eigen[,1,j]),pch=16)
  
  arrows(1+.25,int_95_rand[1,1],1+.25,int_95_rand[2,1],code=3,length=.05,angle=90)
  arrows(1+.25, int_80_rand[1,1],1+.25, int_80_rand[2,1],code=3,length=0,angle=90,lwd=2,col="orange")
  points(1+.25,mean(all_rand_Var_eigen[,1,j]),pch=16)

  arrows(2, HPDinterval(as.mcmc(A6140_subset_Var_eigen[,j]),.95)[1],2, HPDinterval(as.mcmc(A6140_subset_Var_eigen[,j]),.95)[2],code=3,length=.05,angle=90)
  arrows(2, HPDinterval(as.mcmc(A6140_subset_Var_eigen[,j]),.8)[1],2, HPDinterval(as.mcmc(A6140_subset_Var_eigen[,j]),.8)[2],code=3,length=0,angle=90,lwd=2,col="magenta")
  points(2,colMeans(A6140_subset_Var_eigen)[j],pch=16)
  
  #arrows(2.25, HPDinterval(as.mcmc(A6140_subset_rand_Var_eigen[,j]),.95)[1],2.25, HPDinterval(as.mcmc(A6140_subset_rand_Var_eigen[,j]),.95)[2],code=3,length=.05,angle=90)
  #arrows(2.25, HPDinterval(as.mcmc(A6140_subset_rand_Var_eigen[,j]),.8)[1],2.25, HPDinterval(as.mcmc(A6140_subset_rand_Var_eigen[,j]),.8)[2],code=3,length=0,angle=90,lwd=2,col="lightblue")
  #points(2.25,colMeans(A6140_subset_rand_Var_eigen)[j],pch=16)
  
}

dev.off()

################################################################
#   A6140 matrix with random
################################################################

pdf(file='plots/G_mat_A6140_with_Randoms.pdf',h=7,w=4.5)

par(mar=c(5,7,4,2))
vect_Var <- c(2:6,9:12,16:18,23,24,30,1,8,15,22,29,36)
vProb <- .95

plot(c(VCV_mat_A6140[[1]]$G1_mat/2)[vect_Var],c(24:10,6:1),yaxt="n",bty="n",xlim=c(-.08,.19),xlab="Genetic (co-)variances",xaxt="n",type='n',ylab="",cex.lab=1.2)

mtext(side=2,"Phenotypic traits    \n Diagonal                                  Off-diagonal            ",padj=-2,cex=1.2)
lines(c(0,0),c(24.5,8.5))
lines(c(0,0),c(.5,5.5),col="red")

axis(side=1,pos=0)

axis(side=2,at=c(24:10,6:1),labels=c("SF*SB","SF*FS","SF*FB","SF*BS","SF*BF",
                                     "SB*FS","SB*FB","SB*BS","SB*BF","FS*FB",
                                     "FS*BS","FS*BF","FB*BS","FB*BF","BS*BF",
                                     "SF","SB","FS","FB","BS","BF"),las=1)

i=1
temp_95 <- HPDinterval(VCV_mat_A6140[[i]]$VCV_Mat[,1:36]/2,prob=.95)
temp_80 <- HPDinterval(VCV_mat_A6140[[i]]$VCV_Mat[,1:36]/2,prob=.8)

temp_95_rand <- HPDinterval(as.mcmc(VCV_mat_rand[[i]][,1:36]/2),prob=.95)
temp_80_rand <- HPDinterval(as.mcmc(VCV_mat_rand[[i]][,1:36]/2),prob=.8)

arrows(temp_95[vect_Var,1],c(24:10,6:1)+(.25*(i-1)),temp_95[vect_Var,2],c(24:10,6:1)+(.25*(i-1)),code=3,length=.02,angle=90)
arrows(temp_80[vect_Var,1],c(24:10,6:1)+(.25*(i-1)),
       temp_80[vect_Var,2],c(24:10,6:1)+(.25*(i-1)),code=3,length=0,angle=90,lwd=2,col="firebrick3")

points(c(VCV_mat_A6140[[i]]$G1_mat/2)[vect_Var],c(24:10,6:1)+(.25*(i-1)),pch=21,bg="black",cex=.6)


arrows(temp_95_rand[vect_Var,1],c(24:10,6:1)+(.25*(i-1)-.3),temp_95_rand[vect_Var,2],c(24:10,6:1)+(.25*(i-1)-.3),code=3,length=.02,angle=90)
arrows(temp_80_rand[vect_Var,1],c(24:10,6:1)+(.25*(i-1)-.3),
       temp_80_rand[vect_Var,2],c(24:10,6:1)+(.25*(i-1)-.3),code=3,length=0,angle=90,lwd=2,col="orange")

points(colMeans(VCV_mat_rand[[i]]/2)[vect_Var],c(24:10,6:1)+(.25*(i-1)-.3),pch=21,bg="black",cex=.6)

legend(.02,23,c("A6140 estimates","NULL distribution \n of post. means"),ncol=1,c("firebrick3","orange"),bty="n")

dev.off()

###

### A6140 - Subset of 50 lines
i=1
pdf(file='plots/G_mat_A6140_subset.pdf',h=7,w=7.2)

plot(c(VCV_mat_A6140[[1]]$G1_mat/2)[vect_Var],c(24:10,6:1),yaxt="n",bty="n",xlim=c(-.08,.19),xlab="Genetic (co-)variances",xaxt="n",type='n',ylab="",cex.lab=1.2)

mtext(side=2,"Phenotypic traits    \n Diagonal                                  Off-diagonal            ",padj=-2,cex=1.2)
lines(c(0,0),c(24.5,8.5))
lines(c(0,0),c(.5,5.5),col="red")

axis(side=1,pos=0)

axis(side=2,at=c(24:10,6:1),labels=c("SF*SB","SF*FS","SF*FB","SF*BS","SF*BF",
                                     "SB*FS","SB*FB","SB*BS","SB*BF","FS*FB",
                                     "FS*BS","FS*BF","FB*BS","FB*BF","BS*BF",
                                     "SF","SB","FS","FB","BS","BF"),las=1)

temp_95 <- HPDinterval(VCV_mat_A6140[[i]]$VCV_Mat[,1:36]/2,prob=.95)
temp_80 <- HPDinterval(VCV_mat_A6140[[i]]$VCV_Mat[,1:36]/2,prob=.8)

arrows(temp_95[vect_Var,1],c(24:10,6:1),temp_95[vect_Var,2],c(24:10,6:1),code=3,length=.05,angle=90)
arrows(temp_80[vect_Var,1],c(24:10,6:1), temp_80[vect_Var,2],c(24:10,6:1),code=3,length=0,angle=90,lwd=2,col="red")

points(c(VCV_mat_A6140[[i]]$G1_mat/2)[vect_Var],c(24:10,6:1)+(.25*(i-1)),pch=21,bg="black",cex=.6)

#temp_95_rand <- HPDinterval(as.mcmc(A6140_mat_subset_rand[,1:36]),prob=.95)/2
#arrows(temp_95_rand[vect_Var,1],c(24:10,6:1)-.3, temp_95_rand[vect_Var,2],c(24:10,6:1)-.3,code=3,length=.05,angle=90,col='grey')

temp_95_subset <- HPDinterval(as.mcmc(A6140_mat_subset),prob=.95)/2
arrows(temp_95_subset[vect_Var,1],c(24:10,6:1)-.3, temp_95_subset[vect_Var,2],c(24:10,6:1)-.3,code=3,length=.05,angle=90,col='grey')

#temp_80_rand <- HPDinterval(as.mcmc(A6140_mat_subset_rand[,1:36]),prob=.8)/2
#arrows(temp_80_rand[vect_Var,1],c(24:10,6:1)-.3, temp_80_rand[vect_Var,2],c(24:10,6:1)-.3,code=3,length=0,angle=90,lwd=2,col="orange")

temp_80_subset <- HPDinterval(as.mcmc(A6140_mat_subset),prob=.8)/2
arrows(temp_80_subset[vect_Var,1],c(24:10,6:1)-.3, temp_80_subset[vect_Var,2],c(24:10,6:1)-.3,code=3,length=0,angle=90,lwd=2,col="orange")

#vectX_rand <- posterior.mode(as.mcmc(A6140_mat_subset_rand[,1:36]))[vect_Var]/2
vectX_subset=colMedians(A6140_mat_subset)[vect_Var]/2
points(vectX_subset,c(24:10,6:1)-.3,pch=8,col="grey",cex=.8)

legend(.1,23,c(paste0("A6140 (188 RILs)"),paste0("A6140 (60 RILs)")),pch=c(16,8),lwd=1)

dev.off()


################




#### The same plot for the CA populations
i=0
pdf(file='plots/G_mat_CA_with_Randoms.pdf',h=7,w=7.2)

par(mar=c(5,7,4,2))
vect_Var <- c(2:6,9:12,16:18,23,24,30,1,8,15,22,29,36)
vProb <- .95

for(ii in 2:7){

plot(c(VCV_mat[[ii]]$G1_mat/2)[vect_Var],c(24:10,6:1),yaxt="n",bty="n",xlim=c(-.08,.19),xlab="Genetic (co-)variances",xaxt="n",type='n',ylab="",cex.lab=1.2)

mtext(side=2,"Phenotypic traits    \n Diagonal                                  Off-diagonal            ",padj=-2,cex=1.2)
lines(c(0,0),c(24.5,8.5))
lines(c(0,0),c(.5,5.5),col="red")

axis(side=1,pos=0)
axis(side=2,at=c(24:10,6:1),labels=c("SF*SB","SF*FS","SF*FB","SF*BS","SF*BF",
                                     "SB*FS","SB*FB","SB*BS","SB*BF","FS*FB",
                                     "FS*BS","FS*BF","FB*BS","FB*BF","BS*BF",
                                     "SF","SB","FS","FB","BS","BF"),las=1)

temp_95 <- HPDinterval(VCV_mat[[ii]]$VCV_Mat[,1:36]/2,prob=.95)
temp_80 <- HPDinterval(VCV_mat[[ii]]$VCV_Mat[,1:36]/2,prob=.8)

temp_95_rand <- HPDinterval(as.mcmc(VCV_mat_rand[[ii]][,1:36]/2),prob=.95)
temp_80_rand <- HPDinterval(as.mcmc(VCV_mat_rand[[ii]][,1:36]/2),prob=.8)

arrows(temp_95[vect_Var,1],c(24:10,6:1)+(.25*(i-1)),temp_95[vect_Var,2],c(24:10,6:1)+(.25*(i-1)),code=3,length=.02,angle=90)
arrows(temp_80[vect_Var,1],c(24:10,6:1)+(.25*(i-1)),
       temp_80[vect_Var,2],c(24:10,6:1)+(.25*(i-1)),code=3,length=0,angle=90,lwd=2,col="red")

points(c(VCV_mat[[ii]]$G1_mat/2)[vect_Var],c(24:10,6:1)+(.25*(i-1)),pch=21,bg="black",cex=.6)

arrows(temp_95_rand[vect_Var,1],c(24:10,6:1)+(.25*(i-1)-.3),temp_95_rand[vect_Var,2],c(24:10,6:1)+(.25*(i-1)-.3),code=3,length=.02,angle=90,col="grey")
arrows(temp_80_rand[vect_Var,1],c(24:10,6:1)+(.25*(i-1)-.3),
       temp_80_rand[vect_Var,2],c(24:10,6:1)+(.25*(i-1)-.3),code=3,length=0,angle=90,lwd=2,col="orange")

points(colMeans(VCV_mat_rand[[ii]]/2)[vect_Var],c(24:10,6:1)+(.25*(i-1)-.3),pch=8,bg="black",cex=.6)

#legend(.02,23,c("A6140 estimates","NULL distribution \n of post. means"),ncol=1,c(v_col[ii],"orange"),bty="n")

}
dev.off()


vect_x <- c(1,3:5,7:9)
#for(i in 2:5) vect_x <- c(vect_x,(vect_x+12*(i-1)))

temp_95_rand_low=NULL
temp_95_rand_high=NULL
temp_rand_means=NULL
temp_TRUE_means=NULL

for(k in 1:7){
  
      temp_95_rand_low <- rbind(temp_95_rand_low,diag(matrix(HPDinterval(as.mcmc(VCV_mat_rand[[k]][,1:36]/2),prob=.95)[,1],6,6)))
      temp_95_rand_high <- rbind(temp_95_rand_high,diag(matrix(HPDinterval(as.mcmc(VCV_mat_rand[[k]][,1:36]/2),prob=.95)[,2],6,6)))

      temp_rand_means <- rbind(temp_rand_means,diag(matrix(colMeans(VCV_mat_rand[[k]][,1:36]/2),6,6)))
      temp_TRUE_means <- rbind(temp_TRUE_means,diag(VCV_mat[[k]]$G1_mat/2))
      
}


plot(1,1,type="n",xlim=c(0,70),ylim=c(0,.15),ylab="Genetic variance",xlab="",xaxt="n",bty="n")
for(i in 1:6){
  
  vect_x <- c(1,3:5,7:9) + 12*(i-1)
  arrows(vect_x,temp_95_rand_low[,i],vect_x,temp_95_rand_high[,i],code=3,length=.05,angle=90,lty=2)  
  points(vect_x,temp_rand_means[,i],pch=21)  

  points(vect_x,temp_TRUE_means[,i],pch=21,bg=v_col)
}

###Create tables

VCV_mat=c(VCV_mat_A6140,VCV_mat_CA50,VCV_mat_CA100)

k=0
for(vect_pop_id in c("A6140","CA150","CA250","CA350","CA1100","CA2100","CA3100")){
  k=k+1
  write.table(matrix(paste0(round(1000*VCV_mat[[k]]$G1_mat)/1000
                            ,"  [",round(1000*matrix(HPDinterval(VCV_mat[[k]]$VCV_Mat)[1:36,1],ncol=6))/1000,";",
                            round(1000*matrix(HPDinterval(VCV_mat[[k]]$VCV_Mat)[1:36,2],ncol=6))/1000,"]"),ncol=6),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE,file=paste0("Output_files/G_mat_tables/", vect_pop_id,".txt"))
  
}



