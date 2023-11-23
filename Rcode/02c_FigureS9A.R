# This code produce the G matrix for the A6140 population restricted to the last CGE - shared with CA[1-3]100 populations
# Then produce Figure S9A

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


## We will select only A6140 populations from the third CGE

final_A6140_third <-  subset(final_merged,population=='A6140' & substring(final_merged$date_str,1,4)=="2018")
dim(final_A6140_third) # 95
length(unique(final_A6140_third$pop_label)) # 63 lines
mean(table(final_A6140_third$pop_label)) # 1.5 obs per population !

for(j in c('temperature',"rel_humidity","logD")){
  final_A6140_third[,j] <- (final_A6140_third[,j]-mean(final_A6140_third[,j]))/sd(final_A6140_third[,j])
}

vect_P_traits <- c("T12","T13","T21","T23","T31","T32")

plot.acfs <- function(x) {
  n <- dim(x)[2]
  par(mfrow = c(ceiling(n/2), 2), mar = c(3, 2, 3, 0.5))
  for (i in 1:n) {
    acf(x[, i], lag.max = 100, main = colnames(x)[i], ci.type = "ma", ylim = c(-0.15, 0.15))
    grid()
  }
}

nb_trait = length(vect_P_traits)

  
  phen.var = diag(nb_trait) * diag(var(subset(final_A6140_third, select = vect_P_traits)))
  prior_mod <- list(G = list(G1 = list(V = phen.var/3, n = nb_trait), G2 = list(V = phen.var/3, n = nb_trait)), 
                    R = list(V = phen.var/3, n = nb_trait))
  
  model_MCMC <- MCMCglmm(cbind(c(T12, T13, T21, T23, T31, T32)) ~ (temperature+rel_humidity+logD)^3 + trait - 1, random = ~us(trait):pop_label + us(trait):date_str,
                         rcov = ~us(trait):units, 
                         family = rep("gaussian", nb_trait), data = final_A6140_third, prior = prior_mod, verbose = FALSE,nitt=2100000, burnin=100000,thin=2000)
  
  post_dist = posterior.mode(model_MCMC$VCV)
  
  VCV_mat_A6140_third=list(Population = "A6140", N_measurement = nrow(final_A6140_third), G1_mat = matrix(post_dist[1:nb_trait^2], 
                                                                                           nb_trait, nb_trait), G2_mat = matrix(post_dist[(nb_trait^2 + 1):(2 * nb_trait^2)], nb_trait, nb_trait), 
                         R_mat = matrix(post_dist[(2 * nb_trait^2 + 1):(3 * nb_trait^2)], nb_trait, nb_trait), VCV_Mat = model_MCMC$VCV)
  
save(list=ls(),file='Output_files/RData/VCV_A6140_third.RData')


rm(list=ls())
gc()
load('Output_files/RData/VCV_A6140_third.RData')
load('Output_files/RData/VCV_CA100.RData')

### Then we should plot them

#pdf(file='plots/G_mat_CA100.pdf',h=8,w=5.5)
par(mar=c(5,7,4,2))
vect_Var <- c(2:6,9:12,16:18,23,24,30,1,8,15,22,29,36)
vProb <- .95


plot(c(VCV_mat_CA100[[1]]$G1_mat)[vect_Var],c(24:10,6:1),yaxt="n",bty="n",xlim=c(-.25,.50),xlab="Genetic (co-)variances",xaxt="n",type='n',ylab="",cex.lab=1.2)
mtext(side=2,"Phenotypic traits    \n Diagonal                                  Off-diagonal            ",padj=-2,cex=1.2)
lines(c(0,0),c(24.5,8.5))
lines(c(0,0),c(.5,5.5),col="red")

axis(side=1,pos=0)

axis(side=2,at=c(24:10,6:1),labels=c("SF*SB","SF*FS","SF*FB","SF*BS","SF*BF",
                                     "SB*FS","SB*FB","SB*BS","SB*BF","FS*FB",
                                     "FS*BS","FS*BF","FB*BS","FB*BF","BS*BF",
                                     "SF","SB","FS","FB","BS","BF"),las=1)

for(i in 1:3){
  
  temp_95 <- HPDinterval(VCV_mat_CA100[[i]]$VCV_Mat[,1:36],prob=.95)
  arrows(temp_95[vect_Var,1],c(24:10,6:1)+(.25*(i-1)),temp_95[vect_Var,2],c(24:10,6:1)+(.25*(i-1)),code=3,length=.02,angle=90)
  
  temp_80 <- HPDinterval(VCV_mat_CA100[[i]]$VCV_Mat[,1:36],prob=.8)
  arrows(temp_80[vect_Var,1],c(24:10,6:1)+(.25*(i-1)),
         temp_80[vect_Var,2],c(24:10,6:1)+(.25*(i-1)),code=3,length=0,angle=90,lwd=2,col=c("orange","firebrick3","magenta")[i])
  
  points(c(VCV_mat_CA100[[i]]$G1_mat)[vect_Var],c(24:10,6:1)+(.25*(i-1)),pch=21,bg="black",cex=.6)
}


i=4
temp_95 <- HPDinterval(VCV_mat_A6140_third$VCV_Mat[,1:36],prob=.95)
arrows(temp_95[vect_Var,1],c(24:10,6:1)+(.25*(i-1)),temp_95[vect_Var,2],c(24:10,6:1)+(.25*(i-1)),code=3,length=.02,angle=90)

temp_80 <- HPDinterval(VCV_mat_A6140_third$VCV_Mat[,1:36],prob=.8)
arrows(temp_80[vect_Var,1],c(24:10,6:1)+(.25*(i-1)),
       temp_80[vect_Var,2],c(24:10,6:1)+(.25*(i-1)),code=3,length=0,angle=90,lwd=2,col="green")

points(c(VCV_mat_A6140_third$G1_mat)[vect_Var],c(24:10,6:1)+(.25*(i-1)),pch=21,bg="black",cex=.6)

load("Output_files/RData/VCV_A6140.RData")

i=5
temp_95 <- HPDinterval(VCV_mat_A6140[[1]]$VCV_Mat[,1:36],prob=.95)
arrows(temp_95[vect_Var,1],c(24:10,6:1)+(.25*(i-1)),temp_95[vect_Var,2],c(24:10,6:1)+(.25*(i-1)),code=3,length=.02,angle=90)

temp_80 <- HPDinterval(VCV_mat_A6140[[1]]$VCV_Mat[,1:36],prob=.8)
arrows(temp_80[vect_Var,1],c(24:10,6:1)+(.25*(i-1)),
       temp_80[vect_Var,2],c(24:10,6:1)+(.25*(i-1)),code=3,length=0,angle=90,lwd=2,col="green")

points(c(VCV_mat_A6140[[1]]$G1_mat)[vect_Var],c(24:10,6:1)+(.25*(i-1)),pch=21,bg="black",cex=.6)

#dev.off()

## Load the Random G A6140_third


### Sampling in the posterior trace
temp_VCV <- list()
temp_VCV[[1]] <- VCV_mat_A6140[[1]]
temp_VCV[[2]] <- VCV_mat_A6140_third
for(i in 1:3) temp_VCV[[(i+2)]] <- VCV_mat_CA100[[i]]


temp_rand=new.env()
VCV_mat_rand=list()

p=0
for(k in vect_populations[c(1,5:7)]){
  p=p+1
  load(paste0("Output_files/RData/Random_G_Analysis_Cemee_Pop_WI_",k,"_G1_matrices.RData"),envir=temp_rand)
  VCV_mat_rand[[p]]=temp_rand$df_G1
  if(p==1){
    p=p+1
    load("Output_files/RData/Random_G_Analysis_Cemee_Pop_WI_A6140_THIRD_CGE_G1_matrices.RData",envir=temp_rand)
    VCV_mat_rand[[p]]=temp_rand$df_G1
  }
}

all_traces <- array(,c(nrow(temp_VCV[[1]]$VCV_Mat),5))
all_traces_rand <- array(,c(nrow(VCV_mat_rand[[1]]),5))

for(i in 1:nrow(temp_VCV[[1]]$VCV_Mat)){
  for(k in 1:length(temp_VCV)){
    all_traces[i,k] <-  sum(diag(matrix(temp_VCV[[k]]$VCV_Mat[i,1:36],6,6)/2))
    all_traces_rand[i,k] <-  sum(diag(matrix(VCV_mat_rand[[k]][i,1:36],6,6)/2))
  }
}

pdf("plots/FigureS9A.pdf")
par(mfrow=c(1,1))
plot(sum(diag(VCV_mat_A6140[[1]]$G1_mat)/2),xlim=c(1,5.5),ylim=c(0,.65),type='n',bty="n",las=1,ylab="G-matrix trace",xaxt="n",xlab="")
axis(side=1,at=1:5,labels=c("A6140","A6140\nsecond","CA1100","CA2100","CA3100"),las=2)
temp_95 <- HPDinterval(as.mcmc(all_traces),prob=.95)
temp_80 <- HPDinterval(as.mcmc(all_traces),prob=.83)

int_95_rand = HPDinterval(as.mcmc(all_traces_rand/2),prob=.95)	
int_80_rand = HPDinterval(as.mcmc(all_traces_rand/2),prob=.83)	

arrows(1:5,temp_95[,1],1:5,temp_95[,2],code=3,length=.02,angle=90)
arrows(1:5,temp_80[,1],1:5,temp_80[,2],code=3,length=0,angle=90,lwd=2,col="firebrick3")

arrows(1:5+.25,int_95_rand[,1],1:5+.25,int_95_rand[,2],code=3,length=.05,angle=90)
arrows(1:5+.25, int_80_rand[,1],1:5+.25, int_80_rand[,2],code=3,length=0,angle=90,lwd=2,col="orange")
points(1:5+.25,colMeans(all_traces_rand/2),pch=16)

points(1:5,as.numeric(lapply(temp_VCV,function(x){sum(diag(x$G1_mat)/2)})),pch=16)
dev.off()