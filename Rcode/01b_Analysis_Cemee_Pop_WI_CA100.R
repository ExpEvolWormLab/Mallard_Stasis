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


## We will select only CA-100 populations to fit the G separately

final_CA100 <-  subset(final_merged,population%in%paste0("CA",1:3,"100"))

for(j in c('temperature',"rel_humidity","logD")){
final_CA100[,j] <- (final_CA100[,j]-mean(final_CA100[,j]))/sd(final_CA100[,j])
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

VCV_mat_CA100 = NULL
nb_trait = length(vect_P_traits)
k=0

for (i in sort(unique(final_CA100$population))) {

	temp_final = subset(final_CA100,   population == i)

	phen.var = diag(nb_trait) * diag(var(subset(temp_final, select = vect_P_traits)))
	prior_mod <- list(G = list(G1 = list(V = phen.var/3, n = nb_trait), G2 = list(V = phen.var/3, n = nb_trait)), 
		R = list(V = phen.var/3, n = nb_trait))

	model_MCMC <- MCMCglmm(cbind(c(T12, T13, T21, T23, T31, T32)) ~ (temperature+rel_humidity+logD)^3 + trait - 1, random = ~us(trait):pop_label + us(trait):date_str,
		 rcov = ~us(trait):units, 
		 family = rep("gaussian", nb_trait), data = temp_final, prior = prior_mod, verbose = FALSE,nitt=2100000, burnin=100000,thin=2000)

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
	VCV_mat_CA100[[k]]=list(Population = i, N_measurement = nrow(temp_final), G1_mat = matrix(post_dist[1:nb_trait^2], 
		nb_trait, nb_trait), G2_mat = matrix(post_dist[(nb_trait^2 + 1):(2 * nb_trait^2)], nb_trait, nb_trait), 
		R_mat = matrix(post_dist[(2 * nb_trait^2 + 1):(3 * nb_trait^2)], nb_trait, nb_trait), VCV_Mat = model_MCMC$VCV)
}


save(list=ls(),file='Output_files/RData/VCV_CA100.RData')


rm(list=ls())
gc()
load('Output_files/RData/VCV_CA100.RData')
### Then we should plot them

pdf(file='plots/G_mat_CA100.pdf',h=8,w=5.5)
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
dev.off()



