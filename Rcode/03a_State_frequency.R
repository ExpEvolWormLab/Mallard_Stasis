rm(list = ls())
library(data.table)
library(pracma)
library(boot)

# Code for the analysis and the production of Fig. B3

#final_merged =read.table("data/Final_merged_data_NGM.txt",h=TRUE,sep="\t")
state_freq=read.table("data/State_freq_NGM.txt",sep='\t',h=TRUE)

load("Output_files/RData/VCV_A6140.RData")
load("Output_files/RData/VCV_CA50.RData")
load("Output_files/RData/VCV_CA100.RData")
rm(model_MCMC)

final_CA100$year <- substring(final_CA100$date_str,1,4)
final_CA50$is_2012 <- (final_CA50$year=="2012")
final_CA100$is_2012 <- (final_CA100$year=="2012")

final_merged=rbind(final_A6140,final_CA100,final_CA50)

state_freq =subset(state_freq,data_group_name!="B400")
dim(state_freq)
state_final=merge(state_freq, final_merged)


extract_P_inf <- function(vect_exp,nb_sec=3){

	P_mat=c(0, vect_exp[1:3],0, vect_exp[4:6],0)
	P_mat=t(matrix(as.numeric(P_mat),3,3))
	for(i in 1:3) P_mat[i,i]= 0
	for(i in 1:3) P_mat[i,i]= - sum(P_mat[i,])
	

		P_temp=expm(nb_sec*P_mat)

	while(!(sum(round(P_temp[1,],digits=2)==round(P_temp[2,],digits=2))==3 & sum(round(P_temp[1,],digits=2)==round(P_temp[3,],digits=2))==3) & nb_sec<1000){
nb_sec= nb_sec+3
	P_temp=expm(nb_sec*P_mat)
}
	if(nb_sec<100){
		return(expm(nb_sec*P_mat)[1,])

	}else{
		return(c(NA,NA,NA))
		}
}

## Here we use the estimation of the transition rates to model the expected state frequencies
state_from_Stan=apply(as.matrix(round(exp(state_final[, vect_P_traits]),digits=4)),1, extract_P_inf)

vect_exp =as.matrix(round(exp(state_final[, vect_P_traits]),digits=4))[44,]
	P_mat=c(0, vect_exp[1:3],0, vect_exp[4:6],0)
	P_mat=matrix(as.numeric(P_mat),3,3)
	for(i in 1:3) P_mat[i,i]= 0
	for(i in 1:3) P_mat[i,i]= - sum(P_mat[i,])
	
state_from_Stan=t(state_from_Stan)
state_from_Stan_toplot=state_from_Stan[rowSums(state_from_Stan)==1,]
state_final$sex=tstrsplit(state_final$pop_label,"_")[[2]]

### We can plot the observed vs predicted values
pdf("plots/Markov_assumptions_state_freq.pdf")
par(mfrow=c(1,3),bty="n",las=1)

plot(inv.logit(state_final$mean_still)[is.na(state_final$sex)], state_from_Stan[is.na(state_final$sex),1],pch=16,log="",ylab="Frequency of staying still (Markov approximation)",xlab="Observed frequency of staying still",main='Still',xlim=c(0,1),ylim=c(0,1))
abline(a=0,b=1)
plot(inv.logit(state_final$mean_forward)[is.na(state_final$sex)], state_from_Stan[is.na(state_final$sex),2],pch=16,xlim=c(0,1),ylim=c(0,1),ylab="Frequency of moving forward (Markov approximation)",xlab="Observed frequency of moving forward",main='Forward')
abline(a=0,b=1)
plot(inv.logit(state_final$mean_backward)[is.na(state_final$sex)], state_from_Stan[is.na(state_final$sex),3],pch=16,log="",ylab="Frequency of moving backward (Markov approximation)",xlab="Observed frequency of moving backward",main="Backward",xlim=c(0,.2),ylim=c(0,.2))
abline(a=0,b=1)
dev.off()

# Remove duplicates
final_merged <- state_final
rm(state_final,state_freq,final_A6140,final_CA100,final_CA50)
gc()
save(list=ls(),file="Output_files/RData/State_frequencies.RData")

