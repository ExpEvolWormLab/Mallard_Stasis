rm(list=ls())
library(data.table)
library(lme4)
library(plyr)
library(dplyr)
library(ggpubr)
library(RColorBrewer)
library(matrixStats)
 
## File that pre-process the population data

load("Output_files/RData/State_frequencies.RData")

state_freq=read.table("data/State_freq_NGM.txt",sep='\t',h=TRUE)
transition_rates =read.table("data/Final_merged_data_NGM.txt",h=TRUE,sep="\t")

dim(state_freq);dim(transition_rates)
all_data_raw=merge(transition_rates,state_freq[,c("exper_data_id","data_group_name","mean_still","mean_forward","mean_backward")])
dim(all_data_raw)


### Add a new population column
###### Another 2 populations column
all_data_raw$population = as.factor((tstrsplit(all_data_raw$pop_label, split = "L", fixed = TRUE)[[1]]))
all_data_raw$population2 <- as.character(all_data_raw $population)
all_data_raw$population2[all_data_raw$population2%in%c("CA150","CA250","CA350")] <-  "CA50"
all_data_raw$population2[all_data_raw$population2%in%c("CA1100","CA2100","CA3100")] <-  "CA100"
all_data_raw$population2[all_data_raw$population2%in%c("GA150","GA250","GA450")] <-  "GA50"
all_data_raw$population2[all_data_raw$population2%in%vect_WI] <-  "WI"
all_data_raw$population2[!all_data_raw$population2=='WI' & all_data_raw$data_group_name%in%c(paste0("B30",1:9),"B310")]<-  "pops"

all_data_raw$population2 <- as.factor(all_data_raw$population2)
table(all_data_raw$population2)

#### Wild isolates and Monoecious/Dioecious populations
data_WI <- subset(all_data_raw,population2=="WI")
data_MD_herm <- subset(all_data_raw,data_group_name%in%"B310" & population2=="pops" &
tstrsplit(pop_label,"_")[[2]]=="herm")
data_MD_male <- subset(all_data_raw,data_group_name%in%"B310" & population2=="pops" &
                         tstrsplit(pop_label,"_")[[2]]=="male")

data_populations0_herm <- subset(all_data_raw,data_group_name%in%paste0("B30",1:9) & population2=="pops" & tstrsplit(pop_label,"_")[[2]]=="herm")
data_populations0_male <- subset(all_data_raw,data_group_name%in%paste0("B30",1:9) & population2=="pops" & tstrsplit(pop_label,"_")[[2]]=="male")

data_populations_herm=rbind(data_populations0_herm, data_MD_herm, data_WI)
data_populations_male=rbind(data_populations0_male, data_MD_male)

dim(data_populations_herm) #174 measurements, includes WI (32 measurements)
dim(data_populations_male) #136 measurements

data_populations_herm$population=tstrsplit(data_populations_herm$population,"noM")[[1]]
data_populations_male$population=tstrsplit(data_populations_male$population,"noM")[[1]]

write.table(data_populations_herm,"Output_files/txt/data_populations_herm_original_covariates.txt",sep='\t',row.names = FALSE,quote=FALSE)
write.table(data_populations_male,"Output_files/txt/data_populations_male_original_covariates.txt",sep='\t',row.names = FALSE,quote=FALSE)

for(i in c("temperature","rel_humidity","logD")){
  data_populations_herm[,i] <- (data_populations_herm[,i]-data_populations_herm[,i])/sd(data_populations_herm[,i])
  data_populations_male[,i] <- (data_populations_male[,i]-data_populations_male[,i])/sd(data_populations_male[,i])  
}

## A linear model to estimate the mean value for each population

vect_traits_labels <- paste0(rep(c("Still","Forward","Backward"),each=2),'-',c("Forward","Backward","Still","Backward","Still","Forward"))
vect_traits_labels=c(vect_traits_labels,rep("Mean trait frequency",3))
vect_traits_labels_short=c(paste0(rep(c("S","F","B"),each=2),c("F","B","S","B","S","F")),c("Still","Forward","Backward"))

data_populations_male$population <- as.factor(as.character(data_populations_male$population))
data_populations_herm$population <- as.factor(as.character(data_populations_herm$population))


CI_list_male=list()
CI_list_herm=list()
vect_P_traits2=c(vect_P_traits,"mean_still","mean_forward","mean_backward")

#### The plots will be done in another script file
ylim_vect_herm <- cbind(c(-2.2,-2.7,-2.4,-5.5,-1.3,-3.4,-2.5,-1.5,-4),
                        c(-0.6,-1.2,0,-2.9,.8,-1,.5,2.5,-2))
ylim_vect_male <- cbind(c(-2.2,-2.5,-2.4,-5.5,-1.3,-3.4,-2.6,-1,-4),
                   c(0,-.3,0,-2,.8,-.7,.5,2.5,-1))
k=0
for(i in 1:length(vect_P_traits2)){
k=k+1


lmer_mod_herm <- lmer(data_populations_herm[, vect_P_traits2[i]]~  population-1+(1|data_group_name),data= data_populations_herm)
lmer_mod_male <- lmer(data_populations_male[, vect_P_traits2[i]]~  population-1+(1|data_group_name),data= data_populations_male)


temp_CI <- confint(lmer_mod_herm)
temp_CI <- temp_CI[3:nrow(temp_CI),]
CI_list_herm[[k]]=temp_CI

temp_CI <- confint(lmer_mod_male);temp_CI <- temp_CI[3:nrow(temp_CI),]
CI_list_male[[k]]=temp_CI
}

rm(all_data_raw)
save(list=ls(),file="Output_files/RData/Analysis_Cemee_Pop_WI_TR_estimates.RData")


