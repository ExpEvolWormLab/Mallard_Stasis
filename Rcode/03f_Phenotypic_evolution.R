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
library(lme4)
library(MuMIn)

vect_P_traits <- c("T12","T13","T21","T23","T31","T32")

data_populations_herm <- read.table("Output_files/txt/data_populations_herm_original_covariates.txt",sep='\t',h=TRUE)
data_populations_male <- read.table("Output_files/txt/data_populations_male_original_covariates.txt",sep='\t',h=TRUE)

# Remove some populations that will not serve in this analysis (D and M populations essentially)
data_populations_herm <- subset(data_populations_herm,!substring(population,1,2)=="CM" & population2=="pops" & is.na(tstrsplit(data_populations_herm$pop_label,"noML")[[2]]) & !pop_label%in%c("CD580L0_herm","CD680L0_herm","CD780L0_herm","OF5L0_herm"))
data_populations_herm$pop_label=as.factor(as.character(data_populations_herm$pop_label))

data_populations_male <- subset(data_populations_male,!substring(population,1,2)=="CM" & population2=="pops" & is.na(tstrsplit(data_populations_male$pop_label,"noML")[[2]]) & !pop_label%in%c("CD580L0_male","CD680L0_male","CD780L0_male","OF5L0_male"))
data_populations_male$pop_label=as.factor(as.character(data_populations_male$pop_label))


# Attribute the number of generations for each population
gen_df <- data.frame(pop_label=sort(unique(as.character(data_populations_herm$pop_label))),gen=c(0,100,10,30,60,100,10,30,60,100,10,30,60,100,10,40,30,60,100,10,40,60,100,10,40,60,100,10,36,50,5,68,100,10,36,50,5,68,100,10,36,50,5,68,32,66,32,66,32))
gen_df[27:49,]$gen=gen_df[27:49,]$gen+140

data_populations_herm <- merge(data_populations_herm, gen_df)

gen_df_male=gen_df
gen_df_male$pop_label= (paste0(tstrsplit(gen_df_male$pop_label,"_")[[1]],"_male"))

data_populations_male <- merge(data_populations_male, gen_df_male)
dim(data_populations_male)

data_populations_herm$anc <- substring(data_populations_herm$population,1,1)
data_populations_male$anc <- substring(data_populations_male$population,1,1)

data_populations_herm$gen[(data_populations_herm$anc=='C')]=data_populations_herm$gen[(data_populations_herm$anc=='C')]-140
data_populations_male$gen[(data_populations_male$anc=='C')]=data_populations_male$gen[(data_populations_male$anc=='C')]-140


par(mfrow=c(2,3))
results_models_all <- NULL
for(k in 1:2){
  
#  k==1 # herm
#  k==2 # male
  
  if(k==1) data_populations=subset(data_populations_herm,population2=="pops")
  if(k==2) data_populations=subset(data_populations_male,population2=="pops")
  
  
  
  for(i in 1:6){

    
    

    mod_temp <- lmer(data_populations[, vect_P_traits[i]]~   gen + (anc|gen) + (1|data_group_name)  ,data=data_populations)
    mod_temp2 <- lmer(data_populations[, vect_P_traits[i]]~ 1 + (anc|gen) + (1|data_group_name) ,data=data_populations)
    
    results_models_all <- rbind(results_models_all,data.frame(trait= substring(vect_P_traits[i],1,3),Chisq=anova(mod_temp,mod_temp2)$Chisq[2],P=anova(mod_temp,mod_temp2)$Pr[2]))
    
    qqnorm(resid(mod_temp))
    qqline(resid(mod_temp))
    
    
  }
}

results_models_all$Sex <- rep(c("H","M"),each=6)
results_models_all$P_adjust <- p.adjust(results_models_all$P,method="BH")
print(subset(results_models_all,P<0.05))






