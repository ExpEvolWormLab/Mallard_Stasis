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

##  Effect of male frequency on the outbred populations

data_populations_herm <- read.table("Output_files/txt/data_populations_herm_original_covariates.txt",sep='\t',h=TRUE)

male_freq <- read.table("data/population_rep_male_freqs.txt",h=TRUE)
names(male_freq)[1] <- 'population'
head(male_freq)

# Add to data_populations, remove noM populations (no mf estimate)
data_populations_herm  <- merge(male_freq[,c(1,2,5)] ,subset(data_populations_herm,is.na(tstrsplit(data_populations_herm$pop_label,"noM",fixed=TRUE)[[2]])))
dim(data_populations_herm) # 102 27


retained_model <- list()

yl1=c(-3,-4,-3,-6,-2,-4)
yl2=c(0,-1,1,-1,1,0)
### Analyze the effect of male frequency on the transition rate values in the populations

pdf("plots/Male_effect_on_Trait.pdf",w=10)
par(mfrow=c(2,3))
for(i in 1:6){
  
  vect_phenotypes= data_populations_herm[, vect_P_traits[i]]
  # Main model
  temp_mod1 <- lmer(vect_phenotypes ~ (1|population)+(1|date_str)+mf,data= data_populations_herm)
  
  # No male freq. effect
  temp_mod2 <- lmer(vect_phenotypes ~ (1|population)+(1|date_str),data= data_populations_herm)

  plot(vect_phenotypes ~ data_populations_herm$mf,pch=16,type="n",las=1,bty="n",ylab="Log transition rate (Hz)",xlab="Male frequency",main=c("SF","SB","FS","FB","BS","BF")[i],
       ylim=c(yl1[i],yl2[i]))
  points(vect_phenotypes ~ data_populations_herm$mf,pch=16,type="p")
  
  #is mf significant ?
  if(anova(temp_mod1,temp_mod2,test="F")$Pr[2] < 0.05){
    
    retained_model[[i]] <- temp_mod1
    print(r.squaredGLMM(temp_mod1))
    lines(I(c(1:10)/10),predict(temp_mod1,data.frame(mf=I(c(1:10)/10)),re.form=NA),col="red")
    text(0.2,c(-2.5,-4,-2,-2,-1,1)[i],paste("p=",round(anova(temp_mod1,temp_mod2,test="F")$Pr[2],digits=3)))
    
    text(0.2,c(-2.5,-4,-2,-2,-1,1)[i],paste("p=",round(anova(temp_mod1,temp_mod2,test="F")$Pr[2],digits=3)))
    text(0.2,(c(-2.5,-4,-2,-2,-1,1)*1.25)[i],paste("r2=",round(r.squaredGLMM(temp_mod1)[1],digits=3)))
    
  }else{
    #no MF effect
    retained_model[[i]] <- temp_mod2
    print(r.squaredGLMM(temp_mod2))
    
    lines(I(c(1:10)/10),predict(temp_mod2,data.frame(is_pop=FALSE,mf=I(c(1:10)/10)),re.form=NA),col="black")
    
  }
}
dev.off()

