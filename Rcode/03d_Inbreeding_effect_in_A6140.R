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

##  Effect of Inbreeding in A6140
#First we need to get the RILs and populations data

data_populations_herm <- read.table("Output_files/txt/data_populations_herm_original_covariates.txt",sep='\t',h=TRUE)
data_populations_A6140= subset(data_populations_herm,pop_label=="A6140L0_herm")
dim(data_populations_A6140) # 6 26

final_merged =read.table("data/Final_merged_data_NGM.txt",h=TRUE,sep="\t")
final_merged$population = as.factor((tstrsplit(final_merged$pop_label, split = "L", fixed = TRUE)[[1]]))
final_merged=subset(final_merged,population=="A6140" & location_label=='Lisbon')
dim(final_merged) # 370 22

final_merged $year <- substring(final_merged $date_str,1,4)
final_merged$is_2012 <- !(final_merged$year=="2012")
final_merged$population <- as.factor(as.character(final_merged$population))

final_merged=merge(read.table("data/State_freq_NGM.txt",sep='\t',h=TRUE), final_merged)
dim(final_merged) # 370 27

data_populations_A6140$year <- substring(data_populations_A6140$date_str,1,4)
data_populations_A6140$is_2012= (data_populations_A6140$year=="2012")
final_merged$population2="RILS"
dim(final_merged);dim(data_populations_A6140)

final_inbreeding <- rbind(final_merged, data_populations_A6140)
final_inbreeding$is_pop <- final_inbreeding$population2=="pops"

final_RILs <- final_merged
rm(final_merged);gc()

# Normalize environmental covariates
for(i in c("temperature","rel_humidity","logD")){
  final_inbreeding[,i] <- (final_inbreeding[,i]-mean(final_inbreeding[,i]))/sd(final_inbreeding[,i])
}

# Then we need to test the effect of inbreeding.
# We need to fit the mean RILs estimate per population (not as random effect anymore)
table(subset(final_inbreeding,is_pop)$population)
data_to_test <- subset(final_inbreeding,population == 'A6140')
data_to_test$pop_ispop <- paste(data_to_test$is_pop,data_to_test$population)
table(data_to_test$is_pop)

#FALSE  TRUE 
#370    6

summary_results <- NULL
pdf('plots/QQplot_Inbreeding_lmer_only_A6140.pdf')
par(mfrow=c(2,3))
for(i in 1:6){
  
  final_model <- lmer(data_to_test[, vect_P_traits[i]] ~  1  + (temperature+rel_humidity+logD)^3 + is_2012 + is_pop +(1|date_str)+(1|pop_label),data= data_to_test)
  final_model2 <- lmer(data_to_test[, vect_P_traits[i]] ~ 1  + (temperature+rel_humidity+logD)^3 + is_2012 +(1|date_str)+(1|pop_label),data= data_to_test)
  
  
  summary_results <- rbind(summary_results,
                           c(anova(final_model,final_model2)$Chisq[2],anova(final_model,final_model2)$Pr[2]))
  qqnorm(resid(final_model))
  qqline(resid(final_model))
}
dev.off()

summary_results <- as.data.frame(summary_results)
names(summary_results) <- c("Chisq","P")
summary_results$traits <- vect_P_traits
summary_results <- summary_results[,c(3,1,2)]
summary_results$P.adjusted=p.adjust(summary_results$P)
print(summary_results)

write.table(summary_results,'Output_files/txt/Inbreeding_table_only_A6140.txt',row.names=FALSE,quote=FALSE,sep='\t')

