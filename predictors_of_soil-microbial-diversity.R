library(MASS)
library(lme4) 
library(car)
library(glmulti)
library(sjlabelled)

setwd("/ciencia/data/soil/")

############# determine best predictors of diversity ##############
taxa_level5 <- read.csv("taxa-barplots_level-5_NTCrm_lowseqrm.csv", header=TRUE)
dim(taxa_level5) #150 591... taxon columns: C-UU (total in UV)

taxaonly_level5 <- taxa_level5[,3:567]

propEnt <- taxa_level5[,485]/taxa_level5$total_sequences
randomtaxonprop <- taxa_level5[,312]/taxa_level5[,568]

shannonall <- glm(taxa_level5$shannondepth.6667_iter.2 ~ taxa_level5$category + propEnt + taxa_level5$tot_carb * taxa_level5$tot_nitro + taxa_level5$carb_nitro_ratio + taxa_level5$ph * taxa_level5$water_content_soil)
summary(shannonall) 
Anova(shannonall)
Anova(shannonall)[,3]
shannonall$aic


plot(propEnt, taxa_level5$shannondepth.6667_iter.2) #is this driven by that one outlier?
abline(lm(taxa_level5$shannondepth.6667_iter.2 ~ propEnt))
summary(lm(taxa_level5$shannondepth.6667_iter.2 ~ propEnt)) #p=0.00148

#try removing outlier
diversitytest <- read.csv("diversity_enterobacteriaceae_level5.csv", header=TRUE)
newpropEnt <- diversitytest[,1]/diversitytest$total_sequences
summary(lm(diversitytest$shannondepth.6667_iter.2 ~ newpropEnt)) #p=0.493
plot(newpropEnt, diversitytest$shannondepth.6667_iter.2)
abline(lm(diversitytest$shannondepth.6667_iter.2 ~ newpropEnt))

#are other taxa equally good predictors?
shannrand_Anovaresults <- data.frame(Factor=character(), pvalue=numeric(), stringsAsFactors=FALSE)
shannrand_AIC <- matrix(NA, nrow=567, ncol=1)

# index of current taxon output
cur_taxon=1;

for (i in 3:567) {
  prop_randtax <- taxa_level5[,i]/taxa_level5[,568]
  shannonrand <- glm(taxa_level5$shannondepth.6667_iter.2 ~ taxa_level5$category + prop_randtax + taxa_level5$tot_carb * taxa_level5$tot_nitro + taxa_level5$carb_nitro_ratio + taxa_level5$ph * taxa_level5$water_content_soil)
  shannrand_Anovaresults  <- rbind(shannrand_Anovaresults, Anova(shannonrand)[,3]) 
  shannrand_AIC[i,] <- shannonrand$aic
  
  # increment the counter
  cur_taxon=cur_taxon+1;
}

write.table(shannrand_Anovaresults, file="shannon_randtaxaANOVA_pvalues_level5.txt", col.names=F, append=T)  
write.table(shannrand_AIC, file="shannon_randtaxaANOVA_AIC_level5.txt", col.names=F, append=T) #i should cbind this to the previous


#what are the means across groups?
aggregate(taxa_level5$shannondepth.6667_iter.2, list(taxa_level5$category), mean, na.rm=TRUE)
#     adjacent 8.577093
#        bones 8.214101
# bones + scat 8.069607
#         dead 8.482763
#     entrance 8.414987
#       plague 8.568545
#         scat 7.903954


faithpall <- glm(taxa_level5$faithdepth.10000_iter.2 ~ taxa_level5$category + propEnt + taxa_level5$ph * taxa_level5$water_content_soil + taxa_level5$tot_carb * taxa_level5$tot_nitro + taxa_level5$carb_nitro_ratio)
summary(faithpall) 
Anova(faithpall)

plot(propEnt, taxa_level5$faithdepth.10000_iter.2)
abline(lm(taxa_level5$faithdepth.10000_iter.2 ~ propEnt))
summary(lm(taxa_level5$faithdepth.10000_iter.2 ~ propEnt)) #p=0.00261

summary(lm(diversitytest$faithdepth.10000_iter.2 ~ newpropEnt)) #p=0.131
plot(newpropEnt, diversitytest$faithdepth.10000_iter.2)
abline(lm(diversitytest$faithdepth.10000_iter.2 ~ newpropEnt))


#are other taxa better predictors?
faithrand_Anovaresults <- data.frame(Factor=character(), pvalue=numeric(), stringsAsFactors=FALSE)
faithrand_AIC <- matrix(NA, nrow=567, ncol=1)

# index of current taxon output
cur_taxon=1;

for (i in 3:567) {
  prop_randtax <- taxa_level5[,i]/taxa_level5[,568]
  faithrand <- glm(taxa_level5$faithdepth.10000_iter.2 ~ taxa_level5$category + prop_randtax + taxa_level5$tot_carb * taxa_level5$tot_nitro + taxa_level5$carb_nitro_ratio + taxa_level5$ph * taxa_level5$water_content_soil)
  faithrand_Anovaresults  <- rbind(faithrand_Anovaresults, Anova(faithrand)[,3]) 
  faithrand_AIC[i,] <- faithrand$aic
  
  # increment the counter
  cur_taxon=cur_taxon+1;
}

write.table(faithrand_Anovaresults, file="faith_randtaxaANOVA_pvalues_level5.txt", col.names=F, append=T)  
write.table(faithrand_AIC, file="faith_randtaxaANOVA_AIC_level5.txt", col.names=F, append=T) 


#means by category
aggregate(taxa_level5$faithdepth.10000_iter.2, list(taxa_level5$category), mean, na.rm=TRUE)
#     adjacent 87.46738
#        bones 80.59463
# bones + scat 77.07492
#         dead 86.34845
#     entrance 82.25544
#       plague 85.56236
#         scat 74.57382


# plot both measures of diversity together
par(mar = c(6, 4, 0.5, 0.5) + 0.1)
par(mfrow=c(1,2)) 
boxplot(taxa_level5$shannondepth.6667_iter.2 ~ taxa_level5$category, na.rm=TRUE, col=c("gray", "saddlebrown", "wheat2", "brown3", "gray14", "gray32", "yellowgreen"), ylab="Shannon diversity index", varwidth=TRUE, las=2,  xlab = "") 
boxplot(taxa_level5$faithdepth.10000_iter.2 ~ taxa_level5$category, na.rm=TRUE, col=c("gray", "saddlebrown", "wheat2", "brown3", "gray14", "gray32", "yellowgreen"), ylab="Faith phylogenetic diversity", varwidth=TRUE, las=2,  xlab = "") 

