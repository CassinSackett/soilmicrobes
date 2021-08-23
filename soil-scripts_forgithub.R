soil_data <- read.csv("/ciencia/data/soil/soil_data.csv") 

category <- as.factor(category)
levels(category) 
# get rid of the blanks as a level 
sex <- Sex[drop=TRUE]

# boxplots for each measure by location in burrow ("category"; 10 levels). 
boxplot(pH ~ category, col=c("gray", "saddlebrown", "wheat2", "brown3", "gray14", "gray32", "yellowgreen"), ylab="pH", varwidth=TRUE, outline=FALSE) 

#test whether each measure varies significantly by burrow region
tuk <- aov(C.N ~ category)
posthoc <- TukeyHSD(tuk, 'category') #defaults to all


########################################################################
############# determine best predictors of diversity ##############
taxa_level5 <- read.csv("taxa-barplots_level-5_NTCrm_lowseqrm.csv", header=TRUE)
dim(taxa_level5) #150 591... taxon columns: C-UU (total in UV)
colnames(taxa_level5[2])

#subset data to exclude nutrient measures, sample metadata etc. 
taxaonly_level5 <- taxa_level5[,3:567]

#procedures were done on Faith's Phylogenetic Diversity index and can be repeated for other measures of diversity such as the shannon index
hist(taxa_level5$faithdepth.10000_iter.2) 

#what are the means/variation across categories (burrow regions?
aggregate(taxa_level5$faithdepth.10000_iter.2, list(taxa_level5$category), mean, na.rm=TRUE)
boxplot(taxa_level5$faithdepth.10000_iter.2 ~ taxa_level5$category, na.rm=TRUE, col=c("gray", "saddlebrown", "wheat2", "brown3", "gray14", "gray32", "yellowgreen"), ylab="Faith phylogenetic diversity", varwidth=TRUE, las=2,  xlab = "") 

#proportion of total sequences that are Enterobacteraceae
propEnt <- taxa_level5[,485]/taxa_level5$total_sequences

# are the same patterns observed with other taxa?
randomtaxonprop <- taxa_level5[,312]/taxa_level5[,568]
Planococcaceae <- taxa_level5[,277]/taxa_level5[,568]
Micrococcaceae <- taxa_level5[,88]/taxa_level5[,568]

# model ecological diversity among samples
faithall <- glm(taxa_level5$faithdepth.10000_iter.2 ~ taxa_level5$category + propEnt + taxa_level5$tot_carb * taxa_level5$tot_nitro + taxa_level5$carb_nitro_ratio + taxa_level5$ph * taxa_level5$water_content_soil)

summary(faithall) 
faithall$aic

plot(propEnt, taxa_level5$faithdepth.10000_iter.2)  

#are other taxa equally good predictors? Loop through them all and repeat the glm with a different predictor taxon each time.

#create a dataframe to store the p values for each variable
faithrand_Anovaresults <- data.frame(Factor=character(), pvalue=numeric(), stringsAsFactors=FALSE)
#create a matrix for storing AIC values for each taxon's model
faithrand_AIC <- matrix(NA, nrow=567, ncol=1)

# index of current taxon output
cur_taxon=1;

#loop through taxa
for (i in 3:567) {
  prop_randtax <- taxa_level5[,i]/taxa_level5[,568]
  #perform the glm
  faithrand <- glm(taxa_level5$faithdepth.10000_iter.2 ~ taxa_level5$category + prop_randtax + taxa_level5$tot_carb * taxa_level5$tot_nitro + taxa_level5$carb_nitro_ratio + taxa_level5$ph * taxa_level5$water_content_soil)
  #add p values to the dataframe
  faithrand_Anovaresults  <- rbind(faithrand_Anovaresults, Anova(faithrand)[,3]) 
  #add AIC values to the matrix
  faithrand_AIC[i,] <- faithrand$aic
  
  # increment the counter
  cur_taxon=cur_taxon+1;
}

#save dataframe and matrix as files
write.table(faithrand_Anovaresults, file="faith_randtaxaANOVA_pvalues_level5.txt", col.names=F, append=T)  
write.table(faithrand_AIC, file="faith_randtaxaANOVA_AIC_level5.txt", col.names=F, append=T) #i should really cbind this to the previous

# visualize the relationship between relative abundance of particular taxa and phylogenetic diversity
colnames(taxa_level5[27])
WD2101 <- taxa_level5[,364]/taxa_level5[,568]
iii115 <- taxa_level5[,27]/taxa_level5[,568]

plot(WD2101, taxa_level5$faithdepth.10000_iter.2)


########################################################################
############ Does Enterobacteriacaea differ among regions? ############

#verify I have the right columns
colnames(taxa_level5[485]) #Enterobacteriaceae
colnames(taxa_level5[591]) #category

#convert total sequences to relative abundance
propentero <- taxa_level5[,485]/taxa_level5$total_sequences
pctentero <- propentero*100
#where is the proportion of Enterobacteriaceae the highest?
max(pctentero)  

#model dependence of Enterobacteriaceae relative abundance on burrow location
Entero <- glm(taxa_level5[,485]/taxa_level5$total_sequences ~ taxa_level5$category) 
#check p values for each location
summary(Entero)$coefficients[,4] 
#check significance of whole model
Anova(Entero) 

#model dependence of Enterobacteriaceae relative abundance on burrow location and soil chemical properties
Enteroall <- glm(taxa_level5[,485]/taxa_level5$total_sequences ~ taxa_level5$category + taxa_level5$ph + taxa_level5$water_content_soil + taxa_level5$tot_carb * taxa_level5$tot_nitro + taxa_level5$carb_nitro_ratio)
#check significance 
summary(Enteroall) 
Anova(Enteroall) 
# p values
Anova(Enteroall)[,3]

boxplot(taxa_level5[,485]/taxa_level5$total_sequences ~ taxa_level5$category, na.rm=TRUE, col=c("gray", "saddlebrown", "wheat2", "brown3", "gray14", "gray32", "yellowgreen"), ylab="proportion Enterobacteriaceae", varwidth=TRUE,  xlab = "")

# calculate mean proportions by category 
aggregate(taxa_level5[,485]/taxa_level5$total_sequences, list(taxa_level5$category), mean)


####################################################################
#### do other taxa abundances vary by category (burrow region)? ####

#set up a dataframe for storing results
category5_Anovaresults <- data.frame(Factor=character(), pvalue=numeric(), stringsAsFactors=FALSE)

# index of current taxon output.
cur_taxon=1;
# loop over all taxa
for (k in 3:567) #columns 3 to 567 have abundance data
{
  # check predictors of proportion total 
  burr_region5 <- glm(taxa_level5[,k]/taxa_level5$total_sequences ~ taxa_level5$category + taxa_level5$ph + taxa_level5$water_content_soil + taxa_level5$tot_carb * taxa_level5$tot_nitro + taxa_level5$carb_nitro_ratio)
  
  # write p values to the dataframe 
  category5_Anovaresults  <- rbind(category5_Anovaresults, Anova(burr_region5)[,3])
  
  # increment the counter
  cur_taxon=cur_taxon+1;
}

#save the dataframe to a file
write.table(category5_Anovaresults, file="glmANOVA_pvalues_taxalevel5_byregion.txt", col.names=F, append=T)  

# check the direction of effect for highly significant taxa
par(mfrow=c(1,1))
aggregate(taxa_level6[,362]/taxa_level6$total_sequences, list(taxa_level6$category), mean) 
boxplot(taxa_level6[,362]/taxa_level6$total_sequences ~ taxa_level6$category, col=c("gray", "saddlebrown", "wheat2", "brown3", "gray14", "gray32", "yellowgreen"), ylab="proportion unknown AKIW781 (Chloroflexi)", xlab="Burrow Regions", varwidth=TRUE, outline=FALSE)


######################################################################
#### test for association with presence of Yersinia and certain lineages, or presence of Yersinia and diversity #####

# Goal: Loop through all taxa and do a simple test to see if Yersinia presence/absence predicted numbers of that taxon (controlled by "total_sequences" to account for variation among individuals in sequencing depth). 

#get the name of a particular taxon (the header)
colnames(taxa_level5[62]) 

#make the vector a factor
yersiniaYN <- as.factor(taxa_level5$contains_Yersinia)

#make a table to store p-values. Need only 1 col for kruskal-wallis. 567 rows for level 5
diversity_p <- matrix(NA, nrow=567, ncol=1)

# index of current taxon outputed.
cur_taxon=1;
# loop over all taxa
for (k in 3:567) #columns 3 to 567 have abundance data
{
  # check if relationship between proportion total & yersinia pres
  yers_pres <- kruskal.test(taxa_level5[,k]/taxa_level5[,568], yersiniaYN)
  
  # check if proportion total predicts yersinia presence
  yers <- kruskal.test(yersiniaYN ~ taxa_level5[,k]/taxa_level5[,568])
  
  # write p values (coefficients[,4]) to a table
  diversity_p[k, ] <- yers_pres$p.value
  
  # increment the counter
  cur_taxon=cur_taxon+1;
}

write.table(diversity_p, file="wilcox_pvalues_taxalevel5_NTCrm.csv", col.names=F, append=T)

Yersiniapresence <- factor(taxa_level6[,2], levels=c(0,1), labels=c("Absent", "Present"))

# get rid of scientific notation for plot
options(scipen=10)
# plot top results level 5 
boxplot(taxa_level5[,436]/taxa_level5[,568] ~ Yersiniapresence, col="honeydew3", xlab="", ylab="Brucellaceae") 


################# diversity vs Yersinia ###################
#I am adding faithpd at depth 10000 to the taxa file. column name is faithpd_depth10000_iter1 

aggregate(taxa_level6[,569], list(taxa_level5[,2]), mean, na.rm=TRUE)

# check if relationship between faith phylogenetic diversity & yersinia presence
faithpd <- kruskal.test(taxa_level5[,569], taxa_level5[,2])

### ... is yersinia presence predicted by diversity? No ###
faithpd <- kruskal.test(taxa_level5[,2], taxa_level6[,569]) #p=0.48


##############################################################
### Are other taxa significant predictors of Yersinia? ###
tax_yers <- matrix(NA, nrow=567, ncol=1)

# index of current taxon output
cur_taxon=1;

#loop through all taxa
for (i in 3:567) {
  #calculate relative abundance - column 568 is total sequences per sample
  prop_randtax <- taxa_level5[,i]/taxa_level5[,568]
  #test whether relative abundance predicts Yersinia presence
  Yersiniapredictors <- glm(taxa_level5$contains_Yersinia ~ prop_randtax)
  #store p values
  tax_yers[i,] <- Anova(Yersiniapredictors)[,3]
  
  # increment the counter
  cur_taxon=cur_taxon+1;
}

#save the file of p values
write.table(tax_yers, file="taxa_predicting_yersinia-presence_level5.txt", col.names=F, append=T) 


##############################################################
### chi square to determine whether some taxa are more/less abundant across burrow regions than expected ###

##### the below files are from hand-curated dominant taxa lists ##

library(vcd) #for mosaic plots of significance of chi-square
# this mosaic plot with colored cases shows where the observed frequencies deviates from the expected frequencies if the variables were independent. The red cases means that the observed frequencies are smaller than the expected frequencies, whereas the blue cases means that the observed frequencies are larger than the expected frequencies.

dominant_taxa <- read.csv("dominant_taxa_archaearm.csv", header=TRUE)
# build a contingency table with table() function - repeat for other taxonomic levels
phylumtable <- table(dominant_taxa$dominantphylum, dominant_taxa$category)
#perform the chi square test
phylum_test <- chisq.test(phylumtable) 

#create mosaic plot to visualize results
mosaic( ~ dominant_taxa$category + dominant_taxa$dominantphylum, direction=c("v", "h"), data=dominant_taxa, shade=TRUE, na.rm=TRUE) 

# check the number of taxa
phlm <- factor(dominant_taxa$dominantphylum)
levels(phlm)


#########################################################################
####### make individual taxa barplots using means of 2 replicates #######

# will need some colors
library("viridis")

level5_barplot <- read.csv("replicate_mean_proportions.csv", header=TRUE)

## sort the dataset by category ##
level5_barplot <- level5_barplot[order(level5_barplot$category),]

# would like to find most abundant and cbind it to the dataframe
level5_wdominant <- cbind(level5_barplot, max(level5_barplot[,7:570]))

#needs to be a matrix for barplot function
level5_matrix <- as.matrix(t(level5_barplot[7:570]))

#adjust margins bottom, left, top, right
par(mar=c(7,4,1,1))
sample_barplot <- barplot(level5_matrix, beside=FALSE, space=0.1, names.arg=level5_barplot$sample, col=viridis(50), border=NA, main="", xlab="", ylab="", las=2)

######## plot per burrow region ########
barplot5category <- read.csv("replicate_mean_proportions_bycategory.csv", header=TRUE)
barplot5_matrix <- as.matrix(t(barplot5category[7:570]))
par(mar=c(7,4,1,1))
category_barplot <- barplot(barplot5_matrix, beside=FALSE, space=0.1, names.arg=barplot5category$category, col=viridis(100), border=NA, main="", xlab="", ylab="Proportion of total sequences", las=2)




