#### test for association with presence of Yersinia and certain lineages, or presence of Yersinia and diversity #####

#these files have a million columns, so looking at "head" may not be the best idea. they do each have a column "contains_Yersinia".
taxa_level5 <- read.csv("taxa-barplots_level-5_contains-Enterobacteriaceae.csv")
dim(taxa_level5) #183 637... taxon columns: C-UU (total in UV)
taxa_level6 <- read.csv("taxa-abundance_level-6_contains-yersinia.csv")
dim(taxa_level6) #183 1063... taxon columns: C-ALC (total in ALE)

#Loop through all taxa and test whether Yersinia presence/absence predicted numbers of that taxon (controlled by "total_sequences" column). 

#get the name of a particular taxon (the header)
colnames(taxa_level6[710]) 

#test first with one
Brucellaceae <- glm.nb(taxa_level6[,642]/taxa_level6[,993] ~ taxa_level6[,2])
summary(Brucellaceae)
#summary(Brucellaceae)$coefficients[,4]

#make a table to store p-values. Need only 1 col for kruskal-wallis. 992 rows for level 6, 567 for level 5
diversity_p <- matrix(NA, nrow=567, ncol=1)

# index of current taxon outputed.
cur_taxon=1;
# loop over all taxa
for (k in 3:567) #columns 3 to 992 (567) have abundance data
{
  # check if relationship between proportion total & yersinia pres
  yers_pres <- kruskal.test(taxa_level5[,k]/taxa_level5[,568], taxa_level5[,2])
  
  # write p values (coefficients[,4]) to a table
  diversity_p[k, ] <- yers_pres$p.value
  
  # increment the counter
  cur_taxon=cur_taxon+1;
}

write.table(diversity_p, file="wilcox_pvalues_taxalevel5_NTCrm.csv", col.names=F, append=T)


# check if relationship between shannon diversity & yersinia presence
shannon <- kruskal.test(taxa_level6[,995], taxa_level6[,2])

# check if relationship between faith phylogenetic diversity & yersinia presence
faithpd <- kruskal.test(taxa_level6[,994], taxa_level6[,2])


