covariates <- read.csv("covariates.csv", header = TRUE)
pdf("letssee.pdf")
for(a in 1:length(Protein.unique)){
  proID <- box %>% filter(Gene == Protein.unique[a]) %>% select(ID)  
  proID <- as.character(proID$ID)
  
  IDbox <- Master[c(as.character(SLnames$Variable[match(Protein.unique[a], SLnames$Gene)]), proID)]
  
  #Converts percentage dosages to 0, 1, 2 factors
  for(j in 2:ncol(IDbox)){
    test <- as.numeric(as.character(IDbox[,j]))
    
    for(i in 1:length(test)){
      
      if(test[i]>1.5){
        test[i] <- 2
      }
      if(test[i]>0.5 && test[i]<=1.5){
        test[i] <- 1
      }
      if(test[i]<=0.5){
        test[i] <- 0
      }
    }
    
    test <- as.factor(test)
    IDbox[,j] <- test
  }
  
  for(k in 2:ncol(IDbox)){
    IDtwo <- IDbox %>% select(c(colnames(IDbox)[1],colnames(IDbox)[k]))
    colnames(IDtwo) <- c("protein", "Alleles")
    IDthree <- IDtwo %>% mutate(age = agesgenome$age_at_LP)
    IDthree$age <- factor(IDthree$age, levels = c("young", "mid-age", "old"))
    print(ggplot(IDthree, aes(x=age, y = protein, fill=Alleles)) + geom_boxplot(alpha=0.3) + geom_point(size = .8, position = position_dodge(width=0.75)) + labs(y=paste(Protein.unique[a], "(RFUs)")) + ggtitle(paste(colnames(IDbox)[k])))
  }
}
dev.off()
# for(c in 1:nrow(covariates)){
# 
#   if(covariates$age_at_LP[c]>=60){
#     covariates$age_at_LP[c] <- "old"
#   }  
#   if(covariates$age_at_LP[c]<=40){
#     covariates$age_at_LP[c] <- "young"
#   }
#   if(covariates$age_at_LP[c]>40 & covariates$age_at_LP[c]<60){
#     covariates$age_at_LP[c] <- "mid-age"
#   }
# }
ages <- covariates %>% select(patient, age_at_LP)
ages$patient <- as.character(ages$patient)

agesgenome <- merge(ages, Master, by="patient")
agesgenome <- agesgenome[,-c(3)]

