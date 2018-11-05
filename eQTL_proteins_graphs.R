#new.genotypes <- read.csv("new genotype data.csv", header = TRUE)
new.genotypes <- read.table("chrAll transposed.txt", header = TRUE)
#new.eqtls <- read.csv("new eqtl significant.csv", header = TRUE)
new.eqtls <- read.csv("old eqtl significant.csv", header = TRUE)
msdss.score <- read.table("Scores.txt", header = TRUE)
#new.proteins <- read.csv("updated proteins.csv", header = TRUE)
new.proteins <- read.csv("old protein levels.csv", header = TRUE)
slnames <- read.csv("all protein slname.csv", header = TRUE)

colnames(new.genotypes)[1] <- "patient"
colnames(new.proteins)[1] <- "patient"
colnames(new.eqtls)[1] <- "rs"
colnames(slnames)[1] <- "SLnames"
new.genotypes$patient <- as.character(new.genotypes$patient)
msdss.score$patient <- as.character(msdss.score$patient)
new.proteins$patient <- as.character(new.proteins$patient)
msdss.score <- msdss.score %>% select(c("patient", "msdss_last"))

new.master <- merge(msdss.score, new.genotypes, by = "patient")
new.master <- merge(new.master, new.proteins, by = "patient")
new.master.names <- colnames(new.master)

Protein.unique <- as.character(unique(new.eqtls$Gene))

pdf("rerun old eqtls.pdf")
for(a in 1:length(Protein.unique)){
  proID <- new.eqtls %>% filter(Gene == Protein.unique[a]) %>% select(rs)  
  proID <- as.character(proID$rs)
  testingdf <- data.frame(matrix(vector(), 194, 0))
  
  for(i in 1:length(proID)){
    if(proID[i] %in% new.master.names){
      testingdf <- cbind(testingdf, new.master[[proID[i]]])
      colnames(testingdf)[length(testingdf)] <- proID[i]
    }
  }
    for(j in 1:ncol(testingdf)){
      test <- as.numeric(as.character(testingdf[,j]))
      
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
      testingdf[,j] <- test
    }
    
    protein.interested <- as.character(slnames[match(Protein.unique[a], slnames$Gene),1])
    
    protein.int.level <- new.master[[protein.interested]]
    
    protein.geno <- cbind(protein.int.level, testingdf)
    
    for(k in 2:length(protein.geno)){
      protein.geno.two <- as.data.frame(protein.geno[[1]])
      protein.geno.two <- cbind(protein.geno.two, protein.geno[[k]])
      colnames(protein.geno.two) <- c("protein", "Alleles")
      bp <- ggplot(protein.geno.two, aes(x=Alleles, y = protein, fill=Alleles)) + geom_boxplot(alpha=0.3) + geom_point(size = .8, position = position_dodge(width=0.75)) + labs(y=paste(Protein.unique[a], "(RFUs)"))
      print(bp + theme(legend.position = "none") + ggtitle(paste(colnames(protein.geno)[k], " (p = ", round(new.eqtls$p.val[match(colnames(protein.geno)[k], new.eqtls$rs)],6), ")", sep = "")))
     }
  }
dev.off()
