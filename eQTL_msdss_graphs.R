new.genotypes <- read.csv("new genotype data.csv", header = TRUE)
old.genotypes <- read.table("chrAll transposed.txt", header = TRUE)
eqtls <- read.csv("eqtl msdss significant.csv", header = TRUE)
msdss.score <- read.table("Scores.txt", header = TRUE)
new.proteins <- read.csv("updated proteins.csv", header = TRUE)
old.proteins <- read.csv("old protein levels.csv", header = TRUE)
slnames <- read.csv("all protein slname.csv", header = TRUE)

colnames(new.genotypes)[1] <- "patient"
colnames(old.genotypes)[1] <- "patient"
colnames(new.proteins)[1] <- "patient"
colnames(old.proteins)[1] <- "patient"
colnames(eqtls)[1] <- "rs"
colnames(slnames)[1] <- "SLnames"
new.genotypes$patient <- as.character(new.genotypes$patient)
old.genotypes$patient <- as.character(old.genotypes$patient)
msdss.score$patient <- as.character(msdss.score$patient)
new.proteins$patient <- as.character(new.proteins$patient)
old.proteins$patient <- as.character(old.proteins$patient)
msdss.score <- msdss.score %>% select(c("patient", "msdss_last"))

new.master <- merge(msdss.score, new.genotypes, by = "patient")
new.master <- merge(new.master, old.genotypes, by = "patient")
new.master <- merge(new.master, new.proteins, by = "patient")
new.master <- merge(new.master, old.proteins, by = "patient")
new.master.names <- colnames(new.master)

Protein.unique <- as.character(unique(eqtls$Gene))

pdf("rerun msdss eqtls.pdf")
for(a in 1:length(Protein.unique)){
  proID <- eqtls %>% filter(Gene == Protein.unique[a]) %>% select(rs)  
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
  
  protein.geno <- cbind(new.master$msdss_last, testingdf)
  
  for(k in 2:length(protein.geno)){
    protein.geno.two <- as.data.frame(protein.geno[[1]])
    protein.geno.two <- cbind(protein.geno.two, protein.geno[[k]])
    colnames(protein.geno.two) <- c("protein", "Alleles")
    bp <- ggplot(protein.geno.two, aes(x=Alleles, y = protein, fill=Alleles)) + geom_boxplot(alpha=0.3) + geom_point(size = .8, position = position_dodge(width=0.75)) + labs(y="MS-DSS Score")
    print(bp + theme(legend.position = "none") + ggtitle(paste(colnames(protein.geno)[k], " (p = ", round(eqtls$p.val[match(colnames(protein.geno)[k], eqtls$rs)],6), ")", sep = "")))
  }
}
dev.off()
