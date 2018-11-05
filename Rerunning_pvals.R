new.genotypes <- read.csv("new genotype data.csv", header = TRUE)
#new.genotypes <- read.table("chrAll transposed.txt", header = TRUE)
new.eqtls <- read.csv("missing protein eqtls.csv", header = TRUE)
#new.eqtls <- read.csv("old eqtls redo.csv", header = TRUE)
msdss.score <- read.table("Scores.txt", header = TRUE)
new.proteins <- read.csv("updated proteins.csv", header = TRUE)
#new.proteins <- read.csv("old protein levels.csv", header = TRUE)
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

p.df <- data.frame(matrix(vector(), 0, 3))
colnames(p.df) <- c("rs", "gene", "p-val")

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
    if(ncol(testingdf) > 0){
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
          
          testingdf[,j] <- test
        }
        
        protein.interested <- as.character(slnames[match(Protein.unique[a], slnames$Gene),1])
        
        protein.int.level <- new.master[[protein.interested]]
        
        protein.geno <- cbind(protein.int.level, testingdf)
        
        for(k in 2:length(protein.geno)){
          p.df <- p.df %>% bind_rows(c("rs" = colnames(protein.geno)[k], "gene" = Protein.unique[a], "p-val" = cor.test(protein.geno[[k]], protein.geno[[1]])$p.value))
        }
    }
}

write.csv(p.df, "new eqtl p-values.csv")
