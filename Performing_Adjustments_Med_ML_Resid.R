#Importing datasets
adjustingSNPs <- read.csv("significant SNPs for adjustments.csv", header = TRUE)
Master <- read.table("newmaster.txt", header=TRUE)
slnames <- read.csv("all protein slname.csv", header = TRUE)

colnames(adjustingSNPs)[1] <- "rs"
colnames(slnames)[1] <- "Variable"

#Identifying unique proteins
Protein.unique <- as.character(unique(adjustingSNPs$gene))

#Creating empty data frames for adjusted protein levels
AdjustedProteins1 = data.frame(matrix(vector(), 194, 5))
colnames(AdjustedProteins1) <- Protein.unique
AdjustedProteins1log = data.frame(matrix(vector(), 194, 5))
colnames(AdjustedProteins1log) <- Protein.unique
AdjustedProteins2 = data.frame(matrix(vector(), 194, 5))
colnames(AdjustedProteins2) <- Protein.unique

#Making a loop to adjust all protein levels (Current loop only handles three genotypes)
for(a in 1:length(Protein.unique)){
    #Extracting protein levels and genotypes from MASTER dataset based on rs number
    proID <- adjustingSNPs %>% filter(gene == Protein.unique[a]) %>% select(rs)
    proID <- as.character(proID$rs)
    
    #Creating data frame with protein levels and genotypes for a single protein
    uniquegeno <- Master[c(as.character(slnames$Variable[match(Protein.unique[a], slnames$Gene)]), proID)]
    
    #Converts percentage dosages to 0, 1, 2 factors
    for(j in 2:ncol(uniquegeno)){
      test <- as.numeric(as.character(uniquegeno[,j]))
      
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
      uniquegeno[,j] <- test
    }
    
    #Creating smaller data frame with only protein level and a single SNP
    adjustdata <- uniquegeno %>% select(c(colnames(uniquegeno)[1], colnames(uniquegeno)[2]))
   
    #Splitting protein levels based on genotype
    adjustdata.split<-split(adjustdata,adjustdata[[2]])
    
    
    if(length(adjustdata.split)==2){
      #Assigning separate data frames for each genotype (0, 1, 2)
      aa <- adjustdata.split[[1]]
      Aa <- adjustdata.split[[2]]
      
      
      #Creating new dataframes that will store new adjustments
      a.adjusted <- aa
      b.adjusted <- Aa
      
      
      #Performing adjustments by subtrating or adding the absolute difference from the medians
      a.adjusted[[1]] <- a.adjusted[[1]] + (median(b.adjusted[[1]])-median(a.adjusted[[1]]))
      
      #Combining adjusted dataframes together
      rsall <- rbind(a.adjusted, b.adjusted) 
      
      #Reordering adjusted dataframe into original order
      dat <- rsall[sample(nrow(rsall), 194), ]
      rn <- rownames(dat)
      rsall <- dat[order(as.numeric(rn)), ]
      AdjustedProteins1[[a]] <- rsall[[1]]
      
      #Same thing but adjusting the medians with logs
      a.adjusted.log <- aa
      b.adjusted.log <- Aa
      a.adjusted.log[[1]] <- log(a.adjusted.log[[1]])
      b.adjusted.log[[1]] <- log(b.adjusted.log[[1]])
      a.adjusted.log[[1]] <- a.adjusted.log[[1]] + (median(b.adjusted.log[[1]])-median(a.adjusted.log[[1]]))
      rsall.log <- rbind(a.adjusted.log, b.adjusted.log) 
      
      dat <- rsall.log[sample(nrow(rsall.log), 194), ]
      rn <- rownames(dat)
      rsall.log <- dat[order(as.numeric(rn)), ]
      rsall.log[[1]] <- exp(rsall.log[[1]])
      AdjustedProteins1log[[a]] <- rsall.log[[1]]
    }
    
    if(length(adjustdata.split) == 3){
      #Assigning separate data frames for each genotype (0, 1, 2)
      aa <- adjustdata.split[[1]]
      Aa <- adjustdata.split[[2]]
      AA <- adjustdata.split[[3]]
      
      #Creating new dataframes that will store new adjustments
      a.adjusted <- aa
      b.adjusted <- Aa
      c.adjusted <- AA
      
      #Performing adjustments by subtrating or adding the absolute difference from the medians
      a.adjusted[[1]] <- a.adjusted[[1]] + (median(b.adjusted[[1]])-median(a.adjusted[[1]]))
      c.adjusted[[1]] <- c.adjusted[[1]] + (median(b.adjusted[[1]])-median(c.adjusted[[1]]))
      
      #Combining adjusted dataframes together
      rsall <- rbind(a.adjusted, b.adjusted, c.adjusted) 
      
      #Reordering adjusted dataframe into original order
      dat <- rsall[sample(nrow(rsall), 194), ]
      rn <- rownames(dat)
      rsall <- dat[order(as.numeric(rn)), ]
      AdjustedProteins1[[a]] <- rsall[[1]]
      
      #Same thing but adjusting the medians with logs
      a.adjusted.log <- aa
      b.adjusted.log <- Aa
      c.adjusted.log <- AA
      a.adjusted.log[[1]] <- log(a.adjusted.log[[1]])
      b.adjusted.log[[1]] <- log(b.adjusted.log[[1]])
      c.adjusted.log[[1]] <- log(c.adjusted.log[[1]])
      a.adjusted.log[[1]] <- a.adjusted.log[[1]] + (median(b.adjusted.log[[1]])-median(a.adjusted.log[[1]]))
      c.adjusted.log[[1]] <- c.adjusted.log[[1]] + (median(b.adjusted.log[[1]])-median(c.adjusted.log[[1]]))
      rsall.log <- rbind(a.adjusted.log, b.adjusted.log, c.adjusted.log) 
      
      dat <- rsall.log[sample(nrow(rsall.log), 194), ]
      rn <- rownames(dat)
      rsall.log <- dat[order(as.numeric(rn)), ]
      rsall.log[[1]] <- exp(rsall.log[[1]])
      AdjustedProteins1log[[a]] <- rsall.log[[1]]
      
      #New list that stores individual SNP
      n <- as.numeric(as.character(uniquegeno[[2]]))
      
      #Generates a linear model of SNP vs Protein
      mod <- lm(uniquegeno[[1]] ~ uniquegeno[[2]])
      
      #Generates list of residuals for each value
      pro.adj.outcome <- resid(mod)
      
      #Adding adjusted proteins to large data frame
      AdjustedProteins2[[a]] <- pro.adj.outcome
    }
}
