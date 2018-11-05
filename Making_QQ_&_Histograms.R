Genenumber <- read.csv("Genenumber.csv", header=TRUE)
SLpvalues <- Genenumber %>% select(ID) %>% mutate(Pvalues=NA)

SLunique <- as.character(unique(Genenumber$Variable))

for(i in 1:length(SLunique)){
  rsSL <- Genenumber %>% filter(Variable == SLunique[i]) %>% select(ID)

  rsSL <- as.character(rsSL$ID)

  genomeSL <- Master[c(SLunique[i], rsSL)]

  for(j in 1:(ncol(genomeSL)-1)){
    genome <- as.numeric(as.character(genomeSL[,j+1]))
    model <- lm(genomeSL[,1]~genome)
    SLpvalues <- SLpvalues %>% mutate(Pvalues = replace(Pvalues, ID==eval(colnames(genomeSL[j+1])), summary(model)$coefficients[2,4]))
  }
}

finalpvals <- merge(SLpvalues, Genenumber, by="ID") %>% select(ID, Pvalues, Gene)

Geneunique <- as.character(unique(finalpvals$Gene))

pdf("qqplots1.pdf")
par(mfrow=c(2,1))
for(i in 1:length(Geneunique)){
  rsps <- finalpvals %>% filter(Gene == Geneunique[i]) %>% select(Pvalues)
  rsps <- rsps$Pvalues
  
  psort <- sort(rsps)
  neg <- -log10(psort)
  expected <- c(1:length(neg))
  lexp <- -(log10(expected/(length(expected)+1)))
  
  par(family="serif")
  plot(neg~lexp, main=paste(Geneunique[i],"QQ Plot"), xlab = "Expected [-log10(P)]", ylab = "Observed [-log10(P)]", asp=1)
  abline(a=0, b=1, col = "blue")
}
dev.off()


pdf("hist.pdf")
par(mfrow=c(2,1))
for(i in 1:length(Geneunique)){
  rsps <- finalpvals %>% filter(Gene == Geneunique[i]) %>% select(Pvalues)
  rsps <- rsps$Pvalues

  par(family="serif")
  hist(rsps, main=paste(Geneunique[i]), xlab = "P-Values", col = "lightblue")
}
dev.off()

