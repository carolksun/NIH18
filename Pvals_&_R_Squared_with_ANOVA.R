covariates <- read.csv("covariates.csv", header=TRUE)

list <-covariates$race
list<-factor(gsub("White", "white", list))
covariates$race <- list

list <-covariates$family_history_MS
list<-factor(gsub("Yes", "yes", list))
covariates$family_history_MS <- list

anova.pvalues <- data.frame(Status = double(), Sex = double(), Family = double(), Race = double())
covarprotein <- colnames(covariates)
covarprotein <- covarprotein[-c(1,2)]
rownames(anova.pvalues) <- covarprotein

rsquaredvalues <- data.frame(Status = double(), Sex = double(), Family = double(), Race = double(), LPage = double(), AgeSquared = double(), IgG = double(), CRP = double(), Serum = double())
rsquaredvalues[52,8] <- 3
rownames(rsquaredvalues) <- covarprotein

for(j in 3:6){
  for(i in 13:64){
    proteinmod <- lm(covariates[,i] ~ covariates[,j])
    fit <- anova(proteinmod)
    anova.pvalues[i-12,j-2] <- fit$`Pr(>F)`[1]
  }
}

for(j in 3:11){
  for(i in 13:64){
    proteinmod <- lm(covariates[,i] ~ covariates[,j])
    rsquaredvalues[i-12,j-2] <- summary(proteinmod)$r.squared
  }
}

write.csv(anova.pvalues, "ANOVA p-vals.csv")
