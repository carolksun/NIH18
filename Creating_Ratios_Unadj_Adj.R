ratios <- read.csv("ratios.csv", header = TRUE)
new.proteins <- read.csv("updated proteins.csv", header = TRUE)
old.proteins <- read.csv("old protein levels.csv", header = TRUE)
Scores <- read.table("Scores.txt", header = TRUE)
adjusted <- read.csv("newAdjustedProteinsResidual.csv", header = TRUE)

Scores <- Scores %>% select(patient, msdss_last)
colnames(new.proteins)[1] <- "patient"
colnames(old.proteins)[1] <- "patient"

all.proteins <- merge(new.proteins, old.proteins, by = "patient")
all.proteins <- merge(Scores, all.proteins, by = "patient")

all.proteins[[colnames(adjusted)[1]]] <- adjusted[[1]]
all.proteins[[colnames(adjusted)[2]]] <- adjusted[[2]]
all.proteins[[colnames(adjusted)[3]]] <- adjusted[[3]]
all.proteins[[colnames(adjusted)[4]]] <- adjusted[[4]]
all.proteins[[colnames(adjusted)[5]]] <- adjusted[[5]]

protein.ratios <- data.frame(matrix(vector(), 194, 25))
protein.ratios[1] <- all.proteins[1]
protein.ratios[2] <- all.proteins[2]
colnames(protein.ratios)[1] <- "patient"
colnames(protein.ratios)[2] <- "msdss"

for(i in 1:23){
  protein.ratios[i+2] <- all.proteins[[as.character(ratios$v1[i])]] / all.proteins[[as.character(ratios$v2[i])]]
}

rf.protein.ratios <- protein.ratios[-1]
ranger.out <- ranger(msdss~., data=rf.protein.ratios, importance = "permutation", num.trees = 800)
ranger.out <- ranger(msdss~., data=rf.snps.protein.ratios, importance = "permutation", num.trees = 800)
