library("ranger")
library("dplyr")

Scores <- read.table("Scores.txt", header = TRUE)
Master <- read.table("Master.txt", header = TRUE)
ID <- read.csv("List of Patients.csv", header = TRUE)
ID <- as.character(ID$x)

Master2 <- merge(Master, Scores, by="patient")

rfdata <- Master2 %>% select(msdss_last, contains("SL"))

Predicterror = data.frame(matrix(vector(), 1, 3))
colnames(Predicterror) <- c("Trees", "mtry", "error")

Templist <- data.frame(matrix(vector(), 1, 3))
colnames(Templist) <- c("Trees", "mtry", "error")

t=500
while(t <= 1000){
  for(i in 7:30){
  Templist <- data.frame(matrix(vector(), 1, 3))
  colnames(Templist) <- c("Trees", "mtry", "error")
  rfmodelresults <- ranger(msdss_last ~., data = rfdata, num.trees = t, mtry = i)
  Templist$Trees <- t
  Templist$mtry <- i
  Templist$error <- rfmodelresults$prediction.error
  Predicterror <- rbind(Predicterror, Templist)
  }
  t = t+50
}

plot(Predicterror$error~Predicterror$mtry)
Predicterror$Trees <- as.factor(Predicterror$Trees)

Predicterror2 <- Predicterror
Predicterror2$Trees <- as.numeric(as.character(Predicterror$Trees))


for(j in seq(500,1000,50)){
  dferror <- Predicterror %>% filter(Trees == j) %>% select(mtry, error)
  print(ggplot(dferror, aes(y=error, x=mtry)) + geom_line() + ggtitle(paste("Trees =", j)))
}

rfmodelresults <- ranger(msdss_last ~., data = rfdata, num.trees = 800, mtry = 19, importance = "permutation")

train=sample(1:nrow(rfdata),120)
testresults <- randomForest(msdss_last ~., data = rfdata, subset = train, importance = TRUE, mtry = 19)

varImpPlot(testresults)
#mtry = 19, trees = 800
#Read Variable Importance Methods, Splitting Training and Validation, how to show results of random forests
