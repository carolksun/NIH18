protein.ratios <- data.frame(matrix(vector(), 194, 25))
protein.ratios[1] <- all.proteins[1]
protein.ratios[2] <- all.proteins[2]
colnames(protein.ratios)[1] <- "patient"
colnames(protein.ratios)[2] <- "msdss_last"

#data.imp <-data.frame()
data.predictions.tmp.train<-data.frame(row.names=trainID)
data.predictions.tmp.valid<-data.frame(row.names=validationID)
data.OOB <- c()

for(s in seeds){
  set.seed(s)
  ranger.out <- ranger(msdss_last~., data=train, importance = "permutation", num.trees = 800)
  
  #store predictions
  data.predictions.tmp.train <- data.predictions.tmp.train %>% bind_cols(pred = ranger.out$predictions)
  data.predictions.tmp.valid <- data.predictions.tmp.valid %>% bind_cols(pred = predict(ranger.out,validation)$predictions)
  
  #store OOB error
  data.OOB <- c(data.OOB, ranger.out$prediction.error)
  
  # store variable importance
  #data.imp <- data.imp %>% bind_rows(importance(ranger.out))
}

#mean OOB error
mean.data.OOB<- mean(data.OOB,na.rm=T)

#add final data predictions
data.predictions.train <- data.predictions.train %>% bind_cols(mean.pred = rowMeans(data.predictions.tmp.train,na.rm=T))
data.predictions.valid <- data.predictions.valid %>% bind_cols(mean.pred = rowMeans(data.predictions.tmp.valid,na.rm=T))

train.pred.corr <- rowMeans(data.predictions.tmp.train,na.rm=T)
valid.pred.corr <- rowMeans(data.predictions.tmp.valid,na.rm=T)

#correlations between predicted and original MS-DSS
correlation.train <- cor.test(train$msdss_last, train.pred.corr)$estimate
pvalue.train <- cor.test(train$msdss_last, train.pred.corr)$p.value
correlation.valid <- cor.test(validation$msdss_last, valid.pred.corr)$estimate
pvalue.valid <- cor.test(validation$msdss_last, valid.pred.corr)$p.value

#mean variable importance
#data.imp <- data.imp %>% summarise_all(mean)

#RMSE
rmse.train <- rmse(train$msdss_last - train.pred.corr)
rmse.valid <- rmse(validation$msdss_last - valid.pred.corr)

# #create new dataset without least important variable, then rebuild
# #identify least important variable
# var.rm <- names(sort(data.imp))[1]
#   
# ratio.num.den <- c(as.character(ratios$v1[which(var.rm == ratios$dfnames)[[1]]]),as.character(ratios$v2[which(var.rm == ratios$dfnames)[[1]]]))
# #create output table
output<-output %>% bind_rows(c("n.vars"=ncol(train)-1, "OOB.error"=mean.data.OOB,"var.rm.num"=ratio.num.den[1],"var.rm.den"=ratio.num.den[2], "corr.train"=correlation.train, "p.train"=pvalue.train, "corr.valid"=correlation.valid, "p.valid"=pvalue.valid, "RMSE.train"=rmse.train, "RMSE.valid"=rmse.valid))
