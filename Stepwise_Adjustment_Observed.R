#############################################################################################
#         This is the code for model 1: No variables will be adjusted in this model         #
#############################################################################################
library(dplyr)
library(ranger)

#set up training and validation splits
df <- read.csv("Ratios_Unadjusted.csv", header = TRUE)
#balance data for age, sex, race, family history of MS, and MS-DSS*
#SWITCH LINE BELOW WITH OTHER ADJUSTMENT METHODS
adjusted <- read.csv("Ratios_of_Residual.csv", header = TRUE)
proteins.adj <- colnames(adjusted)
adjusted <- adjusted %>% mutate(cohort = df$cohort)

#load data for modeling, test should only contain training data
train <- df %>% filter(cohort == "training")
validation <- df %>% filter(cohort == "validation")
adjusted.train <- adjusted %>% filter(cohort == "training")
adjusted.valid <- adjusted %>% filter(cohort == "validation")

trainID <- as.character(train$patient)
validationID <- as.character(validation$patient)

train <- train %>% select(-c("patient", "cohort"))
validation <- validation %>% select(-c("patient", "cohort"))

#possible seeds
seeds<-c(7, 12, 300, 2000, 999, 345, 89)

# Function that returns Root Mean Squared Error
rmse <- function(error)
{
  sqrt(mean(error^2))
}

#data for final predictions
data.predictions.train<-data.frame(row.names=trainID)
data.predictions.valid<-data.frame(row.names=validationID)

mean.data.OOB<-c()
output<-data.frame("n.vars"=c(),"OOB.error"=c(),"var.rm"=c(), "corr.train"=c(), "p.train"=c(), "corr.valid"=c(), "p.valid"=c(), "RMSE.train"=c(), "RMSE.valid"=c())
k=1

  data.imp <-data.frame()
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
    data.imp <- data.imp %>% bind_rows(importance(ranger.out))
  }
    
  #mean OOB error
  mean.data.OOB<- mean(data.OOB,na.rm=T)
    
  #add final data predictions
  data.predictions.train <- data.predictions.train %>% bind_cols(mean.pred = rowMeans(data.predictions.tmp.train,na.rm=T))
  data.predictions.valid <- data.predictions.valid %>% bind_cols(mean.pred = rowMeans(data.predictions.tmp.valid,na.rm=T))
  
  train.pred.corr <- rowMeans(data.predictions.tmp.train,na.rm=T)
  valid.pred.corr <- rowMeans(data.predictions.tmp.valid,na.rm=T)
  
  #correlations between predicted and original MS-DSS
  correlation.train <- sqrt(summary(lm(train.pred.corr ~ train$msdss_last))$r.squared)
  pvalue.train <- summary(lm(train.pred.corr ~ train$msdss_last))$coefficients[2,4]
  correlation.valid <- sqrt(summary(lm(valid.pred.corr ~ validation$msdss_last))$r.squared)
  pvalue.valid <- summary(lm(valid.pred.corr ~ validation$msdss_last))$coefficients[2,4]
  
  #mean variable importance
  data.imp <- data.imp %>% summarise_all(mean)
  
  #RMSE
  rmse.train <- rmse(train$msdss_last - train.pred.corr)
  rmse.valid <- rmse(validation$msdss_last - valid.pred.corr)
  
  #create new dataset without least important variable, then rebuild
  #identify least important variable
  var.rm <- names(which.min(data.imp))
    
  #create output table
  output<-output %>% bind_rows(c("n.vars"=ncol(train)-1, "OOB.error"=mean.data.OOB,"var.rm"=var.rm,"corr.train"=correlation.train, "p.train"=pvalue.train, "corr.valid"=correlation.valid, "p.valid"=pvalue.valid, "RMSE.train"=rmse.train, "RMSE.valid"=rmse.valid))
    
  while(k < 15){
    if(var.rm %in% proteins.adj){
      train[var.rm] <- adjusted.train[var.rm]
      validation[var.rm] <- adjusted.validation[var.rm]
      
      #data.imp <-data.frame()
      data.predictions.tmp.train<-data.frame(matrix(vector(), 335, 1))
      data.predictions.tmp.valid<-data.frame(matrix(vector(), 165, 1))
      data.OOB <- c()
      
      for(s in seeds){
        set.seed(s)
        ranger.out <- ranger(SS~., data=train, importance = "permutation")
        
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
      correlation.train <- cor.test(train$SS, train.pred.corr)$estimate
      pvalue.train <- cor.test(train$SS, train.pred.corr)$p.value
      correlation.valid <- cor.test(validation$SS, valid.pred.corr)$estimate
      pvalue.valid <- cor.test(validation$SS, valid.pred.corr)$p.value
      
      rsquared <- summary(lm(validation$SS ~ valid.pred.corr))$adj.r.squared
      #mean variable importance
      #data.imp <- data.imp %>% summarise_all(mean)
      
      #RMSE
      rmse.train <- rmse(train$SS - train.pred.corr)
      rmse.valid <- rmse(validation$SS - valid.pred.corr)
      
      # #create new dataset without least important variable, then rebuild
      # #identify least important variable
      var.rm <- names(sort(data.imp))[k]
      
      #create output table
      output<-output %>% bind_rows(c("n.vars"=ncol(train)-1, "OOB.error"=mean.data.OOB,"var.rm"=var.rm, "corr.train"=correlation.train, "p.train"=pvalue.train, "corr.valid"=correlation.valid, "p.valid"=pvalue.valid, "RMSE.train"=rmse.train, "RMSE.valid"=rmse.valid, "r.squared"=rsquared))
      
      if(correlation.valid < max(output$corr.valid.cor)){
        train[var.rm] <- original.train[var.rm]
        validation[var.rm] <- original.validation[var.rm]
        output<-output %>% bind_rows(c("n.vars"=ncol(train)-1, "OOB.error"=0,"var.rm"=var.rm, "corr.train.cor"=0, "p.train"=0, "corr.valid.cor"=0, "p.valid"=0, "RMSE.train"=0, "RMSE.valid"=0,"r.squared"=0))
        k=k+1
        var.rm <- names(sort(data.imp))[k]
      }
      if(correlation.valid >= max(output$corr.valid.cor)){
        k=k+1
        var.rm <- names(sort(data.imp))[k]
      }
    }
    if(var.rm %in% proteins.adj == FALSE){
      output<-output %>% bind_rows(c("n.vars"=ncol(train)-1, "OOB.error"=0,"var.rm"=var.rm, "corr.train.cor"=0, "p.train"=0, "corr.valid.cor"=0, "p.valid"=0, "RMSE.train"=0, "RMSE.valid"=0, "r.squared"=0))
      k=k+1
      var.rm <- names(sort(data.imp))[k]
    }
  }

