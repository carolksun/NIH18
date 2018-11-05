
#library(dplyr)
#library(ranger)
k=1
#set up training and validation splits
df <- read.csv("Ratios_Unadjusted.csv", header = TRUE)
#balance data for age, sex, race, family history of MS, and MS-DSS*
#SWITCH LINE BELOW WITH OTHER ADJUSTMENT METHODS
adjusted <- read.csv("Ratios_of_Residual.csv", header = TRUE)
proteins.adj <- colnames(adjusted)
#adjusted <- adjusted %>% mutate(cohort = df$cohort)

ratios <- read.csv("ratios.csv", header = TRUE)
rationames <- colnames(df)[4:26]
ratios <- ratios %>% mutate(dfnames = rationames)

all.proteins <- read.csv("all proteins.csv", header = TRUE)
#all.proteins.org <- read.csv("all proteins.csv", header = TRUE)
#load data for modeling, test should only contain training data


train <- df %>% filter(cohort == "training")
validation <- df %>% filter(cohort == "validation")
original.train <- df %>% filter(cohort == "training")
original.validation <- df %>% filter(cohort == "validation")
adjusted.train <- adjusted %>% filter(cohort == "training")
adjusted.validation <- adjusted %>% filter(cohort == "validation")

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
  correlation.train <- cor.test(train$msdss_last, train.pred.corr)$estimate
  pvalue.train <- cor.test(train$msdss_last, train.pred.corr)$p.value
  correlation.valid <- cor.test(validation$msdss_last, valid.pred.corr)$estimate
  pvalue.valid <- cor.test(validation$msdss_last, valid.pred.corr)$p.value
  
  #mean variable importance
  data.imp <- data.imp %>% summarise_all(mean)
  
  #RMSE
  rmse.train <- rmse(train$msdss_last - train.pred.corr)
  rmse.valid <- rmse(validation$msdss_last - valid.pred.corr)
  
  # #create new dataset without least important variable, then rebuild
  # #identify least important variable
   var.rm <- names(sort(data.imp))[1]

  #create output table
  output<-output %>% bind_rows(c("n.vars"=ncol(train)-1, "OOB.error"=mean.data.OOB,"var.rm"=var.rm, "corr.train"=correlation.train, "p.train"=pvalue.train, "corr.valid"=correlation.valid, "p.valid"=pvalue.valid, "RMSE.train"=rmse.train, "RMSE.valid"=rmse.valid))
    
  if(var.rm %in% proteins.adj){
    train[var.rm] <- adjusted.train[var.rm]
    validation[var.rm] <- adjusted.validation[var.rm]
    
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
    var.rm <- names(sort(data.imp))[k]
    
    #create output table
    output<-output %>% bind_rows(c("n.vars"=ncol(train)-1, "OOB.error"=mean.data.OOB,"var.rm"=var.rm, "corr.train"=correlation.train, "p.train"=pvalue.train, "corr.valid"=correlation.valid, "p.valid"=pvalue.valid, "RMSE.train"=rmse.train, "RMSE.valid"=rmse.valid))
  }

