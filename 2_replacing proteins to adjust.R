if((ratio.num.den[1] %in% proteins.adj || ratio.num.den[2] %in% proteins.adj) == TRUE){
  if(ratio.num.den[1] %in% proteins.adj){
    all.proteins[ratio.num.den[1]] <- adjusted[ratio.num.den[1]]
  }
  
  if(ratio.num.den[2] %in% proteins.adj){
    all.proteins[ratio.num.den[2]] <- adjusted[ratio.num.den[2]]
  }
  
  for(i in 1:23){
    protein.ratios[i+2] <- all.proteins[[as.character(ratios$v1[i])]] / all.proteins[[as.character(ratios$v2[i])]]
  }
  
  protein.ratios <- protein.ratios %>% mutate(cohort = df$cohort)
  train <- protein.ratios %>% filter(cohort == "training")
  validation <- protein.ratios %>% filter(cohort == "validation")
  
  train <- train %>% select(-c("patient", "cohort"))
  validation <- validation %>% select(-c("patient", "cohort"))
}