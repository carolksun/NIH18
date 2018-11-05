mse <- function(sm) 
  mean(sm$residuals^2)

#random forest
rf.protein.ratios <- protein.ratios[-1]
ranger.out <- ranger(msdss~., data=rf.protein.ratios, importance = "permutation", num.trees = 800)
pred = ranger.out$predictions

#linear regression
lm.out <- lm(msdss~., data=rf.protein.ratios)
pred = predict.lm(lm.out)

#PCR
pcr_model <- pcr(msdss~., data = rf.protein.ratios, scale = TRUE, validation = "CV") 
validationplot(pcr_model, val.type = "MSEP")
pred = predict(pcr_model,rf.protein.ratios,ncomp = 9)

#PLS, ridge, lasso and elastic net regressions
pls.model = plsr(msdss~., data = rf.protein.ratios, validation = "CV")
validationplot(pls.model)
summary(pls.model)
pred = predict(pls.model, ncomp = 12)
