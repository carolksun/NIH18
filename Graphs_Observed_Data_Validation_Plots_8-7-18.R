#load in observed MS-DSS and predicted MS-DSS
msdss <- read.csv("unadjustedratios.csv", header = TRUE)
valid.pred.corr <- read.csv("unadjusted validation prediction.csv", header = T)
valid.pred.corr <- read.csv("training unadj predictions.csv", header = T)

validation <- msdss %>% filter(cohort == "validation")
training <- msdss %>% filter(cohort == "training")
graph2 <- cbind(training$msdss_last, train.pred.corr)
graph2 <- as.data.frame(graph2)
#Creating graph dataset
graph1 <- validation$msdss_last
graph1 <- cbind(graph1, valid.pred.corr$x)
graph1 <- as.data.frame(graph1)

graph1 <- training$msdss_last
graph1 <- cbind(graph1, valid.pred.corr$x)
graph1 <- as.data.frame(graph1)

#graph
ggplot(graph1, aes(x=graph1, y=V2)) + geom_point(size = 2, color = "slateblue") +
  geom_smooth(method=lm, color = "turquoise3") +
  xlab("") + scale_y_continuous("", breaks = waiver(), labels = c(1.5,2.0,2.5,3.0,3.5,4.0), limits = c(1.5,4)) +
  theme_grey() + theme(axis.text=element_text(size=15),
                       axis.title=element_text(size=15),
                       plot.title = element_text(size=20))

#load in predicted MS-DSS for best stepwise residual
bestresidual <- read.table("bestresidual.txt", header = T)
graph2 <- cbind(validation$msdss_last, bestresidual)

#graph
ggplot(graph2, aes(x=validation$msdss_last, y=x)) + geom_point(size = 2, color = "salmon") +
  ggtitle(paste("Residual Validation (r = ", round(cor.test(graph2[[1]], graph2[[2]])$estimate,5), ")", sep="")) +
  geom_smooth(method=lm, color = "violetred1") +
  xlab("Observed MS-DSS") + scale_y_continuous("Predicted MS-DSS", breaks = waiver(), labels = c(1.5,2.0,2.5,3.0,3.5,4.0), limits = c(1.5,4))+
  theme_grey() + theme(axis.text=element_text(size=15),
                       axis.title=element_text(size=15),
                       plot.title = element_text(size=20))

ggplot(graph2, aes(x=V1, y=train.pred.corr)) + geom_point(size = 2, color = "salmon") +
  ggtitle(paste("Residual Validation (r = ", round(cor.test(graph2[[1]], graph2[[2]])$estimate,5), ")", sep="")) +
  geom_smooth(method=lm, color = "violetred1") +
  xlab("Observed MS-DSS") + scale_y_continuous("Predicted MS-DSS", breaks = waiver(), labels = c(1.5,2.0,2.5,3.0,3.5,4.0), limits = c(1.5,4))+
  theme_grey() + theme(axis.text=element_text(size=15),
                       axis.title=element_text(size=15),
                       plot.title = element_text(size=20))

