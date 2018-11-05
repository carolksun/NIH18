library(ggplot2)

#loading dataset
OOB <- read.csv("step simulation OOB.csv", header = T)

#calculating percent decrease
median.percent.decrease <- 100*(OOB$Median-OOB$Unadjusted)/OOB$Unadjusted

#creating groups for the dataset for the graph
effectsize <- c(0.01, 0.11, 0.21, 0.31, 0.41, 0.51, 0.01, 0.11, 0.21, 0.31, 0.41, 0.51, 0.01, 0.11, 0.21, 0.31, 0.41, 0.51,
                0.01, 0.11, 0.21, 0.31, 0.41, 0.51, 0.01, 0.11, 0.21, 0.31, 0.41, 0.51, 0.01, 0.11, 0.21, 0.31, 0.41, 0.51)
allelefreq <- c(0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.45, 0.45, 0.45, 0.45, 0.45, 0.45,
                0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.45, 0.45, 0.45, 0.45, 0.45, 0.45)
samplesize <- c(200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200,
                500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500)

#grouping to make colors for the graph
group <- paste(allelefreq, ".", samplesize, sep = "")
graphing <- cbind(effectsize, group, as.numeric(median.percent.decrease), allelefreq, samplesize)
graphing <- as.data.frame(graphing)

#colors
cols <- c("#2faeae", "#dea9ec", "#80e59b")

#graph
ggplot(data = graphing, aes(x = effectsize, y = median.percent.decrease, group = group)) +
  geom_line(aes(color=allelefreq, linetype=samplesize), size = 1.2) +
  geom_point(aes(color=allelefreq))+
  ylab("Percent Difference of OOB Error") +
  xlab("Effect Size") +
  scale_color_manual(name = "Allele Frequency", values=cols) +
  scale_linetype_manual(name  = "Sample Size", values=c("solid", "dotdash")) +
  ggtitle("Percent Decrease of Out of Bag Error after Median Adjustment") +
  theme_grey() + theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        plot.title = element_text(size=18), legend.position = "none")





#loading dataset
rsquared <- read.csv("step simulation rsquared.csv", header = T)

#calculating percent increase
rsquared.percent.increase <- 100*(rsquared$Median-rsquared$Unadjusted)/rsquared$Unadjusted

#creating groups for the dataset for the graph
effectsize <- c(0.01, 0.11, 0.21, 0.31, 0.41, 0.51, 0.01, 0.11, 0.21, 0.31, 0.41, 0.51, 0.01, 0.11, 0.21, 0.31, 0.41, 0.51,
                0.01, 0.11, 0.21, 0.31, 0.41, 0.51, 0.01, 0.11, 0.21, 0.31, 0.41, 0.51, 0.01, 0.11, 0.21, 0.31, 0.41, 0.51)
allelefreq <- c(0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.45, 0.45, 0.45, 0.45, 0.45, 0.45,
                0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.45, 0.45, 0.45, 0.45, 0.45, 0.45)
samplesize <- c(200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200,
                500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500)

#grouping to make colors for the graph
group <- paste(allelefreq, ".", samplesize, sep = "")
graphing <- cbind(effectsize, group, as.numeric(rsquared.percent.increase), allelefreq, samplesize)
graphing <- as.data.frame(graphing)

#colors
cols <- c("#2faeae", "#dea9ec", "#80e59b")

#graph
ggplot(data = graphing, aes(x = effectsize, y = rsquared.percent.increase, group = group)) +
  geom_line(aes(color=allelefreq, linetype=samplesize), size = 1.2) +
  geom_point(aes(color=allelefreq)) +
  ylab("Percent Difference of Variance Explained") +
  xlab("Effect Size") +
  scale_color_manual(name = "Allele Frequency", values=cols) +
  scale_linetype_manual(name  = "Sample Size", values=c("solid", "dotdash")) +
  ggtitle("Percent Increase of Variance Explained after Median Adjustment") +
  theme_grey() + theme(axis.text=element_text(size=15),
                       axis.title=element_text(size=15),
                       plot.title = element_text(size=18), legend.position = "none")
