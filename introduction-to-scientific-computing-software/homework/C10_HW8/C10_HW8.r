setwd("C:\\Users\\user\\Documents\\code\\introduction-to-scientific-computing-software\\homework\\C10_HW8")
getwd()

#Q1
data1 <- read.csv("HW_Database.csv")
library("Hmisc")
rcorr(as.matrix(data1[,5:8]), type=c("pearson"))

#Q2
summary(lm(NO2~PM25+Temperature+RH, data=data1))

#Q4
dim(data1)
n=39
k=3
round(1-((n-1)/(n-k-1))*(1-0.60),2)