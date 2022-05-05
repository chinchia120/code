setwd("C:\\Users\\user\\Documents\\code\\introduction-to-scientific-computing-software\\homework\\C9_HW7")
getwd()

#Q1
data1 <- read.csv("Tainan_pollution_station.csv")
cor.test(data1$PM25, data1$NO2, method="spearman") 

#Q2
library("Hmisc")
rcorr(as.matrix(data1[,5:8]), type=c("spearman"))

library("ppcor")
pcor.test(data1$PM25, data1$NO2, c(data1$Temperature, data1$RH), method="spearman")
