setwd("C:\\Users\\user\\Documents\\code\\introduction-to-scientific-computing-software\\homework\\C6_HW5")
getwd()

#Q1
data1<-read.csv("C6_HW1.csv")
summary(aov(Yield~ factor(Brand_ID), data=data1)) 

#Q2
data2<-read.csv("C6_HW2.csv")
summary(aov(Speed~ factor(Train_ID), data=data2)) 
library(car)
leveneTest(data2$Speed, data2$Train_ID, center=mean)
TukeyHSD(aov(Speed~ factor(Train_ID), data=data2))