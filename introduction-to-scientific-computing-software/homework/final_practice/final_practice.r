setwd("C:\\Users\\user\\Documents\\code\\introduction-to-scientific-computing-software\\homework\\final_practice")
getwd()

data1 <- read.csv("HW_data_1.csv")
data2 <- read.csv("HW_data_2.csv")

#Q1_1
library("Hmisc")
rcorr(as.matrix(data1[,6:15]), type=c("pearson"))

#Q1_2
data1_model <- lm(PM25_obs~NO2+PM10+O3+Temperature+RH+NDVI+Residential_Area+Majorroad+Waterbody, data=data1)
library("olsrr")
ols_step_both_p(data1_model, penter=0.05, prem=0.1, details=TRUE)

#Q1_3
ols_coll_diag(data1_model) 

#Q2_1
