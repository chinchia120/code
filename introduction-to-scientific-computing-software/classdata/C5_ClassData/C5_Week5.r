install.packages("lattice")
install.packages("nortest")     #Lilliefors(Kolmogorov-Smirnov) Normality Test
install.packages("car") 

#########################
# Set Working Directory #
#########################

# Get your current working directory #
setwd("C:\\Users\\user\\Documents\\code\\introduction-to-scientific-computing-software\\classdata\\C5_ClassData")
getwd()
data1 <- read.csv("C5.csv")  #例題操作資料
data2 <- read.csv("C3.csv")  #課堂演示資料
str(data1)   #Check the variable format
#str(data2)
View(data1)		#Check Dataset
#View(data2)
dim(data1)		#Check Dataset (how many observations and variables)
#dim(data2)


###################
#      資料轉換     #
###################

#########################
#   Test of Normality   #
#########################
#1.shapiro.test
shapiro.test(data2$Rainfall)    #The Shapiro-Wilk Normality Test

#2.lillie.test
library(nortest)
lillie.test(data2$Rainfall)

#3.
#開方根轉換
data2$Rainfall_Change1 <- sqrt(data2$Rainfall)
#對數轉換
data2$Rainfall_Change2 <- log(data2$Rainfall+1)
#View(data2)
#檢定是否為常態分布
lillie.test(data2$Rainfall_Change1)
lillie.test(data2$Rainfall_Change2)

###############
#  histogram  #  
###############
library("lattice")
bartlett.test(NDVI~ T.NT, data=data1)     

hist(data2$Rainfall,main="Rainfall in 2008",xlab="Rainfall(mm)")
hist(data2$Rainfall_Change1,main="Rainfall in 2008",xlab="Rainfall(mm)")
hist(data2$Rainfall_Change2,main="Rainfall in 2008",xlab="Rainfall(mm)")


#########
##例題一##
###################
#Test of Normality#
###################
shapiro.test(data1$NDVI[data1$T.NT==1])   
#The Shapiro-Wilk Normality Test (sample size<50) 
    
lillie.test(data1$NDVI[data1$T.NT==2])
#Lilliefors(Kolmogorov-Smirnov) Normality Test (sample size>50)

###################
#      T-test     #
###################
    ##變方檢定##
 #1. Bartlett Test#
bartlett.test(NDVI~ T.NT, data=data1)           

 #2. Levene's Test#                         
library(car)
leveneTest(data1$NDVI, data1$T.NT, center=mean)

## Two-sample t-test##
t.test(data1$NDVI[data1$T.NT==1], data1$NDVI[data1$T.NT==2], var.equal=TRUE)





#########
##例題二##
###################
#Test of Normality#
###################
shapiro.test(data1$Rainfall[data1$T.NT==1])        
#The Shapiro-Wilk Normality Test (sample size<50)

lillie.test(data1$Rainfall[data1$T.NT==2])
#Lilliefors(Kolmogorov-Smirnov) Normality Test (sample size>50)

###################
#      T-test     #
###################
    ##變方檢定##
 #1. Bartlett Test#
bartlett.test(Rainfall~ T.NT, data=data1) 

 #2. Levene's Test#
leveneTest(data1$Rainfall, data1$T.NT, center=mean)        

## Two-sample t-test##
t.test(data1$Rainfall[data1$T.NT==1], data1$Rainfall[data1$T.NT==2])

#########################################

#Q1
data_weight <- read.csv("Weight.csv")

shapiro.test(data_weight$Weight[data_weight$Brand=="A"])
shapiro.test(data_weight$Weight[data_weight$Brand=="B"])   

bartlett.test(Weight~ Brand, data=data_weight)

t.test(data_weight$Weight[data_weight$Brand=="A"], data_weight$Weight[data_weight$Brand=="B"], var.equal=TRUE)

#Q2
data_height <- read.csv("Height.csv")

shapiro.test(data_height$Height[data_height$class=="A"])
shapiro.test(data_height$Height[data_height$class=="B"]) 

bartlett.test(Height~ class, data=data_height)

t.test(data_height$Height[data_height$class=="A"], data_height$Height[data_height$class=="B"], var.equal=FALSE)

#Q3
data_price <- read.csv("Price.csv")

shapiro.test(data_price$Price)
hist(data_price$Price, main="Price", xlab="Price")

data_price$Price_sqrt <- sqrt(data_price$Price)
shapiro.test(data_price$Price_sqrt)
hist(data_price$Price_sqrt, main="Price in Sqrt", xlab="Price Sqrt")

data_price$Price_log <- log(data_price$Price)
shapiro.test(data_price$Price_log)
hist(data_price$Price_log, main="Price in Log", xlab="Price Log")
