#########################
# Set Working Directory #
#########################

# Get your current working directory #

setwd("C:\\Users\\user\\Documents\\code\\introduction-to-scientific-computing-software\\classdata\\C6_ClassData")
getwd()

data1<-read.csv("C6_1.csv") 
data2<-read.csv("C6_2.csv") 
str(data1)   #Check the variable format
View(data1)		#Check Dataset
dim(data1)		#Check Dataset (how many observations and variables)

###################
#      T-test     #
###################

#1#
## One-way ANOVA##
 #First method#         
aov1 <- aov(yeild~ factor(brand_ID), data=data1) #factor() for categorical variable
summary(aov1)#Check p-value 

 #Second method#
anova(lm(yeild~factor(brand_ID), data=data1))      

#Test for Homogeneity of Variance 
install.packages("car")                          
library(car)
leveneTest(data1$yeild, data1$brand_ID, center=mean)

##Variation equal-ANOVA post-hoc test##
 #Tukey#
TukeyHSD(aov1)    
 
 #Duncan#
#install.packages("laercio") 
library(laercio)  
LDuncan(aov1,"brand_ID") 

####
install.packages("remotes")
library(remotes)
install_version("laercio","1.0.1")

#2#
## One-way ANOVA##
 #First method#         
aov2 <- aov(PM25~ factor(county_ID), data=data2)
summary(aov1)

 #Second method#
anova(lm(PM25~factor(county_ID), data=data2))      

#Test for Homogeneity of Variance 
#install.packages("car")                          
library(car)
leveneTest(data2$PM25, data2$county_ID, center=mean)

##Variation is not equal-ANOVA post-hoc test##
 #Games-Howell#
install.packages("PMCMRplus")
library(PMCMRplus)
res<-gamesHowellTest(aov1)
summary(res)
