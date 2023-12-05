#########################
# Set Working Directory #
#########################

# Get your current working directory #
getwd() # Re-check the path for the working directory

# Change your current working directory #
setwd("C:\\Users\\user\\Documents\\code\\introduction-to-scientific-computing-software\\classdata\\C2_ClassData")
getwd() # Re-check the path for the working directory

#####################
# Install R Package #
#####################

# See which packages are installed on your computer #
library()

# How to install new library #
install.packages("MASS")

# How to load a package #
library("MASS")
require("MASS")
library("MASS2")	#show the error message and stopped

#Get information on a package #

#Information for MASS packages#
library(help=stats)     #Information for stats packages#

#Help menu for a function under stats package#
help(lm)	#help menu for linear regression model (lm funcion)
?lm

######################################
# In-Class Practice (Install Package)#
######################################

# Q1. Install the package of "survival"

install.packages("survival")
library("survival")
# make sure you have internet connectivity

# Q2.  Find a description for the package "survival" #

library(help=survival)

# Extra Information (install package manually)#
# set CRAN mirror
# set repository
# select packages

#############################
# Import an example dataset #
#############################

dataset <- read.csv("Week2.csv")
#dataset <- read.table("Week2.csv",header = TRUE, sep = ",")  #Import dataset  


str(dataset)     #Check the variable format
View(dataset)		#Check Dataset
dim(dataset)		#Check Dataset (how many observations and variables)
colnames(dataset)	#Check variable names
summary(dataset)	#Get the summary statistics about the object


#R的因子(factor)變數是專門用來儲存類別資料的變數，它同時具有字串與整數的特性
#使用as.numeric將資料轉為數字格式
#使用as.character將資料轉為字串格式

dataset$Year_NEW <- as.character(dataset$Year)
dataset$SeasonNum_NEW <- as.character(dataset$SeasonNum)
dataset$Rainfall_TEST2<- as.numeric(as.character(dataset$Rainfall_TEST))

############################### 
# Basic Statistical Test in R #
###############################

#mean
mean(dataset$Rainfall,na.rm=TRUE)


#data:向量,矩陣或數字資料
#na.rm:data: exclude(true)/include(false) missing value during calculations;
#default is false

#max
max(dataset$Rainfall,na.rm = TRUE)

#median
median(dataset$Rainfall,na.rm=TRUE)

#min
min(dataset$Rainfall,na.rm = TRUE)

#variance
var(dataset$Rainfall,na.rm = TRUE)

#Standard Deviation
sd(dataset$Rainfall,na.rm = TRUE)

#range
range(dataset$Rainfall,na.rm = TRUE)

#IQR
IQR(dataset$Rainfall,na.rm = TRUE)

#quantile
quantile(dataset$Rainfall,0.25,na.rm = TRUE)
quantile(dataset$Rainfall,0.50,na.rm = TRUE)
quantile(dataset$Rainfall,0.75,na.rm = TRUE)

fivenum(dataset$Rainfall,na.rm = TRUE) #min 25% 50% 75% max


#generate a histogram/boxplot 

#install a new library "lattice"
install.packages("lattice")

# How to load a package #
library("lattice")
require("lattice")

hist(dataset$Rainfall,main="Rainfall in 2008", xlab = "Rainfall(mm)", ylab = "Count") 

histogram(x= ~Rainfall|SeasonNum_NEW,data =dataset, xlab = "Rainfall(mm)",ylab = "Count",layout=c(3,1)) 

histogram(x= ~Rainfall|SeasonNum_NEW,data =dataset, xlab = "Rainfall(mm)",ylab = "Count",type="count",layout=c(3,1))

############
###Output###
############

write.table(dataset,file="D:\\Week2\\output.csv",sep=",",row.names=F, na = "NA")