#########################
# Set Working Directory #
#########################

# Get your current working directory #
getwd() # Re-check the path for the working directory

# Change your current working directory #
setwd("C:\\Users\\user\\Documents\\code\\introduction-to-scientific-computing-software\\classdata\\C3_ClassData")
getwd() # Re-check the path for the working directory

#############################
# Import an example dataset #
#############################

dataset <- read.csv("C3.csv", encoding = "utf-8")
#dataset <- read.table("C3.csv",header = TRUE, sep = ",")  #Import dataset  


str(dataset)     #Check the variable format
View(dataset)		#Check Dataset
dim(dataset)		#Check Dataset (how many observations and variables)
colnames(dataset)	#Check variable names
summary(dataset)	#Get the summary statistics about the object

######################
# generate a boxplot #
######################

#install.packages("lattice")
library("lattice")
require("lattice")

boxplot(formula = Rainfall~SeasonNum, data = dataset, xlab = "Season", ylab = "Rainfall(mm)", col ="red") 

boxplot(formula = Rainfall~Town_name, data = dataset, xlab = "Town", ylab = "Rainfall(mm)", col = "blue")

dataset2 <- na.omit(dataset)

##Skewness & Kurtosis ##
install.packages("psych")
library(psych)
skew(dataset2$Rainfall)   #Skewness
kurtosi(dataset2$Rainfall)  #Kurtosis

###################
#testing normality#
###################
shapiro.test(dataset2$Rainfall)        #The Shapiro-Wilk Normality Test

install.packages("nortest")     #Lilliefors(Kolmogorov-Smirnov) Normality Test
library(nortest)
lillie.test(dataset2$Rainfall)
