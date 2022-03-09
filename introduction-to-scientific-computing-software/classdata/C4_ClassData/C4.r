#########################
# Set Working Directory #
#########################

# Get your current working directory #

setwd("C:\\Users\\user\\Documents\\code\\introduction-to-scientific-computing-software\\classdata\\C4_ClassData")
getwd()

dataset1<-read.csv("C4_1.csv") 
dataset2<-read.csv("C4_2.csv")
str(dataset1)   #Check the variable format
View(dataset1)		#Check Dataset
dim(dataset1)

###################
#      T-test     #
###################

## One-sample t-test##
#read dataset1
t.test(dataset1$PM25, mu =15) 

##Two-sample paired T-test
#read dataset2
t.test(dataset2$Rainfall_3, dataset2$Rainfall_2,paired=TRUE)
