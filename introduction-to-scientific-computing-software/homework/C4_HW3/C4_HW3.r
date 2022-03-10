setwd("C:\\Users\\user\\Documents\\code\\introduction-to-scientific-computing-software\\homework\\C4_HW3")
getwd()

dataset1 <- read.csv("C4_1_HW.csv") 
dataset2 <- read.csv("C4_2_HW.csv")

t.test(dataset1$材積, mu = 20)

t.test(dataset2$A飼料體重增重, dataset2$B飼料體重增重, paired = TRUE)