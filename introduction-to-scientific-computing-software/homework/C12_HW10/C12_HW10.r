setwd("C:\\Users\\user\\Documents\\code\\introduction-to-scientific-computing-software\\homework\\C12_HW10")
getwd()

#Q1
data1 <- read.csv("Titanic.csv")

data1$Class_fac <- as.factor(data1$Class)
data1$Sex_fac <- as.factor(data1$Sex)
data1$Age_fac <- as.factor(data1$Age)
data1$Survived_fac <- as.factor(data1$Survived)

library("ResourceSelection")
mod_null <- glm(Class ~ 1, family = "binomial", data = data1)
mod_full <- glm(Class ~ 5:8, family = "binomial", data = data1)
step(mod_null, scope = list(lower = mod_null, upper = mod_full), direction = "forward", trace = 1)