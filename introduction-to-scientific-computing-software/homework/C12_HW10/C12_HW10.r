setwd("C:\\Users\\user\\Documents\\code\\introduction-to-scientific-computing-software\\homework\\C12_HW10")
getwd()

#Q1
data1 <- read.csv("Titanic.csv")

library("ResourceSelection")
mod_null <- glm(Class ~ 1, family = "binomial", data = data1)
mod_full <- glm(Class ~ as.factor(Class) + as.factor(Sex) + as.factor(Age) + as.factor(Survived), family = "binomial", data = data1)
step(mod_null, scope = list(lower = mod_null, upper = mod_full), direction = "forward", trace = 1)