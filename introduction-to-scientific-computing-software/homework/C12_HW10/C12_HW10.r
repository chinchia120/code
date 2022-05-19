setwd("C:\\Users\\user\\Documents\\code\\introduction-to-scientific-computing-software\\homework\\C12_HW10")
getwd()

#Q1
data1 <- read.csv("Titanic.csv")

mod_null <- glm(as.factor(Survived) ~ 1, family = "binomial", data = data1)
mod_full <- glm(as.factor(Survived) ~ as.factor(Class) + as.factor(Sex) + as.factor(Age), family = "binomial", data = data1)
mod = step(mod_null, scope = list(lower = mod_null, upper = mod_full), direction = "forward", trace = 1)

#Q2
summary(mod)

#Q3
exp(coef(mod))
exp(confint(mod))