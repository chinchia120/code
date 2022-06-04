setwd("C:\\Users\\user\\Documents\\code\\introduction-to-scientific-computing-software\\homework\\final_practice")
getwd()

data1 <- read.csv("HW_data_1.csv")
data2 <- read.csv("HW_data_2.csv")

#Q1_1
library("Hmisc")
rcorr(as.matrix(data1[,6:15]), type = c("pearson"))

#Q1_2
data1_model <- lm(PM25_obs ~ NO2 + PM10 + O3 + Temperature + RH + NDVI + Majorroad + Waterbody, data = data1)
library("olsrr")
ols_step_both_p(data1_model, penter = 0.05, prem = 0.1, details = TRUE)

#Q1_3
ols_coll_diag(data1_model) 
data1_model <- lm(PM25_obs ~ NO2 + PM10+ NDVI + Majorroad + Waterbody, data = data1)
summary(data1_model)

#Q2_1
data2_mod_null <- glm(as.factor(fracture) ~ 1, family = "binomial", data = data2)
data2_mod_full <- glm(as.factor(fracture) ~ age + as.factor(sex) + weight + height + as.factor(medication) + bmd, family = "binomial", data = data2)
data2_mod = step(data2_mod_null, scope = list(lower = data2_mod_null, upper = data2_mod_full), direction = "both", trace = 1)
summary(data2_mod)

#Q2_2
View(exp(coef(data2_mod_full)))
anova(data2_mod, test = "Chisq")

#Q3_1
data3 <- data.frame(row.names = c("increase", "same"), A = c(32, 13), B = c(27, 8), C = c(19, 5))
library(vcd)
independence_table(data.matrix(data3))
chisq.test(data3, simulate.p.value = FALSE)