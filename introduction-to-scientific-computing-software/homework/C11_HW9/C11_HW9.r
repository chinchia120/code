data(mtcars, package="datasets")
data1 <- mtcars
View(data1)

#Q1
library("Hmisc")
rcorr(as.matrix(data1), type=c("pearson"))

#Q2
data1_model <- lm(mpg~cyl+disp+hp+drat+wt+qsec+vs+am+gear+carb, data=data1)
library("olsrr")
ols_regress(data1_model)

#Q3
ols_step_forward_p(data1_model, penter=0.05, details=TRUE)

#Q4
ols_step_backward_p(data1_model, prem=0.1, details=TRUE)

#Q5
ols_step_both_p(data1_model, penter=0.05, prem=0.1, details=TRUE)

