#　Logistic Regression
########################
### Import packages ####
install.packages("ResourceSelection") # For Hosmer & Lemeshow Test
install.packages("caret") # For Example Database 

library("ResourceSelection")
library("caret")
########################
#　Example 1 - Logistic regression
########################
data(GermanCredit, package = "caret")
ex1 <- GermanCredit

View(ex1)

#範例: 新申請人信用評估資訊
#Class	(Good/Bad)	綜合評價信用好/不好
#Age	年齡
#ForeignWorker	(1/0)	是否為外籍勞工
#Property.RealEstate	(1/0)	有無房地產
#Housing.Own	(1/0)	有無自己的房子
#CreditHistory.Critical	(1/0)	有無信用帳戶
#其他略

#Generalized linear model  #glm(y ~ x, family = data distribution (poission/gaussian/binomial... etc.), data = database)
mod <- glm(Class ~ Age + ForeignWorker + Property.RealEstate 
		+ Housing.Own + CreditHistory.Critical, family="binomial", data = ex1)

summary(mod)

hoslem.test(mod$y, fitted(mod)) #hoslem.test(model$y, fitted(model))
?hoslem.test

#Calculate Odds ratio
exp(coef(mod))
#可以用View函數令其以表格顯示
View(exp(coef(mod)))

#Get 95% Confidence Interval
confint(mod)	#95% CI of Coefficient
exp(confint(mod))	
View(exp(confint(mod)))



########################
#　Example 2 - variable selected logistic regression
########################
data(GermanCredit, package = "caret")
ex2 <- GermanCredit

View(ex2)

mod.null <- glm(Class ~ 1, family="binomial", data = ex2)
mod.full <- glm(Class ~ ., family="binomial", data = ex2)

f.model = step(mod.null, scope =list(lower=mod.null, upper=mod.full), direction = "forward", trace = 1)

summary(f.model)

#b.model = step(mod.full, direction = "backward", trace = 1)
#summary(b.model)

########################
#　Practice - Backward logistic regression
########################





#Appendix
########################
# 其他 Logistic model 係數顯著性測驗
########################
#https://www.r-bloggers.com/evaluating-logistic-regression-models/
########################
#Omnibus6
anova(mod, test ="Chisq")

########################
#Cox & Snell R 平方值, Nagelkerke R 平方值
install.packages("pscl")
library("pscl")

pR2(mod)  # look for 'McFadden'

########################
#Hosmer 和 Lemeshow
install.packages("ResourceSelection")
library("ResourceSelection")

hoslem.test(mod$y, fitted(mod), g=10)

########################
#Wald 值
install.packages("survey")
library("survey")

regTermTest(mod, "ForeignWorker")

