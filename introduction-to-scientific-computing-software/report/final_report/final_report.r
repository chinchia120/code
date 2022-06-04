setwd("C:\\Users\\user\\Documents\\code\\introduction-to-scientific-computing-software\\report\\final_report")
getwd()

data_outcome <- read.csv("平均每戶全年經常性支出.csv", header = TRUE)
data_income <- read.csv("平均每戶全年經常性收入.csv", header = TRUE)
data_stationary <- read.csv("平均每戶書報雜誌文具支出占消費支出比率.csv", header = TRUE)
data_area <- read.csv("平均每人居住面積(坪).csv", header = TRUE)
data_education <- read.csv("15歲以上民間人口之教育程度結構-大專及以上(％).csv", header = TRUE)
data_dominate <- read.csv("平均每戶可支配所得.csv", header = TRUE)
data_spend <- read.csv("平均每戶消費支出.csv", header = TRUE)
data_house <- read.csv("自有住宅率.csv", header = TRUE)
data_taipei <- read.csv("data_Taipei.csv", header = TRUE)
data_2020 <- read.csv("data_2020.csv", header = TRUE)

#敘述統計
str(data_income)
max(subset(data_income, data_income$year == 2020)$ntd, na.rm = TRUE)
min(subset(data_income, data_income$year == 2020)$ntd, na.rm = TRUE)
mean(subset(data_income, data_income$year == 2020)$ntd, na.rm = TRUE)
sd(subset(data_income, data_income$year == 2020)$ntd, na.rm = TRUE)
quantile(subset(data_income, data_income$year == 2020)$ntd, 0.25, na.rm = TRUE)
quantile(subset(data_income, data_income$year == 2020)$ntd, 0.50, na.rm = TRUE)
quantile(subset(data_income, data_income$year == 2020)$ntd, 0.75, na.rm = TRUE)

shapiro.test(subset(data_income, data_income$year == 2020)$ntd)

str(data_outcome)
max(subset(data_outcome, data_outcome$year == 2020)$ntd, na.rm = TRUE)
min(subset(data_outcome, data_outcome$year == 2020)$ntd, na.rm = TRUE)
mean(subset(data_outcome, data_outcome$year == 2020)$ntd, na.rm = TRUE)
sd(subset(data_outcome, data_outcome$year == 2020)$ntd, na.rm = TRUE)
quantile(subset(data_outcome, data_outcome$year == 2020)$ntd, 0.25, na.rm = TRUE)
quantile(subset(data_outcome, data_outcome$year == 2020)$ntd, 0.50, na.rm = TRUE)
quantile(subset(data_outcome, data_outcome$year == 2020)$ntd, 0.75, na.rm = TRUE)

shapiro.test(subset(data_outcome, data_outcome$year == 2020)$ntd)

#盒形圖直方圖
library("lattice")
histogram(x = ~ ntd|city, data = data_dominate, xlab = "NT", ylab = "COUNT", layout = c(12,2))
boxplot(formula = ntd ~ city, data = data_dominate, xlab = "City", ylab = "NT", col = "blue")

#t-test
t.test(subset(data_spend, data_spend$city == "臺北市")$ntd, mu = mean(data_spend$ntd, na.rm = TRUE))

data_spend_city <- read.csv("平均每戶消費支出(分縣市).csv", header = TRUE)
str(data_spend_city)
leveneTest(data_spend_city$ntd, data_spend_city$city, center = mean)
t.test(data_spend_city$ntd[data_spend_city$city == "臺北市"], data_spend_city$ntd[data_spend_city$city == "高雄市"], var.equal = TRUE)

#ANOVA_Q1
data_outcome_city <- read.csv("平均每戶全年經常性支出(分縣市).csv", header = TRUE)
str(data_outcome_city)
aov <- aov(ntd ~ factor(city), data = data_outcome_city)
summary(aov)
library(car)
leveneTest(data_outcome_city$ntd, data_outcome_city$city, center = mean)
library(PMCMRplus)
res <- gamesHowellTest(aov)
summary(res)

#ANOVA_Q2
data_house_city <- read.csv("自有住宅率(分縣市).csv", header = TRUE)
str(data_house_city)
aov <- aov(rate ~ factor(city), data = data_house_city)
summary(aov)
library(car)
leveneTest(data_house_city$rate, data_house_city$city, center = mean)
library(PMCMRplus)
res <- gamesHowellTest(aov)
summary(res)

#ANOVA_Q3
data_outcome_city = subset(data_outcome, data_outcome$city == "臺北市" | city == "臺中市" | city == "高雄市" | city == "花蓮縣")
summary(aov(ntd ~ factor(city), data = data_outcome_city)) 
library(car)
leveneTest(data_outcome_city$ntd, data_outcome_city$city, center = mean)
TukeyHSD(aov(ntd ~ factor(city), data = data_outcome_city))

#ANOVA_Q4
data_stationary_city = subset(data_stationary, data_stationary$city == "臺北市" | city == "臺中市" | city == "高雄市" | city == "花蓮縣")
summary(aov(rate ~ factor(city), data = data_stationary_city)) 
library(car)
leveneTest(data_stationary_city$rate, data_stationary_city$city, center = mean)
TukeyHSD(aov(rate ~ factor(city), data = data_stationary_city))

#連續型
library("Hmisc")
rcorr(as.matrix(data_taipei[,2:9]), type = c("spearman"))
summary(lm(Income ~ Education + Area + Dominate + Income + Outcome + Stationary + Spend + House, data = data_taipei))

#類別型
mod_null <- glm(as.factor(Six_Capital) ~ 1, family = "binomial", data = data_2020)
mod_full <- glm(as.factor(Six_Capital) ~ Education + Area + Dominate + Income + Outcome + Stationary + Spend + House, family = "binomial", data = data_2020)
mod = step(mod_null, scope = list(lower = mod_null, upper = mod_full), direction = "forward", trace = 1)
summary(mod)

View(exp(coef(mod)))
View(exp(confint(mod)))
