setwd("C:\\Users\\user\\Documents\\code\\introduction-to-scientific-computing-software\\homework\\midterm_practice")
getwd()

#Q1
data1 <- read.csv("期中前課堂練習.csv")
data2 <- subset(data1, data1$year=="2009")

summary(subset(data2, data2$city=="桃園市"))
summary(subset(data2, data2$city=="臺北市"))
round(sd(subset(data2, data2$city=="桃園市")$溫度,na.rm = TRUE),4)
round(sd(subset(data2, data2$city=="臺北市")$溫度,na.rm = TRUE),4)

library(psych)
skew(subset(data2, data2$city=="桃園市")$溫度)   
skew(subset(data2, data2$city=="臺北市")$溫度)  
kurtosi(subset(data2, data2$city=="桃園市")$溫度)
kurtosi(subset(data2, data2$city=="臺北市")$溫度)

library(nortest)
lillie.test(subset(data2, data2$city=="桃園市")$溫度)
lillie.test(subset(data2, data2$city=="臺北市")$溫度)

data1$month_ch <- as.character(data1$month)

library("lattice")
histogram(x = ~溫度|month_ch, data = data1, xlab = "Temperature(℃)", ylab = "Count", type = "count", layout = c(4,3))
boxplot(formula = 溫度~month, data = data1, xlab = "month", ylab = "Temperature(℃)", col = "blue")

#Q2
data3 <- read.csv("weight.csv")

t.test(data3$weight_1, mu = 20.80)
t.test(data3$weight_2, mu = 23.35)

t.test(data3$score_1, data3$score_2, paired = TRUE)

#Q3_1
dim(subset(data2, data2$city=="桃園市"))
dim(subset(data2, data2$city=="臺北市"))

library(nortest)
lillie.test(data2$溫度[data2$city=="桃園市"])
lillie.test(data2$溫度[data2$city=="臺北市"])

library(car)
leveneTest(data2$溫度, data2$city=="桃園市", center=mean)
leveneTest(data2$溫度, data2$city=="臺北市", center=mean)

t.test(data2$溫度[data2$city=="桃園市"], data2$溫度[data2$city=="臺北市"], var.equal=FALSE)

#Q3_2
dim(subset(data2, data2$city=="高雄市"))
dim(subset(data2, data2$city=="臺北市"))

library(nortest)
lillie.test(data2$溫度[data2$city=="高雄市"])
lillie.test(data2$溫度[data2$city=="臺北市"])

library(car)
leveneTest(data2$溫度, data2$city=="高雄市", center=mean)
leveneTest(data2$溫度, data2$city=="臺北市", center=mean)

t.test(data2$溫度[data2$city=="高雄市"], data2$溫度[data2$city=="臺北市"], var.equal=FALSE)

#Q4
data4 <- read.csv("cost.csv")

shapiro.test(data4$NT_dollar)
hist(data4$NT_dollar, main="NT_dollar", xlab="NTD")

data4$NT_dollar_sqrt <- sqrt(data4$NT_dollar)
shapiro.test(data4$NT_dollar_sqrt)
hist(data4$NT_dollar_sqrt, main="NT_dollar in sqrt", xlab="NTD")

data4$NT_dollar_log <- log(data4$NT_dollar)
shapiro.test(data4$NT_dollar_log)
hist(data4$NT_dollar_log, main="NT_dollar in log", xlab="NTD")

#Q5
summary(aov(溫度~ factor(year), data=data1)) 

library(car)
leveneTest(data1$溫度, data1$year, center=mean)

library(PMCMRplus)
summary(gamesHowellTest(aov(溫度~ factor(year), data=data1)))

data5 <- subset(data2, data2$city=="臺北市" | data2$city=="新竹縣" | data2$city=="臺南市")
summary(aov(溫度~ factor(city), data=data5))