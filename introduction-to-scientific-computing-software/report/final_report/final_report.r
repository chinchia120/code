setwd("C:\\Users\\user\\Documents\\code\\introduction-to-scientific-computing-software\\report\\final_report")
getwd()

data_outcome <- read.csv("平均每戶全年經常性支出.csv", header = TRUE)
data_income <- read.csv("平均每戶全年經常性收入.csv", header = TRUE)
data_stationary <- read.csv("平均每戶書報雜誌文具支出占消費支出比率.csv", header = TRUE)
data_area <- read.csv("平均每人居住面積(坪).csv", header = TRUE)
data_education <- read.csv("15歲以上民間人口之教育程度結構-大專及以上(％).csv", header = TRUE)

#View(data_outcome)
#View(data_income)
#View(data_stationary)
#View(data_area)
#View(data_education)

#ANOVA_Q3
data_merge1 = subset(data_outcome, data_outcome$city == "臺北市" | city == "臺中市" | city == "高雄市" | city == "花蓮縣")
summary(aov(ntd ~ factor(city), data = data_merge1)) 

#ANOVA_Q3
data_merge2 = subset(data_stationary, data_stationary$city == "臺北市" | city == "臺中市" | city == "高雄市" | city == "花蓮縣")
summary(aov(ntd ~ factor(city), data = data_merge2)) 

#連續型

#類別型
