setwd("C:\\Users\\user\\Documents\\code\\introduction-to-scientific-computing-software\\report\\final_report")
getwd()

data_outcome <- read.csv("平均每戶全年經常性支出.csv")
data_income <- read.csv("平均每戶全年經常性收入.csv")
data_stationary <- read.csv("平均每戶書報雜誌文具支出占消費支出比率.csv")
data_area <- read.csv("平均每人居住面積(坪).csv")
data_education <- read.csv("15歲以上民間人口之教育程度結構-大專及以上(％).csv")

#View(data_outcome)
#View(data_income)
#View(data_stationary)
#View(data_area)
#View(data_education)

#ANOVA_Q2
cor.test(data_outcome$臺北市, data_stationary$臺北市 , method = "spearman")

#連續型

#類別型


