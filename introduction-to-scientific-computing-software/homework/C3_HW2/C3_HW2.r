setwd("C:\\Users\\user\\Documents\\code\\introduction-to-scientific-computing-software\\homework\\C3_HW2")
getwd()

dataset <- read.csv("C3_HW.csv", encoding="utf-8")

library(psych)
skew(dataset$放)   
kurtosi(dataset$放)

library(nortest)
lillie.test(dataset$放)

boxplot(formula = 放~Season, data = dataset, xlab = "Season", ylab = "放()", col ="blue")

boxplot(formula = 放~year, data = subset(dataset, dataset$city=="穝カ"), xlab = "Year", ylab = "放()", col = "blue")


