setwd("C:\\Users\\user\\Documents\\code\\introduction-to-scientific-computing-software\\homework\\C2_HW1")
getwd()

dataset_2009 <- read.csv("C2_HW_2009.csv")
summary(dataset_2009)
round(sd(dataset_2009$放, na.rm = TRUE),2)

dataset_2010 <- read.csv("C2_HW_2010.csv")
summary(dataset_2010)
round(sd(dataset_2010$放, na.rm = TRUE),2)

dataset_2011 <- read.csv("C2_HW_2011.csv")
summary(dataset_2011)
round(sd(dataset_2011$放, na.rm = TRUE),2)

dataset_hist <- read.csv("C2_HW_hist.csv")
histogram(x= ~放|Season, data = dataset_hist, xlab = "Temperature(?J)", ylab = "Count", layout = c(2,2))
