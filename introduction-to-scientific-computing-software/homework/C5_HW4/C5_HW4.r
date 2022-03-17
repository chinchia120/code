setwd("C:\\Users\\user\\Documents\\code\\introduction-to-scientific-computing-software\\homework\\C5_HW4")
getwd()

#Q1
data1 <- read.csv("C5_Q1.csv")

shapiro.test(data1$糤[data1$箎==1])
shapiro.test(data1$糤[data1$箎==2])

bartlett.test(糤~箎, data=data1)

t.test(data1$糤[data1$箎==1], data1$糤[data1$箎==2], var.equal=TRUE)

#Q2
data2 <- read.csv("C5_Q2.csv")

shapiro.test(data2$痜Ω计)
hist(data2$痜Ω计, main="痜Ω计", xlab="痜Ω计")

data2$痜Ω计_sqrt <- sqrt(data2$痜Ω计)
shapiro.test(data2$痜Ω计_sqrt)
hist(data2$痜Ω计_sqrt, main="痜Ω计 in Sqrt", xlab="痜Ω计")
