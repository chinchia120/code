setwd("C:\\Users\\user\\Documents\\code\\introduction-to-scientific-computing-software\\homework\\C5_HW4")
getwd()

#Q1
data1 <- read.csv("C5_Q1.csv")

shapiro.test(data1$??�kdata1$�W??==1])
shapiro.test(data1$??�kdata1$�W??==2])

bartlett.test(??���W??, data=data1)

t.test(data1$??�kdata1$�W??==1], data1$??�kdata1$�W??==2], var.equal=TRUE)

#Q2
data2 <- read.csv("C5_Q2.csv")

shapiro.test(data2$??f�[?)
hist(data2$??f�[?, main="??f�[?", xlab="??f�[?")

data2$??f�[?_sqrt <- sqrt(data2$??f�[?)
shapiro.test(data2$??f�[?_sqrt)
hist(data2$??f�[?_sqrt, main="??f�[? in Sqrt", xlab="??f�[??")
