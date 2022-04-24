#Q1
data1 <- data.frame(blue=45, green=55, red=17)
chisq.test(data1, p=c(1/2, 3/8, 1/8), simulate.p.value=FALSE)

#Q2
data2 <- data.frame(row.names=c("recover", "unrecover"), A=c(5,7), B=c(2,6)) 
library(vcd)
independence_table(data.matrix(data2))
fisher.test(data2, simulate.p.value=FALSE)