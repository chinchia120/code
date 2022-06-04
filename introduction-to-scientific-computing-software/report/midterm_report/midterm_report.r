#Q1
setwd("C:\\Users\\user\\Documents\\code\\introduction-to-scientific-computing-software\\report\\midterm_report")
getwd() 
data <- read.csv("data.csv")
str(data)
max(subset(data,data$Year == 2010)$NTD,na.rm = TRUE)
min(subset(data,data$Year == 2010)$NTD,na.rm = TRUE)
mean(subset(data,data$Year == 2010)$NTD,na.rm=TRUE)
sd(subset(data,data$Year == 2010)$NTD,na.rm = TRUE)
quantile(subset(data,data$Year == 2010)$NTD,0.25,na.rm = TRUE)
quantile(subset(data,data$Year == 2010)$NTD,0.50,na.rm = TRUE)
quantile(subset(data,data$Year == 2010)$NTD,0.75,na.rm = TRUE)
shapiro.test(subset(data,data$Year == 2010)$NTD) 
histogram(x= ~NTD|City,data = data, xlab = "NT",ylab = "Count",layout=c(13,2)) 
boxplot(formula= NTD~City,data =data, xlab = "City", ylab = "NT", col ="blue") 

#Q2
t.test(data_spend$NTD, mu = mean(data$NTD,na.rm=TRUE)) 

#Q3
data5<-read.csv("data5.csv")
aov <- aov(NTD~ factor(City), data=data5) 
summary(aov)#Check p-value 
leveneTest(data5$NTD, data5$City, center=mean)
TukeyHSD(aov)    
LDuncan(aov,"City") 