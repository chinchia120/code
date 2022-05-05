#　Week 9
#　Pearson
#　Spearman
#　Partial correlation
#
install.packages("Hmisc")
#https://cran.r-project.org/web/packages/Hmisc/index.html
#Hmisc: Harrell Miscellaneous
#Contains many functions useful for data analysis, 
#high-level graphics, utility operations, 
#functions for computing sample size and power, 
#importing and annotating datasets, imputing missing values, 
#advanced table making, variable clustering, character string manipulation, 
#conversion of R objects to LaTeX and html code, and recoding variables.

install.packages("ppcor")
#https://cran.r-project.org/web/packages/ppcor/index.html
#Partial and Semi-Partial (Part) Correlation
#Calculates partial and semi-partial (part) correlations along with p-value.

########################
#　Example 1
########################
#　Create example table

Ex1 <- data.frame(坪數=c(32,34,52,48,44,42,64,36,60,40), 
			收入=c(66,78,135,111,84,150,168,102,180,120)) #Column name = value 欄位名稱: 數值

Ex1 # view Ex1 確認Ex1內容

cor.test(Ex1$坪數, Ex1$收入, method="pearson")  #相關性檢定，用Pearson法，檢定Ex1的坪數跟Ex1的收入之間的相關性

########################
#　Exercise 1
########################


########################
#　Exercise 2
########################


########################
#　Example 2
########################
#　Create example table

Ex2 <- data.frame(temp=c(15,14,16,18,25,22,23,20),
			rain=c(23.5,20.4,22.8,25.9,28.9,27.1,26.8,25.2), 
			yield=c(34.0,23.0,33.5,35.8,44.6,41.2,45.4,39.0)) #Column name = value 欄位名稱: 數值

Ex2 # view Ex2 確認Ex2內容

library("Hmisc")

#執行複數變數之間的相關性分析
rcorr(as.matrix(Ex2), type=c("pearson"))
#?rcorr

library("ppcor")

#考量溫度後之雨量與產量之淨相關
pcor.test(Ex2$rain, Ex2$yield, Ex2$temp, method = "pearson")
#?pcor.test

########################
#　Practice 1 - Spearman
########################

########################
#　Practice 2 - Partial correlation
########################

