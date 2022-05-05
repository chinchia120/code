#　Week 8 2022/04/20
#　Chi-square test
#　Fisher’s exact test
#
install.packages("vcd") #安裝vcd套件
########################
#　Example 1
########################
#　Create example table

Ex1 <- data.frame(Phenotype=c("1","1","1","1","2","2","2","2"))
#Column name = value 欄位名稱: 數值
Ex1        # view Ex1 確認Ex1內容

Ex1_Chi <- table(Ex1$Phenotype) 
#generate a contingency table 建立列聯表
Ex1_Chi    # view Ex1_Chi 確認Ex1_Chi內容

chisq.test(Ex1_Chi)  #use “chisq.test()” to obtain the results of chi-square tests 卡方分析，假設每種機率都相同

#Warning message:
#　In chisq.test(Ex1_Chi) : Chi-squared approximation may be incorrect
##because many of the expected values will be very small and therefore the approximations of p may not be right.
##if expected values < 5, this warning will occur
##reference: https://stats.stackexchange.com/questions/81483/warning-in-r-chi-squared-approximation-may-be-incorrect
##
##solve method: add a parameter to calculate simulate p value, simulate.p.value = TRUE
###

chisq.test(Ex1_Chi, simulate.p.value = TRUE) 
# dosen't assign the probability of phenotype 卡方分析，不指定機率(都相同)，計算模擬p值

chisq.test(Ex1_Chi, p = c(3/4, 1/4), simulate.p.value = TRUE) 
# Assign probability. (type1: 3/4, type2: 1/4) 指定機率第一型3/4, 第二型1/4

########################
#　Exercise 1
########################

Exe1_Chi <- data.frame(Gold=10, Red=5, White=45) #直接建立列聯表
Exe1_Chi 

chisq.test(Exe1_Chi, p = c(3/100, 9/100, 88/100), simulate.p.value = TRUE)





########################
#　Example 2
########################
Ex2_Chi <- data.frame(row.names=c("Drug", "No Drug"), Healthy= c(105,110), Unhealthy = c(15,40)) 
#直接建立列聯表
Ex2_Chi

prop.table(Ex2_Chi)  #查看百分比 
prop.table(data.matrix(Ex2_Chi), 2)  #查看各欄百分比 #要先將dataframe格式轉為data.matrix
prop.table(data.matrix(Ex2_Chi), 1)  #查看各列百分比 #要先將dataframe格式轉為data.matrix
# for a matrix 1 indicates rows, 2 indicates columns

chisq.test(Ex2_Chi)

#Dataframe 每一行都可以有不同的資料型態
#Data Matrix 同一行需要統一資料型態

########################
#　Exercise 2
########################

Exe2_Chi <- data.frame(row.names=c("有拜拜", "沒拜拜"), 中獎 = c(38,22), 沒中獎 = c(40,40)) #直接建立列聯表
Exe2_Chi

prop.table(Exe2_Chi)  #查看百分比 
prop.table(data.matrix(Exe2_Chi), 2)  #查看各欄百分比 #要先將dataframe格式轉為data.matrix
prop.table(data.matrix(Exe2_Chi), 1)  #查看各列百分比 #要先將dataframe格式轉為data.matrix

chisq.test(Exe2_Chi)

########################
#　Example 3
########################

Ex3_Chi <- data.frame(row.names=c("改善", "未改善"), A = c(24,12), B = c(30,18), C = c(21,10)) #直接建立列聯表
Ex3_Chi

prop.table(Ex3_Chi)  #查看百分比 
prop.table(data.matrix(Ex3_Chi), 2)  #查看各欄百分比 #要先將dataframe格式轉為data.matrix
prop.table(data.matrix(Ex3_Chi), 1)  #查看各列百分比 #要先將dataframe格式轉為data.matrix

chisq.test(Ex3_Chi)

########################
#　Exercise 3
########################

Exe3_Chi <- data.frame(row.names=c("中獎", "未中獎"), 父 = c(20,10), 母 = c(28,22), 祖父 = c(30,10), 祖母 = c(15,15)) #直接建立列聯表
Exe3_Chi

prop.table(Exe3_Chi)  #查看百分比 
prop.table(data.matrix(Exe3_Chi), 2)  #查看各欄百分比 #要先將dataframe格式轉為data.matrix
prop.table(data.matrix(Exe3_Chi), 1)  #查看各列百分比 #要先將dataframe格式轉為data.matrix

chisq.test(Exe3_Chi)

########################
#　Example 4
########################

Ex4_Fis <- data.frame(row.names=c("Drug", "No Drug"), Healthy= c(21,22), Unhealthy = c(3,8)) #直接建立列聯表
Ex4_Fis

prop.table(Ex4_Fis)  #查看百分比 
prop.table(data.matrix(Ex4_Fis), 2)  #查看各欄百分比 #要先將dataframe格式轉為data.matrix
prop.table(data.matrix(Ex4_Fis), 1)  #查看各列百分比 #要先將dataframe格式轉為data.matrix

chisq.test(Ex4_Fis)

library(vcd) #啟用vcd套件
independence_table(data.matrix(Ex4_Fis)) #建立期望值表

#chisq.test(Ex4_Fis, simulate.p.value = TRUE)

fisher.test(Ex4_Fis)


########################
#　Practice 1 - Chi-squared Test 
########################
data1 <- data.frame(one=20, two=28, three=30, four=15, five=24, six=19)
chisq.test(data1, simulate.p.value=FALSE)

########################
#　Practice 2 - Fisher's exact test 
########################
data2 <- data.frame(row.names=c("yes", "no"), come=c(7,7),  no_come=c(15,2)) 
fisher.test(data2, simulate.p.value=FALSE)