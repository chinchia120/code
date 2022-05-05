#　Week 12 2020/05/21 Thursday
#　Linear regression
########################
#　Example 1 
########################
#fertilizer: 0, 0.5, 1.0, 1.5, 2.0, 2.5
#yield: 10, 18, 32, 48, 55, 62
########################
#　Create example table
Ex1 <- data.frame(fertilizer=c(0, 0.5, 1.0, 1.5, 2.0, 2.5), 
			yield=c(10, 18, 32, 48, 55, 62)) #Column name = value 欄位名稱: 數值

Ex1 # view Ex1 確認Ex1內容

Ex1_lm <- lm(yield ~ fertilizer, data=Ex1) #線性迴歸模型函數lm(y ~ x1 + x2 + ... + xn, data = database)

names(Ex1_lm) #檢視model內含資訊名稱
# [1] "coefficients"   "residuals"     "effects"       "rank"         
# [5] "fitted.values" "assign"        "qr"            "df.residual"  
# [9] "xlevels"       "call"          "terms"         "model"     

Ex1_lm$residuals #查看model殘差

Ex1_predicted <- predict(Ex1_lm, Ex1[1]) #計算推估藉由model推估出的yield: predict(model, database[column number])

rd <- data.frame(realY = Ex1$yield, predictedY = Ex1_predicted, residuals = Ex1$yield-Ex1_predicted) #記錄真實Y值、推估Y值及殘差

View(rd) #展示真實Y值、推估Y值及殘差結果

summary(Ex1_lm) #檢視model概要

confint(Ex1_lm) #show β's 95%CI

########################
#　Exercise 1 
########################
#age: 34,39,44,46,48,51,53,60,61,65,66,67
#cholesterol: 141.4,180.5,178.4,212.0,203.2,224.1,186.0,350.0,286.3,287.6,330.3,371.3
########################
#　Create example table
Exe1 <- data.frame(age=c(34,39,44,46,48,51,53,60,61,65,66,67), 
			cholesterol=c(141.4,180.5,178.4,212.0,203.2,224.1,186.0,350.0,286.3,287.6,330.3,371.3)) #Column name = value 欄位名稱: 數值

Exe1 # view Exe1

Exe1_lm <- lm(cholesterol~ age, data=Exe1) #線性迴歸模型函數

summary(Exe1_lm) #檢視model概要


########################
#　Example 2 
########################
#fertilizer: 0, 0.5, 1.0, 1.5, 2.0, 2.5
#yield: 10, 18, 32, 48, 55, 62
########################
#　Create example table
Ex2 <- data.frame(fertilizer=c(0, 0.5, 1.0, 1.5, 2.0, 2.5), 
			yield=c(10, 18, 32, 48, 55, 62)) #Column name = value 欄位名稱: 數值

Ex2 # view Ex2 確認Ex2內容

Ex2_lm <- lm(yield ~ fertilizer, data=Ex2) #線性迴歸模型函數lm(y ~ x1 + x2 + ... + xn, data = database)
summary(Ex2_lm) #檢視model概要

#變方分析(ANOVA)#其他ANOVA方法請參考第七週課程
anova(lm(yield~fertilizer, data=Ex2))  
confint(Ex2_lm) #show β's 95%CI

########################
#　Exercise 3 
########################
#age: 34,39,44,46,48,51,53,60,61,65,66,67
#cholesterol: 141.4,180.5,178.4,212.0,203.2,224.1,186.0,350.0,286.3,287.6,330.3,371.3
########################
#　Create example table
Exe3 <- data.frame(age=c(34,39,44,46,48,51,53,60,61,65,66,67), 
			cholesterol=c(141.4,180.5,178.4,212.0,203.2,224.1,186.0,350.0,286.3,287.6,330.3,371.3)) #Column name = value 欄位名稱: 數值

Exe3 # view Exe3

Exe3_lm <- lm(cholesterol~ age, data=Exe3) #線性迴歸模型函數

summary(Exe3_lm) #檢視model概要

confint(Exe3_lm) #show β's 95%CI


########################
#　Exercise 4 
########################
#編號 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16
#樹木學分數 35,45,58,43,51,24,49,66,37,40,19,26,54,38,36,38
#型態學分數 56,63,91,77,86,51,89,84,49,81,44,69,81,60,56,75
########################
#　Create example table
Exe4 <- data.frame(編號=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16), 
			樹木學分數=c(35,45,58,43,51,24,49,66,37,40,19,26,54,38,36,38),
			型態學分數=c(56,63,91,77,86,51,89,84,49,81,44,69,81,60,56,75))

Exe4 # view Exe4

Exe4_lm <- lm(樹木學分數 ~ 型態學分數, data=Exe4) #線性迴歸模型函數

anova(lm(樹木學分數 ~ 型態學分數, data=Exe4)) #ANOVA for 各別變數

confint(Exe4_lm) #show β's 95%CI

summary(Exe4_lm) #檢視model概要


########################
#　Example 3 
########################
#運動時數=0,0,0,2,2,2,2,4,4,4
#食物攝取量=2,4,6,2,4,6,8,4,6,8
#體重減輕量=6,2,4,8,6,8,5,11,13,12
########################
#　Create example table
Ex5 <- data.frame(運動時數=c(0,0,0,2,2,2,2,4,4,4), 
			食物攝取量=c(2,4,6,2,4,6,8,4,6,8),
			體重減輕量=c(6,2,4,8,6,8,5,11,13,12))

Ex5 # view Ex5

#install.packages("Hmisc")
library("Hmisc")
#透過相關性檢定確認各個變數間之關聯
rcorr(as.matrix(Ex5 ), type=c("pearson"))

Ex5_lm <- lm(體重減輕量 ~ 運動時數 + 食物攝取量, data=Ex5) #線性迴歸模型函數

summary(Ex5_lm) #檢視model概要

####### Adj R2 #########
r2 = 0.8433 #原始決定係數 R square
k = 2 #變數數量
n = 10 #樣本數
p = k + 1 

Radj = 1-(n-1)/(n-p)*(1-r2) 
Radj 
#####################################
####### 顯著性檢定(整體模型配適度) #########
##原始anova()無法展示整體模型配適度，僅能展示各別變數之配適度，因此改用anova_alt()
##第一次使用anova_alt()請預先執行最下方附錄Anova for full model fitting的定義式

#anova_alt(lm(體重減輕量 ~ 運動時數 + 食物攝取量, data=Ex5)) #ANOVA for 整體模型
anova_alt(Ex5_lm) #ANOVA for 整體模型

#####################################
####### 個別迴歸係數顯著性檢定 ##########
summary(Ex5_lm) #檢視model概要
#####################################

########################
#　Exercise 6 
########################
#平均溫度=25,28,30,26,24,28,29,27,32,26,27,28
#日照=8.6,9,8.7,8.5,7,7.8,10,9.2,8.8,7.5,7.2,8.9
#糖分=13,12.5,11.8,10.4,13.5,9.9,9.5,13.2,10,12.8,11.4,10.9
########################
#　Create example table
Exe6 <- data.frame(Temperature=c(25,28,30,26,24,28,29,27,32,26,27,28),
		Daylight=c(8.6,9,8.7,8.5,7,7.8,10,9.2,8.8,7.5,7.2,8.9),
		Sugar=c(13,12.5,11.8,10.4,13.5,9.9,9.5,13.2,10,12.8,11.4,10.9))

Exe6 # view Exe6

summary(Exe6) # 展示簡易描述性統計

#library("Hmisc")
#透過相關性檢定確認各個變數間之關聯
rcorr(as.matrix(Exe6), type=c("pearson"))

Exe6_lm <- lm(Sugar ~ Temperature + Daylight, data=Exe6) #線性迴歸模型函數

summary(Exe6_lm) #檢視model概要

####### 顯著性檢定(整體模型配適度) #########
#原始anova()無法展示整體模型配適度，僅能展示各別變數之配適度，因此改用anova_alt()
#第一次使用anova_alt()請預先執行最下方附錄Anova for full model fitting的定義式
#summary(model)的最下方也有紀錄簡易版的模型整體配適度結果
anova_alt(Exe6_lm) #ANOVA for 整體模型

####### 個別迴歸係數顯著性檢定 ##########
summary(Exe6_lm) #檢視個別迴歸係數顯著性檢定結果

########################
#　Practice 1
########################



########################
#　Practice 2 
########################




#Appendix
###########################
#常用描述性統計
###########################
#mean() 		平均數
#median() 		中位數
#range() 		全距
#quantile() 	四分位數
#IQR() 		四分位差
#summary() 		描述統計摘要
#sd() 		標準差
#var() 		變異數 
#skewness() 	偏度
#kurtosis() 	峰度
###########################
###########################
#相關性檢定
###########################
#cor.test(Exe1$age, Exe1$cholesterol, method="pearson")
#
#library("Hmisc")
#執行複數變數之間的相關性分析
#rcorr(as.matrix(Ex2), type=c("pearson"))
#install.packages("ppcor")
#library("ppcor")
#考量溫度後之雨量與產量之淨相關
#pcor.test(Ex2$rain, Ex2$yield, Ex2$temp, method = "pearson")
###########################
###########################
#deal with collinearity
###########################
#http://www.sthda.com/english/articles/39-regression-model-diagnostics/160-multicollinearity-essentials-and-vif-in-r/
#install.packages("car")
#library("car")
#vif(model)
#vif(Exe6_lm)

############################
############################
#Anova for full model fitting
############################
#(source: https://community.rstudio.com/t/anova-table-for-full-linear-model/42074/10 created by RussS)
############################Definition function
anova_alt = function (object, reg_collapse=TRUE,...) 
{
  if (length(list(object, ...)) > 1L) 
    return(anova.lmlist(object, ...))
  if (!inherits(object, "lm")) 
    warning("calling anova.lm(<fake-lm-object>) ...")
  w <- object$weights
  ssr <- sum(if (is.null(w)) object$residuals^2 else w * object$residuals^2)
  mss <- sum(if (is.null(w)) object$fitted.values^2 else w * 
               object$fitted.values^2)
  if (ssr < 1e-10 * mss) 
    warning("ANOVA F-tests on an essentially perfect fit are unreliable")
  dfr <- df.residual(object)
  p <- object$rank
  if (p > 0L) {
    p1 <- 1L:p
    comp <- object$effects[p1]
    asgn <- object$assign[stats:::qr.lm(object)$pivot][p1]
    nmeffects <- c("(Intercept)", attr(object$terms, "term.labels"))
    tlabels <- nmeffects[1 + unique(asgn)]
    ss <- c(vapply(split(comp^2, asgn), sum, 1), ssr)
    df <- c(lengths(split(asgn, asgn)), dfr)
    if(reg_collapse){
      if(attr(object$terms, "intercept")){
        collapse_p<-2:(length(ss)-1)
        ss<-c(ss[1],sum(ss[collapse_p]),ss[length(ss)])
        df<-c(df[1],sum(df[collapse_p]),df[length(df)])
        tlabels<-c(tlabels[1],"Source")
      } else{
        collapse_p<-1:(length(ss)-1)
        ss<-c(sum(ss[collapse_p]),ss[length(ss)])
        df<-c(df[1],sum(df[collapse_p]),df[length(df)])
        tlabels<-c("Regression")
      }
    }
  }else {
    ss <- ssr
    df <- dfr
    tlabels <- character()
    if(reg_collapse){
      collapse_p<-1:(length(ss)-1)
      ss<-c(sum(ss[collapse_p]),ss[length(ss)])
      df<-c(df[1],sum(df[collapse_p]),df[length(df)])
    }
  }
  
  ms <- ss/df
  f <- ms/(ssr/dfr)
  P <- pf(f, df, dfr, lower.tail = FALSE)
  table <- data.frame(df, ss, ms, f, P)
  table <- rbind(table, 
                 colSums(table))
  if (attr(object$terms, "intercept")){
   table$ss[nrow(table)]<- table$ss[nrow(table)] - table$ss[1]
  }
  table$ms[nrow(table)]<-table$ss[nrow(table)]/table$df[nrow(table)]
  table[length(P):(length(P)+1), 4:5] <- NA
  dimnames(table) <- list(c(tlabels, "Error","Total"), 
                          c("Df","SS", "MS", "F", 
                            "P"))
  if (attr(object$terms, "intercept")){
    table <- table[-1, ]
    table$MS[nrow(table)]<-table$MS[nrow(table)]*(table$Df[nrow(table)])/(table$Df[nrow(table)]-1)
    table$Df[nrow(table)]<-table$Df[nrow(table)]-1
  }
  structure(table, heading = c("Analysis of Variance Table\n"), 
            class = c("anova", "data.frame"))
}

