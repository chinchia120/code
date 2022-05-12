#　Week 11
#　Variable selection and collinearity test
########################
### Import packages ####

install.packages("descriptr")
install.packages("olsrr")		# advancd way to do stepwise regression
install.packages("Hmisc")
install.packages("MASS")

library("olsrr")
library("Hmisc")

#https://cran.r-project.org/web/packages/olsrr/readme/README.html
#https://cran.r-project.org/web/packages/olsrr/vignettes/intro.html




############################################
#  R Script – olsrr                        
#  further summary of lm by olsrr package  #
############################################
#fertilizer: 0, 0.5, 1.0, 1.5, 2.0, 2.5
#yield: 10, 18, 32, 48, 55, 62
############################################
#　Create example table
f1 <- data.frame(fertilizer = c(0, 0.5, 1.0, 1.5, 2.0, 2.5), 
			yield = c(10, 18, 32, 48, 55, 62))
   # Column name = value 欄位名稱: 數值
f1 # view f1 確認f1內容

ols_regress(yield ~ fertilizer, data = f1) #透過olsrr package的ols_regress函數建立線性模型並檢視更詳盡的model概要




#############################
#　Example 1 - Collinearity 
###########################################################################
#一株穗數=17.9,16.1,15.7,19.3,15.9,15,18.5,15.7,19.6,20.1,21.2,17
#日平均氣溫=25.8,24.3,23.5,26.2,24.9,23.8,25.8,24.2,26.5,25.4,26.6,24.3
#日照時數=5.5,4.9,5.2,7.8,5.6,4.5,7.5,5.6,6.6,7.4,8.2,5.2
#公頃產量=6104,5331,5566,6644,5617,5026,6293,5478,6583,6217,6894,5906
############################################################################
ex1 <- data.frame(一株穗數 = c(17.9,16.1,15.7,19.3,15.9,15,18.5,15.7,19.6,20.1,21.2,17),
		日平均氣溫 = c(25.8,24.3,23.5,26.2,24.9,23.8,25.8,24.2,26.5,25.4,26.6,24.3),
		日照時數 = c(5.5,4.9,5.2,7.8,5.6,4.5,7.5,5.6,6.6,7.4,8.2,5.2),
		公頃產量 = c(6104,5331,5566,6644,5617,5026,6293,5478,6583,6217,6894,5906)) 
View(ex1)

# 透過相關性檢定檢視各個變數間之關聯
rcorr(as.matrix(ex1), type = c("pearson"))

ex1_model <- lm(公頃產量 ~ 一株穗數 + 日平均氣溫 + 日照時數, data = ex1)
# 利用olsrr函數展示共線性檢定結果
ols_coll_diag(ex1_model) 

############################
#　Exercise 1 - Collinearity
############################################################
#平均溫度=25,28,30,26,24,28,29,27,32,26,27,28
#日照=8.6,9,8.7,8.5,7,7.8,10,9.2,8.8,7.5,7.2,8.9
#糖分=13,12.5,11.8,10.4,13.5,9.9,9.5,13.2,10,12.8,11.4,10.9
############################################################
exe1 <- data.frame(平均溫度 = c(25,28,30,26,24,28,29,27,32,26,27,28),
		日照 = c(8.6,9,8.7,8.5,7,7.8,10,9.2,8.8,7.5,7.2,8.9),
		糖分 = c(13,12.5,11.8,10.4,13.5,9.9,9.5,13.2,10,12.8,11.4,10.9))
View(exe1)

rcorr(as.matrix(exe1), type=c("pearson"))
exe1_model <- lm(糖分 ~ 平均溫度 + 日照, data=exe1)
ols_coll_diag(exe1_model)

########################
#　Example 2
########################

data(birthwt, package = "MASS") #調用MASS模組提供的範例檔birthwt
View(birthwt)

#source: https://www.rdocumentation.org/packages/MASS/versions/7.3-50/topics/birthwt
#Format
#This data frame contains the following columns:
#low	 新生兒出生時體重小於2.5公斤(Low)
#    indicator of birth weight less than 2.5 kg.
#age	 新生兒母親生產年齡(MAge)
#    mother's age in years.
#lwt	 新生兒母親產前體重(MWt)
#    mother's weight in pounds at last menstrual period.
#race	 新生兒母親族裔(MRace)
#    mother's race (1 = white, 2 = black, 3 = other).
#smoke 新生兒母親懷孕時是否吸菸(MSmoke)
#    smoking status during pregnancy.
#ptl	 新生兒母親過去早產次數(PreLabours)
#    number of previous premature labours.
#ht	 新生兒家族有無高血壓病史(HPT)
#    history of hypertension.
#ui	 新生兒母親有無子宮煩躁症(Uterine_irr.)
#    presence of uterine irritability.
#ftv	 新生兒產後三個月內就醫次數(FTV)
#    number of physician visits during the first trimester.
#bwt	 新生兒體重(BWt)
#    birth weight in grams.

ex2 <- birthwt
colnames(ex2) <- c("Low","MAge","MWt","MRace",
			"MSmoke","PreLabours","HPT",
			"Uterine_irr.","FTV","BWt")
View(ex2)

rcorr(as.matrix(ex2), type=c("pearson"))

#ols_step_forward_p(model, penter = 0.05) #penter(p enter): p value; variables with p value less than penter will enter into the model
#ols_step_backward_p(model, prem = 0.1) #prem(p remove): p value; variables with p more than prem will be removed from the model.
#ols_step_both_p(model, penter = 0.05, prem = 0.1)

############
# 強制輸入法 #
############
ex2_model <- lm(BWt ~ MAge + MWt + as.factor(MRace) + as.factor(MSmoke) + PreLabours + as.factor(HPT) + as.factor(Uterine_irr.) + FTV, data=ex2)
ols_regress(ex2_model)

############
# 順向進入法 #
############
ex2_model <- lm(BWt ~ MAge + MWt + as.factor(MRace) + as.factor(MSmoke) + PreLabours + as.factor(HPT) + as.factor(Uterine_irr.) + FTV, data=ex2)
ols_step_forward_p(ex2_model, penter = 0.05, details = TRUE)

############
# 反向淘汰法 #
############
ex2_model <- lm(BWt ~ MAge + MWt + as.factor(MRace) + as.factor(MSmoke) + PreLabours + as.factor(HPT) + as.factor(Uterine_irr.) + FTV, data=ex2)
ols_step_backward_p(ex2_model, prem = 0.1, details = TRUE)

############
# 逐步分析法 #
############
ex2_model <- lm(BWt ~ MAge + MWt + as.factor(MRace) + as.factor(MSmoke) + PreLabours + as.factor(HPT) + as.factor(Uterine_irr.) + FTV, data=ex2)
ols_step_both_p(ex2_model, penter = 0.05, prem = 0.1, details = TRUE)



########################
#　Example 2
########################################################################
#一株穗數=17.9,16.1,15.7,19.3,15.9,15,18.5,15.7,19.6,20.1,21.2,17
#日平均氣溫=25.8,24.3,23.5,26.2,24.9,23.8,25.8,24.2,26.5,25.4,26.6,24.3
#日照時數=5.5,4.9,5.2,7.8,5.6,4.5,7.5,5.6,6.6,7.4,8.2,5.2
#公頃產量=6104,5331,5566,6644,5617,5026,6293,5478,6583,6217,6894,5906
########################################################################
exe2 <- data.frame(一株穗數=c(17.9,16.1,15.7,19.3,15.9,15,18.5,15.7,19.6,20.1,21.2,17),
		日平均氣溫=c(25.8,24.3,23.5,26.2,24.9,23.8,25.8,24.2,26.5,25.4,26.6,24.3),
		日照時數=c(5.5,4.9,5.2,7.8,5.6,4.5,7.5,5.6,6.6,7.4,8.2,5.2),
		公頃產量=c(6104,5331,5566,6644,5617,5026,6293,5478,6583,6217,6894,5906)) 
View(exe2)
rcorr(as.matrix(exe2), type=c("pearson"))

#強制輸入法
exe2_model <- lm(公頃產量 ~ 一株穗數 + 日平均氣溫 + 日照時數, data=exe2)
ols_regress(exe2_model)

#順向進入法
exe2_model <- lm(公頃產量 ~ 一株穗數 + 日平均氣溫 + 日照時數, data=exe2)
ols_step_forward_p(exe2_model, penter = 0.05, details = TRUE)

#反向淘汰法
exe2_model <- lm(公頃產量 ~ 一株穗數 + 日平均氣溫 + 日照時數, data=exe2)
ols_step_backward_p(exe2_model, prem = 0.1, details = TRUE)

#逐步分析法
exe2_model <- lm(公頃產量 ~ 一株穗數 + 日平均氣溫 + 日照時數, data=exe2)
ols_step_both_p(exe2_model, penter = 0.05, prem = 0.1, details = TRUE)





########################
#　Practice 1
########################
install.packages("MASS") 
data(road, package = "MASS")
practice_data <- road
View(practice_data)

rcorr(as.matrix(practice_data), type=c("pearson"))



#Appendix
############################
#常用描述性統計
############################
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
############################
#其他參考資料
############################
#https://www.rdocumentation.org/packages/lmSupport/versions/1.05
#https://cran.r-project.org/web/packages/lmSupport/index.html
#https://cran.r-project.org/web/packages/lmSupport/lmSupport.pdf

#https://cran.r-project.org/web/packages/pubh/vignettes/regression.html

#Stepwise by olsrr
#https://www.guru99.com/r-simple-multiple-linear-regression.html#8

#Stepwise by leaps package
#http://www.sthda.com/english/articles/37-model-selection-essentials-in-r/154-stepwise-regression-essentials-in-r/

#Stepwise by lm() function
#http://jackthisisamazing.blogspot.com/2016/12/rstepwise-regression_11.html

