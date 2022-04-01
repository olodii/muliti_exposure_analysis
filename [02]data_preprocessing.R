
#####################################################################################
#라이브러리
library(readxl)    ;library(dplyr)     ;library(ggplot2)
library(mediation) ;library(multilevel);library(bda)
library(gvlma)     ;library(stargazer) ;library(medflex)
library(lubridate) ;library(mgcv)

select <- dplyr::select
#----------------------------------------------------------------------------------#
setwd("C://Users//USER//Desktop//최종")
exam2<- read.csv('raw_data_after_merge.csv')

#이용할 변수만 일부 
exam3<-exam2 %>% dplyr:: select(ID,ag,mbi,medu,wks,dbwt,sex,shs,
                                pm10_1T,pm10_2T,pm10_3T,pm10_preg,pm10_6m,pm10_1y,pm10_2y,
                                no2_1T,no2_2T,no2_3T,no2_preg,no2_6m,no2_1y,no2_2y,
                                pm25_1T,pm25_2T,pm25_3T,pm25_preg,pm25_6m,pm25_1y,pm25_2y,
                                mtemp_1T,mtemp_2T,mtemp_3T,mtemp_preg,mtemp_6m,mtemp_1y.ac,mtemp_2y.ac,mtemp_1y,mtemp_2y,
                                mhumi_1T,mhumi_2T,mhumi_3T,mhumi_preg,mhumi_6m,mhumi_1y.ac,mhumi_2y.ac,mhumi_1y,mhumi_2y,
                                green, green2, c_resi_eq.y, c1_resi_eq, c2_resi_eq,
                                yc,yc1,yc2,ym,ym1,ym2,
                                p.date, p2.date, p3.date, b.date, m6.date, m12.date, m24.date)

nrow(exam3)

exam3<-exam3[complete.cases(exam3 %>% select(green,green2)),];nrow(exam3)
exam3 <-exam3[complete.cases(exam3 %>% select(yc,ym)),];nrow(exam3)
exam3<-exam3[complete.cases(exam3 %>%select(pm10_1T:pm25_2y)),];nrow(exam3) 

exam3 <- exam3[!is.na(exam3$wks),];nrow(exam3)
exam3 <- exam3[!is.na(exam3$dbwt),];nrow(exam3)
exam3<-exam3[complete.cases(exam3 %>% select(age,medu,mbi,shs,sex)),];nrow(exam3)
exam3<-subset(exam3,age!=0) ;nrow(exam3)

exam3.s<-subset(exam3,c_resi_eq.y==2) ;nrow(exam3.s) 
exam3.ns<-subset(exam3,c_resi_eq.y==1) ;nrow(exam3.ns)

naniar::miss_var_summary(exam3.s) %>% View()

library(mice)
cov_meteo <- exam3.s %>% select(mtemp_preg,mtemp_1y, mtemp_2y,
                                mhumi_preg,mhumi_1y, mhumi_2y,
                                pm10_preg, pm10_2y, no2_preg, no2_2y, pm25_preg, pm25_2y)
imputed_cov_met <- mice(cov_meteo, m=5, maxit = 100, method = 'pmm', seed = 500) # numeric value 이므로 pmm 사용
summary(imputed_cov_met)

cov_meteo.im <- cov_meteo
cov_meteo.im[is.na(cov_meteo.im$mtemp_1y),'mtemp_1y'] <- apply(imputed_cov_met$imp$mtemp_1y,1,mean)
cov_meteo.im[is.na(cov_meteo.im$mtemp_2y),'mtemp_2y'] <- apply(imputed_cov_met$imp$mtemp_2y,1,mean)
cov_meteo.im[is.na(cov_meteo.im$mhumi_1y),'mhumi_1y'] <- apply(imputed_cov_met$imp$mhumi_1y,1,mean)
cov_meteo.im[is.na(cov_meteo.im$mhumi_2y),'mhumi_2y'] <- apply(imputed_cov_met$imp$mhumi_2y,1,mean)

exam3.s$mtemp_1y <- cov_meteo.im$mtemp_1y
exam3.s$mtemp_2y <- cov_meteo.im$mtemp_2y
exam3.s$mhumi_1y <- cov_meteo.im$mhumi_1y
exam3.s$mhumi_2y <- cov_meteo.im$mhumi_2y

exam3.s[is.na(exam3.s$mtemp_1y.ac),'mtemp_1y.ac'] <- apply(exam3.s[is.na(exam3.s$mtemp_1y.ac),] %>% select(mtemp_6m,mtemp_1y),1,mean)
exam3.s[is.na(exam3.s$mtemp_2y.ac),'mtemp_2y.ac'] <- apply(exam3.s[is.na(exam3.s$mtemp_2y.ac),] %>% select(mtemp_6m,mtemp_1y,mtemp_2y),1,mean)
exam3.s[is.na(exam3.s$mhumi_1y.ac),'mhumi_1y.ac'] <- apply(exam3.s[is.na(exam3.s$mhumi_1y.ac),] %>% select(mhumi_6m,mhumi_1y),1,mean)
exam3.s[is.na(exam3.s$mhumi_2y.ac),'mhumi_2y.ac'] <- apply(exam3.s[is.na(exam3.s$mhumi_2y.ac),] %>% select(mhumi_6m,mhumi_1y,mhumi_2y),1,mean)

naniar::miss_var_summary(exam3.s) %>% View()

#-------------------------------------------------------------------------
dat <- exam3.s

# 단위 변환
dat %>% select(pm10_1T:pm25_2y) %>% summary()

#NO2 단위 ppm->ppb 변환 
dat$no2_1T.r  =dat$no2_1T*1000
dat$no2_2T.r  =dat$no2_2T*1000
dat$no2_3T.r  =dat$no2_3T*1000
dat$no2_6m.r  =dat$no2_6m*1000
dat$no2_preg.r=dat$no2_preg*1000
dat$no2_1y.r  =dat$no2_1y*1000
dat$no2_2y.r  =dat$no2_2y*1000

#PM10 단위 10ug/m3 증가당, NO2  단위 10ppb 증가당
dat$pm10_1T.u  =dat$pm10_1T/10  ;dat$no2_1T.u  =dat$no2_1T.r/10; dat$pm25_1T.u  =dat$pm25_1T/10 
dat$pm10_2T.u  =dat$pm10_2T/10  ;dat$no2_2T.u  =dat$no2_2T.r/10; dat$pm25_2T.u  =dat$pm25_2T/10 
dat$pm10_3T.u  =dat$pm10_3T/10  ;dat$no2_3T.u  =dat$no2_3T.r/10; dat$pm25_3T.u  =dat$pm25_3T/10 
dat$pm10_preg.u=dat$pm10_preg/10;dat$no2_preg.u=dat$no2_preg.r/10; dat$pm25_preg.u  =dat$pm25_preg/10 
dat$pm10_6m.u  =dat$pm10_6m/10  ;dat$no2_6m.u  =dat$no2_6m.r/10; dat$pm25_6m.u  =dat$pm25_6m/10 
dat$pm10_1yc.u  =dat$pm10_1y/10  ;dat$no2_1yc.u  =dat$no2_1y.r/10; dat$pm25_1yc.u  =dat$pm25_1y/10 
dat$pm10_2yc.u  =dat$pm10_2y/10  ;dat$no2_2yc.u  =dat$no2_2y.r/10; dat$pm25_2yc.u  =dat$pm25_2y/10 


dat %>% select(pm10_1T.u:pm25_2yc.u) %>% summary()

# 타입 변환
dat$medu <- factor(dat$medu)

table(dat$sex)
dat_female <- subset(dat, sex==1);nrow(dat_female)
dat_male <- subset(dat, sex==0);nrow(dat_male)

b1 <- dat %>% select(ID,age,mbi,medu,wks,dbwt,sex,shs,
                     pm10_1T.u,pm10_2T.u,pm10_3T.u,pm10_preg.u,pm10_6m.u,
                     no2_1T.u,no2_2T.u,no2_3T.u,no2_preg.u,no2_6m.u,
                     pm25_1T.u,pm25_2T.u,pm25_3T.u,pm25_preg.u,pm25_6m.u,
                     mtemp_1T,mtemp_2T,mtemp_3T,mtemp_preg,mtemp_6m,
                     mhumi_1T,mhumi_2T,mhumi_3T,mhumi_preg,mhumi_6m,
                     green, green2, c_resi_eq.y, 
                     yc,ym,
                     p.date, p2.date, p3.date, b.date, m6.date)
naniar::miss_var_summary(b1)

b2 <- dat %>% select(ID,age,mbi,medu,wks,dbwt,sex,shs,
                     pm10_1T.u,pm10_2T.u,pm10_3T.u,pm10_preg.u,pm10_6m.u,pm10_1yc.u,pm10_2yc.u,
                     no2_1T.u,no2_2T.u,no2_3T.u,no2_preg.u,no2_6m.u,no2_1yc.u,no2_2yc.u,
                     pm25_1T.u,pm25_2T.u,pm25_3T.u,pm25_preg.u,pm25_6m.u,pm25_1yc.u,pm25_2yc.u,
                     mtemp_1T,mtemp_2T,mtemp_3T,mtemp_preg,mtemp_6m,mtemp_1y,mtemp_2y,
                     mhumi_1T,mhumi_2T,mhumi_3T,mhumi_preg,mhumi_6m,mhumi_1y,mhumi_2y,
                     green, green2, c1_resi_eq,c2_resi_eq,
                     yc,yc1,ym,ym1,
                     p.date, b.date, m12.date,m24.date)

b2.s <- subset(b2, b2$c1_resi_eq == 2)
b2.ns <- subset(b2, b2$c1_resi_eq == 1)
b2.s <- b2.s[complete.cases(b2.s$yc1),]
naniar::miss_var_summary(b2.s)

b3 <- b2.s %>% select(ID,age,mbi,medu,wks,dbwt,sex,shs,
                      pm10_1T.u,pm10_2T.u,pm10_3T.u,pm10_preg.u,pm10_6m.u,pm10_1yc.u,pm10_2yc.u,
                      no2_1T.u,no2_2T.u,no2_3T.u,no2_preg.u,no2_6m.u,no2_1yc.u,no2_2yc.u,
                      pm25_1T.u,pm25_2T.u,pm25_3T.u,pm25_preg.u,pm25_6m.u,pm25_1yc.u,pm25_2yc.u,
                      mtemp_1T,mtemp_2T,mtemp_3T,mtemp_preg,mtemp_6m,mtemp_1y,mtemp_2y,
                      mhumi_1T,mhumi_2T,mhumi_3T,mhumi_preg,mhumi_6m,mhumi_1y,mhumi_2y,
                      green, green2, c2_resi_eq,
                      yc2,ym2,
                      p.date, b.date, m24.date)

b3.s <- subset(b3, b3$c2_resi_eq == 2);nrow(b3.s)
b3.ns <- subset(b3, b3$c2_resi_eq == 1)
b3.s <- b3.s[complete.cases(b3.s$yc2),];nrow(b3.s)
naniar::miss_var_summary(b3.s)

Data <- NULL
Data[[1]] <- b1
Data[[2]] <- b2.s
Data[[3]] <- b3.s

save(Data, file='data_formation.RData')