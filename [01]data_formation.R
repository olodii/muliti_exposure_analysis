##################################
### Data     #####################
##################################

setwd("C://Users//USER//Desktop//최종")
#----------------------------------------------------------------------------------#
#라이브러리
library(readxl)    ;library(dplyr)     ;library(ggplot2)
library(mediation) ;library(multilevel);library(bda)
library(gvlma)     ;library(stargazer) ;library(medflex)
library(lubridate) ;library(mgcv)

select <- dplyr::select
#----------------------------------------------------------------------------------#
trap<-as.data.frame(read_excel("air pollution_PM10_NO2.xls",sheet=1))

pm10 <- trap %>% dplyr::select(pm10_6m,pm10_1y,pm10_2y)
pm10_pst <- data.frame(pm10_pst=apply(pm10,1,mean))
no2 <- trap %>% dplyr::select(no2_6m,no2_1y,no2_2y)
no2_pst <- data.frame(no2_pst=apply(no2,1,mean))

trap <- cbind(trap, pm10_pst, no2_pst)

trap %>% names()

pm25<-as.data.frame(read_excel("pm25_cmaq.xlsx",sheet=1)) %>% select(id,m1,m2,m3,pregnancy,y1,y2,y3) 
colnames(pm25) <- c('id','pm25_1T','pm25_2T','pm25_3T','pm25_preg','pm25_6m','pm25_1y','pm25_2y')
post_pm25 <- pm25 %>% select(c('pm25_6m','pm25_1y','pm25_2y'))
pm25$pm25_pst <- apply(post_pm25,1,mean)

traps <- merge(pm25, trap, all.x=TRUE, by.x='id',by.y='ID')
colnames(traps)
#----------------------------------------------------------------------------------#
lee.data<-as.data.frame(read_excel("green_and_development.xlsx",sheet=1))

names(trap)[1]    ="ID"
names(lee.data)[2]="ID"

lee.data$dob=dmy(lee.data$dob)

lee.data$b_year =year(lee.data$dob)
lee.data$b_month=month(lee.data$dob)
lee.data$b_day  =day(lee.data$dob)

#Season variable (Winter vs other season)
lee.data$season=factor(ifelse(lee.data$b_month==12 | lee.data$b_month<=2,0,1))

bout<-read_excel("..//data/measurement.xlsx",sheet=6) %>% dplyr:: select(id,d_del_ty,d_del_date_y:d_ga_d)

names(bout)[1]="ID"

#TRAPs와 merge
exam1<-merge(lee.data,traps,by.x="ID",by.y='id',all.x=T)
exam1<-merge(exam1,bout,by="ID",all.x=T)

exam1$b.date <- ymd(paste0(exam1$d_del_date_y,'-',exam1$d_del_date_m, '-', exam1$d_del_date_d))
str(exam1$d_ga_wks)
exam1$p.days <- as.numeric(exam1$d_ga_wks)*7 + exam1$d_ga_d


exam1$p.date <- exam1$b.date - days(exam1$p.days)
exam1$p2.date <- exam1$p.date + weeks(12)
exam1$p3.date <- exam1$p2.date + weeks(15)
exam1$m6.date <- exam1$b.date + months(6)
exam1$m12.date <- exam1$b.date + months(12)
exam1$m24.date <- exam1$b.date + months(24)

#-----------------------------------------------------------------------#
# 기상청 변수 연계, 기온 & 습도
meteo<-read.csv("wea_2000_2016.csv")
meteo$key=paste0(meteo$ddate,"-",meteo$sido)

a <-read_excel("..//_data//examination.xlsx",sheet=1) %>% dplyr::select(id,y_yr,y_mo,y_day)
a1<-read_excel("..//data//examination.xlsx",sheet=2) %>% dplyr::select(id,y1_yr,y1_mo,y1_day)
a2<-read_excel("..//data//examination.xlsx",sheet=3) %>% dplyr::select(id,y2_yr,y2_mo,y2_day)

a$date1=with(a ,paste0(y_yr ,"-",ifelse(y_mo<10 ,paste0(0,y_mo) ,y_mo),"-" ,ifelse(y_day<10 ,paste0(0,y_day) ,y_day)))
a1$date2=with(a1,paste0(y1_yr,"-",ifelse(y1_mo<10,paste0(0,y1_mo),y1_mo),"-",ifelse(y1_day<10,paste0(0,y1_day),y1_day)))
a2$date3=with(a2,paste0(y2_yr,"-",ifelse(y2_mo<10,paste0(0,y2_mo),y2_mo),"-",ifelse(y2_day<10,paste0(0,y2_day),y2_day)))

a$bayley_6m_date =as.Date(ifelse(a$date1=="NA-NA-NA" ,NA,a$date1))
a1$bayley_1y_date=as.Date(ifelse(a1$date2=="NA-NA-NA",NA,a1$date2))
a2$bayley_2y_date=as.Date(ifelse(a2$date3=="NA-NA-NA",NA,a2$date3))

m <-merge(a ,a1,by="id",all.x=T)
m2<-merge(m ,a2,by="id",all.x=T)

head(m2, 10)

m4<-m2 %>% select(id,date1,date2,date3)
m4 <- unique(m4)

sum(is.na(lee.data$yc)) 
sum(is.na(lee.data$yc1)) 
sum(is.na(lee.data$yc2))

sum(is.na(exam1$p.date))
exam1.r <- exam1[!is.na(exam1$p.date),]

M1 <- m4[!is.na(m4$bayley_6m_date),];nrow(M1)
M2 <- m4[!is.na(m4$bayley_1y_date),];nrow(M2)
M3 <- m4[!is.na(m4$bayley_2y_date),];nrow(M3)

m_id1<-merge(M1,M2,by="id",all.x=T)
m_id2<-merge(m_id1,M3,by="id",all.x=T)

m_id4  <-merge(m_id2,exam1,by.x="id",by.y="ID",all.x=T)
m_id4  <-m_id4[!is.na(m_id4$p.date),]

nrow(m_id4)

m_id4$arecode<-ifelse(m_id4$region=="EW",11,ifelse(m_id4$region=="US",26,44))

m_id4<-m_id4[!duplicated(m_id4),]
#---------------------------------------------------------------------------------------------------------#
meteo.list=NULL
for(i in 1:nrow(m_id4)){
  s<-subset(m_id4,id==m_id4$id[i])
  df<-data.frame(id=s$id,
                 date=seq(s$p.date,s$bayley_6m_date,1),
                 areacode=s$arecode,
                 sggcode =s$code)
  df$key=paste0(df$date,"-",df$areacode)

  df$T1  =ifelse(df$date<s$p2.date,1,0)
  df$T2  =ifelse(df$date>=s$p2.date & df$date<s$p3.date,1,0)
  df$T3  =ifelse(df$date>=s$p3.date & df$date<s$b.date,1,0)
  df$preg=ifelse(df$date<s$b.date,1,0) 
  
  df$m6 = ifelse(df$date>=s$b.date & df$date<s$bayley_6m_date,1,0)
  
  df2<-merge(df ,meteo,by="key",all.x=T)
  df2$SIGUNGU_DATE=with(df2,paste0(ddate,"-",sggcode))

  meteo_exposure<-c(apply(subset(df2,T1==1)   %>% select(meantemp,meanhumi),2,mean,na.rm=T),
                    apply(subset(df2,T2==1)   %>% select(meantemp,meanhumi),2,mean,na.rm=T),
                    apply(subset(df2,T3==1)   %>% select(meantemp,meanhumi),2,mean,na.rm=T),
                    apply(subset(df2,preg==1) %>% select(meantemp,meanhumi),2,mean,na.rm=T),
                    apply(subset(df2,m6==1)   %>% select(meantemp,meanhumi),2,mean,na.rm=T))
  
  
  names(meteo_exposure)=paste0(rep(c("mtemp_","mhumi_"),5),rep(c("1T","2T","3T","preg","6m"),each=2))
  s2<-cbind(s,t(as.data.frame(meteo_exposure))) %>% select(-region)
  
  meteo.list[[i]]<-s2
  print(i)
  }

meteo.df1<-do.call(rbind,meteo.list)

exam2<-merge(exam1,meteo.df1 %>% dplyr::select(id,mtemp_1T:mhumi_6m),by.x="ID",by.y="id",all.x=T)

mdu<-read_excel("..//data//questionnaire.xlsx",sheet=1) %>% dplyr::select(id, m_ma_edu)
sum(is.na(mdu))
exam2 <- merge(exam2, mdu, by.x="ID", by.y="id", all.x=T)

move1 <- read_excel("..//data//questionnaire.xlsx", sheet=1) %>% dplyr::select(id, c_resi_eq)
move2 <- read_excel("..//data//questionnaire.xlsx", sheet=2) %>% dplyr::select(id, c1_resi_eq)
move3 <- read_excel("..//data//questionnaire.xlsx", sheet=3) %>% dplyr::select(id, c2_resi_eq)
move <- merge(move1, move2, by="id", all.x=T)
move <- merge(move, move3, by="id", all.x=T)

exam2<-merge(exam2,move,by.x="ID",by.y="id",all.x=T)

write.csv(exam2,"C://Users//USER//Desktop//최종//raw_data_after_merge.csv")
