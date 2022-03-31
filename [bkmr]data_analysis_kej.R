###################################################################################################
##########################Bayesian kernel machine regression#######################################
###################################################################################################
library(bkmr);library(ggplot2)
library(dplyr);library(writexl)
library(cvTools);library(foreach)
library(gam);library(mediation)
library(readxl);library(mgcv);library(dplyr)
library(caret)

#setwd("C://Users//EJ//Desktop//ewha")
#setwd("\Users\olodi\Downloads")
setwd("C://Users//USER//Desktop//kej//univ//2020-2학기 연구//분석//나//결과//최종//추가분석")
select <- dplyr::select

load('data_formation_kej.RData')

######dat[[i]] : i=1 control, i=2 prenatal~6m, i=3 birth to 1y, i=4 birth to 2y##########
# preganacy, 6m 631명
dat <- Data[[1]]
exp <- dat %>% dplyr::select(pm10_1T.u:pm25_6m.u)

head(exp)
par(mfrow=c(1,1)); boxplot(exp)

#Correlation matrix
x11();par(mfrow=c(1,1)); corrplot::corrplot(cor(exp),method="square",type="upper",addCoef.col = TRUE,tl.col="black")

# 다른변수
cov <- dat %>% dplyr::select(mtemp_1T:mhumi_6m,d_ma_age,mbmi,d_b_wt,d_ga_wks,edu_university,shs,sex,lngreen,green234_10000)

# 건강영향
nd <- dat %>% dplyr::select(y_cogqo,y_motqo)

dat2 <- cbind(exp, cov, nd)

dat2$edu_university <- as.numeric(dat2$edu_university)
dat2$shs <- as.numeric(dat2$shs)
dat2$sex <- as.numeric(dat2$sex)
dat2$d_ma_age <- as.numeric(dat2$d_ma_age)
dat2$d_b_wt <- as.numeric(dat2$d_b_wt *0.001)
dat2$d_ga_wks <- as.numeric(dat2$d_ga_wks)
dat2$y_cogqo <- as.numeric(dat2$y_cogqo)
dat2$y_motqo <- as.numeric(dat2$y_motqo)

which(names(dat2)==c("pm25_6m.u"))
#표준화
dat2[,1:15] <- scale(dat2[,1:15], center=T, scale=T)
d_dat<-dat2

D1 <- d_dat %>% dplyr::select(y_cogqo, y_motqo, pm10_preg.u, no2_preg.u, pm25_preg.u,
                             mtemp_preg, mhumi_preg, d_ma_age, d_ga_wks, d_b_wt, mbmi, edu_university:sex,
                             lngreen, green234_10000)
D2 <- d_dat %>% dplyr::select(y_cogqo, y_motqo, pm10_1T.u, no2_1T.u, pm25_1T.u,
                             mtemp_1T, mhumi_1T, d_ma_age, d_ga_wks, d_b_wt, mbmi, edu_university.0:sex.1,
                             lngreen, green234_10000)
D3 <- d_dat %>% dplyr::select(y_cogqo, y_motqo, pm10_2T.u, no2_2T.u, pm25_2T.u,
                             mtemp_2T, mhumi_2T, d_ma_age, d_ga_wks, d_b_wt, mbmi, edu_university.0:sex.1,
                             lngreen, green234_10000)
D4 <- d_dat %>% dplyr::select(y_cogqo, y_motqo, pm10_3T.u, no2_3T.u, pm25_3T.u,
                             mtemp_3T, mhumi_3T, d_ma_age, d_ga_wks, d_b_wt, mbmi, edu_university.0:sex.1,
                             lngreen, green234_10000)
D5 <- d_dat %>% dplyr::select(y_cogqo, y_motqo, pm10_6m.u, no2_6m.u, pm25_6m.u,
                             mtemp_6m, mhumi_6m, d_ma_age, d_ga_wks, d_b_wt, mbmi, edu_university.0:sex.1,
                             lngreen, green234_10000)


# preganacy, 1y 407명
dat <- Data[[2]]

exp <- dat %>% dplyr::select(pm10_1T.u:pm25_1yc.u)

# 보정변수
cov <- dat %>% dplyr::select(mtemp_1T:mhumi_1y, d_ma_age,mbmi,d_b_wt,d_ga_wks,edu_university,shs,sex,lngreen,green234_10000)

# 건강영향
nd <- dat %>% dplyr::select(y1_cogqo,y1_motqo)

dat2 <- cbind(exp, cov, nd)

dat2$edu_university <- as.numeric(dat2$edu_university)
dat2$shs <- as.numeric(dat2$shs)
dat2$sex <- as.numeric(dat2$sex)
dat2$d_ma_age <- as.numeric(dat2$d_ma_age)
dat2$d_b_wt <- as.numeric(dat2$d_b_wt *0.001)
dat2$d_ga_wks <- as.numeric(dat2$d_ga_wks)
dat2$y1_cogqo <- as.numeric(dat2$y1_cogqo)
dat2$y1_motqo <- as.numeric(dat2$y1_motqo)

which(names(dat2)==c("pm25_1yc.u"))
#표준화
dat2[,1:20] <- scale(dat2[,1:20], center=T, scale=T)
d_dat<-dat2

D7 <- d_dat %>% dplyr::select(y1_cogqo, y1_motqo, pm10_1yc.u, no2_1yc.u, pm25_1yc.u,
                             mtemp_1y, mhumi_1y, d_ma_age, d_ga_wks, d_b_wt, mbmi, edu_university.0:sex.1, 
                             lngreen, green234_10000)
D8 <- d_dat %>% dplyr::select(y1_cogqo, y1_motqo, pm10_preg.u, no2_preg.u, 
                             mtemp_preg, mhumi_preg, d_ma_age, d_ga_wks, d_b_wt, mbmi, edu_university.0:sex.1,
                             lngreen, green234_10000)


# preganacy, 2y 258명
dat <- Data[[3]]

exp <- dat %>% dplyr::select(pm10_1T.u:pm25_2yc.u)

# 보정변수
cov <- dat %>% dplyr::select(mtemp_1T:mhumi_2y,mbmi,d_ma_age,d_b_wt,d_ga_wks,edu_university,shs,sex,lngreen,green234_10000)

# 건강영향
nd <- dat %>% dplyr::select(y2_cogqo,y2_motqo)

dat2 <- cbind(exp, cov, nd)

dat2$edu_university <- as.numeric(dat2$edu_university)
dat2$shs <- as.numeric(dat2$shs)
dat2$sex <- as.numeric(dat2$sex)
dat2$d_ma_age <- as.numeric(dat2$d_ma_age)
dat2$d_b_wt <- as.numeric(dat2$d_b_wt *0.001)
dat2$d_ga_wks <- as.numeric(dat2$d_ga_wks)
dat2$y1_cogqo <- as.numeric(dat2$y1_cogqo)
dat2$y1_motqo <- as.numeric(dat2$y1_motqo)

which(names(dat2)==c("pm25_2yc.u"))
#표준화
dat2[,1:21] <- scale(dat2[,1:21], center=T, scale=T)
d_dat<-dat2

D9 <- d_dat %>% dplyr::select(y2_cogqo, y2_motqo, pm10_2yc.u, no2_2yc.u, pm25_2yc.u,
                             mtemp_2y, mhumi_2y, d_ma_age, d_ga_wks, d_b_wt, mbmi, edu_university.0:sex.1,
                             lngreen, green234_10000)
D10 <- d_dat %>% dplyr::select(y2_cogqo, y2_motqo, pm10_preg.u, no2_preg.u, 
                             mtemp_preg, mhumi_preg, d_ma_age, d_ga_wks, d_b_wt, mbmi, edu_university.0:sex.1,
                             lngreen, green234_10000)


#-----------------------------------------------------------------------------------
panel.lm <- function(x, y, col=par("col"), bg=NA, pch=par("pch"),
                     cex=1, col.smooth="black", ...) {
  points(x, y, pch=pch, col=col, bg=bg, cex=cex)
  abline(stats::lm(y~x), col=col.smooth, ...)
} 

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
} 

panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
} 

x11();pairs(D1 %>% dplyr::select(pm10_preg.u, no2_preg.u, mtemp_preg, lngreen, y_cogqo, y_motqo), lower.panel = panel.lm,upper.panel = panel.cor, diag.panel = panel.hist)

pairs(D2 %>% dplyr::select(pm10_1T.u, no2_1T.u, mtemp_1T, lngreen, y_cogqo, y_motqo), lower.panel = panel.lm,upper.panel = panel.cor, diag.panel = panel.hist)

pairs(D3 %>% dplyr::select(pm10_2T.u, no2_2T.u, mtemp_2T, lngreen, y_cogqo, y_motqo), lower.panel = panel.lm,upper.panel = panel.cor, diag.panel = panel.hist)

pairs(D4 %>% dplyr::select(pm10_3T.u, no2_3T.u, mtemp_3T, lngreen, y_cogqo, y_motqo), lower.panel = panel.lm,upper.panel = panel.cor, diag.panel = panel.hist)


pairs(D6 %>% dplyr::select(pm10_6m.u, no2_6m.u, mtemp_6m, lngreen, y_cogqo, y_motqo), lower.panel = panel.lm,upper.panel = panel.cor, diag.panel = panel.hist)

pairs(D7 %>% dplyr::select(pm10_1yc.u, no2_1yc.u, mtemp_1y.ac, lngreen, y1_cogqo, y1_motqo), lower.panel = panel.lm,upper.panel = panel.cor, diag.panel = panel.hist)

pairs(D8 %>% dplyr::select(pm10_preg.u, no2_preg.u, mtemp_preg, lngreen, y1_cogqo, y1_motqo), lower.panel = panel.lm,upper.panel = panel.cor, diag.panel = panel.hist)

pairs(D9 %>% dplyr::select(pm10_2yc.u, no2_2yc.u, mtemp_2y, lngreen, y2_cogqo, y2_motqo), lower.panel = panel.lm,upper.panel = panel.cor, diag.panel = panel.hist)

pairs(D10 %>% dplyr::select(pm10_preg.u, no2_preg.u, mtemp_preg.ac, lngreen, y2_cogqo, y2_motqo), lower.panel = panel.lm,upper.panel = panel.cor, diag.panel = panel.hist)
#-------------------------------------------------------------------------------------------------------------------------#

# split train set:test set =1:1 >>> 1000회 반복.
# BKMR, mediation 둘다 몬테카를로 기법을 사용하기 때문에 모델의 안정화를 위해 각각 1000번의 부트스트래핑 시행
# >>> {BMKR*1000, mediation*1000}*1000회 반복

#-------------------------------------------------------------------------------------------------------------------------#
# Data <- Dx 데이터 할당
names(Data)[1:7] <- c("MDI","PDI","pm10","no2","pm25","temp","humi")
res_name <- data.frame(c("acme.est", "acme.se", "ade.est", "ade.se","tote.est", "tote.se", "prop.est", "prop.se"))
pip_name <-data.frame(c("pm10.pip", "no2.pip","pm25.pip"))

risk.s <- NULL
risk.o <- NULL
h_hat <- NULL
RES <- NULL
PIP <- NULL

for (k in 1:10){
  rn <- createDataPartition(y = Data$lngreen, p = 0.5, list = F)   # row 색인을 위해 list = F 로 설정
  
  train <- Data[rn,]            # 50%의 train data (랜덤 샘플링된 행 번호를 색인)
  test <- Data[-rn,]            # 나머지 50%의 test data
  
  # 모델 훈련
  Z_train <- train[, 3:5] # exposures
  X_train <- train[, 6:14] # covariates
  y_train <- train[, 1] # MDI
  
  # 예측
  bkmrmod <- kmbayes(y=y_train, Z=Z_train, X=X_train, iter =100, family = "gaussian", est.h = TRUE, varsel=TRUE)
  summary(bkmrmod)
  
  pip <- data.frame(summary(bkmrmod)$PIP)
  
  risks.singvar <- SingVarRiskSummaries(fit = bkmrmod, 
                                        qs.diff = c(0.25, 0.75), 
                                        q.fixed = c(0.25, 0.50, 0.75),
                                        method = "exact")
  risk.s[[k]] <- data.frame(risks.singvar)[,3:4]
  
  risks.overall <- OverallRiskSummaries(fit = bkmrmod, qs=seq(0.10, 0.90, by =0.05),
                                        q.fixed=0.5, method ="exact")
  
  risk.o[[k]] <- data.frame(risks.overall)[,2:3]
  
  # 성능 평가
  Z_test <- test[, 3:5]
  X_test <- test[, 6:14]
  y_test <- data.frame(test[, 1])
  names(y_test) <- c("MDI")
  
  h_pred <- ComputePostmeanHnew(bkmrmod, Znew = Z_test, method = "exact")
  mixture <- data.frame(h_pred$postmean)
  names(mixture) <- c("mixture")
  h_hat[[k]] <- mixture
  
  df <- cbind(y_test,mixture,X_test,test %>% select(lngreen,green234_10000))
  
  # X -> M
  model.m1<-lm(mixture~lngreen+temp+humi,data=df)
  model.m2<-lm(mixture~green234_10000+temp+humi,data=df)
  
  # X -> Y
  model.y1<-lm(MDI~lngreen+mixture+temp+humi+d_ma_age+mbmi+d_ga_wks+d_b_wt+edu_university+shs+sex,data=df)
  model.y2<-lm(MDI~green234_10000+mixture+temp+humi+d_ma_age+mbmi+d_ga_wks+d_b_wt+edu_university+shs+sex,data=df)
  
  # causal mediation
  outcome1 <- mediate(model.m1, model.y1, sims=100, boot=TRUE, treat="lngreen", mediator="mixture")
  outcome2 <- mediate(model.m2, model.y2, sims=100, boot=TRUE, treat="green234_10000", mediator="mixture")
  
  #summary(outcome)
  cme <- data.frame(summary(outcome1)$d.avg)
  cme.se <- data.frame((outcome1$d.avg.ci[2]-outcome1$d.avg.ci[1])/3.92)
  de <- data.frame(summary(outcome1)$z.avg)
  de.se <- data.frame((outcome1$z.avg.ci[2]-outcome1$z.avg.ci[1])/3.92)
  tt <- data.frame(summary(outcome1)$tau.coef)
  tt.se <- data.frame((outcome1$tau.ci[2]-outcome1$tau.ci[1])/3.92)
  pp <- data.frame(summary(outcome1)$n.avg)
  pp.se <- data.frame((outcome1$n.avg.ci[2]-outcome1$n.avg.ci[1])/3.92)
  
  rs1 <- t(data.frame(cbind(cme, cme.se, de, de.se, tt, tt.se, pp, pp.se)))
  
  cme <- data.frame(summary(outcome2)$d.avg)
  cme.se <- data.frame((outcome2$d.avg.ci[2]-outcome2$d.avg.ci[1])/3.92)
  de <- data.frame(summary(outcome2)$z.avg)
  de.se <- data.frame((outcome2$z.avg.ci[2]-outcome2$z.avg.ci[1])/3.92)
  tt <- data.frame(summary(outcome2)$tau.coef)
  tt.se <- data.frame((outcome2$tau.ci[2]-outcome2$tau.ci[1])/3.92)
  pp <- data.frame(summary(outcome2)$n.avg)
  pp.se <- data.frame((outcome2$n.avg.ci[2]-outcome2$n.avg.ci[1])/3.92)
  
  rs2 <- t(data.frame(cbind(cme, cme.se, de, de.se, tt, tt.se, pp, pp.se)))
  
  pm10.pip <- pip[1,]
  no2.pip <- pip[2,]
  pm25.pip <- pip[3,]
  
  pip <- t(data.frame(cbind(pm10.pip, no2.pip, pm25.pip)))
  res <- data.frame(cbind(rs1, rs2))
  RES[[k]] <- res 
  PIP[[k]] <- pip
  
  print(k)}

######################
Data <- D1

ov_hat <- do.call(cbind, h_hat);risks.s <- do.call(cbind, risk.s);risks.o <- do.call(cbind, risk.o)
write.csv(result, "mdi1T_1000times_bkmr,cme(SE)_results.xlsx");
write.csv(ov_hat,"Mpredicted_oveff1.csv");write.csv(risks.s,"Msingle_oveff1.csv");write.csv(risks.o,"MOver_oveff1.csv")
######################
Data <- D2

ov_hat <- do.call(cbind, h_hat);risks.s <- do.call(cbind, risk.s);risks.o <- do.call(cbind, risk.o)
write.csv(result, "mdi2T_1000times_bkmr,cme(SE)_results.xlsx")
write.csv(ov_hat,"Mpredicted_oveff2.csv");write.csv(risks.s,"Msingle_oveff2.csv");write.csv(risks.o,"MOver_oveff2.csv")
######################
Data <- D3

ov_hat <- do.call(cbind, h_hat);risks.s <- do.call(cbind, risk.s);risks.o <- do.call(cbind, risk.o)
write.csv(result, "mdi3T_1000times_bkmr,cme(SE)_results.xlsx")
write.csv(ov_hat,"Mpredicted_oveff3.csv");write.csv(risks.s,"Msingle_oveff3.csv");write.csv(risks.o,"MOver_oveff3.csv")
######################
Data <- D4

ov_hat <- do.call(cbind, h_hat);risks.s <- do.call(cbind, risk.s);risks.o <- do.call(cbind, risk.o)
write.csv(result, "mdipreg_1000times_bkmr,cme(SE)_results.xlsx")
write.csv(ov_hat,"Mpredicted_oveff4.csv");write.csv(risks.s,"Msingle_oveff4.csv");write.csv(risks.o,"MOver_oveff4.csv")
######################
Data <- D5

ov_hat <- do.call(cbind, h_hat);risks.s <- do.call(cbind, risk.s);risks.o <- do.call(cbind, risk.o)
write.csv(result, "mdi6m.ac_1000times_bkmr,cme(SE)_results.xlsx")
write.csv(ov_hat,"Mpredicted_oveff5.csv");write.csv(risks.s,"Msingle_oveff5.csv");write.csv(risks.o,"MOver_oveff5.csv")
######################
Data <- D6

ov_hat <- do.call(cbind, h_hat);risks.s <- do.call(cbind, risk.s);risks.o <- do.call(cbind, risk.o)
write.csv(result, "mdi6m_1000times_bkmr,cme(SE)_results.xlsx")
write.csv(ov_hat,"Mpredicted_oveff6.csv");write.csv(risks.s,"Msingle_oveff6.csv");write.csv(risks.o,"MOver_oveff6.csv")
######################
Data <- D7

ov_hat <- do.call(cbind, h_hat);risks.s <- do.call(cbind, risk.s);risks.o <- do.call(cbind, risk.o)
write.csv(result, "mdi1y.ac_1000times_bkmr,cme(SE)_results.xlsx")
write.csv(ov_hat,"Mpredicted_oveff7.csv");write.csv(risks.s,"Msingle_oveff7.csv");write.csv(risks.o,"MOver_oveff7.csv")

######################
Data <- D8

ov_hat <- do.call(cbind, h_hat);risks.s <- do.call(cbind, risk.s);risks.o <- do.call(cbind, risk.o)
write.csv(result, "mdi1y_1000times_bkmr,cme(SE)_results.xlsx")
write.csv(ov_hat,"Mpredicted_oveff8.csv");write.csv(risks.s,"Msingle_oveff8.csv");write.csv(risks.o,"MOver_oveff8.csv")

######################
Data <- D9

ov_hat <- do.call(cbind, h_hat);risks.s <- do.call(cbind, risk.s);risks.o <- do.call(cbind, risk.o)
write.csv(result, "mdi2y.ac_1000times_bkmr,cme(SE)_results.xlsx")
write.csv(ov_hat,"Mpredicted_oveff9.csv");write.csv(risks.s,"Msingle_oveff9.csv");write.csv(risks.o,"MOver_oveff9.csv")

######################
Data <- D10

ov_hat <- do.call(cbind, h_hat);risks.s <- do.call(cbind, risk.s);risks.o <- do.call(cbind, risk.o)
write.csv(result, "mdi2y_1000times_bkmr,cme(SE)_results.xlsx")
write.csv(ov_hat,"Mpredicted_oveff10.csv");write.csv(risks.s,"Msingle_oveff10.csv");write.csv(risks.o,"MOver_oveff10.csv")

#-------------------------------------------------------------------------------------------------------------------------#
library(caret)
set.seed(2021)

names(Data)[1:6] <- c("MDI","PDI","pm10","no2","temp","humi")
result <- data.frame(c("acme.est", "acme.se", "ade.est", "ade.se","tote.est", "tote.se", "prop.est", "prop.se", "pm10.pip", "no2.pip","temp.pip"))
risk.s <- NULL
risk.o <- NULL
h_hat <- NULL

for (k in 1:1000){
  rn <- createDataPartition(y = Data$lngreen, p = 0.5, list = F)   # row 색인을 위해 list = F 로 설정
  
  train <- Data[rn,]            # 50%의 train data (랜덤 샘플링된 행 번호를 색인)
  test <- Data[-rn,]            # 나머지 50%의 test data
  
  # 모델 훈련
  Z_train <- train[, 3:5]
  X_train <- train[, 6:14]
  y_train <- train[, 2] # PDI
  
  # 예측
  bkmrmod <- kmbayes(y=y_train, Z=Z_train, X=X_train, iter =1000, family = "gaussian", est.h = TRUE, varsel=TRUE)
  #summary(bkmrmod)
  
  #h_hat <- data.frame(apply(bkmrmod$h.hat,2,mean))
  pip <- data.frame(summary(bkmrmod)$PIP)
  
  risks.singvar <- SingVarRiskSummaries(fit = bkmrmod, 
                                        qs.diff = c(0.25, 0.75), 
                                        q.fixed = c(0.25, 0.50, 0.75),
                                        method = "exact")
  risk.s[[k]] <- data.frame(risks.singvar)[,3:4]
  
  risks.overall <- OverallRiskSummaries(fit = bkmrmod, qs=seq(0.10, 0.90, by =0.05),
                                        q.fixed=0.5, method ="exact")
  
  risk.o[[k]] <- data.frame(risks.overall)[,2:3]
  
  # 성능 평가      
  Z_test <- test[, 3:5]
  X_test <- test[, 6:14]
  y_test <- data.frame(test[,2])
  names(y_test) <- c("PDI")
  
  h_pred <- ComputePostmeanHnew(bkmrmod, Znew = Z_test, method = "exact")
  mixture <- data.frame(h_pred$postmean)
  names(mixture) <- c("mixture")
  h_hat[[k]] <- mixture
  
  df <- cbind(y_test,mixture,X_test)
  
  # X -> M
  model.m<-gam(mixture~lngreen+humi,data=df)
  
  # X -> Y
  model.y<-gam(PDI~lngreen+mixture+d_ma_age+mbmi+d_ga_wks+d_b_wt+m_ma_edu+shs+sex+humi,data=df)
  
  # causal mediation
  med <- c("mixture")
  outcome <- mediate(model.m, model.y, sims=1000, boot=TRUE, treat="lngreen", mediator=med, treat.value=5, control.value=3)
  
  #summary(outcome)
  cme <- data.frame(summary(outcome)$d.avg)
  cme.se <- data.frame((outcome$d.avg.ci[2]-outcome$d.avg.ci[1])/3.92)
  de <- data.frame(summary(outcome)$z.avg)
  de.se <- data.frame((outcome$z.avg.ci[2]-outcome$z.avg.ci[1])/3.92)
  tt <- data.frame(summary(outcome)$tau.coef)
  tt.se <- data.frame((outcome$tau.ci[2]-outcome$tau.ci[1])/3.92)
  pp <- data.frame(summary(outcome)$n.avg)
  pp.se <- data.frame((outcome$n.avg.ci[2]-outcome$n.avg.ci[1])/3.92)
  pm10.pip <- pip[1,]
  no2.pip <- pip[2,]
  temp.pip <- pip[3,]
  
  rs <- t(data.frame(cbind(cme, cme.se, de, de.se, tt, tt.se, pp, pp.se, pm10.pip, no2.pip, temp.pip)))
  result <- cbind(result, rs)
  
  print(k)}


######################
Data <- D1

ov_hat <- do.call(cbind, h_hat);risks.s <- do.call(cbind, risk.s);risks.o <- do.call(cbind, risk.o)
write.csv(result, "pdi1T_1000times_bkmr,cme(SE)_results.xlsx");
write.csv(ov_hat,"Ppredicted_oveff1.csv");write.csv(risks.s,"Psingle_oveff1.csv");write.csv(risks.o,"POver_oveff1.csv")
######################
Data <- D2

ov_hat <- do.call(cbind, h_hat);risks.s <- do.call(cbind, risk.s);risks.o <- do.call(cbind, risk.o)
write.csv(result, "pdi2T_1000times_bkmr,cme(SE)_results.xlsx");
write.csv(ov_hat,"Ppredicted_oveff2.csv");write.csv(risks.s,"Psingle_oveff2.csv");write.csv(risks.o,"POver_oveff2.csv")
######################
Data <- D3

ov_hat <- do.call(cbind, h_hat);risks.s <- do.call(cbind, risk.s);risks.o <- do.call(cbind, risk.o)
write.csv(result, "pdi3T_1000times_bkmr,cme(SE)_results.xlsx");
write.csv(ov_hat,"Ppredicted_oveff3.csv");write.csv(risks.s,"Psingle_oveff3.csv");write.csv(risks.o,"POver_oveff3.csv")
######################
Data <- D4

ov_hat <- do.call(cbind, h_hat);risks.s <- do.call(cbind, risk.s);risks.o <- do.call(cbind, risk.o)
write.csv(result, "pdipreg_1000times_bkmr,cme(SE)_results.xlsx");
write.csv(ov_hat,"Ppredicted_oveff4.csv");write.csv(risks.s,"Psingle_oveff4.csv");write.csv(risks.o,"POver_oveff4.csv")
######################
Data <- D5

ov_hat <- do.call(cbind, h_hat);risks.s <- do.call(cbind, risk.s);risks.o <- do.call(cbind, risk.o)
write.csv(result, "pdi6m.ac_1000times_bkmr,cme(SE)_results.xlsx");
write.csv(ov_hat,"Ppredicted_oveff5.csv");write.csv(risks.s,"Psingle_oveff5.csv");write.csv(risks.o,"POver_oveff5.csv")
######################
Data <- D6

ov_hat <- do.call(cbind, h_hat);risks.s <- do.call(cbind, risk.s);risks.o <- do.call(cbind, risk.o)
write.csv(result, "pdi6m_1000times_bkmr,cme(SE)_results.xlsx");
write.csv(ov_hat,"Ppredicted_oveff6.csv");write.csv(risks.s,"Psingle_oveff6.csv");write.csv(risks.o,"POver_oveff6.csv")
######################
Data <- D7

ov_hat <- do.call(cbind, h_hat);risks.s <- do.call(cbind, risk.s);risks.o <- do.call(cbind, risk.o)
write.csv(result, "pdi1y.ac_1000times_bkmr,cme(SE)_results.xlsx");
write.csv(ov_hat,"Ppredicted_oveff7.csv");write.csv(risks.s,"Psingle_oveff7.csv");write.csv(risks.o,"POver_oveff7.csv")
######################
Data <- D8

ov_hat <- do.call(cbind, h_hat);risks.s <- do.call(cbind, risk.s);risks.o <- do.call(cbind, risk.o)
write.csv(result, "pdi1y_1000times_bkmr,cme(SE)_results.xlsx");
write.csv(ov_hat,"Ppredicted_oveff8.csv");write.csv(risks.s,"Psingle_oveff8.csv");write.csv(risks.o,"POver_oveff8.csv")
######################
Data <- D9

ov_hat <- do.call(cbind, h_hat);risks.s <- do.call(cbind, risk.s);risks.o <- do.call(cbind, risk.o)
write.csv(result, "pdi2y.ac_1000times_bkmr,cme(SE)_results.xlsx");
write.csv(ov_hat,"Ppredicted_oveff9.csv");write.csv(risks.s,"Psingle_oveff9.csv");write.csv(risks.o,"POver_oveff9.csv")
######################
Data <- D10

ov_hat <- do.call(cbind, h_hat);risks.s <- do.call(cbind, risk.s);risks.o <- do.call(cbind, risk.o)
write.csv(result, "pdi2y_1000times_bkmr,cme(SE)_results.xlsx");
write.csv(ov_hat,"Ppredicted_oveff10.csv");write.csv(risks.s,"Psingle_oveff10.csv");write.csv(risks.o,"POver_oveff10.csv")
