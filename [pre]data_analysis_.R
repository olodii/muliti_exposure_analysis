###################################################################################################
###################################### BKMR -> CME ################################################
###################################################################################################
library(bkmr);library(ggplot2)
library(dplyr);library(writexl)
library(cvTools);library(foreach)
library(gam);library(mediation)
library(readxl);library(mgcv);library(dplyr)
library(caret)

setwd("C://Users//USER//Desktop//최종")
select <- dplyr::select

load('data_formation_fin.RData')
#-------------------------------------------------------------------------------------------------------------------------#

# split train set:test set =1:1 >>> 1000회 반복.
# BKMR, mediation 둘다 몬테카를로 기법을 사용하기 때문에 모델의 안정화를 위해 각각 1000번의 부트스트래핑 시행
# >>> {BMKR*1000, mediation*1000}*1000회 반복

#-------------------------------------------------------------------------------------------------------------------------#
# Prenatal period: 631명 D1~D5
pre <- NULL
pre[[1]] <- D1
pre[[2]] <- D2
pre[[3]] <- D3
pre[[4]] <- D4
pre[[5]] <- D5
#-------------------------------------------------------------------------------------------------------------------------#
# Data <- Dx 데이터 할당
res_name <- data.frame(c("acme.est", "acme.se", "ade.est", "ade.se","tote.est", "tote.se", "prop.est", "prop.se"))
pip_name <-data.frame(c("pm10.pip", "no2.pip","pm25.pip"))

risk.s <- NULL
risk.o <- NULL
h_hat <- NULL
cme <- NULL
pip <- NULL

all.risk.s <- NULL
all.risk.o <- NULL
all.pip <- NULL
all.pred <- NULL
all.cme <- NULL
all.cme2 <- NULL

for (m in 1:length(pre)){
  # prenatal period iterate
  data <- pre[[m]]
  period <- unlist(strsplit(names(data)[6], split='_'))[2]
  
  # colnames는 동일하게 할당
  names(data)[1:7] <- c("MDI","PDI","pm10","no2","pm25","temp","humi")
  
  for (k in 1:1000){
    
    # split train/test -> 이 과정이 1000번 반복되면서 다른 결과 산출
    rn <- createDataPartition(y = data$lngreen, p = 0.5, list = F)   # row 색인을 위해 list = F 로 설정
    
    train <- data[rn,] # 50%의 train data (랜덤 샘플링된 행 번호를 색인)
    test <- data[-rn,] # 나머지 50%의 test data
    
    #동일한 exp, cov, out 개수와 순서
    Z_train <- train[, 3:5] # exposures
    X_train <- train[, 6:14] # covariates
    y_train <- train[, 1] # MDI
    
    # 1. BKMR
    # training
    bkmrmod <- kmbayes(y=y_train, Z=Z_train, X=X_train, iter =1000, family = "gaussian",
                       est.h = TRUE, varsel=TRUE, verbose = FALSE)
    
    # bkmr result: feature importance
    ps <- data.frame(PIPs=ExtractPIPs(bkmrmod)[,2])
    rownames(ps) <- c('pm10','no2','pm25')
    
    # bkmr result: single exposure interaction
    risk.sing <- SingVarRiskSummaries(fit = bkmrmod, 
                                          qs.diff = c(0.25, 0.75), 
                                          q.fixed = c(0.25, 0.50, 0.75),
                                          method = "exact")
    risk.sing <- data.frame(risk.sing)[,3:4]

    risk.s[[k]] <- risk.sing
    
    # bkmr result: single exposure interaction
    risk.over <- OverallRiskSummaries(fit = bkmrmod, qs=seq(0.10, 0.90, by =0.05),
                                          q.fixed=0.5, method ="exact")
    rownames(risk.over) <- risk.over$quantile
    risk.over <- data.frame(risk.over)[,2:3] 
    risk.o[[k]] <- risk.over
    
    # prediction of personal exposure level
    Z_test <- test[, 3:5]
    X_test <- test[, 6:14]
    y_test <- data.frame(test[, 1])
    names(y_test) <- c("MDI")
    
    h_pred <- ComputePostmeanHnew(bkmrmod, Znew = Z_test, method = "exact")
    mixture <- data.frame(h_pred$postmean)
    names(mixture) <- c("mixture")
    h_hat[[k]] <- mixture
    
    #2. Causal mediation
    df <- cbind(y_test,mixture,X_test,test %>% select(lngreen,green234_10000))
    
    # X -> M
    model.m1<-lm(mixture~lngreen+temp+humi,data=df)
    model.m2<-lm(mixture~green234_10000+temp+humi,data=df)
    
    # X -> Y
    model.y1<-lm(MDI~lngreen+mixture+temp+humi+d_ma_age+mbmi+d_ga_wks+d_b_wt+edu_university+shs+sex,data=df)
    model.y2<-lm(MDI~green234_10000+mixture+temp+humi+d_ma_age+mbmi+d_ga_wks+d_b_wt+edu_university+shs+sex,data=df)
    
    # causal mediation
    outcome1 <- mediate(model.m1, model.y1, sims=1000, boot=TRUE, treat="lngreen", mediator="mixture")
    outcome2 <- mediate(model.m2, model.y2, sims=1000, boot=TRUE, treat="green234_10000", mediator="mixture")
    
    #summary(outcome)
    rs1 <- t(data.frame(summary(outcome1)$d.avg,
                        (outcome1$d.avg.ci[2]-outcome1$d.avg.ci[1])/3.92,
                        summary(outcome1)$z.avg,
                        (outcome1$z.avg.ci[2]-outcome1$z.avg.ci[1])/3.92,
                        summary(outcome1)$tau.coef,
                        (outcome1$tau.ci[2]-outcome1$tau.ci[1])/3.92,
                        summary(outcome1)$n.avg,
                        (outcome1$n.avg.ci[2]-outcome1$n.avg.ci[1])/3.92))
    rownames(rs1) <- c('cme','cme.se','de','de.se','tt','tt.se','pp','pp.se')

    
    rs2 <- t(data.frame(summary(outcome2)$d.avg,
                        (outcome2$d.avg.ci[2]-outcome2$d.avg.ci[1])/3.92,
                        summary(outcome2)$z.avg,
                        (outcome2$z.avg.ci[2]-outcome2$z.avg.ci[1])/3.92,
                        summary(outcome2)$tau.coef,
                        (outcome2$tau.ci[2]-outcome2$tau.ci[1])/3.92,
                        summary(outcome2)$n.avg,
                        (outcome2$n.avg.ci[2]-outcome2$n.avg.ci[1])/3.92))
    rownames(rs2) <- c('cme','cme.se','de','de.se','tt','tt.se','pp','pp.se')
    
    res <- data.frame(cbind(rs1, rs2))
    colnames(res) <- c('lng200m','g10km')
    
    cme[[k]] <- res 
    pip[[k]] <- ps
    
    print(paste0('-----In ',period,': The ',k,'th iteration is done-------'))
    
    }
  
  risks.s <- do.call(cbind, risk.s) 
  risks.o <- do.call(cbind, risk.o)
  pips <- do.call(cbind, pip)
  preds <- do.call(cbind, h_hat)
  cmes <- do.call(cbind, cme)
  
  #### mean, 95% CI for 1000 iteration
  # average & 95% percentile of Risk ratio b/w SingleVariables
  avg.risks.s <- data.frame(round(apply(risks.s[,colnames(risks.s)=='est'],1,mean),3),
                            round(apply(risks.s[,colnames(risks.s)=='est'],1,quantile, probs=0.025),3),
                            round(apply(risks.s[,colnames(risks.s)=='est'],1,quantile, probs=0.975),3))
  
  avg.risks.s$sig <- ifelse(avg.risks.s[,2]<0 & avg.risks.s[,3]>0, '', '*') 
  
  avg.risks.s$full <- paste0(avg.risks.s[,1], ' (',
                             avg.risks.s[,2], ', ',
                             avg.risks.s[,3], ')')
  
  colnames(avg.risks.s) <- c(paste('rsksing',period,'avg',sep='_'),
                             paste('rsksing',period,'lci',sep='_'),
                             paste('rsksing',period,'uci',sep='_'),
                             'sig',
                             paste('rsksing',period, sep='_'))
                             
  rownames(avg.risks.s) <- c('0.25_pm10','0.25_no2','0.25_pm25',
                             '0.50_pm10','0.50_no2','0.50_pm25',
                             '0.75_pm10','0.75_no2','0.75_pm25')
  
  print(paste('In SingVarRiskSummaries, there is',
              length(subset(avg.risks.s, sig =='*')),
              'significant association in the results.', sep=" "))
  
  #print(avg.risks.s[,1:4])
  
  # average & 95% percentile of Risk ratio of overall effect
  avg.risks.o <- data.frame(round(apply(risks.o[,colnames(risks.o)=='est'],1,mean),3),
                            round(apply(risks.o[,colnames(risks.o)=='est'],1,quantile, probs=0.025),3),
                            round(apply(risks.o[,colnames(risks.o)=='est'],1,quantile, probs=0.975),3))
  
  avg.risks.o$sig <- ifelse(avg.risks.o[,2]<0 & avg.risks.o[,3]>0, '', '*')
  
  avg.risks.o$full <- paste0(avg.risks.o[,1], ' (',
                             avg.risks.o[,2], ', ',
                             avg.risks.o[,3], ')')
  
  colnames(avg.risks.o) <- c(paste('rskover',period,'avg',sep='_'),
                             paste('rskover',period,'lci',sep='_'),
                             paste('rskover',period,'uci',sep='_'),
                             'sig',
                             paste('rskover',period, sep='_'))
  
  rownames(avg.risks.o) <- seq(from=0.10, to=0.90, by=0.05)
  
  print(paste('In OverallRiskSummaries, there is',
              length(subset(avg.risks.o, sig =='*')),
              'significant association in the results.', sep=" "))
  
  print(avg.risks.o[,1:4])
  
  # average & 95% percentile of PIP
  av.pp <- data.frame(avg=apply(pips,2,mean),
                      std=apply(pips,2,sd))
  avg.pips <- data.frame(round(apply(av.pp,2,mean),3))
  rownames(avg.pips) <- c('pips_avg','pips_std')
  colnames(avg.pips) <- paste('pre',period, sep='_')
  
  # average of personal exposure levels
  av.prd <- data.frame(avg=apply(preds,2,mean),
                       std=apply(preds,2,sd))
  avg.preds <- data.frame(round(apply(av.prd,2,mean),3))
  rownames(avg.preds) <- c('mix_avg','mix_std')
  colnames(avg.preds) <- paste('pre',period, sep='_')
  
  ## average of causal mediation results
  # lngreen 200m buffer
  avg.cme <- data.frame(round(apply(cmes[,colnames(cmes)=='lng200m'],1,mean),3),
                        round(apply(cmes[,colnames(cmes)=='lng200m'],1,quantile, probs=0.025),3),
                        round(apply(cmes[,colnames(cmes)=='lng200m'],1,quantile, probs=0.975),3))
  
  avg.cme$sig <- ifelse(avg.cme[,2]<0 & avg.cme[,3]>0, '', '*')
  
  avg.cme$full <- paste0(avg.cme[,1], ' (',
                         avg.cme[,2], ', ',
                         avg.cme[,3], ')')
  
  colnames(avg.cme) <- c(paste('cme',period,'avg',sep='_'),
                         paste('cme',period,'lci',sep='_'),
                         paste('cme',period,'uci',sep='_'),
                         'sig',
                         paste('cme',period, sep='_'))
 
  message <- paste('lngreen 200m: There is significant mediation effect in', paste('pre',period,sep='_'), sep=' ')
  message2 <-paste('lngreen 200m: There is no mediation effect in', paste('pre',period,sep='_'), sep=' ')
  ifelse((avg.cme['tt','sig'] == '*') & (avg.cme['cme','sig'] == '*'), message, message2)

  # green234 10000m buffer
  avg.cme2 <- data.frame(round(apply(cmes[,colnames(cmes)=='g10km'],1,mean),3),
                         round(apply(cmes[,colnames(cmes)=='g10km'],1,quantile, probs=0.025),3),
                         round(apply(cmes[,colnames(cmes)=='g10km'],1,quantile, probs=0.975),3))
  
  avg.cme2$sig <- ifelse(avg.cme2[,2]<0 & avg.cme2[,3]>0, '', '*')
  
  avg.cme2$full <- paste0(avg.cme2[,1], ' (',
                          avg.cme2[,2], ', ',
                          avg.cme2[,3], ')')
  
  colnames(avg.cme2) <- c(paste('cme',period,'avg',sep='_'),
                          paste('cme',period,'lci',sep='_'),
                          paste('cme',period,'uci',sep='_'),
                          'sig',
                          paste('cme',period,sep='_'))
  
  message <- paste('green 10km: There is significant mediation effect in', paste('pre',period,sep='_'), sep=' ')
  message2 <-paste('green 10km: There is no mediation effect in', paste('pre',period,sep='_'), sep=' ')
  ifelse((avg.cme2['tt','sig'] == '*') & (avg.cme2['cme','sig'] == '*'), message, message2)
  
  # final results
  all.risk.s[[m]] <- avg.risks.s[,5]
  all.risk.o[[m]] <- avg.risks.o[,5]
  all.pip[[m]] <- avg.pips[,1]
  all.pred[[m]] <- avg.preds[,1]
  all.cme[[m]] <- avg.cme[,5]
  all.cme2[[m]] <- avg.cme2[,5]

  }

all.RISK.s <- do.call(cbind, all.risk.s)
rownames(all.RISK.s) <- c('0.25_pm10','0.25_no2','0.25_pm25',
                          '0.50_pm10','0.50_no2','0.50_pm25',
                          '0.75_pm10','0.75_no2','0.75_pm25')

all.RISK.o <- do.call(cbind, all.risk.o)
rownames(all.RISK.o) <- seq(from=0.10, to=0.90, by=0.05)

all.PIP <- do.call(cbind, all.pip)
rownames(all.PIP) <- c('pip_avg','pip_std')

all.PRED <- do.call(cbind, all.pred)
rownames(all.PRED) <- c('mix_avg','mix_std')

all.CME <- do.call(cbind, all.cme)
rownames(all.CME) <- c('cme','cme.se','de','de.se','tt','tt.se','pp','pp.se')




#######################R E S U L T#########################################
# MDI, prenatal period
write.csv(all.RISK.s, "MDI_PRE_bkmr_single.csv", row.names=TRUE)
write.csv(all.RISK.o,"MDI_PRE_bkmr_overall.csv", row.names=TRUE)
write.csv(all.PIP,"MDI_PRE_bkmr_pip.csv", row.names=TRUE)
write.csv(all.PRED,"MDI_PRE_bkmr_mixture.csv", row.names=TRUE)
write.csv(all.CME,"MDI_PRE_causal_mediation.csv", row.names=TRUE)

# PDI, prenatal period
write.csv(all.RISK.s, "PDI_PRE_bkmr_single.csv", row.names=TRUE)
write.csv(all.RISK.o,"PDI_PRE_bkmr_overall.csv", row.names=TRUE)
write.csv(all.PIP,"PDI_PRE_bkmr_pip.csv", row.names=TRUE)
write.csv(all.PRED,"PDI_PRE_bkmr_mixture.csv", row.names=TRUE)
write.csv(all.CME,"PDI_PRE_causal_mediation.csv", row.names=TRUE)

