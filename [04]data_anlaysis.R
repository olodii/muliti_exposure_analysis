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

load('data_formation.RData')
#-------------------------------------------------------------------------------------------------------------------------#

# split train set:test set =1:1 >>> 1000회 반복.
# BKMR, mediation 둘다 몬테카를로 기법을 사용하기 때문에 모델의 안정화를 위해 각각 1000번의 부트스트래핑 시행
# >>> {BMKR*1000, mediation*1000}*1000회 반복

#-------------------------------------------------------------------------------------------------------------------------#
pre <- NULL
pre[[1]] <- D1
pre[[2]] <- D2
pre[[3]] <- D3
pre[[4]] <- D4
pre[[5]] <- D5
#-------------------------------------------------------------------------------------------------------------------------#
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
all.cme3 <- NULL
all.cme4 <- NULL

for (m in 1:length(pre)){
  # prenatal period iterate
  data <- pre[[m]]
  period <- unlist(strsplit(names(data)[6], split='_'))[2]
  
  # colnames는 동일하게 할당
  names(data)[1:7] <- c("MDI","PDI","pm10","no2","pm25","temp","humi")
  
  for (k in 1:1000){
    
    rn <- createDataPartition(y = data$green, p = 0.5, list = F)
    
    train <- data[rn,] # 50%의 train data
    test <- data[-rn,] # 나머지 50% test data
    
    #동일한 exp, cov, out 개수와 순서
    Z_train <- train[, 3:7] # exposures
    X_train <- train[, 8:14] # covariates
    y_train <- train[, 1] 
    
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
    Z_test <- test[, 3:7]
    X_test <- test[, 8:14]
    y_test <- data.frame(test[, 1])
    names(y_test) <- c("MDI")
    
    h_pred <- ComputePostmeanHnew(bkmrmod, Znew = Z_test, method = "exact")
    mixture <- data.frame(h_pred$postmean)
    names(mixture) <- c("mixture")
    h_hat[[k]] <- mixture
    
    #2. Causal mediation
    df <- cbind(y_test,mixture,X_test,test %>% select(temp,humi,green,green2))
    
    # X -> M
    model.m1<-lm(mixture~lng100m+temp+humi,data=df)
    model.m2<-lm(mixture~lng200m+temp+humi,data=df)
    model.m3<-lm(mixture~lng300m+temp+humi,data=df)
    model.m4<-lm(mixture~lng500m+temp+humi,data=df)
    
    # X -> Y
    model.y1<-lm(MDI~lng100m+mixture+temp+humi+age+mbi+wks+dbwt+medu+shs+sex,data=df)
    model.y2<-lm(MDI~lng200m+mixture+temp+humi+age+mbi+wks+dbwt+medu+shs+sex,data=df)
    model.y3<-lm(MDI~lng300m+mixture+temp+humi+age+mbi+wks+dbwt+medu+shs+sex,data=df)
    model.y4<-lm(MDI~lng500m+mixture+temp+humi+age+mbi+wks+dbwt+medu+shs+sex,data=df)
    
    # causal mediation
    outcome1 <- mediate(model.m1, model.y1, sims=1000, boot=TRUE, treat="lng100m", mediator="mixture")
    outcome2 <- mediate(model.m2, model.y2, sims=1000, boot=TRUE, treat="lng200m", mediator="mixture")
    outcome3 <- mediate(model.m2, model.y2, sims=1000, boot=TRUE, treat="lng300m", mediator="mixture")
    outcome4 <- mediate(model.m2, model.y2, sims=1000, boot=TRUE, treat="lng500m", mediator="mixture")
    
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
    
    rs3 <- t(data.frame(summary(outcome3)$d.avg,
                        (outcome3$d.avg.ci[2]-outcome3$d.avg.ci[1])/3.92,
                        summary(outcome3)$z.avg,
                        (outcome3$z.avg.ci[2]-outcome3$z.avg.ci[1])/3.92,
                        summary(outcome3)$tau.coef,
                        (outcome3$tau.ci[2]-outcome3$tau.ci[1])/3.92,
                        summary(outcome3)$n.avg,
                        (outcome3$n.avg.ci[2]-outcome3$n.avg.ci[1])/3.92))
    rownames(rs3) <- c('cme','cme.se','de','de.se','tt','tt.se','pp','pp.se')
    
    rs4 <- t(data.frame(summary(outcome4)$d.avg,
                        (outcome4$d.avg.ci[2]-outcome4$d.avg.ci[1])/3.92,
                        summary(outcome4)$z.avg,
                        (outcome1$z.avg.ci[2]-outcome4$z.avg.ci[1])/3.92,
                        summary(outcome4)$tau.coef,
                        (outcome4$tau.ci[2]-outcome4$tau.ci[1])/3.92,
                        summary(outcome4)$n.avg,
                        (outcome4$n.avg.ci[2]-outcome4$n.avg.ci[1])/3.92))
    rownames(rs4) <- c('cme','cme.se','de','de.se','tt','tt.se','pp','pp.se')
    
    res <- data.frame(cbind(rs1, rs2,rs3,rs4))
    colnames(res) <- c('100m','200m','300m','500m')
    
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
                             
  rownames(avg.risks.s) <- c('0.25_pm10','0.25_no2','0.25_pm25','0.25_temp','0.25_humi',
                             '0.50_pm10','0.50_no2','0.50_pm25','0.50_temp','0.50_humi',
                             '0.75_pm10','0.75_no2','0.75_pm25','0.75_temp','0.75_humi')
  
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
  avg.pips <- t(data.frame(pips_avg=round(apply(pips,1,mean),3)))
  rownames(avg.pips) <- c('pips_avg')
  colnames(avg.pips) <- paste('pre',period, colnames(avg.pips), sep='_')
  
  # average of personal exposure levels
  av.prd <- data.frame(avg=apply(preds,2,mean),
                       std=apply(preds,2,sd))
  avg.preds <- data.frame(round(apply(av.prd,2,mean),3))
  rownames(avg.preds) <- c('mix_avg','mix_std')
  colnames(avg.preds) <- paste('pre',period, sep='_')
  
  ## average of causal mediation results
  # green 100m buffer
  avg.cme <- data.frame(round(apply(cmes[,colnames(cmes)=='lng100m'],1,mean),3),
                        round(apply(cmes[,colnames(cmes)=='lng100m'],1,quantile, probs=0.025),3),
                        round(apply(cmes[,colnames(cmes)=='lng100m'],1,quantile, probs=0.975),3))
  
  avg.cme$sig <- ifelse(avg.cme[,2]<0 & avg.cme[,3]>0, '', '*')
  
  avg.cme$full <- paste0(avg.cme[,1], ' (',
                         avg.cme[,2], ', ',
                         avg.cme[,3], ')')
  
  colnames(avg.cme) <- c(paste('cme',period,'avg',sep='_'),
                         paste('cme',period,'lci',sep='_'),
                         paste('cme',period,'uci',sep='_'),
                         'sig',
                         paste('cme',period, sep='_'))
  
  message <- paste('green 100m: There is significant mediation effect in', paste('pre',period,sep='_'), sep=' ')
  message2 <-paste('green 100m: There is no mediation effect in', paste('pre',period,sep='_'), sep=' ')
  ifelse((avg.cme['tt','sig'] == '*') & (avg.cme['cme','sig'] == '*'), message, message2)
  
  # lngreen 200m buffer
  avg.cme2 <- data.frame(round(apply(cmes[,colnames(cmes)=='lng200m'],1,mean),3),
                         round(apply(cmes[,colnames(cmes)=='lng200m'],1,quantile, probs=0.025),3),
                         round(apply(cmes[,colnames(cmes)=='lng200m'],1,quantile, probs=0.975),3))
  
  avg.cme2$sig <- ifelse(avg.cme2[,2]<0 & avg.cme2[,3]>0, '', '*')
  
  avg.cme2$full <- paste0(avg.cme2[,1], ' (',
                          avg.cme2[,2], ', ',
                          avg.cme2[,3], ')')
  
  colnames(avg.cme2) <- c(paste('cme',period,'avg',sep='_'),
                          paste('cme',period,'lci',sep='_'),
                          paste('cme',period,'uci',sep='_'),
                          'sig',
                          paste('cme',period, sep='_'))
 
  message <- paste('green 200m: There is significant mediation effect in', paste('pre',period,sep='_'), sep=' ')
  message2 <-paste('green 200m: There is no mediation effect in', paste('pre',period,sep='_'), sep=' ')
  ifelse((avg.cme2['tt','sig'] == '*') & (avg.cme2['cme','sig'] == '*'), message, message2)

  # 300m buffer
  avg.cme3 <- data.frame(round(apply(cmes[,colnames(cmes)=='lng300m'],1,mean),3),
                         round(apply(cmes[,colnames(cmes)=='lng300m'],1,quantile, probs=0.025),3),
                         round(apply(cmes[,colnames(cmes)=='lng300m'],1,quantile, probs=0.975),3))
  
  avg.cme3$sig <- ifelse(avg.cme3[,2]<0 & avg.cme3[,3]>0, '', '*')
  
  avg.cme3$full <- paste0(avg.cme3[,1], ' (',
                          avg.cme3[,2], ', ',
                          avg.cme3[,3], ')')
  
  colnames(avg.cme3) <- c(paste('cme',period,'avg',sep='_'),
                          paste('cme',period,'lci',sep='_'),
                          paste('cme',period,'uci',sep='_'),
                          'sig',
                          paste('cme',period,sep='_'))
  
  message <- paste('green 300m: There is significant mediation effect in', paste('pre',period,sep='_'), sep=' ')
  message2 <-paste('green 300m: There is no mediation effect in', paste('pre',period,sep='_'), sep=' ')
  ifelse((avg.cme3['tt','sig'] == '*') & (avg.cme3['cme','sig'] == '*'), message, message2)
  
  # 500m buffer
  avg.cme4 <- data.frame(round(apply(cmes[,colnames(cmes)=='300m'],1,mean),3),
                         round(apply(cmes[,colnames(cmes)=='300m'],1,quantile, probs=0.025),3),
                         round(apply(cmes[,colnames(cmes)=='300m'],1,quantile, probs=0.975),3))
  
  avg.cme4$sig <- ifelse(avg.cme4[,2]<0 & avg.cme4[,3]>0, '', '*')
  
  avg.cme4$full <- paste0(avg.cme4[,1], ' (',
                          avg.cme4[,2], ', ',
                          avg.cme4[,3], ')')
  
  colnames(avg.cme4) <- c(paste('cme',period,'avg',sep='_'),
                          paste('cme',period,'lci',sep='_'),
                          paste('cme',period,'uci',sep='_'),
                          'sig',
                          paste('cme',period,sep='_'))
  
  message <- paste('green 300m: There is significant mediation effect in', paste('pre',period,sep='_'), sep=' ')
  message2 <-paste('green 300m: There is no mediation effect in', paste('pre',period,sep='_'), sep=' ')
  ifelse((avg.cme4['tt','sig'] == '*') & (avg.cme4['cme','sig'] == '*'), message, message2)
  
  # final results
  all.risk.s[[m]] <- avg.risks.s[,5]
  all.risk.o[[m]] <- avg.risks.o[,5]
  all.pip[[m]] <- avg.pips
  all.pred[[m]] <- avg.preds[,1]
  all.cme[[m]] <- avg.cme[,5]
  all.cme2[[m]] <- avg.cme2[,5]

  }

all.RISK.s <- do.call(cbind, all.risk.s)
rownames(all.RISK.s) <- c('0.25_pm10','0.25_no2','0.25_pm25','0.25_temp','0.25_humi',
                          '0.50_pm10','0.50_no2','0.50_pm25','0.50_temp','0.50_humi',
                          '0.75_pm10','0.75_no2','0.75_pm25','0.75_temp','0.75_humi')

all.RISK.o <- do.call(cbind, all.risk.o)
rownames(all.RISK.o) <- seq(from=0.10, to=0.90, by=0.05)

all.PIP <- do.call(rbind, all.pip)
all.PIP <- t(data.frame(pips_avg=round(apply(all.PIP,2,mean),3),
                        pips_std=round(apply(all.PIP,2,sd),3)))

all.PRED <- do.call(cbind, all.pred)
rownames(all.PRED) <- c('mix_avg','mix_std')

all.CME <- do.call(cbind, all.cme)
rownames(all.CME) <- c('cme','cme.se','de','de.se','tt','tt.se','pp','pp.se')


all.CME2 <- do.call(cbind, all.cme2)
rownames(all.CME2) <- c('cme','cme.se','de','de.se','tt','tt.se','pp','pp.se')

all.CME3 <- do.call(cbind, all.cme3)
rownames(all.CME3) <- c('cme','cme.se','de','de.se','tt','tt.se','pp','pp.se')

all.CME4 <- do.call(cbind, all.cme4)
rownames(all.CME4) <- c('cme','cme.se','de','de.se','tt','tt.se','pp','pp.se')


#######################R E S U L T#########################################
# MDI, prenatal period
write.csv(all.RISK.s, "MDI_PRE_bkmr_single.csv", row.names=TRUE)
write.csv(all.RISK.o,"MDI_PRE_bkmr_overall.csv", row.names=TRUE)
write.csv(all.PIP,"MDI_PRE_bkmr_pip.csv", row.names=TRUE)
write.csv(all.PRED,"MDI_PRE_bkmr_mixture.csv", row.names=TRUE)
write.csv(all.CME,"MDI_PRE_lng100m_causal_mediation.csv", row.names=TRUE)
write.csv(all.CME2,"MDI_PRE_lng200m_causal_mediation.csv", row.names=TRUE)
write.csv(all.CME3,"MDI_PRE_lng300m_causal_mediation.csv", row.names=TRUE)
write.csv(all.CME4,"MDI_PRE_lng500m_causal_mediation.csv", row.names=TRUE)

# PDI, prenatal period
write.csv(all.RISK.s, "PDI_PRE_bkmr_single.csv", row.names=TRUE)
write.csv(all.RISK.o,"PDI_PRE_bkmr_overall.csv", row.names=TRUE)
write.csv(all.PIP,"PDI_PRE_bkmr_pip.csv", row.names=TRUE)
write.csv(all.PRED,"PDI_PRE_bkmr_mixture.csv", row.names=TRUE)
write.csv(all.CME,"PDI_PRE_lng100m_causal_mediation.csv", row.names=TRUE)
write.csv(all.CME2,"PDI_PRE_lng200m_causal_mediation.csv", row.names=TRUE)
write.csv(all.CME3,"PDI_PRE_lng300m_causal_mediation.csv", row.names=TRUE)
write.csv(all.CME4,"PDI_PRE_lng500m_causal_mediation.csv", row.names=TRUE)
