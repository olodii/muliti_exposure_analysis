###################################################################################################
##########################Bayesian kernel machine regression#######################################
###################################################################################################
library(bkmr);library(ggplot2)
library(dplyr);library(writexl)
library(cvTools);library(foreach)
library(gam);library(mediation)
library(readxl);library(mgcv);library(dplyr)
library(caret)

setwd("C://Users//user//Desktop//최종")
select <- dplyr::select

load('data_formation.RData')

####<- Data[[1]]
exp <- dat %>% dplyr::select(pm10_1T.u:pm25_6m.u)
cov <- dat %>% dplyr::select(mtemp_1T:mhumi_6m, age,mbi,dbwt,wks,medu,sex,green:green4)
nd <- dat %>% dplyr::select(yc,ym)

dat2<- cbind(exp, cov, nd)

head(exp)
par(mfrow=c(1,1)); boxplot(exp)

#Correlation matrix
x11();par(mfrow=c(1,1)); corrplot::corrplot(cor(exp),method="square",type="upper",addCoef.col = TRUE,tl.col="black")


dat2$lng100m <- log1p(dat2$green)
dat2$lng200m <- log1p(dat2$green2)
dat2$lng300m <- log(dat2$green3)
dat2$lng500m <- log(dat2$green3)
dat2 %>% select(lng100m:lng500m) %>% summary()

dat2$medu <- as.factor(dat2$medu)
dat2$shs <- as.factor(dat2$shs)
dat2$sex <- as.factor(dat2$sex)
dat2$age<- as.factor(dat2$age)

dat2$dbwt - as.numeric(dat2$dbwt *0.001)
dat2$wks <- as.numeric(dat2$wks)
dat2$yc <- as.numeric(dat2$yc)
dat2$ym <- as.numeric(dat2$ym)

which(names(dat2)==c("mhumi_6m"))
dat2[,1:25] <- scale(dat2[,1:25], center=T, scale=T)

dummy <- dummyVars(" ~ .",data=dat2)
d_dat <- data.frame(predict(dummy, newdata=dat2))

D1 <- d_dat %>% dplyr::select(yc, ym, pm10_preg.u, no2_preg.u, pm25_preg.u,
                             mtemp_preg, mhumi_preg, age, wks, dbwt, mbi, medu.0:sex.1,
                             lng100m:lng500m)
D2 <- d_dat %>% dplyr::select(y_cogc, ym, pm10_1T.u, no2_1T.u, pm25_1T.u,
                             mtemp_1T, mhumi_1T, d_ma_age, age, wks, dbwt, mbi, medu.0:sex.1,
                             lng100m:lng500m)
D3 <- d_dat %>% dplyr::select(yc, ym, pm10_2T.u, no2_2T.u, pm25_2T.u,
                             mtemp_2T, mhumi_2T, d_ma_age, age, wks, dbwt, mbi, medu.0:sex.1,
                             lng100m:lng500m)
D4 <- d_dat %>% dplyr::select(yc, ym10_3T.u, no2_3T.u, pm25_3T.u,
                             mtemp_3T, mhumi_3T, d_ma_age, d_gaage, wks, dbwt, mbi, medu.0:sex.1,
                             lng100m:lng500m)
D5 <- d_dat %>% dplyr::select(y_cogcmotm10_6m.u, no2_6m.u, pm25_6m.u,
                             mtemp_6m, mhumi_6m, d_ma_age, d_gaage, wks, dbwt, mbi, medu.0:sex.1,
                             lng100m:lng500m)

#
dat<- Data[[2]]

exp <- dat %>% dplyr::select(pm10_1T.u:pm25_1yc.u)
cov <- dat %>% dplyr::select(mtemp_1T:mhumi_1y, age,mbi,dbwt,wks,medu,sex,green:green4)
nd <- dat %>% dplyr::select(yc1,ym1)

dat2<- cbind(exp, cov, nd)

dat2$lng100m <- log1p(dat2$green)
dat2$lng200m <- log1p(dat2$green2)
dat2$lng300m <- log(dat2$green3)
dat2$lng500m <- log(dat2$green4)

dat2$medu <- as.factor(dat2$medu)
dat2$shs <- as.factor(dat2$shs)
dat2$sex <- as.factor(dat2$sex)

dat2$age <- as.numeric(dat2$age)
dat2$dbwt <- as.numeric(dat2$dbwt *0.001)
dat2$wks <- as.numeric(dat2$wks)
dat2$yc1 <- as.numeric(dat2$yc1)
dat2$ym1 <- as.numeric(dat2$ym1)

which(names(dat2)==c("mhumi_1y"))
dat2[,1:33] <- scale(dat2[,1:33], center=T, scale=T)
dummy <- dummyVars(" ~ .",data=dat2)
d_dat <- data.frame(predict(dummy, newdata=dat2))

D6 <- d_dat %>% dplyr::select(y1_coc1_moc1, pm10_1yc.u, no2_1yc.u, pm25_1yc.u,
                             mtemp_1y, mhumi_1y, age, wks, dbwt, mbi, medu.0:sex.1,
                             lng100m:lng500m)
D7 <- d_dat %>% dplyr::select(yc1, yc110_preg.u, no2_preg.u, pm25_preg.u,
                             mtemp_preg, mhumi_preg, d_mage, wks, dbwt, mbi, medu.0:sex.1,
                             lng100m:lng500m)


# pr<- Data[[3]]

exp <- dat %>% dplyr::select(pm10_1T.u:pm25_2yc.u)
cov <- dat %>% dplyr::select(mtemp_1T:mhumi_2y, age,mbi,dbwt,wks,medu,sex,green:green4)
nd <- dat %>% dplyr::select(yc2,ym2)

dat2 <- cbind(exp, cov, nd)

dat2$lng100m <- log1p(dat2$green)
dat2$lng200m <- log1p(dat2$green2)
dat2$lng300m <- log(dat2$green300)
dat2$lng500m <- log(dat2$green500)


dat2$medu <- as.factor(dat2$medu)
dat2$shs <- as.factor(dat2$shs)
dat2$sex <- as.factor(dat2$sex)

dat2$age <- as.numeric(dat2$age)
dat2$dbwt <- as.numeric(dat2$dbwt *0.001)
dat2$wks <- as.numeric(dat2$wks)
dat2$yc2 <- as.numeric(dat2$yc2)
dat2$ym2 <- as.numeric(dat2$ym2)

which(names(dat2)==c("mhumi_2y"))
dat2[,1:35] <- scale(dat2[,1:35], center=T, scale=T)

dummy <- dummyVars(" ~ .",data=dat2)
d_dat <- data.frame(predict(dummy, newdata=dat2))

D8 <- d_dat %>% dplyr::select(y2_coc2_mom2, pm10_2yc.u, no2_2yc.u, pm25_2yc.u,
                             mtemp_2y, mhumi_2y, age, wks, dbwt, mbi, medu.0:sex.1,
                             lng100m:lng500m)
D9 <- d_dat %>% dplyr::select(yc2, ym210_preg.u, no2_preg.u, pm25_preg.u,
                             mtemp_preg, mhumi_preg, d_maage, wks, dbwt, mbi, medu.0:sex.1,
                             lng100m:lng500m)


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

x11();pairs(D1 %>% dplyr::select(pm10_preg.u, no2_preg.u, mtemp_preg, lngreen100m,yc, ym),
            lower.panel = panel.lm, upper.panel = panel.cor, diag.panel = panel.hist)