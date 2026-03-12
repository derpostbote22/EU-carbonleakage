

# --------------------------------
# 0. Preparations
# --------------------------------

#clear work
rm(list =ls())

#import
library(lmtest)
library(tseries)
library(dplyr)
library(orcutt)
library(HoRM)
library(prais)
library(ggplot2)
library(ggfortify)
library(forecast)
library(stargazer)
library(tikzDevice)
library(cowplot)

# --------------------------------
# 1. Import data and preprocess
# --------------------------------

#open dataset
df <- read.csv("C:/Users/jt3022/OneDrive - Imperial College London/Output/Private and Former Research/TFM paper/Data analysis ground work/New filing system/2 Data processing & R/conference_base.csv")

#turn fp and months into factor
df$fp <- factor (df$fp)
df$fac_mon <- factor(df$month)

# turn important variables (NET?) into smaller portions
df$NET_adj <- df$NET/1000000
df$NET_adj2 <- df$NET/1000000000

# --------------------------------
# 2. Run regressions
# --------------------------------

#2a. OLS regression model
lrmod1 <- lm(NET_adj ~ euindout2 + lEUA:fp + fac_mon, data = df)
summary(lrmod1)

#set time series 
NET_ts <- ts(df$NET_adj, start = c(2005,5), end = c(2021,12), frequency = 12)

#draw data against fitted values
tsfit1 <- ts(fitted(lrmod1), start = c(2005,5), end = c(2021,12), frequency = 12)
NET_with_lrmod1 <- cbind(NET_ts, tsfit1)
autoplot(NET_with_lrmod1) + 
  ylab("Net Aluminium Imports (???MM)") + xlab("Time") +
  #ggtitle("EU Net Aluminium Imports (April 2005 - December 2021)") +  
  scale_color_manual(labels = c("Original data", "Fitted linear model"), values=c("black","blue")) +
  theme(legend.title = element_blank()) +
  theme(text=element_text(size=17)) +
  theme(legend.position = c(0.2, 0.85),
        legend.background = element_rect(fill = "white"))

#residuals
tsres1 = NET_ts - tsfit1
df$zeros <- replicate(200, 0)
zerross <- ts(df$zeros, start = c(2005,5), end = c(2021,12), frequency = 12)
resid1 <- cbind(tsres1, zerross)
autoplot(resid1) +
  #ggtitle("Residuals of linear regression") +
  xlab ("Time") +
  ylab("Residuals (???MM)") + 
  #ggtitle("EU Net Aluminium Imports (April 2005 - December 2021)") +  
  scale_color_manual(labels = c("Original data", "Fitted linear model"), values=c("black","blue")) +
  theme(legend.title = element_blank()) +
  theme(text=element_text(size=17)) +
  theme(legend.position = c(0.2, 0.85),
        legend.background = element_rect(fill = "white"))

#ACF and PACF of residuals
require(gridExtra)
plot1 <- ggAcf(tsres1, lag.max = 36) + ggtitle("Autocorrelation function for residuals of linear regression") +xlab ("lags") + theme(axis.title.y = element_blank()) + ylim(c(-0.3,0.7))
plot2 <- ggPacf(tsres1, lag.max = 36) + ggtitle("Partial autocorrelation function for residuals of linear regression") +xlab ("lags") + theme(axis.title.y = element_blank()) + ylim(c(-0.3,0.7))
grid.arrange(plot1, plot2, nrow = 2)


#bgtest
teststat = bgtest(lrmod1, order=1)
teststat$statistic


#set time series 
NET_ts <- ts(df$NET_adj, start = c(2005,5), end = c(2021,12), frequency = 12)

#2b. Cochrane Orcutt
coch1 = cochrane.orcutt(lrmod1)
summary(coch1)
coch1$rho

#2c. Prais Winsten
df$mtime <- 1:nrow(df)
pw1 <- prais_winsten(NET_adj ~ lEUA:fp + euindout2 + fac_mon, data = df, index = "mtime")
pw1$rho
summary(pw1)

#2d. Hildreth-Lu
#define hildreth lu modelling function
hildreth.lu.func<- function(r, model){
  x1 <- model.matrix(model)[,2]
  x2 <- model.matrix(model)[,14]
  x3 <- model.matrix(model)[,15]
  x4 <- model.matrix(model)[,16]
  x5 <- model.matrix(model)[,17]
  xm2 <- model.matrix(model)[,3]
  xm3 <- model.matrix(model)[,4]
  xm4 <- model.matrix(model)[,5]
  xm5 <- model.matrix(model)[,6]
  xm6 <- model.matrix(model)[,7]
  xm7 <- model.matrix(model)[,8]
  xm8 <- model.matrix(model)[,9]
  xm9 <- model.matrix(model)[,10]
  xm10 <- model.matrix(model)[,11]
  xm11 <- model.matrix(model)[,12]
  xm12 <- model.matrix(model)[,13]
  y <- model.response(model.frame(model))
  n <- length(y)
  m <- ncol(model.matrix(model))
  t <- 2:n
  y <- y[t]-r*y[t-1]
  x1 <- x1[t]-r*x1[t-1]
  x2 <- x2[t]-r*x2[t-1]
  x3 <- x3[t]-r*x3[t-1]
  x4 <- x4[t]-r*x4[t-1]
  x5 <- x5[t]-r*x5[t-1]
  xm2 <- xm2[t-1]
  xm3 <- xm3[t-1]
  xm4 <- xm4[t-1]
  xm5 <- xm5[t-1]
  xm6 <- xm6[t-1]
  xm7 <- xm7[t-1]
  xm8 <- xm8[t-1]
  xm9 <- xm9[t-1]
  xm10 <- xm10[t-1]
  xm11 <- xm11[t-1]
  xm12 <- xm12[t-1]
  return(lm(y~x1 + x2 + x3 + x4 + x5 + xm2 + xm3 + xm4 + xm5 + xm6 + xm7 + xm8 + xm9 + xm10 + xm11 + xm12))
}

#optimise rho
r <- c(seq(0,1.05, by= 0.01))
tab <- data.frame("rho" = r, "SSE" = sapply(r, function(i){deviance(hildreth.lu.func(i, lrmod1))}))
optrho <- which.min(round(tab, 4)[,2])

#optimise Breusch-Godfrey
tab2 <- data.frame("rho" = r, "SSE" = sapply(r, function(i){bgtest(hildreth.lu.func(i, lrmod1))$statistic}))
opt2rho <- which.min(round(tab2, 4)[,2])

#optimise Breusch-Godfrey(p-value)
tab3 <- data.frame("rho" = r, "SSE" = sapply(r, function(i){bgtest(hildreth.lu.func(i, lrmod1))$p.value}))
opt3rho <- which.min(round(tab, 4)[,2])

#optimise Box-Ljung
tab4 <- data.frame("rho" = r, "SSE" = sapply(r, function(i){Box.test(resid((hildreth.lu.func(i, lrmod1))), type="Ljung",lag=20,fitdf=1)$p.value}))
opt4rho <- which.min(round(tab, 4)[,2])

#2e. First Differences
#transform variables to differences
dfd <- data.frame(NET = diff(df$NET_adj),
                  lEUA = diff(df$lEUA),
                  euindout2 = diff(df$euindout2),
                  fp = df$fp[2:200])

# --------------------------------
# 3. Visualisations
# --------------------------------

#visualise rho finding
plot(tab$SSE ~ tab$rho , type = "l")
#abline(v = tab[tab$SSE==min(tab$SSE),"rho"], lty = 3)
abline(v=0.6, lty = 3)
round(tab, 4)[optrho,]

#visualise other test statistic
tab <- data.frame("rho" = r, "SSE" = sapply(r, function(i){deviance(hildreth.lu.func(i, lrmod1))}))
plot(tab3$SSE ~ tab3$rho , type = "l")
abline(v = tab3[tab3$SSE==min(tab3$SSE),"rho"], lty = 3)
round(tab3, 4)[opt3rho,]

#preprocess for final plot
tab$SSE <- tab$SSE/1000000

#final plot
par(mar = c(5, 4, 4, 6) + 0.3)
plot(tab$SSE ~ tab$rho , type = "l", xlab=expression(paste(delta)), ylab="", col="blue")
abline(v = tab[tab$SSE==min(tab$SSE),"rho"], lty = 3)
mtext(expression("Sum of Squared Errors (SSE) x10"^~6),side=2,line=2.5)
round(tab, 4)[optrho,]
par(new = TRUE)
plot(tab2$SSE ~ tab2$rho, type = "l", axes = FALSE, bty = "n", xlab = "", ylab = "", col="red", lty=2)
abline(v = tab[tab2$SSE==min(tab2$SSE),"rho"], lty = 3)
mtext("Lagrange Multiplier (LM)",side=4,line=2.5) #col="red"
axis(side=4, at = pretty(range(tab2$SSE)))
abline(v = 0, lty = 3)
abline(v = 1, lty = 3)
mtext('HL/CO/PW', side=3, line=0.25, at=0.85)
mtext(expression('LM'[min]), side=3, line=0.25, at=0.63)
mtext('OLS', side=3, line=0.25, at=0)
mtext('FD', side=3, line=0.25, at=1)
legend(0.05, 30, legend=c("SSE", "LM"),
       col=c("blue", "red"), lty=1:2, cex=0.8)
hlm <- hildreth.lu.func(0.6, lrmod1)
summary(hlm)



#fit linear model
dmod1 <- lm(NET ~ lEUA:fp + euindout2, data = dfd)
summary(dmod1)

# --------------------------------
# 4. Summarise regression results
# --------------------------------

#USe Hildreth Lu function for all five regressions

#OLS
OLS <- hildreth.lu.func(0, lrmod1)
summary(OLS)

#COM
coch1$rho
COM <- hildreth.lu.func(0.847, lrmod1)
summary(COM)

#PWE
summary(pw1)
PWE <- hildreth.lu.func(0.843, lrmod1)
summary(PWE)

#HLM
HLM <- hildreth.lu.func(0.86, lrmod1)
summary(HLM)

#CPH
CPH <- hildreth.lu.func(0.85, lrmod1)
summary(CPH)

#FDE
FDE <- hildreth.lu.func(1, lrmod1)
summary(FDE)

#LMmin
LMmin <- hildreth.lu.func(0.63, lrmod1)
summary(LMmin)

Box.test(resid(LMmin), type="Ljung",lag=20,fitdf=1)

# --------------------------------
# 5. Create summary table
# --------------------------------

#make results into Latex table using stargazer
library(stargazer)
stargazer(OLS, LMmin, CPH, FDE, title="Results", align=TRUE)

stargazer(lrmod1, coch1, pw1$model, hlm, dmod1, type = "html", out = "fit_lm.html")

# --------------------------------
# 6. Autoregression checks
# --------------------------------

bgtest(CPH)
bgtest(OLS)
bgtest(LMmin)
bgtest(PWE)
bgtest(COM)
bgtest(HLM)
bgtest(FDE)

dwtest(CPH)
dwtest(OLS)
dwtest(LMmin)
dwtest(PWE)
dwtest(COM)
dwtest(HLM)
dwtest(FDE)

Box.test(resid(CPH), type="Ljung",lag=20,fitdf=1)
Box.test(resid(OLS), type="Ljung",lag=20,fitdf=1)
Box.test(resid(LMmin), type="Ljung",lag=20,fitdf=1)
Box.test(resid(PWE), type="Ljung",lag=20,fitdf=1)
Box.test(resid(COM), type="Ljung",lag=20,fitdf=1)
Box.test(resid(HLM), type="Ljung",lag=20,fitdf=1)
Box.test(resid(FDE), type="Ljung",lag=20,fitdf=1)



dwtest(LMmin)
dwtest(lrmod1)

Box.test(resid(LMmin), type="Ljung",lag=20,fitdf=1)$p.value

dwtest(PWE)
dwtest(HLM)


#compare performance of procedures
bgtest(coch1)
bgtest(COM)

summary(coch1)
summary(COM)

#check performance of LMmin
bgtest(LMmin)
dwtest(LMmin)
Box.test(resid(LMmin), type="Ljung",lag=1,fitdf=1)




