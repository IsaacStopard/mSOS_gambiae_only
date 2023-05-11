##Analyzing Survival##
##Load all required packages
library(survival) #survival analysis
library(eha) #event history analysis
library(reshape) #data restructuring
library(cluster) #cluster analysis
library(lattice) #graphics
library(foreign) #read in SPSS, Stata, text data
library(mstate) #multistate and competing risk analysis
library(TraMineR) #sequence analysis
library(Hmisc) #descriptive statistics

names(surv.temp)


##recode variables as factors, as necessary
temp.surv<-surv.temp
survtemp$expt<-as.factor(survtemp$expt)
survtemp$temp<-as.factor(survtemp$temp)


#rename factor levels
levels(survtemp$expt)[levels(survtemp$expt)=="1"]<-"Block 1"
levels(survtemp$expt)[levels(survtemp$expt)=="2"]<-"Block 2"

levels(survtemp$temp)[levels(survtemp$temp)=="21"]<-"21º C"
levels(survtemp$temp)[levels(survtemp$temp)=="24"]<-"24º C"
levels(survtemp$temp)[levels(survtemp$temp)=="27"]<-"27º C"
levels(survtemp$temp)[levels(survtemp$temp)=="30"]<-"30º C"
levels(survtemp$temp)[levels(survtemp$temp)=="32"]<-"32º C"
levels(survtemp$temp)[levels(survtemp$temp)=="34"]<-"34º C"


#KM Estimations

.Survfit<-survfit(Surv(day,status)~1,conf.type="log",conf.int=0.95,type="kaplan-meier",error="greenwood",data=survtemp)
summary(.Survfit)

#create new file 
survtemp.km.1<-survfit(Surv(day,status)~temp,type="kaplan-meier",data=survtemp,se.fit=T,conf.int=0.95,error="greenwood")
survtemp.km.2<-survfit(Surv(day,status)~temp+expt,type="kaplan-meier",data=survtemp,se.fit=T,conf.int=0.95,error="greenwood")
summary(survtemp.km.1)
summary(survtemp.km.2)

survblock1<-subset(survtemp,expt=="Block 1")
survblock2<-subset(survtemp,expt=="Block 2")
surv21<-subset(survtemp,temp=="21º C")
surv24<-subset(survtemp,temp=="24º C")
surv27<-subset(survtemp,temp=="27º C")
surv30<-subset(survtemp,temp=="30º C")
surv32<-subset(survtemp,temp=="32º C")
surv34<-subset(survtemp,temp=="34º C")

surv21.1<-subset(survblock1,temp=="21º C")
surv24.1<-subset(survblock1,temp=="24º C")
surv27.1<-subset(survblock1,temp=="27º C")
surv30.1<-subset(survblock1,temp=="30º C")
surv32.1<-subset(survblock1,temp=="32º C")
surv34.1<-subset(survblock1,temp=="34º C")

surv21.2<-subset(survblock2,temp=="21º C")
surv24.2<-subset(survblock2,temp=="24º C")
surv27.2<-subset(survblock2,temp=="27º C")
surv30.2<-subset(survblock2,temp=="30º C")
surv32.2<-subset(survblock2,temp=="32º C")
surv34.2<-subset(survblock2,temp=="34º C")


survtemp.km.fit.block1<-survfit(Surv(day,status)~temp,type="kaplan-meier",data=survblock1,se.fit=T,conf.int=0.95,error="greenwood")
survtemp.km.fit.block2<-survfit(Surv(day,status)~temp,type="kaplan-meier",data=survblock2,se.fit=T,conf.int=0.95,error="greenwood")

survtemp.km.fit.21<-survfit(Surv(day,status)~1,type="kaplan-meier",data=surv21,se.fit=T,conf.int=0.95,error="greenwood")
survtemp.km.fit.24<-survfit(Surv(day,status)~1,type="kaplan-meier",data=surv24,se.fit=T,conf.int=0.95,error="greenwood")
survtemp.km.fit.27<-survfit(Surv(day,status)~1,type="kaplan-meier",data=surv27,se.fit=T,conf.int=0.95,error="greenwood")
survtemp.km.fit.30<-survfit(Surv(day,status)~1,type="kaplan-meier",data=surv30,se.fit=T,conf.int=0.95,error="greenwood")
survtemp.km.fit.32<-survfit(Surv(day,status)~1,type="kaplan-meier",data=surv32,se.fit=T,conf.int=0.95,error="greenwood")
survtemp.km.fit.34<-survfit(Surv(day,status)~1,type="kaplan-meier",data=surv34,se.fit=T,conf.int=0.95,error="greenwood")


survtemp.km.fit.21.1<-survfit(Surv(day,status)~1,type="kaplan-meier",data=surv21.1,se.fit=T,conf.int=0.95,error="greenwood")
survtemp.km.fit.21.2<-survfit(Surv(day,status)~1,type="kaplan-meier",data=surv21.2,se.fit=T,conf.int=0.95,error="greenwood")
survtemp.km.fit.24.1<-survfit(Surv(day,status)~1,type="kaplan-meier",data=surv24.1,se.fit=T,conf.int=0.95,error="greenwood")
survtemp.km.fit.24.2<-survfit(Surv(day,status)~1,type="kaplan-meier",data=surv24.2,se.fit=T,conf.int=0.95,error="greenwood")
survtemp.km.fit.27.1<-survfit(Surv(day,status)~1,type="kaplan-meier",data=surv27.1,se.fit=T,conf.int=0.95,error="greenwood")
survtemp.km.fit.27.2<-survfit(Surv(day,status)~1,type="kaplan-meier",data=surv27.2,se.fit=T,conf.int=0.95,error="greenwood")
survtemp.km.fit.30.1<-survfit(Surv(day,status)~1,type="kaplan-meier",data=surv30.1,se.fit=T,conf.int=0.95,error="greenwood")
survtemp.km.fit.30.2<-survfit(Surv(day,status)~1,type="kaplan-meier",data=surv30.2,se.fit=T,conf.int=0.95,error="greenwood")
survtemp.km.fit.32.1<-survfit(Surv(day,status)~1,type="kaplan-meier",data=surv32.1,se.fit=T,conf.int=0.95,error="greenwood")
survtemp.km.fit.32.2<-survfit(Surv(day,status)~1,type="kaplan-meier",data=surv32.2,se.fit=T,conf.int=0.95,error="greenwood")
survtemp.km.fit.34.1<-survfit(Surv(day,status)~1,type="kaplan-meier",data=surv34.1,se.fit=T,conf.int=0.95,error="greenwood")
survtemp.km.fit.34.2<-survfit(Surv(day,status)~1,type="kaplan-meier",data=surv34.2,se.fit=T,conf.int=0.95,error="greenwood")



library(flexsurv)
gomp21b1<-flexsurvreg(Surv(day,status)~1,data=surv21.1,dist="gompertz")
weibull21b1<-flexsurvreg(Surv(day,status)~1,data=surv21.1,dist="weibull")
exp21b1<-flexsurvreg(Surv(day,status)~1,data=surv21.1,dist="exp")
loglog21b1<-flexsurvreg(Surv(day,status)~1,data=surv21.1,dist="llogis")
logn21b1<-flexsurvreg(Surv(day,status)~1,data=surv21.1,dist="lnorm")
AIC(gomp21b1,weibull21b1,exp21b1,loglog21b1,logn21b1)
df      AIC
gomp21b1     2 1880.692
weibull21b1  2 1894.909
exp21b1      1 1909.012
loglog21b1   2 1906.421
logn21b1     2 1916.377

gomp21b1
est      L95%     U95%     se     
shape  0.06449  0.04202  0.08696  0.01146
rate   0.01290  0.00980  0.01698  0.00181

N = 609,  Events: 201,  Censored: 408
Total time at risk: 8494
Log-likelihood = -938.3462, df = 2
AIC = 1880.692

##Getting median survival times##
median.gompertz<-function(shape,rate){qgompertz(.5,shape=shape,rate=rate)}
summary(gomp30b1,fn=median.gompertz,t=1,B=10000)


gomp24b1<-flexsurvreg(Surv(day,status)~1,data=surv24.1,dist="gompertz")
weibull24b1<-flexsurvreg(Surv(day,status)~1,data=surv24.1,dist="weibull")
exp24b1<-flexsurvreg(Surv(day,status)~1,data=surv24.1,dist="exp")
loglog24b1<-flexsurvreg(Surv(day,status)~1,data=surv24.1,dist="llogis")
logn24b1<-flexsurvreg(Surv(day,status)~1,data=surv24.1,dist="lnorm")
AIC(gomp24b1,weibull24b1,exp24b1,loglog24b1,logn24b1)
df      AIC
gomp24b1     2 2114.072
weibull24b1  2 2136.338
exp24b1      1 2146.343
loglog24b1   2 2154.105
logn24b1     2 2164.537

gomp24b1
est      L95%     U95%     se     
shape  0.06416  0.04312  0.08520  0.01073
rate   0.01500  0.01162  0.01935  0.00195

N = 629,  Events: 233,  Censored: 396
Total time at risk: 8541
Log-likelihood = -1055.036, df = 2
AIC = 2114.072

gomp27b1<-flexsurvreg(Surv(day,status)~1,data=surv27.1,dist="gompertz")
weibull27b1<-flexsurvreg(Surv(day,status)~1,data=surv27.1,dist="weibull")
exp27b1<-flexsurvreg(Surv(day,status)~1,data=surv27.1,dist="exp")
loglog27b1<-flexsurvreg(Surv(day,status)~1,data=surv27.1,dist="llogis")
logn27b1<-flexsurvreg(Surv(day,status)~1,data=surv27.1,dist="lnorm")
AIC(gomp27b1,weibull27b1,exp27b1,loglog27b1,logn27b1)
df      AIC
gomp27b1     2 1900.502
weibull27b1  2 1922.424
exp27b1      1 1948.976
loglog27b1   2 1939.657
logn27b1     2 1953.082

gomp27b1
est      L95%     U95%     se     
shape  0.09284  0.06785  0.11784  0.01275
rate   0.01436  0.01096  0.01882  0.00198

N = 611,  Events: 217,  Censored: 394
Total time at risk: 7087
Log-likelihood = -948.2509, df = 2
AIC = 1900.502

gomp30b1<-flexsurvreg(Surv(day,status)~1,data=surv30.1,dist="gompertz")
weibull30b1<-flexsurvreg(Surv(day,status)~1,data=surv30.1,dist="weibull")
exp30b1<-flexsurvreg(Surv(day,status)~1,data=surv30.1,dist="exp")
loglog30b1<-flexsurvreg(Surv(day,status)~1,data=surv30.1,dist="llogis")
logn30b1<-flexsurvreg(Surv(day,status)~1,data=surv30.1,dist="lnorm")
AIC(gomp30b1,weibull30b1,exp30b1,loglog30b1,logn30b1)
df      AIC
gomp30b1     2 1963.330
weibull30b1  2 1967.464
exp30b1      1 1983.718
loglog30b1   2 1975.313
logn30b1     2 1977.917

gomp30b1
est     L95%    U95%    se    
shape  0.0746  0.0442  0.1050  0.0155
rate   0.0236  0.0184  0.0303  0.0030

N = 648,  Events: 232,  Censored: 416
Total time at risk: 6110
Log-likelihood = -979.6651, df = 2
AIC = 1963.33



gomp32b1<-flexsurvreg(Surv(day,status)~1,data=surv32.1,dist="gompertz")
weibull32b1<-flexsurvreg(Surv(day,status)~1,data=surv32.1,dist="weibull")
exp32b1<-flexsurvreg(Surv(day,status)~1,data=surv32.1,dist="exp")
loglog32b1<-flexsurvreg(Surv(day,status)~1,data=surv32.1,dist="llogis")
logn32b1<-flexsurvreg(Surv(day,status)~1,data=surv32.1,dist="lnorm")
AIC(gomp32b1,weibull32b1,exp32b1,loglog32b1,logn32b1)
df      AIC
gomp32b1     2 1402.918
weibull32b1  2 1406.576
exp32b1      1 1425.592
loglog32b1   2 1413.972
logn32b1     2 1416.228

gomp32b1
shape  0.1006  0.0617  0.1396  0.0199
rate   0.0257  0.0192  0.0343  0.0038

N = 476,  Events: 174,  Censored: 302
Total time at risk: 3827
Log-likelihood = -699.4591, df = 2
AIC = 1402.918

gomp34b1<-flexsurvreg(Surv(day,status)~1,data=surv34.1,dist="gompertz")
weibull34b1<-flexsurvreg(Surv(day,status)~1,data=surv34.1,dist="weibull")
exp34b1<-flexsurvreg(Surv(day,status)~1,data=surv34.1,dist="exp")
loglog34b1<-flexsurvreg(Surv(day,status)~1,data=surv34.1,dist="llogis")
logn34b1<-flexsurvreg(Surv(day,status)~1,data=surv34.1,dist="lnorm")
AIC(gomp34b1,weibull34b1,exp34b1,loglog34b1,logn34b1)

df      AIC
gomp34b1     2 1968.127
weibull34b1  2 1975.173
exp34b1      1 2000.823
loglog34b1   2 1988.319
logn34b1     2 1988.353

gomp34b1
est      L95%     U95%     se     
shape  0.10440  0.07035  0.13844  0.01737
rate   0.03148  0.02483  0.03990  0.00381

N = 627,  Events: 256,  Censored: 371
Total time at risk: 4671
Log-likelihood = -982.0636, df = 2
AIC = 1968.127


gomp21b2<-flexsurvreg(Surv(day,status)~1,data=surv21.2,dist="gompertz")
weibull21b2<-flexsurvreg(Surv(day,status)~1,data=surv21.2,dist="weibull")
exp21b2<-flexsurvreg(Surv(day,status)~1,data=surv21.2,dist="exp")
loglog21b2<-flexsurvreg(Surv(day,status)~1,data=surv21.2,dist="llogis")
logn21b2<-flexsurvreg(Surv(day,status)~1,data=surv21.2,dist="lnorm")
AIC(gomp21b2,weibull21b2,exp21b2,loglog21b2,logn21b2)
df      AIC
gomp21b2     2 1531.252
weibull21b2  2 1552.045
exp21b2      1 1592.608
loglog21b2   2 1563.134
logn21b2     2 1582.569

gomp21b2
Estimates: 
  est      L95%     U95%     se     
shape  0.10014  0.07590  0.12438  0.01237
rate   0.00627  0.00447  0.00878  0.00108

N = 612,  Events: 158,  Censored: 454
Total time at risk: 8921
Log-likelihood = -763.6259, df = 2
AIC = 1531.252


gomp24b2<-flexsurvreg(Surv(day,status)~1,data=surv24.2,dist="gompertz")
weibull24b2<-flexsurvreg(Surv(day,status)~1,data=surv24.2,dist="weibull")
exp24b2<-flexsurvreg(Surv(day,status)~1,data=surv24.2,dist="exp")
loglog24b2<-flexsurvreg(Surv(day,status)~1,data=surv24.2,dist="llogis")
logn24b2<-flexsurvreg(Surv(day,status)~1,data=surv24.2,dist="lnorm")
AIC(gomp24b2,weibull24b2,exp24b2,loglog24b2,logn24b2)
df      AIC
gomp24b2     2 1825.788
weibull24b2  2 1842.794
exp24b2      1 1888.428
loglog24b2   2 1855.735
logn24b2     2 1877.503

gomp24b2
est      L95%     U95%     se     
shape  0.09291  0.07071  0.11511  0.01133
rate   0.00867  0.00644  0.01168  0.00132

N = 613,  Events: 196,  Censored: 417
Total time at risk: 8870
Log-likelihood = -910.8938, df = 2
AIC = 1825.788

gomp27b2<-flexsurvreg(Surv(day,status)~1,data=surv27.2,dist="gompertz")
weibull27b2<-flexsurvreg(Surv(day,status)~1,data=surv27.2,dist="weibull")
exp27b2<-flexsurvreg(Surv(day,status)~1,data=surv27.2,dist="exp")
loglog27b2<-flexsurvreg(Surv(day,status)~1,data=surv27.2,dist="llogis")
logn27b2<-flexsurvreg(Surv(day,status)~1,data=surv27.2,dist="lnorm")
AIC(gomp27b2,weibull27b2,exp27b2,loglog27b2,logn27b2)
df      AIC
gomp27b2     2 2022.034
weibull27b2  2 2059.906
exp27b2      1 2107.706
loglog27b2   2 2086.606
logn27b2     2 2113.219

gomp27b2
est      L95%     U95%     se     
shape  0.09950  0.07918  0.11982  0.01037
rate   0.01034  0.00789  0.01356  0.00143

N = 638,  Events: 229,  Censored: 409
Total time at risk: 8361
Log-likelihood = -1009.017, df = 2
AIC = 2022.034


gomp30b2<-flexsurvreg(Surv(day,status)~1,data=surv30.2,dist="gompertz")
weibull30b2<-flexsurvreg(Surv(day,status)~1,data=surv30.2,dist="weibull")
exp30b2<-flexsurvreg(Surv(day,status)~1,data=surv30.2,dist="exp")
loglog30b2<-flexsurvreg(Surv(day,status)~1,data=surv30.2,dist="llogis")
logn30b2<-flexsurvreg(Surv(day,status)~1,data=surv30.2,dist="lnorm")
AIC(gomp30b2,weibull30b2,exp30b2,loglog30b2,logn30b2)
df      AIC
gomp30b2     2 1881.492
weibull30b2  2 1899.545
exp30b2      1 1933.508
loglog30b2   2 1915.426
logn30b2     2 1929.366

gomp30b2
est      L95%     U95%     se     
shape  0.10397  0.07694  0.13100  0.01379
rate   0.01520  0.01162  0.01989  0.00209

N = 628,  Events: 219,  Censored: 409
Total time at risk: 6627
Log-likelihood = -938.746, df = 2
AIC = 1881.492

gomp32b2<-flexsurvreg(Surv(day,status)~1,data=surv32.2,dist="gompertz")
weibull32b2<-flexsurvreg(Surv(day,status)~1,data=surv32.2,dist="weibull")
exp32b2<-flexsurvreg(Surv(day,status)~1,data=surv32.2,dist="exp")
loglog32b2<-flexsurvreg(Surv(day,status)~1,data=surv32.2,dist="llogis")
logn32b2<-flexsurvreg(Surv(day,status)~1,data=surv32.2,dist="lnorm")
AIC(gomp32b2,weibull32b2,exp32b2,loglog32b2,logn32b2)
df      AIC
gomp32b2     2 2049.869
weibull32b2  2 2079.875
exp32b2      1 2119.524
loglog32b2   2 2104.020
logn32b2     2 2122.457

gomp32b2
est      L95%     U95%     se     
shape  0.11702  0.09047  0.14357  0.01355
rate   0.01611  0.01243  0.02088  0.00213

N = 640,  Events: 248,  Censored: 392
Total time at risk: 6520
Log-likelihood = -1022.935, df = 2
AIC = 2049.869


gomp34b2<-flexsurvreg(Surv(day,status)~1,data=surv34.2,dist="gompertz")
weibull34b2<-flexsurvreg(Surv(day,status)~1,data=surv34.2,dist="weibull")
exp34b2<-flexsurvreg(Surv(day,status)~1,data=surv34.2,dist="exp")
loglog34b2<-flexsurvreg(Surv(day,status)~1,data=surv34.2,dist="llogis")
logn34b2<-flexsurvreg(Surv(day,status)~1,data=surv34.2,dist="lnorm")
AIC(gomp34b2,weibull34b2,exp34b2,loglog34b2,logn34b2)
df      AIC
gomp34b2     2 2253.437
weibull34b2  2 2303.081
exp34b2      1 2388.637
loglog34b2   2 2349.205
logn34b2     2 2371.465

gomp34b2
Estimates: 
  est      L95%     U95%     se     
shape  0.15880  0.13249  0.18510  0.01342
rate   0.01724  0.01345  0.02210  0.00219

N = 630,  Events: 302,  Censored: 328
Total time at risk: 5778
Log-likelihood = -1124.719, df = 2
AIC = 2253.437

##Working with Full Dataset that is separated by individuals
head(spz.data)
names(spz.data)
spzdata<-as.data.frame(spz.data)

#change expt to factor 
spz.data$expt<-as.factor(spz.data$expt)
#change temp to factor
spz.data$temp<-as.factor(spz.data$temp)
#change cup to factor
spz.data$cup<-as.factor(spz.data$cup)
#change cup.id to factor
spz.data$cup.id<-as.factor(spz.data$cup.id)

#spz prevalence means by day
library(plyr)
prevplyr<-ddply(spz.data,c("expt","temp","day"),summarise, N=length(spz.prevalence),
                mean=mean(spz.prevalence),sd=sd(spz.prevalence),se=sd/sqrt(N))
prevplyr2<-ddply(spz.data,c("temp","day"),summarise, N=length(spz.prevalence),
                 mean=mean(spz.prevalence),sd=sd(spz.prevalence),se=sd/sqrt(N))

#subset by temperature and block

T21<-subset(prevplyr,temp=="21")
T21.1<-subset(T21, expt=="0")
T21.2<-subset(T21, expt=="1")
T24<-subset(prevplyr,temp=="24")
T24.1<-subset(T24, expt=="0")
T24.2<-subset(T24, expt=="1")
T27<-subset(prevplyr,temp=="27")
T27.1<-subset(T27,expt=="0")
T27.2<-subset(T27,expt=="1")
T30<-subset(prevplyr,temp=="30")
T30.1<-subset(T30,expt=="0")
T30.2<-subset(T30,expt=="1")
T32<-subset(prevplyr,temp=="32")
T32.1<-subset(T32,expt=="0")
T32.2<-subset(T32,expt=="1")
T34<-subset(prevplyr,temp=="34")
T34.1<-subset(T34,expt=="0")
T34.2<-subset(T34,expt=="1")

##FOR THOSE THAT CAN BE FIT SIGMOIDALLY###
g=asymptote
k=rate
t=inflection point


#for all 21C
fitmodel21<-nls(mean~g/(1+exp(-k*(day-t))),data=T21,start=list(g=.4,k=2,t=12))
model: mean ~ g/(1 + exp(-k * (day - t)))
data: T21
g       k       t 
0.4799  0.6217 14.6124 
residual sum-of-squares: 0.132

Number of iterations to convergence: 10 
Achieved convergence tolerance: 4.391e-06

#for 21C block 1
fitmodel21.1<-nls(mean~g/(1+exp(-k*(day-t))),data=T21.1,start=list(g=.4,k=2,t=12))
model: mean ~ g/(1 + exp(-k * (day - t)))
data: T21.1
g       k       t 
0.3966  0.9207 13.9720 
residual sum-of-squares: 0.02976

Number of iterations to convergence: 6 
Achieved convergence tolerance: 5.451e-06

#for 21C block 2
fitmodel21.2<-nls(mean~g/(1+exp(-k*(day-t))),data=T21.2,start=list(g=.4,k=2,t=12))
model: mean ~ g/(1 + exp(-k * (day - t)))
data: T21.2
g       k       t 
0.6374  0.4073 16.1279 
residual sum-of-squares: 0.04465

Number of iterations to convergence: 10 
Achieved convergence tolerance: 3.606e-06

##plotting the data + curves
op<-par(mfrow=c(1,1),mar=c(5,5,3,2))
plab<-"Proportion Infectious"; tlab="Day Post Infection"
plot(0,ylim=c(0,1),xlim=c(5,25),xlab=tlab,ylab=plab,main = "21 C")
points(mean~day,data=T21.2,col="blue",pch=16,cex=1.7)
points(mean~day,data=T21.3,col="red",pch=16,cex=1.7)
curve(coef(fitmodel21)[1]/(1+exp(-coef(fitmodel21)[2]*(x-coef(fitmodel21)[3]))),from=5,to=25,add=TRUE,lwd=2)
curve(coef(fitmodel21)[1]/(1+exp(-coef(fitmodel21)[2]*(x-coef(fitmodel21)[3]))),from=5,to=25,add=TRUE,lwd=4)
with (data = T21.1, expr = errbar(day, mean, mean+se, 
                                  mean-se, add=T, pch=16, col="blue",cap=.05,lwd=2,cex=1.7))
with(data=T21.2,expr=errbar(day, mean,mean+se,mean-se,add=T,pch=16,col="firebrick1",cap=.05,lwd=2,cex=1.7))


#for all 24C
fitmodel24<-nls(mean~g/(1+exp(-k*(day-t))),data=T24,start=list(g=.4,k=2,t=12))

Parameters:
  Estimate Std. Error t value Pr(>|t|)    
g  0.60077    0.02679  22.426 3.76e-16 ***
  k  2.02607    0.76478   2.649    0.015 *  
  t 11.50812    0.22585  50.954  < 2e-16 ***
  ---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.1051 on 21 degrees of freedom

Number of iterations to convergence: 7 
Achieved convergence tolerance: 9.316e-06

#for 24C block 1
fitmodel24.1<-nls(mean~g/(1+exp(-k*(day-t))),data=T24.1,start=list(g=.4,k=2,t=12))
model: mean ~ g/(1 + exp(-k * (day - t)))
data: T21.1

Parameters:
  Estimate Std. Error t value Pr(>|t|)    
g  0.51127    0.02781  18.381 1.91e-08 ***
  k  3.51482    2.22668   1.579    0.149    
t 11.34207    0.26783  42.349 1.14e-11 ***
  ---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.08286 on 9 degrees of freedom

Number of iterations to convergence: 9 
Achieved convergence tolerance: 1.47e-06

#for 24C block 2
fitmodel24.2<-nls(mean~g/(1+exp(-k*(day-t))),data=T24.2,start=list(g=.4,k=2,t=12))

Estimate Std. Error t value Pr(>|t|)    
g  0.70364    0.02174   32.37 1.26e-10 ***
  k  1.55048    0.34528    4.49  0.00151 ** 
  t 11.66065    0.16608   70.21 1.22e-13 ***
  ---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.05541 on 9 degrees of freedom

Number of iterations to convergence: 8 
Achieved convergence tolerance: 2.354e-06

#for all 27C
fitmodel27<-nls(mean~g/(1+exp(-k*(day-t))),data=T27,start=list(g=.4,k=2,t=12))
Parameters:
  Estimate Std. Error t value Pr(>|t|)    
g  0.59520    0.01798  33.099   <2e-16 ***
  k  4.96957    1.83895   2.702    0.013 *  
  t 10.50947    0.18618  56.447   <2e-16 ***
  ---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.07616 on 22 degrees of freedom

Number of iterations to convergence: 7 
Achieved convergence tolerance: 1.021e-06

#for 27C block 1
fitmodel27.1<-nlsLM(mean~g/(1+exp(-k*(day-t))),data=T27.1,start=list(g=.6,k=2,t=10))
model: mean ~ g/(1 + exp(-k * (day - t)))

Parameters:
  Estimate Std. Error t value Pr(>|t|)    
g 5.733e-01  1.817e-02  31.550 2.41e-11 ***
  k 1.827e+01  1.266e+06   0.000    1.000    
t 1.089e+01  7.283e+03   0.001    0.999    
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.05746 on 10 degrees of freedom

Number of iterations to convergence: 16 
Achieved convergence tolerance: 1.49e-08

#for 27C block 2
fitmodel27.2<-nls(mean~g/(1+exp(-k*(day-t))),data=T27.2,start=list(g=.4,k=2,t=12))

Parameters:
  Estimate Std. Error t value Pr(>|t|)    
g  0.62255    0.03245  19.188 1.31e-08 ***
  k  4.85184    3.77834   1.284    0.231    
t 10.37154    0.31944  32.468 1.23e-10 ***
  ---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.09158 on 9 degrees of freedom

Number of iterations to convergence: 10 
Achieved convergence tolerance: 1.571e-06

#Truncated Data for 30C, 32C, 34C
library(plyr)
truncplyr<-ddply(trunc.dat.final,c("expt","temp","day"),summarise, N=length(spz.prev),
                 mean=mean(spz.prev),sd=sd(spz.prev),se=sd/sqrt(N))

#subset by temperature and block
T27trunc<-subset(truncplyr,temp=="27")
T27.1.t<-subset(T27trunc,expt=="0")

T30trunc<-subset(truncplyr,temp=="30")
T30.1.t<-subset(T30trunc,expt=="0")
T30.2.t<-subset(T30trunc,expt=="1")

T32trunc<-subset(truncplyr,temp=="32")
T32.1.t<-subset(T32trunc,expt=="0")
T32.2.t<-subset(T32trunc,expt=="1")

T34trunc<-subset(truncplyr,temp=="34")
T34.1.t<-subset(T34trunc,expt=="0")
T34.2.t<-subset(T34trunc,expt=="1")

#30C Both Blocks

fitmodel30.t<-nls(mean~g/(1+exp(-k*(day-t))),data=T30trunc,start=list(g=.4,k=2,t=10))
Parameters:
  Estimate Std. Error t value Pr(>|t|)    
g  0.43589    0.01464  29.766 2.66e-10 ***
  k  3.11206    0.68774   4.525  0.00144 ** 
  t  8.28692    0.09002  92.056 1.07e-14 ***
  ---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.03155 on 9 degrees of freedom

Number of iterations to convergence: 11 
Achieved convergence tolerance: 3.486e-06

#30C Block 1
fitmodel30.t1<-nls(mean~g/(1+exp(-k*(day-t))),data=T30.1.t,start=list(g=.4,k=2,t=10))
Parameters:
  Estimate Std. Error t value Pr(>|t|)    
g  0.42488    0.01149  36.990 3.19e-06 ***
  k  4.79988    1.87954   2.554   0.0631 .  
t  8.27130    0.11667  70.892 2.37e-07 ***
  ---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.0198 on 4 degrees of freedom

Number of iterations to convergence: 9 
Achieved convergence tolerance: 1.056e-06

#30C Block 2
fitmodel30.t2<-nls(mean~g/(1+exp(-k*(day-t))),data=T30.2.t,start=list(g=.4,k=2,t=10))
Parameters:
  Estimate Std. Error t value Pr(>|t|)    
g  0.45300    0.03253  13.924  0.00512 ** 
  k  2.31897    0.79583   2.914  0.10036    
t  8.28246    0.17779  46.584  0.00046 ***
  ---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.04069 on 2 degrees of freedom

Number of iterations to convergence: 16 
Achieved convergence tolerance: 7.138e-06

##Piecewise NLS models for 30,32,34

##30C##
##First you have to divide into "breaks" - I did this in Excel and entered each break in as a
####.csv file. There should be 2 breaks per temp x block combination

# 30C Block 1: BREAK 2: DAY 10+ #
Nonlinear regression model
model: spz ~ b * (exp(k2 * day))

k2        b 
-0.08134  1.04538 
residual sum-of-squares: 0.02433

Number of iterations to convergence: 3 
Achieved convergence tolerance: 1.49e-08

#30C Block 2: BREAK 2: DAY 12+ #
Nonlinear regression model
model: spz ~ x * (exp(k1 * (time)))

k1        x 
-0.09545  1.30647 
residual sum-of-squares: 0.00764

Number of iterations to convergence: 9 
Achieved convergence tolerance: 1.49e-08

#32C Both Blocks
fitmodel32.t<-nls(mean~g/(1+exp(-k*(day-t))),data=T32trunc,start=list(g=.4,k=2,t=10))
Parameters:
  Estimate Std. Error t value Pr(>|t|)    
g  0.40747    0.01846   22.07 4.40e-11 ***
  k  2.98967    0.85654    3.49  0.00446 ** 
  t  7.65590    0.13318   57.49 5.07e-16 ***
  ---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.04662 on 12 degrees of freedom

Number of iterations to convergence: 21 
Achieved convergence tolerance: 8.273e-06

#32C Block 1
fitmodel32.1.t<-nls(mean~g/(1+exp(-k*(day-t))),data=T32.1.t,start=list(g=.4,k=2,t=10))
Parameters:
  Estimate Std. Error t value Pr(>|t|)    
g  0.40186    0.02714  14.806 5.97e-06 ***
  k  2.36673    0.97595   2.425   0.0515 .  
t  7.65150    0.21870  34.986 3.63e-08 ***
  ---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.04982 on 6 degrees of freedom

Number of iterations to convergence: 44 
Achieved convergence tolerance: 9.247e-06

#32C Block 2
fitmodel32.2.t<-nls(mean~g/(1+exp(-k*(day-t))),data=T32.2.t,start=list(g=.4,k=2,t=10))
Parameters:
  Estimate Std. Error t value Pr(>|t|)    
g   0.4153     0.0295  14.078 0.000776 ***
  k   4.0420     2.2403   1.804 0.168957    
t   7.6743     0.2157  35.577 4.88e-05 ***
  ---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.05026 on 3 degrees of freedom

Number of iterations to convergence: 13 
Achieved convergence tolerance: 7.995e-06

#32C Block 1: BREAK 2: DAY 12+
Nonlinear regression model
model: spz ~ b * (exp(k2 * day))

k2        b 
-0.08134  1.04538 
residual sum-of-squares: 0.02433

Number of iterations to convergence: 3 
Achieved convergence tolerance: 1.49e-08

#32C Block 2 BREAK 2: DAY 12+
Nonlinear regression model
model: spz ~ x * (exp(k2 * (day)))

x      k2 
2.6239 -0.1782 
residual sum-of-squares: 0.0006344

Number of iterations to convergence: 4 
Achieved convergence tolerance: 5.199e-06

#34C Both Blocks
fitmodel34.t<-nls(mean~g/(1+exp(-k*(day-t))),data=T34trunc,start=list(g=.2,k=1,t=7))
Parameters:
  Estimate Std. Error t value Pr(>|t|)    
g  0.36295    0.01837  19.753 1.01e-08 ***
  k  4.01189    3.12565   1.284    0.231    
t  7.02849    0.07927  88.663 1.50e-14 ***
  ---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.03473 on 9 degrees of freedom

Number of iterations to convergence: 10 
Achieved convergence tolerance: 1.849e-06

#34C Block 1
fitmodel34.1.t<-nls(mean~g/(1+exp(-k*(day-t))),data=T34.1.t,start=list(g=.4,k=4,t=7))
Parameters:
  Estimate Std. Error t value Pr(>|t|)    
g  0.33259    0.01494  22.264 2.41e-05 ***
  k  3.09997    1.24739   2.485   0.0678 .  
t  6.94403    0.09421  73.711 2.03e-07 ***
  ---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.02296 on 4 degrees of freedom

Number of iterations to convergence: 12 
Achieved convergence tolerance: 7.188e-06

#34C Block 2
fitmodel34.2.t<-nls(mean~g/(1+exp(-k*(day-t))),data=T34.2.t,start=list(g=.4,k=4,t=7))
Parameters:
  Estimate Std. Error t value Pr(>|t|)    
g 0.411187   0.001248   329.4 9.22e-06 ***
  k 5.037714   0.381735    13.2  0.00569 ** 
  t 7.097291   0.008438   841.1 1.41e-06 ***
  ---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.001306 on 2 degrees of freedom

Number of iterations to convergence: 7 
Achieved convergence tolerance: 2.991e-06


#34 Block 1 BREAK 2: Day 10+
Nonlinear regression model
model: spz ~ x * (exp(k2 * (day)))

x      k2 
1.8530 -0.1663 
residual sum-of-squares: 0.006396

Number of iterations to convergence: 2 
Achieved convergence tolerance: 5.636e-06


#34 Block 2 BREAK 2: Day 8+
Nonlinear regression model
model: spz ~ x * (exp(k2 * (day)))
data: t34b3break2
x      k2 
3.7855 -0.2324 
residual sum-of-squares: 0.01257

Number of iterations to convergence: 8 
Achieved convergence tolerance: 1.491e-06
