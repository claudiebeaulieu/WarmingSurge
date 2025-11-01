#####################################
###### Temperature Anomalies Analysis
#####################################

library(WeightedPortTest) # for the portmanteau test
library(EnvCpt) # for changepoint detection
library(changepoint) # for changepoint detection

# load data
load('temperature_anomalies.RData')
# use Tanom_annual_df, matrix, 6 columns, first year
# at the time of analysis, the Japan Met dataset had not been updated and is not included in the paper

# setup for PORTMANTEAU results - Table S1
names = c("NASA","Japan Met","HadCRUT","NOAA","Berkeley")
FGstatdisc = matrix(data=NA,nrow=4,ncol=5,dimnames=list(c("IID","GlobalAR(1)","GlobalAR(4)","ChangingAR(1)"),names))
FGstatcont = matrix(data=NA,nrow=4,ncol=5,dimnames=list(c("IID","GlobalAR(1)","GlobalAR(4)","ChangingAR(1)"),names))
LAG=20

# Each of the below sets of code fits a single model to all the datasets, extracts the changepoints,
# saves the results, then calculates the residuals from the fit and performs a portmanteau test


############### Discontinuous model IID - Figure 2 in the paper
trend=list()
mresiduals=list()
fits=list()
for(i in 2:6){
  data=Tanom_annual_df[,c(1,i)][!is.na(Tanom_annual_df[,i]),]
  n=nrow(data)
  trend[[i]]=EnvCpt:::cpt.reg(cbind(data[,2],rep(1,n),data[,1]),penalty="BIC",method="PELT",minseglen=10)
  plot(trend[[i]],main=names(Tanom_annual_df)[i])
  abline(v=cpts(trend[[i]]),col='blue')
  cpts(trend[[i]])=data[cpts(trend[[i]]),1]  # put cpts in terms of year
}
cptstrend=lapply(trend[-1],FUN=function(x){cpts(x)})

# calculate iid fit residuals and perform portmanteau test; row 1 is for iid discontinuous fit
for (i in 1:5){
  times=Tanom_annual_df[!is.na(Tanom_annual_df[,(i+1)]),1]
  if(i==2) times=c(times,2023)
  n=length(times)
  y=Tanom_annual_df[(times-1849),(i+1)]
  changetimes=c(times[1]-1,cptstrend[i][[1]],2023)
  ns=diff(changetimes)
  nchanges=length(changetimes)-1
  X=matrix(nrow=n,ncol=2*nchanges,0)
  lastchange=times[1]
  for(j in 1:nchanges){
    X[(changetimes[j]-times[1]+2):(changetimes[j+1]-times[1]+1),(2*j-1)]=rep(1,ns[j])
    X[(changetimes[j]-times[1]+2):(changetimes[j+1]-times[1]+1),(2*j)]=(changetimes[j]-times[1]+2):(changetimes[j+1]-times[1]+1)
  }
  lmfit=lm(y~0+X)
  fits[[i+1]] = lmfit$fitted.values
  FGstatdisc[1,i] = Weighted.Box.test(lmfit$resid,lag=LAG,type="Ljung")$p.value
  mresiduals[[i+1]] = lmfit$resid
}

# trend contains results for all datasets
save(trend,cptstrend,mresiduals,fits,file="./Results/resultstrend.Rdata")


############# Discontinuous model AR(1) 

source('./MethodCode/PELTtrendARp.R') # replacement trendARp function
trendar=list()
mresiduals=list()
fits=list()
dates=list()

for(i in 2:6){
  data=Tanom_annual_df[,c(1,i)][!is.na(Tanom_annual_df[,i]),]
  n=nrow(data)
  #Model fit
  trendar[[i]]=PELT.trendARp(data[,2],p=1,pen=5*log(n),minseglen=10)
  fittrend = fit.trendARp(data[,2],trendar[[i]],p=1,dates=data[,1],plot=T,add.ar=F,
                              title=names(Tanom_annual_df)[i])# get fit without AR - to visualize trend segments
  fits[[i]] = fittrend$fit
  dates[[i]] = fittrend$dates
  fittrendAR = fit.trendARp(data[,2],trendar[[i]],p=1,dates=data[,1],plot=T,add.ar=T,
                                title=names(Tanom_annual_df)[i])#get fit with AR - to compute residuals
  trendar[[i]]=data[trendar[[i]],1]  # put cpts in terms of year
  cat(paste(names(Tanom_annual_df)[i],':TrendAR1 \n'))
  print(fittrend$coeffs)
  #Compute residuals
  resid = data[,2] - fittrendAR$fit
  mresiduals[[i]] = resid
  FGstatdisc[4,i-1] = Weighted.Box.test(mresiduals[[i]],lag=LAG,type="Ljung",fitdf=1)$p.value
}
save(trendar,mresiduals,fits,dates,file="./Results/resultstrendar1.Rdata")



############ Discontinuous model with fixed AR(1)

source('./MethodCode/PELTtrendARp.R') # replacement trendAR1 function
trendFIXar=list()
mresiduals=list()
fits=list()
dates=list()
ar=c(NA,0.5367,0.4449,0.3669,0.4178,0.5402) # identified via EM from starting value of zero

for(i in 2:6){
  data=Tanom_annual_df[,c(1,i)][!is.na(Tanom_annual_df[,i]),]
  n=nrow(data)
  #Model fit
  trendFIXar[[i]]=PELT.trendFIXARp(data[,2],pen=5*log(n),arp=ar[i],minseglen=10)
  fittrend = fit.trendFIXARp(data[,2],trendFIXar[[i]],dates=data[,1],arp=ar[i],plot=T,add.ar=F,
                             title=names(Tanom_annual_df)[i])
  fits[[i]] = fittrend$fit
  dates[[i]] = fittrend$dates
  fittrendAR = fit.trendFIXARp(data[,2],trendFIXar[[i]],dates=data[,1],arp=ar[i],plot=F,fit=T,
                               add.ar=T,title=names(Tanom_annual_df)[i]) # save ar fit
  trendFIXar[[i]]=data[trendFIXar[[i]],1]  # put cpts in terms of year
  cat(paste(names(Tanom_annual_df)[i],':TrendFIXAR1 \n'))
  print(fittrendAR$coeffs)
  #Compute residuals
  resid = data[,2] - fittrendAR$fit
  mresiduals[[i]] = resid
  FGstatdisc[2,i-1]=Weighted.Box.test(mresiduals[[i]],lag=LAG,type="Ljung",fitdf=1)$p.value
}
save(trendFIXar,mresiduals,fits,dates,file="./Results/resultstrendFIXar1.Rdata")


############ Discontinuous model with fixed AR(4)

source('./MethodCode/PELTtrendARp.R') # replacement trendARp function
trendFIXarp=list()
mresiduals=list()
fits=list()
dates=list()
arp=matrix(c(rep(NA,4),
             0.5347,-0.1079,0.0991,0.1268,
             0.5046,-0.1608,0.0256,0.1594,
             0.4174,-0.1455,0.0442,0.0859,
             0.4697,-0.1384,0.04536,0.0632,
             0.5603,-0.1294,0.0446,0.1827),ncol=4,nrow=6,byrow=T)

for(i in 2:6){
  data=Tanom_annual_df[,c(1,i)][!is.na(Tanom_annual_df[,i]),]
  n=nrow(data)
  #Model fit
  trendFIXarp[[i]]=PELT.trendFIXARp(data[,2],pen=5*log(n),arp=arp[i,],minseglen=10)
  fittrend = fit.trendFIXARp(data[,2],trendFIXarp[[i]],dates=data[,1],arp=arp[i,],plot=T,add.ar=F,fit=F,
                             title=names(Tanom_annual_df)[i])# get fit without AR - to visualize trend segments
  fits[[i]] = fittrend$fit
  dates[[i]] = fittrend$dates
  fittrendAR = fit.trendFIXARp(data[,2],trendFIXarp[[i]],dates=data[,1],arp=arp[i,],plot=F,add.ar=T,fit=T,
                             title=names(Tanom_annual_df)[i])
  trendFIXarp[[i]]=data[trendFIXarp[[i]],1]  # put cpts in terms of year
  cat(paste(names(Tanom_annual_df)[i],':TrendAR4 \n'))
  print(fittrend$coeffs)
  #Compute residuals
  resid = data[,2] - fittrendAR$fit
  mresiduals[[i]] = resid
  FGstatdisc[3,i-1] = Weighted.Box.test(mresiduals[[i]],lag=LAG,type="Ljung",fitdf=4)$p.value
}
# trendFIXarp contains results for all datasets
save(trendFIXarp,mresiduals,fits,dates,file="./Results/resultstrendFIXar4.Rdata")



########### Continuous model IID

source("./MethodCode/optpelt2b.R") ##code for fitting continuous piecewise linear
source("./MethodCode/CROPS_optp.R") ##code for using CROPS
source("./MethodCode/FCPS.R") ##code for fitting function given changepoints

sigsquared=c(NA,0.0133,0.0114,0.00951,0.0100,0.0152)
jointrend=list()
varjt=list()
mresiduals=list()

for(i in 2:6){
  data=Tanom_annual_df[,c(1,i)][!is.na(Tanom_annual_df[,i]),]
  n=nrow(data)
  jointrend[[i]]=optp(data[,2],beta=3*log(n),sigsquared[i])
  functionfit=FCPS(data[,2],jointrend[[i]][[2]],sigsquared[i],beta=3*log(n))
  varjt[[i]]=var(data[,2]-functionfit[[4]])
  jointrend[[i]][[2]]=data[jointrend[[i]][[2]],1] # put cpts in terms of year

  #Compute residuals
  resid = data[,2] - functionfit[[4]]
  mresiduals[[i]] = resid
  FGstatcont[1,i-1] = Weighted.Box.test(mresiduals[[i]],lag=LAG,type="Ljung",fitdf=0)$p.value
}
cptsjointrend=lapply(jointrend[-1],FUN=function(x){x[[2]]})

save(jointrend,cptsjointrend,file="./Results/resultsjointrend.Rdata")



########### Continuous model AR(1)

source('./MethodCode/PELTtrendARpJOIN.R')
trendarjoin=list()
mresiduals=list()
fits=list()
dates=list()

for(i in 2:6){
  data=Tanom_annual_df[,c(1,i)][!is.na(Tanom_annual_df[,i]),]
  n=nrow(data)
  #Model fit
  trendarjoin[[i]]=PELT.trendARpJOIN(data[,2],p=1,pen=4*log(n),minseglen=10)
  fittrend = fit.trendARpJOIN(data[,2],trendarjoin[[i]],p=1,dates=data[,1],plot=T,add.ar=F,fit=T,
                              title=names(Tanom_annual_df)[i],pred=F)# get fit without AR - to visualize trend segments
  fits[[i]] = fittrend$fit
  dates[[i]] = fittrend$dates
  fittrendAR = fit.trendARpJOIN(data[,2],trendarjoin[[i]],p=1,dates=data[,1],plot=F,add.ar=T,fit=T,
                                title=names(Tanom_annual_df)[i]) #get fit with AR - to compute residuals below
  trendarjoin[[i]]=data[trendarjoin[[i]],1]  # put cpts in terms of year
  cat(paste(names(Tanom_annual_df)[i],':TrendAR1join \n'))
  print(fittrend$coeffs)
  #Compute residuals
  resid = data[,2] - fittrendAR$fit
  mresiduals[[i]] = resid
  FGstatcont[4,i-1]=Weighted.Box.test(resid,lag=LAG,type="Ljung",fitdf=1)$p.value
}
save(trendarjoin,mresiduals,fits,dates,file="./Results/resultstrendar1join.Rdata")



############ Continuous model fixed AR(1)

source('./MethodCode/PELTtrendARpJOIN.R')
trendFIXarjoin=list()
mresiduals=list()
fits=list()
dates=list()
ar=c(NA,0.5347,0.4683,0.2875,0.5052,0.5384)

for(i in 2:6){
  data=Tanom_annual_df[,c(1,i)][!is.na(Tanom_annual_df[,i]),]
  n=nrow(data)
  #Model fit
  trendFIXarjoin[[i]]=PELT.trendFIXARpJOIN(data[,2],pen=4*log(n),arp=ar[i],minseglen=10)
  fittrend = fit.trendFIXARpJOIN(data[,2],trendFIXarjoin[[i]],arp=ar[i],dates=data[,1],plot=T,add.ar=F,
                                 title=names(Tanom_annual_df)[i])
  fits[[i]] = fittrend$fit
  dates[[i]] = fittrend$dates
  fittrendAR = fit.trendFIXARpJOIN(data[,2],trendFIXarjoin[[i]],arp=ar[i],dates=data[,1],plot=F,fit=T,
                      add.ar=T,title=names(Tanom_annual_df)[i]) # save AR fit
  trendFIXarjoin[[i]]=data[trendFIXarjoin[[i]],1]  # put cpts in terms of year
  cat(paste(names(Tanom_annual_df)[i],':TrendFIXAR1join \n'))
  print(fittrendAR$coeffs)
  #Compute residuals
  resid = data[,2] - fittrendAR$fit
  mresiduals[[i]] = resid
  FGstatcont[2,i-1]=Weighted.Box.test(mresiduals[[i]],lag=LAG,type="Ljung",fitdf=1)$p.value
}
save(trendFIXarjoin,mresiduals,fits,dates,file="./Results/resultstrendFIXar1join.Rdata")



########### Continuous model fixed AR(4)

source('./MethodCode/PELTtrendARpJOIN.R')
trendFIXar4join=list()
mresiduals=list()
fits=list()
dates=list()

arp=matrix(c(rep(NA,4),
             0.5325,-0.0977,0.0907,0.1310,
             0.5066,-0.1437,0.0294,0.1963,
             0.3331,-0.1916,-0.0195,0.0407,
             0.5918,-0.0937,0.0909,0.1605,
             0.5616,-0.1148,0.0391,0.1877),ncol=4,nrow=6,byrow=T)

for(i in 2:6){
  data=Tanom_annual_df[,c(1,i)][!is.na(Tanom_annual_df[,i]),]
  n=nrow(data)
  #Model fit
  trendFIXar4join[[i]]=PELT.trendFIXARpJOIN(data[,2],pen=4*log(n),arp=arp[i,],minseglen=10)
  fittrend= fit.trendFIXARpJOIN(data[,2],trendFIXar4join[[i]],arp=arp[i,],dates=data[,1],plot=T,add.ar=F,
                                title=names(Tanom_annual_df)[i])
  fits[[i]] = fittrend$fit
  dates[[i]] = fittrend$dates
  fittrendAR = fit.trendFIXARpJOIN(data[,2],trendFIXar4join[[i]],arp=arp[i,],dates=data[,1],plot=F,fit=T,
                                   add.ar=T,title=names(Tanom_annual_df)[i]) # save AR fit
  trendFIXar4join[[i]]=data[trendFIXar4join[[i]],1]  # put cpts in terms of year
  cat(paste(names(Tanom_annual_df)[i],':TrendFIXAR4join \n'))
  print(fittrend$coeffs)
  #Compute residuals
  resid = data[,2] - fittrendAR$fit
  mresiduals[[i]] = resid
  FGstatcont[3,i-1]=Weighted.Box.test(mresiduals[[i]],lag=LAG,type="Ljung",fitdf=4)$p.value
}

save(trendFIXar4join,mresiduals,fits,dates,file="./Results/resultstrendFIXar4join.Rdata")

save(FGstatcont,FGstatdisc,file='./Results/FGstats.Rdata')


########### Continuous model AR(1) withholding 2023 

source('./MethodCode/PELTtrendARpJOIN.R')
trendarjoinm23=list()
predtrendarjoinm23=list()
fits=list()
dates=list()

datam23=Tanom_annual_df[-nrow(Tanom_annual_df),] # remove 2023 measurement
for(i in 2:6){
  data=datam23[,c(1,i)][!is.na(datam23[,i]),]
  n=nrow(data)
  
  trendarjoinm23[[i]]=PELT.trendARpJOIN(data[,2],p=1,pen=4*log(n),minseglen=10)
  fittrendARm23 = fit.trendARpJOIN(data[,2],trendarjoinm23[[i]],p=1,dates=data[,1],plot=F,add.ar=T,fit=T,
                                title=names(Tanom_annual_df)[i]) #get fit with AR
  fits[[i]] = fittrendARm23$fit
  dates[[i]] = fittrendARm23$dates
  predtrendarjoinm23[[i]]=fit.trendARpJOIN(data[,2],trendarjoinm23[[i]],p=1,dates=data[,1],plot=F,fit=T,add.ar=T,
                   title=names(Tanom_annual_df)[i],pred=T) # get 1 step ahead prediction
}
preds=NULL
sepred=NULL
for(i in 2:6){
  preds[i]=predtrendarjoinm23[[i]][[1]]$pred
  sepred[i]=predtrendarjoinm23[[i]][[1]]$se
}
predtrendarjoinm23=rbind(preds,sepred)
rm(preds,sepred)

Z=(Tanom_annual_df[nrow(Tanom_annual_df),]-predtrendarjoinm23[1,])/predtrendarjoinm23[2,]
pnorm(q=as.numeric(Z[2:6]))
# 0.9925837 0.9983212 0.9962647 0.9982221 0.9932941

save(predtrendarjoinm23,fits,dates,file="./Results/predtrendar1joinm2023.Rdata")

