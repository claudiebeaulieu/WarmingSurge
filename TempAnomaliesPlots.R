#Figures of different models fits considered in the main paper and supplementary

# Load dataset-updated for 2023 data except for Japan Met
load("temperature_anomalies.RData")

# Initialise figure params
cols = c("darkblue","orange","red","darkgrey")
names = c("NASA","HadCRUT","NOAA","Berkeley")
index = c(2,4,5,6)

##### 
# Continuous and discontinuous fits with changing AR (Figure 1 Main)
pdf("./Results/Fits_changingAR.pdf")
#Joinpin- changing AR(1)
load("./Results/resultstrendar1join.Rdata")
matplot(Tanom_annual_df$year,Tanom_annual_df[,index],type="l",col=cols,lwd=0.5,lty=1,xlab="Year",ylab="Anomaly (°C)",xlim=c(1850,2023),main="a)")
legend("topleft",legend=names,col=cols,lwd=2,cex=0.75)
for (i in 1:4){
  lines(dates[[index[i]]],fits[[index[i]]],col=cols[i],lwd=2)
}
#Discontinuous - changing AR(1)
load("./Results/resultstrendar1.Rdata")
matplot(Tanom_annual_df$year,Tanom_annual_df[,index],type="l",col=cols,lwd=0.5,lty=1,xlab="Year",ylab="Anomaly (°C)",xlim=c(1850,2023),main="b)")
legend("topleft",legend=names,col=cols,lwd=2,cex=0.75)
for (i in 1:4){
  lines(dates[[index[i]]],fits[[index[i]]],col=cols[i],lwd=2)
}
dev.off()

##### 
# Discontinuous fits IID (Figure 2 Main & Figure S4 Supp)
pdf("./Results/Fits_discontinous_IID.pdf")
load("./Results/resultstrend.Rdata")
for (i in 1:4){
  plot(Tanom_annual_df$year,Tanom_annual_df[,index[i]],ylab="Anomaly",xlab="",type="l",main=names[i])
  dates = Tanom_annual_df$year[!is.na(Tanom_annual_df[,index[i]])]
  lines(dates,fits[[index[i]]],col="red",lwd=2)
}
dev.off()


##### 
# Trend AR1Join fit without 2023 + 2023 predictions (Figure 3a)
pdf("./Results/2023_predicted.pdf")
# Fits without 2023 and predictions
load("./Results/predtrendar1joinm2023.Rdata")
matplot(Tanom_annual_df$year,Tanom_annual_df[,index],type="l",col=cols,lwd=0.5,lty=1,xlab="Year",ylab="Anomaly (°C)",xlim=c(1850,2023),main="")
legend("topleft",legend=c(names,"2023 observed","2023 predicted"),col=c(cols,"black","black"),lwd=c(2,2,2,2,2,2),lty=c(1,1,1,1,NA,NA),pch=c(NA,NA,NA,NA,21,15),cex=1)
for (i in 1:4){
  lines(dates[[index[i]]],fits[[index[i]]],col=cols[i],lwd=2)
  points(2023,predtrendarjoinm23[1,index[i]],pch=15,col=cols[i],lwd=2,cex=1)
  points(2023,Tanom_annual_df[174,index[i]],pch=21,col="black",bg=cols[i],lwd=2,cex=1)
}
# 2023 predictions + prediction intervals (Figure 3b)
plot(1,type="n",xlab="",ylab="",main="",xlim=c(1,4),ylim=c(0.5,1.5),xaxt="n",yaxt="n")
axis(side=4)
axis(1, at = c(1,2,3,4), labels=names)
for (i in 1:4){
  points(i,predtrendarjoinm23[1,index[i]],pch=15,col=cols[i],lwd=3,cex=2)
  points(i,Tanom_annual_df[174,index[i]],pch=21,col="black",bg=cols[i],lwd=3,cex=2)
  arrows(i,predtrendarjoinm23[1,index[i]]+2*predtrendarjoinm23[2,index[i]], i, predtrendarjoinm23[1,index[i]]-2*predtrendarjoinm23[2,index[i]], angle=90, code=3,length=0.1,col=cols[i])
}
dev.off()



#####
# Continuous and discontinuous fits with global AR(1) (Figure S1)
pdf("./Results/Fits_GlobalAR1.pdf")
#Joinpin- global AR(1)
load("./Results/resultstrendFIXar1join.Rdata")
matplot(Tanom_annual_df$year,Tanom_annual_df[,index],type="l",col=cols,lwd=0.5,lty=1,xlab="Year",ylab="Anomaly (°C)",xlim=c(1850,2023),main="a)")
legend("topleft",legend=names,col=cols,lwd=2,cex=0.75)
for (i in 1:4){
  lines(dates[[index[i]]],fits[[index[i]]],col=cols[i],lwd=2)
}
#Discontinuous - global AR(1)
load("./Results/resultstrendFIXar1.Rdata")
matplot(Tanom_annual_df$year,Tanom_annual_df[,index],type="l",col=cols,lwd=0.5,lty=1,xlab="Year",ylab="Anomaly (°C)",xlim=c(1850,2023),main="b)")
legend("topleft",legend=names,col=cols,lwd=2,cex=0.75)
for (i in 1:4){
  lines(dates[[index[i]]],fits[[index[i]]],col=cols[i],lwd=2)
}
dev.off()



#####
# Residuals of discontinuous model with fixed AR(1) (Figure S2)
pdf("./Results/pacf_res_GlobalAR1.pdf")
load("./Results/resultstrendFIXar1.Rdata")
for (i in 1:4){
  pacf(mresiduals[[index[i]]],lag.max=10,main=names[i])
}
dev.off()



#####
# Continuous and discontinuous fits with global AR(4) (Figure S3)
pdf("./Results/Fits_GlobalAR4.pdf")
#Joinpin- global AR(4)
load("./Results/resultstrendFIXar4join.Rdata")
matplot(Tanom_annual_df$year,Tanom_annual_df[,index],type="l",col=cols,lwd=0.5,lty=1,xlab="Year",ylab="Anomaly (°C)",xlim=c(1850,2023),main="a)")
legend("topleft",legend=names,col=cols,lwd=2,cex=0.75)
for (i in 1:4){
  lines(dates[[index[i]]],fits[[index[i]]],col=cols[i],lwd=2)
}
#Discontinuous - global AR(4)
load("./Results/resultstrendFIXar4.Rdata")
matplot(Tanom_annual_df$year,Tanom_annual_df[,index],type="l",col=cols,lwd=0.5,lty=1,xlab="Year",ylab="Anomaly (°C)",xlim=c(1850,2023),main="b)")
legend("topleft",legend=names,col=cols,lwd=2,cex=0.75)
for (i in 1:4){
  lines(dates[[index[i]]],fits[[index[i]]],col=cols[i],lwd=2)
}
dev.off()



#####
# Residuals of discontinuous model IID (Figure S5)
pdf("./Results/pacf_res_discontinous_IID.pdf")
load("./Results/resultstrend.Rdata")
for (i in 1:4){
  pacf(mresiduals[[index[i]]],lag.max=10,main=names[i])
}
dev.off()


##### Code for figure 4 in Main and Figures S6-S10 is in VantageSimulation.R

