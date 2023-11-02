## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this code or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .


library(yuima)
source("732A51_BioinformaticsHT2023_Lecture05code_YuimaSimF.R")

## we want to spiral in towards the optimum, cf. with ODE system for mean
A <-c("-(X1)-0.5*X2","-(X2)+0.5*X1")
S <- t(chol(rbind(c(0.4,0.15),c(0.15,0.4)))) ## one provides Sigma, not Sigma%*%t(Sigma)
yuima_OU_2d <- yuima::setModel(drift = A, diffusion = S,
state.variable=c("X1","X2"),solve.variable=c("X1","X2") )

step<-0.0001
time<-7.5 ## convergence to optimum is rapid

vylim<-c(-2.5,2.5)
vxlim<-c(-2.5,2.5)

X0<-c(0,0)
mSimData<-f_sim_yuima(yuima_OU_2d,X0,step=step,time=time)
png("mvOU_yuima_2d_X000.png");plot(mSimData[2,],mSimData[3,],pch=19,cex=0.1,xlab="X1",ylab="X2",main="",cex.lab=1.75,cex.axis=1.5,ylim=vylim,xlim=vxlim);points(0,0,col="red",cex=2,pch=19);dev.off()

X0<-c(5,5)
mSimData<-f_sim_yuima(yuima_OU_2d,X0,step=step,time=time)
png("mvOU_yuima_2d_X055.png");plot(mSimData[2,],mSimData[3,],pch=19,cex=0.1,xlab="X1",ylab="X2",main="",cex.lab=1.75,cex.axis=1.5,ylim=vylim,xlim=vxlim);points(0,0,col="red",cex=2,pch=19);dev.off()

X0<-c(-5,-5)
mSimData<-f_sim_yuima(yuima_OU_2d,X0,step=step,time=time)
png("mvOU_yuima_2d_X0m5m5.png");plot(mSimData[2,],mSimData[3,],pch=19,cex=0.1,xlab="X1",ylab="X2",main="",cex.lab=1.75,cex.axis=1.5,ylim=vylim,xlim=vxlim);points(0,0,col="red",cex=2,pch=19);dev.off()

X0<-c(5,-5)
mSimData<-f_sim_yuima(yuima_OU_2d,X0,step=step,time=time)
png("mvOU_yuima_2d_X05m5.png");plot(mSimData[2,],mSimData[3,],pch=19,cex=0.1,xlab="X1",ylab="X2",main="",cex.lab=1.75,cex.axis=1.5,ylim=vylim,xlim=vxlim);points(0,0,col="red",cex=2,pch=19);dev.off()





