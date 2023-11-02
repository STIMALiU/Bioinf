## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this code or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .


library(yuima)
source("732A51_BioinformaticsHT2023_Lecture05code_YuimaSimF.R")

A <- c(0,0)
S <- t(chol(rbind(c(1,-0.5),c(-0.5,1)))) ## one provides Sigma, not Sigma%*%t(Sigma)
yuima_BM_2d <- yuima::setModel(drift = A, diffusion = S,
state.variable=c("X1","X2"),solve.variable=c("X1","X2") )
X0<-c(0,0)

vylim<-c(-50,50)
vxlim<-c(-50,50)

mSimData<-f_sim_yuima(yuima_BM_2d,X0)
png("mvBM_yuima_2d_1.png");plot(mSimData[2,],mSimData[3,],pch=19,cex=0.1,xlab="X1",ylab="X2",main="",cex.lab=1.75,cex.axis=1.5,ylim=vylim,xlim=vxlim);points(0,0,col="red",cex=2,pch=19);dev.off()

mSimData<-f_sim_yuima(yuima_BM_2d,X0)
png("mvBM_yuima_2d_2.png");plot(mSimData[2,],mSimData[3,],pch=19,cex=0.1,xlab="X1",ylab="X2",main="",cex.lab=1.75,cex.axis=1.5,ylim=vylim,xlim=vxlim);points(0,0,col="red",cex=2,pch=19);dev.off()

mSimData<-f_sim_yuima(yuima_BM_2d,X0)
png("mvBM_yuima_2d_3.png");plot(mSimData[2,],mSimData[3,],pch=19,cex=0.1,xlab="X1",ylab="X2",main="",cex.lab=1.75,cex.axis=1.5,ylim=vylim,xlim=vxlim);points(0,0,col="red",cex=2,pch=19);dev.off()

mSimData<-f_sim_yuima(yuima_BM_2d,X0)
png("mvBM_yuima_2d_4.png");plot(mSimData[2,],mSimData[3,],pch=19,cex=0.1,xlab="X1",ylab="X2",main="",cex.lab=1.75,cex.axis=1.5,ylim=vylim,xlim=vxlim);points(0,0,col="red",cex=2,pch=19);dev.off()





