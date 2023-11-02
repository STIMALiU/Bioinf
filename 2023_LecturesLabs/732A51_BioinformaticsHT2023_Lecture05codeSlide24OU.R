## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this code or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .


library(yuima)
A <- c("-X")
S <- 1
yuima_OU_1d <- yuima::setModel(drift = A, diffusion = S,
state.variable=c("X"),solve.variable=c("X") )

step<-0.001
time<-100
numsteps<-time/step
samp_scheme<-yuima::setSampling(Terminal=time, n=numsteps)
simulobj_yuima <- yuima::setYuima(model=yuima_OU_1d, sampling=samp_scheme)

vylim=c(-5.5,10.5)


X0<-c(0)
simulobj_yuima_res_X00 <- yuima::simulate(simulobj_yuima, xinit=X0,space.discretized=TRUE)
png("OU_yuima_1d_X00.png");plot(simulobj_yuima_res_X00,lwd=3,cex=2,cex.axis=1.5,ylim=vylim,ylab="X(t)");dev.off()

X0<-c(5)
simulobj_yuima_res_X05 <- yuima::simulate(simulobj_yuima, xinit=X0,space.discretized=TRUE)
png("OU_yuima_1d_X05.png");plot(simulobj_yuima_res_X05,lwd=3,cex=2,cex.axis=1.5,ylim=vylim,ylab="X(t)");dev.off()

X0<-c(-5)
simulobj_yuima_res_X0m5 <- yuima::simulate(simulobj_yuima, xinit=X0,space.discretized=TRUE)
png("OU_yuima_1d_X0m5.png");plot(simulobj_yuima_res_X0m5,lwd=3,cex=2,cex.axis=1.5,ylim=vylim,ylab="X(t)");dev.off()

X0<-c(10)
simulobj_yuima_res_X010 <- yuima::simulate(simulobj_yuima, xinit=X0,space.discretized=TRUE)
png("OU_yuima_1d_X010.png");plot(simulobj_yuima_res_X010,lwd=3,cex=2,cex.axis=1.5,ylim=vylim,ylab="X(t)");dev.off()



