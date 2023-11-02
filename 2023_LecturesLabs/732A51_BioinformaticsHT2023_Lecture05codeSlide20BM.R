## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this code or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .


library(yuima)
A <- 0
S <- 1
yuima_BM_1d <- yuima::setModel(drift = A, diffusion = S,
state.variable=c("X(t)"),solve.variable=c("X(t)") )
X0<-c(0)
step<-0.01
time<-1000
numsteps<-time/step
samp_scheme<-yuima::setSampling(Terminal=time, n=numsteps)
simulobj_yuima <- yuima::setYuima(model=yuima_BM_1d, sampling=samp_scheme)
simulobj_yuima_res <- yuima::simulate(simulobj_yuima, xinit=X0,space.discretized=TRUE)
plot(simulobj_yuima_res,lwd=3,cex=2,cex.axis=1.5,ylim=c(-60,60))


