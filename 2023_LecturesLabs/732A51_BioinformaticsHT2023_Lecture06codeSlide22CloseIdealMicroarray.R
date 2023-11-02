## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this code or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .


yx<-function(x,sdscale){abs(sapply(x,function(u,sdscale){rnorm(1,u,u/sdscale)},sdscale=sdscale,simplify=TRUE))};

plotlim<-12;sdscale<-4.5;plotlim2<-3

Cy3<-rexp(10000);
Cy5<-yx(Cy3,sdscale)
M<-log2(Cy5/Cy3);A<-(log2(Cy5)+log2(Cy3))/2

png("IdealMicroarray.png", width = 3*480, height = 480)
par(mfrow=c(1,3))
plot(Cy3,Cy5,pch=19,xlim=c(0,plotlim),ylim=c(0,plotlim),font.lab=2,xlab="Cy3",ylab="Cy5",cex.lab=1.5,cex.axis=1.5,font.axis=2);abline(a=0,b=0.5,lwd=2);abline(a=0,b=2,lwd=2);abline(a=0,b=1,lwd=2)
plot(log2(Cy3),log2(Cy5),pch=19,font.lab=2,xlab="log(Cy3)",ylab="log(Cy5)",cex.lab=1.5,cex.axis=1.5,font.axis=2);abline(a=-2,b=1,lwd=2);abline(a=2,b=1,lwd=2);abline(a=0,b=1,lwd=2)
plot(A,M,pch=19,font.lab=2,xlab="(logCy5+logCy3)/2",ylab="log(Cy5/Cy3)",xlim=c(-4*plotlim2,plotlim2),ylim=c(-plotlim2,plotlim2),cex.lab=1.5,cex.axis=1.5,font.axis=2);abline(a=2,b=0,lwd=2);abline(a=-2,b=0,lwd=2);abline(a=0,b=0,lwd=2)
dev.off()