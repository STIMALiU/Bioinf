## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this code or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .

## modified solution of Owlright from
##https://stackoverflow.com/questions/14196551/add-isoclines-and-or-direction-field-to-plot

library(deSolve) 
library(phaseR)

## ======================================================================
## edit these two variables to control the aesthetics of the flow field
ff_points<-15 ## controls the density of points for the flow field
ff_frac<-1.1 ## controls length of arrows for the flow field
## ======================================================================

model.OU <- function(t, y, parameters){
  with(as.list(parameters),{

    OU1<-y[1] 
    OU2<-y[2]
    dOU1 <- (-1)*( a11*(OU1-theta1) + a12*(OU2-theta2))
    dOU2 <- (-1)*( a21*(OU1-theta1) + a22*(OU2-theta2))

    list(c(dOU1,dOU2))
  })
}



P1<-cbind(c(1,1),c(-1,5))
P1com<-cbind(c(1+1i,5-1i),c(-1-1i,5+1i))
mL1L2<-matrix(NA,ncol=2,nrow=10)
mL1L2<-rbind(
c(-0.01,-0.05), #1
c(0.01,0.05), #2
c(-0.01,0.05), #3
c(NA,NA), #4
c(NA,NA), #5
c(-0.01,-0.01), #6
c(0.01,0.01), #7
c(-0.01+0.05i,-0.01-0.05i), #8
c(0.01+0.05i,0.01-0.05i), #9
c(0.02i,-0.02i)  #10
)
invP1<-solve(P1)
invP1com<-solve(P1com)

lA<-vector("list",nrow(mL1L2))
lA<-sapply(1:nrow(mL1L2),function(i,mL1L2,P1,invP1,P1com,invP1com){
    if(i>7){
	res<-P1com%*%diag(mL1L2[i,])%*%invP1com
    }else{
	res<-P1%*%diag(mL1L2[i,])%*%invP1
    }
    Re(res)
},mL1L2=mL1L2,P1=P1,invP1=invP1,P1com=P1com,invP1com=invP1com,simplify=FALSE)
lA[[4]]<-rbind(c(0.02,0.01),c(0,0.02))
lA[[5]]<-rbind(c(-0.02,0.01),c(0,-0.02))
mXY_lim<-rbind(
c(-5,5,-5,5),#1
c(-2,2,-2,2),#2
c(-5,5,-5,5),#3
c(-2,2,-2,2),#4
c(-5,5,-5,5),#5
c(-5,5,-5,5),#6
c(-2,2,-2,2),#7
c(-10,10,-10,10),#8
c(-2.5,2.5,-2.5,2.5),#9
c(-5,5,-5,5)#10
)

m_y0<-rbind(
c(-1,1,1,-1), #1
c(-1,1,1,-1), #2
c(-1,1,1,-1), #3
c(-1,1,1,-1), #4
c(-1,1,1,-1), #5
c(-1,1,1,-1), #6
c(-1,1,1,-1), #7
c(-0.1,0.1,0.1,-0.1), #8
c(-1,1,1,-1), #9
c(-1,-1,0.5,0.5) #10
)

sapply(1:length(lA),function(i,lA,mXY_lim,m_y0){
	mA<-lA[[i]]
	params.OU<-c(a11=mA[1,1],a12=mA[1,2],a21=mA[2,1],a22=mA[2,2],theta1=0,theta2=0)
	data.OU<-as.data.frame(lsoda(c(OU1=0,OU2=0),seq(1,250,by=0.01), model.OU, params.OU))
	#v_xlim<-c(-5,5);v_ylim<-c(-5,5) #c(-0,1,0.1)
	v_xlim<-mXY_lim[i,1:2];v_ylim<-mXY_lim[i,3:4];
	png(paste0("PPODE_setup_",i,".png"))

    # plot the trajectories of the system
    plot(data.OU$OU1, data.OU$OU2, type="l", col="blue", main="",
        xlab="", ylab="", xlim=v_xlim, ylim=v_ylim)
    #add Nullclines
    nullclines(model.OU, xlim=v_xlim,ylim=v_ylim, parameters=params.OU, system="two.dim", col=c("blue","blue"), add=TRUE,add.legend=FALSE,lwd=2)
    flowField(model.OU, xlim=v_xlim,ylim=v_ylim, parameters=params.OU, system="two.dim", add=TRUE,add.legend=FALSE,frac=ff_frac,points=ff_points)#, col=c("black"))
    trajectory(model.OU,tlim=c(0,1000),y0=rbind(m_y0[i,1:2],m_y0[i,3:4]), parameters=params.OU, system="two.dim",frac=ff_frac)
    dev.off()
},lA=lA,mXY_lim=mXY_lim,m_y0=m_y0,simplify=FALSE
)




