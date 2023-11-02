## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this code or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .


f_sim_yuima<-function(yuima_obj,X0,step=0.01,time=1000){

    numsteps<-time/step
    samp_scheme<-yuima::setSampling(Terminal=time, n=numsteps)
    simulobj_yuima <- yuima::setYuima(model=yuima_obj, sampling=samp_scheme)
    simulobj_yuima_res <- yuima::simulate(simulobj_yuima, xinit=X0,space.discretized=TRUE)

    ## We want an (X1,X2) picture, yuima will draw two graphs (t,X1(t)), (t,X2(t))
    ## extract data from yuima object
    time_grid_length<-length(simulobj_yuima_res@sampling@grid[[1]])
    sdedim<-length(X0)
    
    sim_data<-NA
    sim_data<-simulobj_yuima_res@data@original.data@.Data
    time_data<-simulobj_yuima_res@sampling@grid[[1]]
    timepoints<-length(sim_data)/sdedim
    mSimData<-rbind(time_data[1:timepoints],matrix(sim_data,ncol=timepoints,byrow=TRUE))
    mSimData
}

