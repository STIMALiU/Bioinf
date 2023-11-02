## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this code or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .

library(igraph)


f_plot5stateDNAHMM<-function(png_filename,v_blue_edge,v_blue_emit){
    HMMadjm<-matrix(0,ncol=12,nrow=12)
    HMMadjm[1,c(2,7)]<-1

    ## P protein-coding region 
     HMMadjm[2,c(3,8)]<-1
     HMMadjm[3,c(4,9)]<-1
     HMMadjm[4,c(5,10)]<-1
     HMMadjm[5,c(6,11)]<-1
     HMMadjm[6,c(12)]<-1

    ## N non-coding region 
     HMMadjm[7,c(3,8)]<-1
     HMMadjm[8,c(4,9)]<-1
     HMMadjm[9,c(5,10)]<-1
     HMMadjm[10,c(6,11)]<-1
     HMMadjm[11,c(12)]<-1

    gHMM <- igraph::graph_from_adjacency_matrix( HMMadjm, mode="directed" ) 
    coords<-rbind(
    c(2,7.5), ##S
    c(5,8),c(8,8),c(11,8),c(14,8),c(17,8), ## P
    c(5,7),c(8,7),c(11,7),c(14,7),c(17,7), ## N
    c(20,7.5) ##E
    )

    v_transition_probs<-c(
	"1/2","1/2",
	"2/3","1/3",
	"2/3","1/3",
	"2/3","1/3",
	"2/3","1/3",
	"1",
	"1/3","2/3",
	"1/3","2/3",
	"1/3","2/3",
	"1/3","2/3",
	"1"
    )

    m_edge_lab<-rbind(
	c(-1+1/18,-0.5), #"1/2"
	c(-1+1/18,0.5), #"1/2"
	c(-1/2,1-1/12), #"2/3"
	c(-0.69,4/12), #"1/3"
	c(-5/32,1-1/12), #"2/3"
	c(-0.32,3/12), #"1/3"
	c(5/32,1-1/12), #"2/3"
	c(0,3/12), #"1/3"
	c(1/2,1-1/12), #"2/3"
	c(0,-3/12), #"1/3"
	c(1-1/12,-0.5), #"1"
	c(-0.69,-4/12), #"1/3"
	c(-1/2,-1+1/12), #"2/3" 
	c(-0.32,-3/12), #"1/3"
	c(-5/32,-1+1/12), #"2/3"
	c(0.33,3/12), #"1/3"
	c(5/32,-1+1/12), #"2/3"
	c(0.33,-3/12), #"1/3"
	c(1/2,-1+1/12), #"2/3"
	c(1-1/12,0.5) #"1" 
    )

    v_edge_lab_x<-m_edge_lab[,1]
    v_edge_lab_y<-m_edge_lab[,2]

    v_edge_color<-rep("gray",20)
    v_edge_color[v_blue_edge]<-"blue"
    #https://stackoverflow.com/questions/60311851/position-edge-labels-on-igraph-plot-in-r

    png(png_filename)
	plot(gHMM, vertex.size = 15, edge.width = 3, vertex.color = "coral",main="",vertex.label=c("S",rep("P",5),rep("N",5),"E"),layout=coords,
	vertex.label.cex=2,vertex.label.font=2,vertex.label.col="black",
	edge.label=v_transition_probs,edge.label.cex=2,edge.label.font=2,edge.label.col="black",
	edge.label.x = v_edge_lab_x,edge.label.y =v_edge_lab_y,edge.color=v_edge_color)

	v_emit_cols<-rep("black",40)
	v_emit_cols[v_blue_emit]<-"blue"
	j<-1

	mtext('0.25', side=1, line=0, at=-0.7,font=2,cex=1.25,col=v_emit_cols[j]);j<-j+1
	mtext('0.25', side=1, line=1, at=-0.7,font=2,cex=1.25,col=v_emit_cols[j]);j<-j+1
	mtext('0.25', side=1, line=2, at=-0.7,font=2,cex=1.25,col=v_emit_cols[j]);j<-j+1
	mtext('0.25', side=1, line=3, at=-0.7,font=2,cex=1.25,col=v_emit_cols[j]);j<-j+1

	mtext('0.25', side=1, line=0, at=-0.33,font=2,cex=1.25,col=v_emit_cols[j]);j<-j+1
	mtext('0.25', side=1, line=1, at=-0.33,font=2,cex=1.25,col=v_emit_cols[j]);j<-j+1
	mtext('0.25', side=1, line=2, at=-0.33,font=2,cex=1.25,col=v_emit_cols[j]);j<-j+1
	mtext('0.25', side=1, line=3, at=-0.33,font=2,cex=1.25,col=v_emit_cols[j]);j<-j+1

	mtext('0.25', side=1, line=0, at=0,font=2,cex=1.25,col=v_emit_cols[j]);j<-j+1
	mtext('0.25', side=1, line=1, at=0,font=2,cex=1.25,col=v_emit_cols[j]);j<-j+1
	mtext('0.25', side=1, line=2, at=0,font=2,cex=1.25,col=v_emit_cols[j]);j<-j+1
	mtext('0.25', side=1, line=3, at=0,font=2,cex=1.25,col=v_emit_cols[j]);j<-j+1

	mtext('0.25', side=1, line=0, at=0.33,font=2,cex=1.25,col=v_emit_cols[j]);j<-j+1
	mtext('0.25', side=1, line=1, at=0.33,font=2,cex=1.25,col=v_emit_cols[j]);j<-j+1
	mtext('0.25', side=1, line=2, at=0.33,font=2,cex=1.25,col=v_emit_cols[j]);j<-j+1
	mtext('0.25', side=1, line=3, at=0.33,font=2,cex=1.25,col=v_emit_cols[j]);j<-j+1

	mtext('0.25', side=1, line=0, at=0.7,font=2,cex=1.25,col=v_emit_cols[j]);j<-j+1
	mtext('0.25', side=1, line=1, at=0.7,font=2,cex=1.25,col=v_emit_cols[j]);j<-j+1
	mtext('0.25', side=1, line=2, at=0.7,font=2,cex=1.25,col=v_emit_cols[j]);j<-j+1
	mtext('0.25', side=1, line=3, at=0.7,font=2,cex=1.25,col=v_emit_cols[j]);j<-j+1

	mtext('0.2', side=3, line=-0.5, at=-0.7,font=2,cex=1.25,col=v_emit_cols[j]);j<-j+1
	mtext('0.3', side=3, line=0.5, at=-0.7,font=2,cex=1.25,col=v_emit_cols[j]);j<-j+1
	mtext('0.3', side=3, line=1.5, at=-0.7,font=2,cex=1.25,col=v_emit_cols[j]);j<-j+1
	mtext('0.2', side=3, line=2.5, at=-0.7,font=2,cex=1.25,col=v_emit_cols[j]);j<-j+1

	mtext('0.2', side=3, line=-0.5, at=-0.33,font=2,cex=1.25,col=v_emit_cols[j]);j<-j+1
	mtext('0.3', side=3, line=0.5, at=-0.33,font=2,cex=1.25,col=v_emit_cols[j]);j<-j+1
	mtext('0.3', side=3, line=1.5, at=-0.33,font=2,cex=1.25,col=v_emit_cols[j]);j<-j+1
	mtext('0.2', side=3, line=2.5, at=-0.33,font=2,cex=1.25,col=v_emit_cols[j]);j<-j+1

	mtext('0.2', side=3, line=-0.5, at=0,font=2,cex=1.25,col=v_emit_cols[j]);j<-j+1
	mtext('0.3', side=3, line=0.5, at=0,font=2,cex=1.25,col=v_emit_cols[j]);j<-j+1
	mtext('0.3', side=3, line=1.5, at=0,font=2,cex=1.25,col=v_emit_cols[j]);j<-j+1
	mtext('0.2', side=3, line=2.5, at=0,font=2,cex=1.25,col=v_emit_cols[j]);j<-j+1

	mtext('0.2', side=3, line=-0.5, at=0.33,font=2,cex=1.25,col=v_emit_cols[j]);j<-j+1
	mtext('0.3', side=3, line=0.5, at=0.33,font=2,cex=1.25,col=v_emit_cols[j]);j<-j+1
	mtext('0.3', side=3, line=1.5, at=0.33,font=2,cex=1.25,col=v_emit_cols[j]);j<-j+1
	mtext('0.2', side=3, line=2.5, at=0.33,font=2,cex=1.25,col=v_emit_cols[j]);j<-j+1
	
	mtext('0.2', side=3, line=-0.5, at=0.7,font=2,cex=1.25,col=v_emit_cols[j]);j<-j+1
	mtext('0.3', side=3, line=0.5, at=0.7,font=2,cex=1.25,col=v_emit_cols[j]);j<-j+1
	mtext('0.3', side=3, line=1.5, at=0.7,font=2,cex=1.25,col=v_emit_cols[j]);j<-j+1
	mtext('0.2', side=3, line=2.5, at=0.7,font=2,cex=1.25,col=v_emit_cols[j]);j<-j+1
	
    dev.off()
}

#f_plot5stateDNAHMM("HMMDNAgraph.png",c(1,3,5,8,19,20),c(22,25,32,13,20))
## A,C,C,G,T
f_plot5stateDNAHMM("HMMDNAgraph.png",c(1,3,5,8,19,20),c(24,27,31,15,20))
f_plot5stateDNAHMM("HMMDNAgraph_2.png",c(2,13,15,17,18,11),c(1,6,10,15,37))
f_plot5stateDNAHMM("HMMDNAgraph_3.png",c(1,4,14,8,18,11),c(24,6,31,15,37))

