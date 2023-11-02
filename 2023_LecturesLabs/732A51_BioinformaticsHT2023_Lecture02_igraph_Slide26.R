## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this code or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .


library(igraph)

MCadjm<-matrix(1,ncol=6,nrow=6)
MCadjm[,1]<-0
MCadjm[1,6]<-0
MCadjm[6,]<-0

gMC <- igraph::graph_from_adjacency_matrix( MCadjm, mode="directed" ) 
coords<-rbind(c(2,7.5),c(3,8),c(4,8),c(3,7),c(4,7),c(5,7.5))

#https://stackoverflow.com/questions/52410379/r-igraph-orientation-of-edge-beginning-and-ending-from-same-vertex
#Note that you need an angle for all edges, not just the loops, even though this only applies to loops. by G5W
png("MCDNAgraph.png")
plot(gMC, vertex.size = 15, edge.width = 3, vertex.color = "coral",main="",vertex.label=c("S","A","C","G","T","E"),layout=coords,
edge.loop.angle=c(
0,0,0,0,
pi,0,0,0,
0,0,0, ## this is the loop
0,0,0,0,0,
pi,0,0,0,
0,0,0, ## this is the loop
0
),vertex.label.cex=2,vertex.label.font=2,vertex.label.col="black"   
# https://stackoverflow.com/questions/13438013/variable-vertex-font-size-in-igraph
# one can control fonts for each vertex also
)
dev.off()

