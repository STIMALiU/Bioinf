## This file accompanies the course 732A51 Bioinformatics 

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this code or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .


library(igraph)
library(ape)
## based on igraph's examples and 
## https://web.stanford.edu/class/bios221/NetworkSlides.html

adjm_1<-matrix(sample(0:1,100,replace=TRUE, prob=c(0.9,0.1)),nc=10)
g1<-igraph::graph_from_adjacency_matrix(adjm_1)

adjm_2<-apply((adjm_1+t(adjm_1)),c(1,2),function(x){min(x,1)})
g2<-igraph::graph_from_adjacency_matrix(adjm_1,mode="max") 
g3<-igraph::graph_from_adjacency_matrix(adjm_2,mode="undirected") 



phyltree<-rtree(10)
graph_phyltree<-igraph::graph_from_edgelist(phyltree$edge,directed=TRUE)

par(mfrow=c(2,3))
plot(g1,vertex.size=25,edge.width=5,vertex.color="coral")
plot(g2,vertex.size=25,edge.width=5,vertex.color="coral")
plot(g3,vertex.size=25,edge.width=5,vertex.color="coral")
plot(phyltree)
plot(graph_phyltree,vertex.size=10,edge.width=1,vertex.color="coral")


g_full<-igraph::make_full_graph(10)
g_bipartite<-igraph::make_full_bipartite_graph(5,5)
g_ring<-igraph::make_ring(10)
m_edges_path<-cbind(1:9,2:10)
g_path<-igraph::graph_from_edgelist(m_edges_path,directed=FALSE)

png("Typical_graphs.png")
par(mfrow=c(2,2))
plot(g_full,vertex.size=25,edge.width=4,vertex.color="coral",main="Full graph")
plot(g_bipartite,vertex.size=25,edge.width=4, vertex.color="coral",main="Full bipartite graph")
plot(g_ring,vertex.size=25,edge.width=4,vertex.color="coral",main="Ring graph (cycle)")
plot(g_path,vertex.size=10,edge.width=4,vertex.color="coral",main="Path")
dev.off()

