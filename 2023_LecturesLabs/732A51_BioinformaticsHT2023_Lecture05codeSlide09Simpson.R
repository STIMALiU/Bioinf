## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this code or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .


## This code replicates the idea behind Figure E1 in
## K. Bartoszek, J. Pienaar, P. Mostad, S. Andersson, and T. F. Hansen, 2012, 
## A phylogenetic comparative method for studying multivariate adaptation. 
## J. Theor. Biol 314:204-215.
## Requires mvsLOUCH >= 2.7.6 to run correctly

library(mvSLOUCH)
library(ape)
library(TreeSim)



f_rstartree<-function(n,branch_len=1){
    phyltree<-stree(n)
    phyltree$edge.length<-rep(branch_len,nrow(phyltree$edge))
    phyltree
}

f_YuleHeightFixed<-function(n,height=1){
    phyltree<-TreeSim::sim.bd.taxa(n,1,1,0)[[1]]
    phyltree$root.edge<-NULL
    phyltree$edge.length<-height*phyltree$edge.length/max(ape::node.depth.edgelength(phyltree))
    phyltree
}

num_traits<-2

v_clustsize<-c(30,30,30)
joining_branchlengths<-c(10,NA)
phylo_phyltree<-mvSLOUCH::simulate_clustered_phylogeny(v_clustsize,joining_branchlengths=joining_branchlengths,f_simclustphyl=f_YuleHeightFixed,joiningphyl=f_rstartree,b_change_joining_branches=TRUE) #,branch_len=joining_branchlengths[1])

v_regimes<-rep("reg0",sum(sapply(phylo_phyltree$edges_clusters,function(x){length(x)},simplify=TRUE)))
v_regimes[phylo_phyltree$edges_clusters[[2]]]<-"reg1"
v_regimes[phylo_phyltree$edges_clusters[[3]]]<-"reg2"
v_regimes[phylo_phyltree$edges_clusters[[4]]]<-"reg3"


mPsi<-cbind(c(0,0),c(5,5),c(10,10),c(15,15))
colnames(mPsi)<-c("reg0","reg1","reg2","reg3")
OUOUparameters<-list(vY0=matrix(0,nrow=num_traits,ncol=1),
A=cbind(c(4,0.4),c(3,10)),mPsi=mPsi,Syy=rbind(c(0.5,-1.5),c(0,1)))
OUOUdata<-mvSLOUCH::simulOUCHProcPhylTree(phylo_phyltree,OUOUparameters,regimes=v_regimes)



png("ThreeClustphyltree.png")
ape::plot.phylo(phylo_phyltree,show.tip.label=FALSE,type="phylogram",direction="rightwards",edge.width=3,x.lim=c(-0.05,joining_branchlengths[1]+2))
points(11.5,15,pch=15,cex=2)
points(11.5,45,pch=17,cex=2)
points(11.5,75,pch=19,cex=2)
dev.off()
png("SimpsonsParadox.png")
vpch<-rep(15,length(phylo_phyltree$tip.label))
vpch[phylo_phyltree$tips_clusters$cluster_2]<-17
vpch[phylo_phyltree$tips_clusters$cluster_3]<-19
plot(OUOUdata[,1],OUOUdata[,2],pch=vpch,cex=1.25,main="",xlab="trait 1",ylab="trait 2",cex.lab=1.5,font=2)
dev.off()

