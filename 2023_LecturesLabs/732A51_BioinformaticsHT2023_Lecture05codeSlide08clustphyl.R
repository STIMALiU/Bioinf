## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this code or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .


## This code replicates the idea behind Figures 5,6,7 in
## J. Felsenstein, 1985. Phylogenies and the comparative method. Am. Nat. 125:1–15.
## Requires mvsLOUCH >= 2.7.6 to run correctly

library(mvSLOUCH)
library(ape)

num_traits<-2

f_rstartree<-function(n,branch_len=1){
    phyltree<-stree(n)
    phyltree$edge.length<-rep(branch_len,nrow(phyltree$edge))
    phyltree
}

v_clustsize<-c(20,20)
joining_branchlengths<-c(0.3,NA)
phylo_phyltree<-mvSLOUCH::simulate_clustered_phylogeny(v_clustsize,joining_branchlengths=joining_branchlengths,f_simclustphyl=f_rstartree,joiningphyl="sim.bd.taxa_Yule1",b_change_joining_branches=TRUE)

v_regimes<-rep("reg0",sum(sapply(phylo_phyltree$edges_clusters,function(x){length(x)},simplify=TRUE)))
v_regimes[phylo_phyltree$edges_clusters[[2]]]<-"reg1"
v_regimes[phylo_phyltree$edges_clusters[[3]]]<-"reg2"
mPsi<-cbind(c(0,0),c(1,1),c(1.75,1.75))
colnames(mPsi)<-c("reg0","reg1","reg2")
OUOUparameters<-list(vY0=matrix(0,nrow=num_traits,ncol=1),A=diag(2,2,2),mPsi=mPsi,Syy=diag(0.5,2,2))
OUOUdata<-mvSLOUCH::simulOUCHProcPhylTree(phylo_phyltree,OUOUparameters,regimes=v_regimes)

png("TwoClustphyltree.png")
ape::plot.phylo(phylo_phyltree,show.tip.label=FALSE,type="cladogram",direction="leftwards",edge.width=3,x.lim=c(-0.05,1.3))
points(-0.05,10,pch=15,cex=2)
points(-0.05,30,pch=17,cex=2)
dev.off()
png("DataClustsHidden.png")
plot(OUOUdata[,1],OUOUdata[,2],pch=19,cex=1,main="",xlab="",ylab="")
regmodel<-lm(OUOUdata[,2]~OUOUdata[,1])
abline(regmodel,lwd=2)
legend("bottomright",legend=paste0("y=",round(regmodel$coefficients[[2]],3),"x","+",round(regmodel$coefficients[[1]],3)),cex=1.5,bty="n")
dev.off()
png("DataClustsShown.png")
vpch<-rep(15,length(phylo_phyltree$tip.label))
vpch[phylo_phyltree$tips_clusters$cluster_2]<-17
plot(OUOUdata[,1],OUOUdata[,2],pch=vpch,cex=1,main="",xlab="",ylab="")
regmodel1<-lm(OUOUdata[phylo_phyltree$tips_clusters$cluster_1,2]~OUOUdata[phylo_phyltree$tips_clusters$cluster_1,1])
regmodel2<-lm(OUOUdata[phylo_phyltree$tips_clusters$cluster_2,2]~OUOUdata[phylo_phyltree$tips_clusters$cluster_2,1])
abline(regmodel1,lwd=2)
abline(regmodel2,lwd=2)
legend("bottomright",legend=c(paste0("y=",round(regmodel1$coefficients[[2]],3),"x","+",round(regmodel1$coefficients[[1]],3)),paste0("y=",round(regmodel2$coefficients[[2]],3),"x","+",round(regmodel2$coefficients[[1]],3))),cex=1.5,bty="n",pch=c(15,17))
dev.off()

