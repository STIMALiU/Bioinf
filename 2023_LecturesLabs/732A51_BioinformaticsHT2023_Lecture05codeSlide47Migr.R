## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this code or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .


## This code replicates the idea behind Figure 2 in
## K. Bartoszek, S. Glemin, I. Kaj, and M. Lascoux, 2017,
## The Ornstein–Uhlenbeck process with migration: evolution with interactions, 
## J. Theor. Biol 429:35-45.

library(pcmabc)


phyltree<-list("edge"=rbind(c(4,1),c(4,5),c(5,2),c(5,3)),"tip.label"=c("t1","t2","t3"),"edge.length"=c(0.1,0.05,0.05,0.05), "Nnode"=2)
class(phyltree)<-"phylo"


## define a variation around the YUIMA model for OU with different optima 
## for different branches
simulate_mvsl_sde<-function(time,params,X0,step){
	 maxval<-switch(X0[2],3,10)	    
	 U<-runif(1,min=0,max=maxval)
	 A<-paste0("(-1.5*(x1-",U,"))")
         yuima.OU <- yuima::setModel(drift = A, diffusion = 0.1, state.variable=c("x1"),solve.variable=c("x1") )
         res<-pcmabc::simulate_sde_on_branch(time,yuima.OU,X0[-2],step)
         rbind(res,X0[2]+1)
}
     
step<-0.00001 
sde.params=list(s11=0.01)
simres<-pcmabc::simulate_phenotype_on_tree(phyltree, simulate_mvsl_sde, simul.params=NULL, X0=c(0,1), step)

simres2<-simres
simres2$root.branch.phenotype<-simres2$root.branch.phenotype[-3,,drop=FALSE]
simres2$phenotype<-sapply(simres2$phenotype,function(x){x[-3,,drop=FALSE]},simplify=FALSE)

png("MigrationPhyl.png")
pcmabc::draw_phylproc(simres2)
num_lines<-40
v_sample_index<-1:(min(ncol(simres2$phenotype[[1]]),ncol(simres2$phenotype[[2]])))
n<-length(v_sample_index)
step_n<-as.integer(n/num_lines)
v_sample_index<-seq(from=1,by=step_n,length.out=num_lines)
apply(rbind(simres2$phenotype[[1]][,v_sample_index],simres2$phenotype[[2]][,v_sample_index]),2,function(x){
    lines(list(y=c(x[1],x[3]),x=c(x[2],x[4])))
})


num_lines<-20
v_sample_index<-1:(as.integer(min(ncol(simres2$phenotype[[3]]),ncol(simres2$phenotype[[4]])))/2)
n<-length(v_sample_index)
step_n<-as.integer(n/num_lines)
v_sample_index<-seq(from=1,by=step_n,length.out=num_lines)
m_migration<-rbind(simres2$phenotype[[3]][,v_sample_index],simres2$phenotype[[4]][,v_sample_index])
m_migration[c(1,3),]<-m_migration[c(1,3),]+0.05
apply(m_migration,2,function(x){
    lines(list(y=c(x[1],x[3]),x=c(x[2],x[4])))
})
dev.off()



