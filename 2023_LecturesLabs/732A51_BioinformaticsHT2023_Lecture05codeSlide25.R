## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this code or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .


## This code replicates the idea behind Figure 1 in
## K. Bartoszek, S. Glemin, I. Kaj, and M. Lascoux, 2017,
## The Ornstein–Uhlenbeck process with migration: evolution with interactions, 
## J. Theor. Biol 429:35-45.

library(TreeSim)
library(mvSLOUCH)

phyltree<-TreeSim::sim.bd.taxa(n=50, numbsim=1, lambda=1, mu=0)[[1]]

## The below lines are needed as mvSLOUCH was being upgraded
## and in the new version a different tree format is used
## we also want the tree to be of height 1

if (packageVersion("mvSLOUCH")<2){
    library(ouch)
    print("Using old version of mvSLOUCH. Please consider upgrading if version greater or equal than 2.0.0 is available on CRAN.")
    phyltree<-ape2ouch(phyltree,scale=1) ## changing to OUCH format and rescaling tree to height 1
}else{
    library(ape)
    print("Using new version of mvSLOUCH")
    tree_height<-max(ape::node.depth.edgelength(phyltree))
    phyltree$edge.length<-phyltree$edge.length/tree_height ## rescaling tree to height 1
}

### Define SDE parameters to be able to simulate data under the OUOU model.
OUOUparameters<-list(vY0=matrix(0,1,1),A=matrix(9,1,1),mPsi=matrix(2,1,1),Syy=matrix(1,1,1))

### Now simulate the trait data keeping the whole trajectory
OUOUdata<-mvSLOUCH::simulOUCHProcPhylTree(phyltree,OUOUparameters,regimes=NULL,NULL,fullTrajectory=TRUE)

### Define Brownian motion parameters to be able to simulate data 
### under the Brownian motion model.
BMparameters<-list(vX0=matrix(0,1,1),Sxx=matrix(1,1,1))

### Now simulate the trait data keeping the whole trajectory
BMdata<-mvSLOUCH::simulBMProcPhylTree(phyltree,X0=BMparameters$vX0,Sigma=BMparameters$Sxx,fullTrajectory=TRUE)

## And now do the plotting
png("YuleTree50.png")
plot(phyltree,show.tip.label = FALSE,direction="downwards");print("Press enter to continue");readline()
dev.off()
png("YOU50ex.png")
mvSLOUCH::drawPhylProcess(PhylTraitProcess=OUOUdata);print("Press enter to continue");readline()
dev.off()
png("YBM50ex.png")
mvSLOUCH::drawPhylProcess(PhylTraitProcess=BMdata)
dev.off()



