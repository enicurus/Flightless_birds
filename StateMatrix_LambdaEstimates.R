#Coding states from list files


require(phytools)
require(plyr)
require(geiger)
library(reshape)
require(ggplot2)




stateMatrix<-function(tree,flightless,simultaneous){
	tree1<-tree[[1]]
	treeTips<-tree1$tip.label
	flight<-rep(0,length(treeTips))
	molt<-rep(0,length(treeTips))
	flightless_index<-unlist(lapply(flightless,grep,treeTips))
	flight[flightless_index]<-1
	molt_index<-unlist(lapply(simultaneous,grep,treeTips))
	molt[molt_index]<-1
	out<-data.frame(treeTips,flight,molt)
	return(out)
	}
	
	
##Fossil node01 0 Struthioniformes_Struthionidae_Struthio_camelus Passeriformes-oscines_Passeridae_Passer_domesticus

F<-read.csv("~/Dropbox/Flightless_project/Trees/All_Birds/Flightless_taxa.csv")
flightless<-F$Flightless_taxa

S<-read.csv("Dropbox/Flightless_project/Simultaneous_species.csv",header=FALSE)
simultaneous<-S$V1


Burleigh<-read.tree("~/Dropbox/Flightless_project/Trees/All_Birds/Burleigh/Burleigh.tre")
JG<-read.nexus("~/Dropbox/Flightless_project/Trees/All_Birds/Jetz_gene.nex")
JGF<-read.nexus("/Users/ryanterrill/Dropbox/Flightless_project/Trees/All_Birds/JetzGeneFlightless/JGF_final.nexus")

JA<-read.nexus("~/Dropbox/Flightless_project/Trees/All_Birds/Jetz_all.nex")
JAF<-read.nexus("~/Dropbox/Flightless_project/Trees/All_Birds/JetzAllFlightless.nex")

Bstates<-stateMatrix(Burleigh,flightless,simultaneous)
JGstates<-stateMatrix(JG,flightless,simultaneous)
JGFstates<-stateMatrix(JGF,flightless,simultaneous)
JAstates<-stateMatrix(JA,flightless,simultaneous)
JAFstates<-stateMatrix(JAF,flightless,simultaneous)




BurFlight<-matrix(Bstates$flight);rownames(BurFlight)<-Bstates$treeTips
JGFlight<-matrix(JGstates$flight);rownames(JGFlight)<-JGstates$treeTips
JGFFlight<-matrix(JGFstates$flight);rownames(JGFFlight)<-JGFstates$treeTips
JAFlight<-matrix(JAstates$flight);rownames(JAFlight)<-JAstates$treeTips
JAFFlight<-matrix(JAFstates$flight);rownames(JAFFlight)<-JAFstates$treeTips




###run this over the trees and get lambda and aic for all 

BurBI<-llply(Burleigh,multi2di,.progress="text")
JGBI<-llply(JG,multi2di,.progress="text")
JGFBI<-llply(JGF,multi2di,.progress="text")
JABI<-llply(JA,multi2di,.progress="text")
JAFBI<-llply(JAF,multi2di,.progress="text")




BurleighFlightLambda<-llply(BurBI,fitDiscrete,BurFlight,transform="lambda",.progress="text")
JGFlightLambda<-llply(JGBI,fitDiscrete,JGFlight,transform="lambda",.progress="text")
JGFFlightLambda<-llply(JGFBI,fitDiscrete,JGFFlight,transform="lambda",.progress="text")
JAFlightLambda<-llply(JABI,fitDiscrete,JAFlight,transform="lambda",.progress="text")
JAFFlightLambda<-llply(JAFBI,fitDiscrete,JAFFlight,transform="lambda",.progress="text")


BFLam<-matrix(nrow=length(BurleighFlightLambda),ncol=2)

for(i in 1:length(BurleighFlightLambda)){
	BFLam[i,1]<-BurleighFlightLambda[[i]]$opt$lam
	BFLam[i,2]<-BurleighFlightLambda[[i]]$opt$aic
	}
colnames(BFLam)<-c("Lambda","AIC")

write.csv(BFLam,"~/Dropbox/Flightless_project/Results/BurleighLambda.csv")


JGFLam<-matrix(nrow=length(JGFlightLambda),ncol=2)

for(i in 1:length(JGFlightLambda)){
	JGFLam[i,1]<-JGFlightLambda[[i]]$opt$lam
	JGFLam[i,2]<-JGFlightLambda[[i]]$opt$aic
	}
colnames(JGFLam)<-c("Lambda","AIC")

write.csv(JGFLam,"~/Dropbox/Flightless_project/Results/JGLambda.csv")


JGFFLam<-matrix(nrow=length(JGFFlightLambda),ncol=2)

for(i in 1:length(JGFFlightLambda)){
	JGFFLam[i,1]<-JGFFlightLambda[[i]]$opt$lam
	JGFFLam[i,2]<-JGFFlightLambda[[i]]$opt$aic
	}
colnames(JGFFLam)<-c("Lambda","AIC")

write.csv(JGFFLam,"~/Dropbox/Flightless_project/Results/JGFLambda.csv")
 


JAFLam<-matrix(nrow=length(JAFlightLambda),ncol=2)

for(i in 1:length(JAFlightLambda)){
	JAFLam[i,1]<-JAFlightLambda[[i]]$opt$lam
	JAFLam[i,2]<-JAFlightLambda[[i]]$opt$aic
	}
colnames(JAFLam)<-c("Lambda","AIC")

write.csv(JAFLam,"~/Dropbox/Flightless_project/Results/JALambda.csv")



JAFFLam<-matrix(nrow=length(JAFFlightLambda),ncol=2)

for(i in 1:length(JAFFlightLambda)){
	JAFFLam[i,1]<-JAFFlightLambda[[i]]$opt$lam
	JAFFLam[i,2]<-JAFFlightLambda[[i]]$opt$aic
	}
colnames(JAFFLam)<-c("Lambda","AIC")

write.csv(JAFFLam,"~/Dropbox/Flightless_project/Results/JAFLambda.csv")


save(BurleighFlightLambda,file="~/Dropbox/Flightless_project/Results/fitDiscrete_Results/Burleigh_lambda.RData")
save(JGFlightLambda,file="~/Dropbox/Flightless_project/Results/fitDiscrete_Results/Jetz_Gene_lambda.RData")
save(JGFFlightLambda,file="~/Dropbox/Flightless_project/Results/fitDiscrete_Results/Jetz_Gene_Flightless_lambda.RData")
save(JAFlightLambda,file="~/Dropbox/Flightless_project/Results/fitDiscrete_Results/Jetz_All_lambda.RData")
save(JAFFlightLambda,file="~/Dropbox/Flightless_project/Results/fitDiscrete_Results/Jetz_All_Flightless_lambda.RData")

n<-rep(NA,(1000-101))
burLambda<-c(BFLam[,1],n)
treeLambdas<-data.frame(burLambda,JGFLam[,1],JGFFLam[,1],JAFLam[,1],JAFFLam[,1])
colnames(treeLambdas)<-c("Burleigh","Jetz Gene","Jetz Gene +\n Flightless","Jetz All","Jetz All +\n Flightless")
colnames(treeLambdas)<-c("Burleigh","Jetz Gene","Jetz Gene +\n Flightless","Jetz All","Jetz All +\n Flightless")




tl.melt<-melt(treeLambdas)
ggplot(data=tl.melt,aes(x=variable,y=value))+geom_boxplot()+theme_bw()+theme(axis.title.x=element_blank())+ylab("Lambda - Flightlessness")
ggsave("/Users/ryanterrill/Dropbox/Flightless_project/Manuscript/Lambda_trees.pdf")


#### transition parameters
JAFFLam<-matrix(nrow=length(JAFFlightLambda),ncol=2)

BrownPar<-matrix(ncol=5,nrow=1000)
for(i in 1:100){
	BrownPar[i,1]<-BurleighFlightLambda[[i]]$opt$q12
	for (j in 101:1000){BrownPar[i,1]<NA
	}
	}
for(i in 1:1000){
	BrownPar[i,2]<-JGFlightLambda[[i]]$opt$q12
	BrownPar[i,3]<-JGFFlightLambda[[i]]$opt$q12
	BrownPar[i,4]<-JAFlightLambda[[i]]$opt$q12
	BrownPar[i,5]<-JAFFlightLambda[[i]]$opt$q12	
}

colnames(BrownPar)<-c("Burleigh","Jetz Gene","Jetz Gene +\n Flightless","Jetz All","Jetz All +\n Flightless")

bp.melt<-melt(BrownPar)
ggplot(data=bp.melt,aes(x=X2,y=log(value)))+geom_boxplot()+theme_bw()+theme(axis.title.x=element_blank())+ylab("Brownian Rate Parameter - Loss of FLight")
ggsave("/Users/ryanterrill/Dropbox/Flightless_project/Manuscript/Brownian_rate_parameter_trees.pdf")


###Also do this in phylosig to see if results are the same




