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

S<-read.csv("~/Dropbox/Flightless_project/Simultaneous_species.csv",header=FALSE)
simultaneous<-S$V1


Burleigh<-read.tree("~/Dropbox/Flightless_project/Trees/All_Birds/Burleigh/Burleigh.tre")
JG<-read.nexus("~/Dropbox/Flightless_project/Trees/All_Birds/Jetz_gene.nex")
JGF<-read.nexus("/Users/ryanterrill/Dropbox/Flightless_project/Trees/All_Birds/JGF.nexus")

JA<-read.nexus("~/Dropbox/Flightless_project/Trees/All_Birds/Jetz_all.nex")
JAF<-read.nexus("~/Dropbox/Flightless_project/Trees/All_Birds/Jetz_all_Flightless.nexus")

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



#Equal rates models
BurleighFlightLambda_ER<-llply(BurBI,fitDiscrete,BurFlight,transform="lambda",model="ER".progress="text")
JGFlightLambda_ER<-llply(JGBI,fitDiscrete,JGFlight,transform="lambda",model="ER".progress="text")
JGFFlightLambda_ER<-llply(JGFBI,fitDiscrete,JGFFlight,transform="lambda",model="ER".progress="text")
JAFlightLambda_ER<-llply(JABI,fitDiscrete,JAFlight,transform="lambda",model="ER".progress="text")
JAFFlightLambda_ER<-llply(JAFBI,fitDiscrete,JAFFlight,transform="lambda",model="ER".progress="text")

#ARD Models
BurleighFlightLambda_ARD<-llply(BurBI,fitDiscrete,BurFlight,transform="lambda",model="ARD",.progress="text")
JGFlightLambda_ARD<-llply(JGBI,fitDiscrete,JGFlight,transform="lambda",model="ARD",.progress="text")
JGFFlightLambda_ARD<-llply(JGFBI,fitDiscrete,JGFFlight,transform="lambda",model="ARD",.progress="text")
JAFlightLambda_ARD<-llply(JABI,fitDiscrete,JAFlight,transform="lambda",model="ARD",.progress="text")
JAFFlightLambda_ARD<-llply(JAFBI,fitDiscrete,JAFFlight,transform="lambda",model="ARD",.progress="text")


#Check for and save P values

##Also need to etxract AIC and P values for ARD models
##IF ARD is better, USE THESE PARAMETERS, not ER

#Burleigh - lambda for flightlessness, All Rates Different model

BFLam_ER<-matrix(nrow=length(BurleighFlightLambda_ER),ncol=2)

for(i in 1:length(BurleighFlightLambda_ER)){
	BFLam_ER[i,1]<-BurleighFlightLambda_ER[[i]]$opt$lam
	BFLam_ER[i,2]<-BurleighFlightLambda_ER[[i]]$opt$aic
	}
colnames(BFLam_ER)<-c("Lambda","AIC")

write.csv(BFLam_ER,"~/Dropbox/Flightless_project/Results/BurleighLambda_ER.csv")

#Burleigh - lambda for flightlessness, All Rates Different model

BFLam_ARD<-matrix(nrow=length(BurleighFlightLambda_ARD),ncol=2)

for(i in 1:length(BurleighFlightLambda_ARD)){
	BFLam_ARD[i,1]<-BurleighFlightLambda_ARD[[i]]$opt$lam
	BFLam_ARD[i,2]<-BurleighFlightLambda_ARD[[i]]$opt$aic
	}
colnames(BFLam_ARD)<-c("Lambda","AIC")

write.csv(BFLam_ARD,"~/Dropbox/Flightless_project/Results/BurleighLambda_ARD.csv")


#Jetz Gene-only - lambda for flightlessness,Equal Rates model

JGFLam_ER<-matrix(nrow=length(JGFlightLambda_ER),ncol=2)

for(i in 1:length(JGFlightLambda_ER)){
	JGFLam_ER[i,1]<-JGFlightLambda_ER[[i]]$opt$lam
	JGFLam_ER[i,2]<-JGFlightLambda_ER[[i]]$opt$aic
	}
colnames(JGFLam_ER)<-c("Lambda","AIC")

write.csv(JGFLam_ER,"~/Dropbox/Flightless_project/Results/JGLambda_ER.csv")


#Jetz Gene-only - lambda for flightlessness, All Rates Different model
JGFLam_ARD<-matrix(nrow=length(JGFlightLambda_ARD),ncol=2)

for(i in 1:length(JGFlightLambda_ARD)){
	JGFLam_ARD[i,1]<-JGFlightLambda_ARD[[i]]$opt$lam
	JGFLam_ARD[i,2]<-JGFlightLambda_ARD[[i]]$opt$aic
	}
colnames(JGFLam_ARD)<-c("Lambda","AIC")

write.csv(JGFLam_ARD,"~/Dropbox/Flightless_project/Results/JGLambda_ARD.csv")




#Jetz Gene-only + flightless - lambda for flightlessness, Equal Rates model

JGFFLam_ER<-matrix(nrow=length(JGFFlightLambda_ER),ncol=2)

for(i in 1:length(JGFFlightLambda_ER)){
	JGFFLam_ER[i,1]<-JGFFlightLambda_ER[[i]]$opt$lam
	JGFFLam_ER[i,2]<-JGFFlightLambda_ER[[i]]$opt$aic
	}
colnames(JGFFLam_ER)<-c("Lambda","AIC")

write.csv(JGFFLam_ER,"~/Dropbox/Flightless_project/Results/JGFLambda_ER.csv")


#Jetz Gene-only + flightless - lambda for flightlessness, All Rates Different model

JGFFLam_ARD<-matrix(nrow=length(JGFFlightLambda_ARD),ncol=2)

for(i in 1:length(JGFFlightLambda_ARD)){
	JGFFLam_ARD[i,1]<-JGFFlightLambda_ARD[[i]]$opt$lam
	JGFFLam_ARD[i,2]<-JGFFlightLambda_ARD[[i]]$opt$aic
	}
colnames(JGFFLam_ARD)<-c("Lambda","AIC")

write.csv(JGFFLam_ARD,"~/Dropbox/Flightless_project/Results/JGFLambda_ARD.csv")
 

#Jetz all - lambda for flightlessness, Equal rates model

JAFLam_ER<-matrix(nrow=length(JAFlightLambda_ER),ncol=2)

for(i in 1:length(JAFlightLambda_ER)){
	JAFLam_ER[i,1]<-JAFlightLambda_ER[[i]]$opt$lam
	JAFLam_ER[i,2]<-JAFlightLambda_ER[[i]]$opt$aic
	}
colnames(JAFLam_ER)<-c("Lambda","AIC")

write.csv(JAFLam_ER,"~/Dropbox/Flightless_project/Results/JALambda_ER.csv")

#Jetz all - lambda for flightlessness, All Rates Different model

JAFLam_ER<-matrix(nrow=length(JAFlightLambda_ER),ncol=2)

for(i in 1:length(JAFlightLambda_ER)){
	JAFLam_ER[i,1]<-JAFlightLambda_ER[[i]]$opt$lam
	JAFLam_ER[i,2]<-JAFlightLambda_ER[[i]]$opt$aic
	}
colnames(JAFLam_ER)<-c("Lambda","AIC")

write.csv(JAFLam_ER,"~/Dropbox/Flightless_project/Results/JALambda_ER.csv")


#Jetz all + flightless- lambda for flightlessness, Equal Rates model

JAFFLam<-matrix(nrow=length(JAFFlightLambda),ncol=2)

for(i in 1:length(JAFFlightLambda)){
	JAFFLam[i,1]<-JAFFlightLambda[[i]]$opt$lam
	JAFFLam[i,2]<-JAFFlightLambda[[i]]$opt$aic
	}
colnames(JAFFLam)<-c("Lambda","AIC")

write.csv(JAFFLam,"~/Dropbox/Flightless_project/Results/JAFLambda_ER.csv")

#Jetz all + flightless - lambda for flightlessness, All Rates Different model

JAFFLam_ARD<-matrix(nrow=length(JAFFlightLambda_ARD),ncol=2)

for(i in 1:length(JAFFlightLambda_ARD)){
	JAFFLam[i,1]<-JAFFlightLambda_ARD[[i]]$opt$lam
	JAFFLam[i,2]<-JAFFlightLambda_ARD[[i]]$opt$aic
	}
colnames(JAFFLam_ARD)<-c("Lambda","AIC")

write.csv(JAFFLam_ARD,"~/Dropbox/Flightless_project/Results/JAFLambda_ARD.csv")

# Save these on an external hard drive

save(BurleighFlightLambda,file="~/Dropbox/Flightless_project/Results/fitDiscrete_Results/Burleigh_lambda.RData")
save(JGFlightLambda,file="~/Dropbox/Flightless_project/Results/fitDiscrete_Results/Jetz_Gene_lambda.RData")
save(JGFFlightLambda,file="~/Dropbox/Flightless_project/Results/fitDiscrete_Results/Jetz_Gene_Flightless_lambda.RData")
save(JAFlightLambda,file="~/Dropbox/Flightless_project/Results/fitDiscrete_Results/Jetz_All_lambda.RData")
save(JAFFlightLambda,file="~/Dropbox/Flightless_project/Results/fitDiscrete_Results/Jetz_All_Flightless_lambda.RData")

n<-rep(NA,(1000-101))
burLambda<-c(BFLam[,2],n)
treeLambdas<-data.frame(burLambda,JGFLam[,2],JGFFLam[,1],JAFLam[,2],JAFFLam[,2])
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




