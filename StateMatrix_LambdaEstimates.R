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
BurleighFlightLambda_ER<-llply(BurBI,fitDiscrete,BurFlight,transform="lambda",model="ER",.progress="text")
JGFlightLambda_ER<-llply(JGBI,fitDiscrete,JGFlight,transform="lambda",model="ER",.progress="text")
JGFFlightLambda_ER<-llply(JGFBI,fitDiscrete,JGFFlight,transform="lambda",model="ER",.progress="text")
JAFlightLambda_ER<-llply(JABI,fitDiscrete,JAFlight,transform="lambda",model="ER",.progress="text")
JAFFlightLambda_ER<-llply(JAFBI,fitDiscrete,JAFFlight,transform="lambda",model="ER",.progress="text")

#ARD Models
BurleighFlightLambda_ARD<-llply(BurBI,fitDiscrete,BurFlight,transform="lambda",model="ARD",.progress="text")
JGFlightLambda_ARD<-llply(JGBI,fitDiscrete,JGFlight,transform="lambda",model="ARD",.progress="text")
JGFFlightLambda_ARD<-llply(JGFBI,fitDiscrete,JGFFlight,transform="lambda",model="ARD",.progress="text")
JAFlightLambda_ARD<-llply(JABI,fitDiscrete,JAFlight,transform="lambda",model="ARD",.progress="text")
JAFFlightLambda_ARD<-llply(JAFBI,fitDiscrete,JAFFlight,transform="lambda",model="ARD",.progress="text")


#Save lambda and AIC values for each model for each tree 


BFLam_ER<-matrix(nrow=length(BurleighFlightLambda_ER),ncol=5)

for(i in 1:length(BurleighFlightLambda_ER)){
	BFLam_ER[i,1]<-BurleighFlightLambda_ER[[i]]$opt$lam
	BFLam_ER[i,2]<-BurleighFlightLambda_ER[[i]]$opt$aic
	BFLam_ER[i,3]<-BurleighFlightLambda_ER[[i]]$opt$q12
	BFLam_ER[i,4]<-BurleighFlightLambda_ER[[i]]$opt$q21
	BFLam_ER[i,5]<-BurleighFlightLambda_ER[[i]]$opt$lnL
	}
colnames(BFLam_ER)<-c("Lambda","AIC")

write.csv(BFLam_ER,"~/Dropbox/Flightless_project/Results/BurleighLambda_ER.csv")

#Burleigh - lambda for flightlessness, All Rates Different model

BFLam_ARD<-matrix(nrow=length(BurleighFlightLambda_ARD),ncol=5)

for(i in 1:length(BurleighFlightLambda_ARD)){
	BFLam_ARD[i,1]<-BurleighFlightLambda_ARD[[i]]$opt$lam
	BFLam_ARD[i,2]<-BurleighFlightLambda_ARD[[i]]$opt$aic
	BFLam_ARD[i,3]<-BurleighFlightLambda_ARD[[i]]$opt$q12
	BFLam_ARD[i,4]<-BurleighFlightLambda_ARD[[i]]$opt$q21
	BFLam_ARD[i,5]<-BurleighFlightLambda_ARD[[i]]$opt$lnL
	}
colnames(BFLam_ARD)<-c("Lambda","AIC")

write.csv(BFLam_ARD,"~/Dropbox/Flightless_project/Results/BurleighLambda_ARD.csv")


#Jetz Gene-only - lambda for flightlessness,Equal Rates model

JGFLam_ER<-matrix(nrow=length(JGFlightLambda_ER),ncol=5)

for(i in 1:length(JGFlightLambda_ER)){
	JGFLam_ER[i,1]<-JGFlightLambda_ER[[i]]$opt$lam
	JGFLam_ER[i,2]<-JGFlightLambda_ER[[i]]$opt$aic
	JGFLam_ER[i,3]<-JGFlightLambda_ER[[i]]$opt$q12
	JGFLam_ER[i,4]<-JGFlightLambda_ER[[i]]$opt$q21
	JGFLam_ER[i,5]<-JGFlightLambda_ER[[i]]$opt$lnL
	}
colnames(JGFLam_ER)<-c("Lambda","AIC")

write.csv(JGFLam_ER,"~/Dropbox/Flightless_project/Results/JGLambda_ER.csv")


#Jetz Gene-only - lambda for flightlessness, All Rates Different model
JGFLam_ARD<-matrix(nrow=length(JGFlightLambda_ARD),ncol=5)

for(i in 1:length(JGFlightLambda_ARD)){
	JGFLam_ARD[i,1]<-JGFlightLambda_ARD[[i]]$opt$lam
	JGFLam_ARD[i,2]<-JGFlightLambda_ARD[[i]]$opt$aic
	JGFLam_ARD[i,3]<-JGFlightLambda_ARD[[i]]$opt$q12
	JGFLam_ARD[i,4]<-JGFlightLambda_ARD[[i]]$opt$q21
	JGFLam_ARD[i,5]<-JGFlightLambda_ARD[[i]]$opt$lnL
	}
colnames(JGFLam_ARD)<-c("Lambda","AIC")

write.csv(JGFLam_ARD,"~/Dropbox/Flightless_project/Results/JGLambda_ARD.csv")




#Jetz Gene-only + flightless - lambda for flightlessness, Equal Rates model

JGFFLam_ER<-matrix(nrow=length(JGFFlightLambda_ER),ncol=5)

for(i in 1:length(JGFFlightLambda_ER)){
	JGFFLam_ER[i,1]<-JGFFlightLambda_ER[[i]]$opt$lam
	JGFFLam_ER[i,2]<-JGFFlightLambda_ER[[i]]$opt$aic
	JGFFLam_ER[i,3]<-JGFFlightLambda_ER[[i]]$opt$q12
	JGFFLam_ER[i,4]<-JGFFlightLambda_ER[[i]]$opt$q21
	JGFFLam_ER[i,5]<-JGFFlightLambda_ER[[i]]$opt$lnL
	}
colnames(JGFFLam_ER)<-c("Lambda","AIC")

write.csv(JGFFLam_ER,"~/Dropbox/Flightless_project/Results/JGFLambda_ER.csv")


#Jetz Gene-only + flightless - lambda for flightlessness, All Rates Different model

JGFFLam_ARD<-matrix(nrow=length(JGFFlightLambda_ARD),ncol=5)

for(i in 1:length(JGFFlightLambda_ARD)){
	JGFFLam_ARD[i,1]<-JGFFlightLambda_ARD[[i]]$opt$lam
	JGFFLam_ARD[i,2]<-JGFFlightLambda_ARD[[i]]$opt$aic
	JGFFLam_ARD[i,3]<-JGFFlightLambda_ARD[[i]]$opt$q12
	JGFFLam_ARD[i,4]<-JGFFlightLambda_ARD[[i]]$opt$q21
	JGFFLam_ARD[i,5]<-JGFFlightLambda_ARD[[i]]$opt$lnL
	}
colnames(JGFFLam_ARD)<-c("Lambda","AIC")

write.csv(JGFFLam_ARD,"~/Dropbox/Flightless_project/Results/JGFLambda_ARD.csv")
 

#Jetz all - lambda for flightlessness, Equal rates model

JAFLam_ER<-matrix(nrow=length(JAFlightLambda_ER),ncol=5)

for(i in 1:length(JAFlightLambda_ER)){
	JAFLam_ER[i,1]<-JAFlightLambda_ER[[i]]$opt$lam
	JAFLam_ER[i,2]<-JAFlightLambda_ER[[i]]$opt$aic
	JAFLam_ER[i,3]<-JAFlightLambda_ER[[i]]$opt$q12
	JAFLam_ER[i,4]<-JAFlightLambda_ER[[i]]$opt$q21
	JAFLam_ER[i,5]<-JAFlightLambda_ER[[i]]$opt$lnL
	}
colnames(JAFLam_ER)<-c("Lambda","AIC")

write.csv(JAFLam_ER,"~/Dropbox/Flightless_project/Results/JALambda_ER.csv")

#Jetz all - lambda for flightlessness, All Rates Different model

JAFLam_ARD<-matrix(nrow=length(JAFlightLambda_ARD),ncol=5)

for(i in 1:length(JAFlightLambda_ARD)){
	JAFLam_ARD[i,1]<-JAFlightLambda_ARD[[i]]$opt$lam
	JAFLam_ARD[i,2]<-JAFlightLambda_ARD[[i]]$opt$aic
	JAFLam_ARD[i,3]<-JAFlightLambda_ARD[[i]]$opt$q12
	JAFLam_ARD[i,4]<-JAFlightLambda_ARD[[i]]$opt$q21
	JAFLam_ARD[i,5]<-JAFlightLambda_ARD[[i]]$opt$lnL
	
	}
colnames(JAFLam_ARD)<-c("Lambda","AIC")

write.csv(JAFLam_ARD,"~/Dropbox/Flightless_project/Results/JALambda_ARD.csv")


#Jetz all + flightless- lambda for flightlessness, Equal Rates model

JAFFLam_ER<-matrix(nrow=length(JAFFlightLambda_ER),ncol=5)

for(i in 1:length(JAFFlightLambda_ER)){
	JAFFLam_ER[i,1]<-JAFFlightLambda_ER[[i]]$opt$lam
	JAFFLam_ER[i,2]<-JAFFlightLambda_ER[[i]]$opt$aic
	JAFFLam_ER[i,3]<-JAFFlightLambda_ER[[i]]$opt$q12
	JAFFLam_ER[i,4]<-JAFFlightLambda_ER[[i]]$opt$q21
	JAFFLam_ER[i,5]<-JAFFlightLambda_ER[[i]]$opt$lnL
	}
colnames(JAFFLam_ER)<-c("Lambda","AIC")

write.csv(JAFFLam_ER,"~/Dropbox/Flightless_project/Results/JAFLambda_ER.csv")

#Jetz all + flightless - lambda for flightlessness, All Rates Different model
#Add all these parameters to the others

JAFFLam_ARD<-matrix(nrow=length(JAFFlightLambda_ARD),ncol=5)

for(i in 1:length(JAFFlightLambda_ARD)){
	JAFFLam_ARD[i,1]<-JAFFlightLambda_ARD[[i]]$opt$lam
	JAFFLam_ARD[i,2]<-JAFFlightLambda_ARD[[i]]$opt$aic
	JAFFLam_ARD[i,3]<-JAFFlightLambda_ARD[[i]]$opt$q12
	JAFFLam_ARD[i,4]<-JAFFlightLambda_ARD[[i]]$opt$q21
	JAFFLam_ARD[i,5]<-JAFFlightLambda_ARD[[i]]$opt$lnL
	}
colnames(JAFFLam_ARD)<-c("Lambda","AIC","q12","q21","lnL")

write.csv(JAFFLam_ARD,"~/Dropbox/Flightless_project/Results/JAFLambda_ARD.csv")


### Load transition parameters from files

B_ER<-read.csv("/Users/ryanterrill/Dropbox/Flightless_project/Results/BurleighLambda_ER.csv")
B_ARD<-read.csv("/Users/ryanterrill/Dropbox/Flightless_project/Results/BurleighLambda_ARD.csv")
JG_ER<-read.csv("/Users/ryanterrill/Dropbox/Flightless_project/Results/JGLambda_ER.csv")
JG_ARD<-read.csv("/Users/ryanterrill/Dropbox/Flightless_project/Results/JGLambda_ARD.csv")
JGF_ER<-read.csv("/Users/ryanterrill/Dropbox/Flightless_project/Results/JGFLambda_ER.csv")
JGF_ARD<-read.csv("/Users/ryanterrill/Dropbox/Flightless_project/Results/JGFLambda_ARD.csv")
JA_ER<-read.csv("/Users/ryanterrill/Dropbox/Flightless_project/Results/JALambda_ER.csv")
JA_ARD<-read.csv("/Users/ryanterrill/Dropbox/Flightless_project/Results/JALambda_ARD.csv")
JAF_ER<-read.csv("/Users/ryanterrill/Dropbox/Flightless_project/Results/JAFLambda_ER.csv")
JAF_ARD<-read.csv("/Users/ryanterrill/Dropbox/Flightless_project/Results/JAFLambda_ARD.csv")

### Compare ER to ARD models from AIC - present lambda from best (ARD) model in paper - but give lambda for ER and ARD, as well as AIC for each in supplemental
##Need to add q12, q21, and logLik to these matrices

Burleigh_ER<-data.frame(B_ER$Lambda,B_ER$AIC,"Burleigh","ER");names(Burleigh_ER)<-c("Lambda","AIC","Tree","Model")
Burleigh_ARD<-data.frame(B_ARD$Lambda,B_ARD$AIC,"Burleigh","ARD");names(Burleigh_ARD)<-c("Lambda","AIC","Tree","Model")
JetzGene_ER<-data.frame(JG_ER$Lambda,JG_ER$AIC,"Jetz_gene","ER");names(JetzGene_ER)<-c("Lambda","AIC","Tree","Model")
JetzGene_ARD<-data.frame(JG_ARD$Lambda,JG_ARD$AIC,"Jetz_gene","ARD");names(JetzGene_ARD)<-c("Lambda","AIC","Tree","Model")
JetzGeneFlightless_ER<-data.frame(JGF_ER$Lambda,JGF_ER$AIC,"Jetz_gene_flightless","ER");names(JetzGeneFlightless_ER)<-c("Lambda","AIC","Tree","Model")
JetzGeneFlightless_ARD<-data.frame(JGF_ARD$Lambda,JGF_ARD$AIC,"Jetz_gene_flightless","ARD");names(JetzGeneFlightless_ARD)<-c("Lambda","AIC","Tree","Model")
JetzAll_ER<-data.frame(JA_ER$Lambda,JA_ER$AIC,"Jetz_All","ER");names(JetzAll_ER)<-c("Lambda","AIC","Tree","Model")
JetzAll_ARD<-data.frame(JA_ARD$Lambda,JA_ARD$AIC,"Jetz_All","ARD");names(JetzAll_ARD)<-c("Lambda","AIC","Tree","Model")
JetzAllFlightless_ER<-data.frame(JAF_ER$Lambda,JAF_ER$AIC,"Jetz_All_flightless","ER");names(JetzAllFlightless_ER)<-c("Lambda","AIC","Tree","Model")
JetzAllFlightless_ARD<-data.frame(JAF_ARD$Lambda,JAF_ARD$AIC,"Jetz_All_flightless","ARD");names(JetzAllFlightless_ARD)<-c("Lambda","AIC","Tree","Model")

Lambdas<-rbind(Burleigh_ER,Burleigh_ARD,JetzGene_ER,JetzGene_ARD,JetzGeneFlightless_ER,JetzGeneFlightless_ARD,JetzAll_ER,JetzAll_ARD,JetzAllFlightless_ER,JetzAllFlightless_ARD)

ERLambdas<-Lambdas[Lambdas$Model=="ER",]
ARDLambdas<-Lambdas[Lambdas$Model=="ARD",]

ggplot(data=Lambdas,aes(y=Lambda,x=Tree,fill=Model))+geom_boxplot()+theme_bw()+scale_fill_grey()

ggplot(data=Lambdas,aes(y=AIC,x=Tree,fill=Model))+geom_boxplot()+theme_bw()+scale_fill_grey()

btest<-rbind(Burleigh_ER,Burleigh_ARD)
btest<-btest[,-btest$Lambda]
index<-1:101
index<-c(index,index)
btest$index<-index
bt.m<-melt(btest,var="Tree","Model","index")
ggplot(data=btest,aes(y=AIC,x=Model))+geom_boxplot(alpha=.1,outlier.shape=NA)+theme_bw()+geom_line(aes(group=index),colour="darkgray",alpha=.5,lwd=.5)+geom_point(cex=.5)

index<-c(1:101,1:101,1:1000,1:1000,1:1000,1:1000,1:1000,1:1000,1:1000,1:1000)
Lambdas$index<-index

#plot AIC values
pdf("~/Dropbox/Flightless_project/Results/AIC_values.pdf")
ggplot(data=Lambdas,aes(y=AIC,x=Model))+facet_wrap(~Tree,scales="free",nrow=1)+theme_bw()+geom_line(aes(group=index),colour="black",alpha=.1,lwd=.015)+geom_point(cex=.5,alpha=.2)+theme(strip.text.x = element_text(size=7))
dev.off()

#plot lambdas under ER
pdf("~/Dropbox/Flightless_project/Results/ERLambdas.pdf")
ggplot(data=ERLambdas,aes(y=Lambda,x=Tree,fill=Model))+geom_boxplot(outlier.shape=NA,fill="transparent")+theme_bw()+geom_jitter(alpha=.1)+scale_fill_grey()
dev.off()

#plot lambdas under ARD
pdf("~/Dropbox/Flightless_project/Results/ARDLambdas.pdf")
ggplot(data=ARDLambdas,aes(y=Lambda,x=Tree,fill=Model))+geom_boxplot(outlier.shape=NA,fill="transparent")+theme_bw()+geom_jitter(alpha=.1)+scale_fill_grey()
dev.off()



##USE AIC values to calculate likelihood ratio for ARD vs ER models

Burleigh_model_likelihood<-exp((B_ARD$AIC-B_ER$AIC)/2)
pdf("/Users/ryanterrill/Dropbox/Flightless_project/Results/Burleigh_Model_Likelihood")
hist(Burleigh_model_likelihood,breaks=50)
dev.off()


JG_ER_likelihood<-exp((JG_ARD$AIC-JG_ER$AIC)/2)
#remove one outlier at 2x10^64
JG_ER_likelihood[JG_ER_likelihood>1000]<-NA
pdf("/Users/ryanterrill/Dropbox/Flightless_project/Results/JG_Model_Likelihood")
hist(JG_ER_likelihood,breaks=100)
dev.off()

JGF_ER_likelihood<-exp((JGF_ARD$AIC-JGF_ER$AIC)/2)
pdf("/Users/ryanterrill/Dropbox/Flightless_project/Results/JGF_Model_Likelihood")
hist(JGF_ER_likelihood,breaks=100)
dev.off()

JA_model_likelihood

JAF_model_likelihood


##Make two boxplots - one with all Lambdas - both models for each tree - one with all AIC values




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




