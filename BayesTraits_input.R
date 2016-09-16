### Import BayesTraits output


bBT<-read.csv("~/Dropbox/Flightless_project/Trees/All_Birds/BurleighDepOutput.csv")
jgBT<-read.csv("~/Dropbox/Flightless_project/Trees/All_Birds/JetzGeneOutput.csv")
jaBT<-read.csv("~/Dropbox/Flightless_project/Trees/All_Birds/JetzAllOutput.csv")
jgfBT<-read.csv("/Users/ryanterrill/Dropbox/Flightless_project/Trees/All_Birds/JGF_Output.csv")
jafBT<-read.csv("~/Dropbox/Flightless_project/Trees/All_Birds/JetzAllFlightlessOutput.csv")


#get and arrange character strings


model.string.table<-function(BT,tree_name){
tmp<-BT$Model.string
tab<-table(tmp)
Tab<-data.frame(tab,tree=tree_name)
names(Tab)<-c("model.string","count","tree")
Tab<-Tab[order(-Tab$count),]
Tab<-Tab[1:10,]
out=Tab
}




Bur.ms<-model.string.table(bBT,"Burleigh")
JG.ms<-model.string.table(jgBT,"JetzGene")
JA.ms<-model.string.table(jaBT,"JetzAll")
JGF.ms<-model.string.table(jgfBT,"JetzGeneFlightless")
JAF.ms<-model.string.table(jafBT,"JetzAllFlightless")

model.strings<-rbind(Bur.ms,JG.ms,JA.ms,JGF.ms,JAF.ms)


#save a table with occurence of top 10 model strings
write.csv(model.strings,"~/Dropbox/Flightless_project/Results/Model.strings.csv")

mod.string<-read.csv("~/Dropbox/Flightless_project/Results/Model.strings.csv",stringsAsFactors=FALSE)

####need to run BayesTraits with independent constraint to calculate Bayes Facotr support for dependent


#split up model strings to visualize individual transition models
library(stringr)

ms.split<-str_split_fixed(mod.string$model.string," ",8)
trans.13<-data.frame(as.character(ms.split[,2]),mod.string$count,mod.string$tree,transition="Q13");names(trans.13)<-c("model","count","tree","transition")
trans.24<-data.frame(as.character(ms.split[,4]),mod.string$count,mod.string$tree,transition="Q24");names(trans.24)<-c("model","count","tree","transition")


trans.mod<-function(BT,tree_name){
	mod.string<-BT$Model.string
	ms.split<-str_split_fixed(mod.string," ",8)
	q13<-as.character(ms.split[,2])
	q24<-as.character(ms.split[,4])
	out<-data.frame(q13,q24,tree=tree_name)
	return(out)	
}

Bur.trans.mod<-trans.mod(bBT,"Burleigh")
JG.trans.mod<-trans.mod(jgBT,"JetzGene")
JA.trans.mod<-trans.mod(jaBT,"JetzAll")
JGF.trans.mod<-trans.mod(jgfBT,"JetzGeneFlightless")
JAF.trans.mod<-trans.mod(jafBT,"JetzAllFlightless")

transition.models<-rbind(Bur.trans.mod,JG.trans.mod,JA.trans.mod,JGF.trans.mod,JAF.trans.mod)

write.csv(transition.models,"~/Dropbox/Flightless_project/Results/transition.models.csv")

# get and arrange transition rates
# loss of flight preceded by sequential molt: Q13
# loss of flight preceded by simultaneous molt: Q24


b.13<-data.frame(bBT$q13,rep("q13",nrow(bBT)),rep("Burleigh",nrow(bBT)));names(b.13)<-c("rate","parameter","tree")
b.24<-data.frame(bBT$q24,rep("q24",nrow(bBT)),rep("Burleigh",nrow(bBT)));names(b.24)<-c("rate","parameter","tree")
JG.13<-data.frame(jgBT$q13,rep("q13",nrow(jgBT)),rep("Jetz Gene",nrow(jgBT)));names(JG.13)<-c("rate","parameter","tree")
JG.24<-data.frame(jgBT$q24,rep("q24",nrow(jgBT)),rep("Jetz Gene",nrow(jgBT)));names(JG.24)<-c("rate","parameter","tree")
JA.13<-data.frame(jaBT$q13,rep("q13",nrow(jaBT)),rep("Jetz All",nrow(jaBT)));names(JA.13)<-c("rate","parameter","tree")
JA.24<-data.frame(jaBT$q24,rep("q24",nrow(jaBT)),rep("Jetz All",nrow(jaBT)));names(JA.24)<-c("rate","parameter","tree")
JGF.13<-data.frame(jgfBT$q13,rep("q13",nrow(jgfBT)),rep("Jetz Gene + Flightless",nrow(jgfBT)));names(JGF.13)<-c("rate","parameter","tree")
JGF.24<-data.frame(jgfBT$q24,rep("q24",nrow(jgfBT)),rep("Jetz Gene + Flightless",nrow(jgfBT)));names(JGF.24)<-c("rate","parameter","tree")
JAF.13<-data.frame(jafBT$q13,rep("q13",nrow(jafBT)),rep("Jetz All + Flightless",nrow(jafBT)));names(JAF.13)<-c("rate","parameter","tree")
JAF.24<-data.frame(jafBT$q24,rep("q24",nrow(jafBT)),rep("Jetz All + Flightless",nrow(jafBT)));names(JAF.24)<-c("rate","parameter","tree")




trans<-rbind(b.13,b.24,JG.13,JG.24,JA.13,JA.24,JGF.13,JGF.24,JAF.13,JAF.24)


write.csv(trans,file="~/Dropbox/Flightless_project/Trees/All_Birds/BayesTraitsTransitions.csv")


##plot transition rates

trans$tree<-factor(trans$tree, levels=c("Burleigh","Jetz Gene","Jetz Gene + Flightless","Jetz All","Jetz All + Flightless"))

ggplot(data=trans,aes(x=log(rate),fill=parameter))+geom_histogram(bins=100)+facet_wrap(~tree,scales="free",ncol=1)+theme_bw()

ggsave("/Users/ryanterrill/Dropbox/Flightless_project/Results/posteriorDensitiesColor.pdf")

ggplot(data=trans,aes(x=log(rate),fill=parameter))+geom_histogram(bins=100)+facet_wrap(~tree,scales="free",ncol=1)+theme_bw()+scale_fill_manual(values=c("gray","black"))+theme(strip.background=element_rect(fill="white",color="white"))+theme(strip.text.x=element_text(size=12))+theme(panel.border=element_rect(colour="transparent"))+theme(axis.line.y=element_line(colour="black",size=.5,linetype="solid"))+theme(axis.line.x=element_line(colour="black",size=.5,linetype="solid"))+theme(panel.grid.major = element_line(colour = "transparent"))

ggsave("/Users/ryanterrill/Dropbox/Flightless_project/Manuscript/posteriorDensitiesBW.pdf")
ggsave("/Users/ryanterrill/Dropbox/Flightless_project/Manuscript/posteriorDensitiesBW.png")

ggplot(data=trans,aes(x=parameter,y=rate))+geom_boxplot()+theme_bw()+facet_wrap(~tree,scales="free_y",nrow=1)+theme_bw()








