graft.replace<-function(backBoneTree,cladeTree){
	require(phytools)
	require(plyr)
		c<-cladeTree[[1]]$tip.label    
		ca<-getMRCA(backBoneTree[[1]],c)
		cat("Calculating donor clade branching times","\n")
		bc<-list()
			for(i in 1:length(cladeTree)){
				bc[[i]]<-branching.times(cladeTree[[i]])[1]
				}
		bc<-unlist(bc)
		bb<-list()
		bbTipNums<-list()
			for(i in 1:length(backBoneTree)){	
				bbTipNums[[i]]<-unlist(lapply(c,grep,backBoneTree[[i]]$tip.label))
				cat("step 1:",i,"of",length(backBoneTree),"\r")
				}	
			cat("\n")
		cat("Calculating receiving clade branching times","\n")
		bbtTips<-list()
		
			for(i in 1:length(backBoneTree)){		
				bb[[i]]<-branching.times(backBoneTree[[i]])
				cat("step 1:",i,"of",length(backBoneTree),"\r")
				}
				cat("\n")
			for(i in 1:length(backBoneTree)){	
				bbtTips[[i]]<-bb[[i]][grep(ca,names(bb[[i]]))]
				cat("step 2:",i,"of",length(backBoneTree),"\r")
				}
				cat("\n")
		cat("Synching branching times","\n")	
			for (i in 1:length(backBoneTree)){
				cladeTree[[i]]$edge.length<-cladeTree[[i]]$edge.length*(bbtTips[[i]]/bc[i])
				
			}	
		o<-list()
		for(i in 1:length(backBoneTree)){
			o[[i]]<-paste(backBoneTree[[i]]$tip.label[bbTipNums[[i]]],"_old")
		
			}
		t<-list()  
			for (i in 1:length(backBoneTree)){ 
				t[[i]]<-backBoneTree[[i]]
				t[[i]]$tip.label[bbTipNums[[i]]]<-o[[i]]
				}
		class(t)<-"multiPhylo"
		cat("Binding trees","\n")
		boundTree<-list()
			for(i in 1:length(backBoneTree)){
				boundTree[[i]]<-bind.tree(t[[i]],cladeTree[[i]],where=ca)
			}
		class(boundTree)<-"multiPhylo"
		cat("Dropping duplicated tips","\n")
		dt<-list()
			for(f in 1:length(backBoneTree)){
				dt[[f]]<-drop.tip(boundTree[[f]],o[[i]])
				cat(i,"of",length(backBoneTree),"\r")
			}
		cat("\n")	
		cat("checking tree quality","\n")	
		dup<-0
		ult<-0
	for(i in 1:length(dt)){
		if (is.ultrametric((dt[[i]]))==FALSE)ult<-1	
		if (is.element(TRUE,duplicated(dt[[i]]$tip.label))==TRUE)dup<-1
		cat(i,"of",length(backBoneTree),"\r")
		}
		cat("\n")
		if(dup==1){cat("Duplicate tips in tree :-(","\n")}
		if(ult==1){cat("Tree being returned is not ultrametric :-(","\n")}
	cat("done.","\n")	
class(dt)<-"multiPhylo"	
return(dt)
}
		
		
		
		
##drop flying non-gene tips from JGF
JG.f<-JGstates[JGstates$flight==0,]
JGF.f<-JGFstates[JGFstates$flight==0,]
JGF.drop<-setdiff(JGF.f$treeTips,JG.f$treeTips)
JGF<-llply(JGF,drop.tip,JGF.drop,.progress="text")
class(JGF)<-"multiPhylo"
write.nexus(JGF,file="/Users/ryanterrill/Dropbox/Flightless_project/Trees/All_Birds/JGF.nexus")

