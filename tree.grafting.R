##Code to bind all subclade trees


### first two functions required for fam.replace

clade.depth<-function(backBoneTree,cladeTree){
	require(phytools)
	cladetips<-cladeTree$tip.label                                                    #get clade tips                                                       
	dropTips<-intersect(backBoneTree$tip.label,cladetips)  	                          #make list of tips to drop from backbone tree so you can leave one tip that is guaranteed to be in the backbone tree
	m<-findMRCA(backBoneTree,dropTips)												#get mrca of clade in backbone tree
	t<-extract.clade(backBoneTree,m)												#extract clade tree form backbone tree
	d<-branching.times(t)[1]
	return(d)
	}




replace.clade<-function(backBoneTree,cladeTree,cladeDepth){
	require(phytools)
	cladetips<-cladeTree$tip.label                                                    #get clade tips                                                       
	dropTips<-intersect(backBoneTree$tip.label,cladetips)  	                          #make list of tips to drop from backbone tree so you can leave one tip that is guaranteed to be in the backbone tree
	depth<-clade.depth(backBoneTree,cladeTree)														#find height of clade tree
	cladeTree$edge.length<-cladeTree$edge.length*depth									#multiply branches of input clade tree by clade tree in backbone
	backBoneTree<-drop.tip(backBoneTree,dropTips[2:length(dropTips)])                 #leave the first tip in the backbone tree
	treeLength<-branching.times(cladeTree)[1]						                  # find total depth of Clade tree
	bindtip<-which(backBoneTree$tip.label==dropTips[1])                               #find the tip number
	backBoneTree$tip.label[bindtip]<-"Clade" 									      #rename the tip so it can be dropped later
	t<-match(which(backBoneTree$tip.label=="Clade"),backBoneTree$edge[,2])            #find edge length of recieving brach                            
	backBoneTree$edge.length[t] <-backBoneTree$edge.length[t]-treeLength              #cut receiving tip by length of clade tree		
	bound.tree<-bind.tree(backBoneTree,cladeTree,where=bindtip)
	return(bound.tree)
	}
	
	
#drop outgroups first	
	
fam.replace<-function(backBoneTree,cladeTree){

	out<-list()
	cd<-clade.depth(backBoneTree[[1]],cladeTree[[1]])
	for(i in 1:length(backBoneTree)){
		out[[i]]<-replace.clade(backBoneTree[[i]],cladeTree[[i]],cd)
		cat(i,"_")
		}
	class(out)<-"multiPhylo"	
	return(out)	
	}	
			
			
###Trying with bind.tree			

graft.replace<-function(backBoneTree,cladeTree){
	out<-list()
	c<-cladeTree[[1]]$tip.label      ###get names of clade tree
	ca<-llply(backBoneTree,getMRCA,c,.progress="text")     ###Find MRCA of clade in backbone
	f<-llply(c,grep,backBoneTree[[1]]$tip.label)	
	bbt<-unlist(f)
	o<-paste(backBoneTree[[1]]$tip.label[bbt],"_old")
		t<-list()    ### Rename tips to be removed in old tree
		for (i in 1:1000){ 
			t[[i]]<-backBoneTree[[i]]
			t[[i]]$tip.label[bbt]<-o	
					}
		class(t)<-"multiPhylo"
		bt<-list()
	for(j in 1:length(backBoneTree)){
			bt[[j]]<-bind.tree(backBoneTree[[j]],cladeTree[[j]],where=ca[[j]])
			cat(j,"\n")
			}
		class(bt)<-"multiPhylo"	
		
	dt<-lapply(bt,drop.tip,tip=o)		###drop old clade in backbone
	class(dt)<-"multiPhylo"
	return(dt)
	}
			### graft new clade at base of old clade

	make tree multiPhylo
	###return tree
	
	
	#### when this is done save all the trees as jpg plots to look at them

btest<-bind.tree(backBoneTree[[1]],cladeTree[[1]],where=ca[[1]])
	
t$tip.label[bbt]

	
backBoneTree[[1]]$tip.label[bbt]	
paste(backBoneTree[[1]]$tip.label[bbt],"old")
	
test<-backBoneTree[[900]]		
test$tip.label[bbt]<-paste(test$tip.label[bbt],"old")

	
backBoneTree<-JG
cladeTree<-PodJGF	
	
	
	
	###get MRCA of bound Podicipedidae clade (so we don't drop the outgroups twice

PodfJAF<-llply(PodJAF,drop.tip,c("Phoenicopterus_ruber","Grus_japonensis"),.progress="text") 

podTip<-PodJAF[[1]]$tip.label
podCut<-llply(JAF,getMRCA,podTip,.progress="text")


PhoenNode<-llply(JAF,getMRCA,Phoen,.progress="text")


getDescendants<-function(tree,node,curr=NULL){
  if(is.null(curr)) curr<-vector()
  daughters<-tree$edge[which(tree$edge[,1]==node),2]
  curr<-c(curr,daughters)
  w<-which(daughters>=length(tree$tip))
  if(length(w)>0) for(i in 1:length(w)) 
    curr<-getDescendants(tree,daughters[w[i]],curr)
  return(curr)
}



JAF.D<-list()
for(i in 1:length(JAF)){
JAF.D[[i]]<-drop.tip(JAF[[i]],getDescendants(JAF[[i]],podCut[[i]]))
}

class(JAF.D)<-"multiPhylo"


JAF.P<-list()
for (i in 1:length(JAF)){
	
	JAF.P[[i]]<-bind.tree(JAF.D[[i]],PodJAF[[i]],where=PhoenNode[[i]])
	}
class(JAF.P)<-"multiPhylo"
	






#remember to only use 1000 subclade tree and remove the outgroups!

### Jetz Gene

##Clades to replace
require(plyr)
require(phytools)
JG<-read.nexus("~/Dropbox/Flightless_project/Trees/All_Birds/Jetz_gene.nex")


	
# Anatide
	Ana_JG<-read.nexus("~/Dropbox/Flightless_project/Trees/Anatidae/JetzGeneFlightless/AnatidaeJetzGeneFlightless.nexus.t")
	AnaJG<-list()
	AnaJG<-Ana_JG[1:1000]
	AnaJGF<-llply(AnaJG,drop.tip,c("Dendrocygna_javanica","Anser_albifrons"),.progress="text")
	class(AnaJGF)<-"multiPhylo"


	JGF.1<-fam.replace(JG,AnaJGF)
	write.nexus(JGF.1,file="~/Dropbox/Flightless_project/Trees/JGF.1.nexus")
	
# Rallidae
JFG.1<-read.tree("/Users/ryanterrill/Dropbox/Flightless_project/Trees/All_Birds/All_Birds_JGF_Ana.tree")
	RGF<-read.nexus("/Users/ryanterrill/Dropbox/Flightless_project/Trees/Rallidae/Rallidae_gene_plus_flightless/Rallidae_JetzGeneFlightless.nexus.t")
	RalJGF<-list()
	RalJGF<-RGF[1:1000]
	RalJGF<-llply(RalJGF,drop.tip,c("Grus_canadensis","Podiceps_cristatus"),.progress="text")
	class(RalJGF)<-"multiPhylo"

	JGF.2<-fam.replace(JGF.1,RalJGF)
	write.nexus(JGF.2,file="~/Dropbox/Flightless_project/Trees/JGF.2.nexus")

	
### 8 June 16 at this point###

# Anserinae
	Ans_JGF<-read.nexus("~/Dropbox/Flightless_project/Trees/Anseridae/JetzGeneFlightless/Anserinae_JetzGeneFLightless.nexus.t")
	AnsJGF<-list()
	AnsJGF<-Ans_JGF[1:1000]
	AnsJGF<-llply(AnsJGF,drop.tip,c("Anas_platyrhynchos","Dendrocygna_javanica"),.progress="text")
	class(AnsJGF)<-"multiPhylo"

	JGF.3<-fam.replace(JGF.2,AnsJGF)
		write.nexus(JGF.3,file="~/Dropbox/Flightless_project/Trees/JGF.3.nexus")

	
### Here on 9 June


# Ciconiidae
	Cic_JGF<-read.nexus("/Users/ryanterrill/Dropbox/Flightless_project/Trees/Ciconiidae/Jetz_Gene_Flightless/Ciconiidae_JetzGeneFlightless.nexus.t")
	CicJGF<-list()
	CicJGF<-Cic_JGF[1:1000]
	CicJGF<-llply(CicJGF,drop.tip,c("Gavia_immer","Theristicus_caudatus"),.progress="text")
	class(CicJGF)<-"multiPhylo"

	JGF.4<-fam.replace(JGF.3,CicJGF)
		write.nexus(JGF.4,file="~/Dropbox/Flightless_project/Trees/JGF.4.nexus")


# Columbidae
	Col_JGF<-read.nexus("/Users/ryanterrill/Dropbox/Flightless_project/Trees/Columbidae/JetzGeneFlightless/Columbidae_JetzGeneFlightless.nexus.t")
	ColJGF<-list()
	ColJGF<-Col_JGF[1:1000]
	ColJGF<-llply(ColJGF,drop.tip,c("Pterocles_gutturalis","Colius_striatus"),.progress="text")
	class(ColJGF)<-"multiPhylo"

	JGF.5<-fam.replace(JGF.4,ColJGF)
		write.nexus(JGF.5,file="~/Dropbox/Flightless_project/Trees/JGF.5.nexus")


# Emberizidae
	Emb_JGF<-read.nexus("~/Dropbox/Flightless_project/Trees/Emberizidae/JetzGeneFlightless/Emberizidae_p20b_JetzGeneFlightless.nexus.t")
	EmbJGF<-list()
	EmbJGF<-Emb_JGF[1:1000]
	EmbJGF<-llply(EmbJGF,drop.tip,c("Parula_pitiayumi","Coereba_flaveola"),.progress="text")
	class(EmbJGF)<-"multiPhylo"

	JGF.6<-fam.replace(JGF.5,EmbJGF)
		write.nexus(JGF.6,file="~/Dropbox/Flightless_project/Trees/JGF.6.nexus")

# Falconidae
	Fal_JGF<-read.nexus("~/Dropbox/Flightless_project/Trees/Falconidae/JetzAllFlightless/Falconidae_JetzAllFlightless.nexus.t")
	FalJGF<-list()
	FalJGF<-Fal_JGF[1:1000]
	FalJGF<-llply(FalJGF,drop.tip,c("Tyto_alba","Parula_pitiayumi"),.progress="text")
	class(FalJGF)<-"multiPhylo"

	JGF.7<-fam.replace(JGF.6,FalJGF)
		write.nexus(JGF.7,file="~/Dropbox/Flightless_project/Trees/JGF.7.nexus")


# Gruidae
	Gru_JGF<-read.nexus("~/Dropbox/Flightless_project/Trees/Gruidae/JetzGeneFlightless/Gruidae_JetzGeneFlightless.nexus.t")
	GruJGF<-list()
	GruJGF<-Gru_JGF[1:1000]
	GruJGF<-llply(GruJGF,drop.tip,c("Psophia_viridis","Fulica_americana"),.progress="text")
	class(GruJGF)<-"multiPhylo"

	JGF.8<-fam.replace(JGF.7,GruJGF)
		write.nexus(JGF.8,file="~/Dropbox/Flightless_project/Trees/JGF.8.nexus")

	
# Megapodidae
	Meg_JGF<-read.nexus("~/Dropbox/Flightless_project/Trees/Megapodidae/JetzGeneFlightless/Megapodidae_JetzGeneFlightless.nexus.t")
	MegJGF<-list()
	MegJGF<-Meg_JGF[1:1000]
	MegJGF<-llply(MegJGF,drop.tip,c("Chauna_torquata","Gallus_gallus"),.progress="text")
	class(MegJGF)<-"multiPhylo"

	JGF.9<-fam.replace(JGF.8,MegJGF)
		write.nexus(JGF.9,file="~/Dropbox/Flightless_project/Trees/JGF.9.nexus")

	
# Podicipedidae
	Pod_JGF<-read.nexus("~/Dropbox/Flightless_project/Trees/Podicipedidae/JetzGeneFlightless/Podicipedidae_JetzGene.nexus.t")
	PodJGF<-list()
	PodJGF<-Pod_JGF[1:1000]
	PodJGF<-llply(PodJGF,drop.tip,c("Phoenicopterus_ruber","Grus_japonensis"),.progress="text")
	class(PodJGF)<-"multiPhylo"

	JGF.10<-fam.replace(JGF.9,PodJGF)
		write.nexus(JGF.10,file="~/Dropbox/Flightless_project/Trees/JGF.10.nexus")

####use bind.tree to bind this as sister to Phoenicopteridae

	
	JGF.9<-llply(JGF.9,drop.tip,c("Phoenicopterus_ruber","Grus_japonensis"),.progress="text")
	jtest<-JGF.9[1:2];ptest<-PodJGF[1:2]
	test<-fam.replace(jtest,ptest)
	

#backbone trees

backbone_JGF<-read.nexus("~/Dropbox/Flightless_project/Trees/Backbone/backbone.nexus.t")
backboneJGF<-list()
backboneJGF<-backbone_JGF[1:1000]
class(backboneJGF)<-"multiPhylo"


# Upupa Antaios
	#extract Upupidae clade	
	
	Phoeniculidae<-read.csv("~/Dropbox/Flightless_project/Trees/Phoeniculidae.csv",header=F)

	Upu_JGF_nodes<-findMRCA(backboneJGF[[1]],c(as.vector(Phoeniculidae$V1),"Upupa_epops","Upupa_antaios"))	
	
	Upu_JGF_nodes<-unlist(Upu_JGF_nodes)
	
	Upu_JGF<-llply(backboneJGF,extract.clade,284,.progress="text")
	
	class(Upu_JGF)<-"multiPhylo"
	
	JGF.11<-fam.replace(JGF.10,Upu_JGF)
	
# Rhynochetos orarius
	# extract Rhynochetidae
	
		Rhy_JGF_Nodes<-findMRCA(backboneJGF[[1]],c("Eurypyga_helias","Rhynochetos_jubatus","Rhynochetos_orarius"))	

		Rhy_JGF_nodes<-unlist(Rhy_JGF_Nodes)
		
		Rhy_JGF<-llply(backboneJGF,extract.clade,344,.progress="text")
		
		class(Rhy_JGF)<-"multiPhylo"
		
		JGF.12<-fam.replace(JGF.11,Rhy_JGF)
		
	
#Acanthisittidae
	#extract Acanthisittidae
	
	Aca_Nodes<-findMRCA(backboneJGF[[1]],c("Dendroscansor_decurvirostris","Traversia_lyalli","Pachyplichas_yaldwyni","Acanthisitta_chloris","Xenicus_gilviventris"))	

	Aca_JGF<-llply(backboneJGF,extract.clade,200,.progress="text")
	
	class(Aca_JGF)<-"multiPhylo"
	
	JGF.13<-fam.replace(JGF.12,Aca_JGF)

### Jetz All

##Clades to replace
# Rallidae
# Anatide
# Anserinae
# Ciconiidae
# Columbidae
# Emberizidae
# Falconidae
# Gruidae
# Megapodidae
# Phalacracoraciidae
# Podicipedidae
# Species to graft to backbone Tree - Acanthisittids, Rhynochetos orarius, Upupa Antaios,






### Jetz All

##Clades to replace
require(plyr)
require(phytools)
JG<-read.tree("~/Dropbox/Flightless_project/Trees/HackettStage1Full_1.tre")


# Rallidae
	Ral_JG<-read.nexus("~/Dropbox/Flightless_project/Trees/Rallidae/Rallidae_gene_plus_flightless/Rallidae_JetzGeneFlightless.nexus.t")
	RalJG<-list()
	RalJG<-Ral_JG[1:1000]
	RalJG<-llply(RalJG,drop.tip,c("Grus_canadensis","Podiceps_cristatus"),.progress="text")
	class(RalJG)<-"multiPhylo"

	JG.1<-fam.replace(JG,RalJG)
# Anatide
	Ana_JG<-read.nexus("~/Dropbox/Flightless_project/Trees/Anatidae/JetzGeneFlightless/AnatidaeJetzGeneFlightless.nexus.t")
	AnaJG<-list()
	AnaJG<-Ana_JG[1:1000]
	AnaJG<-llply(AnaJG,drop.tip,c("Anser_albifrons_anser","Dendrocygna_javanica"),.progress="text")
	class(AnaJG)<-"multiPhylo"

	JG.2<-fam.replace(JG.1,AnaJG)
	
# Anserinae
	Ans_JG<-read.nexus("~/Dropbox/Flightless_project/Trees/Anseridae/JetzGeneFlightless/Anserinae_JetzGeneFLightless.nexus.t")
	AnsJG<-list()
	AnsJG<-Ans_JG[1:1000]
	AnsJG<-llply(AnsJG,drop.tip,c("Anas_platyrhynchos","Dendrocygna_javanica"),.progress="text")
	class(AnsJG)<-"multiPhylo"

	JG.3<-fam.replace(JG.2,AnsJG)
	
# Ciconiidae
	Cic_JG<-read.nexus("~/Dropbox/Flightless_project/Trees/Ciconiidae/Jetz_Gene_Flightless/Ciconiidae_JetzGeneFlightless.nexus.t")
	CicJAF<-list()
	CicJAF<-Cic_JAF[1:1000]
	CicJAF<-llply(CicJAF,drop.tip,c("Gavia_immer","Theristicus_caudatus"),.progress="text")
	class(CicJAF)<-"multiPhylo"

	JAF.1<-fam.replace(JAF,CicJAF)
# Columbidae
	Col_JG<-read.nexus("~/Dropbox/Flightless_project/Trees/Columbidae/JetzGeneFlightless/Columbidae_JetzGeneFlightless.nexus.t")
	ColJG<-list()
	ColJG<-Col_JG[1:1000]
	ColJG<-llply(ColJG,drop.tip,c("Pterocles_gutturalis","Colius_striatus"),.progress="text")
	
	JG.5<-fam.replace(JG.4,ColJG)
	
# Emberizidae
	Emb_JG<-read.nexus("~/Dropbox/Flightless_project/Trees/Emberizidae/JetzGeneFlightless/Emberizidae_p20b_JetzGeneFlightless.nexus.t")
	EmbJG<-list()
	EmbJG<-Emb_JG[1:1000]
	EmbJG<-llply(EmbJG,drop.tip,c("Parula_pitiayumi","Coereba_flaveola"),.progress="text")
	JG.6<-fam.replace(JG.5,EmbJG)

# Falconidae
	Fal_JG<-read.nexus("~/Dropbox/Flightless_project/Trees/Falconidae/JetzGeneFlightless/Falconidae_JetzGeneFLightless.nexus.t")
	FalJG<-list()
	FalJG<-Fal_JG[1:1000]
	FalJG<-llply(FalJG,drop.tip,c("Tyto_alba","Parula_pitiayumi"),.progress="text")
	
	JG.7<-fam.replace(JG.6,FalJG)
	
# Gruidae
	Gru_JG<-read.nexus("~/Dropbox/Flightless_project/Trees/Gruidae/JetzGeneFlightless/Gruidae_JetzGeneFlightless.nexus.t")
	GruJG<-list()
	GruJG<-Gru_JG[1:1000]
	GruJG<-llply(GruJG,drop.tip,c("Psophia_viridis","Fulica_americana"),.progress="text")
	
	JG.8<-fam.replace(JG.7,FalJG)

# Megapodidae
	Meg_JG<-read.nexus("~/Dropbox/Flightless_project/Trees/Megapodidae/JetzGeneFlightless/Megapodidae_JetzGeneFlightless.nexus.t")
	MegJG<-list()
	MegJG<-Meg_JG[1:1000]
	MegJG<-llply(MegJG,drop.tip,c("Chauna_torquata","Gallus_gallus"),.progress="text")
	
	JG.9<-fam.replace(JG.8,MegJG)

# Podicipedidae
JG.9<-read.nexus("~/Dropbox/Flightless_project/Trees/All_Birds/JGF_final.nexus")
	Pod_JG<-read.nexus("~/Dropbox/Flightless_project/Trees/Podicipedidae/JetzGeneFlightless/Podicipedidae_JetzGene.nexus.t")
	PodJG<-list()
	PodJG<-Pod_JG[1:1000]
	PodJG<-llply(PodJG,drop.tip,c("Phoenicopterus_ruber","Grus_japonensis"),.progress="text")
	
	JG.10<-fam.replace(JG.9,PodJG)



	
# Upupa Antaios
	#extract Upupidae clade	
	
	Phoeniculidae<-read.csv("~/Dropbox/Flightless_project/Trees/Phoeniculidae.csv",header=F)

	Upu_Ja_nodes<-findMRCA(backboneJA[[1]],c(as.vector(Phoeniculidae$V1),"Upupa_epops","Upupa_antaios"))	
	
	Upu_Ja_nodes<-unlist(Upu_Ja_nodes)
	
	Upu_JA<-llply(backboneJA,extract.clade,284,.progress="text")
	
	class(Upu_JA)<-"multiPhylo"
	
	JG.11<-fam.replace(JA.10,Upu_JA)
	
# Rhynochetos orarius
	# extract Rhynochetidae
	
		Rhy_JA_Nodes<-findMRCA(backboneJA[[1]],c("Eurypyga_helias","Rhynochetos_jubatus","Rhynochetos_orarius"))	

		Rhy_JA_nodes<-unlist(Rhy_JA_Nodes)
		
		Rhy_JA<-llply(backboneJA,extract.clade,344,.progress="text")
		
		class(Rhy_JA)<-"multiPhylo"
		
		JG.12<-fam.replace(JG.11,Rhy_JA)
		
	
#Acanthisittidae
	#extract Acanthisittidae
	
	Aca_Nodes<-findMRCA(backboneJA[[1]],c("Dendroscansor_decurvirostris","Traversia_lyalli","Pachyplichas_yaldwyni","Acanthisitta_chloris","Xenicus_gilviventris"))	

	Aca_JA<-llply(backboneJA,extract.clade,200,.progress="text")
	
	class(Aca_JA)<-"multiPhylo"
	
	Jg.13<-fam.replace(JG.12,Aca_JA)
 
 
 #Done! Now run that shit through BayesTraits!
