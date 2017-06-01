library(phytools)

read.tree("~/Dropbox/new_nototheniod/Final_summary_tree_runs/cyronoto.phy")->treefornodenumbers
treefornodenumbers->tree1
#plot(treefornodenumbers)
#nodelabels()
setwd("~/Dropbox/new_nototheniod/No_model")
resfn = "nototreeDECA_138"
library(BioGeoBEARS)

 source("http://phylo.wdfiles.com/local--files/biogeobears/cladoRcpp.R") # (needed now that traits model added; source FIRST!)
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_add_fossils_randomly_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_basics_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_calc_transition_matrices_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_classes_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_detection_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_DNA_cladogenesis_sim_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_extract_Qmat_COOmat_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_generics_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_models_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_on_multiple_trees_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_plots_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_readwrite_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_simulate_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_SSEsim_makePlots_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_SSEsim_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_stochastic_mapping_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_stratified_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_univ_model_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/calc_uppass_probs_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/calc_loglike_sp_v01.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/get_stratified_subbranch_top_downpass_likelihoods_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/runBSM_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/stochastic_map_given_inputs.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/summarize_BSM_tables_v1.R")
source("http://phylo.wdfiles.com/local--files/biogeobears/BioGeoBEARS_traits_v1.R") # added traits model
calc_loglike_sp = compiler::cmpfun(calc_loglike_sp_prebyte)    # crucial to fix bug in uppass calculations
calc_independent_likelihoods_on_each_branch = compiler::cmpfun(calc_independent_likelihoods_on_each_branch_prebyte)
    # slight speedup hopefully
 load(resfn)
 areas<-c("PSA","SGSS","KI","HA","SA", "NZ","AU") 
  rcpp_areas_list_to_states_list(areas=areas, maxareas=7, include_null_range=TRUE)->states

mal<-vector(length=128)
for (i in 1:128)
{
	3%in%states[[i]]->mal[i]
	
}
which(mal=="TRUE")->allthrees


  ##note node = c("taxonA", "taxonB")
 get.split<-function(node,length, string_tree, string_file)
 {
 	HA<-matrix(nrow=length)
 	 	SA<-matrix(nrow=length)
 	 		timer<-matrix(nrow=length)
 for(i in 1:length)
 {
 	paste(string_tree,i,sep="")->trfn
 	read.tree(trfn)->tree
 	findMRCA(tree,node, type="node")->rownum
 paste(string_file,i,sep="")->bioresult	
 resfn = bioresult	
  res$ML_marginal_prob_each_state_at_branch_top_AT_node->probtable
  
  probtable[rownum,]->row_of_interest
  row_of_interest[5]->HA[i]
  1-row_of_interest[5]->SA[i]
  branching.times(tree)->ll
  which(names(ll)== rownum)->part_of_vector
(ll[[part_of_vector]])->timer[i]
   #print(timer)
 }
  cbind("SA", "HA", "time")->header
  cbind(SA,HA,timer)->output
  colnames(output)<-header
  return(output)
   }  
   
    get.split.general<-function(node,length, string_tree, string_file,column)
 {
 	HA<-matrix(nrow=length)
 	 	SA<-matrix(nrow=length)
 	 		timer<-matrix(nrow=length)
 for(i in 1:length)
 {
 	paste(string_tree,i,sep="")->trfn
 	read.tree(trfn)->tree
 	findMRCA(tree,node, type="node")->rownum
 paste(string_file,i,sep="")->bioresult	
 resfn = bioresult	
  res$ML_marginal_prob_each_state_at_branch_top_AT_node->probtable
  
  probtable[rownum,]->row_of_interest
  sum(row_of_interest[column])->HA[i]
  1-sum(row_of_interest[column])->SA[i]
  branching.times(tree)->ll
  which(names(ll)== rownum)->part_of_vector
(ll[[part_of_vector]])->timer[i]
   #print(timer)
 }
  cbind("Everything_Else", "Row_Interest", "time")->header
  cbind(SA,HA,timer)->output
  colnames(output)<-header
  return(output)
   }  
   


get.likelihoods<-function(string_fileA, length)
{
	
	fitA<-matrix(nrow=length)
				
for (i in 1:length)
{
	 paste(string_fileA,i,sep="")->bioresult	
 	resfn = bioresult	
get_LnL_from_BioGeoBEARS_results_object(res)-> fitA[i]

	
}
return(fitA)
}
weight<-function(AIC1, AIC2)
{	
	exp((-1*abs(AIC1-AIC2))/2)->weight
	return(weight)	
}
AICw<-function(AICvector)
{
	length(AICvector)->n
	min(AICvector)->min
	weights<-vector(length=n)
	for(i in 1:n)
	{
		AICvector[i]->currentmodel
		weight(min,currentmodel)->weights[i]
	}
	sum(weights)->denom
	AICweights<-vector(length=n)
		for(i in 1:n)
	{
		weights[i]->currentweight	
		currentweight/denom->AICweights[i]	
	}
	#cbind(AICvector, AICweights)-> output
	return(AICweights)
	
}

  model.average<-function(node,length, string_tree, string_fileA, string_fileB, string_fileC, string_fileD, p1,p2,p3,p4)
  {
  	

  	
  	get.split(node, length, string_tree, string_fileA)->cladeA
  	get.split(node, length, string_tree, string_fileB)->cladeB
  	get.split(node, length, string_tree, string_fileC)->cladeC
  	get.split(node, length, string_tree, string_fileD)->cladeD
  
 
  	
get.likelihoods(string_fileA, length)->likeA
get.likelihoods(string_fileB, length)->likeB
get.likelihoods(string_fileC, length)->likeC
get.likelihoods(string_fileD, length)->likeD
  	
-2* likeA +2*p1->model1_AIC
-2* likeB +2*p2->model2_AIC
-2* likeC +2*p3->model3_AIC
-2* likeD +2*p4->model4_AIC

cbind(model1_AIC, model2_AIC, model3_AIC, model4_AIC)->AICS

#return(AICS)

weights<-matrix(nrow=length, ncol=4)
for (i in 1:length)
{
AICS[i,]->vector
AICw(vector)->weights[i,]
	
}
cbind((cladeA[,1]*weights[,1]),(cladeA[,2]*weights[,1]))->cladaA.a
	cbind((cladeB[,1]*weights[,2]),(cladeB[,2]*weights[,2]))->cladaB.a
		cbind((cladeC[,1]*weights[,3]),(cladeC[,2]*weights[,3]))->cladaC.a
			cbind((cladeD[,1]*weights[,4]),(cladeD[,2]*weights[,4]))->cladaD.a

cladaA.a+ cladaB.a+ cladaC.a+ cladaD.a->modelav
cbind(modelav, cladeA[,3])-> modelav2
  cbind("SA", "HA", "Time")->header
  colnames(modelav2)<-header

return(modelav2)

  	
  } 
  
  #  Usage()
#model.average("Notothenia_coriicepsB", "Paranotothenia_magellanicaA", length=5, "nototree_", "nototreeDEC_", "nototreeDECJ_", "REALnototreeDECA_", "nototreeDECA_", 2,3,3,4)->test
#model.average("Notothenia_coriicepsB", "Paranotothenia_magellanicaA", length=5, "nototree_", "nototreeDEC_", "nototreeDECJ_", "REALnototreeDECA_", "nototreeDECA_", 2,3,3,4)->test

#median(test[,1])

 model.average.general<-function(node,length, string_tree, string_fileA, string_fileB, string_fileC, string_fileD, p1,p2,p3,p4, column)
  {
  	

  	
  	get.split.general(node, length, string_tree, string_fileA, column)->cladeA
  	  	get.split.general(node, length, string_tree, string_fileB, column)->cladeB
  	  	get.split.general(node, length, string_tree, string_fileC, column)->cladeC
  	  	get.split.general(node, length, string_tree, string_fileD, column)->cladeD
  
 
  	
get.likelihoods(string_fileA, length)->likeA
get.likelihoods(string_fileB, length)->likeB
get.likelihoods(string_fileC, length)->likeC
get.likelihoods(string_fileD, length)->likeD
  	
-2* likeA +2*p1->model1_AIC
-2* likeB +2*p2->model2_AIC
-2* likeC +2*p3->model3_AIC
-2* likeD +2*p4->model4_AIC

cbind(model1_AIC, model2_AIC, model3_AIC, model4_AIC)->AICS

#return(AICS)

weights<-matrix(nrow=length, ncol=4)
for (i in 1:length)
{
AICS[i,]->vector
AICw(vector)->weights[i,]
	
}
cbind((cladeA[,1]*weights[,1]),(cladeA[,2]*weights[,1]))->cladaA.a
	cbind((cladeB[,1]*weights[,2]),(cladeB[,2]*weights[,2]))->cladaB.a
		cbind((cladeC[,1]*weights[,3]),(cladeC[,2]*weights[,3]))->cladaC.a
			cbind((cladeD[,1]*weights[,4]),(cladeD[,2]*weights[,4]))->cladaD.a

cladaA.a+ cladaB.a+ cladaC.a+ cladaD.a->modelav
cbind(modelav, cladeA[,3])-> modelav2
  cbind("SA", "HA", "Time")->header
  colnames(modelav2)<-header

return(modelav2)

  	
  } 


################################################################
################################
####Giant Heatmap Time...
########################################
############################################################################

###core function to get single value per million years

record.parse<-function(node.output)
{
	node.output[which(node.output[,3]<30 & node.output[,3]>29),]->sub29
		length(dim(sub29))->d29
		if (d29==2){
			mean(sub29[,2])->sub29
			}else sub29[2]->sub29
				
	node.output[which(node.output[,3]<29 & node.output[,3]>28),]->sub28
		length(dim(sub28))->d28
		if (d28==2){
			mean(sub28[,2])->sub28
			}else sub28[2]->sub28
	node.output[which(node.output[,3]<28 & node.output[,3]>27),]->sub27
		length(dim(sub27))->d27
		if (d27==2){
			mean(sub27[,2])->sub27
			}else sub27[2]->sub27
	node.output[which(node.output[,3]<27 & node.output[,3]>26),]->sub26
		length(dim(sub26))->d26
		if (d26==2){
			mean(sub26[,2])->sub26
			}else sub26[2]->sub26
	node.output[which(node.output[,3]<26 & node.output[,3]>25),]->sub25
			length(dim(sub25))->d25
		if (d25==2){
			mean(sub25[,2])->sub25
			}else sub25[2]->sub25
	node.output[which(node.output[,3]<24 & node.output[,3]>23),]->sub24
			length(dim(sub24))->d24
		if (d24==2){
			mean(sub24[,2])->sub24
			}else sub24[2]->sub24
	node.output[which(node.output[,3]<23 & node.output[,3]>22),]->sub23
			length(dim(sub23))->d23
		if (d23==2){
			mean(sub23[,2])->sub23
			}else sub23[2]->sub23
	node.output[which(node.output[,3]<22 & node.output[,3]>21),]->sub22
			length(dim(sub22))->d22
		if (d22==2){
			mean(sub22[,2])->sub22
			}else sub22[2]->sub22
	node.output[which(node.output[,3]<21 & node.output[,3]>20),]->sub21
			length(dim(sub21))->d21
		if (d21==2){
			mean(sub21[,2])->sub21
			}else sub21[2]->sub21
	node.output[which(node.output[,3]<20 & node.output[,3]>19),]->sub20
			length(dim(sub20))->d20
		if (d20==2){
			mean(sub20[,2])->sub20
			}else sub20[2]->sub20
	node.output[which(node.output[,3]<19 & node.output[,3]>18),]->sub19
			length(dim(sub19))->d19
		if (d19==2){
			mean(sub19[,2])->sub19
			}else sub19[2]->sub19
	node.output[which(node.output[,3]<18 & node.output[,3]>17),]->sub18
			length(dim(sub18))->d18
		if (d18==2){
			mean(sub18[,2])->sub18
			}else sub18[2]->sub18
	node.output[which(node.output[,3]<17 & node.output[,3]>16),]->sub17
			length(dim(sub17))->d17
		if (d17==2){
			mean(sub17[,2])->sub17
			}else sub17[2]->sub17
	node.output[which(node.output[,3]<16 & node.output[,3]>15),]->sub16		
			length(dim(sub16))->d16
		if (d16==2){
			mean(sub16[,2])->sub16
			}else sub16[2]->sub16				
	node.output[which(node.output[,3]<15 & node.output[,3]>14),]->sub15
			length(dim(sub15))->d15
		if (d15==2){
			mean(sub15[,2])->sub15
			}else sub15[2]->sub15	
	node.output[which(node.output[,3]<14 & node.output[,3]>13),]->sub14
			length(dim(sub14))->d14
		if (d14==2){
			mean(sub14[,2])->sub14
			}else sub14[2]->sub14
	node.output[which(node.output[,3]<13 & node.output[,3]>12),]->sub13
			length(dim(sub13))->d13
			if (d13==2){
			mean(sub13[,2])->sub13
			}else sub13[2]->sub13
	node.output[which(node.output[,3]<12 & node.output[,3]>11),]->sub12
			length(dim(sub12))->d12
		if (d12==2){
			mean(sub12[,2])->sub12
			}else sub12[2]->sub12
	node.output[which(node.output[,3]<11 & node.output[,3]>10),]->sub11
			length(dim(sub11))->d11
		if (d11==2){
			mean(sub11[,2])->sub11
			}else sub11[2]->sub11
	node.output[which(node.output[,3]<10 & node.output[,3]>9),]->sub10
			length(dim(sub10))->d10
		if (d10==2){
			mean(sub10[,2])->sub10
			}else sub10[2]->sub10
	node.output[which(node.output[,3]<9 & node.output[,3]>8),]->sub9
			length(dim(sub9))->d9
		if (d9==2){
			mean(sub9[,2])->sub9
			}else sub9[2]->sub9
	node.output[which(node.output[,3]<8 & node.output[,3]>7),]->sub8
			length(dim(sub8))->d8
		if (d8==2){
			mean(sub8[,2])->sub8
			}else sub8[2]->sub8
	node.output[which(node.output[,3]<7 & node.output[,3]>6),]->sub7
			length(dim(sub7))->d7
		if (d7==2){
			mean(sub7[,2])->sub7
			}else sub7[2]->sub7
	node.output[which(node.output[,3]<6 & node.output[,3]>5),]->sub6
			length(dim(sub6))->d6
		if (d6==2){
			mean(sub6[,2])->sub6
			}else sub6[2]->sub6
	node.output[which(node.output[,3]<5 & node.output[,3]>4),]->sub5
			length(dim(sub5))->d5
		if (d5==2){
			mean(sub5[,2])->sub5
			}else sub5[2]->sub5
	node.output[which(node.output[,3]<4 & node.output[,3]>3),]->sub4
			length(dim(sub4))->d4
		if (d4==2){
			mean(sub4[,2])->sub4
			}else sub4[2]->sub4
	node.output[which(node.output[,3]<3 & node.output[,3]>2),]->sub3
			length(dim(sub3))->d3
		if (d3==2){
			mean(sub3[,2])->sub3
			}else sub3[2]->sub3
	node.output[which(node.output[,3]<2 & node.output[,3]>1),]->sub2
			length(dim(sub2))->d2
		if (d2==2){
			mean(sub2[,2])->sub2
			}else sub2[2]->sub2
	node.output[which(node.output[,3]<1 & node.output[,3]>0),]->sub1
			length(dim(sub1))->d1
		if (d1==1){
			mean(sub1[,2])->sub1
			}else sub1[2]->sub1

c(sub29,sub28,sub27,sub26,sub25,sub24,sub23,sub22,sub21,sub20,sub19,sub18,sub17,sub16,sub15,sub14,sub13,sub12,sub11,sub10,sub9,sub8,sub7,sub6,sub5,sub4,sub3,sub2,sub1)->rowdata
rowdata[is.na(rowdata)]<-0


return(rowdata)
										
}

create.map.data<-function(guidetree,rootnode, maxnode, length, string_tree, string_fileA, string_fileB, string_fileC, string_fileD, p1,p2,p3,p4, filename)
{
	maxnode-rootnode+1->total.nodes
	rows.of.doom<-matrix(ncol=29)
	for (i in rootnode:maxnode)
	{
	getDescendants(guidetree, i)->desc
	as.vector(na.omit(guidetree$tip.label[desc]))->node
	model.average(node, length, string_tree, string_fileA, string_fileB, string_fileC, string_fileD, p1,p2,p3,p4)->test	
	record.parse(test)->rows.of.doomtemp
	rbind(rows.of.doom, rows.of.doomtemp)->rows.of.doom
	}
	as.data.frame(rows.of.doom)-> rows.of.doom
	print(dim(rows.of.doom))
	colnames(rows.of.doom)<-c("29","28","27","26","25","24","23","22","21","20","19","18","17","16","15","14","13","12","11","10","9","8","7","6","5","4","3","2","1")
	rownames(rows.of.doom)<-c("ignore",rootnode:maxnode)
	write.csv(rows.of.doom, file=filename)
	return(rows.of.doom)
	
}


create.map.data.general<-function(guidetree,rootnode, maxnode, length, string_tree, string_fileA, string_fileB, string_fileC, string_fileD, p1,p2,p3,p4, filename, column)
{
	maxnode-rootnode+1->total.nodes
	rows.of.doom<-matrix(ncol=29)
	for (i in rootnode:maxnode)
	{
	getDescendants(guidetree, i)->desc
	as.vector(na.omit(guidetree$tip.label[desc]))->node
	model.average.general(node, length, string_tree, string_fileA, string_fileB, string_fileC, string_fileD, p1,p2,p3,p4,column)->test	
	record.parse(test)->rows.of.doomtemp
	rbind(rows.of.doom, rows.of.doomtemp)->rows.of.doom
	}
	as.data.frame(rows.of.doom)-> rows.of.doom
	print(dim(rows.of.doom))
	colnames(rows.of.doom)<-c("29","28","27","26","25","24","23","22","21","20","19","18","17","16","15","14","13","12","11","10","9","8","7","6","5","4","3","2","1")
	rownames(rows.of.doom)<-c("ignore",rootnode:maxnode)
	write.csv(rows.of.doom, file=filename)
	return(rows.of.doom)
	
}

#create.map.data(tree1, 85,167,length=5, "nototree_", "nototreeDEC_", "nototreeDECJ_", "REALnototreeDECA_", "nototreeDECA_", 2,3,3,4, filename="test_5reps")->test

#create.map.data.general(tree1, 85,167,length=500, "nototree_", "nototreeDEC_", "nototreeDECJ_", "REALnototreeDECA_", "nototreeDECA_", 2,3,3,4, filename="ANYwithHA_500reps", column=mal)->test3
# following code limits the lowest and highest color to 5%, and 95% of your range, respectively

colors = c(seq(0,300,length=100),seq(301,600,length=100),seq(601,1000,length=100))
my_palette <- colorRampPalette(c("gray", "blue", "orange"))(n = 299)
#heatmap.2(HA, scale="none", dendrogram="none", Colv=FALSE, col=my_palette )



###### All you need to plot once this is run Multiply each by 1000
#### Just HA
read.csv("~/Dropbox/new_nototheniod/No_model/final_HA_endemic")->HA
HA[,2:30]->HA
as.matrix(HA)->HA
 colnames(HA)<-rev(c(1:29))
 rownames(HA)<-c(85:167)
 heatmap.2(HA, scale="none", dendrogram="none", Colv=FALSE,Rowv=FALSE, symm=F,symkey=F,symbreaks=T,col=my_palette, breaks=colors)
#### any HA
read.csv("~/Dropbox/new_nototheniod/No_model/ANYwithHA_500reps")->anyHA
anyHA[2:84,2:30]->aHA
as.matrix(aHA)->aHA2
colnames(aHA2)<-rev(c(1:29)) 
rownames(aHA2)<-c(85:167)
heatmap.2(aHA2, scale="none", dendrogram="none", Colv=FALSE,Rowv=FALSE, col=my_palette )

###PSA
read.csv("~/Dropbox/new_nototheniod/No_model/Correct_PSA+HA_500reps")->PSAHA
PSAHA[2:84,2:30]->aPSAHA
 as.matrix(aPSAHA)->ap2
 colnames(ap2)<-rev(c(1:29))
 rownames(ap2)<-c(85:167)
 heatmap.2(ap2, scale="none", dendrogram="none", Colv=FALSE,Rowv=FALSE, col=my_palette )








#############Get each node's state
create.node.data.general<-function(guidetree,rootnode, maxnode, length, string_tree, string_fileA, string_fileB, string_fileC, string_fileD, p1,p2,p3,p4, filename, column)
{
	maxnode-rootnode+1->total.nodes
	rows.of.doom<-0
	for (i in rootnode:maxnode)
	{
	getDescendants(guidetree, i)->desc
	as.vector(na.omit(guidetree$tip.label[desc]))->node
	model.average.general(node, length, string_tree, string_fileA, string_fileB, string_fileC, string_fileD, p1,p2,p3,p4,column)->test
	median(test[,2])-> test2	
	c(rows.of.doom, test2)->rows.of.doom
	}
	as.data.frame(rows.of.doom)-> rows.of.doom
	print(dim(rows.of.doom))
	rownames(rows.of.doom)<-c("ignore",rootnode:maxnode)
	write.csv(rows.of.doom, file=filename)
	return(rows.of.doom)
	
}

create.node.data.general(tree1, 85,167,length=500, "nototree_", "nototreeDEC_", "nototreeDECJ_", "REALnototreeDECA_", "nototreeDECA_", 2,3,3,4, filename="Acrosstime_Just_HA", column=5)->test3

create.node.data.general(tree1, 85,167,length=500, "nototree_", "nototreeDEC_", "nototreeDECJ_", "REALnototreeDECA_", "nototreeDECA_", 2,3,3,4, filename="Acrosstime_Any_subwith_HA", column=mal)->test4