### Four clades 
library(geiger)
tree<- read.tree("~/Dropbox/Final_nototheniod_Biogeography/Transitions/Notothen_tree.tre"); ## all tree
read.tree("~/Dropbox/Final_nototheniod_Biogeography/making_files/newick.trees")->trees
length(trees)->loop






Patago<-c("Patagonotothen_cornucolaC","Patagonotothen_elegansA", "Patagonotothen_guntheriD", "Patagonotothen_longipesB", "Patagonotothen_ramsayiP",  "Patagonotothen_simaA","Patagonotothen_squamicepsA","Patagonotothen_tessellataB","Patagonotothen_wiltoniA") 
name<-"Harpagifer_bispinisB"
###Get Harpagifer
get.age<-function(tree,name)
{
n<-length(tree$tip.label)
ee<-setNames(tree$edge.length[sapply(1:n,function(x,y) which(y==x),y=tree$edge[,2])],tree$tip.label)

ee[which(names(ee)==name)]->age
return(age)
}


harp<-matrix(nrow=loop)
for(i in 1:loop)
{
	get.age(trees[[i]],name)->harp[i]
	
}

harp

hist(harp[,1])

##Next two!
name<-"Champsocephalus_esoxB"
name<-"Paranotothenia_magellanicaA"


name<-"Champsocephalus_esoxB"
###Get Champsocephalus
get.age<-function(tree,name)
{
  n<-length(tree$tip.label)
  ee<-setNames(tree$edge.length[sapply(1:n,function(x,y) which(y==x),y=tree$edge[,2])],tree$tip.label)
  
  ee[which(names(ee)==name)]->age
  return(age)
}

###Get Champsocephalus

champ<-matrix(nrow=loop)
for(i in 1:loop)
{
  get.age(trees[[i]],name)->champ[i]
  
}

#harp

hist(champ[,1])

###Get dissostichus

name<-"Dissostichus_eleginoidesC"
get.age<-function(tree,name)
{
  n<-length(tree$tip.label)
  ee<-setNames(tree$edge.length[sapply(1:n,function(x,y) which(y==x),y=tree$edge[,2])],tree$tip.label)
  
  ee[which(names(ee)==name)]->age
  return(age)
}

dis<-matrix(nrow=loop)
for(i in 1:loop)
{
  get.age(trees[[i]],name)->dis[i]
  
}



hist(dis[,1], breaks=24)


name<-"Paranotothenia_magellanicaA"
###Get Paranotothenia
get.age<-function(tree,name)
{
  n<-length(tree$tip.label)
  ee<-setNames(tree$edge.length[sapply(1:n,function(x,y) which(y==x),y=tree$edge[,2])],tree$tip.label)
  
  ee[which(names(ee)==name)]->age
  return(age)
}

###Get Paranotothenia

para<-matrix(nrow=loop)
for(i in 1:loop)
{
  get.age(trees[[i]],name)->para[i]
  
}

para

hist(para[,1])

###Patago attempt

name<-"Paranotothenia_magellanicaA"
###Get Paranotothenia
get.age<-function(tree,name)
{
  n<-length(tree$tip.label)
  ee<-setNames(tree$edge.length[sapply(1:n,function(x,y) which(y==x),y=tree$edge[,2])],tree$tip.label)
  
  ee[which(names(ee)==name)]->age
  return(age)
}

###Get Paranotothenia

para<-matrix(nrow=loop)
for(i in 1:loop)
{
  get.age(trees[[i]],name)->para[i]
  
}

para

hist(para[,1])


Patago = c("Patagonotothen_cornucolaC","Patagonotothen_elegansA", "Patagonotothen_guntheriD", "Patagonotothen_longipesB", "Patagonotothen_ramsayiP",  "Patagonotothen_simaA","Patagonotothen_squamicepsA","Patagonotothen_tessellataB","Patagonotothen_wiltoniA")
mrca(tree)["Patagonotothen_longipesB", "Patagonotothen_tessellataB"]

#Patago <- node.leaves(tree, 102)

Patago <- matrix(nrow=loop)
for(i in 1:loop)
{
  trees[[i]]->tree2
  mrca(tree2)["Patagonotothen_longipesB", "Patagonotothen_tessellataB"]-> temp
  which(tree2$tip.label=="Patagonotothen_longipesB") -> test
Patago[i]<- dist.nodes(tree2)[temp,test]

  }
rbind(Patago[,1],para[,1],champ[,1],harp[,1])->all
