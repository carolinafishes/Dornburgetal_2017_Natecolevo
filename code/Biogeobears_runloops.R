setwd("/Users/admadornburg/Dropbox/new_nototheniod/Final_summary_tree_runs")

##loop DEC
read.tree("~/Dropbox/notothen_biogeography/making_files/newick.trees")->trees
length(trees)->loop



geogfn = "~/Dropbox/notothen_biogeography/notothenioid_temp.data"
multmatrix = "~/Dropbox/notothen_biogeography/making_files/noto.rate.matrix.mat.txt"
for (i in 1:1000)
{
	paste("~/Dropbox/notothen_biogeography/No_model/nototree_",i,sep="")->trfn
	paste("nototreeDEC_",i,sep="")->res_file
	Bio_dec(trfn, geo=geogfn,res_file,7)
}

###now DECJ
DECFILE<-"~/Dropbox/notothen_biogeography/No_model/nototreeDEC_1"

for (i in 1:loop)
{
	paste("~/Dropbox/notothen_biogeography/No_model/nototree_",i,sep="")->trfn
	paste("nototreeDECJ_",i,sep="")->res_file
	Bio_decj(trfn, geogfn,res_file,7,DECFILE)
} 



#single


setwd("/Users/admadornburg/Dropbox/new_nototheniod/Final_summary_tree_runs")
"/Users/admadornburg/Dropbox/new_nototheniod/Final_summary_tree_runs/Notothen_tree.tre"->trfn
geogfn = "/Users/admadornburg/Dropbox/new_nototheniod/Final_summary_tree_runs/notothenioid_temp.data"


Bio_dec(trfn, geo=geogfn,res_file="nototreeDEC",7)
DECFILE<-"/Users/admadornburg/Dropbox/new_nototheniod/Final_summary_tree_runs/nototreeDEC"
Bio_decj(trfn, geo=geogfn,res_file="nototreeDECJ",7,DECFILE)
Bio_deca(trfn, geo=geogfn,res_file="nototreeDECA",7,DECFILE)
Bio_decaj(trfn, geo=geogfn,res_file="nototreeDECAJ",7,DECFILE)
