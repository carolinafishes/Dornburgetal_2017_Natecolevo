library(lattice)
##############
##############part 1 comparing DEC and DEC J for each
##############


###+J versus Regular DEC

ilikebio<-function(directory,string,number)
{
likelihoods<-matrix(nrow=length(number), ncol=1)
for (i in 1:number)
{
paste(directory, string,i,sep="")-> resfn	

    load(resfn)
    likelihoods [i]<- get_LnL_from_BioGeoBEARS_results_object(res)

}
return(likelihoods)
	
}

###phase 2 use likelihoods from 1 folder and compare AIC scores, weights
stringA<-"nototreeDEC_"
stringB<-"nototreeDECJ_"
stringC<-"REALnototreeDECA_"
stringD<-"nototreeDECA_"

model_AICS<-function(directory,stringA,stringB,number,p1,p2)
{

	ilikebio(directory,stringA,number)->likemodel_1
	ilikebio(directory,stringB,number)->likemodel_2
	##change
	-2*likemodel_1+2*p1->model1_AIC
	-2*likemodel_2+2*p2->model2_AIC
	cbind(model1_AIC, model2_AIC)->aics
	
	abs(model1_AIC-model2_AIC)->DeltaAIC
	
	mean(likemodel_1)->l1
	mean(likemodel_2)->l2
	sd(likemodel_1)->sdl1
	sd(likemodel_2)->sdl2
	
	mean(model1_AIC)->a1
	mean(model2_AIC)->a2
	sd(model1_AIC)->sda1
	sd(model2_AIC)->sda2
	
	mean(DeltaAIC)->da
	sd(DeltaAIC)->sdda
	
	c(l1,sdl1,l2,sdl2,a1,sda1,a2,sda2,da,sdda)->output
	c("likelhood model 1","SD model 1","likelhood model 2","SD model 2","AIC model 1","SD AIC model 1","AIC model 2","SD AIC model 2","delta AIC", "SD delta aic")->legend
	rbind(legend, output)->results
	hist(DeltaAIC, xlim=c(0,max(DeltaAIC)))
	arrows(4,0,x1=4,y1=100, col="blue")
	return(results)


}

##summary
setwd("~/Dropbox/New_Nototheniod/No_model/")
directory<-("~/Dropbox/New_Nototheniod/No_model/")


number<-(500)
#Nomodel J vs DEC




##############part2
#####compare all models
##############
#rbind(directoryA,directoryA1,directoryA2,directoryA3,directoryA4)->directories
directoryA

rbind(stringB,stringC)->stringsB



##size corrected weights, read in
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
p1<-2
p2<-3
p3<-4
#note this function assumes parameter count of 2,3,3,4 have the same number of paramters 
mini_AICS<-function(directory,stringA,stringB, stringC,stringD, number,p1,p2,p3)
{

#likeliA<-matrix(nrow=number, ncol=columns)
#likeliB<-matrix(nrow=number, ncol=columns)

	ilikebio(directory,stringA,number)->likeliA
		ilikebio(directory,stringB,number)->likeliB
			ilikebio(directory,stringC,number)->likeliC
				ilikebio(directory,stringD,number)->likeliD
	

	##change
	-2* likeliA+2*p1->model1_AIC
	-2*likeliB+2*p2->model2_AIC
	-2*likeliC+2*p2->model3_AIC
	-2*likeliD+2*p3->model4_AIC
	cbind(model1_AIC, model2_AIC, model3_AIC,model4_AIC)->all_AIC_scores



AICweightdistro<-matrix(nrow=number, ncol=4)
for (i in 1:number)
{
AICw(all_AIC_scores[i,])->AICweightdistro[i,]
}
#namings<-c(stringsA,stringsB)
AICweightdistro->weights
return(weights)
}




#note this function assumes all results in stringA or stringB have the same number of paramters 
mega_AICS<-function(directories,stringsA,stringsB,number,p1,p2)
{
length(directories)->columns
likeliA<-matrix(nrow=number, ncol=columns)
likeliB<-matrix(nrow=number, ncol=columns)
for (i in 1:columns)
{

	ilikebio(directories[i],stringsA[i],number)->likeliA[,i]
	}
for (i in 1:columns)
{

	ilikebio(directories[i],stringsB[i],number)->likeliB[,i]
	}	
	##change
	-2* likeliA+2*p1->model1_AIC
	-2*likeliB+2*p2->model2_AIC
	cbind(model1_AIC, model2_AIC)->all_AIC_scores



AICweightdistro<-matrix(nrow=number, ncol=columns*2)
for (i in 1:number)
{
AICw(all_AIC_scores[i,])->AICweightdistro[i,]
}
#namings<-c(stringsA,stringsB)
AICweightdistro->weights
return(weights)
}

#plot it using lattice!
rbind(as.numeric(weights[,1]) ,as.numeric(weights[,2]) ,as.numeric(weights[,3]) ,as.numeric(weights[,4]))->stack
weight_histo <- data.frame(gp = factor(rep(paste('model', 1:4, sep = ''), each = 
 1000)), x = stack) 

 histogram( ~ as.numeric(weights[,1]) +as.numeric(weights[,2]) +as.numeric(weights[,3]) +as.numeric(weights[,4]),layout=c(1,4), data = weight_histo) 

























