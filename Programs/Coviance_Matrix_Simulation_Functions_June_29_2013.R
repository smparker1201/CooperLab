#!/usr/bin/env Rscript

#Here are all the functions necessary for the simulation, with comments.
#In summary, the simulation code makes a bunch of simulated datasets with parameters that you can set; it chooses some of the genes in the simulated data to be in the 
#gene set; and it tests the significance of the gene set in the simulated data using both the GSEA and Siegel-Tukey (SiT) algorithms

#You will need to make sure you have downloaded this special R library (MASS); it contains the multivariate normal distribution function
library(MASS) 

#Gene_Set_Analysis_Simulated_Data is the main function.
#sigma is the correlation matrix for you to input. Make sure it's dimensions are equal to the number of genes you specify in the input. If you don't specify a correlation matrix (ie let sigma=""), my original simulation method will run.
#sigma_method specifies how to use the given correlation matrix. sigma_method="best_updown_separate" takes the maximum correlation in the matrix and gives that correlation to every gene in the gene set (actually, the genes going up have that correlation with the other genes going up and the genes going down have that correlation with the other genes going down). sigma_method="best_updown_together" gives the best correlation to all genes in the gene set regardless of direction. That's all the methods I coded in - any input besides that will just leave the correlations as input, essentially choosing random rows and columns in the correlation matrix to be the gene set correlations. This is just the first way I could think of to do this - feel free to manipulate this part in particular, or let me know your idea and I'll code it for you.
#nsim is the number of simulated datasets to make and then evaluate. I'd recommend nsim=10000 for a good simulation if you can let it run that long (took my computer about 8 hours or so). nsim=1000 may be ok though.
#ngenes is the number of genes in each of the simulated data sets
#nsamples is the number of experimental and control samples in each of the simulated data sets (so nsamples=6 means 6 control and 6 experimental)
#nperm is the number of permuted data sets tried in the p value calculation of both the GSEA and Siegel-Tukey algorithms.
#control_variance is the variance in the control expression values - log fold change between experimental and control is just mean(experimental)-mean(control). This code fixes mean(control) to be zero for each gene, but with a variance that is chosen here.
#top_X, of_Y and direction indicate the ranks of the genes that are in the gene set. For example: if top_X=10, of_Y=50 and direction="top", the code makes sure that 10 genes are in the gene set and all of them are in the top 50 by rank (ranked by fold change).
#direction can take values of "top", "bottom" and "both". For "both", you take a double set; for example 10 of the top 50 as well as 10 of the bottom 50.
#output_file_name is simply what you want to name your file with the simulation results.
#write.data is a TRUE or FALSE; if true, it saves all the simulated gene expression values - which is (2*nsamples)*(ngenes)*(nsim) entries. The file is often huge, so the default is false.
Gene_Set_Analysis_Simulated_Data<-function(sigma=readMatrix(),sigma_method="best_updown_together",covariance_method="good",nsim=1000,ngenes=1000,nsamples=6,nperm=1000,control_variance=1,top_X=10,of_Y=50,direction="top",output_file_name="sub1togtopcontvar1",write.data=FALSE){

  #Rename the output file, if no name is given.
  print(sigma_method)
  if(output_file_name==""){
    output_file_name=paste("Results_nsim_",as.character(nsim),"_ngenes_",as.character(ngenes),"_nsamples_",as.character(nsamples),"_control_variance_",as.character(control_variance),"_top_",as.character(top_X),"_of_",as.character(of_Y),"_direction_",direction,".txt",sep="")
  }
  

#Start loop over simulations, define variables
  results<-matrix(0,nsim,10) #for Siegel-Tukey and GSEA (GSEA listed first), keep score, normalized score, p value, FDR, FDR non-normalized
  nullNormGSEA<-matrix(0,nsim,nperm)
  nullNormSiT<-matrix(0,nsim,nperm)
  if(write.data){
	allDataUsed<-matrix(0,ngenes,2*nsamples*nsim)
  }

  for(i in 1:nsim){
    print(paste("Simulation Number ",i,sep=""))
    #choose fold changes for all genes (log values, normal dist, mean 0, var 1 - pretty standard assumptions)
    fc<-rnorm(ngenes,0,1)
    #knowing those fold changes, select which genes will be in your gene set
    fcOrderTop<-order(fc,decreasing=T) #the indices of the biggest positive fc down to the indices of the most negative fc
    indexTop<-sample(fcOrderTop[1:of_Y],top_X,replace=FALSE) #draws some of the top fold changes - for example, it could be 10 of the top 50 fold change genes.
    fcOrderBottom<-order(fc,decreasing=F)
    indexBottom<-sample(fcOrderBottom[1:of_Y],top_X,replace=FALSE)
    if(direction=="both"){
      geneSetIndex<-c(indexTop,indexBottom)
    }
    if(direction=="top"){
      geneSetIndex<-indexTop
    }
    if(direction=="bottom"){
      geneSetIndex<-indexBottom
    }


    #define the covariance matrix - if no input for it, it runs the way I initially defined things.
	#NOTE: Originally this had if(sigma="") and it gave an error so I changed it -Sarah
    if(sigma==""){
      sigma<-matrix(0,ngenes,ngenes)
      diag(sigma)<-1 #every gene correlates with itself
      for(j in 1:(length(geneSetIndex)-1)){
        for(k in (j+1):length(geneSetIndex)){
          sigma[geneSetIndex[j],geneSetIndex[k]]<-1 #note that in my original simulations, I didn't separate the correlation of the up and down genes!
          sigma[geneSetIndex[k],geneSetIndex[j]]<-1
        }
      }
    }
	
	#NOTE: I also did a '==' edit here -Sarah
    if(sigma!="" & sigma_method=="best_updown_separated"){
      bestValue<-max(sigma)
      if(length(indexTop)>0){
        for(j in 1:(length(indexTop)-1)){
          for(k in j:length(indexTop)){
            sigma[indexTop[j],indexTop[k]]<-bestValue #this code sets the genes to have "bestValue" correlation with themselves as well. I think that's a good idea but maybe you should think about it too to make sure.
            sigma[indexTop[k],indexTop[j]]<-bestValue
          }
        }
      }
      if(length(indexBottom)>0){
        for(j in 1:(length(indexBottom)-1)){
          for(k in j:length(indexBottom)){
            sigma[indexBottom[j],indexBottom[k]]<-bestValue #this code sets the genes to have "bestValue" correlation with themselves as well. I think that's a good idea but maybe you should think about it too to make sure.
            sigma[indexBottom[k],indexBottom[j]]<-bestValue
          }
        }
      }
    }
    if(sigma!="" & sigma_method=="best_updown_together"){
      bestValue<-max(sigma)
      for(j in 1:(length(geneSetIndex)-1)){
        for(k in j:length(geneSetIndex)){
          sigma[geneSetIndex[j],geneSetIndex[k]]<-bestValue #this code sets the genes to have "bestValue" correlation with themselves as well. I think that's a good idea but maybe you should think about it too to make sure.
          sigma[geneSetIndex[k],geneSetIndex[j]]<-bestValue
        }
      }
    }
	


    #use mvrnorm to create the "experimental samples" - representing log2 gene expression values
    experimentalSamples<-t(mvrnorm(nsamples,fc,sigma,tol=2)) #the output of mvrnorm comes out transposed (rows/columns switched); the command t() transposes it back
    #subtract exactly the right number so that all the means are equal the fold changes drawn above; otherwise the ranks computed in "fcOrderTop", etc. would be somewhat off
    diff<-fc-apply(experimentalSamples,1,mean)
    experimentalSamples<-apply(experimentalSamples,2,"+",diff)
    #create control sample values
    controlSamples<-rnorm(nsamples*ngenes,0,control_variance)
    dim(controlSamples)<-c(ngenes,nsamples)
    diff<-0-apply(controlSamples,1,mean)
    controlSamples<-apply(controlSamples,2,"+",diff)



    #run GSEA, Siegel-Tukey on the data. Notice we're running both methods on the same sets of simulated data (to make it a fair comparison)
    STRanks<-Siegel_Tukey_Rank_Producer(ngenes)
    fcAll<-apply(experimentalSamples,1,mean)-apply(controlSamples,1,mean) #this will equal the variable "fc" computed above
    fcSimAll<-Shuffled_Data_Maker(cbind(experimentalSamples,controlSamples),nsamples,nsamples,nperm)
    scoreDataGSEA<-GSEA_Score_Data(fcAll,geneSetIndex)
    scoreDataSiT<-Siegel_Tukey_Score_Data(fcAll,geneSetIndex,STRanks)
    permGSEA<-GSEA_Permutations_P_Value(fcSimAll,geneSetIndex,scoreDataGSEA)
    permSiT<-Siegel_Tukey_Permutations_P_Value(fcSimAll,geneSetIndex,scoreDataSiT,STRanks)
    nullNormGSEA[i,]<-permGSEA$scores
    nullNormSiT[i,]<-permSiT$scores
    results[i,c(1:3,6:8)]<-c(scoreDataGSEA,permGSEA$NES,permGSEA$p_value,scoreDataSiT,permSiT$NES,permSiT$p_value)
    if(write.data){
      allDataUsed[,(2*nsamples*(i-1)+1):(2*nsamples*i)]<-cbind(controlSamples,experimentalSamples)
    }
}
  
  
  #now that we have computed scores and p values for nsim number of simulations, use these values to correct for multiple hypothesis testing and to get normalized scores
  #this is for the GSEA results:
  dim(nullNormGSEA)<-NULL
  NESUp<-nullNormGSEA[which(nullNormGSEA>=0)] #notice for the GSEA data, we follow their methodology and divide the positive score gene sets from the negative
  NESDown<-nullNormGSEA[which(nullNormGSEA<0)]
  totUp<-length(NESUp)
  totDown<-length(NESDown)
  indexUp<-which(results[,2]>=0)
  indexDown<-which(results[,2]<0)
  for(geneSetCount in indexUp){
    results[geneSetCount,4]<-sum(NESUp>=results[geneSetCount,2])/totUp #this is a partial result on the way to computing the false discovery rate
  }
  for(geneSetCount in indexDown){
    results[geneSetCount,4]<-sum(NESDown<=results[geneSetCount,2])/totDown #this is a partial result on the way to computing the false discovery rate
  }
  
  indexUpGeneSets<-which(results[,2]>=0)
  indexDownGeneSets<-which(results[,2]<0)
  NESUpGeneSets<-results[indexUpGeneSets,2]
  NESDownGeneSets<-results[indexDownGeneSets,2]
  for(geneSetCount in indexUpGeneSets){
	results[geneSetCount,5]<-results[geneSetCount,4]/(sum(NESUpGeneSets>=results[geneSetCount,2])/length(indexUpGeneSets)) #this is the false discovery rate (p value corrected for multiple hypothesis testing), using the GSEA method.
  }
  for(geneSetCount in indexDownGeneSets){
	results[geneSetCount,5]<-results[geneSetCount,4]/(sum(NESDownGeneSets<=results[geneSetCount,2])/length(indexDownGeneSets)) #this is the false discovery rate
  }

  #Now the same process for the Siegel-Tukey (SiT) results, except that all the Siegel-Tukey scores will be positive:
  dim(nullNormSiT)<-NULL
  NESUp<-nullNormSiT[which(nullNormSiT>=0)]
  totUp<-length(NESUp)
  indexUp<-which(results[,7]>=0)
  for(geneSetCount in indexUp){
    results[geneSetCount,9]<-sum(NESUp<=results[geneSetCount,7])/totUp #changed direction
  }
  indexUpGeneSets<-which(results[,7]>=0)
  NESUpGeneSets<-results[indexUpGeneSets,7]
  for(geneSetCount in indexUpGeneSets){
	results[geneSetCount,10]<-results[geneSetCount,9]/(sum(NESUpGeneSets<=results[geneSetCount,7])/length(indexUpGeneSets)) #changed direction
  }

  #Write results to file
  if(write.data){
    write.table(allDataUsed,paste("all_data_used_control_first_",output_file_name,sep=""),sep="\t",quote=F,row.names=F,col.names=F)
  }

  colnames(results)<-c("ES_GSEA","NES_GSEA","pval_GSEA","FDR_Non_Norm_GSEA","FDR_GSEA","ES_SiT","NES_SiT","pval_SiT","FDR_Non_Norm_SiT","FDR_SiT")
  results[which(results[,5]>1),5]<-1 #sometimes the false discovery rate is greater than one; this has the same meaning as being equal to one.
  results[which(results[,10]>1),10]<-1
  results<-results[,c(1:3,5:8,10)] #I don't need the non-normalized false discovery rates anymore.
  write.table(results,output_file_name,sep="\t",quote=F,row.names=F,col.names=T)

}


#These are the other functions that are called by the main function. At one point, I used these functions on data and got exactly the same result as running the GSEA program.

Siegel_Tukey_Rank_Producer<-function(numGenes){
  #The Siegel-Tukey test assigns ranks in a different and important way (ie, top is 1, bottom is 2, etc.). For speed, this is only called once
  a<-1:numGenes
  b<-c(1,rep(c(-1,-1,1,1),ceiling(length(a)/4)))[1:numGenes]
  bTop<-a[which(b==1)]
  bBottom<-rev(a[which(b== -1)])
  ranks<-c(bTop,bBottom)
  return(ranks)
}

Shuffled_Data_Maker<-function(data,n1,n2,nsim){
  #This makes all of the permuted data in the GSEA method. data is the data, n1 and n2 are the number of experimental and control samples respectively.
  sampleBinary<-matrix(0,nsim,n1+n2)
  sampleBinary[,1:n1]<-1
  sampleBinary<-apply(sampleBinary,1,sample)
  
  mean1<-data %*% sampleBinary #matrix multiplication
  mean2<-data %*% (1-sampleBinary)
  
  fcSim<-mean1-mean2

  return(fcSim)
}

GSEA_Score_Data<-function(measureStat,geneSetIndices){
  #Computes the GSEA score. measureStat is the statistic used to rank (here, fold change) and geneSetIndices indicate which genes are in the gene set.
  ranksAll<-order(order(measureStat,decreasing=T),decreasing=F)
  ranksGeneSet<-ranksAll[geneSetIndices]
  pHit<-rep(0,length(ranksAll))
  pMiss<-rep(0,length(ranksAll))
  pHit[ranksGeneSet]<-abs(measureStat[geneSetIndices])
  pHit<-pHit/sum(pHit)
  pMiss[setdiff(1:length(ranksAll),ranksGeneSet)]<-1/(length(ranksAll)-length(ranksGeneSet))

  pHit<-cumsum(pHit)
  pMiss<-cumsum(pMiss)

  pDiff<-pHit-pMiss

  ES<-ifelse(max(pDiff)> -1*min(pDiff),max(pDiff),min(pDiff))
  return(ES)
}


GSEA_Permutations_P_Value<-function(fcSim,geneSetIndices,ES){
  #Computes the p value in the GSEA method.
  scores<-signif(apply(fcSim,2,GSEA_Score_Data,geneSetIndices),5)
  ESUp<-scores[which(scores>=0)]  
  ESDown<-scores[which(scores<0)]
  
  if(ES>=0){
    p_value<-sum(ESUp>=ES)/length(ESUp)
    NES<-ES/mean(ESUp)
  }
  if(ES<0){
    p_value<-sum(ESDown<=ES)/length(ESDown)
    NES<-ES/mean(abs(ESDown))
  }
  NESUp<-ESUp/mean(ESUp)
  NESDown<-ESDown/mean(abs(ESDown))
  return(list("p_value"=p_value,"scores"=c(NESUp,NESDown),"NES"=NES))
}


Siegel_Tukey_Score_Data<-function(measureStat,geneSetIndices,STRanks){
  #Computes the Siegel-Tukey score. STRanks is the ranks from the function Siegel_Tukey_Rank_Producer.
  ranksAll<-order(order(measureStat,decreasing=T),decreasing=F)
  ranksGeneSet<-ranksAll[geneSetIndices]
  STRanksGeneSet<-STRanks[ranksGeneSet]
  n<-length(geneSetIndices)
  STScore<-sum(STRanksGeneSet)-n*(n+1)/2
  return(STScore)
}


Siegel_Tukey_Permutations_P_Value<-function(fcSim,geneSetIndices,ES,STRanks){
  #Computes the Siegel-Tukey method p value. Differs from GSEA in that all scores are positive.
  scores<-apply(fcSim,2,Siegel_Tukey_Score_Data,geneSetIndices,STRanks)
  ESUp<-scores[which(scores>=0)]  
  p_value<-sum(ESUp<=ES)/length(ESUp)
  NES<-ES/mean(ESUp)
  NESUp<-ESUp/mean(ESUp)
  return(list("p_value"=p_value,"scores"=NESUp,"NES"=NES))
}

readMatrix<-function(){

	print("Read In Matrix")
	myFile=read.table("/home/smparker/Research/Cooper/Matrices/ToProcess/mat1.csv",sep=",")
	matrix=as.matrix(myFile,nrow=1000,ncol=1000)

	origEig = do.call(c,eigen(matrix,only.values=TRUE))
	
	cholError = any(origEig<=0)
	print(cholError)


	iter <- 0
		value=0
	while (cholError) {

		iter <- iter + 1
		cat("iteration ", iter, "\n")

		for(i in 1:1000){
			matrix[i,i]=matrix[i,i]+.00001
			}

		# replace -ve eigen values with small +ve number
		value=value+.00001
		newEig= eigen(matrix,only.values=TRUE)
		newEig=do.call(c,newEig)
		newEig=array(newEig)
		# try chol again
		
		cholError = any(newEig<=0)
	
		}
		print('VALUE ADDED:')
		print(value)
		return(matrix)


	}

