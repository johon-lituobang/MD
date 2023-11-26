if (!require("foreach")) install.packages("foreach")
library(foreach)
if (!require("doParallel")) install.packages("doParallel")
library(doParallel)
#require randtoolbox for random number generations
if (!require("randtoolbox")) install.packages("randtoolbox")
library(randtoolbox)
if (!require("Rcpp")) install.packages("Rcpp")
library(Rcpp)
if (!require("Rfast")) install.packages("Rfast")
library(Rfast)
if (!require("matrixStats")) install.packages("matrixStats")
library(matrixStats)






raw1 <- read.csv(("GSE147507_RawReadCounts_Ferret.csv"))
#colnames(raw1)[1] <- 'NAME'
# --- mean value
meta1 <- read.csv(("metadata.csv"))
#colnames(meta1)[1] <- 'NAME'


orderpairGSE147507<-orderpair(rawdata = raw1,metadata = meta1,paired=FALSE,ci = FALSE)

orderpairGSE147507A<-orderpairGSE147507$result
write.csv(orderpairGSE147507A,paste("orderpairGSE147507_day.csv", sep = ","), row.names = FALSE)
raw1 <- read.csv(("AEC_SARSCoV2_17AAG_readcounts.csv"))
#colnames(raw1)[1] <- 'NAME'
# --- mean value
meta1 <- read.csv(("metadata_AEC_SARSCoV2_17AAG.csv"))
#colnames(meta1)[1] <- 'NAME'

orderpairAEC_SARSCoV2<-orderpair(rawdata = raw1,metadata = meta1,paired=FALSE,ci = FALSE)

orderpairAEC_SARSCoV2<-orderpairAEC_SARSCoV2$result
write.csv(orderpairAEC_SARSCoV2,paste("orderpairAEC_SARSCoV2_day.csv", sep = ","), row.names = FALSE)

raw1 <- read.csv(("AEC_SARSCoV2_17AAG_readcounts.csv"))
#colnames(raw1)[1] <- 'NAME'
# --- mean value
meta1 <- read.csv(("metadata_AEC_SARSCoV2_17AAGmock.csv"))
#colnames(meta1)[1] <- 'NAME'

orderpairAEC_SARSCoV2<-orderpair(rawdata = raw1,metadata = meta1,paired=FALSE,ci = FALSE)

orderpairAEC_SARSCoV2<-orderpairAEC_SARSCoV2$result
write.csv(orderpairAEC_SARSCoV2,paste("orderpairAEC_SARSCoV2_day_mock.csv", sep = ","), row.names = FALSE)

orderpair<- function(rawdata,metadata,phylogenetictree=NULL,paired=FALSE,ci = FALSE,alpha=0.05,nboot=100,robust = 1){
  metadata$Order2 <- with(metadata,paste0(Subgroup,Order))
  metadataorder <- split(metadata,metadata$Subgroup)
  pairs0 <- data.frame()
  results0 <- data.frame()
  for (i in 1:length(metadataorder)) {
    cnt1 <- length(unique(metadataorder[[i]]$Order2))
    metadata1=metadataorder[[i]]
    ordernumber=as.numeric(length(unique(metadataorder[[i]]$Order)))
    orderlist=sort(unique(metadataorder[[i]]$Order),decreasing = FALSE)
    for (j in 1:(ordernumber-1)) {
      dissimilarity4<-dissimilarity(rawdata = raw1,metadata = metadata1,Order1 = toString(orderlist[j]), Order2 = toString(orderlist[j+1]),phylogenetictree=phylogenetictree,paired=paired,ci = ci,alpha=alpha,nboot=nboot,robust = robust)
      pair1 <- paste0(unique(metadataorder[[i]]$Subgroup),orderlist[j])
      pair2 <- paste0(unique(metadataorder[[i]]$Subgroup),orderlist[j+1])
      #pairs <- data.frame(cbind(pair1,pair2,meandivergens$estimate,sparsitydivergens$estimate,variancedivergens$estimaste,meandivergens$ci,sparsitydivergens$ci,variancedivergens$ci))
      pairs<-c(pair1=pair1,
               pair2=pair2,
               dissimilarity=(dissimilarity4$results))
      pairs0 <- rbind(pairs0,pairs)
      resul1<-c(pair1=pair1,
                pair2=pair2,
                rawdata=dissimilarity4$rawdata)
      resultlist<-list(resul1)
      names(resultlist) <- c(paste("A", i,j, sep = ""))
      results0 <- append(results0,resultlist)
    }
  }
  All<-list(result=pairs0,rawdata=results0)
  return(All)#,'orderpair_batchhuman.csv',row.names = FALSE)
}


moments<-function (x){
  n<-length(x)
  mean1=mean(x)
  var1<-sum((x-mean1)^2)/n
  tm1<-sum((x-mean1)^3)/n
  fm1<-sum((x-mean1)^4)/n
  listall<-list(mean=mean1,var=var1,tm=tm1,fm=fm1)
  (listall)
}
unbiasedmoments<-function (x){
  n<-length(x)
  moments1<-moments(x)
  m1<-moments1$mean
  var1<-moments1$var
  tm1<-moments1$tm
  fm1<-moments1$fm
  var2<-var1*(n/(n-1))
  tm2<-(tm1)*(n^2/((n-1)*(n-2)))
  ufm1<--3*var1^2*(2*n-3)*n/((n-1)*(n-2)*(n-3))+(n^2-2*n+3)*fm1*n/((n-1)*(n-2)*(n-3))
  listall<-c(mean=m1,var=var2,tm=tm2,fm=ufm1)
  (listall)
}
standardizedmoments<-function (x,releaseall=FALSE){
  unbiasedmoments1<-unbiasedmoments(x)
  all1<-c(mean=unbiasedmoments1[1],var=sqrt(unbiasedmoments1[2]),skew=unbiasedmoments1[3]/((unbiasedmoments1[2])^(3/2)),kurt=unbiasedmoments1[4]/((unbiasedmoments1[2])^(2)))
  if(releaseall){
    c(unbiasedmoments1,all1)
  }else{
    (all1)
  }
}

findmoments<- function(selected1=NULL,selected2=NULL,paired=FALSE){
  if(paired){
    selected3<-selected2-selected1
    meandi <- data.frame(apply(selected3,1,FUN = mean))
    meancombined3 <-data.frame(apply(meandi,1,abs))
    
    sd1 <- data.frame(apply(selected1,1,FUN = sd))
    sd2 <- data.frame(apply(selected2,1,FUN = sd))
    
    sdall <- cbind(sd1,sd2)
    sdcombined1 <-data.frame(apply(sdall,1,max))
    sdcombined2 <-data.frame(apply(sdall,1,min))
    
    sdcombined3 <-sdcombined1-sdcombined2
    sddi<- sd2-sd1
    
    select <- cbind(selected1,selected2)
    sumall<- data.frame(apply(select,1,FUN = sum))
    
    selected1zero<-rowSums(selected1==0)/(dim(selected1)[2])
    selected2zero<-rowSums(selected2==0)/(dim(selected2)[2])
  } else{
    mean1 <- data.frame(apply(selected1,1,FUN = mean))
    mean2 <- data.frame(apply(selected2,1,FUN = mean))
    
    meanall <- cbind(mean1,mean2)
    meancombined1 <-data.frame(apply(meanall,1,max))
    meancombined2 <-data.frame(apply(meanall,1,min))

    meancombined3 <-meancombined1-meancombined2
    meandi<- mean2-mean1
    
    sd1 <- data.frame(apply(selected1,1,FUN = sd))
    sd2 <- data.frame(apply(selected2,1,FUN = sd))
    
    sdall <- cbind(sd1,sd2)
    sdcombined1 <-data.frame(apply(sdall,1,max))
    sdcombined2 <-data.frame(apply(sdall,1,min))
    
    sdcombined3 <-sdcombined1-sdcombined2
    sddi<- sd2-sd1
    
    select <- cbind(selected1,selected2)
    sumall<- data.frame(apply(select,1,FUN = sum))
    
    selected1zero<-rowSums(selected1==0)/(dim(selected1)[2])
    selected2zero<-rowSums(selected2==0)/(dim(selected2)[2])
    
  }
  results<-list(meanmean=mean(unlist(meanall)),sdsd=mean(unlist(sdall)),sumall=sumall,meancombined3=meancombined3,meandi=meandi,sdcombined3=sdcombined3,sddi=sddi,selected1zero=selected1zero,selected2zero=selected2zero)
  return(results)
}
#treetoweight(phylogenetictree)
#matrixtoweight(correlationmatrix)
#phyloadjust<-treetoweight(phylotree)



dissimilarity<- function(rawdata,metadata,phylogenetictree=FALSE,Order1,Order2,paired=NULL,ci=NULL,alpha=NULL,nboot=NULL,robust=NULL){
  # rawdata<-raw1
  # metadata<-meta1
  # full=TRUE
  # paired=TRUE
  # Group1="faeces-amox-1"
  # Group2="faeces-amox-2"
  # NAME1=NULL
  # Treatment1=NULL
  # Subgroup1=NULL
  # Order1=NULL
  # NAME2=NULL
  # Treatment2=NULL
  # Subgroup2=NULL
  # Order2=NULL
  # phylo=TRUE
  # robust = 1
  rawdata[is.na(rawdata)] = 0
  if (!require("foreach")) install.packages("foreach")
  library(foreach)
  if (!require("doParallel")) install.packages("doParallel")
  library(doParallel)
  #require randtoolbox for random number generations
  if (!require("randtoolbox")) install.packages("randtoolbox")
  library(randtoolbox)
  if (!require("Rcpp")) install.packages("Rcpp")
  library(Rcpp)
  if (!require("Rfast")) install.packages("Rfast")
  library(Rfast)
  if (!require("matrixStats")) install.packages("matrixStats")
  library(matrixStats)
  indx1 <- metadata$NAME[metadata$Order==Order1]
  indx2 <- metadata$NAME[metadata$Order==Order2]

  selected1 <- rawdata[,indx1]
  selected2 <- rawdata[,indx2]
  selected1 <- data.frame(selected1)
  selected2 <- data.frame(selected2)
  results<-findmoments(selected1=selected1,selected2=selected2,paired=paired)
  
  
  meanmean=results$meanmean
  sdsd=results$sdsd
  sumall=results$sumall
  meancombined3=results$meancombined3
  meandi=results$meandi
  sdcombined3=results$sdcombined3
  sddi=results$sddi
  selected1zero=results$selected1zero
  selected2zero=results$selected2zero
  
  mean_dissimilarity <- wilcox.test(x=unlist(meancombined3),conf.int = TRUE)$estimate

  mean_dissimilarity_direct <- wilcox.test(x=unlist(meandi),conf.int = TRUE)$estimate

  standardized_mean_dissimilarity<-mean_dissimilarity/meanmean
  standardized_mean_dissimilarity_direct<-mean_dissimilarity_direct/meanmean

  sd_dissimilarity <- wilcox.test(x=unlist(sdcombined3),conf.int = TRUE)$estimate

  sd_dissimilarity_direct <- wilcox.test(x=unlist(sddi),conf.int = TRUE)$estimate

  standardized_sd_dissimilarity<-sd_dissimilarity/sdsd
  standardized_sd_dissimilarity_direct<-sd_dissimilarity_direct/sdsd

  sparsity_dissimilarity<-sum(meancombined3*abs(selected1zero-selected2zero))#*2/(sparsity_all[2]+sparsity_all[1])

  standardized_sparsity_dissimilarity<-sum(meancombined3*abs(selected1zero-selected2zero))*2/(sum(unlist(sumall)))

  sparsity_dissimilarity_direct<-sum(meandi*abs(selected1zero-selected2zero))#*2/(sparsity_all[2]+sparsity_all[1])

  standardized_sparsity_dissimilarity_direct<-sum(meandi*abs(selected1zero-selected2zero))*2/(sum(unlist(sumall)))
  
  R=50
  
  indices_list <- replicate(R, list(indices1 = sample(1:(dim(selected1)[1]), replace = TRUE)), simplify = FALSE)
  
  numCores <- detectCores()-4 # Detect the number of available cores
  cl <- makeCluster(numCores) # Create a cluster with the number of cores
  registerDoParallel(cl) # Register the parallel backend
  
  bootresults1 <- foreach(indices = indices_list, .packages = c("stats")) %dopar% {
    selected1_resampled <- selected1[indices$indices1,]
    selected2_resampled <- selected2[indices$indices1,]
    
    findmoments<- function(selected1=NULL,selected2=NULL,paired=FALSE){
      if(paired){
        selected3<-selected2-selected1
        meandi <- data.frame(apply(selected3,1,FUN = mean))
        meancombined3 <-data.frame(apply(meandi,1,abs))
        
        sd1 <- data.frame(apply(selected1,1,FUN = sd))
        sd2 <- data.frame(apply(selected2,1,FUN = sd))
        
        sdall <- cbind(sd1,sd2)
        sdcombined1 <-data.frame(apply(sdall,1,max))
        sdcombined2 <-data.frame(apply(sdall,1,min))
        
        sdcombined3 <-sdcombined1-sdcombined2
        sddi<- sd2-sd1
        
        select <- cbind(selected1,selected2)
        sumall<- data.frame(apply(select,1,FUN = sum))
        
        selected1zero<-rowSums(selected1==0)/(dim(selected1)[2])
        selected2zero<-rowSums(selected2==0)/(dim(selected2)[2])
      } else{
        mean1 <- data.frame(apply(selected1,1,FUN = mean))
        mean2 <- data.frame(apply(selected2,1,FUN = mean))
        
        meanall <- cbind(mean1,mean2)
        meancombined1 <-data.frame(apply(meanall,1,max))
        meancombined2 <-data.frame(apply(meanall,1,min))
        
        meancombined3 <-meancombined1-meancombined2
        meandi<- mean2-mean1
        
        sd1 <- data.frame(apply(selected1,1,FUN = sd))
        sd2 <- data.frame(apply(selected2,1,FUN = sd))
        
        sdall <- cbind(sd1,sd2)
        sdcombined1 <-data.frame(apply(sdall,1,max))
        sdcombined2 <-data.frame(apply(sdall,1,min))
        
        sdcombined3 <-sdcombined1-sdcombined2
        sddi<- sd2-sd1
        
        select <- cbind(selected1,selected2)
        sumall<- data.frame(apply(select,1,FUN = sum))
        
        selected1zero<-rowSums(selected1==0)/(dim(selected1)[2])
        selected2zero<-rowSums(selected2==0)/(dim(selected2)[2])
        
      }
      results<-list(meanmean=mean(unlist(meanall)),sdsd=mean(unlist(sdall)),sumall=sumall,meancombined3=meancombined3,meandi=meandi,sdcombined3=sdcombined3,sddi=sddi,selected1zero=selected1zero,selected2zero=selected2zero)
      return(results)
    }
    
    results<-findmoments(selected1=selected1_resampled,selected2=selected2_resampled,paired=paired)
    
    
    meanmean=results$meanmean
    sdsd=results$sdsd
    sumall=results$sumall
    meancombined3=results$meancombined3
    meandi=results$meandi
    sdcombined3=results$sdcombined3
    sddi=results$sddi
    selected1zero=results$selected1zero
    selected2zero=results$selected2zero
    
    mean_dissimilarity <- wilcox.test(x=unlist(meancombined3),conf.int = TRUE)$estimate
    
    mean_dissimilarity_direct <- wilcox.test(x=unlist(meandi),conf.int = TRUE)$estimate
    
    standardized_mean_dissimilarity<-mean_dissimilarity/meanmean
    standardized_mean_dissimilarity_direct<-mean_dissimilarity_direct/meanmean
    
    sd_dissimilarity <- wilcox.test(x=unlist(sdcombined3),conf.int = TRUE)$estimate
    
    sd_dissimilarity_direct <- wilcox.test(x=unlist(sddi),conf.int = TRUE)$estimate
    
    standardized_sd_dissimilarity<-sd_dissimilarity/sdsd
    standardized_sd_dissimilarity_direct<-sd_dissimilarity_direct/sdsd
    
    sparsity_dissimilarity<-sum(meancombined3*abs(selected1zero-selected2zero))#*2/(sparsity_all[2]+sparsity_all[1])
    
    standardized_sparsity_dissimilarity<-sum(meancombined3*abs(selected1zero-selected2zero))*2/(sum(unlist(sumall)))
    
    sparsity_dissimilarity_direct<-sum(meandi*abs(selected1zero-selected2zero))#*2/(sparsity_all[2]+sparsity_all[1])
    
    standardized_sparsity_dissimilarity_direct<-sum(meandi*abs(selected1zero-selected2zero))*2/(sum(unlist(sumall)))
    
    
    results2<-c(mean_dissimilarity=mean_dissimilarity,
                   mean_dissimilarity_direct=mean_dissimilarity_direct,standardized_mean_dissimilarity=standardized_mean_dissimilarity,
                   standardized_mean_dissimilarity_direct=standardized_mean_dissimilarity_direct,
                   sd_dissimilarity=sd_dissimilarity,
                   sd_dissimilarity_direct=sd_dissimilarity_direct,standardized_sd_dissimilarity=standardized_sd_dissimilarity,
                   standardized_sd_dissimilarity_direct=standardized_sd_dissimilarity_direct,
                   
                   sparsity_dissimilarity=sparsity_dissimilarity,
                   standardized_sparsity_dissimilarity=standardized_sparsity_dissimilarity,sparsity_dissimilarity_direct=sparsity_dissimilarity_direct,
                   standardized_sparsity_dissimilarity_direct=standardized_sparsity_dissimilarity_direct)
    results2
  }
  
  
  stopCluster(cl)
  registerDoSEQ()
  
  bootresults2 <- data.frame(bootresults1)
  
  ci0<-c()
  for(i in 1:12){
    ci1<-quantile(unlist(bootresults2[i,]), probs = c(0.05, 0.95))
    ci0<-rbind(ci0,ci1)
  }
  
  
  results2<-list(mean_dissimilarity=mean_dissimilarity,mean_dissimilarityci1=(unlist(ci0[1,1])),mean_dissimilarityci2=(unlist(ci0[1,2])),
                 mean_dissimilarity_direct=mean_dissimilarity_direct,mean_dissimilarity_directci1=unlist(ci0[2,1]),mean_dissimilarity_directci2=unlist(ci0[2,2]),
                 standardized_mean_dissimilarity=standardized_mean_dissimilarity,standardized_mean_dissimilarityci1=unlist(ci0[3,1]),standardized_mean_dissimilarityci2=unlist(ci0[3,2]),
                 standardized_mean_dissimilarity_direct=standardized_mean_dissimilarity_direct,standardized_mean_dissimilarity_directci1=unlist(ci0[4,1]),standardized_mean_dissimilarity_directci2=unlist(ci0[4,2]),
                 sd_dissimilarity=sd_dissimilarity,sd_dissimilarityci1=unlist(ci0[5,1]),sd_dissimilarityci2=unlist(ci0[5,2]),
                 sd_dissimilarity_direct=sd_dissimilarity_direct,sd_dissimilarity_directci1=unlist(ci0[6,1]),sd_dissimilarity_directci2=unlist(ci0[6,2]),
                 standardized_sd_dissimilarity=standardized_sd_dissimilarity,standardized_sd_dissimilarityci1=unlist(ci0[7,1]),standardized_sd_dissimilarityci2=unlist(ci0[7,2]),
                 standardized_sd_dissimilarity_direct=standardized_sd_dissimilarity_direct,standardized_sd_dissimilarity_directci1=unlist(ci0[8,1]),standardized_sd_dissimilarity_directci2=unlist(ci0[8,2]),
                 
                 sparsity_dissimilarity=sparsity_dissimilarity,sparsity_dissimilarityci1=unlist(ci0[9,1]),sparsity_dissimilarityci2=unlist(ci0[9,2]),
                 standardized_sparsity_dissimilarity=standardized_sparsity_dissimilarity,standardized_sparsity_dissimilarityci1=unlist(ci0[10,1]),standardized_sparsity_dissimilarityci2=unlist(ci0[10,2]),
                 sparsity_dissimilarity_direct=sparsity_dissimilarity_direct,sparsity_dissimilarity_directci1=unlist(ci0[11,1]),sparsity_dissimilarity_directci2=unlist(ci0[11,2]),
                 standardized_sparsity_dissimilarity_direct=standardized_sparsity_dissimilarity_direct,standardized_sparsity_dissimilarity_directci1=unlist(ci0[12,1]),standardized_sparsity_dissimilarity_directci2=unlist(ci0[12,2])
                 )
  
  
  rawdata<-list(meanmean=results$meanmean,
                 sdsd=results$sdsd,
                
                 sumall=as.data.frame(results$sumall),
                 meancombined3=as.data.frame(results$meancombined3),
                 meandi=as.data.frame(results$meandi),
                 sdcombined3=as.data.frame(results$sdcombined3),
                 sddi=as.data.frame(results$sddi),
                
                selected1zero=as.data.frame(results$selected1zero),
                selected2zero=as.data.frame(results$selected2zero)
                 )
  Results<-list(results=(results2),rawdata=rawdata)
  return(Results)
}

