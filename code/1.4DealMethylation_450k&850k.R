library(minfi)
library(wateRmelon)
path<-"F:/"
filenames<-list.files(path)

for(i in seq(1,length(filenames),2)){
  file<-paste(path,filenames[i],sep="/")
  print(file)
  data1 <- read.metharray(file, extended = TRUE, verbose = TRUE) # extended needed for get NBeads
  dim(getRed(data1))
  dim(getGreen(data1))
  manifest <- getManifest(data1)
  manifest
  head(getProbeInfo(manifest))
  head(getRed(data1))
  
  ###sample filter by detectionP
  detP <- detectionP(data1)
  Pfailed <- detP>0.01
  colMeans(Pfailed) # Fraction of failed positions per sample
  sum(colMeans(Pfailed)>0.01)
  if(sum(colMeans(Pfailed)>0.01)>0){
    data1_cutSam <- data1
    result_name<-substr(filenames[i],1,10)
    result_file<-paste(paste("F:/GSE87571_result",result_name,sep="/"),"_cutSample.txt",sep="")
  }else{
    cutSample <- colMeans(Pfailed)>0.01
    data1_cutSam <- data1[,!cutSample]
    result_name<-substr(filenames[i],1,10)
    result_file<-paste(paste("F:/result",result_name,sep="/"),".txt",sep="")
  }
  ###find failed probe by detectionP
  detP <- detectionP(data1_cutSam)
  Pfailed <- detP>0.01
  sum(rowMeans(Pfailed)>0.1) # How many positions failed in >10% of samples?
  failedProbesP <-rownames(Pfailed)[rowMeans(Pfailed)>0.1]
  
  ###find failed probe by Nbeads use wateRmelon packages
  
  beadcount <- beadcount(data1_cutSam)
  NBfailed <- is.na(beadcount)
  sum(rowMeans(NBfailed)>0.05)
  failedProbesNB <-rownames(NBfailed)[rowMeans(NBfailed)>0.05]
  
  ###get the list of failed probes
  failedProbes <- union(failedProbesP, failedProbesNB) 
  
  ###preprocess and get beta value
  Mset <-preprocessIllumina(data1_cutSam, bg.correct=TRUE, normalize = c("controls"), reference =1) 
  Mset
  
  ###remove failed probes
  Mset_cutProbe <- Mset[!rownames(Mset)%in% failedProbes,] 
  head(getMeth(Mset_cutProbe))
  
  ###get beta value
  RSet <- ratioConvert(Mset_cutProbe, what = "both", keepCN = TRUE)
  RSet
  beta <- getBeta(RSet)
  
  write.table(beta,file=result_file,sep="\t",row.names=TRUE,col.names=TRUE)
}