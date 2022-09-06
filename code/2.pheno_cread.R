gse<-"GSE45291"
exprSet<-read.delim(paste(gse,"txt",sep = "."))
data<-getGEO(gse,getGPL = F)
pheno<-pData(data[[1]])

unique(pheno$`clinical stage:ch1`)
#move<-pheno$geo_accession[which(pheno$source_name_ch1=="synovial tissue from osteoarthritic joint")]
move<-pheno$geo_accession[which(pheno$`clinical stage:ch1`=="Relative in preclinical autoimmune stage" )]

exp<-exprSet[,-which(match(names(exprSet),move)!="NA")]
phenotype<-data.frame(gsm=pheno$geo_accession,pheno=pheno$`clinical stage:ch1`)
phenotype$phenotype<-2
phenotype$phenotype[which(phenotype$pheno=="Healthy relatives")]<-0
phenotype$phenotype[which(phenotype$pheno=="RA patient, less than one year of evolution untreated")]<-1
phenodata<-data.frame(gsm=phenotype$gsm,phenotype=phenotype$phenotype)[-which(phenotype$phenotype==2),]
dir.create(gse)
pheno_file<-paste(gse,paste("phenotype","csv",sep = "."),sep = "/")
expr_file<-paste(gse,paste("exprSet","txt",sep = "."),sep = "/")
names(phenodata)<-c("sample","status")
write.table(phenodata,pheno_file,row.names = F,col.names =T,quote = F,sep = ",")
write.table(exp,expr_file,row.names = F,quote = F,sep = "\t")

length(which(phenodata$phenotype==0))
length(which(phenodata$phenotype==1))

###
gse<-"GSE93777"
exp<-read.delim(paste(gse,"txt",sep = "."))
path<-"E:\\"
dir.create(gse)
expr_file<-paste(gse,paste("exprSet","txt",sep = "."),sep = "/")
write.table(exp,expr_file,row.names = F,quote = F,sep = "\t")

phe<-paste(path,paste(paste(gse,"1",sep = "_"),"txt",sep = "."),sep="\\")
phenodata<-read.table(phe)
names(phenodata)<-c("sample","status")
pheno_file<-paste(gse,paste("phenotype","csv",sep = "."),sep = "/")
write.table(phenodata,pheno_file,row.names = F,col.names =T,quote = F,sep = ",")

