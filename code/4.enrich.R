require(DOSE)
require(clusterProfiler)
require(org.Hs.eg.db)
library(ggplot2)
library(stringr)
library(GOplot)
mRNAlist<-list.files("F:\\datasetlist\\")
path<-"F:\\dataset\\"
for (i in 1:length(mRNAlist)) {
  ID<-mRNAlist[i]
  datapath<-str_c(path,ID,"\\nrDEG_",ID,".txt")
  print(str_c(i,ID,"start"))
  if(file.exists(datapath)){  
    DEGfile<-read.delim(datapath)
    pCutoff<-0.05
    resultpCutoff<-0.05
    genelist<-DEGfile[which(DEGfile$P.Value<pCutoff),1] 
    gene<-as.character(genelist) 
    #symbol2id
    ids <- bitr(gene,fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = "org.Hs.eg.db")
    genes = ids[,2]
    #GO
    GO_enrich<-enrichGO(gene=genes,
                        OrgDb = "org.Hs.eg.db",
                        keyType = "ENTREZID",
                        ont="ALL",
                        pvalueCutoff = 1,
                        qvalueCutoff = 0.5,
                        readable = T) 
    go_all<-data.frame(GO_enrich@result)
    if(nrow(go_all)==0){
      print(str_c(ID,"not enrich"))
      KEGG_enrich <- enrichKEGG(gene = genes,
                                organism = 'hsa',
                                pAdjustMethod = "BH",
                                pvalueCutoff = 1,
                                qvalueCutoff = 0.5
      )
      kegg_symbol <- setReadable(KEGG_enrich,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
      kegg_all<-data.frame(ONTOLOGY="KEGG",kegg_symbol@result)
      kegg_allresult<-data.frame(kegg_all[,1:5],signif(kegg_all[,6:8],4),kegg_all[,9:10])
      keggresult_path<-str_c(path,ID,"\\KEGGenrich_", ID, ".txt")
      write.table(kegg_allresult,keggresult_path,row.names = F,quote = F,sep = "\t")
      print(str_c(ID," KEGG done"))
    }else{
      go_allresult<-data.frame(go_all[,1:5],signif(go_all[,6:8],4),go_all[,9:10])
      result_path<-str_c(path,ID,"\\GOenrich_", ID, ".txt")
      write.table(go_allresult,result_path,row.names = F,quote = F,sep = "\t")
      print("GO done")
      KEGG_enrich <- enrichKEGG(gene = genes,
                                organism = 'hsa',
                                pAdjustMethod = "BH",
                                pvalueCutoff = 1,
                                qvalueCutoff = 0.5)
      kegg_symbol <- setReadable(KEGG_enrich,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
      kegg_all<-data.frame(ONTOLOGY="KEGG",kegg_symbol@result)
      kegg_allresult<-data.frame(kegg_all[,1:5],signif(kegg_all[,6:8],4),kegg_all[,9:10])
      keggresult_path<-str_c(path,ID,"\\KEGGenrich_", ID, ".txt")
      write.table(kegg_allresult,keggresult_path,row.names = F,quote = F,sep = "\t")
      print("KEGG done")
    }
  }else{
    print(str_c(ID," not nrDEG"))
  }
}

  
  


