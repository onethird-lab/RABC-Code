RABC_goplot <- function(ID,shownum,width,height,path,picpath) {

  #ggplot bar
  datapath<-str_c(path,ID,"/GOenrich_", ID, ".txt")
  go_all<-read.delim(datapath)
  go_MF<-go_all[go_all$ONTOLOGY=="MF",][1:shownum,]
  go_CC<-go_all[go_all$ONTOLOGY=="CC",][1:shownum,]
  go_BP<-go_all[go_all$ONTOLOGY=="BP",][1:shownum,]
  go_enrich_df<-rbind(go_BP,go_CC,go_MF)
  go_enrich_df$Description <- factor(go_enrich_df$Description,levels=unique(go_enrich_df$Description))
  go_enrich_df$type<-factor(c(rep("biological process", shownum), rep("cellular component", shownum),rep("molecular function",shownum)),levels=c("molecular function", "cellular component", "biological process"))
  go_enrich_df$number <- factor(rev(1:nrow(go_enrich_df)))
  go_enrich_df <- na.omit(go_enrich_df)
  CPCOLS <- c("#3d405b", "#81b29a", "#f2cc8f")
  p1 <- ggplot(data=go_enrich_df, aes(x=Description, y=Count, fill=type)) +
    geom_bar(stat="identity", width=0.8) + coord_flip() +
    scale_fill_manual(values = CPCOLS) + theme_test() +
    xlab("GO term") +
    theme(axis.text=element_text(face = "bold", color="gray50")) +
    labs(title = "RABC Enriched GO Terms ",ID)
  bar_path<-str_c(picpath,"/gobar_",ID,".png")
  ggsave(plot = p1, filename = bar_path, width = 15, height =10 )
  #
  go_enrich_df$genenum = as.numeric(str_split_fixed(go_enrich_df$GeneRatio,'/',2)[,2])
  p2<-ggplot(go_enrich_df,aes(Count/genenum,Description))+
    geom_point(aes(size=Count,color=p.adjust))+
    scale_colour_gradient(low="#718355",high="red",
                          guide = guide_colorbar(reverse = TRUE))+
    labs(x="GeneRatio",title = str_c("RABC Enriched GO_",ID))+
    theme_bw()+
    theme(axis.text = element_text(size=12),
          axis.title = element_text(size = 14))+
    facet_grid(ONTOLOGY~., scale="free")
  gobubble_path<-str_c(picpath,"/gobubble_",ID,".png")
  ggsave(plot = p2, filename = gobubble_path, width = 18, height =10 )
  }

require(DOSE)
require(clusterProfiler)
require(org.Hs.eg.db)
library(ggplot2)
library(stringr)
library(GOplot)
args=commandArgs(T)
RABC_goplot(ID= args[1], shownum=args[2],width=args[3],height=args[4],path=args[5], picpath=args[6])
