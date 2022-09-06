RABC_keggplot <- function(ID,shownum,width,height,path,picpath) {
  print(path)
  #ggplot bar
  datapath<-str_c(path,ID,"/KEGGenrich_", ID, ".txt")
  print(datapath)
  go_all<-read.delim(datapath)
  go_enrich_df<-go_all[go_all$ONTOLOGY=="KEGG",][1:shownum,]
  go_enrich_df$Description <- factor(go_enrich_df$Description,levels=unique(go_enrich_df$Description))
  
  go_enrich_df$number <- factor(rev(1:nrow(go_enrich_df)))
  go_enrich_df <- na.omit(go_enrich_df)
  p1 <- ggplot(data=go_enrich_df, aes(x=Description, y=Count)) +
    geom_bar(stat="identity", width=0.8, fill="#d77a61") + coord_flip() +
    theme_test() +
    xlab("KEGG term") +
    theme(axis.text=element_text(face = "bold", color="gray50")) +
    labs(title = "RABC Enriched KEGG Terms ",ID)
  bar_path<-str_c(picpath,"/keggbar_",ID,".png")
  ggsave(plot = p1, filename = bar_path, width = 15, height =10 )
  #
  go_enrich_df$genenum = as.numeric(str_split_fixed(go_enrich_df$GeneRatio,'/',2)[,2])
  p2<-ggplot(go_enrich_df,aes(Count/genenum,Description))+
    geom_point(aes(size=Count,color=p.adjust))+
    scale_colour_gradient(low="#718355",high="red",
                          guide = guide_colorbar(reverse = TRUE))+
    labs(x="GeneRatio",title = str_c("RABC Enriched KEGG_",ID))+
    theme_bw()+
    theme(axis.text = element_text(size=12),
          axis.title = element_text(size = 14))+
    facet_grid(ONTOLOGY~., scale="free")
  gobubble_path<-str_c(picpath,"/keggbubble_",ID,".png")
  ggsave(plot = p2, filename = gobubble_path, width = 18, height =10 )
}

library(ggplot2)
library(stringr)

args=commandArgs(T)
RABC_keggplot(ID= args[1], shownum=args[2],width=args[3],height=args[4],path=args[5], picpath=args[6])
