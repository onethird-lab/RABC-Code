FUN_volcano <- function(ID, adjP = F, pCutoff = 0.05, logFCcutoff = 1,datapath,picpath) {
  pCutoff<-as.numeric(pCutoff)
  logFCcutoff<-as.numeric(logFCcutoff)
  res <- fread(datapath,data.table = F)
  if (adjP == F) {
    y <- "P.Value"
    ylab <- bquote(~-Log[10]~italic("P-value"))
    #print("p")
  } else {
    y <- "adj.P.Val"
    ylab <- bquote(~-Log[10]~italic("adj.P-value"))
    #print("adj")
  }

  p <- EnhancedVolcano(res,
                       lab = rep("", nrow(res)),
                       x = "logFC",
                       y = y,
                       ylab = ylab,
                       title = ID,
                       pCutoff = pCutoff,
                       FCcutoff = logFCcutoff,
                       transcriptPointSize = 2
  ) +
    theme(title = element_blank(), plot.subtitle = element_blank(), plot.caption = element_blank(), legend.position = "none")
  ggsave(picpath, width = 5, height = 5)
  return(p)
}
library(tidyverse)
library(data.table)
library(EnhancedVolcano)
args=commandArgs(T)
FUN_volcano(ID= args[1], adjP=args[2], pCutoff=args[3], logFCcutoff=args[4],datapath=args[5],picpath = args[6])
