## ID <- "GSE120178"
## gene1 <- "TP53"
## gene2 <- "JAK2"
## method <- "pearson/spearman"

FUN_corrplot <- function(ID, gene1, gene2, method = "pearson", addNormal = F,datapath,picpath) {
  proj <- ID
  exprSet <- fread(str_c(datapath,"/exprSet.txt"), data.table = F) %>% column_to_rownames(names(.)[1])
  pd <- read.csv(str_c(datapath,"/phenotype.csv"))
  pd$status <- ifelse(pd$status == 1, "RA", "Normal")
  sel <- which(pd$status == "RA")
  if (!gene1 %in% rownames(exprSet)) {
    stop(gene1, " is not in the study!")
  }
  if (!gene2 %in% rownames(exprSet)) {
    stop(gene2, " is not in the study!")
  }
  if (addNormal == F) {
    temp <- exprSet[c(gene1, gene2), sel]
  } else {
    temp <- exprSet[c(gene1, gene2), ]
  }
  df <- temp %>%
    t() %>%
    as.data.frame()
  p <- ggplot(data = df, aes(x = df[, 1], y = df[, 2])) +
    labs(x = gene1, y = gene2) +
    geom_point(color = "#C45A44", size = 2) +
    geom_smooth(method = lm, lwd = 1.5, se = T, color = "#C45A44") +
    geom_rug(color = "#C45A44") +
    ggtitle(ID) +
    theme_classic2() +
    theme(
      plot.title = element_text(size = 18, colour = "black"),
      axis.title = element_text(size = 18),
      axis.text.y = element_text(size = 12, colour = "black"),
      axis.text.x = element_text(size = 12, colour = "black")
    ) +
    stat_cor(method = method, aes(x = df[, 1], y = df[, 2]), size = 5.5)
  ggsave(plot = p, filename = picpath, width = 4.5, height = 4)
  #return(p)
}
library(tidyverse)
library(data.table)
library(magrittr)
library(tibble)
library(ggpubr)
#FUN_corrplot(ID = "GSE120178", gene1 = "TP53", gene2 = "JAK2")
args=commandArgs(T)
FUN_corrplot(ID = args[1], gene1 = args[2],gene2 =args[3],method =args[4],addNormal =args[5],datapath=args[6],picpath = args[7])
