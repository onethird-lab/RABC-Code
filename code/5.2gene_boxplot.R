
FUN_box <- function(ID, gene, method= "t.test",datapath,picpath) {
  proj <- ID
  exprSet <- fread(str_c(datapath,"/exprSet.txt"), data.table = F) %>% column_to_rownames(names(.)[1])
  pd <- read.csv(str_c(datapath, "/phenotype.csv"))
  pd$status <- ifelse(pd$status == 1, "RA", "Normal")
  group_list <- factor(pd$status, levels = c("RA", "Normal"))
  if(!gene %in% rownames(exprSet)) {
    stop(str_c(gene, " is not in the study!"))
  }
  temp <- exprSet[gene, ]
  df <- exprSet %>%
    t() %>%
    as.data.frame()
  df$type <- factor(group_list, levels = c("RA", "Normal"))
  p <- ggplot(df, aes(x = type, y = df[, gene], fill = type)) +
    stat_boxplot(geom = "errorbar", width = 0.2) +
    geom_boxplot(
      width = 0.5, size = 0.8,
      position = position_dodge(width = 0.5),
      outlier.size = 0
    ) +
    geom_jitter(width = 0.05, size = 1, alpha = 0.6) +
    scale_fill_manual(values = c("#D62B2B", "#7498C4")) +
    theme_bw() +
    geom_signif(
      comparisons = list(c("RA", "Normal")), map_signif_level = T,
      test = method, textsize = 9, vjust = 0.6, color = "black"
    ) + # 或者换成t检验，就是t.test
    guides(fill = FALSE) +
    xlab(NULL) +
    ylab(str_c("Expression of ", gene)) +
    ggtitle(ID) +
    theme(
      plot.title = element_text(size = 18, colour = "black"),
      axis.text.y = element_text(size = 12, colour = "black"),
      axis.text.x = element_text(size = 18, colour = "black"),
      axis.title.y = element_text(size = 18, colour = "black"),
      panel.grid = element_blank()
    )
  ggsave(plot = p, filename = picpath, width = 4.5, height = 4)
  return(p)
}
library(tidyverse)
library(data.table)
library(magrittr)
library(tibble)
library(ggsignif)
args=commandArgs(T)
FUN_box(ID = args[1],gene=args[2],method =args[3],datapath=args[4],picpath=args[5])
