library(tidyverse)
library(data.table)
library(magrittr)
library(tibble)
library(limma)
setwd("E:\\RABC-dataset\\RABC-Methylation\\RA-healthy")
l <- list.dirs()
l<-l[c(4,5)]
l <- l[l %>% str_detect("GSE")]
FUN <- function(dir) {
  proj <- dir %>% str_extract("GSE\\d+")
  expr <- fread(str_c(dir, "/exprSet.txt"), data.table = F) %>% column_to_rownames(names(.)[1])
  pd <- read.csv(str_c(dir, "/phenotype.csv"))
  pd$status <- ifelse(pd$status == 1, "RA", "Normal")
  group_list <- factor(pd$status, levels = c("RA", "Normal"))
  
  design <- model.matrix(~ 0 + group_list)
  colnames(design) <- levels(group_list)
  fit <- lmFit(expr, design)
  contrast.matrix <- makeContrasts(RA - Normal,
                                   levels = design
  )
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  tempOutput <- topTable(fit2, coef = 1, n = Inf)
  nrDEG <- na.omit(tempOutput)
  nrDEG <- nrDEG %>%
    rownames_to_column("symbol") %>%
    distinct()
  head(nrDEG)
  fwrite(nrDEG, file = str_c("nrDEG_", proj, ".txt"), row.names = F, quote = F, sep = "\t")
}
res <- lapply(l, FUN = FUN)
dir<-choose.dir()
FUN("E:\\RABC-dataset\\RABC-aff\\renew\\geneSymbol\\RA-healthy\\GSE45291")
