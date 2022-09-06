## library the packages
library(tidyverse)
library(tibble)
library(GEOquery)
library(affy)
library(limma)
library(idmap3)

## deal the illumina microarray data
gse<-"GSE49604"
gpl<-"GPL10558"
x <- read.ilmn("GSE49604_GPL8432_non-normalized.txt",expr = "SAMPLE",probeid="ID_REF")
y <- neqc(x)
EListRaw<-x
# data were normalized uniformly using the RMA method
BgCorrect <- backgroundCorrect(EListRaw, method = "normexp", normexp.method = "rma", offset = 50)
BGandNormalized <- normalizeBetweenArrays(BgCorrect, method = "quantile")
EList <- BGandNormalized
dim(EList$E)
averEList <- avereps(EList, ID = EList$genes$ProbeName)
expy<-y$E
exp<-averEList$E
dim(exp)
dim(expy)
GSE<-getGEO(gse,getGPL = F)
pheno<-GSE$`GSE49604-GPL8432_series_matrix.txt.gz`@assayData$exprs
dim(pheno)
colnames(exp) <- colnames(pheno) 
#head(exp)
exprSet <- exp %>%
  as.data.frame() %>%
  rownames_to_column("probe_id")
dim(exprSet)
write.table(exprSet,paste("renew/probe",paste(gse,"txt",sep = "."),sep = "/"),row.names = F,sep = "\t",quote = F)

## annotate the probes to gene symbol according to platform-specific probe annotation files.
file <- "E:\\renew_data\\data\\GSE49604_family.soft.gz"
anno <- getGEO(filename = file, getGPL = F)
anno <- anno@gpls$GPL10558@dataTable@table
rm('probe2symbol_df')
probe2symbol_df <- anno %>% dplyr::select(c("ID", "Symbol"))
names(probe2symbol_df) <- c("probe_id", "symbol")
head(probe2symbol_df)
probe2symbol_df <- probe2symbol_df %>%
  dplyr::filter(!symbol %>% str_detect("///")) %>%
  dplyr::filter(symbol %>% str_detect("\\w+"))
head(probe2symbol_df)
exprSet_gene <- exprSet %>%
  dplyr::inner_join(probe2symbol_df,
                    by = "probe_id"
  ) %>%
  dplyr::select(-"probe_id") %>%
  dplyr::select("symbol", dplyr::everything())
names.col <- names(exprSet_gene)
exprSet <- aggregate(exprSet_gene[, -1],
                     by = list(exprSet_gene$symbol),
                     FUN = median
) %>% as.data.frame()
names(exprSet) <- names.col
head(exprSet)
dim(exprSet)
