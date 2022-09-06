## library the packages
library(tidyverse)
library(tibble)
library(GEOquery)
library(affy)
library(limma)
library(idmap3)

# deal the agilent microarray raw data
gse<-"GSE78068"
gpl<-"GPL6480"

# start
rawname<-paste(gse,"RAW.tar",sep ="_")
rawpath<-"E:\\raw"
file <- paste(rawpath,rawname,sep = "\\")
exdir <- file %>% stringr::str_extract("GSE\\d+")
untar(file, exdir = exdir)
gz <- list.files(path = exdir, pattern = "*.gz$")
sapply(paste(exdir, gz, sep = "/"), gunzip)
files <- list.files(exdir, pattern = "GSM", full.names = T)
RawData <- read.maimages(files = files, source = "agilent", green.only = T)
EListRaw <- RawData
unlink(exdir, recursive = TRUE)

# data were normalized uniformly using the RMA method.
BgCorrect <- backgroundCorrect(EListRaw, method = "normexp", normexp.method = "rma", offset = 50)##数据标准�?
BGandNormalized <- normalizeBetweenArrays(BgCorrect, method = "quantile")
EList <- BGandNormalized
averEList <- avereps(EList, ID = EList$genes$ProbeName)
exp <- averEList$E
row.names(exp) <- averEList$genes$ProbeName
colnames(exp) <- colnames(exp) %>% str_extract("GSM\\d+")
#head(exp)
exprSet <- exp %>%
  as.data.frame() %>%
  rownames_to_column("probe_id")
print(dim(exprSet)[2])
head(exprSet)
write.table(exprSet,paste("renew/probe",paste(gse,"txt",sep = "."),sep = "/"),row.names = F,sep = "\t",quote = F)

## annotate the probes to gene symbol according to platform-specific probe annotation files.
file <- "GSE78068_family.soft.gz"
anno <- getGEO(filename = file, getGPL = F)
anno <- anno@gpls$GPL6480@dataTable@table
probe2symbol_df <- anno %>% dplyr::select(c("ID", "Gene Symbol"))
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
