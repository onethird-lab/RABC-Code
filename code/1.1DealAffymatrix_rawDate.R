## library the packages
library(tidyverse)
library(tibble)
library(GEOquery)
library(affy)
library(limma)
library(idmap3)

## deal the affymatrix microarray raw data
gse<-"GSE10500"
gpl<-"GPL8300"

# start
rawname<-paste(gse,"RAW.tar",sep ="_")
rawpath<-"E:\\raw"
file <- paste(rawpath,rawname,sep = "\\")
exdir <- rawname %>% stringr::str_extract("GSE\\d+")
untar(file, exdir = exdir)
gz <- list.files(path = exdir, pattern = "*.gz$")
sapply(paste(exdir, gz, sep = "/"), gunzip)##if you are first use gunzip  you need to install.packages("R.utils")
# read data
cels <- affy::list.celfiles(exdir, full.names = TRUE)
affy_data <- affy::ReadAffy(filenames = cels)
unlink(exdir, recursive = TRUE)
# data were normalized uniformly using the RMA method
affy_rma <- affy::rma(affy_data)
exprSet <- exprs(affy_rma) %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "probe_id")
exprSet$probe_id <- as.character(exprSet$probe_id)
names(exprSet)[-1] <- names(exprSet)[-1] %>% str_extract("GSM\\d+")
print(dim(exprSet)[2])
head(exprSet)
write.table(exprSet,paste("renew/probe",paste(gse,"txt",sep = "."),sep = "/"),row.names = F,sep = "\t",quote = F)

## annotate the probes to gene symbol according to platform-specific probe annotation files.
file <- "GSE10500_family.soft.gz"
anno <- getGEO(filename = file, getGPL = F)
anno <- anno@gpls$GPL8300@dataTable@table
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

