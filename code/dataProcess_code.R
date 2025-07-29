
## GSE75010
remove(list = ls())
library(AnnoProbe) 
library(GEOquery) 
geo_data <- getGEO(GEO = "GSE75010", destdir = './raw_data/GSE75010/', getGPL = F)
geo_data0 <- geo_data[[1]]
exp <- as.data.frame(Biobase::exprs(geo_data0)) 
gpl_number <- geo_data0@annotation
ids <- AnnoProbe::idmap(gpl_number)
exp$probe_id <- rownames(exp)
exp1 <- dplyr::inner_join(ids,exp)
library(limma)
exp2 <- exp1[,-1]
exp3 <- as.data.frame(limma::avereps(exp2[,-1],ID = exp2$symbol))
saveRDS(exp3,"./raw_data/GSE75010/GSE75010_exp.rds")
cli <- Biobase::pData(geo_data0)
saveRDS(cli,"./raw_data/GSE75010/GSE75010_cli.rds")



### GSE25906
remove(list = ls())
library(AnnoProbe) 
library(GEOquery) 
geo_data <- getGEO(GEO = "GSE25906", destdir = './raw_data/GSE25906/', getGPL = F)
geo_data0 <- geo_data[[1]]
exp <- as.data.frame(Biobase::exprs(geo_data0)) 
gpl_number <- geo_data0@annotation
ids <- AnnoProbe::idmap(gpl_number)
exp$probe_id <- rownames(exp)
exp1 <- dplyr::inner_join(ids,exp)
library(limma)
exp2 <- exp1[,-1]
exp3 <- as.data.frame(limma::avereps(exp2[,-1],ID = exp2$symbol))
saveRDS(exp3,"./raw_data/GSE25906/GSE25906_exp.rds")
cli <- Biobase::pData(geo_data0)
saveRDS(cli,"./raw_data/GSE25906/GSE25906_cli.rds")
cli <- readRDS("./raw_data/GSE25906/GSE25906_cli.rds")
use_clin <- cli %>% dplyr::select(geo_accession,"classification:ch1","gender:ch1","gestational age:ch1")
colnames(use_clin) <- c("sample","status","gender","MMage")
use_clin1 <- use_clin
use_clin1$MMage <- as.numeric(gsub(".*: ", "", use_clin1$MMage))
use_clin1 <- use_clin1 %>% dplyr::mutate(status = if_else(status == "preeclamptic","PE","NonPE"))
saveRDS(use_clin1,"./use_data/GSE25906_clin.rds")




## GSE48424
library(AnnoProbe)
library(GEOquery) 
geo_data <- getGEO(GEO = "GSE48424", destdir = './raw_data/GSE48424/', getGPL = F)
geo_data0 <- geo_data[[1]]
exp <- as.data.frame(Biobase::exprs(geo_data0)) 
ids <- fread("./raw_data/GSE48424/GPL6480-9577.txt",data.table = F)
ids1 <- ids[,c(1,7)]
colnames(ids1) <- c("probe_id","symbol")
exp$probe_id <- rownames(exp)
exp1 <- dplyr::inner_join(ids1,exp)
library(limma)
exp2 <- exp1[,-1]
exp3 <- as.data.frame(limma::avereps(exp2[,-1],ID = exp2$symbol))
saveRDS(exp3,"./raw_data/GSE48424/GSE48424_exp.rds")
cli <- Biobase::pData(geo_data0)
saveRDS(cli,"./raw_data/GSE48424/GSE48424_cli.rds")


