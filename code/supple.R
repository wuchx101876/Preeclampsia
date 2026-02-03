

remove(list = ls())

## GSE149437

library(clusterProfiler) 
library(org.Hs.eg.db) 
library(stats)
library(data.table) 
library(dplyr)
library(AnnoProbe)
library(GEOquery) 

geo_data <- getGEO(
  filename  = "./raw_data/GSE149437/GSE149437_series_matrix.txt.gz",
  GSEMatrix = TRUE,
  getGPL    = FALSE
)


expr  <- exprs(geo_data)   
pd    <- pData(geo_data)   


ids <- fread("./GPL28460.txt",data.table = F)


gene.df <- bitr(ids$ENTREZ_GENE_ID, fromType = "ENTREZID",
                toType = c("SYMBOL","ENTREZID"), 
                OrgDb = org.Hs.eg.db) 
gene.df$ENTREZID <- as.integer(gene.df$ENTREZID)

colnames(gene.df)

ids <- ids %>% dplyr::select(ID,ENTREZ_GENE_ID)
colnames(ids)[2] <- "ENTREZID"
use_id <- dplyr::inner_join(ids,gene.df)
use_id <- use_id %>% dplyr::select(ID,SYMBOL)

exp <- data.frame(expr)
exp$ID <- rownames(exp)

exp1 <- dplyr::inner_join(use_id,exp)

library(limma)
exp2 <- exp1[,-1]
exp3 <- as.data.frame(limma::avereps(exp2[,-1],ID = exp2$SYMBOL))
saveRDS(exp3,"./raw_data/GSE149437/GSE149437_exp.rds")

cli <- pd %>% dplyr::select(geo_accession,characteristics_ch1.3)
table(cli$characteristics_ch1.3)
colnames(cli) <- c("sample","type")
saveRDS(cli,"./raw_data/GSE149437/GSE149437_cli.rds")



## AUC
remove(list = ls())

exp <- readRDS("./raw_data/GSE149437/GSE149437_exp.rds")
BACH1 <- data.frame(t(exp["BACH1",]))
BACH1$sample <- rownames(BACH1)
cli <- readRDS("./raw_data/GSE149437/GSE149437_cli.rds")

use_data <- dplyr::inner_join(BACH1,cli)



df <- readRDS("./raw_data/GSE149437/use_data.rds")

library(ggplot2)
library(ggpubr)
library(ggprism)


ggplot(data = df, aes(x = type, y = BACH1 ,fill = type)) +
  geom_boxplot() +
  stat_compare_means(size = 5, comparisons = list(c("Normal","PE")),label = "p.format",
                     method = "wilcox.test") +
  
  theme_prism(border = TRUE) +
  scale_x_discrete(labels = c("Normal","PE")) +
  theme(axis.title.x = element_blank())+
  scale_fill_manual(values = c("#66C2A5","#FC8D62"))




df <- use_data
df$status <- ifelse(df$type == "PE",1,0)

library(pROC)
roc_obj <- roc(
  response = df$status,  
  predictor = df$BACH1,  
  levels = c(0, 1),  
  direction = "auto"
)


p <- ggroc(roc_obj, color = "red", size = 1.2) + 
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), 
               color = "grey", linetype = "dashed") +  
  labs(title = paste0("BACH1"," AUC = ", 
                      round(auc(roc_obj), 3)," (GSE149437)"), 
       x = "1 - Specificity",
       y = "Sensitivity") +
  ggprism::theme_prism(border = T)+
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  ) + coord_equal() 

p






