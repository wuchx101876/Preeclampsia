

remove(list = ls())
setwd("~/BioXCG/income/Service/2024/2024.04.22/")

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



remove(list = ls())

exp <- readRDS("./raw_data/GSE48424/GSE48424_exp.rds")
cli <- readRDS("./raw_data/GSE48424/GSE48424_cli.rds")
cli1 <- cli %>% dplyr::select(geo_accession,`disease status:ch1`)
colnames(cli1) <- c("sample","type")
use_data <- data.frame(t(exp["BACH1",]))
use_data$sample <- rownames(use_data)

use_data1 <- dplyr::inner_join(use_data,cli1)
use_data1$type <- factor(use_data1$type,levels = names(table(use_data1$type)))

ggplot(data=use_data1,aes(x=type,y=BACH1,fill=type))+
  geom_boxplot()+
  stat_compare_means(aes(label = ..p.signif..),method = "wilcox.test",hide.ns = F)+
  theme_prism()+
  theme(legend.position = "none")+
  labs(y="BACH1 expression",x= NULL,title = "GSE48424")+
  scale_x_discrete(labels = c("Normal","PE"))


remove(list = ls())
model_data <- readRDS("./result/ML/ML_usedata.rds")

exp <- readRDS("./raw_data/GSE48424/GSE48424_exp.rds")
use_exp <- as.data.frame(t(exp[colnames(model_data)[3:13],]))
use_exp$sample <- rownames(use_exp)

cli <- readRDS("./raw_data/GSE48424/GSE48424_cli.rds")
cli1 <- cli %>% dplyr::select(geo_accession,`disease status:ch1`)
colnames(cli1) <- c("sample","status")

use_clin <- cli1[,1:2]
use_data <- dplyr::inner_join(use_clin,use_exp)
use_data <- use_data %>% dplyr::mutate(status = ifelse(status == "healthy",0,1))
saveRDS(use_data,"./result/Validation/GSE48424_model_data.rds")


### 

use_data <- readRDS("./result/Validation/GSE48424_model_data.rds")

for (i in colnames(use_data)[3:13]) {
  
  df <- use_data %>% dplyr::select(status,i)
  colnames(df) <- c("status","gene_expression")

  roc_obj <- roc(
    response = df$status,  
    predictor = df$gene_expression,  
    levels = c(0, 1),  
    direction = "auto"
  )
  

  p <- ggroc(roc_obj, color = "red", size = 1.2) + 
    geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), 
                 color = "grey", linetype = "dashed") +  
    labs(title = paste0(i," AUC = ", 
                        round(auc(roc_obj), 3)), 
         x = "1 - Specificity",
         y = "Sensitivity") +
    ggprism::theme_prism(border = T)+
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    ) + coord_equal() 
  
  
  pdf(file.path(paste0("./result/BACH1/",i,".pdf")),
      width = 8, height = 6)
  print(p, newpage = FALSE)
  dev.off()
  
}

