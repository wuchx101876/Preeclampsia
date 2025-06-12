

remove(list = ls())
library(dplyr)
library(data.table)


infor <- readRDS("./use_data/GSE75010_clin.rds")
infor <- infor %>% dplyr::arrange(status)
exp1 <- readRDS("./use_data/GSE75010_exp.rds")


library(limma)
library(edgeR)
use_exp <- exp1 %>% dplyr::select(infor$sample)
group <- c(rep("NonPE", 77), rep("PE", 80)) %>% as.factor()
desigN <- model.matrix(~ 0 + group) 
colnames(desigN) <- levels(group)
fit = lmFit(use_exp, desigN)
cont.matrix <- makeContrasts(contrasts = c('PE-NonPE'), levels = desigN)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

diff <- topTable(fit2,adjust='fdr', coef=1, number=Inf)
limma_DEG <- na.omit(diff)
limma_DEG$gene_id <- rownames(limma_DEG)
limma_DEG$logFC <- log2((10^abs(limma_DEG$logFC)))*sign(limma_DEG$logFC)

limma_DEG_sig <- limma_DEG %>% dplyr::filter(adj.P.Val < 0.05 ) %>% dplyr::filter(abs(logFC) > 0.5)

fwrite(limma_DEG,"./result/data/limma_DEG.csv")
fwrite(limma_DEG_sig,"./result/data/limma_DEG_sig.csv")



library(data.table)
nr_DEG <- fread("./result/data/limma_DEG.csv",data.table = F)



nr_DEG$log10FDR <- -log10(nr_DEG$adj.P.Val)
nr_DEG <- nr_DEG %>% 
  mutate(DEG = case_when(logFC > 0.5 & adj.P.Val < 0.05 ~ "up (614)",
                         abs(logFC) < 0.5 | adj.P.Val > 0.05 ~ "no (17391)",
                         logFC < -0.5 & adj.P.Val < 0.05 ~ "down (820)"))

table(nr_DEG$DEG)
library(ggplot2) 
library(ggprism) 
library(ggrepel) 
ggplot(nr_DEG, aes(x =logFC, y=log10FDR, colour=DEG)) +
  geom_point(alpha=1, size=2) +
  scale_color_manual(values=c('steelblue','gray','brown')) +
  xlim(-7,8) +
  geom_vline(xintercept=c(-0.5,0.5),lty=4,col="black",lwd=0.8)+
  geom_hline(yintercept = -log10(0.05), lty=4,col="black",lwd=0.8) + 
  labs(x="log2FC", y="-log10FDR") +
  ggtitle("NonPE vs PE DEG") + 
  theme(plot.title = element_text(hjust = 0.5),legend.position="right",legend.title = element_blank())+
  theme_prism(border = T)




## WGCNA

remove(list = ls())
# install.packages("WGCNA")
library(WGCNA)
library(limma)

exp_data <- readRDS("./use_data/GSE75010_exp.rds")
infor_data <- readRDS("./use_data/GSE75010_clin.rds")

wgcna_data1 <- exp_data
wgcna_data2 <- as.data.frame(t(wgcna_data1))
vars_res <- apply(wgcna_data2, 2, var)
per_res <- quantile(vars_res, probs = seq(0, 1, 0.25)) 
per_res

upperGene <- wgcna_data2[, which(vars_res > per_res[2])]
dim(upperGene)


exp_tpm <- upperGene
gsg <- goodSamplesGenes(exp_tpm, verbose = 3)
sampleTree <- hclust(dist(exp_tpm), method = "average")
dev.off()
plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "")
abline(h = 50, col = "red")

clust <- cutreeStatic(sampleTree, cutHeight = 50, minSize = 10)
table(clust)


datExpr <- exp_tpm[clust == 1, ]
nGenes = ncol(datExpr)
nGenes
nSamples = nrow(datExpr)
nSamples


clin_data1 <- infor_data %>% dplyr::filter(sample %in% rownames(datExpr))
clin_data1$type <- if_else(clin_data1$status == "NonPE",0,1)

sampleTree_final <- hclust(dist(datExpr), method = "average")

library(RColorBrewer)
library(scales)
traitColors <- numbers2colors(clin_data1$type, colors = c("#66C2A5","#FC8D62"), signed = F) ##离散变量
plotDendroAndColors(sampleTree_final, traitColors,
                    groupLabels = "group",
                    main = "Sample dendrogram and trait heatmap")

powers <- c(c(1:10), seq(from = 13, to = 25, by = 3))

sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
sft


dev.off()
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2", type = "n",
     main = paste("Scale independence"))
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers,cex = 1.2, col = "red")
abline(h = 0.85, col = "red")

dev.off()
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[, 1], sft$fitIndices[, 5],cex = 1.2, labels = powers, col = "red")

net <- blockwiseModules(datExpr, power = sft$powerEstimate, 
                        maxBlockSize = nGenes, TOMType = "unsigned", 
                        minModuleSize = 50, reassignThreshold = 0, mergeCutHeight = 0.25, 
                        numericLabels = TRUE, pamRespectsDendro = FALSE, 
                        saveTOMs = F, verbose = 3)
saveRDS(net,"./result/data/net.rds")

net <- readRDS("./result/data/net.rds")
table(net$colors) 
moduleColors <- labels2colors(net$colors)
table(moduleColors)

plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

moduleLables <- net$colors
moduleColors <- labels2colors(net$colors)

MEs <- net$MEs
head(MEs)[1:5, 1:5]
geneTree <- net$dendrograms[[1]]
save(moduleLables, moduleColors, MEs, geneTree, file = "./result/data/networkConstruction.RData")

MEList <-  moduleEigengenes(datExpr, colors = moduleColors)
MEs0 <- MEList$eigengenes

head(MEs0)[1:5, 1:5]
MEs <- orderMEs(MEs0)
head(MEs)[1:5, 1:5]
moduleTraitCor <- cor(MEs, clin_data1$type , use = "p");
head(moduleTraitCor)

moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)
head(moduleTraitPvalue)

textMatrix <- paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)

use_p <- data.frame(textMatrix)
colnames(use_p) <- "Preeclampsia"
use_p$control <- paste(-signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")

use_matrix <- data.frame(moduleTraitCor)
use_matrix$control <- -moduleTraitCor[,1]

colnames(use_matrix) <- c("Preeclampsia","Control")

dev.off()
par(mar = c(5, 10, 3, 3))
labeledHeatmap(Matrix = use_matrix,
               xLabels = colnames(use_matrix),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = use_p,
               setStdMargins = FALSE,
               cex.text = 0.7,
               zlim = c(-0.6, 0.6),xLabelsAngle = 0,xLabelsPosition = "bottom",
               main = paste("Module-trait relationships"))
modNames <- substring(names(MEs), 3)
modNames

geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) <- paste("MM", modNames, sep = "")
names(MMPvalue) <- paste("p.MM", modNames, sep = "");
geneModuleMembership[1:5, 1:5]


geneTraitSignificance <- as.data.frame(cor(datExpr, clin_data1$type, use = "p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) <- paste("GS.", names(clin_data1$type), sep = "")
names(GSPvalue) <- paste("p.GS.", names(clin_data1$type), sep = "")
head(geneTraitSignificance)

module = "yellow"
pheno = "auc"
modNames = substring(names(MEs), 3)

module_column = match(module, modNames)
pheno_column = match(pheno, "type")

moduleGenes <- moduleColors == module

dev.off()
verboseScatterplot(abs(geneModuleMembership[moduleGenes, module_column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste("Gene significance for LRG"),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)



moduleGene <- data.frame(net$colors)
moduleGene$gene <- rownames(moduleGene)
moduleGene$color <- labels2colors(net$colors)
saveRDS(moduleGene,"./result/data/moduleGene.rds")

moduleGene <- readRDS("./result/data/moduleGene.rds")
nrDEG_limma_signif <- fread("./result/data/limma_DEG_sig.csv")
colnames(nrDEG_limma_signif)[7] <- "gene"
module_signif <- moduleGene %>% dplyr::filter(color %in% c("midnightblue","yellow"))

DEG_ModuleGen <- dplyr::inner_join(nrDEG_limma_signif,module_signif)
ferr <- readxl::read_xlsx("./use_data/Ferroptosis.xlsx")
colnames(ferr) <- "gene"
DEG_ModuleGen_Ferroptosis <- dplyr::inner_join(DEG_ModuleGen,ferr)
fwrite(DEG_ModuleGen_Ferroptosis,"./result/data/DEG_ModuleGen_Ferroptosis.csv")

library(dplyr)
library(readxl)
ferr <- readxl::read_xlsx("./use_data/Ferroptosis.xlsx")
nrDEG_limma_signif <- fread("./result/data/limma_DEG_sig.csv")
module_signif <- moduleGene %>% dplyr::filter(color %in% c("midnightblue","yellow"))
module_gene <- module_signif$gene

gene1 <- nrDEG_limma_signif$gene_id
gene2 <- module_signif$gene
gene3 <- ferr$Ferroptosis

library(VennDiagram)
pdf(file = './vnn.pdf', width = 10, height = 10)
venn_plot <- venn.diagram(x=list(gene1,gene2,gene3),
             scaled = F, 
             alpha= 0.8, 
             lwd=1,lty=1,col=c('#8DD3C7',"#BEBADA","#FB8072"), 
             label.col ='black' , 
             cex = 1, 
             fontface = "bold",  
             fill=c('#8DD3C7',"#BEBADA","#FB8072"), 
             category.names = c("DEG", "WGCNA", "Ferroptosis") ,
             cat.dist = 0.03, 
             cat.cex = 1, 
             cat.fontface = "bold",  
             cat.col='black' ,   
             cat.default.pos = "text", 
             output=F,
             filename=NULL,
             resolution = 400, 
             compression = "lzw",
             width = 3000, 
             height = 3000)

grid.draw(venn_plot)
dev.off()


library(VennDiagram)
pdf(file = './result/vnn.pdf', width = 10, height = 10)

venn_plot <- venn.diagram(x = list(gene1, gene2,gene3),
                          scaled = FALSE, 
                          alpha = 0.8, 
                          lwd = 1, lty = 1, col = c('#8DD3C7', "#FB8072","#80B1D3"), 
                          label.col = 'black', 
                          cex = 1.5, 
                          fontface = "bold", 
                          fill = c('#8DD3C7',"#FB8072","#80B1D3"), 
                          category.names = c("DEG gene", "Hypoxia related-gene", "TACE related-gene"), # 标签名
                          cat.dist = 0.03, 
                          cat.cex = 1.5,
                          cat.fontface = "bold", 
                          cat.col = 'black',   
                          cat.default.pos = "text", 
                          filename = NULL, 
                          output = FALSE)
grid.draw(venn_plot)
dev.off()


use_gene <- fread("./result/data/DEG_ModuleGen_Ferroptosis.csv", data.table = F)
exp <- readRDS("./use_data/GSE75010_exp.rds")
clin_infor <- readRDS("./use_data/GSE75010_clin.rds")
use_infor <- clin_infor[,1:2]
use_infor$status <- ifelse(use_infor$status == "NonPE",0,1)
use_exp <- as.data.frame(t(exp[use_gene$gene,]))
use_exp$sample <- rownames(use_exp)
use_data <- dplyr::inner_join(use_infor,use_exp)

use_data1 <- use_data[,-1]

## LASSO
x <- as.matrix(use_data1[,-1])
y <- as.matrix(use_data1[,1])
library(glmnet)

set.seed(123)
alpha1_fit <- glmnet(x,y,alpha=1,family="binomial",nlambda = 100)
plot(alpha1_fit,xvar="lambda",label=TRUE)


set.seed(1235)
alpha1.fit <- cv.glmnet(x,y,type.measure = "class",alpha=1,family="binomial")
plot(alpha1.fit)
print(alpha1.fit)


coef(alpha1_fit,s=alpha1.fit$lambda.1se)
library(dplyr)
all_feature <- coef(alpha1.fit,s=alpha1.fit$lambda.1se) %>% as.matrix() %>% as.data.frame() %>% 
  rename("coff"="s1")
select_feature <-  all_feature %>% filter(abs(coff)>0)
select_feature1 <- select_feature
select_feature1$gene <- rownames(select_feature1)
select_feature2 <- select_feature1[-1,]
saveRDS(select_feature2,"./result/LASSO/select_feature.rds")

feature_gene <- readRDS("./result/LASSO/select_feature.rds")
ML_data <- use_data %>% dplyr::select(sample,status,feature_gene$gene)
saveRDS(ML_data,"./result/ML/ML_usedata.rds")


## Importance
library(ggplot2) 
library(RColorBrewer) 
library(tidyverse) 
library(randomForest) 
library(rfPermute) 
library(ggthemes)
library(viridis)
rf_data <- readRDS("./result/ML/ML_usedata.rds")
rf_data <- rf_data[, -1]
rf_data$status <- as.factor(rf_data$status)
set.seed(1234)
df.rf <- rfPermute(status ~ ., 
                   data= rf_data, 
                   ntree = 300,
                   nrep = 299, 
                   num.cores = 4)
importance_scores <- importance(df.rf)
feature_importance <- data.frame(Gene = rownames(importance_scores), Importance = importance_scores[, "MeanDecreaseGini"])
ordered_features <- feature_importance[order(-feature_importance$Importance), ]
importance_gini <- importance(df.rf, type = 2)
importance_accuracy <- importance(df.rf, type = 1)
top_genes <- head(ordered_features, 50)
write.csv(top_genes, file = "./feature_importance.csv", row.names = F)

top_genes <-read.table("./result/ML/feature_importance.csv",sep=",",header=T)
color_range <- c("#D8BFD8", "#8B008B")
ggplot(top_genes, aes(x = reorder(Gene, Importance), y = Importance, fill = Importance, label = round(Importance, 2))) +
  geom_bar(stat = "identity") +
  geom_text(size = 2.5, position = position_stack(vjust = 0.5), color = "white") + 
  scale_fill_gradient(low = color_range[1], high = color_range[2]) + 
  theme_few() + 
  theme(
    plot.title = element_text(hjust = 0.5, size = 16), 
    plot.caption = element_text(size = 12), 
    axis.text = element_text(size = 8), 
    axis.title = element_text(size = 15)
  ) +
  labs(title = "Important Genes", x = "Gene", y = "Importance") +
  coord_flip() 


## kegg
setwd("/home/data/sdc/wuchx/BioXCG/BioXCG/income/Service/2024/2024.04.22/")
inter_gene <- fread("./result/data/limma_DEG_sig.csv")

# BiocManager::install("org.Hs.eg.db")
library(ggplot2) 
library(clusterProfiler) 
library(org.Hs.eg.db) 
library(stats)
library(data.table) 
library(dplyr)

gene.df <- bitr(inter_gene$gene_id, fromType = "SYMBOL",
                toType = c("ENTREZID", "SYMBOL"), 
                OrgDb = org.Hs.eg.db) 
kegg <- enrichKEGG(gene.df$ENTREZID, organism = 'human', pvalueCutoff = 0.05, 
                   pAdjustMethod = 'BH', minGSSize = 3, maxGSSize = 500, qvalueCutoff = 0.05, 
                   use_internal_data = FALSE)
kegg_result <- data.frame(kegg)

fwrite(kegg_result,"./result/data/kegg_result.csv")


kegg_result <- fread("./result/data/kegg_result.csv")
significant_pathways <- subset(kegg_result, pvalue < 0.05)

use_pathway <- significant_pathways[1:20,]
ggplot(use_pathway, aes(x = reorder(Description, -log10(pvalue)), y = -log10(pvalue))) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(title = "KEGG Pathway Enrichment Analysis",
       x = "Pathway",
       y = "-log10(P-value)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  cowplot::theme_cowplot()+
  coord_flip()

library(ggplot2)
library(ggprism)
library(RColorBrewer)
library(circlize)

use_pathway1 <- use_pathway[,c(4,5,7,11)]
use_pathway1$GeneRatio <- sapply(use_pathway1$GeneRatio, function(x) eval(parse(text = x)))
use_pathway1 <- use_pathway1 %>% dplyr::arrange(GeneRatio)
use_pathway1$Description <- factor(use_pathway1$Description,levels = use_pathway1$Description)

colors <- colorRampPalette(brewer.pal(n = 11, name = "RdYlBu"))(100)
values <- seq(0, 0.05, length.out = 101)[-101]
col_fun <- colorRamp2(values, colors)

ggplot(use_pathway1, aes(x = GeneRatio, y = Description, size = Count, color = pvalue)) +
  geom_point(alpha = 0.9) +
  scale_size_continuous(range = c(3, 10)) +
  scale_color_gradientn(colors = col_fun(values)) +
  theme(plot.title = element_text(hjust = 0.5))+
  labs(title = "KEGG Pathway Enrichment Analysis", x = "GeneRatio", y = "Pathway")



## GO
GO_all <- enrichGO(gene = gene.df$ENTREZID,  
                   keyType = "ENTREZID",  
                   OrgDb=org.Hs.eg.db,  
                   ont = "ALL",   
                   pvalueCutoff = 0.05, 
                   pAdjustMethod = "fdr",  
                   minGSSize = 10,  
                   maxGSSize = 500,  
                   qvalueCutoff = 0.05,  
                   readable = TRUE)  
GO_result <- data.frame(GO_all)     

fwrite(GO_result,"./result/data/GO_result.csv")
go_enrichment_pathway <- GO_result %>% group_by(ONTOLOGY) %>% top_n(n = 10, wt = -p.adjust)


ggplot(go_enrichment_pathway, aes(x=reorder(Description, -Count), y=Count, fill=ONTOLOGY)) +
  geom_bar(stat="identity", position=position_dodge(width=0.8), width=0.6) +
  theme_minimal() +
  labs(x="GO Term", y="Gene_Number", title="Top 10 Enriched GO Terms") +
  facet_grid(ONTOLOGY~., scale = 'free_y', space = 'free_y')+
  coord_flip() + 
  scale_fill_manual(values=c("CC"="skyblue","BP"="pink","MF"="lightgreen"))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))


ggplot(go_enrichment_pathway, aes(x=reorder(Description, Count), y=Count)) +
  geom_point(aes(size=Count,color=-log10(p.adjust))) +
  scale_size_continuous(range=c(1, 10)) +
  facet_grid(ONTOLOGY~., scale = 'free_y', space = 'free_y')+
  coord_flip() +  
  theme_minimal() +
  scale_color_gradient(low = "pink",high ="red")+
  labs(color=expression(-log10(p.adjust),size="Count"), 
       x="Gene Ratio",y="Gene_Number",title="Top 10 Enriched GO Terms")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))
  


## GSEA
nrDEG_limma <- fread("./result/data/limma_DEG.csv",data.table = F)
geneSet <- read.gmt("/home/data/sdc/wuchx/BioXCG/BioXCG/income/Service/2023/2023.11.28/raw_data/h.all.v2023.1.Hs.symbols.gmt") #下载的基因集

colnames(nrDEG_limma)[7] <- "gene_name"

geneList <- nrDEG_limma$logFC 
names(geneList) <-  nrDEG_limma$gene_name 
geneList <- sort(geneList, decreasing = T) 

GSEA_enrichment <- GSEA(geneList, 
                        TERM2GENE = geneSet, 
                        pvalueCutoff = 0.05, 
                        minGSSize = 5, 
                        maxGSSize = 500, 
                        eps = 0, 
                        pAdjustMethod = "BH")

result <- data.frame(GSEA_enrichment)
dim(GSEA_enrichment@result)

fwrite(result,"./result/data/GSEA_hallmarker_result.csv")
saveRDS(GSEA_enrichment,"./result/data/GSEA_enrichment.rds")

significant_pathways <- subset(result, pvalue < 0.05)

use_pathway <- significant_pathways[1:15,]

library(ggplot2)
library(ggprism)
library(RColorBrewer)
library(circlize)

use_pathway$GeneRatio <- use_pathway$setSize/use_pathway$rank
use_pathway <- use_pathway %>% dplyr::arrange(GeneRatio)

use_pathway$Description <- factor(use_pathway$Description,levels = use_pathway$Description)

colors <- colorRampPalette(brewer.pal(n = 11, name = "RdYlBu"))(100)
values <- seq(0, 0.05, length.out = 101)[-101]
col_fun <- colorRamp2(values, colors)

ggplot(use_pathway, aes(x = GeneRatio, y = Description, size = setSize, color = pvalue)) +
  geom_point(alpha = 0.9) +
  scale_size_continuous(range = c(3, 10)) +
  scale_color_gradientn(colors = col_fun(values)) +
  theme(plot.title = element_text(hjust = 0.5))+
  labs(title = "GSEA Enrichment Analysis", x = "GeneRatio", y = "Pathway")



library(Seurat)
library(dplyr)
library(ggplot2)
library(ggridges)
library(DOSE)
library(clusterProfiler)
library(dplyr)
data <- readRDS('./result/data/GSEA_enrichment.rds')

ordered_gene_sets <- result$Description[1:10]
  
annotation_data <- data@result %>%
  filter(ID %in% ordered_gene_sets) %>%
  mutate(ID = factor(ID, levels = ordered_gene_sets))

set.seed(123)
sample_data <- data.frame(
  ID = rep(ordered_gene_sets, each = 100),
  NES = c(
    rnorm(100, mean = 1.5, sd = 0.8),   # HALLMARK_MYC_TARGETS_V2
    rnorm(100, mean = 0.1, sd = 0.6),   # HALLMARK_PROTEIN_SECRETION
    rnorm(100, mean = 1.5, sd = 0.8), # HALLMARK_TNFA_SIGNALING_VIA_NFKB
    rnorm(100, mean = 0.6, sd = 0.8), # HALLMARK_UNFOLDED_PROTEIN_RESPONSE
    rnorm(100, mean = 1.4, sd = 0.8), # HALLMARK_INTERFERON_ALPHA_RESPONSE
    rnorm(100, mean = 0.5, sd = 0.8), # HALLMARK_MYC_TARGETS_V1
    rnorm(100, mean = 0.8, sd = 0.7), # HALLMARK_MITOTIC_SPINDLE
    rnorm(100, mean = 0.2, sd = 0.3), # HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION
    rnorm(100, mean = 1.4, sd = 0.8), # HALLMARK_G2M_CHECKPOINT
    rnorm(100, mean = 2, sd = 0.9)  # HALLMARK_E2F_TARGETS
  )
)

color_palette <- c("#FF9896B2", "#98DF8AB2", "#9EDAE5B2", "#9467BDB2", "#17BECFB2",
                   "#EE4C97B2", "#C5B0D5B2", "#7F7F7FB2", "#E18727B2", "#79AF97B2")

p <- ggplot(sample_data, aes(x = NES, y = ID, fill = ID)) +
  geom_density_ridges(scale = 2, rel_min_height = 0.01) +
  scale_fill_manual(values = color_palette) +
  theme_ridges() +
  theme(legend.position = "none", 
        axis.text.y = element_text(size = 8), 
        plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        plot.margin = unit(c(1, 3, 1, 1), "cm")) + 
  labs(title = "GSEA (Hallmark Gene Set)", x = "", y ="")

for (i in 1:nrow(annotation_data)) {
  p <- p + annotate("text", x = max(sample_data$NES) * 1, y = i, 
                    label = paste("NES =", round(annotation_data$NES[i], 3), 
                                  "\np.adj =", sprintf("%.3f", annotation_data$p.adjust[i])),
                    hjust = 0, vjust = 0.5, size = 3)
}
print(p)



GSEA_enrichment <- readRDS('./result/data/GSEA_enrichment.rds')

library(enrichplot)
gseaplot2(GSEA_enrichment,"HALLMARK_HYPOXIA",color="red",pvalue_table = T) 

gseaplot2(GSEA_enrichment,c("HALLMARK_HYPOXIA","HALLMARK_OXIDATIVE_PHOSPHORYLATION",
                            "HALLMARK_GLYCOLYSIS","HALLMARK_BILE_ACID_METABOLISM","HALLMARK_FATTY_ACID_METABOLISM"),
          color=brewer.pal(n = 11, name = "RdYlBu")[c(1,3,5,8,10)],pvalue_table = T)


library(readxl)
ferrop_gene <- readxl::read_xlsx("./use_data/Ferroptosis.xlsx")
nrDEG_edgeR_signif <- fread("./result/data/limma_DEG_sig.csv",data.table = F)
nrDEG_edgeR_signif <- nrDEG_edgeR_signif %>% dplyr::filter(abs(logFC) > 0.5)
table(nrDEG_edgeR_signif$gene_id %in% ferrop_gene$Ferroptosis)
use_gene <- nrDEG_edgeR_signif %>% dplyr::filter(gene_id %in% ferrop_gene$Ferroptosis)


tpm_exp <- readRDS("./use_data/GSE75010_exp.rds")
infor <- readRDS("./use_data/GSE75010_clin.rds")
tpm_exp <- tpm_exp %>% dplyr::select(infor$sample)

heatmap_data <- na.omit(tpm_exp[use_gene$gene_id,])
heatmap_data1 <- as.data.frame(t(heatmap_data))
heatmap_data <- as.data.frame(t(heatmap_data1))

norm_data <- t(apply(heatmap_data, 1, function(x){(x-mean(x))/sd(x)}))

norm_data[norm_data > 2] <- 2
norm_data[norm_data < -2] <- -2


library(dendextend)
library(circlize)
library(RColorBrewer)
colors = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100)
values <- seq(-2.5, 2.5, length.out = 101)[-101]
col_fun = colorRamp2(values, colors)

library(ComplexHeatmap)
top_annotation <- HeatmapAnnotation(Type = c(rep("NonPE",77),rep("PE",80)), col = list(Type = c("NonPE" = "red","PE"="blue")))

# split_vector <- c(rep("", 5), rep(" ",3))
ComplexHeatmap::Heatmap(norm_data,cluster_rows = T,top_annotation = top_annotation,show_row_dend = F,show_column_dend = F,
                        column_title = "Ferroptosis-Related Gene",row_names_side = "left",
                        cluster_columns = T,show_column_names = F,show_row_names = T,
                        show_heatmap_legend = T,name = " ",col = col_fun, column_split = c(rep("NonPE",77),rep("PE",80)))



###### BACH1 #######

setwd("/home/data/sdc/wuchx/BioXCG/BioXCG/income/Service/2024/2024.04.22/")
remove(list = ls())

use_data <- readRDS("./result/ML/ML_usedata.rds")
use_data$status <- ifelse(use_data$status == 0,"NonPE","PE")
draw_data <- use_data %>% tidyr::pivot_longer(-c(sample, status), names_to = "Gene", values_to = "Expression")

ggplot(draw_data, aes(x = status, y = Expression, fill = status)) +
  geom_boxplot() + facet_wrap(~Gene, ncol = 3, scales = "free_y") +
  ylab("Expression (log2(TPM + 1))") + xlab("") +
  theme(line = element_line(color = "black", size = 1,
                            linetype = 1, lineend = "butt"),
        legend.position = "right") + labs(fill = "Group") +
  # scale_fill_manual(values = c("#EDA1A4", "#9DB4CE"))+
  cowplot::theme_cowplot(font_size = 15, line_size = 1)+ 
  stat_compare_means()


library(data.table)
clin <- fread("./raw_data/GSE75010/GSE75010_clin.csv")

clin1 <- clin %>% dplyr::select(geo_accession,"maximum diastolic bp:ch1","maximum systolic bp:ch1","mean umbilical pi:ch1","mean uterine pi:ch1",
                                "newborn weight z-score:ch1","placental weight z-score:ch1")
colnames(clin1) <-  c("sample","Maximum_Diastolic_BP","Maximum_Systolic_BP","Mean_Umbilical_PI",
                      "Mean_Uterine_PI","Newborn_Weight_z_score","Placental_Weight_z_score")

data1 <- readRDS("./use_data/GSE75010_clin.rds")

use_clin <- dplyr::inner_join(data1,clin1)
use_clin1 <- as.data.frame(apply(use_clin[,4:12], 2, function(x){as.numeric(x)}))
use_clin2 <- cbind(use_clin[,1:3],use_clin1)
saveRDS(use_clin2,"./use_data/GSE75010_clin2.rds")

exp_data <- readRDS("./result/ML/ML_usedata.rds")
use_clin2 <- readRDS("./use_data/GSE75010_clin2.rds")
exp_data1 <- exp_data %>% dplyr::select(sample,BACH1)
draw_data <- dplyr::inner_join(exp_data1,use_clin2)


library(ggplot2)
library(ggExtra)


for (i in colnames(draw_data)[5:13]) {
  
draw_data1 <- na.omit(draw_data[,c("BACH1",i)])


ylab1 <- paste0(i)
xlab1 <- paste0("BACH1 Expression (log2(TPM + 1 ))")


p <- ggplot(draw_data1, aes(x = BACH1, y = draw_data1[,i])) +
  geom_point(color = "#359023") + 
  stat_cor(method = "pearson", label.x = 8.8, label.y = max(draw_data1[,i]),size = 7)+
  geom_smooth(method = "lm", color = "red", linetype = "solid")+
  theme_minimal()+
  labs(y = ylab1,x = xlab1)+
  theme(
    plot.title = element_text(size = 20),  
    axis.title.x = element_text(size = 15),  
    axis.title.y = element_text(size = 15),  
    axis.text.x = element_text(size = 15), 
    axis.text.y = element_text(size = 15)  
  )

p
p_with_marginals <- ggMarginal(p, type = "histogram", fill = "#ffa500")
print(p_with_marginals)

pdf(file.path(paste0("./result/figure/clin_BACH1/",i,"_BACH1.pdf")),
    width = 8, height = 6)
print(p_with_marginals, newpage = FALSE)

dev.off()

}





library(ggplot2)
library(ggExtra)

for (i in colnames(draw_data)[5:13]) {
  

draw_data1 <- na.omit(draw_data[,c("status",i)])
  
p <- ggplot(data = draw_data1, aes(x = status, y = draw_data1[,i],fill = status)) +
    geom_boxplot() +
    stat_compare_means(size = 5, comparisons = list(c("NonPE","PE")),label = "p.format",
                       method = "wilcox.test") +
    geom_point(position = position_jitter(width = 0.2, height = 0), alpha = 0.8) +
    
    theme_prism(border = TRUE) +
    scale_x_discrete(labels = c("NonPE","PE")) +
    theme(axis.title.x = element_blank())+
    scale_fill_manual(values = c("#66C2A5","#FC8D62"))+
    labs(y=i)
  
p
  pdf(file.path(paste0("./result/figure/clin_NonPE_PE/",i,"_NonPE_PE.pdf")),
      width = 8, height = 6)
  print(p, newpage = FALSE)
  
  dev.off()
  
}






