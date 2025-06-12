
remove(list = ls())
library(dplyr)
library(data.table)
library(Matrix)
library(Seurat)

gene1 <- fread("./raw_data/GSE173193_RAW/GSM5261699_P1_features.tsv.gz",data.table = F,header = F)
barcode1 <- fread("./raw_data/GSE173193_RAW/GSM5261699_P1_barcodes.tsv.gz",data.table = F,header = F)
data1 <- readMM("./raw_data/GSE173193_RAW/GSM5261699_P1_matrix.mtx.gz")
colnames(data1) <- barcode1$V1 
rownames(data1) <- gene1$V2
data1 <- data.frame(data1)

gene2 <- fread("./raw_data/GSE173193_RAW/GSM5261700_P2_features.tsv.gz",data.table = F,header = F)
barcode2 <- fread("./raw_data/GSE173193_RAW/GSM5261700_P2_barcodes.tsv.gz",data.table = F,header = F)
data2 <- readMM("./raw_data/GSE173193_RAW/GSM5261700_P2_matrix.mtx.gz")
colnames(data2) <- barcode2$V1 
rownames(data2) <- gene2$V2
data2 <- data.frame(data2)

seurat1 <- CreateSeuratObject(counts = data1, project = "P1", min.cells = 3, min.features = 200)
seurat2 <- CreateSeuratObject(counts = data2, project = "P2", min.cells = 3, min.features = 200)

combined <- merge(seurat1, y = list(seurat2), 
                  add.cell.ids = c("P1", "P2"), 
                  project = "GSE173193")

pbmc <- JoinLayers(combined)

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt < 10)
dim(pbmc)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(pbmc)

# PCA
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc)) 
ElbowPlot(pbmc)
pbmc <- FindNeighbors(pbmc, dims = 1:15)
pbmc <- FindClusters(pbmc, resolution = 0.05)

pbmc <- RunUMAP(pbmc, dims = 1:15) 
DimPlot(pbmc, reduction = "umap", pt.size = 0.5,label = F) 

FeaturePlot(pbmc, features = c("BACH1"),order = T,pt.size = 0.5)

# pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# use_marker <- pbmc.markers %>% dplyr::filter(avg_log2FC > 2)
# table(use_marker$cluster)
# fwrite(use_marker,"./result/use_marker.csv")

DimPlot(pbmc, reduction = "umap", pt.size = 0.5,label = T) 
new.cluster.ids <- c("Villous cytotrophoblast (VCT)", "Neutrophils", "Syncytiotrophoblast (SCT)", 
                     "Macrophages", "Extravillous trophoblasts (EVT)", "T cells / NK cells", "Erythroid progenitors",
                     "Monocyte", "Erythroid progenitors")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
pbmc@meta.data$cell_type <- Idents(pbmc)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5)

new_order <- c(
  "Villous cytotrophoblast (VCT)",
  "Extravillous trophoblasts (EVT)",
  "Syncytiotrophoblast (SCT)",
  "Erythroid progenitors",
  "Neutrophils",
  "Macrophages",
  "Monocyte",
  "T cells / NK cells"
)

pbmc$cell_type <- factor(Idents(pbmc), levels = new_order)
Idents(pbmc) <- pbmc$cell_type 

saveRDS(pbmc, file = "./result/单细胞/final.rds")


