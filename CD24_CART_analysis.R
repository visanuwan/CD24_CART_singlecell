library(Seurat)
library(SeuratDisk)
library(MAST)
library(dplyr)
library(ggplot2)

#################################################
###### Load and QC data
#################################################

expression_matrix <- ReadMtx(cells = "Normal_1_barcodes.tsv.gz", features = "Normal_1_features.tsv.gz", mtx = "Normal_1_matrix.mtx.gz")
Normal_1 <- CreateSeuratObject(counts = expression_matrix, project = "Normal_1", min.features = 200)
Normal_1[["percent.mt"]] <- PercentageFeatureSet(wt_normal, pattern = "^mt-")
Normal_1 <- subset(Normal_1, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & nCount_RNA < 70000 & percent.mt < 10)

expression_matrix <- ReadMtx(cells = "Normal_2_barcodes.tsv.gz", features = "Normal_2_features.tsv.gz", mtx = "Normal_2_matrix.mtx.gz")
Normal_2 <- CreateSeuratObject(counts = expression_matrix, project = "Normal_2", min.features = 200)
Normal_2[["percent.mt"]] <- PercentageFeatureSet(wt_normal, pattern = "^mt-")
Normal_2 <- subset(Normal_2, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & nCount_RNA < 70000 & percent.mt < 10)

expression_matrix <- ReadMtx(cells = "Normal_3_barcodes.tsv.gz", features = "Normal_3_features.tsv.gz", mtx = "Normal_3_matrix.mtx.gz")
Normal_3 <- CreateSeuratObject(counts = expression_matrix, project = "Normal_3", min.features = 200)
Normal_3[["percent.mt"]] <- PercentageFeatureSet(wt_normal, pattern = "^mt-")
Normal_3 <- subset(Normal_3, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & nCount_RNA < 70000 & percent.mt < 10)

expression_matrix <- ReadMtx(cells = "PBS_1_barcodes.tsv.gz", features = "PBS_1_features.tsv.gz", mtx = "PBS_1_matrix.mtx.gz")
PBS_1 <- CreateSeuratObject(counts = expression_matrix, project = "PBS_1", min.features = 200)
PBS_1[["percent.mt"]] <- PercentageFeatureSet(wt_normal, pattern = "^mt-")
PBS_1 <- subset(PBS_1, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & nCount_RNA < 70000 & percent.mt < 10)

expression_matrix <- ReadMtx(cells = "PBS_2_barcodes.tsv.gz", features = "PBS_2_features.tsv.gz", mtx = "PBS_2_matrix.mtx.gz")
PBS_2 <- CreateSeuratObject(counts = expression_matrix, project = "PBS_2", min.features = 200)
PBS_2[["percent.mt"]] <- PercentageFeatureSet(wt_normal, pattern = "^mt-")
PBS_2 <- subset(PBS_2, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & nCount_RNA < 70000 & percent.mt < 10)

expression_matrix <- ReadMtx(cells = "PBS_3_barcodes.tsv.gz", features = "PBS_3_features.tsv.gz", mtx = "PBS_3_matrix.mtx.gz")
PBS_3 <- CreateSeuratObject(counts = expression_matrix, project = "PBS_3", min.features = 200)
PBS_3[["percent.mt"]] <- PercentageFeatureSet(wt_normal, pattern = "^mt-")
PBS_3 <- subset(PBS_3, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & nCount_RNA < 70000 & percent.mt < 10)

expression_matrix <- ReadMtx(cells = "MOCK_1_barcodes.tsv.gz", features = "MOCK_1_features.tsv.gz", mtx = "MOCK_1_matrix.mtx.gz")
MOCK_1 <- CreateSeuratObject(counts = expression_matrix, project = "MOCK_1", min.features = 200)
MOCK_1[["percent.mt"]] <- PercentageFeatureSet(wt_normal, pattern = "^mt-")
MOCK_1 <- subset(MOCK_1, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & nCount_RNA < 70000 & percent.mt < 10)

expression_matrix <- ReadMtx(cells = "MOCK_2_barcodes.tsv.gz", features = "MOCK_2_features.tsv.gz", mtx = "MOCK_2_matrix.mtx.gz")
MOCK_2 <- CreateSeuratObject(counts = expression_matrix, project = "MOCK_2", min.features = 200)
MOCK_2[["percent.mt"]] <- PercentageFeatureSet(wt_normal, pattern = "^mt-")
MOCK_2 <- subset(MOCK_2, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & nCount_RNA < 70000 & percent.mt < 10)

expression_matrix <- ReadMtx(cells = "MOCK_3_barcodes.tsv.gz", features = "MOCK_3_features.tsv.gz", mtx = "MOCK_3_matrix.mtx.gz")
MOCK_3 <- CreateSeuratObject(counts = expression_matrix, project = "MOCK_3", min.features = 200)
MOCK_3[["percent.mt"]] <- PercentageFeatureSet(wt_normal, pattern = "^mt-")
MOCK_3 <- subset(MOCK_3, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & nCount_RNA < 70000 & percent.mt < 10)

expression_matrix <- ReadMtx(cells = "CD24_1_barcodes.tsv.gz", features = "CD24_1_features.tsv.gz", mtx = "CD24_1_matrix.mtx.gz")
CD24_1 <- CreateSeuratObject(counts = expression_matrix, project = "CD24_1", min.features = 200)
CD24_1[["percent.mt"]] <- PercentageFeatureSet(wt_normal, pattern = "^mt-")
CD24_1 <- subset(CD24_1, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & nCount_RNA < 70000 & percent.mt < 10)

expression_matrix <- ReadMtx(cells = "CD24_2_barcodes.tsv.gz", features = "CD24_2_features.tsv.gz", mtx = "CD24_2_matrix.mtx.gz")
CD24_2 <- CreateSeuratObject(counts = expression_matrix, project = "CD24_2", min.features = 200)
CD24_2[["percent.mt"]] <- PercentageFeatureSet(wt_normal, pattern = "^mt-")
CD24_2 <- subset(CD24_2, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & nCount_RNA < 70000 & percent.mt < 10)

expression_matrix <- ReadMtx(cells = "CD24_3_barcodes.tsv.gz", features = "CD24_3_features.tsv.gz", mtx = "CD24_3_matrix.mtx.gz")
CD24_3 <- CreateSeuratObject(counts = expression_matrix, project = "CD24_3", min.features = 200)
CD24_3[["percent.mt"]] <- PercentageFeatureSet(wt_normal, pattern = "^mt-")
CD24_3 <- subset(CD24_3, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & nCount_RNA < 70000 & percent.mt < 10)



#################################################
###### Normalize and prepare data for integration
#################################################

cart.list <- list(Normal_1, Normal_2, Normal_3, PBS_1, PBS_2, PBS_3, MOCK_1, MOCK_2, MOCK_3, CD24_1, CD24_2, CD24_3)
cart.list <- lapply(X = cart.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = cart.list)
cart.list <- lapply(X = cart.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

#################################################
###### Perform integration
#################################################

cart.anchors <- FindIntegrationAnchors(object.list = cart.list, anchor.features = features, reduction = "rpca")
cart.combined <- IntegrateData(anchorset = cart.anchors)
DefaultAssay(cart.combined) <- "integrated"

#################################################
###### Perform the standard workflow for visualization and clustering
#################################################

cart.combined <- ScaleData(cart.combined, verbose = FALSE)
cart.combined <- RunPCA(cart.combined, npcs = 50, verbose = FALSE)
cart.combined <- RunUMAP(cart.combined, reduction = "pca", dims = 1:30, min.dist = 0.4, n.neighbors = 20)
cart.combined <- FindNeighbors(cart.combined, reduction = "pca", dims = 1:30)
cart.combined <- FindClusters(cart.combined, resolution = 0.2)

cart.combined@meta.data <- cart.combined@meta.data %>% mutate(combined.ident = orig.ident)
cart.combined@meta.data[cart.combined@meta.data$orig.ident == 'Normal_1',]$combined.ident <- 'Normal'
cart.combined@meta.data[cart.combined@meta.data$orig.ident == 'Normal_2',]$combined.ident <- 'Normal'
cart.combined@meta.data[cart.combined@meta.data$orig.ident == 'Normal_3',]$combined.ident <- 'Normal'
cart.combined@meta.data[cart.combined@meta.data$orig.ident == 'PBS_1',]$combined.ident <- 'PBS'
cart.combined@meta.data[cart.combined@meta.data$orig.ident == 'PBS_2',]$combined.ident <- 'PBS'
cart.combined@meta.data[cart.combined@meta.data$orig.ident == 'PBS_3',]$combined.ident <- 'PBS'
cart.combined@meta.data[cart.combined@meta.data$orig.ident == 'MOCK_1',]$combined.ident <- 'MOCK'
cart.combined@meta.data[cart.combined@meta.data$orig.ident == 'MOCK_2',]$combined.ident <- 'MOCK'
cart.combined@meta.data[cart.combined@meta.data$orig.ident == 'MOCK_3',]$combined.ident <- 'Mock'
cart.combined@meta.data[cart.combined@meta.data$orig.ident == 'CD24_1',]$combined.ident <- 'CD24'
cart.combined@meta.data[cart.combined@meta.data$orig.ident == 'CD24_2',]$combined.ident <- 'CD24'
cart.combined@meta.data[cart.combined@meta.data$orig.ident == 'CD24_3',]$combined.ident <- 'CD24'

cart.combined$combined.ident <- factor(x = cart.combined$combined.ident, levels = c("Normal", "PBS", "MOCK", "CD24"))

DimPlot(cart.combined, reduction = "umap",label = T, label.size = 5, repel = F)
DimPlot(cart.combined, reduction = "umap",label = T, label.size = 5, repel = F, split.by = 'orig.ident') + NoLegend()
DimPlot(cart.combined, reduction = "umap",label = T, label.size = 5, repel = F, group.by = 'orig.ident') + NoLegend()

#################################################
###### Cell type annotation
#################################################

DefaultAssay(cart.combined) <- "RNA"

# Neutrophils
FeaturePlot(cart.combined, features = c("Ly6g", "Mmp8", "Cxcr2"))

# Myelocytes
FeaturePlot(cart.combined, features = c("Fcnb", "Ltf", "Lcn2"))

# Myeloblasts
FeaturePlot(cart.combined, features = c("Elane", "Mpo", "Ctsg", "Ms4a3"))

# Monoblasts
FeaturePlot(cart.combined, features = c("F13a1", "Irf8", "Ly86"))

# Basophils
FeaturePlot(cart.combined, features = c("Ccl3", "Fcer1a", "Ms4a2", "Hdc", "Cyp4f18", "Csf1", "Ccl4", "Hgf"))

# DCs
FeaturePlot(cart.combined, features = c("Itgax", "Siglech"))

# Monocytes
FeaturePlot(cart.combined, features = c("S100a4", "Pld4", "Csf1r"))

# Macrophages
FeaturePlot(cart.combined, features = c("Adgre1"))

# NK and T cells
FeaturePlot(cart.combined, features = c("Gzma", "Klra4", "Cd3d", "Cd3e"))

# Pro-B cells
FeaturePlot(cart.combined, features = c("Vpreb3", "Akap12"))

# Pre-B cells
FeaturePlot(cart.combined, features = c("Cd74", "H2Ab1"))

# Erythroblasts
FeaturePlot(cart.combined, features = c("Car2", "Gypa", "Prdx2", "Alas2", "Slc4a1", "Hemgn", "Spta1", "Slc25a37"))



# Assign cell type to clusters
cart.combined@active.ident <- factor(cart.combined@active.ident,
                                     levels = c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16'))
cart.combined$cell_types <- Idents(cart.combined)
cart.combined@meta.data <- cart.combined@meta.data %>% mutate(name_cell_types = cell_types)
levels(cart.combined$name_cell_types)[levels(cart.combined$name_cell_types) == '1'] <- 'Neutrophils-I'
levels(cart.combined$name_cell_types)[levels(cart.combined$name_cell_types) == '2'] <- 'Neutrophils-II'
levels(cart.combined$name_cell_types)[levels(cart.combined$name_cell_types) == '3'] <- 'Myelocytes-I'
levels(cart.combined$name_cell_types)[levels(cart.combined$name_cell_types) == '4'] <- 'Myelocytes-II'
levels(cart.combined$name_cell_types)[levels(cart.combined$name_cell_types) == '5'] <- 'Myeloblasts'
levels(cart.combined$name_cell_types)[levels(cart.combined$name_cell_types) == '6'] <- 'Monoblasts-I'
levels(cart.combined$name_cell_types)[levels(cart.combined$name_cell_types) == '7'] <- 'Monoblasts-II'
levels(cart.combined$name_cell_types)[levels(cart.combined$name_cell_types) == '8'] <- 'Basophils'
levels(cart.combined$name_cell_types)[levels(cart.combined$name_cell_types) == '9'] <- 'DCs'
levels(cart.combined$name_cell_types)[levels(cart.combined$name_cell_types) == '10'] <- 'Monocytes'
levels(cart.combined$name_cell_types)[levels(cart.combined$name_cell_types) == '11'] <- 'Macrophages'
levels(cart.combined$name_cell_types)[levels(cart.combined$name_cell_types) == '12'] <- 'NK/T'
levels(cart.combined$name_cell_types)[levels(cart.combined$name_cell_types) == '13'] <- 'Pro-B'
levels(cart.combined$name_cell_types)[levels(cart.combined$name_cell_types) == '14'] <- 'Pre-B-I'
levels(cart.combined$name_cell_types)[levels(cart.combined$name_cell_types) == '15'] <- 'Pre-B-II'
levels(cart.combined$name_cell_types)[levels(cart.combined$name_cell_types) == '16'] <- 'Erythroblasts'


#################################################
###### Cell-type specific markers
#################################################

DefaultAssay(cart.combined) <- "RNA"
all.markers <- FindAllMarkers(object = cart.combined, test.use = "MAST")
top20.markers <- all.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

#################################################
###### Subcluster analysis
#################################################

# Macrophages
Idents(cart.combined) <- "name_cell_types"
cart.combined.mcp <- subset(cart.combined, idents = "Macrophages")
cart.combined.mcp <- RunPCA(cart.combined.mcp, npcs = 50, verbose = TRUE)
cart.combined.mcp <- FindNeighbors(cart.combined.mcp, reduction = "pca", dims = 1:50)
cart.combined.mcp <- FindClusters(cart.combined.mcp, graph.name = 'integrated_snn',resolution = 0.85, verbose = TRUE)
cart.combined.mcp <- RunUMAP(cart.combined.mcp, reduction = "pca", dims = 1:50, min.dist = 0.6, n.neighbors = 6)
DimPlot(cart.combined.mcp, reduction = "umap",label = T, label.size = 4, repel = F, split.by = 'combined.ident')
