library(Seurat)
library(SeuratDisk)
library(patchwork)
library(dplyr)
library(MAST)

## treat_01
import_h5 <- Read10X_h5("C:/Users/sunfumou/Desktop/Emory CART/01d28.h5")
sobj <- CreateSeuratObject(counts = import_h5$`Gene Expression`, assay = 'RNA', project = 'CD24', min.features = 200)
sobj[['percent.mt']] <- PercentageFeatureSet(sobj, pattern = "^MT-")
treat_01 <- subset(sobj, subset = percent.mt < 10 & nCount_RNA < 70000)

treat_01$sample_id <- '01d28'
treat_01$state <- 'after-treatment'
treat_01$pair <- '01'
treat_01$comb_id <- 'treat_01'
treat_01$group <- 'long-relapse'
treat_01$comb_group <- 'long-relapse after-treatment'

## treat_19
import_h5 <- Read10X_h5("C:/Users/sunfumou/Desktop/Emory CART/19d28.h5")
sobj <- CreateSeuratObject(counts = import_h5$`Gene Expression`, assay = 'RNA', project = 'CD24', min.features = 200)
sobj[['percent.mt']] <- PercentageFeatureSet(sobj, pattern = "^MT-")
treat_19 <- subset(sobj, subset = percent.mt < 10 & nCount_RNA < 70000)

treat_19$sample_id <- '19d28'
treat_19$state <- 'after-treatment'
treat_19$pair <- '19'
treat_19$comb_id <- 'treat_19'
treat_19$group <- 'long-relapse'
treat_19$comb_group <- 'long-relapse after-treatment'

## treat_33
import_h5 <- Read10X_h5("C:/Users/sunfumou/Desktop/Emory CART/33d28.h5")
sobj <- CreateSeuratObject(counts = import_h5$`Gene Expression`, assay = 'RNA', project = 'CD24', min.features = 200)
sobj[['percent.mt']] <- PercentageFeatureSet(sobj, pattern = "^MT-")
treat_33 <- subset(sobj, subset = percent.mt < 10 & nCount_RNA < 70000)

treat_33$sample_id <- '33d28'
treat_33$state <- 'after-treatment'
treat_33$pair <- '33'
treat_33$comb_id <- 'treat_33'
treat_33$group <- 'long-relapse'
treat_33$comb_group <- 'long-relapse after-treatment'

## treat_16
import_h5 <- Read10X_h5("C:/Users/sunfumou/Desktop/Emory CART/16d28.h5")
sobj <- CreateSeuratObject(counts = import_h5$`Gene Expression`, assay = 'RNA', project = 'CD24', min.features = 200)
sobj[['percent.mt']] <- PercentageFeatureSet(sobj, pattern = "^MT-")
treat_16 <- subset(sobj, subset = percent.mt < 10 & nCount_RNA < 70000)

treat_16$sample_id <- '16d28'
treat_16$state <- 'after-treatment'
treat_16$pair <- '16'
treat_16$comb_id <- 'treat_16'
treat_16$group <- 'short-relapse'
treat_16$comb_group <- 'short-relapse after-treatment'

## treat_20
import_h5 <- Read10X_h5("C:/Users/sunfumou/Desktop/Emory CART/20d28.h5")
sobj <- CreateSeuratObject(counts = import_h5$`Gene Expression`, assay = 'RNA', project = 'CD24', min.features = 200)
sobj[['percent.mt']] <- PercentageFeatureSet(sobj, pattern = "^MT-")
treat_20 <- subset(sobj, subset = percent.mt < 10 & nCount_RNA < 70000)

treat_20$sample_id <- '20d28'
treat_20$state <- 'after-treatment'
treat_20$pair <- '20'
treat_20$comb_id <- 'treat_20'
treat_20$group <- 'short-relapse'
treat_20$comb_group <- 'short-relapse after-treatment'

## treat_27
import_h5 <- Read10X_h5("C:/Users/sunfumou/Desktop/Emory CART/27d28.h5")
sobj <- CreateSeuratObject(counts = import_h5$`Gene Expression`, assay = 'RNA', project = 'CD24', min.features = 200)
sobj[['percent.mt']] <- PercentageFeatureSet(sobj, pattern = "^MT-")
treat_27 <- subset(sobj, subset = percent.mt < 10 & nCount_RNA < 70000)

treat_27$sample_id <- '27d28'
treat_27$state <- 'after-treatment'
treat_27$pair <- '27'
treat_27$comb_id <- 'treat_27'
treat_27$group <- 'short-relapse'
treat_27$comb_group <- 'short-relapse after-treatment'

## pre_01
import_h5 <- Read10X_h5("C:/Users/sunfumou/Desktop/Emory CART/01pre.h5")
sobj <- CreateSeuratObject(counts = import_h5$`Gene Expression`, assay = 'RNA', project = 'CD24', min.features = 200)
sobj[['percent.mt']] <- PercentageFeatureSet(sobj, pattern = "^MT-")
pre_01 <- subset(sobj, subset = percent.mt < 10 & nCount_RNA < 70000)

pre_01$sample_id <- '01pre'
pre_01$state <- 'pre-treatment'
pre_01$pair <- '01'
pre_01$comb_id <- 'pre_01'
pre_01$group <- 'long-relapse'
pre_01$comb_group <- 'long-relapse pre-treatment'

## pre_16
import_h5 <- Read10X_h5("C:/Users/sunfumou/Desktop/Emory CART/16pre.h5")
sobj <- CreateSeuratObject(counts = import_h5$`Gene Expression`, assay = 'RNA', project = 'CD24', min.features = 200)
sobj[['percent.mt']] <- PercentageFeatureSet(sobj, pattern = "^MT-")
pre_16 <- subset(sobj, subset = percent.mt < 10 & nCount_RNA < 70000)

pre_16$sample_id <- '16pre'
pre_16$state <- 'pre-treatment'
pre_16$pair <- '16'
pre_16$comb_id <- 'pre_16'
pre_16$group <- 'short-relapse'
pre_16$comb_group <- 'short-relapse pre-treatment'

## pre_19
import_h5 <- Read10X_h5("C:/Users/sunfumou/Desktop/Emory CART/19pre.h5")
sobj <- CreateSeuratObject(counts = import_h5$`Gene Expression`, assay = 'RNA', project = 'CD24', min.features = 200)
sobj[['percent.mt']] <- PercentageFeatureSet(sobj, pattern = "^MT-")
pre_19 <- subset(sobj, subset = percent.mt < 10 & nCount_RNA < 70000)

pre_19$sample_id <- '19pre'
pre_19$state <- 'pre-treatment'
pre_19$pair <- '19'
pre_19$comb_id <- 'pre_19'
pre_19$group <- 'long-relapse'
pre_19$comb_group <- 'long-relapse pre-treatment'

## pre_32
import_h5 <- Read10X_h5("C:/Users/sunfumou/Desktop/Emory CART/32pre.h5")
sobj <- CreateSeuratObject(counts = import_h5$`Gene Expression`, assay = 'RNA', project = 'CD24', min.features = 200)
sobj[['percent.mt']] <- PercentageFeatureSet(sobj, pattern = "^MT-")
pre_32 <- subset(sobj, subset = percent.mt < 10 & nCount_RNA < 70000)

pre_32$sample_id <- '32pre'
pre_32$state <- 'pre-treatment'
pre_32$pair <- '32'
pre_32$comb_id <- 'pre_32'
pre_32$group <- 'short-relapse'
pre_32$comb_group <- 'short-relapse pre-treatment'

## pre_33
import_h5 <- Read10X_h5("C:/Users/sunfumou/Desktop/Emory CART/33pre.h5")
sobj <- CreateSeuratObject(counts = import_h5$`Gene Expression`, assay = 'RNA', project = 'CD24', min.features = 200)
sobj[['percent.mt']] <- PercentageFeatureSet(sobj, pattern = "^MT-")
pre_33 <- subset(sobj, subset = percent.mt < 10 & nCount_RNA < 70000)

pre_33$sample_id <- '33pre'
pre_33$state <- 'pre-treatment'
pre_33$pair <- '33'
pre_33$comb_id <- 'pre_33'
pre_33$group <- 'long-relapse'
pre_33$comb_group <- 'long-relapse pre-treatment'


####################
#### integrate data
####################

# list of Seurat object
cart.list <- list(pre_01, treat_01,
                  pre_16, treat_16,
                  pre_19, treat_19,
                  pre_32, treat_20,
                  pre_33, treat_27,
                  treat_33)

# normalize and identify variable features for each dataset independently
cart.list <- lapply(X = cart.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

## fasta integration using rpca
# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
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

cart.combined <- ScaleData(cart.combined, verbose = FALSE)
cart.combined <- RunPCA(cart.combined, npcs = 50, verbose = FALSE)
cart.combined <- FindNeighbors(cart.combined, reduction = "pca", dims = 1:30)
cart.combined <- FindClusters(cart.combined, resolution = 0.5)
cart.combined <- RunUMAP(cart.combined, reduction = "pca", dims = 1:30, min.dist = 0.58, n.neighbors = 20)

Idents(cart.combined) <- 'state'
Idents(cart.combined) <- factor(Idents(cart.combined), levels = c('pre-treatment',
                                                                  'after-treatment'))
cart.combined$state <- factor(x = cart.combined$state, levels = c('pre-treatment',
                                                                  'after-treatment'))
levels(cart.combined)
Idents(cart.combined) <- 'seurat_clusters'

DimPlot(cart.combined, reduction = "umap", label = T, label.size = 4,  repel = F)

DimPlot(cart.combined, reduction = "umap", label = T, label.size = 3,  repel = F,
               split.by = 'state', ncol = 2)




# SaveH5Seurat(cart.combined, filename = 'results/integrated_pre_treat')

SaveH5Seurat(cart.combined, filename = paste0(Project(object = cart.combined), ".h5seurat"),
             overwrite = FALSE,
             verbose = TRUE,)

# top20.markers
cart.combined <- LoadH5Seurat("C:/Users/sunfumou/Desktop/scRNASeq/Emory_pre_d28.h5seurat")
all.markers <- FindAllMarkers(object = cart.combined, test.use = 'MAST', min.pct = 0.1, only.pos = TRUE)
top20.markers <- all.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
write.csv(top20.markers, file = 'top20.markers.csv', row.names = FALSE)

rep_cell_counts <- group_by(cart.combined@meta.data, seurat_clusters, comb_id) %>% summarise(count = n())
write.table(rep_cell_counts,file = "Cell counts.csv",sep = ",",row.names = FALSE)

####################
#### select data
####################

sub.cart.combined <- subset(cart.combined, idents = c('3'))
DimPlot(sub.cart.combined, reduction = "umap", label = T, label.size = 3,  repel = F,
        split.by = 'state', ncol = 2)

####################
#### find DEG
####################
Idents(sub.cart.combined) <- 'state'

markers <- FindMarkers(sub.cart.combined, assay = 'RNA', slot = 'data',
                         ident.1 = 'after-treatment', ident.2 = 'pre-treatment',
                         test.use = 'MAST',
                         min.pct = 0.1, only.pos = FALSE, verbose = TRUE)
markers$gene <- rownames(markers_2)

write.csv(markers, file = 'M.csv', row.names = FALSE) ## write CSV



DefaultAssay(sub.cart.combined) <- 'RNA'
VlnPlot(sub.cart.combined, features = c("CD24"), split.by = 'state')


exprs <- data.frame(FetchData(sub.cart.combined, c("CD24")))
vln.markers <- exprs
write.csv(vln.markers, file = 'vln.csv')
write.csv(sub.cart.combined@meta.data, file = 'mata.csv')



####################
#### subset
####################

sub_cluster <- subset(cart.combined, idents = c('4','10','11','15','16'))
sub_cluster <- RunPCA(sub_cluster, npcs = 50, verbose = FALSE)
sub_cluster <- FindNeighbors(sub_cluster, reduction = "pca", dims = 1:50)
sub_cluster <- FindClusters(sub_cluster, graph.name = 'integrated_snn', resolution = 0.2) # resolution: number of clusters
sub_cluster <- RunUMAP(sub_cluster, reduction = "pca", dims = 1:50, min.dist = 0.01, n.neighbors = 10) # min.dist: cluster tightness, n.neighbors: structure of cluster
DimPlot(sub_cluster, reduction = "umap",label = T, label.size = 4, repel = F, split.by = 'state')

##Gene expression
DefaultAssay(cart.combined) <- 'RNA'
FeaturePlot(cart.combined, features = c('CD24'), split.by = 'state')
VlnPlot(cart.combined, features = c("CD24"),split.by = 'state')
vlnplot_data <- as.data.frame(vlnplot)


##Differemt Gene
DefaultAssay(sub_cluster) <- "RNA"
all.markers <- FindAllMarkers(object = sub_cluster, test.use = 'MAST')
top20.markers <- all.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
write.table(all.markers,file = "DE.csv",sep = ",",row.names = FALSE)


####################
#### combine clusters
####################

# change all ids to the same cluster id : ('1','4','25) -> ('1','1','1')

new.cluster.ids <- c('1','1','2','2','3','1','5','1','2','4','3','3','6','5','5','3',
                     '3','7','5','8','1','9','10','2','5','8','11','2')

names(new.cluster.ids) <- levels(cart.combined)
cart.combined <- RenameIdents(cart.combined, new.cluster.ids)

DimPlot(cart.combined, reduction = "umap", label = T, label.size = 3,  repel = F)

DimPlot(cart.combined, reduction = "umap", label = T, label.size = 3,  repel = F, split.by = 'state', ncol = 2)

