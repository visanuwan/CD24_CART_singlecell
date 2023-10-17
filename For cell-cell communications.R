library(Seurat)
library(ggplot2)
library(dplyr)
library(MAST)
suppressMessages({
  library(celltalker)
  library(Seurat)
  suppressWarnings(
    library(SeuratData)
  )
  library(tidyverse)
})

library(CellChat)
library(NMF)
library(tidyverse)
library(ggalluvial)
library(patchwork)
options(stringsAsFactors = FALSE)
dput(unique(combined_cluster$combined.ident))
Idents(combined_cluster)<- "orig.ident"
combined_cluster$combined.ident<- as.character(combined_cluster$combined.ident)
combined_cluster$name_cell_types<- as.character(combined_cluster$name_cell_types)
Idents(combined_cluster)

normal <- combined_cluster[,combined_cluster$combined.ident %in% c("normal")]
pbs <- combined_cluster[,combined_cluster$combined.ident %in% c("pbs")]
mock <- combined_cluster[,combined_cluster$combined.ident %in% c("mock")]
cd24<- combined_cluster[,combined_cluster$combined.ident %in% c("cd24")]

b<- combined_cluster[,combined_cluster$orig.ident %in% c("cd24_3")]
dput(unique(b$cell_types))


##cellchat
Idents(normal)
b1<-createCellChat(normal@assays$RNA@data, meta = normal@meta.data, group.by = "name_cell_types")
b2<-createCellChat(pbs@assays$RNA@data, meta = pbs@meta.data, group.by = "name_cell_types")
b3<-createCellChat(mock@assays$RNA@data, meta = mock@meta.data, group.by = "name_cell_types")
b4<-createCellChat(cd24@assays$RNA@data, meta = cd24@meta.data, group.by = "name_cell_types")
save(b1, b2,b3,b4, file ="1nor_2pbs_3moc_4cd14.rda")



##for normal
cellchat<-b1
cellchat@DB<- CellChatDB.mouse
cellchat<- subsetData(cellchat)
cellchat<- identifyOverExpressedGenes(cellchat)
cellchat<- identifyOverExpressedInteractions(cellchat)
cellchat<- computeCommunProb(cellchat, raw.use =TRUE, population.size = TRUE)
cellchat<-computeCommunProbPathway(cellchat)
cellchat<- aggregateNet(cellchat)
cellchat<- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

b1<- cellchat
saveRDS(b1,"normal.cellchat.rds")

##for pbs
cellchat<- b2
cellchat@DB<- CellChatDB.mouse
cellchat<- subsetData(cellchat)
cellchat<- identifyOverExpressedGenes(cellchat)
cellchat<- identifyOverExpressedInteractions(cellchat)
cellchat<- computeCommunProb(cellchat, raw.use =TRUE, population.size = TRUE)
cellchat<-computeCommunProbPathway(cellchat)
cellchat<- aggregateNet(cellchat)
cellchat<- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

b2<- cellchat
saveRDS(b2,"pbs.cellchat.rds")

##for mock
cellchat<- b3
cellchat@DB<- CellChatDB.mouse
cellchat<- subsetData(cellchat)
cellchat<- identifyOverExpressedGenes(cellchat)
cellchat<- identifyOverExpressedInteractions(cellchat)
cellchat<- computeCommunProb(cellchat, raw.use =TRUE, population.size = TRUE)
cellchat<-computeCommunProbPathway(cellchat)
cellchat<- aggregateNet(cellchat)
cellchat<- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

b3<- cellchat
saveRDS(b3,"mock.cellchat.rds")

##for cd24
cellchat<- b4
cellchat@DB<- CellChatDB.mouse
cellchat<- subsetData(cellchat)
cellchat<- identifyOverExpressedGenes(cellchat)
cellchat<- identifyOverExpressedInteractions(cellchat)
cellchat<- computeCommunProb(cellchat, raw.use =TRUE, population.size = TRUE)
cellchat<-computeCommunProbPathway(cellchat)
cellchat<- aggregateNet(cellchat)
cellchat<- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

b4<- cellchat
saveRDS(b4,"cd24.cellchat.rds")
##in one group
cellchat<- b132
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)

netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions") 
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

h1<- netVisual_heatmap(cellchat)
h2<- netVisual_heatmap(cellchat, measure = "weight")
h1+h2
# Access all the signaling pathways showing significant communications
pathways.show.all <- cellchat@netP$pathways
# check the order of cell identity to set suitable vertex.receiver

levels(cellchat@idents)
vertex.receiver = seq(1,4)
for (i in 1:length(pathways.show.all)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
  ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
}



##differential analysis

gg1<- compareInteractions(cellchat, show.legend = F, group = c(1,2,3,4), measure = "count")
gg2<- compareInteractions(cellchat, show.legend = F, group = c(1,2,3,4), measure = "weight")
p<- gg1+gg2
p
dput(unique(mock$name_cell_types))

cco.list<- list(moc=b3, cd=b4)
cellchat<- mergeCellChat(cco.list, add.names= names(cco.list), cell.prefix= TRUE)



par(mfrow=c(1,2))
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

par(mfrow=c(1,1))
h1<- netVisual_heatmap(cellchat)
h2<- netVisual_heatmap(cellchat, measure = "weight")
h1+h2
