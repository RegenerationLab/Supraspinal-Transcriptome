# This script creates Figure 6, the creation of a Dimplot and Dotplot
# with marker genes for the final clusters.

# Input for this script is LUN1_Final

saveRDS(LUN_Final, "Z:/Blackmore_Lab/scRNAseq_Analyses/1_StndID/LuN/LuN_Final.rds")
LUN_Final <- readRDS("Z:/Blackmore_Lab/scRNAseq_Analyses/1_StndID/LUN/LUN_Final.rds")

######################################################################################
# This is a review of how clusters were assigned and ordered and the source file
# for the dotplot

#Map seurat_clusters to manual.cluster.1

LUN_Final_Map1 <- read.csv("Z:/Blackmore_Lab/scRNAseq_Analyses/1_StndID/LUN/Fig3_Final/LUN_Final_Map1.csv",
                           sep = ",")

Idents(LUN_Final) <- "seurat_clusters"
levels(LUN_Final)

LUN_Final@meta.data$manual.clusters.1 <- plyr::mapvalues(
  x = LUN_Final@meta.data$seurat_clusters,
  from = LUN_Final_Map1$seurat_clusters,
  to = LUN_Final_Map1$manual.cluster.1)

#Order Seurat Clusters in a rostral-caudal sequence

LUN_Final@meta.data$seurat_clusters <- factor(
  x = LUN_Final@meta.data$seurat_clusters,
  levels = LUN_Final_Map1$seurat_clusters)

#Order manual.clusters.1 in a rostral-caudal sequence

LUN_Final_Order1 <- read.csv("Z:/Blackmore_Lab/scRNAseq_Analyses/1_StndID/LUN/Fig3_Final/LUN_Final_Order1.csv",
                             sep = ",")

LUN_Final@meta.data$manual.clusters.1 <- factor(
  x = LUN_Final@meta.data$manual.clusters.1,
  levels = LUN_Final_Order1$Order)

#####################################################################################


#6a - a large labeled overview of the clusters


DimPlot(LUN_Final, group.by = 'manual.clusters.1', label = FALSE, 
        label.size = 40, repel = TRUE, 
        pt.size = 4) + 
  theme(text = element_text(size = 100),
        axis.text = element_text(size = 70)) +
  NoLegend()  + NoAxes() #2000x3000



#6b - Dotplot with marker genes

# Read csv with selected marker genes

LUN_Final_DotPlot2 <- read.csv("Z:/Blackmore_Lab/scRNAseq_Analyses/1_StndID/LUN/Fig3_Final/LUN_Final_DotPlot2.csv",
                               sep = ",")

DefaultAssay(LUN_Final) <- "RNA"
Idents(LUN_Final) <- "manual.clusters.1"

DotPlot(LUN_Final, features = LUN_Final_DotPlot2$Marker, dot.scale = 6, dot.min = 0) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 20)) #1000 x 1500 then x2 in Adobe