# Load libraries
library('Seurat')
library('ggplot2')
library('dplyr')

#########   Input Data    ###########################################

Sample1.data <- Read10X(data.dir = "Z:/Blackmore_Lab/scRNAseq_Analyses/1_StndID/LU/1_SourceFiles/LU2")
Sample1 <- CreateSeuratObject(counts = Sample1.data, project = "Sample1", 
                              min.cells = 3, min.features = 1000)
Sample1 
# Initially 2710
VlnPlot(Sample1, features = c("nFeature_RNA"))
Sample1 <- subset(Sample1, subset = nFeature_RNA > 1000 & nFeature_RNA < 6000)
Sample1 
# Now 2654

Sample2.data <- Read10X(data.dir = "Z:/Blackmore_Lab/scRNAseq_Analyses/1_StndID/LU/1_SourceFiles/LU3")
Sample2 <- CreateSeuratObject(counts = Sample2.data, project = "Sample2", 
                              min.cells = 3, min.features = 1000)
Sample2 
# Initially 2770
VlnPlot(Sample2, features = c("nFeature_RNA"))
Sample2 <- subset(Sample2, subset = nFeature_RNA > 1000 & nFeature_RNA < 8000)
Sample2 
# Now 2608


Sample3.data <- Read10X(data.dir = "Z:/Blackmore_Lab/scRNAseq_Analyses/1_StndID/LU/1_SourceFiles/LU4")
Sample3 <- CreateSeuratObject(counts = Sample3.data, project = "Sample3", 
                              min.cells = 3, min.features = 1000)
Sample3
#Initially 3579
VlnPlot(LUN1, features = c("nFeature_RNA"))
Sample3 <- subset(Sample3, subset = nFeature_RNA > 1000 & nFeature_RNA < 7000)
Sample3
#Now 3518

# Save RDS files
saveRDS(Sample1, "Z:/Blackmore_Lab/scRNAseq_Analyses/1_StndID/LU/Sample1_1000_6000.rds")
saveRDS(Sample2, "Z:/Blackmore_Lab/scRNAseq_Analyses/1_StndID/LU/Sample2_1000_8000.rds")
saveRDS(Sample3, "Z:/Blackmore_Lab/scRNAseq_Analyses/1_StndID/LU/Sample3_1000_7000.rds")


#To open the individual sample files if needed
Sample1 <- readRDS("Z:/Blackmore_Lab/scRNAseq_Analyses/1_StndID/LU/Sample1_1000_6000.rds")
Sample2 <- readRDS("Z:/Blackmore_Lab/scRNAseq_Analyses/1_StndID/LU/Sample2_1000_8000.rds")
Sample3 <- readRDS("Z:/Blackmore_Lab/scRNAseq_Analyses/1_StndID/LU/Sample3_1000_7000.rds")

############# Merge Data #####################
#Step 2 - merge the three datasets

LUN_merge.list <- list (Sample1, Sample2, Sample3)

LUN_merge.list <- lapply(X=LUN_merge.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = LUN_merge.list)

LUN_merge.anchors <- FindIntegrationAnchors(object.list = LUN_merge.list,
                                            anchor.features = features)

LUN <- IntegrateData(anchorset = LUN_merge.anchors)

DefaultAssay(LUN) <- "integrated"

LUN <- ScaleData(LUN, verbose = FALSE)

LUN <- RunPCA(LUN, npcs = 50, verbose = FALSE)

LUN <- RunUMAP(LUN, reduction = "pca", dims = 1:30)

LUN <- FindNeighbors(LUN, reduction = "pca", dims = 1:30)

LUN <- FindClusters(LUN, resolution = 0.5)

DefaultAssay(LUN) <- "RNA"

LUN

#Indicates 8780


# Some visualizations for initial assessment of clustering and sample bias


DimPlot(LUN, group.by = 'seurat_clusters', split.by = 'orig.ident',
        label = TRUE, label.size = 5) + NoLegend() # 1500 x 500

saveRDS(LuN, "Z:/Blackmore_Lab/scRNAseq_Analyses/1_StndID/LuN/LuN.rds")
LUN <- readRDS("Z:/Blackmore_Lab/scRNAseq_Analyses/1_StndID/LUN/LUN.rds")

###############################################################################
#                             Save and Restart Point
###############################################################################



#Remove the clearly stressed, dying, or artificially low cell count clusters.

FeaturePlot(LUN, c("Creb5"), label = TRUE, 
            label.size = 5, repel = TRUE, order = TRUE) #500 x 500

FeaturePlot(LUN, c("Atf3"), label = TRUE, 
            label.size = 5, repel = TRUE, order = TRUE) #500 x 500

FeaturePlot(LUN, c("nFeature_RNA"), label = TRUE, 
            label.size = 5, repel = TRUE, order = TRUE) #500 x 500

VlnPlot(LUN, "nFeature_RNA", pt.size = 0) + NoLegend() +
  theme(axis.text.x = element_text(angle = 0, hjust = .5)) +
  theme(text = element_text(size = 20)) #1500 x 500

####Based on the above QC metrics, cluster 6, 7, 20, and 13 are removed

Idents(LUN) <- "seurat_clusters"

rm(LUN1)

LUN1 <- subset(x = LUN, idents = c("6", "7", "20", "13"), invert = TRUE)

DimPlot(LUN1, group.by = 'seurat_clusters', split.by = 'orig.ident', label = TRUE, label.size = 3) + NoLegend()

FeaturePlot(LUN1, c("Daam2"), 
            label = TRUE, label.size = 3,
            repel = FALSE, order = TRUE, ) #500 x 1000

############  Renormalize the new dataset that lacks the dying and low-count clusters cells   ################

DefaultAssay(LUN1) <- "RNA"

LUN1 <- NormalizeData(LUN1, normalization.method = "LogNormalize", scale.factor = 10000)

LUN1 <- FindVariableFeatures(LUN1, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(LUN1)

LUN1 <- ScaleData(LUN1, features = all.genes)

LUN1 <- RunPCA(LUN1, features = VariableFeatures(object = LUN1), npcs = 50)

LUN1 <- FindNeighbors(LUN1, dims = 1:30)
LUN1 <- FindClusters(LUN1, resolution = 0.5)

LUN1 <- RunUMAP(LUN1, dims = 1:30)

LUN1 #7748

# This gets me to a pretty good looking Umap.

# Various visualizations to decide whether clustering is appropriate

FeaturePlot(LUN1, c("Gad2"), 
            label = TRUE, label.size = 3,
            repel = FALSE, order = TRUE, ) #500 x 1000

DimPlot(LUN1, group.by = 'seurat_clusters',  label = TRUE, label.size = 3) + NoLegend()


saveRDS(LUN, "Z:/Blackmore_Lab/scRNAseq_Analyses/1_StndID/LUN/LUN.rds")
saveRDS(LUN1, "Z:/Blackmore_Lab/scRNAseq_Analyses/1_StndID/LU/LUN1.rds")
LUN1 <- readRDS("Z:/Blackmore_Lab/scRNAseq_Analyses/1_StndID/LUN/LUN1.rds")


############################################################

##Launch point for Figure 1, Basic Characterization###

##Also Launch for Figure 2, CST identification ####

############################################################


### Identify CST and also establish markers to exclude non-
#CST cortical cells and determine their number

# Make violin plots of Satb2, Camk2a, and Crym and Fezf2 

Idents(LUN1) <- "seurat_clusters"

VlnPlot(LUN1, "Satb2", pt.size = 0) + NoLegend() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  theme(text = element_text(size = 15),
        axis.text = element_text(size = 15)) #Sized for display

VlnPlot(LUN1, "Camk2a", pt.size = 0) + NoLegend() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  theme(text = element_text(size = 15),
        axis.text = element_text(size = 15)) #Sized for display

VlnPlot(LUN1, "Crym", pt.size = 0) + NoLegend() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  theme(text = element_text(size = 15),
        axis.text = element_text(size = 15)) #Sized for display

VlnPlot(LUN1, "Fezf2", pt.size = 0) + NoLegend() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  theme(text = element_text(size = 15),
        axis.text = element_text(size = 15)) #Sized for display

VlnPlot(LUN1, "Slc30a3", pt.size = 0) + NoLegend() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  theme(text = element_text(size = 15),
        axis.text = element_text(size = 15)) #Sized for display



VlnPlot(LUN1, "nFeature_RNA", pt.size = 0) + NoLegend() +
  theme(axis.text.x = element_text(angle = 0, hjust = .5)) +
  theme(text = element_text(size = 20)) #1500 x 500




### Create lists of differentially expressed genes that mark putative
# CST clusters and any non-layer V cortical neurons

Idents(LUN1) <- "seurat_clusters"

LUN1_CST_Markers <- FindMarkers(LUN1,
                                ident.1 = c("0","2"))

write.csv(x = LUN1_CST_Markers, 
          file = "Z:/Blackmore_Lab/scRNAseq_Analyses/1_StndID/LUN/Fig2_CST_ID/LUN1_CST_Markers.csv",
          quote = FALSE)




LUN1_13_Markers <- FindMarkers(LUN1,
                               ident.1 = c("13"))
write.csv(x = LUN1_13_Markers, 
          file = "Z:/Blackmore_Lab/scRNAseq_Analyses/1_StndID/LUN/Fig2_CST_ID/LUN1_13_Markers.csv",
          quote = FALSE)



LUN1_CSTvs13_Markers <- FindMarkers(LUN1,
                                    ident.1 = c("13"),
                                    ident.2 = c("0", "1"))

write.csv(x = LUN1_CSTvs13_Markers, 
          file = "Z:/Blackmore_Lab/scRNAseq_Analyses/1_StndID/LUN/Fig2_CST_ID/LUN1_CSTvs13_Markers.csv",
          quote = FALSE)


#Check these markers on a dot plot: they are pan-cortex, layer V, then intracortical

DotPlot(LUN1, features = c("Mef2c", "Grin2a", "Camk2a", "Ptk2b", "Satb2", 
                           "Crym", "Fezf2", "Bcl6", "Slc30a3"), dot.scale = 4, dot.min = 0) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 0, hjust = .5),
        axis.text = element_text(size = 10)) #Sized for display


# Based on these plots, cluster 0 and 1 are likely CST and 13 is non-CST cortical



#Clean up the environment

rm(LUN1_13_Markers)
rm(LUN1_CST_Markers)
rm(LUN1_CSTvs13_Markers)


saveRDS(LUN1, "Z:/Blackmore_Lab/scRNAseq_Analyses/1_StndID/LUN/LUN1.rds")
LUN1 <- readRDS("Z:/Blackmore_Lab/scRNAseq_Analyses/1_StndID/LUN/LUN1.rds")

###############################################################################
#                             Save and Restart Point
###############################################################################






################################################################################################

#Launch Point for Figure 2, CST Identification

################################################################################################


# Now subcluster the subcortical neurons, renormalize, and assign identity and order

Idents(LUN1) <- "seurat_clusters"


# Note that cluster 13 has high levels of slc30a3, which is an intracortical marker

LUN_Final <- subset(LUN1, idents = c("13"), invert = TRUE)


#Read a dotplot with standard markers

LUN_Final_DotPlot1 <- read.csv("Z:/Blackmore_Lab/scRNAseq_Analyses/1_StndID/LUN/Fig3_Final/LUN_Final_DotPlot1.csv",
                               sep = ",")

Idents(LUN_Final)  <- "seurat_clusters"

DefaultAssay(LUN_Final) <- "RNA"

DotPlot(LUN_Final, features = LUN_Final_DotPlot1$Marker, dot.scale = 9, dot.min = 0) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 0, hjust = .5),
        axis.text = element_text(size = 10)) #Sized for display


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

Idents(LUN_Final) <- "seurat_clusters"
levels(LUN_Final)


#check that dotplots still work
DefaultAssay(LUN_Final) <- "RNA"

DotPlot(LUN_Final, features = LUN_Final_DotPlot1$Marker, dot.scale = 9, dot.min = 0) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 0, hjust = .5),
        axis.text = element_text(size = 10)) #Sized for display

FeaturePlot(LUN_Final, c("Gad2"), 
            label = TRUE, label.size = 3,
            repel = FALSE, order = TRUE, ) #500 x 1000

DimPlot(LUN_Final, group.by = 'seurat_clusters',  label = TRUE, label.size = 3) + NoLegend()

#Yes, still looks good


#Order manual.clusters.1 in a rostral-caudal sequence

LUN_Final_Order1 <- read.csv("Z:/Blackmore_Lab/scRNAseq_Analyses/1_StndID/LUN/Fig3_Final/LUN_Final_Order1.csv",
                             sep = ",")

LUN_Final@meta.data$manual.clusters.1 <- factor(
  x = LUN_Final@meta.data$manual.clusters.1,
  levels = LUN_Final_Order1$Order)

Idents(LUN_Final) <- "manual.clusters.1"
levels(LUN_Final)


#Check identities with Dotplots

LUN_Final_DotPlot2 <- read.csv("Z:/Blackmore_Lab/scRNAseq_Analyses/1_StndID/LUN/Fig3_Final/LUN_Final_DotPlot2.csv",
                               sep = ",")


DefaultAssay(LUN_Final) <- "RNA"


DotPlot(LUN_Final, features = LUN_Final_DotPlot1$Marker, dot.scale = 3, dot.min = 0) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 10)) #Sized for display


DotPlot(LUN_Final, features = LUN_Final_DotPlot2$Marker, dot.scale = 4, dot.min = 0) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 10)) #Sized for display



Idents(LUN_Final) <- "manual.clusters.1"
DotPlot(LUN_Final, features = LUN_Final_DotPlot2$Marker, dot.scale = 4, dot.min = 0) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 10)) #Sized for display












#Read a dotplot with standard markers

LUN_Final_DotPlot1 <- read.csv("Z:/Blackmore_Lab/scRNAseq_Analyses/1_StndID/LUN/Fig3_SC/LUN_Final_DotPlot1.csv",
                               sep = ",")

Idents(LUN_Final)  <- "seurat_clusters"
DefaultAssay(LUN_Final) <- "RNA"

DotPlot(LUN_Final, features = LUN_Final_DotPlot1$Marker, dot.scale = 9, dot.min = 0) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 0, hjust = .5),
        axis.text = element_text(size = 10)) #Sized for display

FeaturePlot(LUN_Final, c("Crym"), 
            label = TRUE, label.size = 3,
            repel = FALSE, order = FALSE, ) #500 x 1000

LUN_Final_11 <- FindMarkers(LUN_Final,
                            ident.1 = c("11"))

write.csv(x = LUN_Final_11, 
          file = "Z:/Blackmore_Lab/scRNAseq_Analyses/1_StndID/LUN/Fig3_SC/LUN_Final_11.csv",
          quote = FALSE)



#Map seurat_clusters to manual.cluster.1

LUN_Final_Map1 <- read.csv("Z:/Blackmore_Lab/scRNAseq_Analyses/1_StndID/LUN/Fig3_Final/LUN_Final_Map1.csv",
                           sep = ",")

LUN_Final@meta.data$manual.clusters.1 <- plyr::mapvalues(
  x = LUN_Final@meta.data$seurat_clusters,
  from = LUN_Final_Map1$seurat_clusters,
  to = LUN_Final_Map1$manual.cluster.1)

#Order Seurat Clusters in a rostral-caudal sequence

LUN_Final@meta.data$seurat_clusters <- factor(
  x = LUN_Final@meta.data$seurat_clusters,
  levels = LUN_Final_Map1$seurat_clusters)

Idents(LUN_Final) <- "seurat_clusters"

DotPlot(LUN_Final, features = LUN_Final_DotPlot1$Marker, dot.scale = 4, dot.min = 0) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 0, hjust = .5),
        axis.text = element_text(size = 10)) #Sized for display

#Order manual.clusters.1 in a rostral-caudal sequence

LUN_Final_Order1 <- read.csv("Z:/Blackmore_Lab/scRNAseq_Analyses/1_StndID/LUN/Fig3_Final/LUN_Final_Order1.csv",
                             sep = ",")

LUN_Final@meta.data$manual.clusters.1 <- factor(
  x = LUN_Final@meta.data$manual.clusters.1,
  levels = LUN_Final_Order1$Order)



Idents(LUN_Final) <- "manual.clusters.1"

DotPlot(LUN_Final, features = LUN_Final_DotPlot1$Marker, dot.scale = 4, dot.min = 0) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 10)) #Sized for display

LUN_Final_DotPlot2 <- read.csv("Z:/Blackmore_Lab/scRNAseq_Analyses/1_StndID/LUN/Fig3_Final/LUN_Final_DotPlot2.csv",
                               sep = ",")


DotPlot(LUN_Final, features = LUN_Final_DotPlot2$Marker, dot.scale = 4, dot.min = 0) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 10)) #Sized for display


saveRDS(LUN_Final, "Z:/Blackmore_Lab/scRNAseq_Analyses/1_StndID/LuN/LuN_Final.rds")
LUN_Final <- readRDS("Z:/Blackmore_Lab/scRNAseq_Analyses/1_StndID/LUN/LUN_Final.rds")

###############################################################################
#                             Save and Restart Point
###############################################################################

#Export a file that lists the average expression of all transcripts for each cluster

DefaultAssay(LUN_Final) <- "RNA"

avgexp <- AverageExpression(LUN_Final, assay = 'RNA', 
                            features = rownames(LUN_Final@assays$RNA@counts),
                            group.by = c('manual.clusters.1'))
avgexp$RNA

write.csv(file = "Z:/Blackmore_Lab/scRNAseq_Analyses/1_StndID/LuN/Fig6_Lists/AvgExpression.csv", x = avgexp)





#Identify marker genes for manual.clusters.1 and export to csv

LU1_SC_all <- FindAllMarkers(LU1_SC, only.pos = TRUE)

write.csv(x = LU1_SC_all, 
          file = "Z:/Blackmore_Lab/scRNAseq_Analyses/1_StndID/LU/Fig3_SC/LU1_SC_all.csv",
          quote = FALSE)

LU1_SC_DotPlot2 <- read.csv("Z:/Blackmore_Lab/scRNAseq_Analyses/1_StndID/LU/Fig3_SC/LU1_SC_Dotplot2.csv",
                            sep = ",")

DotPlot(LU1_SC, features = LU1_SC_DotPlot2$Marker, dot.scale = 4, dot.min = 0) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 0, hjust = .5),
        axis.text = element_text(size = 10)) #Sized for display

DotPlot(LUN_Final, features = "Lhx3", dot.scale = 4, dot.min = 0) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 0, hjust = .5),
        axis.text = element_text(size = 10)) #Sized for display


DimPlot(LUN_Final, group.by = 'manual.clusters.1', label = TRUE, 
        label.size = 15, repel = TRUE, 
        pt.size = 4) + 
  theme(text = element_text(size = 70),
        axis.text = element_text(size = 70))  #3000x3000


#### Determine the number of nuclei in each cluster

Idents(LUN_Final) <- "manual.clusters.1"
cell.num <- table(Idents(LUN_Final))
cell.num

write.csv(x = cell.num, 
          file = "Z:/Blackmore_Lab/scRNAseq_Analyses/1_StndID/LUN/Fig4_HB/Counts.csv",
          quote = FALSE)



LUN_Final

Idents(LUN1) <- "seurat_clusters"
cell.num <- table(Idents(LUN1))
cell.num

LUN1

### Find DEGs in CST versus all others

Idents(LUN_Final) <- "manual.clusters.1"
levels(LUN_Final)
LUN_CST_Markers <- FindMarkers(LUN_Final,
                               ident.1 = c("CST"),
                               ident.2 = c("RN", "HB-1-Daam2",
                                           "HB-2-Pard3b",
                                           "HB-3-Col19a1",
                                           "HB-4-Nox4",
                                           "HB-5-Other"))


write.csv(x = LUN_CST_Markers, 
          file = "Z:/Blackmore_Lab/scRNAseq_Analyses/1_StndID/LUN/Fig2_CST_ID/LUN_CST_Markers.csv",
          quote = FALSE)



######Final save of all seurat objects

saveRDS(Sample1, "Z:/Blackmore_Lab/scRNAseq_Analyses/1_StndID/LUn/Sample1_1000_6000.rds")
saveRDS(Sample2, "Z:/Blackmore_Lab/scRNAseq_Analyses/1_StndID/LUN/Sample2_1000_8000.rds")
saveRDS(Sample3, "Z:/Blackmore_Lab/scRNAseq_Analyses/1_StndID/LUN/Sample3_1000_7000.rds")

saveRDS(LUN, "Z:/Blackmore_Lab/scRNAseq_Analyses/1_StndID/LUN/LUN.rds")
saveRDS(LUN1, "Z:/Blackmore_Lab/scRNAseq_Analyses/1_StndID/LUN/LUN1.rds")
saveRDS (LUN_Final, "Z:/Blackmore_Lab/scRNAseq_Analyses/1_StndID/LUN/LUN_Final.rds")






