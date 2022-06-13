# This script creates Figure 1
# The input for this script is LUN1

# Load libraries
library('Seurat')
library('ggplot2')
library('plyr')
library('dplyr')

#Input data






#1c - a large labeled overview of the clusters


DimPlot(LUN1, group.by = 'seurat_clusters', label = FALSE, 
        label.size = 40, repel = TRUE, 
        pt.size = 4) + 
  theme(text = element_text(size = 100),
        axis.text = element_text(size = 70)) +
  NoLegend()  + NoAxes() #3000x3000


#1d - showing consistency across samples

Idents(LUN1) <- "seurat_clusters"

DimPlot(LUN1, group.by = 'seurat_clusters', label = FALSE, 
        label.size = 15, repel = TRUE, split.by = "orig.ident", 
        pt.size = 2) + 
  theme(text = element_text(size = 40),
        axis.text = element_text(size = 40)) +
  NoLegend() + NoAxes() #3000x1000

#Create a dotplot with markers for off-target cells (based on Tuszynski contamination)

Fig1markers <- read.csv("Z:/Blackmore_Lab/scRNAseq_Analyses/1_StndID/LUN/Fig1_InitChar/Fig1markers.csv",
                        sep=",")

DefaultAssay(LUN1) <- "RNA"

DotPlot(LUN1, features = Fig1markers$Gene, dot.scale = 6, dot.min = 0) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1),
        axis.text = element_text(size = 20)) #1000 x 1500 then x2 in Adobe

## Create a series of violin plots for neurotransmitter types

VlnPlot(LUN1, "Slc17a7", pt.size = 0) + NoLegend() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  theme(text = element_text(size = 30),
        axis.text = element_text(size = 30)) #1000 x 600

VlnPlot(LUN1, "Slc17a6", pt.size = 0) + NoLegend() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  theme(text = element_text(size = 30),
        axis.text = element_text(size = 30)) #1000 x 600

VlnPlot(LUN1, "Gad2", pt.size = 0) + NoLegend() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  theme(text = element_text(size = 30),
        axis.text = element_text(size = 30)) #1000 x 600

VlnPlot(LUN1, "Tph2", pt.size = 0) + NoLegend() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  theme(text = element_text(size = 30),
        axis.text = element_text(size = 30)) #1000 x 600

VlnPlot(LUN1, "Slc6a2", pt.size = 0) + NoLegend() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  theme(text = element_text(size = 30),
        axis.text = element_text(size = 30)) #1000 x 600





# Make feature plots

FeaturePlot(LUN1, c("Slc17a7"), label = FALSE, pt.size = 1,
            label.size = 5, repel = TRUE, order = FALSE) +
  theme(text = element_text(size = 25),
        axis.text = element_text(size = 40)) +
  NoLegend() + NoAxes() #700 x 1000

FeaturePlot(LUN1, c("Slc17a6"), label = FALSE, pt.size = 1,
            label.size = 5, repel = TRUE, order = FALSE) +
  theme(text = element_text(size = 25),
        axis.text = element_text(size = 40)) +
  NoLegend() + NoAxes() #700 x 1000

FeaturePlot(LUN1, c("Gad2"), label = FALSE, pt.size = 1,
            label.size = 5, repel = TRUE, order = FALSE) +
  theme(text = element_text(size = 25),
        axis.text = element_text(size = 40)) +
  NoLegend() + NoAxes() #700 x 1000

FeaturePlot(LUN1, c("Tph2"), label = FALSE, pt.size = 1,
            label.size = 5, repel = TRUE, order = FALSE) +
  theme(text = element_text(size = 25),
        axis.text = element_text(size = 40)) +
  NoLegend() + NoAxes() #700 x 1000

FeaturePlot(LUN1, c("Slc6a2"), label = FALSE, pt.size = 1,
            label.size = 5, repel = TRUE, order = FALSE) +
  theme(text = element_text(size = 25),
        axis.text = element_text(size = 40)) +
  NoLegend() + NoAxes() #700 x 1000




