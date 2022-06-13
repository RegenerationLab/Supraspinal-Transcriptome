# This file is to load and visualize dot plot of interesting genes

# Load libraries
library('Seurat')
library('ggplot2')
library('plyr')
library('dplyr')

# Input file LUN_Final

LUN_Final <- readRDS("Z:/Blackmore_Lab/scRNAseq_Analyses/1_StndID/LU/LUN_Final.rds")


#Export a file that lists the average expression of all transcripts for each cluster

DefaultAssay(LUN_Final) <- "RNA"

avgexp <- AverageExpression(LUN_Final, assay = 'RNA', 
                            features = rownames(LUN_Final@assays$RNA@counts),
                            group.by = c('manual.clusters.1'))
avgexp$RNA

write.csv(file = "Z:/Blackmore_Lab/scRNAseq_Analyses/1_StndID/LuN/Fig9_Lists/AvgExpression.csv", x = avgexp)





#Input list of TFs

TFs <- read.csv("Z:/Blackmore_Lab/scRNAseq_Analyses/1_StndID/LUN/Fig9_Lists/Tfs_v1.csv",
                sep = ",")


DotPlot(LUN_Final, features = TFs$TFs, dot.scale = 5, dot.min = 0.01) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 10)) #600 x 1200 for v1



#Input list of receptors

Receptors <- read.csv("Z:/Blackmore_Lab/scRNAseq_Analyses/1_StndID/LUN/Fig9_Lists/Receptors_v1.csv",
                      sep = ",")


Idents(LUN_Final)  <- "manual.clusters.1"

DefaultAssay(LUN_Final) <- "RNA"

DotPlot(LUN_Final, features = Receptors$Receptor, dot.scale = 5, dot.min = 0.01) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 10)) #600 x 1200


# Adhesive receptors

Receptors3 <- read.csv("Z:/Blackmore_Lab/scRNAseq_Analyses/1_StndID/LUN/Fig6_Lists/Receptors3.csv",
                       sep = ",")


Idents(LUN_Final)  <- "manual.clusters.1"

DefaultAssay(LUN_Final) <- "RNA"

DotPlot(LUN_Final, features = Receptors3$Category, dot.scale = 4, dot.min = 0) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 10)) #Sized for display

# Regeneration Associated Genes

RAGs <- read.csv("Z:/Blackmore_Lab/scRNAseq_Analyses/1_StndID/LUN/Fig6_Lists/RAGs.csv",
                 sep = ",")


Idents(LUN_Final)  <- "manual.clusters.1"

DefaultAssay(LUN_Final) <- "RNA"

DotPlot(LUN_Final, features = RAGs$Gene, dot.scale = 7, dot.min = 0) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 10)) #Sized for display

#################################################################################



VlnPlot(LUN_Final, "Hk2", pt.size = 0) + NoLegend() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(text = element_text(size = 20), axis.text = element_text(size = 10)) #1000 x 300

VlnPlot(LUN_Final, "Klf7", pt.size = 0) + NoLegend() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(text = element_text(size = 20), axis.text = element_text(size = 10)) #1000 x 300

VlnPlot(LUN_Final, "", pt.size = 0) + NoLegend() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(text = element_text(size = 20), axis.text = element_text(size = 10)) #1000 x 300


VlnPlot(LUN_Final, "Plxna2", pt.size = 0) + NoLegend() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(text = element_text(size = 20), axis.text = element_text(size = 10)) #1000 x 300

VlnPlot(LUN_Final, "Sox11", pt.size = 0) + NoLegend() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(text = element_text(size = 20), axis.text = element_text(size = 10)) #1000 x 300


Idents(LUN_Final) <- "manual.clusters.1"
VlnPlot(LUN_Final, "Ttc6", pt.size = 0) + NoLegend() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(text = element_text(size = 20), axis.text = element_text(size = 10)) #1000 x 300


VlnPlot(LUN_Final, "Sox14", pt.size = 0) + NoLegend() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(text = element_text(size = 20), axis.text = element_text(size = 10)) #1000 x 300

VlnPlot(LUN_Final, "Plppr1", pt.size = 0) + NoLegend() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(text = element_text(size = 20), axis.text = element_text(size = 10)) #1000 x 300




FeaturePlot(LUN_Final, c("Ttc6"), 
            label = FALSE, label.size = 4, 
            repel = TRUE, order = TRUE, cols = c("grey", "red")) + NoAxes() #1500 x 500


FeaturePlot(LUN_Final, c("Sox14"), 
            label = TRUE, label.size = 4,
            repel = TRUE, order = TRUE, ) + NoAxes() #1500 x 500


FeaturePlot(LUN_Final, c("Rtn4r"), 
            label = TRUE, label.size = 4,
            repel = TRUE, order = TRUE, ) + NoAxes() #1500 x 500

FeaturePlot(LUN_Final, c("Klf7"), 
            label = TRUE, label.size = 4,
            repel = TRUE, order = TRUE, ) + NoAxes() #1500 x 500

FeaturePlot(LUN_Final, c("Klf6"), 
            label = TRUE, label.size = 4,
            repel = TRUE, order = TRUE, ) + NoAxes() #1500 x 500



FeaturePlot(LUN_Final, c("Inpp5k"), 
            label = TRUE, label.size = 4,
            repel = TRUE, order = TRUE, ) + NoAxes() #1500 x 500


FeaturePlot(LUN_Final, c("Nrp1"), 
            label = TRUE, label.size = 4,
            repel = TRUE, order = TRUE, ) + NoAxes() #1500 x 500

FeaturePlot(LUN_Final, c("Plxna2"), 
            label = TRUE, label.size = 4,
            repel = TRUE, order = TRUE, ) + NoAxes() #1500 x 500

FeaturePlot(LUN_Final, c("Prdm6"), 
            label = FALSE, label.size = 4,
            repel = TRUE, order = TRUE, cols = c("grey", "red")) + NoAxes() #1500 x 500


FeaturePlot(LUN_Final, c("Ttc6"), 
            label = FALSE, label.size = 4,
            repel = TRUE, order = TRUE, cols = c("grey", "blue")) + NoAxes() #1500 x 500
