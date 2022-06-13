#The purpose of this script is to generate violin plots for Figure 4 (initial cluster ID)

## This script starts from LUN_Final and an initial ordering determined by 
# LUN_Final_Map1 (for seurat_clusters) and LUN_Final_Order1 for (manual.clusters.1)

# Load libraries
library('Seurat')
library('ggplot2')
library('plyr')
library('dplyr')

#Load relevant data
LUN_Final <- readRDS("Z:/Blackmore_Lab/7_Single_Cell_Manuscript/Github/LUN_Final.rds")

#violin and feature plots for 4A,E,I,M,Q,U

Idents(LUN_Final) <- "seurat_clusters"

# Sim1
VlnPlot(LUN_Final, "Sim1", pt.size = 0) + NoLegend() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 15)) #700 x 300

FeaturePlot(LUN_Final, c("Sim1"), 
            label = FALSE, label.size = 4, 
            repel = TRUE, order = FALSE) +  #300 x 300
  theme(text = element_text(size = 15),
        axis.text = element_text(size = 15))

# Plagl1
VlnPlot(LUN_Final, "Plagl1", pt.size = 0) + NoLegend() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 15)) #700 x 300

FeaturePlot(LUN_Final, c("Plagl1"), 
            label = FALSE, label.size = 4, 
            repel = TRUE, order = FALSE) +  #300 x 300
  theme(text = element_text(size = 15),
        axis.text = element_text(size = 15))


# ttc6
VlnPlot(LUN_Final, "Ttc6", pt.size = 0) + NoLegend() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 15)) #700 x 300

FeaturePlot(LUN_Final, c("Ttc6"), 
            label = FALSE, label.size = 4, 
            repel = TRUE, order = FALSE) +  #300 x 300
  theme(text = element_text(size = 15),
        axis.text = element_text(size = 15))


# Cartpt
VlnPlot(LUN_Final, "Cartpt", pt.size = 0) + NoLegend() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 15)) #700 x 300

FeaturePlot(LUN_Final, c("Cartpt"), 
            label = FALSE, label.size = 4, 
            repel = TRUE, order = TRUE) +  #300 x 300
  theme(text = element_text(size = 15),
        axis.text = element_text(size = 15))

# Slc6a2
VlnPlot(LUN_Final, "Slc6a2", pt.size = 0) + NoLegend() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 15)) #700 x 300

FeaturePlot(LUN_Final, c("Slc6a2"), 
            label = FALSE, label.size = 4, 
            repel = TRUE, order = FALSE) +  #300 x 300
  theme(text = element_text(size = 15),
        axis.text = element_text(size = 15))


# Prdm6
VlnPlot(LUN_Final, "Prdm6", pt.size = 0) + NoLegend() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 15)) #700 x 300

FeaturePlot(LUN_Final, c("Prdm6"), 
            label = FALSE, label.size = 4, 
            repel = TRUE, order = FALSE) +  #300 x 300
  theme(text = element_text(size = 15),
        axis.text = element_text(size = 15))












#3a Dotplot of Hox genes
Idents(LUN_Final) <- "seurat_clusters"

DotPlot(LUN_Final, features = c("Hoxd3", "Hoxb3", "Hoxa3",
                                "Hoxd4", "Hoxc4"),
        dot.scale = 9, dot.min = 0) + coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
        axis.text.y = element_text(size = 30)) #1000 x 800,


#A dimplot to help order the dotplot (put adjacent clusters in order)
DimPlot(LUN_Final, group.by = 'seurat_clusters', label = TRUE, 
        label.size = 5, repel = TRUE, 
        pt.size = .5) + 
  theme(text = element_text(size = 10),
        axis.text = element_text(size = 10)) +
  NoLegend()  #Sized for display

#3b Feature plot of Hox 
FeaturePlot(LUN_Final, c("Hoxd3"), 
            label = TRUE, label.size = 8, 
            repel = TRUE, order = TRUE, pt.size = 2) + 
  theme(text = element_text(size = 30),
        axis.text = element_text(size = 30)) +
  NoLegend()  #Sized for display#1000 x 800


FeaturePlot(LUN_Final, c("Lhx4"), 
            label = TRUE, label.size = 8, 
            repel = TRUE, order = TRUE, pt.size = 2) + 
  theme(text = element_text(size = 30),
        axis.text = element_text(size = 30)) +
  NoLegend()  #Sized for display#1000 x 800

FeaturePlot(LUN_Final, c("Vsx2"), 
            label = TRUE, label.size = 8, 
            repel = TRUE, order = TRUE, pt.size = 2) + 
  theme(text = element_text(size = 30),
        axis.text = element_text(size = 30)) +
  NoLegend()  #Sized for display#1000 x 800


FeaturePlot(LUN_Final, c("Daam2"), 
            label = TRUE, label.size = 8, 
            repel = TRUE, order = TRUE, pt.size = 1) + 
  theme(text = element_text(size = 30),
        axis.text = element_text(size = 30)) +
  NoLegend()  #Sized for display#1000 x 800


FeaturePlot(LUN_Final, c("Nox4"), 
            label = TRUE, label.size = 8, 
            repel = TRUE, order = FALSE, pt.size = 1) + 
  theme(text = element_text(size = 30),
        axis.text = element_text(size = 30)) +
  NoLegend()  #Sized for display#1000 x 800

FeaturePlot(LUN_Final, c("Pard3b"), 
            label = TRUE, label.size = 8, 
            repel = TRUE, order = FALSE, pt.size = 1) + 
  theme(text = element_text(size = 30),
        axis.text = element_text(size = 30)) +
  NoLegend()  #Sized for display#1000 x 800

FeaturePlot(LUN_Final, c("Col19a1"), 
            label = TRUE, label.size = 8, 
            repel = TRUE, order = FALSE, pt.size = 1) + 
  theme(text = element_text(size = 30),
        axis.text = element_text(size = 30)) +
  NoLegend()  #Sized for display#1000 x 800

FeaturePlot(LUN_Final, c("Kit"), 
            label = TRUE, label.size = 8, 
            repel = TRUE, order = TRUE, pt.size = 1) + 
  theme(text = element_text(size = 30),
        axis.text = element_text(size = 30)) +
  NoLegend()  #Sized for display#1000 x 800

FeaturePlot(LUN_Final, c("Gad2"), 
            label = TRUE, label.size = 8, 
            repel = TRUE, order = TRUE, pt.size = 1) + 
  theme(text = element_text(size = 30),
        axis.text = element_text(size = 30)) +
  NoLegend()  #Sized for display#1000 x 800


# Lhx4
VlnPlot(LUN_Final, "Tph2", pt.size = 0) + NoLegend() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 15)) #1000 x 300

FeaturePlot(LUN_Final, c("Lhx4"), 
            label = FALSE, label.size = 4, 
            repel = TRUE, order = TRUE) +  #300 x 300
  theme(text = element_text(size = 15),
        axis.text = element_text(size = 15))

# Vsx2
VlnPlot(LUN_Final, "Vsx2", pt.size = 0) + NoLegend() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 15)) #700 x 300

FeaturePlot(LUN_Final, c("Vsx2"), 
            label = FALSE, label.size = 4, 
            repel = TRUE, order = TRUE) +  #300 x 300
  theme(text = element_text(size = 15),
        axis.text = element_text(size = 15))

# Pard3b
VlnPlot(LUN_Final, "Pard3b", pt.size = 0) + NoLegend() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 15)) #700 x 300

FeaturePlot(LUN_Final, c("Pard3b"), 
            label = FALSE, label.size = 4, 
            repel = TRUE, order = FALSE) +  #300 x 300
  theme(text = element_text(size = 15),
        axis.text = element_text(size = 15))

# Cartpt
VlnPlot(LUN_Final, "Cartpt", pt.size = 0) + NoLegend() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 15)) #700 x 300

FeaturePlot(LUN_Final, c("Cartpt"), 
            label = FALSE, label.size = 4, 
            repel = TRUE, order = FALSE) +  #300 x 300
  theme(text = element_text(size = 15),
        axis.text = element_text(size = 15))

# Slc6a2
VlnPlot(LUN_Final, "Slc6a2", pt.size = 0) + NoLegend() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 15)) #700 x 300

FeaturePlot(LUN_Final, c("Slc6a2"), 
            label = FALSE, label.size = 4, 
            repel = TRUE, order = FALSE) +  #300 x 300
  theme(text = element_text(size = 15),
        axis.text = element_text(size = 15))

# Prdm6
VlnPlot(LUN_Final, "Prdm6", pt.size = 0) + NoLegend() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 15)) #700 x 300

FeaturePlot(LUN_Final, c("Prdm6"), 
            label = FALSE, label.size = 4, 
            repel = TRUE, order = FALSE) +  #300 x 300
  theme(text = element_text(size = 15),
        axis.text = element_text(size = 15))

# Lhx4
VlnPlot(LUN_Final, "Lhx4", pt.size = 0) + NoLegend() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 15)) #700 x 300

FeaturePlot(LUN_Final, c("Lhx4"), 
            label = FALSE, label.size = 4, 
            repel = TRUE, order = TRUE) +  #300 x 300
  theme(text = element_text(size = 15),
        axis.text = element_text(size = 15))

# Lhx4
VlnPlot(LUN_Final, "Nkx2-2", pt.size = 0) + NoLegend() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 15)) #700 x 300

FeaturePlot(LUN_Final, c("Lhx4"), 
            label = FALSE, label.size = 4, 
            repel = TRUE, order = TRUE) +  #300 x 300
  theme(text = element_text(size = 15),
        axis.text = element_text(size = 15))




levels(LUN_Final)

# Identify markers for the seurat clusters
#Naming scheme is to give file name, then Markers, then
#slot for grouping, then what is compared (all)

LUN_Final_Markers_Seurat_All <- FindAllMarkers(LUN_Final, only.pos = TRUE)

write.csv(x = LUN_Final_Markers_Seurat_All, 
          file = "Z:/Blackmore_Lab/Zac/1SingleCell/1/Results/LUN_Final_Markers_Seurat_All.csv",
          quote = FALSE)


#Assign LUN_Final seurat clusters to anatomical names

#Make a "map" file and save it in Lists folder.

LUN_Final_Map1 <- read.csv("Z:/Blackmore_Lab/Zac/1SingleCell/1/Inputs/LUN_Final_Map1.csv",
                           sep=",")

LUN_Final@meta.data$manual.clusters <- plyr::mapvalues(
  x = LUN_Final@meta.data$seurat_clusters,
  from = LUN_Final_Map1$Cluster,
  to = LUN_Final_Map1$ID1)

levels(LUN_Final$manual.clusters)

LUN_Final_Order <- read.csv("Z:/Blackmore_Lab/Zac/1SingleCell/1/Inputs/LUN_Final_Order.csv",
                            sep=",")

LUN_Final@meta.data$manual.clusters <- factor(
  x = LUN_Final@meta.data$manual.clusters,
  levels = LUN_Final_Order$ID2)

levels(LUN_Final$manual.clusters)

Idents(LUN_Final) <- "manual.clusters"

DotPlot(LUN_Final, features = FinalDP$Gene, dot.min = 0, dot.scale = 5) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

FeaturePlot(LUN_Final, c("Kit"), 
            label = TRUE, label.size = 3, 
            repel = TRUE, order = TRUE) #500 x 500

#Now that the names are logically ordered, make a new diff expression file
#with the proper names. Call this "LUN_Final_Manual_all" 

Idents(LUN_Final) <- "manual.clusters"

LUN_Final_Manual_all <- FindAllMarkers(LUN_Final, only.pos = TRUE)

write.csv(x = LUN_Final_Manual_all, 
          file = "Z:/Blackmore_Lab/Zac/1SingleCell/1/Results/LUN_Final_Manual_all.csv",
          quote = FALSE)


#Test the specificity of genes by making dot plots using a "test" file

Test <- read.csv("Z:/Blackmore_Lab/Zac/1SingleCell/1/Inputs/Test.csv",
                 sep=",")

DotPlot(LUN_Final, features = Test$Genes, dot.min = 0, dot.scale = 5) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

LUN_Final_Final <- read.csv("Z:/Blackmore_Lab/Zac/1SingleCell/1/Inputs/LUN_Final_Final1.csv",
                            sep=",")

DotPlot(LUN_Final, features = LUN_Final_Final$Gene, dot.min = 0, dot.scale = 5) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


FeaturePlot(LUN_Final, c("Pard3b"), 
            label = TRUE, label.size = 3, 
            repel = FALSE, order = TRUE) #500 x 500

FeaturePlot(LUN_Final, c("Col19a1"), 
            label = TRUE, label.size = 3, 
            repel = FALSE, order = TRUE) #500 x 500

FeaturePlot(LUN_Final, c("Kit"), 
            label = TRUE, label.size = 3, 
            repel = FALSE, order = TRUE) #500 x 500

FeaturePlot(LUN_Final, c("Dsg2"), 
            label = FALSE, label.size = 3, 
            repel = FALSE, order = TRUE) #500 x 500

FeaturePlot(LUN_Final, c("Lhx4"), 
            label = FALSE, label.size = 3, 
            repel = FALSE, order = TRUE) #500 x 500

DotPlot(LUN_Final, features = AllTF1$TF, dot.min = 0, dot.scale = 5) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))