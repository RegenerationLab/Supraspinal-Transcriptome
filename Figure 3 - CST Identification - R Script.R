#The purpose of this script is to establish CST identity in the UNINJURED dataset


#Positive markers for layer V cortical / CST: Crym, Fezf2, Bcl11b, 
# Bcl6, Pdlim1, Cacna1h, Mylip, Kcng1
#markers I considered but didn't use:(Slco2a1, Crim1, Gng7, Cthrc1)

#Positive markers for Cortex: Satb2, Iptka

#Negative markers for CST: Sulf1, Rims3, Hpgd (thalamus projecting, didn't work)


#Load relevant data
LUN1 <- readRDS("Z:/Blackmore_Lab/scRNAseq_Analyses/1_StndID/LUN/LUN1.rds")



# Make a dimplot to remind me of which clusters are which

Idents(LUN1) <- "seurat_clusters"
DimPlot(LUN1, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()


# Example Feature plots


FeaturePlot(LUN1, c("Bcl11b"), label = TRUE, pt.size = 1,
            label.size = 10, repel = TRUE, order = FALSE) +
  theme(text = element_text(size = 25),
        axis.text = element_text(size = 40)) +
  NoLegend() + NoAxes() #1000 x 1000

FeaturePlot(LUN1, c("Crym"), label = TRUE, pt.size = 1,
            label.size = 10, repel = TRUE, order = FALSE) +
  theme(text = element_text(size = 25),
        axis.text = element_text(size = 40)) +
  NoLegend() + NoAxes() #1000 x 1000

FeaturePlot(LUN1, c("Fezf2"), label = TRUE, pt.size = 1,
            label.size = 10, repel = TRUE, order = TRUE) +
  theme(text = element_text(size = 25),
        axis.text = element_text(size = 40)) +
  NoLegend() + NoAxes() #1000 x 1000



FeaturePlot(LUN1, c("Bcl11b"), 
            label = TRUE, label.size = 5, split.by = "orig.ident",
            repel = TRUE, order = TRUE) + NoLegend() #1000 x 500 then 2x in adobe

FeaturePlot(LUN1, c("Crym"), 
            label = TRUE, label.size = 5, split.by = "orig.ident",
            repel = TRUE, order = TRUE, ) #1000 x 500 then 2x in adobe


FeaturePlot(LUN1, c("Fezf2"), 
            label = TRUE, label.size = 5, split.by = "orig.ident",
            repel = TRUE, order = TRUE, ) #1000 x 500 then 2x in adobe


#Example Violin Plots

VlnPlot(LUN1, "Bcl11b", pt.size = 0) + NoLegend() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0)) +
  theme(text = element_text(size = 20), 
        axis.text = element_text(size = 20)) #1000 x 300

VlnPlot(LUN1, "Crym", pt.size = 0) + NoLegend() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0)) +
  theme(text = element_text(size = 20), axis.text = element_text(size = 20)) #1000 x 300

VlnPlot(LUN1, "Fezf2", pt.size = 0) + NoLegend() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0)) +
  theme(text = element_text(size = 20), axis.text = element_text(size = 20)) #1000 x 300

VlnPlot(LUN1, "Bcl6", pt.size = 0) + NoLegend() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0)) +
  theme(text = element_text(size = 20), axis.text = element_text(size = 20)) #1000 x 300

VlnPlot(LUN1, "Pdlim1", pt.size = 0) + NoLegend() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0)) +
  theme(text = element_text(size = 20)) #1000 x 300

VlnPlot(LUN1, "Cacna1h", pt.size = 0) + NoLegend() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0)) +
  theme(text = element_text(size = 20)) #1000 x 300

VlnPlot(LUN1, "Mylip", pt.size = 0) + NoLegend() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0)) +
  theme(text = element_text(size = 20)) #1000 x 300

VlnPlot(LUN1, "Kcng1", pt.size = 0) + NoLegend() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0)) +
  theme(text = element_text(size = 20)) #1000 x 300






