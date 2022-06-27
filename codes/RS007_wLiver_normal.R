# 2021 June 15 created 
# script is the Seurat4 analysis of the RS007

setwd("D:/pCloud Sync/Columbia/SchwabeLab/Aveline_HSCnHCC/github/HSCinHCC/") #ajay

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

# *****************************************read the data
data <- Read10X_h5("data/RS007_wLiver_normal/RS007_filtered_gene_bc_matrices_h5.h5")

rs007 <- CreateSeuratObject(counts = data, project = "wLiver_normal", min.cells = 3, min.features = 200)
rs007
rs007[["percent.mt"]] <- PercentageFeatureSet(rs007, pattern = "^mt-")
rs007[["percent.ribo"]] <- PercentageFeatureSet(rs007, pattern = "^Rpl|^Rps|Mrps|^Mrpl")#

# Visualize QC metrics as a violin plot
VlnPlot(rs007, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo"), ncol = 4)
FeatureScatter(rs007, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(rs007, feature1 = "nFeature_RNA", feature2 = "percent.mt")
FeatureScatter(rs007, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p1 <- FeatureScatter(rs007, feature1 = "nCount_RNA", feature2 = "percent.mt") + NoLegend()
p2 <- FeatureScatter(rs007, feature1 = "nFeature_RNA", feature2 = "percent.mt") + NoLegend()
p3 <- FeatureScatter(rs007, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
wrap_plots(list(p1,p2,p3),nrow = 2,ncol = 2)
# keep only the required cells

# Normalize and cluster seurat object

DimPlot(rs007,label = T)
DimPlot(rs007,label = F)


# selecting cell groups
# https://satijalab.org/seurat/articles/visualization_vignette.html
plot <- DimPlot(rs007, reduction = "umap")
select.cells <- CellSelector(plot = plot)
head(select.cells)
Idents(rs007, cells = select.cells) <- 21
Idents(rs007)
DimPlot(rs007,label = T)

FeaturePlot(object = rs007, features = c("percent.mt"))
VlnPlot(rs007, features = c("percent.mt"))
FeaturePlot(object = rs007, features = c("percent.ribo"))
FeaturePlot(object = rs007, features = c("nCount_RNA"))
FeaturePlot(object = rs007, features = c("nFeature_RNA"))





# ******************************************************************saved object
rs007 <- readRDS(file = "results/RS007/seurat4_rs007_wliver.rds")
DimPlot(rs007,label = T)
DimPlot(rs007,label = T, group.by = "Cell.types")
DimPlot(rs007,label = T, group.by = "Cell.clusters")

# get only the clean clusters
rs007 <- SetIdent(rs007, value = "Cell.types")
req_subset <- subset(rs007, idents = c("Multi"), invert=T)
req_subset
DimPlot(req_subset,label = T, group.by = "ident",pt.size = 1)
DimPlot(req_subset,label = T, group.by = "Cell.types")
DimPlot(req_subset,label = T, group.by = "Cell.clusters")
DimPlot(req_subset,label = T, group.by = "Cell.typesFine")



#******************************************************#******************************************************
#******************************************************PLOTS
#******************************************************#******************************************************
library(RColorBrewer)
library("ggsci")
library("scales")


# CELL TYPES UMAP
data_req <- req_subset
DimPlot(object = data_req, group.by = "ident", label = T)
levels(data_req@active.ident)

DimPlot(object = data_req, group.by = "ident", label = F, cols = cols, pt.size = 2) 

DimPlot(object = data_req, group.by = "Cell.types", label = T, cols = cols, pt.size = 1)  + labs(title = "Mouse normal whole-liver")
ggsave(filename = paste0("results/RS007/figs/mus_RS007_wliver_normal_celltypes_label.pdf"))
DimPlot(object = data_req, group.by = "Cell.types", label = F, cols = cols, pt.size = 2) + NoLegend() + labs(title = "Mouse normal whole-liver")
ggsave(filename = paste0("results/RS007/figs/mus_RS007_wliver_normal_celltypes_noLabel.pdf"))


# CELL TYPES DOT PLOT
markerGenes <- c("Alb","Serpina1a",
                 # 'Epcam', 'Krt7',
                 "Lrat","Rgs5","Timp1","Col1a1", #,"Acta2","Lox"
                 "Kdr","Aqp1",
                 # 'Upk3b','Msln',
                 "Il7r","Cd3d",
                 "Nkg7",
                 "Top2a","Mki67",
                 "Cd79a","Ebf1",
                 "Mzb1","Jchain",
                 "Hbb-bt","Snca",
                 "Lgals3","Ccr2",
                 "Clec4f","Vsig4",
                 "Siglech","Klk1")               
DotPlot(data_req, features = markerGenes, group.by = "Cell.types") + RotatedAxis() & scale_color_distiller(palette = "RdBu")
ggsave(filename = "results/RS007/figs/mus_RS007_wliver_normal_celltypes_dotplot_markers.pdf", width = 9, height = 4)



