# 2022 April 28 created 
# this is the integration of nuc-seq human liver data (normal, nafld, HCC, tumor adjacent) after cellbender


library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

setwd("D:/pCloud Sync/Columbia/SchwabeLab/Aveline_HSCnHCC/github/HSCinHCC")

# the  data
markerGenes <- c("ALB","APOE",
                 "ANXA4","EPCAM", "CFTR",#"KRT7","KRT8",'KRT9','KRT19','SOX9',
                 "COL4A1","LAMA2","CCBE1","HGF", #CALD1
                 "PTPRB","PECAM1","FLT1",
                 'UPK3B','MSLN',
                 "IL7R","THEMIS","FYN",#"AOAH",
                 "NKG7",
                 "TOP2A","MKI67",
                 "EBF1","IGKC","IGHA1",
                 "MZB1","JCHAIN",
                 "HBB","HBA1","SNCA",
                 "LGALS3","CCR2",
                 # "Clec4f","Vsig4",
                 "CD163","VSIG4",'DMXL2' )

# markerGenes <- c("MYH11","ACTG2","ACTA2",
#                  "RGS5","RERGL","COLEC11","LRAT","HGF",
#                  "LUM","COL1A1","COL3A1","LOX","TIMP1","DCN",
#                  "DPT","GPX3","MFAP4","GSN", "MGP","ELN","ATP1A2","DCLK1","CLIC5","NAKIN3","FLT1","SVEP1", #,"Gsn", "Mfap4"
#                  "KDR","AQP1","VWF",
#                  "ENG","PECAM1","RAMP3","INMT",#Portal Endo
#                  "STAB1",#CV LSEC
#                  "SPARCL1","CLEC14A",#PERIPortal  LSEC
#                  "TOP2A","MKI67","CENPF",
#                  'UPK3B','MSLN', 
#                  # 'EPCAM', 'KRT7',"ALB","SERPINA1A",
#                  "ALB","APOE","APOB","CPS1",
#                  "ANXA4","EPCAM", "CFTR","KRT7","KRT8",'KRT19','SOX9',
#                  "IL7R","CD8A","CD3D","TRAC","FOXP3","THEMIS",
#                  "NKG7","GNLY",
#                  "MS4A1","CD79A","IGKC",
#                  "JCHAIN",
#                  "HBB","HBA1","SNCA",
#                  "CLEC4F","VSIG4","LGALS3", #Kupffer cells
#                  "ADGRE1","CD68","CD163","ITGB2","FCGR2A","ITGAM",
#                  "FCGR3A","MS4A7", # FCGR3A + Monocytes
#                  "CD14","LYZ", # LyZ monocytes
#                  "FCER1A", "CST3") # Dendritic cells

# **********************************Reading multiple samples
# all cell bender samples
loc <- "data/hum_nucseq_RS_NH_Andrews/cellBender"
dir <- list.dirs(path = loc,full.names = F, recursive = F)
dir
i=dir[1]
seurat.list <- vector(mode = "list", length = length(dir)) 
names(seurat.list) <- dir
seurat.list
# name.list <- rep(NA, length(dir))

for (i in 1:length(dir)) {
  if (file.exists(file.path(loc, dir[i], "CB_raw_feature_bc_matrix_filtered.h5"))){
    sampleName <- dir[i]
    
    # name.list[i] <- sampleName
    print(sampleName)
    dat_filt <- Read10X_h5(file.path(loc, sampleName, "CB_raw_feature_bc_matrix_filtered.h5"))
    dim(dat_filt)
    
    seurat <- CreateSeuratObject(counts = dat_filt, project = sampleName, min.cells = 3, min.features = 200)
    seurat
    seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")
    seurat[["percent.ribo"]] <- PercentageFeatureSet(seurat, pattern = "^RPL|^RPS|MRPS|^MRPL")#
    seurat <- subset(seurat, subset = percent.mt < 20 )#& nFeature_RNA > 200 & nFeature_RNA < 6500 & nCount_RNA < 40000 )#&
    seurat.list[i] <- seurat
    # names(seurat.list[i]) <- sampleName
  }
}
length(seurat.list)
seurat.list[sapply(seurat.list, is.null)] <- NULL
seurat.list


# *************************************************************
# combining all the samples # *************************************************************
tmp <- names(seurat.list)
tmp 
# required samples
req_samples <-tmp #[-which(tmp %in% c("S41","S42","S43","S44"))] 
req_sample_pos <- which(tmp %in% req_samples)

req.list <- seurat.list[req_sample_pos]
names(req.list)

req.list <- lapply(X = req.list, FUN = function(x) {
  print(unique(x@meta.data$orig.ident))
  # x <- subset(x, subset = percent.mt < 40 )#& nFeature_RNA > 200 & nFeature_RNA < 6500 & nCount_RNA < 40000 )#&
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, verbose = FALSE, nfeatures = 3000)
})

names(req.list)

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = req.list, nfeatures = 3000)

# scale and compute pca for RPCA integration method
req.list <- lapply(X = req.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})
# perform RPCA integration of all samples
data.anchors <- FindIntegrationAnchors(object.list = req.list, reduction = "rpca", anchor.features = features) # 
# create an 'integrated' data assay
data.combined_all <- IntegrateData(anchorset = data.anchors) #, dims = 1:50

DefaultAssay(data.combined_all) <- "integrated"
# Run the standard workflow for visualization and clustering
data.combined_all <- ScaleData(data.combined_all, verbose = FALSE)
data.combined_all <- RunPCA(data.combined_all, npcs = 50, verbose = FALSE)
ElbowPlot(data.combined_all, ndims = 50)
ndim <- 14
data.combined_all <- RunUMAP(data.combined_all, reduction = "pca", dims = 1:ndim)
data.combined_all <- FindNeighbors(data.combined_all, reduction = "pca", dims = 1:ndim)
data.combined_all <- FindClusters(data.combined_all)#, resolution = 0.5
DimPlot(data.combined_all, label = T)
DimPlot(data.combined_all,group.by = "orig.ident", label = F,shuffle = T)
FeaturePlot(object = data.combined_all, features = c("percent.mt"))
FeaturePlot(object = data.combined_all, features = c("percent.ribo"))
FeaturePlot(object = data.combined_all, features = c("nCount_RNA"))
FeaturePlot(object = data.combined_all, features = c("nFeature_RNA"))
VlnPlot(object = data.combined_all, features = c("percent.mt"))
VlnPlot(object = data.combined_all, features = c("percent.ribo"))
VlnPlot(object = data.combined_all, features = c("nCount_RNA"))
VlnPlot(object = data.combined_all, features = c("nFeature_RNA"))

DefaultAssay(data.combined_all) <- "RNA"
DotPlot(data.combined_all, features = markerGenes) + RotatedAxis()
DotPlot(data.combined_all, features = markerGenes, cluster.idents = T) + RotatedAxis()

FeaturePlot(object = data.combined_all, features = c("HGF","COL1A1","COLEC11","MYH11"), max.cutoff = "q99", pt.size = 0.5,order = T)
FeaturePlot(object = data.combined_all, features = c("VCAN","DCN","RGS5","LUM"), max.cutoff = "q99", pt.size = 0.5,order = T)
FeaturePlot(object = data.combined_all, features = c("PDGFRA","PDGFRB","LRAT","DES"), max.cutoff = "q99", pt.size = 0.5,order = T)
date()

table(data.combined_all@meta.data$orig.ident)
# saving the cluster numbers
# data.combined_all[["Cell.clusters"]] <- Idents(object = data.combined_all)

#*# selecting cell groups from the original clusters, manually assinging to new clusters based on gene expression and subclusters and reordering the cluster numbers
DimPlot(data.combined_all, reduction = "umap",label = T,pt.size = 1, shuffle = T,group.by = "Cell.clusters")
# CellsByIdentities(data.combined_all, idents = 1, cells = NULL, return.null = FALSE) #getting cells by their identities

# data.combined_all <- SetIdent(data.combined_all, value = "Cell.clusters")
identSel <- c(22,3)
DimPlot(data.combined_all, reduction = "umap", cells = unlist(CellsByIdentities(data.combined_all, idents = identSel)))
FeaturePlot(data.combined_all,features = "MKI67",cells = unlist(CellsByIdentities(data.combined_all, idents = identSel)))
FeaturePlot(data.combined_all,features = "IL7R",cells = unlist(CellsByIdentities(data.combined_all, idents = identSel)))
FeaturePlot(data.combined_all,features = "THEMIS",cells = unlist(CellsByIdentities(data.combined_all, idents = identSel)))
FeaturePlot(data.combined_all,features = "CPS1",cells = unlist(CellsByIdentities(data.combined_all, idents = identSel)))
FeaturePlot(data.combined_all,features = "CAMK4",cells = unlist(CellsByIdentities(data.combined_all, idents = identSel)))
FeaturePlot(object = data.combined_all, features = c("percent.mt"), cells = unlist(CellsByIdentities(data.combined_all, idents = identSel)))
FeaturePlot(object = data.combined_all, features = c("percent.ribo"), cells = unlist(CellsByIdentities(data.combined_all, idents = identSel)))
FeaturePlot(object = data.combined_all, features = c("nCount_RNA"), cells = unlist(CellsByIdentities(data.combined_all, idents = identSel)))
FeaturePlot(object = data.combined_all, features = c("nFeature_RNA"), cells = unlist(CellsByIdentities(data.combined_all, idents = identSel)))

# manually assigning separately clustering cells to new clusters
# https://satijalab.org/seurat/articles/visualization_vignette.html
plot <- DimPlot(data.combined_all, reduction = "umap", cells = unlist(CellsByIdentities(data.combined_all, idents = identSel)))
select.cells <- CellSelector(plot = plot)
head(select.cells)
Idents(data.combined_all, cells = select.cells) <- 35
# Idents(data.combined_all)
DimPlot(data.combined_all,label = T, pt.size = 1)
table((data.combined_all@active.ident))
levels(data.combined_all@active.ident)
data.combined_all@active.ident <- factor(data.combined_all@active.ident, levels = as.character(c(0:35)))
DimPlot(data.combined_all,label = T, pt.size = 1)
data.combined_all[["Cell.clusters.new"]] <- Idents(object = data.combined_all)

# get cluster markers
data.combined_all <- SetIdent(data.combined_all, value = "Cell.clusters.new")
DimPlot(data.combined_all,label = T, pt.size = 1)
DefaultAssay(data.combined_all) <- "RNA"
DotPlot(data.combined_all, features = markerGenes) + RotatedAxis()
table(data.combined_all@active.ident)

# get the marker genes for each cluster
markers <- FindAllMarkers(data.combined_all,only.pos = T, logfc.threshold = 0.25,max.cells.per.ident = 1000)
markers <- markers[order(markers$cluster,markers$avg_log2FC,decreasing = T),]
head(markers,20)
# write.csv(markers, file = "results/hum_wLiv/markers_humLivNucseq_Norm_Nafld_HCC_tumorAdj_clusters.csv")


# giving cell type annotation
DimPlot(data.combined_all,label = T)
DimPlot(data.combined_all,label = T,group.by = "Cell.clusters.new")
data.combined_all <- SetIdent(data.combined_all, value = "Cell.clusters.new")
DimPlot(data.combined_all,label = T)
data.combined_all <- RenameIdents(object = data.combined_all, 
                        '0' = "Hepatocytes",
                        '1' = "Hepatocytes",
                        '2' = "Hepatocytes",
                        '4' = "Hepatocytes",
                        '5' = "Hepatocytes",
                        '6' = "Hepatocytes",
                        '7' = "Hepatocytes",
                        '8' = "Hepatocytes",
                        '9' = "Hepatocytes",
                        '13' = "Hepatocytes",
                        '16' = "Hepatocytes",
                        '17' = "Hepatocytes",
                        '18' = "Hepatocytes",
                        '20' = "Hepatocytes",
                        '28' = "Hepatocytes",
                        '29' = "Hepatocytes",
                        '30' = "Hepatocytes",
                        '34' = "Hepatocytes",
                        
                        '19' = "Cycling",
                        '33' = "Cycling",
                        '34' = "Cycling",
                        '35' = "Cycling",
                        
                        '11' = "Cholangiocytes",
                        '26' = "Cholangiocytes",
                        '27' = "Cholangiocytes",
                        
                        '12' = "Fibroblasts",
                        '21' = "Fibroblasts",
                        
                        '14' = "Endothelial",
                        '15' = "Endothelial",
                        
                        '10' = "Myeloid",
                        
                        '23' = "B cell",
                        '24' = "Plasma",
                        
                        '3' = "T_NK_NKT",
                        
                        '22' = "Mixed",
                        '25' = "Mixed",
                        '31' = "Mixed",
                        '32' = "Mixed")






data.combined_all[["Cell.types"]] <- Idents(object = data.combined_all)
DimPlot(data.combined_all,label = T, group.by = "Cell.types")
table(data.combined_all@meta.data$Cell.types)
table(data.combined_all@meta.data$Cell.clusters)
# table(data.combined_all@meta.data$cell.type)

DimPlot(data.combined_all, reduction = "umap",group.by = "Cell.types", label = F)
DimPlot(data.combined_all, reduction = "umap",group.by = "Cell.clusters.new", label = F)

DotPlot(data.combined_all, features = markerGenes) + RotatedAxis()

# giving disease type annotations
tmp <- data.combined_all[["orig.ident"]]
table(tmp$orig.ident)

tmppos <- grep("TST",tmp$orig.ident)
table(tmp$orig.ident[tmppos])
tmp$orig.ident[tmppos] <- "Normal"

tmppos <- grep("nafld",tmp$orig.ident)
table(tmp$orig.ident[tmppos])
tmp$orig.ident[tmppos] <- "NAFLD"

tmppos <- which(tmp$orig.ident %in% c("S41","S43"))
table(tmp$orig.ident[tmppos])
tmp$orig.ident[tmppos] <- "Adj.tumor"

tmppos <- which(tmp$orig.ident %in% c("S42","S44"))
table(tmp$orig.ident[tmppos])
tmp$orig.ident[tmppos] <- "Tumor"

table(tmp$orig.ident)

data.combined_all[["DiseaseTypes"]] <- tmp$orig.ident
DimPlot(data.combined_all, reduction = "umap",group.by = "DiseaseTypes", label = F, shuffle = T)

# giving disease type annotations -simple
tmp <- data.combined_all[["orig.ident"]]
table(tmp$orig.ident)

tmppos <- grep("TST",tmp$orig.ident)
table(tmp$orig.ident[tmppos])
tmp$orig.ident[tmppos] <- "Normal"

tmppos <- grep("nafld",tmp$orig.ident)
table(tmp$orig.ident[tmppos])
tmp$orig.ident[tmppos] <- "Fibrotic"

tmppos <- which(tmp$orig.ident %in% c("S41","S43"))
table(tmp$orig.ident[tmppos])
tmp$orig.ident[tmppos] <- "Fibrotic"

tmppos <- which(tmp$orig.ident %in% c("S42","S44"))
table(tmp$orig.ident[tmppos])
tmp$orig.ident[tmppos] <- "Tumor"

table(tmp$orig.ident)

data.combined_all[["DiseaseTypesSimple"]] <- tmp$orig.ident
DimPlot(data.combined_all, reduction = "umap",group.by = "DiseaseTypesSimple", label = F, shuffle = T)




# saveRDS(data.combined_all, file = "results/hum_wLiv/seuratObj_humLivNucseq_Norm_Nafld_HCC_tumorAdj_mito20.rds")


# #******************************************************#******************************************************
# #*************************************************** reclustering T-NK-NKT to identify doublets
# #******************************************************#******************************************************
# data.combined_all <- readRDS(file = "results/hum_wLiv/seuratObj_humLivNucseq_Norm_Nafld_HCC_tumorAdj_mito20.rds")
# DimPlot(data.combined_all, label = T)
# DimPlot(data.combined_all, reduction = "umap",group.by = "Cell.clusters.new", label = T)
# 
# data.combined_all <- SetIdent(data.combined_all, value = "Cell.clusters.new")
# req_subset <- subset(data.combined_all, idents = c("3", "18","22", "25"), invert=F)
# req_subset
# DimPlot(req_subset, label = T)
# # RECLUSTERING
# DefaultAssay(req_subset) <- "integrated"
# req_subset <- ScaleData(req_subset, features = VariableFeatures(object = req_subset))
# req_subset <- RunPCA(req_subset, features = VariableFeatures(object = req_subset))
# ElbowPlot(req_subset, ndims = 50)
# ndim=10
# req_subset <- FindNeighbors(req_subset, dims = 1:ndim)
# req_subset <- FindClusters(req_subset)
# req_subset <- RunUMAP(req_subset, reduction = "pca", dims = 1:ndim)
# DimPlot(req_subset,label = T, pt.size = 3)
# DimPlot(req_subset,label = T, group.by = "orig.ident")
# DimPlot(req_subset,label = T, group.by = "Cell.clusters.new")
# 
# DefaultAssay(req_subset) <- "RNA"
# DotPlot(req_subset, features = markerGenes) + RotatedAxis()
# 
# req_subset[["Cell.clusters.TNKNKT"]] <- Idents(object = req_subset)
# DimPlot(req_subset,label = T,group.by = "Cell.clusters.TNKNKT")
# req_subset <- SetIdent(req_subset, value = "Cell.clusters.TNKNKT")
# DimPlot(req_subset,label = T)
# req_subset <- RenameIdents(object = req_subset, 
#                            '0' = "T_NK_NKT",
#                            '1' = "T_NK_NKT",
#                            '2' = "T_NK_NKT",
#                            '3' = "T_NK_NKT",
#                            '5' = "T_NK_NKT",
#                            # '8' = "T_NK_NKT",
#                            '9' = "T_NK_NKT",
#                            
#                            '14' = "Cycling",
#                            
#                            '6' = "Hepatocytes",
#                            # '7' = "Hepatocytes",
#                            '10' = "Hepatocytes",
#                            '13' = "Hepatocytes",
#                            
#                            
#                            '7' = "Mixed",
#                            '8' = "Mixed",
#                            '11' = "Mixed",
#                            '12' = "Mixed",
#                            '15' = "Mixed",
#                            '4' = "Mixed")
# 
# 
# 
# 
# 
# 
# req_subset[["Cell.types.TNKNKT"]] <- Idents(object = req_subset)
# DimPlot(req_subset,label = T, group.by = "Cell.types.TNKNKT")
# DimPlot(req_subset,label = T, group.by = "Cell.clusters.new")


# # GETTING THIS ANNOTATION BACK TO ORIGINAL DATA
# DimPlot(data.combined_all, label = T)
# data.combined_all[["Cell.types.new"]] <-data.combined_all[["Cell.types"]]
# DimPlot(data.combined_all, reduction = "umap",group.by = "Cell.types.new", label = T)
# data.combined_all <- SetIdent(data.combined_all, value = "Cell.types.new")
# DimPlot(data.combined_all, label = T)
# # T_NK_NKT
# tmppos <- which(req_subset@meta.data$Cell.types.TNKNKT=="T_NK_NKT")
# tmpCells <- rownames(req_subset@meta.data)[tmppos]
# Idents(object = data.combined_all, cells = tmpCells) <- "T_NK_NKT"
# DimPlot(data.combined_all,label = T)
# # Cycling
# tmppos <- which(req_subset@meta.data$Cell.types.TNKNKT=="Cycling")
# tmpCells <- rownames(req_subset@meta.data)[tmppos]
# Idents(object = data.combined_all, cells = tmpCells) <- "Cycling"
# DimPlot(data.combined_all,label = T)
# # Mixed
# tmppos <- which(req_subset@meta.data$Cell.types.TNKNKT=="Mixed")
# tmpCells <- rownames(req_subset@meta.data)[tmppos]
# Idents(object = data.combined_all, cells = tmpCells) <- "Mixed"
# DimPlot(data.combined_all,label = T)
# # Hepatocytes
# tmppos <- which(req_subset@meta.data$Cell.types.TNKNKT=="Hepatocytes")
# tmpCells <- rownames(req_subset@meta.data)[tmppos]
# Idents(object = data.combined_all, cells = tmpCells) <- "Hepatocytes"
# DimPlot(data.combined_all,label = T)
# 
# # correcting Hepatocytes to Epithelial_Hep_tumor
# tmppos <- which(data.combined_all@meta.data$Cell.types.new=="Hepatocytes")
# tmpCells <- rownames(data.combined_all@meta.data)[tmppos]
# Idents(object = data.combined_all, cells = tmpCells) <- "Epithelial_Hep_tumor"
# DimPlot(data.combined_all,label = T)
# # # writing to a new metadata column
# data.combined_all[["Cell.types.new"]] <- Idents(object = data.combined_all)
# DimPlot(data.combined_all,label = T, group.by = "Cell.types.new")
# data.combined_all[["Cell.types.new"]] <- NULL #no longer used

# saveRDS(data.combined_all, file = "results/hum_wLiv/seuratObj_humLivNucseq_Norm_Nafld_HCC_tumorAdj_mito20.rds")

#******************************************************#******************************************************
#*************************************************** reclustering Fibroblasts to identify cy-myHSC
#******************************************************#******************************************************
data.combined_all <- readRDS(file = "results/hum_wLiv/seuratObj_humLivNucseq_Norm_Nafld_HCC_tumorAdj_mito20.rds")
DimPlot(data.combined_all, label = T)
DimPlot(data.combined_all, reduction = "umap",group.by = "Cell.types", label = T)

data.combined_all <- SetIdent(data.combined_all, value = "Cell.types")
req_subset <- subset(data.combined_all, idents = c("Fibroblasts"), invert=F)
req_subset
DimPlot(req_subset, label = T)
# RECLUSTERING
DefaultAssay(req_subset) <- "integrated"
req_subset <- ScaleData(req_subset, features = VariableFeatures(object = req_subset))
req_subset <- RunPCA(req_subset, features = VariableFeatures(object = req_subset))
ElbowPlot(req_subset, ndims = 50)
ndim=20
req_subset <- FindNeighbors(req_subset, dims = 1:ndim)
req_subset <- FindClusters(req_subset)
req_subset <- RunUMAP(req_subset, reduction = "pca", dims = 1:ndim)
DimPlot(req_subset,label = T, pt.size = 3)
DimPlot(req_subset,label = T, group.by = "orig.ident")
DimPlot(req_subset,label = T, group.by = "Cell.clusters.new")
DimPlot(req_subset, reduction = "umap",group.by = "DiseaseTypes", label = F,shuffle = T)
DimPlot(req_subset, reduction = "umap",group.by = "DiseaseTypesSimple", label = F,shuffle = T)

DefaultAssay(req_subset) <- "RNA"
DotPlot(req_subset, features = markerGenes) + RotatedAxis()
FeaturePlot(object = req_subset, features = c("HGF","COL1A1","COLEC11","MYH11"), max.cutoff = "q99", pt.size = 0.5,order = T)
FeaturePlot(object = req_subset, features = c("VCAN","DCN","RGS5","LUM"), max.cutoff = "q99", pt.size = 0.5,order = T)
FeaturePlot(object = req_subset, features = c("PDGFRA","PDGFRB","LRAT","DES"), max.cutoff = "q99", pt.size = 0.5,order = T)
FeaturePlot(object = req_subset, features = c("ACTA2","CTHRC1","LOX","TIMP1"), max.cutoff = "q99", pt.size = 0.5,order = T)
FeaturePlot(object = req_subset, features = c("TGFBR1","COL3A1","COL15A1","ACTG1"), max.cutoff = "q99", pt.size = 0.5,order = T)
FeaturePlot(object = req_subset, features = c("ROBO2","RELN","HHIP","COL25A1"), max.cutoff = "q99", pt.size = 0.5,order = T)
FeaturePlot(object = req_subset, features = c("DDR1","HAS3","HAS1","HAS2"), max.cutoff = "q99", pt.size = 0.5,order = T)
FeaturePlot(object = req_subset, features = c("LMNA","MMP19","MMP2","MMP16"), max.cutoff = "q99", pt.size = 0.5,order = T)
FeaturePlot(object = req_subset, features = c("MMP14","MMP24","XIST","C1S"), max.cutoff = "q99", pt.size = 0.5,order = T)
FeaturePlot(object = req_subset, features = c("HGF","MKI67","TOP2A","CENPF"), max.cutoff = "q99", pt.size = 0.5,order = T)

# manually assigning separately clustering cells to new clusters
# https://satijalab.org/seurat/articles/visualization_vignette.html
req_subset <- SetIdent(req_subset, value = "Cell.clusters.Fib")
DimPlot(req_subset,label = T)
identSel <- c(1,3,5,0,7,8,9,2)
plot <- DimPlot(req_subset, reduction = "umap", cells = unlist(CellsByIdentities(req_subset, idents = identSel)))
select.cells <- CellSelector(plot = plot)
head(select.cells)
Idents(req_subset, cells = select.cells) <- 17
# Idents(req_subset)
DimPlot(req_subset,label = T, pt.size = 1)
table((req_subset@active.ident))
levels(req_subset@active.ident)
req_subset@active.ident <- factor(req_subset@active.ident, levels = as.character(c(0:17)))
DimPlot(req_subset,label = T, pt.size = 1)
req_subset[["Cell.clusters.Fib"]] <- Idents(object = req_subset)

# checking the markers
DotPlot(req_subset, features = markerGenes,group.by = "Cell.clusters.Fib") + RotatedAxis()
# req_subset[["Cell.clusters.Fib"]] <- Idents(object = req_subset)
DimPlot(req_subset,label = T,group.by = "Cell.clusters.Fib")

# getting marker genes
req_subset <- SetIdent(req_subset, value = "Cell.clusters.Fib")
DimPlot(req_subset,label = T)
marker <- FindAllMarkers(req_subset,only.pos = T)
marker <- marker[order(marker$cluster, marker$avg_log2FC, decreasing = T),]
head(marker,20)
tail(marker,20)
# write.csv(marker, file = "results/hum_wLiv/markers_Fib_clusters.csv")


# settling cell types
req_subset <- SetIdent(req_subset, value = "Cell.clusters.Fib")
DimPlot(req_subset,label = T)
req_subset <- RenameIdents(object = req_subset, 
                           '0' = "HSC",
                           '1' = "HSC",
                           '3' = "HSC",
                           '5' = "HSC",
                           '2' = "HSC",
                           '7' = "HSC",
                           '8' = "HSC",
                           '9' = "HSC",
                           
                           
                           '6' = "VSMC",
                           '15' = "Mesothelial",
                           
                           '10' = "Mixed",
                           '11' = "Mixed",
                           '12' = "Mixed",
                           '13' = "Mixed",
                           '14' = "Mixed",
                           '16' = "Mixed",
                           '17' = "Mixed",
                           '4' = "Mixed")






req_subset[["Cell.types.Fib"]] <- Idents(object = req_subset)
DimPlot(req_subset,label = T, group.by = "Cell.types.Fib")
DimPlot(req_subset,label = T, group.by = "Cell.clusters.Fib")

# saveRDS(req_subset, file = "results/hum_wLiv/seuratObj_humLivNucseq_Norm_Nafld_HCC_tumorAdj_mito20_Fib.rds")

req_subset <- readRDS(file = "results/hum_wLiv/seuratObj_humLivNucseq_Norm_Nafld_HCC_tumorAdj_mito20_Fib.rds")
DimPlot(req_subset, label = T)
DimPlot(req_subset, label = T,group.by = "Cell.clusters.Fib")
DimPlot(req_subset,group.by = "orig.ident", label = F,shuffle = T)
DimPlot(req_subset, reduction = "umap",group.by = "Cell.types", label = T)
DimPlot(req_subset, reduction = "umap",group.by = "Cell.types.Fib", label = T)
DimPlot(req_subset, reduction = "umap",group.by = "Cell.clusters", label = T)
DimPlot(req_subset, reduction = "umap",group.by = "Cell.clusters.new", label = T)
DimPlot(req_subset, reduction = "umap",group.by = "DiseaseTypes", label = F,shuffle = T)
DimPlot(req_subset, reduction = "umap",group.by = "DiseaseTypesSimple", label = F,shuffle = T)

# ************************************************** getting the cy-my signature from the mouse dataset
DefaultAssay(req_subset) <- "RNA"
DotPlot(req_subset, features = markerGenes) + RotatedAxis()
library(readxl)
dataset <- data.frame(read_excel("results/cymySig_secreated.xlsx", skip = 1))
head(dataset)

cySig <- na.omit(dataset$cySig2)
mySig <- na.omit(dataset$mySgi2)

# Mouse to human gene conversion using homologus genes
homologs_human_mouse <- readRDS(file="D:/pCloud Sync/Columbia/CommonResource/Data/homologs_mouse_human_geneSymbol.rds")
head(homologs_human_mouse)

mus2hum <- cySig
tmppos <- match(mus2hum, homologs_human_mouse[,"MGI.symbol"])
humGenes <- homologs_human_mouse[tmppos,"HGNC.symbol"]
cysig_hum <- na.omit(humGenes)

mus2hum <- mySig
tmppos <- match(mus2hum, homologs_human_mouse[,"MGI.symbol"])
humGenes <- homologs_human_mouse[tmppos,"HGNC.symbol"]
mySig_hum <- na.omit(humGenes)

ptsize <- 2
sig <- cysig_hum
sigName <- "cySig"
callName <- paste0(sigName,"1")
req_subset <- AddModuleScore(req_subset,features = list(sig), name = sigName)
FeaturePlot(req_subset, callName, max.cutoff = "q99") + scale_color_viridis_c(option = "viridis")
FeaturePlot(req_subset, callName,pt.size = ptsize,order = T) + scale_color_viridis_c(option = "viridis")+labs(title = "Human cy signature")

sig <- mySig_hum
sigName <- "mySig"
callName <- paste0(sigName,"1")
req_subset <- AddModuleScore(req_subset,features = list(sig), name = sigName)
FeaturePlot(req_subset, callName,pt.size = ptsize,order = T) + scale_color_viridis_c(option = "viridis")+labs(title = "Human my signature")

VlnPlot(req_subset, features = c("mySig1", "cySig1"),group.by = "Cell.clusters.Fib", sort = F)

# quantify the cy-my score average in each cluster
req_subset <- SetIdent(req_subset, value = "Cell.clusters.Fib")
df <- data.frame(clust=req_subset@active.ident, cy=req_subset@meta.data$cySig1, my=req_subset@meta.data$mySig1)
head(df)
ggplot(df, aes(x=clust, y=cy))+geom_violin()

library(tidyverse)
# function to compute the median cy-my score  in each cluster
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(median = median(x[[col]], na.rm=TRUE) )
  } 
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("median" = varname))
  return(data_sum)
}

tmp <- data_summary(df, varname = "cy", groupnames = "clust")
tmp
cymyscore <- tmp[,c(1,2)]
tmp <- data_summary(df, varname = "my", groupnames = "clust")
tmp
cymyscore <- cbind(cymyscore, my=tmp[,c(2)], cy_my=(cymyscore[,2]-tmp[,c(2)]))
cymyscore

req_subset <- SetIdent(req_subset, value = "Cell.clusters.Fib")
DimPlot(req_subset,label = T)
req_subset <- RenameIdents(object = req_subset, 
                           '0' = "myHSC",
                           '1' = "cyHSC",
                           '3' = "cyHSC",
                           '5' = "cyHSC",
                           '2' = "myHSC",
                           '7' = "myHSC",
                           '8' = "myHSC",
                           '9' = "myHSC",
                           
                           '6' = "VSMC",
                           '15' = "Mesothelial",
                           
                           '10' = "Mixed",
                           '11' = "Mixed",
                           '12' = "Mixed",
                           '13' = "Mixed",
                           '14' = "Mixed",
                           '16' = "Mixed",
                           '17' = "Mixed",
                           '4' = "Mixed")

req_subset[["Cell.types.Fib.cymy"]] <- Idents(object = req_subset)
DimPlot(req_subset,label = T, group.by = "Cell.types.Fib.cymy")

VlnPlot(req_subset, features = c("mySig1", "cySig1"),group.by = "Cell.types.Fib.cymy", sort = F)
DotPlot(req_subset, features = markerGenes,group.by = "Cell.types.Fib.cymy") + RotatedAxis()

# saveRDS(req_subset, file = "results/hum_wLiv/seuratObj_humLivNucseq_Norm_Nafld_HCC_tumorAdj_mito20_Fib.rds")

req_subset <- readRDS(file = "results/hum_wLiv/seuratObj_humLivNucseq_Norm_Nafld_HCC_tumorAdj_mito20_Fib.rds")



# **************************************************************************************
# # ****************************** subset the HSC population 
req_subset <- readRDS(file = "results/hum_wLiv/seuratObj_humLivNucseq_Norm_Nafld_HCC_tumorAdj_mito20_Fib.rds")
req_subset <- SetIdent(req_subset, value = "Cell.types.Fib")
DimPlot(req_subset,label = T)
req_subset_HSC <- subset(req_subset, idents = c("HSC"))
DimPlot(req_subset_HSC,label = T, group.by = "ident")
DimPlot(req_subset_HSC,label = T, group.by = "Cell.types.Fib.cymy") +xlim(-9,5) #+ylim(-4,8)
DimPlot(req_subset_HSC,label = T, group.by = "Cell.clusters.Fib") +xlim(-9,5)

# cy-my cluster selection
sig <- toupper(unlist(strsplit("Acta2, Col1a1, Lox, Timp1",split = ", ")))
sigName <- "Myofib"
callName <- paste0(sigName,"1")
req_subset_HSC <- AddModuleScore(req_subset_HSC,features = list(sig), name = sigName)
FeaturePlot(req_subset_HSC, callName,pt.size = 2,order = T) + scale_color_viridis_c(option = "viridis")+labs(title = paste0("Human ",sigName," signature"))

VlnPlot(req_subset_HSC, features = c("Myofib1"),group.by = "Cell.clusters.Fib", sort = T)

req_subset_HSC <- readRDS(file = "results/hum_wLiv/seuratObj_humLivNucseq_Norm_Nafld_HCC_tumorAdj_mito20_Fib_HSCsubset.rds")

DimPlot(req_subset_HSC,label = T, pt.size = 3)  +xlim(-9,5)



# **************************************************************************************
# ******************************get the fibroblast annotations back into whole dataset
req_subset <- readRDS(file = "results/hum_wLiv/seuratObj_humLivNucseq_Norm_Nafld_HCC_tumorAdj_mito20_Fib.rds")
DimPlot(req_subset,label = T)
table(req_subset@meta.data$Cell.types.Fib)
data.combined_all <- readRDS(file = "results/hum_wLiv/seuratObj_humLivNucseq_Norm_Nafld_HCC_tumorAdj_mito20.rds")
DimPlot(data.combined_all, label = T)
table(data.combined_all@meta.data$Cell.types)


# writing the cell types annotation
df <- data.frame(cells=row.names(data.combined_all@meta.data), samples=data.combined_all@meta.data$orig.ident, clusters=data.combined_all@meta.data$Cell.clusters.new, Cell.types=data.combined_all@meta.data$Cell.types, DiseaseTypes=data.combined_all@meta.data$DiseaseTypes, DiseaseTypesSimple=data.combined_all@meta.data$DiseaseTypesSimple)
head(df)

df.fib <- data.frame(cells=row.names(req_subset@meta.data),clusters.Fib=req_subset@meta.data$Cell.clusters.Fib,Cell.types.Fib.subtypes=req_subset@meta.data$Cell.types.Fib, Cell.types.Fib.cymy=req_subset@meta.data$Cell.types.Fib.cymy)
head(df.fib)
table(df.fib$Cell.types.Fib)
table((df.fib$Cell.types.Fib.cymy))

df$clusters.Fib <- as.character(df$clusters)
df$Cell.types.Fib <-as.character(df$Cell.types)
df$Cell.types.Fib.subtypes <-as.character(df$Cell.types)
df$Cell.types.Fib.cymy <- as.character(df$Cell.types)
head(df)

tmppos <- match(df.fib$cells, df$cells)
head(df.fib)
head(df[tmppos,])
df$clusters.Fib[tmppos] <- paste0("fib_",df.fib$clusters.Fib)
df$Cell.types.Fib.subtypes[tmppos] <- as.character(df.fib$Cell.types.Fib.subtypes)
df$Cell.types.Fib.cymy[tmppos] <- as.character(df.fib$Cell.types.Fib.cymy)

# get the fibroblast mixed back to w-liv annotation
tmppos <- which(df$Cell.types.Fib.subtypes=="Mixed")
table(df$Cell.types.Fib.subtypes)
table(df$Cell.types.Fib.subtypes[tmppos])
table(df$Cell.types.Fib)
table(df$Cell.types.Fib[tmppos])
df$Cell.types.Fib[tmppos] <- df$Cell.types.Fib.subtypes[tmppos]
table(df$Cell.types.Fib)
table(df$Cell.types.Fib[tmppos])
head(df)
head(df[tmppos,])


# write.csv(df, file = "results/hum_wLiv/humLivNucseq_Norm_Nafld_HCC_tumorAdj_annotation.csv")

df <- read.csv(file = "results/hum_wLiv/humLivNucseq_Norm_Nafld_HCC_tumorAdj_annotation.csv", header = T, row.names = 1)
head(df)
table(df$samples)
table((df$clusters))
table((df$Cell.types))
table(df$DiseaseTypes)
table((df$DiseaseTypesSimple))
table(df$Cell.types.Fib)
table(df$Cell.types.Fib.subtypes)
table((df$Cell.types.Fib.cymy))


data.combined_all[["Cell.types.Fib"]] <- df$Cell.types.Fib
DimPlot(data.combined_all, label = T,group.by = "Cell.types.Fib")

data.combined_all[["Cell.types.Fib.subtypes"]] <- df$Cell.types.Fib.subtypes
DimPlot(data.combined_all, label = T,group.by = "Cell.types.Fib.subtypes")

data.combined_all[["Cell.types.Fib.cymy"]] <- df$Cell.types.Fib.cymy
DimPlot(data.combined_all, label = T,group.by = "Cell.types.Fib.cymy")


# quantify the Fibroblast cells in tumor and adj.tumor
table(df$samples)
tmppos <- which(df$samples =="S41")
df_sample_adj <- df[tmppos,]
table(df_sample_adj$DiseaseTypes)
table(df_sample_adj$Cell.types)
table(df_sample_adj$Cell.types.Fib)
table(df_sample_adj$Cell.types.Fib.cymy)

tmppos <- which(df$samples =="S42")
df_sample_tumor <- df[tmppos,]
table(df_sample_tumor$DiseaseTypes)
table(df_sample_tumor$Cell.types)
table(df_sample_tumor$Cell.types.Fib)
table(df_sample_tumor$Cell.types.Fib.cymy)

table(df$samples)
tmppos <- which(df$samples =="S43")
df_sample_adj <- df[tmppos,]
table(df_sample_adj$DiseaseTypes)
table(df_sample_adj$Cell.types)
table(df_sample_adj$Cell.types.Fib)
table(df_sample_adj$Cell.types.Fib.cymy)

tmppos <- which(df$samples =="S44")
df_sample_tumor <- df[tmppos,]
table(df_sample_tumor$DiseaseTypes)
table(df_sample_tumor$Cell.types)
table(df_sample_tumor$Cell.types.Fib)
table(df_sample_tumor$Cell.types.Fib.cymy)



# saveRDS(data.combined_all, file = "results/hum_wLiv/seuratObj_humLivNucseq_Norm_Nafld_HCC_tumorAdj_mito20.rds")


# *************************************************************
# get the required input data for cellphoneDB
# *************************************************************
data.combined_all <- readRDS(file = "results/hum_wLiv/seuratObj_humLivNucseq_Norm_Nafld_HCC_tumorAdj_mito20.rds")
DimPlot(data.combined_all, label = T)
DimPlot(data.combined_all, reduction = "umap",group.by = "DiseaseTypesSimple", label = F,shuffle = T)
data.combined_all <- SetIdent(data.combined_all, value = "DiseaseTypesSimple")
DimPlot(data.combined_all, label = T)
req_subset <- subset(data.combined_all, idents = c("Fibrotic"),invert=F)
req_subset
DimPlot(req_subset, label = T)
req_subset <- SetIdent(req_subset, value = "Cell.types.Fib.subtypes")
DimPlot(req_subset, label = T)
req_subset <- subset(req_subset, idents = c("Mixed"),invert=T)
DimPlot(req_subset, label = T)
DimPlot(req_subset, label = F)


data_req <- req_subset
DimPlot(data_req,label = T)#, group.by = "Cell.types.Fib")

dim(data_req@assays$RNA@counts)
dim(data_req@assays$RNA@data)
dim(data_req@assays$integrated@data)

table(data_req@meta.data$Cell.types.Fib.subtypes)
table(data_req@meta.data$Cell.clusters)
DimPlot(data_req, label = T, group.by = "Cell.clusters")
data_req <- SetIdent(data_req, value = "Cell.clusters")
DimPlot(data_req, label = T)
# Downsample the number of cells per identity class
data_req <- subset(x = data_req, downsample = 400)
dim(data_req@assays$RNA@counts)

table(data_req@meta.data$Cell.types.Fib.subtypes)
DimPlot(data_req, label = T)
DimPlot(data_req, label = T,group.by = "DiseaseTypesSimple")
DimPlot(data_req, label = T,group.by = "Cell.types.Fib.subtypes")


# getting the normalized counts
rawCounts <- as.matrix(data_req@assays$RNA@counts)
rawCounts[1:5,1:5]
dim(rawCounts)

# getting the cell types
cellLabels <- data.frame(cells=colnames(data_req@assays$RNA@counts), cellType=data_req@meta.data$Cell.types.Fib.subtypes)
head(cellLabels)
dim(cellLabels)
table(cellLabels$cellType)
all(colnames(rawCounts)==as.character(cellLabels[,"cells"]))

#normalizing for library depth
temp <- rawCounts 
dim(temp)
temp[1:5,1:5]
temp_norm_mat <- apply(temp, 2, function(x) (x/sum(x))*10000)
colSums(temp_norm_mat)[1:10]
rowSums(temp_norm_mat)[1:10]
remove(temp)

# normazlied counts data
cpdbNormData <- temp_norm_mat
dim(cpdbNormData)
# write.table(t(c("Genes", colnames(cpdbNormData))),file = "cpdb/hum_wLiv/data_cellphoneDB_countNorm_hum_wLiv.txt", quote = FALSE, sep = "\t", append = FALSE, col.names = FALSE, row.names = FALSE)
# write.table(cpdbNormData,file = "cpdb/hum_wLiv/data_cellphoneDB_countNorm_hum_wLiv.txt", quote = FALSE, sep = "\t", append = TRUE, col.names = FALSE)



# #getting the metadata; all cell types
cellLabels_req <- cellLabels 
metaData <- data.frame(Cell=cellLabels_req[,"cells"],cell_type=cellLabels_req[,"cellType"])
all(metaData[,"Cell"]==colnames(cpdbNormData))
head(metaData)
dim(metaData)
table(metaData[,"cell_type"])
# write.table(metaData, file = "cpdb/hum_wLiv/data_cellphoneDB_metaData_hum_wLiv_cellTypes.txt", sep = '\t', quote = F, row.names = F, col.names = T)



# *************************************************************
# using saved data
# *************************************************************
data.combined_all <- readRDS(file = "results/hum_wLiv/seuratObj_humLivNucseq_Norm_Nafld_HCC_tumorAdj_mito20.rds")
DimPlot(data.combined_all, label = T)
DimPlot(data.combined_all, label = T,group.by = "Cell.clusters")
DimPlot(data.combined_all,group.by = "orig.ident", label = F,shuffle = T)
DimPlot(data.combined_all, reduction = "umap",group.by = "Cell.types", label = T)
# DimPlot(data.combined_all, reduction = "umap",group.by = "Cell.types.new", label = T)
DimPlot(data.combined_all, reduction = "umap",group.by = "Cell.clusters", label = T)
DimPlot(data.combined_all, reduction = "umap",group.by = "Cell.clusters.new", label = T)
DimPlot(data.combined_all, reduction = "umap",group.by = "DiseaseTypes", label = F,shuffle = T)
DimPlot(data.combined_all, reduction = "umap",group.by = "DiseaseTypesSimple", label = F,shuffle = T)

DefaultAssay(data.combined_all) <- "RNA"
DotPlot(data.combined_all, features = markerGenes) + RotatedAxis()
DotPlot(data.combined_all, features = markerGenes, cluster.idents = T) + RotatedAxis()

DotPlot(data.combined_all, features = markerGenes,group.by = "Cell.clusters.new") + RotatedAxis()
DotPlot(data.combined_all, features = markerGenes, cluster.idents = T) + RotatedAxis()

# Fibroblasts
req_subset <- readRDS(file = "results/hum_wLiv/seuratObj_humLivNucseq_Norm_Nafld_HCC_tumorAdj_mito20_Fib.rds")
DimPlot(req_subset, label = T)
DimPlot(req_subset, label = T,group.by = "Cell.clusters")
DimPlot(req_subset,group.by = "orig.ident", label = F,shuffle = T)
DimPlot(req_subset, reduction = "umap",group.by = "Cell.types", label = T)
# DimPlot(req_subset, reduction = "umap",group.by = "Cell.types.new", label = T)
DimPlot(req_subset, reduction = "umap",group.by = "Cell.clusters", label = T)
DimPlot(req_subset, reduction = "umap",group.by = "Cell.clusters.new", label = T)
DimPlot(req_subset, reduction = "umap",group.by = "Cell.clusters.Fib", label = T)
DimPlot(req_subset, reduction = "umap",group.by = "Cell.types.Fib", label = T)
DimPlot(req_subset, reduction = "umap",group.by = "Cell.types.Fib.cymy", label = T)
DimPlot(req_subset, reduction = "umap",group.by = "DiseaseTypes", label = F,shuffle = T)
DimPlot(req_subset, reduction = "umap",group.by = "DiseaseTypesSimple", label = F,shuffle = T)

# HSC
req_subset_HSC <- readRDS(file = "results/hum_wLiv/seuratObj_humLivNucseq_Norm_Nafld_HCC_tumorAdj_mito20_Fib_HSCsubset.rds")

DimPlot(req_subset_HSC,label = T, pt.size = 3)  +xlim(-9,5)
DimPlot(req_subset_HSC,label = T, group.by = "orig.ident") +xlim(-9,5)
DimPlot(req_subset_HSC,label = T, group.by = "Cell.types.Fib.cymy") +xlim(-9,5)
DimPlot(req_subset_HSC,label = T, group.by = "Cell.clusters.Fib") +xlim(-9,5)




#******************************************************#******************************************************
#******************************************************PLOTS
#******************************************************#******************************************************
library(RColorBrewer)
library("ggsci")
library("scales")
library(ggplot2)

data.combined_all <- readRDS(file = "results/hum_wLiv/seuratObj_humLivNucseq_Norm_Nafld_HCC_tumorAdj_mito20.rds")
DimPlot(data.combined_all,label = T, group.by = "Cell.clusters.new")
DefaultAssay(data.combined_all) <- "RNA"
# ggsave(filename = paste0("results/hum_wLiv/figs//hum_wLiv_cellClusters_label.pdf"), width = 8, height = 8)

DimPlot(object = data.combined_all, group.by = "orig.ident", label = F, shuffle = T, pt.size = 1)  + labs(title = "Human whole-liver")+  scale_color_aaas()
ggsave(filename = paste0("results/hum_wLiv/figs//hum_wLiv_samples.pdf"), width = 8, height = 8)
DimPlot(object = data.combined_all, group.by = "DiseaseTypes", label = F, shuffle = T, pt.size = 1)  + labs(title = "Human whole-liver")+  scale_color_aaas()
ggsave(filename = paste0("results/hum_wLiv/figs//hum_wLiv_disease.pdf"), width = 8, height = 8)
DimPlot(object = data.combined_all, group.by = "DiseaseTypesSimple", label = F, shuffle = T, pt.size = 1)  + labs(title = "Human whole-liver")+  scale_color_aaas()
ggsave(filename = paste0("results/hum_wLiv/figs//hum_wLiv_diseaseSimple.pdf"), width = 8, height = 8)

markerGenes <- c("ALB","APOE","CPS1",
                 "IL7R","THEMIS","NKG7", #,"Trac" #TCells
                 "CD163","FCGR2A", #"Trem2", #myeloid
                 'ANXA4', 'CFTR',
                 "HGF","COLEC11",#"COL1A1", #"Col1a1", #"Lox","Col3a1",
                 "DCN",#"PDGFRA","PDGFRB",
                 "KDR","ENG",
                 "MS4A1","BLK",
                 "JCHAIN","IGKC",#,"Jchain","Siglech",
                 "TOP2A","MKI67")

DotPlot(data.combined_all, features = markerGenes, group.by = "Cell.clusters.new", cluster.idents = F) + RotatedAxis()+ scale_color_viridis_c() #, cols = "RdBu"
ggsave(filename = paste0("results/hum_wLiv/figs//hum_wLiv_all_clusters_dotplot.pdf"))#, width = 6, height = 6)


# get only the clean clusters
data.combined_all <- SetIdent(data.combined_all, value = "Cell.types")
DimPlot(data.combined_all,label = T)
req_subset <- subset(data.combined_all, idents = c("Mixed"), invert=T)
req_subset
DimPlot(req_subset,label = T, group.by = "ident",pt.size = 1)
DimPlot(req_subset,label = T, group.by = "Cell.types")
DimPlot(req_subset,label = T, group.by = "Cell.clusters.new")
# DimPlot(req_subset,label = T, group.by = "Cell.types.cymy")


# CELL TYPES UMAP
data_req <- req_subset
remove(req_subset)
DimPlot(object = data_req, label = T)
levels(data_req@active.ident)


cols <- c('#00468BFF','#CC99FFFF','#9900CCFF','#1B1919FF','#CCFF00FF','#FF7F00','#42B540FF','#B2df8a','#AD002AFF')
DimPlot(object = data_req, group.by = "ident", label = F, cols = cols, pt.size = 1) 

DimPlot(object = data_req, group.by = "Cell.types", label = T, cols = cols, pt.size = 1)  + labs(title = "Human whole-liver")
ggsave(filename = paste0("results/hum_wLiv/figs//hum_wLiv_celltypes_label.pdf"), width = 8, height = 8)
DimPlot(object = data_req, group.by = "Cell.types", label = F, cols = cols, pt.size = 1)  + labs(title = "Human whole-liver")+ NoLegend()
ggsave(filename = paste0("results/hum_wLiv/figs//hum_wLiv_celltypes_noLabel.pdf"), width = 8, height = 8)


genes <- c(unlist(strsplit("COL1A1, HGF",split = ", "))) #Hgf, Col1a1, 
# markers required
ptsize <- 1
genes <- intersect(genes, rownames(data_req@assays$RNA@data))
FeaturePlot(object = data_req, features = c(genes[2]), max.cutoff = "q99", pt.size = ptsize,order = F) + scale_color_viridis_c()

for(gene in genes){
  FeaturePlot(object = data_req, features = c(gene), max.cutoff = "q99", pt.size = ptsize,order = F) + scale_color_viridis_c()
  ggsave(filename = paste0("results/hum_wLiv/figs/hum_wLiv_gene_",gene,".pdf"), width = 8, height = 8)
}

# violin plots
data_req <- SetIdent(data_req, value = "Cell.types.Fib.subtypes")
DimPlot(data_req,label = T)
gene <- "HGF"
VlnPlot(data_req, features = gene, idents = c("HSC","Endothelial")) #+geom_boxplot()
ggsave(filename = paste0("results/hum_wLiv/figs/hum_wLiv_vln_gene_",gene,".pdf"), width = 8, height = 8)



# plots for subsets of samples and cell types
ptsize <- 2
celltype <- "Hepatocytes"
DimPlot(data_req,label = T)
tmp <- subset(data_req, idents = celltype) 
DimPlot(tmp,label = T)+xlim(-12,8)+ylim(-8,6)
tmp <- SetIdent(tmp, value = "DiseaseTypesSimple")
DimPlot(tmp,label = T)+xlim(-12,8)+ylim(-8,6)


gene <- "DDR1"
sampleSel <- c("Normal")
# DimPlot(tmp, reduction = "umap", cells = unlist(CellsByIdentities(tmp, idents = sampleSel)))+xlim(-12,8)+ylim(-8,6)
# FeaturePlot(tmp,features = gene,cells = unlist(CellsByIdentities(tmp, idents = sampleSel)))+xlim(-12,8)+ylim(-8,6)
FeaturePlot(object = tmp, features = c(gene), max.cutoff = "q99", pt.size = ptsize,order = T,cells = unlist(CellsByIdentities(tmp, idents = sampleSel)))+xlim(-12,8)+ylim(-8,6) + scale_color_viridis_c()
ggsave(filename = paste0("results/hum_wLiv/figs/hum_wLiv_",sampleSel,"_",celltype,"_gene_",gene,".pdf"), width = 8, height = 8)

sampleSel <- c("Fibrotic")
FeaturePlot(object = tmp, features = c(gene), max.cutoff = "q99", pt.size = ptsize,order = T,cells = unlist(CellsByIdentities(tmp, idents = sampleSel)))+xlim(-12,8)+ylim(-8,6) + scale_color_viridis_c()
ggsave(filename = paste0("results/hum_wLiv/figs/hum_wLiv_",sampleSel,"_",celltype,"_gene_",gene,".pdf"), width = 8, height = 8)


sampleSel <- c("Tumor")
FeaturePlot(object = tmp, features = c(gene), max.cutoff = "q99", pt.size = ptsize,order = T,cells = unlist(CellsByIdentities(tmp, idents = sampleSel)))+xlim(-12,8)+ylim(-8,6) + scale_color_viridis_c()
ggsave(filename = paste0("results/hum_wLiv/figs/hum_wLiv_",sampleSel,"_",celltype,"_gene_",gene,".pdf"), width = 8, height = 8)


VlnPlot(tmp, features = gene)
ggsave(filename = paste0("results/hum_wLiv/figs/hum_wLiv_vln_",celltype,"_gene_",gene,".pdf"), width = 8, height = 8)


data_req <- SetIdent(data_req, value = "DiseaseTypes")
DimPlot(data_req,label = T)
tmp <- subset(data_req, idents = c("Adj.tumor","Tumor")) 
DimPlot(tmp,label = T)#+xlim(-12,8)+ylim(-8,6)

sampleSel <- c("Adj.tumor")
DimPlot(object = tmp, group.by = "Cell.types", label = T, cols = cols, pt.size = 1,cells = unlist(CellsByIdentities(tmp, idents = sampleSel)))  + labs(title = "Human whole-liver")
ggsave(filename = paste0("results/hum_wLiv/figs//hum_wLiv_",sampleSel,"_celltypes_label.pdf"), width = 8, height = 8)
DimPlot(object = tmp, group.by = "Cell.types", label = F, cols = cols, pt.size = 1,cells = unlist(CellsByIdentities(tmp, idents = sampleSel)))  + labs(title = "Human whole-liver")+ NoLegend()
ggsave(filename = paste0("results/hum_wLiv/figs//hum_wLiv_",sampleSel,"_celltypes_noLabel.pdf"), width = 8, height = 8)


sampleSel <- c("Tumor")
DimPlot(object = tmp, group.by = "Cell.types", label = T, cols = cols, pt.size = 1,cells = unlist(CellsByIdentities(tmp, idents = sampleSel)))  + labs(title = "Human whole-liver")
ggsave(filename = paste0("results/hum_wLiv/figs//hum_wLiv_",sampleSel,"_celltypes_label.pdf"), width = 8, height = 8)
DimPlot(object = tmp, group.by = "Cell.types", label = F, cols = cols, pt.size = 1,cells = unlist(CellsByIdentities(tmp, idents = sampleSel)))  + labs(title = "Human whole-liver")+ NoLegend()
ggsave(filename = paste0("results/hum_wLiv/figs//hum_wLiv_",sampleSel,"_celltypes_noLabel.pdf"), width = 8, height = 8)


# ********************************************************** plots for Fib part
req_subset <- readRDS(file = "results/hum_wLiv/seuratObj_humLivNucseq_Norm_Nafld_HCC_tumorAdj_mito20_Fib.rds")
req_subset <- SetIdent(req_subset, value = "Cell.types.Fib")
DimPlot(req_subset,label = T)
DefaultAssay(req_subset) <- "RNA"

data_req <- req_subset
DefaultAssay(data_req) <- "RNA"
levels(data_req@active.ident)
ptsize <- 3

# samples present
DimPlot(object = data_req, group.by = "ident", label = F, shuffle = T, pt.size = ptsize) + labs(title = "Human whole-liver Fib")+  scale_color_aaas()
ggsave(filename = paste0("results/hum_wLiv/figs//hum_wLiv_Fib_celltypes.pdf"), width = 8, height = 8)



# ********************************************************** plots for HSC part
req_subset_HSC <- readRDS(file = "results/hum_wLiv/seuratObj_humLivNucseq_Norm_Nafld_HCC_tumorAdj_mito20_Fib_HSCsubset.rds")
DimPlot(req_subset_HSC,label = T, pt.size = 3)  +xlim(-9,5)
DefaultAssay(req_subset_HSC) <- "RNA"

data_req <- req_subset_HSC
DefaultAssay(data_req) <- "RNA"
data_req@meta.data$Cell.types.Fib.cymy <- factor(data_req@meta.data$Cell.types.Fib.cymy, levels<-c("cyHSC","myHSC"))
data_req <- SetIdent(data_req,value = "Cell.types.Fib.cymy")
DimPlot(object = data_req, label = T)+xlim(-9,5)
levels(data_req@active.ident)
ordered(data_req@active.ident)
ptsize <- 3

# removing tumor samples
table((data_req@meta.data$orig.ident))
data_req <- SetIdent(data_req,value = "orig.ident")
DimPlot(object = data_req, label = T)+xlim(-9,5)
data_req <- subset(data_req, idents = c("S42","S44"),invert=T)
DimPlot(object = data_req, label = T)+xlim(-9,5)
DimPlot(object = data_req, label = T,group.by = "DiseaseTypes")+xlim(-9,5)

# samples present
DimPlot(object = data_req, group.by = "orig.ident", label = F, shuffle = T, pt.size = ptsize) +xlim(-9,5) + labs(title = "Human whole-liver HSC")+  scale_color_aaas()
ggsave(filename = paste0("results/hum_wLiv/figs//hum_wLiv_HSC_samples.pdf"), width = 8, height = 8)
DimPlot(object = data_req, group.by = "DiseaseTypes", label = F, shuffle = T, pt.size = ptsize) +xlim(-9,5) + labs(title = "Human whole-liver HSC")+  scale_color_aaas()
ggsave(filename = paste0("results/hum_wLiv/figs//hum_wLiv_HSC_disease.pdf"), width = 8, height = 8)
DimPlot(object = data_req, group.by = "DiseaseTypesSimple", label = F, shuffle = T, pt.size = ptsize) +xlim(-9,5) + labs(title = "Human whole-liver HSC")+  scale_color_aaas()
ggsave(filename = paste0("results/hum_wLiv/figs//hum_wLiv_HSC_diseaseSimple.pdf"), width = 8, height = 8)


# cy-my signature
FeaturePlot(data_req, "cySig1", max.cutoff = "q99",pt.size = ptsize,order = T) +xlim(-9,5)+ scale_color_viridis_c(option = "viridis") +labs(title = "Human whole-liver HSC cy signature")
ggsave(filename = paste0("results/hum_wLiv/figs/hum_wLiv_HSC_signature_cy.pdf"), width = 8, height = 8)

FeaturePlot(data_req, "mySig1", max.cutoff = "q99",pt.size = ptsize,order = T) +xlim(-9,5)+ scale_color_viridis_c(option = "viridis") +labs(title = "Human whole-liver HSC my signature")
ggsave(filename = paste0("results/hum_wLiv/figs/hum_wLiv_HSC_signature_my.pdf"), width = 8, height = 8)

cols_cymy <- c("#FF9900","#000099")
DimPlot(data_req, reduction = "umap",group.by = "Cell.types.Fib.cymy", pt.size = ptsize,cols = cols_cymy, label = F)+xlim(-9,5) +labs(title = "Human whole-liver HSC cy-my clusters")
ggsave(filename = paste0("results/hum_wLiv/figs/hum_wLiv_HSC_signature_clust_cymy.pdf"), width = 8, height = 8)


# getting the cy-my for normal and fibrotic
data_req <- SetIdent(data_req, value = "DiseaseTypesSimple")
DimPlot(data_req,label = T) +xlim(-9,5)
tmpObj <- data_req #subset(data_req, idents = c("Normal","Fibrotic")) 
DimPlot(tmpObj,label = T)+xlim(-9,5)

sampleSel <- c("Normal")
FeaturePlot(tmpObj, "cySig1", max.cutoff = "q99",pt.size = ptsize,order = T,cells = unlist(CellsByIdentities(tmpObj, idents = sampleSel))) +xlim(-9,5)+ scale_color_viridis_c(option = "viridis") +labs(title = "Human whole-liver HSC Normal cy signature")
ggsave(filename = paste0("results/hum_wLiv/figs//hum_wLiv_HSC_signature_cy_",sampleSel,".pdf"), width = 8, height = 8)
FeaturePlot(tmpObj, "mySig1", max.cutoff = "q99",pt.size = ptsize,order = T,cells = unlist(CellsByIdentities(tmpObj, idents = sampleSel))) +xlim(-9,5)+ scale_color_viridis_c(option = "viridis") +labs(title = "Human whole-liver HSC Normal my signature")
ggsave(filename = paste0("results/hum_wLiv/figs//hum_wLiv_HSC_signature_my_",sampleSel,".pdf"), width = 8, height = 8)


sampleSel <- c("Fibrotic")
FeaturePlot(tmpObj, "cySig1", max.cutoff = "q99",pt.size = ptsize,order = T,cells = unlist(CellsByIdentities(tmpObj, idents = sampleSel))) +xlim(-9,5)+ scale_color_viridis_c(option = "viridis") +labs(title = "Human whole-liver HSC Fibrotic cy signature")
ggsave(filename = paste0("results/hum_wLiv/figs//hum_wLiv_HSC_signature_cy_",sampleSel,".pdf"), width = 8, height = 8)
FeaturePlot(tmpObj, "mySig1", max.cutoff = "q99",pt.size = ptsize,order = T,cells = unlist(CellsByIdentities(tmpObj, idents = sampleSel))) +xlim(-9,5)+ scale_color_viridis_c(option = "viridis") +labs(title = "Human whole-liver HSC Fibrotic my signature")
ggsave(filename = paste0("results/hum_wLiv/figs//hum_wLiv_HSC_signature_my_",sampleSel,".pdf"), width = 8, height = 8)


cols <- c("red","blue")
gene <-"cySig1"
p <- VlnPlot(tmpObj, features = gene) 
df <- p$data
colnames(df)[1] <- "gene"
head(df)
tmppos <- which(df$ident=="Normal")
tmpx <- df$gene[tmppos]
tmpy <- df$gene[-tmppos]
tmp <- wilcox.test(x=tmpx, y=tmpy)
tmp$p.value
ggplot(data = df, aes(x=ident, y=gene)) +
  geom_violin(aes(fill=ident),scale = "width") + 
  geom_boxplot(width=0.1) +scale_fill_manual(values = cols) +
  theme_classic()+
  ylab(gene)+
  ggtitle(paste0("Human whole-liver HSC ",gene," MW p ",scientific(tmp$p.value,digits = 3)))

# VlnPlot(tmpObj, features = "cySig1",cols = cols,pt.size = 0) +
#   geom_boxplot(width=0.1, fill="white") + #scale_fill_manual(values = c("white","white"))  +
#   labs(title = "Human whole-liver HSC cy signature")
ggsave(filename = paste0("results/hum_wLiv/figs//hum_wLiv_HSC_signature_cy_vln.pdf"))

gene <-"mySig1"
p <- VlnPlot(tmpObj, features = gene) 
df <- p$data
colnames(df)[1] <- "gene"
head(df)
tmppos <- which(df$ident=="Normal")
tmpx <- df$gene[tmppos]
tmpy <- df$gene[-tmppos]
tmp <- wilcox.test(x=tmpx, y=tmpy)
tmp$p.value
ggplot(data = df, aes(x=ident, y=gene)) +
  geom_violin(aes(fill=ident),scale = "width") + 
  geom_boxplot(width=0.1) +scale_fill_manual(values = cols) +
  theme_classic()+
  ylab(gene)+
  ggtitle(paste0("Human whole-liver HSC ",gene," MW p ",scientific(tmp$p.value,digits = 3)))
ggsave(filename = paste0("results/hum_wLiv/figs//hum_wLiv_HSC_signature_my_vln.pdf"))



# markers required
genes <- c(unlist(strsplit("COL1A1, HGF",split = ", "))) #Hgf, Col1a1, 
genes <- intersect(genes, rownames(data_req@assays$RNA@data))
FeaturePlot(object = data_req, features = c(genes[2]), max.cutoff = "q99", pt.size = ptsize,order = F) +xlim(-9,5)+ scale_color_viridis_c()

for(gene in genes){
  FeaturePlot(object = data_req, features = c(gene), max.cutoff = "q99", pt.size = ptsize,order = F) +xlim(-9,5)+ scale_color_viridis_c()
  ggsave(filename = paste0("results/hum_wLiv/figs/hum_wLiv_HSC_gene_",gene,".pdf"), width = 8, height = 8)
}

# violin plots
gene <- "COL1A1"
p <- VlnPlot(data_req, features = gene, group.by = "Cell.types.Fib.cymy", cols = cols_cymy) #+geom_boxplot()
df <- p$data
colnames(df)[1] <- "gene"
head(df)
tmppos <- which(df$ident=="cyHSC")
tmpx <- df$gene[tmppos]
tmpy <- df$gene[-tmppos]
tmp <- wilcox.test(x=tmpx, y=tmpy)
tmp$p.value

ggplot(data = df, aes(x=ident, y=gene)) +
  geom_violin(aes(fill=ident),scale = "width") + 
  geom_boxplot(width=0.1) +scale_fill_manual(values = cols_cymy) +
  theme_classic()+
  ylab(gene)+
  ggtitle(paste0("Hum ",gene," Mann Whitney p ",scientific(tmp$p.value,digits = 3)))
ggsave(filename = paste0("results/hum_wLiv/figs/hum_wLiv_HSC_signature_clust_cymy_vln_",gene,".pdf"),width = 5)


gene <- "HGF" 
p <- VlnPlot(data_req, features = gene, group.by = "Cell.types.Fib.cymy", cols = cols_cymy) #+geom_boxplot()
df <- p$data
colnames(df)[1] <- "gene"
head(df)
tmppos <- which(df$ident=="cyHSC")
tmpx <- df$gene[tmppos]
tmpy <- df$gene[-tmppos]
tmp <- wilcox.test(x=tmpx, y=tmpy)
tmp$p.value

ggplot(data = df, aes(x=ident, y=gene)) +
  geom_violin(aes(fill=ident),scale = "width") +
  geom_boxplot(width=0.1) +scale_fill_manual(values = cols_cymy) +
  theme_classic()+
  ylab(gene)+
  ggtitle(paste0("Hum ",gene," Mann Whitney p ",scientific(tmp$p.value,digits = 3)))
ggsave(filename = paste0("results/hum_wLiv/figs/hum_wLiv_HSC_signature_clust_cymy_vln_",gene,".pdf"),width = 5)

# get the cy my correlation plots
library(ggpubr)
df <- data.frame(cyHSC=data_req@meta.data$cySig1, myHSC=data_req@meta.data$mySig1, score=data_req@meta.data$mySig1-data_req@meta.data$cySig1)
head(df)
ggscatter(df, x = "cyHSC", y = "myHSC",
          add = "reg.line",                                 # Add regression line
          conf.int = TRUE,                                  # Add confidence interval
          add.params = list(color = "black", fill = "lightgray"),
          color = "score")+
  stat_cor(method = "pearson", label.x = 0.35, label.y = .7)+  # Add correlation coefficient
  gradient_color(c("blue", "yellow", "red"))
ggsave(filename = paste0("results/hum_wLiv/figs/hum_wLiv_HSC_signature_clust_cymy_correln.pdf"))

