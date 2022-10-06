# 2021 Nov 01 created
# this is the read and analysis of the mus wLiver CCl4 datasets


setwd("D:/pCloud Sync/Columbia/SchwabeLab/Aveline_HSCnHCC/github/HSCinHCC/")

library(Seurat)
library(ggplot2)
library(Matrix)
source("D:/pCloud Sync/Columbia/CommonResource/Codes/commonFunctions.R")


markerGenes <- c("Myh11","Actg2", #VSMC
                 "Rgs5", #Pericytes
                 "Lrat","Dcn","Colec11","Pdgfra","Pdgfrb", #HSC
                 "Acta2","Lox","Col3a1","Col1a1", #Act-HSC
                 "Mgp","Eln","Mfap4","Dpt","Gsn", #PF markers
                 "Kdr","Aqp1", #Endothelial
                 "Top2a","Mki67", #Cycling
                 'Upk3b','Msln', #Mesothelial
                 'Epcam', 'Krt7', #Cholagiocytes
                 "Il7r","Cd3d","Trac", #T-cells
                 "Nkg7","Gzma", #NK-T
                 "Clec4f","Vsig4",#Kupffer cells
                 "Lgals3","Ccr2", "Lyz2", #Macrophages
                 "Klrd1","Mycl","Cd74","Itgax", # Monocytes
                 "Cxcr2", "S100a9","Ngp", #Neutorphils
                 "Bpgm","Hba-a1", #Erythrocyte
                 "Cd79a","Ebf1","Ms4a1", #B-cell
                 "Jchain",          #Plasma
                 "Siglech","Cox6a2", #Siglech
                 "Alb","Serpina1a", #Hepatocytes
                 "Apoa1","Serpina3k",
                 "Abcc2","Sema4g")

# *****************************************************************************
# ****************************************Reading all data individually
# *****************************************************************************
# Reading dataset1
# i=1
data <- Read10X_h5("Data/RS025_whole_liver_CCl4/RS025_filtered_gene_bc_matrices_h5.h5")
dim(data)
# data$prefix
sampleName <-  'RS025' #gsub("data/","",data.dirs[i])
# get seurat obj and get results
dataset1 <- CreateSeuratObject(counts = data, project = sampleName, min.cells = 3, min.features = 200)
dataset1
# pdf(file = paste0("results/reports/sample_",sampleName,".pdf"))
# plot.new()
# text(x=.5, y=0.5,paste0("Sample: ",sampleName))
# text(x=.5, y=0.4,paste0("Dimensions: ", paste(dim(dataset1), collapse = "x")))
dataset1[["percent.mt"]] <- PercentageFeatureSet(dataset1, pattern = "^mt-")
dataset1[["percent.ribo"]] <- PercentageFeatureSet(dataset1, pattern = "^Rpl|^Rps|^Mrps|^Mrpl")#
dataset1 <- subset(dataset1, subset = percent.mt < 20 & nFeature_RNA > 200 & nFeature_RNA < 4000 & nCount_RNA < 20000 )#&
VlnPlot(dataset1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo"), ncol = 4)
FeatureScatter(dataset1, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(dataset1, feature1 = "nFeature_RNA", feature2 = "percent.mt")
FeatureScatter(dataset1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dataset1 <- NormalizeData(object = dataset1)
dataset1 <- FindVariableFeatures(object = dataset1, nfeatures = 3000)
dataset1 <- ScaleData(object = dataset1)
dataset1 <- RunPCA(object = dataset1)
dataset1 <- FindNeighbors(object = dataset1, dims = 1:30)
dataset1 <- RunUMAP(object = dataset1, reduction = "pca", dims = 1:30)
dataset1 <- FindClusters(object = dataset1)

FeaturePlot(object = dataset1, features = c("percent.mt"))
FeaturePlot(object = dataset1, features = c("percent.ribo"))
FeaturePlot(object = dataset1, features = c("nCount_RNA"))
FeaturePlot(object = dataset1, features = c("nFeature_RNA"))
DimPlot(object = dataset1, label = T)
DotPlot(dataset1, features = markerGenes,cluster.idents = T) + RotatedAxis()
# dev.off()
# dev.off()

# Reading dataset2
data <- Read10X_h5("Data/RS045_whole_liver_CCl4/RS045_filtered_feature_bc_matrix.h5")
dim(data)
# data$prefix
sampleName <- 'RS045' #gsub("data/","",data.dirs[i])
# get seurat obj and get results
dataset2 <- CreateSeuratObject(counts = data, project = sampleName, min.cells = 3, min.features = 200)
dataset2
# pdf(file = paste0("results/reports/sample_",sampleName,".pdf"))
# plot.new()
# text(x=.5, y=0.5,paste0("Sample: ",sampleName))
# text(x=.5, y=0.4,paste0("Dimensions: ", paste(dim(dataset2), collapse = "x")))
dataset2[["percent.mt"]] <- PercentageFeatureSet(dataset2, pattern = "^mt-")
dataset2[["percent.ribo"]] <- PercentageFeatureSet(dataset2, pattern = "^Rpl|^Rps|^Mrps|^Mrpl")#
dataset2 <- subset(dataset2, subset = percent.mt < 20 & nFeature_RNA > 200 & nFeature_RNA < 4000 & nCount_RNA < 20000 )#&
VlnPlot(dataset2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo"), ncol = 4)
FeatureScatter(dataset2, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(dataset2, feature1 = "nFeature_RNA", feature2 = "percent.mt")
FeatureScatter(dataset2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dataset2 <- NormalizeData(object = dataset2)
dataset2 <- FindVariableFeatures(object = dataset2, nfeatures = 3000)
dataset2 <- ScaleData(object = dataset2)
dataset2 <- RunPCA(object = dataset2)
dataset2 <- FindNeighbors(object = dataset2, dims = 1:30)
dataset2 <- RunUMAP(object = dataset2, reduction = "pca", dims = 1:30)
dataset2 <- FindClusters(object = dataset2)

FeaturePlot(object = dataset2, features = c("percent.mt"))
FeaturePlot(object = dataset2, features = c("percent.ribo"))
FeaturePlot(object = dataset2, features = c("nCount_RNA"))
FeaturePlot(object = dataset2, features = c("nFeature_RNA"))
DimPlot(object = dataset2, label = T)
DotPlot(dataset2, features = markerGenes,cluster.idents = T) + RotatedAxis()
# dev.off()
# dev.off()

# Reading dataset3
data <- Read10X_h5("Data/RS046_whole_liver_CCl4/RS046_filtered_feature_bc_matrix.h5")
dim(data)
# data$prefix
sampleName <- 'RS046' #gsub("data/","",data.dirs[i])
# get seurat obj and get results
dataset3 <- CreateSeuratObject(counts = data, project = sampleName, min.cells = 3, min.features = 200)
dataset3
# pdf(file = paste0("results/reports/sample_",sampleName,".pdf"))
# plot.new()
# text(x=.5, y=0.5,paste0("Sample: ",sampleName))
# text(x=.5, y=0.4,paste0("Dimensions: ", paste(dim(dataset3), collapse = "x")))
dataset3[["percent.mt"]] <- PercentageFeatureSet(dataset3, pattern = "^mt-")
dataset3[["percent.ribo"]] <- PercentageFeatureSet(dataset3, pattern = "^Rpl|^Rps|^Mrps|^Mrpl")#
dataset3 <- subset(dataset3, subset = percent.mt < 20 & nFeature_RNA > 200 & nFeature_RNA < 4000 & nCount_RNA < 20000 )#&
VlnPlot(dataset3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo"), ncol = 4)
FeatureScatter(dataset3, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(dataset3, feature1 = "nFeature_RNA", feature2 = "percent.mt")
FeatureScatter(dataset3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dataset3 <- NormalizeData(object = dataset3)
dataset3 <- FindVariableFeatures(object = dataset3, nfeatures = 3000)
dataset3 <- ScaleData(object = dataset3)
dataset3 <- RunPCA(object = dataset3)
dataset3 <- FindNeighbors(object = dataset3, dims = 1:30)
dataset3 <- RunUMAP(object = dataset3, reduction = "pca", dims = 1:30)
dataset3 <- FindClusters(object = dataset3)

FeaturePlot(object = dataset3, features = c("percent.mt"))
FeaturePlot(object = dataset3, features = c("percent.ribo"))
FeaturePlot(object = dataset3, features = c("nCount_RNA"))
FeaturePlot(object = dataset3, features = c("nFeature_RNA"))
DimPlot(object = dataset3, label = T)
DotPlot(dataset3, features = markerGenes,cluster.idents = T) + RotatedAxis()
# dev.off()
# dev.off()




# *****************************************************************************
# **************************************** integrating
# *****************************************************************************
# integrating
data.list <- list(dataset1,dataset2,dataset3)#,dataset4,dataset5)#,dataset6,dataset7)
names(data.list) <- c('RS025','RS045','RS046')#,'RS032','RS007') #gsub("data/","",data.dirs)
names(data.list) 


# rPCA based integration
# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = data.list,nfeatures = 3000)
# perform integration of all samples
data.anchors <- FindIntegrationAnchors(object.list = data.list, reduction = "rpca", anchor.features = features) #
date()


# create an 'integrated' data assay
data.combined <- IntegrateData(anchorset = data.anchors)#, dims = 1:10, k.weight = 50)
# specify that we will perform downstream analysis on the corrected data note that the original
# unmodified data still resides in the 'RNA' assay



DefaultAssay(data.combined) <- "integrated"
# Run the standard workflow for visualization and clustering
data.combined <- ScaleData(data.combined, verbose = FALSE)
data.combined <- RunPCA(data.combined, npcs = 50, verbose = FALSE)
ElbowPlot(data.combined, ndims = 50)
ndim <- 22
data.combined <- FindNeighbors(data.combined, reduction = "pca", dims = 1:ndim)
data.combined <- FindClusters(data.combined)#, resolution = 0.5
data.combined <- RunUMAP(data.combined, reduction = "pca", dims = 1:ndim)
DimPlot(data.combined, label = T)
DimPlot(data.combined,group.by = "orig.ident", label = F,shuffle = T)


DefaultAssay(data.combined) <- "RNA"
DotPlot(data.combined, features = markerGenes,cluster.idents = T) + RotatedAxis()
DotPlot(data.combined, features = markerGenes) + RotatedAxis()
dim(data.combined@assays$RNA@data)
dim(data.combined@assays$integrated@data)
DimPlot(data.combined, label = T)

# selecting cell groups
DimPlot(data.combined, reduction = "umap",label = T,pt.size = 1)
# https://satijalab.org/seurat/articles/visualization_vignette.html
plot <- DimPlot(data.combined, reduction = "umap")
select.cells <- CellSelector(plot = plot)
head(select.cells)
Idents(data.combined, cells = select.cells) <- 12
# Idents(data.combined)
DimPlot(data.combined,label = T, pt.size = 1)
table((data.combined@active.ident))
levels(data.combined@active.ident)
data.combined@active.ident <- factor(data.combined@active.ident, levels = c("0",  "1",  "2",  "3",  "4",  "5",  "6",  "7",  "8",  "9",  "10", "11", "12", "13", "14", "15", "16", "17", "18", "19","20","21","22","23","24","25", "26","27","28","29","30"))
DimPlot(data.combined,label = T, pt.size = 1)

DotPlot(data.combined, features = markerGenes,cluster.idents = T) + RotatedAxis()
DotPlot(data.combined, features = markerGenes) + RotatedAxis()

# get cluster markers
DefaultAssay(data.combined) <- "RNA"
markers <- FindAllMarkers(data.combined,only.pos = T, logfc.threshold = 0.25)
head(markers,20)
# write.csv(markers, file = "results/wLiv_CCl4/markers_mus_wLiv_CCl4_clusters.csv")

# giving cell types
DimPlot(data.combined,label = T)
# data.combined[["Cell.clusters"]] <- Idents(object = data.combined)
DimPlot(data.combined,label = T,group.by = "Cell.clusters")
data.combined <- SetIdent(data.combined, value = "Cell.clusters")
data.combined <- RenameIdents(object = data.combined, 
                        '24' = "Cholangiocytes",
                        
                        '0' = "Hepatocytes",
                        '1' = "Hepatocytes",
                        '3' = "Hepatocytes",
                        '5' = "Hepatocytes",
                        '10' = "Hepatocytes",
                        '13' = "Hepatocytes",
                        '18' = "Hepatocytes",
                        
                       
                        '12' = "HSC",
                        '25' = "Mesothelial",
                        
                        '4' = "Endothelial",
                        '15' = "Endothelial",
                        '22' = "Endothelial",
                        '27' = "Endothelial",
                        
                        '8' = "Macrophages",
                        '14' = "Kupffer cell",
                        '16' = "Monocytes",
                        
                        '2' = "B cell",
                        '26' = "Plasma cell",
                       
                        '6' = "T_NKT_NK cell",
                        '17' = "T_NKT_NK cell",
                        
                        '9' = "Erythroid",
                        '19' = "pDC",
                        '11' = "Neutrophils",
                        
                        '7' = "Mixed",
                        '20' = "Mixed",
                        '21' = "Mixed",
                        '23' = "Mixed",
                        '21' = "Mixed",
                        '28' = "Mixed",
                        '29' = "Mixed",
                        '30' = "Mixed")


data.combined[["Cell.types"]] <- Idents(object = data.combined)
DimPlot(data.combined,label = T, group.by = "Cell.types")

table(data.combined@meta.data$Cell.types)
table(data.combined@meta.data$Cell.clusters)
table(data.combined@meta.data$cell.type)

DimPlot(data.combined, reduction = "umap",group.by = "Cell.types", label = F)
DimPlot(data.combined, reduction = "umap",group.by = "Cell.clusters", label = F)
data.combined <- SetIdent(data.combined, value = "Cell.clusters")


# saveRDS(data.combined, file = "results/wLiv_CCl4/seuratRPCAmerge_mus_wLiv_CCl4_Samples.rds")

# *****************************************************************************
# **************************************** working with integrated data
# *****************************************************************************
data.combined <- readRDS(file = "results/wLiv_CCl4/seuratRPCAmerge_mus_wLiv_CCl4_Samples.rds")
DimPlot(data.combined, label = T)
DimPlot(data.combined, reduction = "umap",group.by = "Cell.types", label = T)
DimPlot(data.combined, reduction = "umap",group.by = "Cell.clusters", label = T)
DimPlot(data.combined,group.by = "orig.ident", label = F,shuffle = T)
DimPlot(data.combined, label = T, split.by = "orig.ident")

DefaultAssay(data.combined) <- "RNA"
DotPlot(data.combined, features = markerGenes,cluster.idents = T, group.by = "Cell.clusters") + RotatedAxis()
DotPlot(data.combined, features = markerGenes, group.by = "Cell.clusters") + RotatedAxis()
DotPlot(data.combined, features = markerGenes,cluster.idents = T, group.by = "Cell.types") + RotatedAxis()
DotPlot(data.combined, features = markerGenes, group.by = "Cell.types") + RotatedAxis()


# ************************************************** cy-my signature analysis
# this part is performed after identifying the cy-my signature using script 'createSignature.R'
# cy-my signature
library(readxl)
dataset <- data.frame(read_excel("results/cymySig_secreated.xlsx", skip = 1))
head(dataset)

cySig <- na.omit(dataset$cySig2)
mySig <- na.omit(dataset$mySgi2)


# reclustering HSC for cy my
data.combined <- SetIdent(data.combined, value = "Cell.types")
req_subset <- subset(data.combined, idents = c("HSC"), invert=F)
req_subset
DimPlot(req_subset,label = T, group.by = "ident",pt.size = 1)
DimPlot(req_subset,label = T, group.by = "Cell.types")
DimPlot(req_subset,label = T, group.by = "Cell.clusters")
# DimPlot(req_subset,label = T, group.by = "Cell.typesFine")

DefaultAssay(req_subset) <- "integrated"
req_subset <- NormalizeData(req_subset, normalization.method = "LogNormalize", scale.factor = 10000)
req_subset <- FindVariableFeatures(req_subset, selection.method = "vst", nfeatures = 2000)
# all.genes <- rownames(req_subset)
req_subset <- ScaleData(req_subset, features = VariableFeatures(object = req_subset), vars.to.regress = c("nCount_RNA", "percent.mt"))
req_subset <- RunPCA(req_subset, features = VariableFeatures(object = req_subset))
ElbowPlot(req_subset, ndims = 50)
ndim=10
req_subset <- FindNeighbors(req_subset, dims = 1:ndim)
req_subset <- FindClusters(req_subset)
req_subset <- RunUMAP(req_subset, reduction = "pca", dims = 1:ndim,min.dist = 0.9, spread = 1.745)
DimPlot(req_subset,label = T, pt.size = 3)
DimPlot(req_subset,label = T, group.by = "orig.ident")


#*# selecting cell groups from a cluster and reordering the cluster numbers
DimPlot(req_subset, reduction = "umap",label = T,pt.size = 1)
# https://satijalab.org/seurat/articles/visualization_vignette.html
plot <- DimPlot(req_subset, reduction = "umap")
select.cells <- CellSelector(plot = plot)
head(select.cells)
Idents(req_subset, cells = select.cells) <- 6
# Idents(req_subset)
DimPlot(req_subset,label = T, pt.size = 1)
table((req_subset@active.ident))
levels(req_subset@active.ident)
req_subset@active.ident <- factor(req_subset@active.ident, levels = c("0",  "1",  "2",  "3",  "4",  "5",  "6"))#,  "7",  "8",  "9",  "10", "11", "12", "13", "14", "15", "16", "17", "18", "19","20","21","22","23","24","25", "26","27","28","29","30"))
DimPlot(req_subset,label = T, pt.size = 3)


sig <- cySig
sigName <- "cySig"
callName <- paste0(sigName,"1")
req_subset <- AddModuleScore(req_subset,features = list(sig), name = sigName)
FeaturePlot(req_subset, callName, max.cutoff = "q99", pt.size = 3) + scale_color_viridis_c(option = "viridis")
VlnPlot(req_subset, callName, sort = T)
# VlnPlot(req_subset, callName, group.by = "Cell.clusters", sort = T)

sig <- mySig
sigName <- "mySig"
callName <- paste0(sigName,"1")
req_subset <- AddModuleScore(req_subset,features = list(sig), name = sigName)
FeaturePlot(req_subset, callName, max.cutoff = "q99", pt.size = 3) + scale_color_viridis_c(option = "viridis")
VlnPlot(req_subset, callName, sort = T)
# VlnPlot(req_subset, callName, group.by = "Cell.clusters", sort = T)

VlnPlot(req_subset, features = c("cySig1","mySig1"), sort = F)
DotPlot(req_subset, features = c("cySig1","mySig1")) + RotatedAxis()
FeaturePlot(req_subset, "Lrat", max.cutoff = "q99", pt.size = 3) + scale_color_viridis_c(option = "viridis")

# quantify the cy-my score median in different clusters
df <- data.frame(clust=req_subset@active.ident, cy=req_subset@meta.data$cySig1, my=req_subset@meta.data$mySig1)
head(df)
ggplot(df, aes(x=clust, y=cy))+geom_violin()

library(tidyverse)
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

# req_subset[["Cell.clusters.HSC"]] <- Idents(object = req_subset)
DimPlot(req_subset, reduction = "umap",group.by = "Cell.clusters.HSC", label = T)
req_subset <- SetIdent(req_subset, value = "Cell.clusters.HSC")
req_subset <- RenameIdents(object = req_subset,
                             '2' = "cyHSC",
                             '3' = "cyHSC",
                             '4' = "cyHSC",
                             '5' = "cyHSC",
                             '6' = "cyHSC",
                             
                           
                             '0' = "myHSC",
                             '1' = "myHSC")


req_subset[["Cell.types.cymy"]] <- Idents(object = req_subset)
DimPlot(req_subset,label = T, group.by = "Cell.types.cymy")
table(req_subset@active.ident)

cols_cymy <- c("#FF9900","#000099")
DimPlot(req_subset, reduction = "umap",group.by = "Cell.types.cymy", pt.size = 2,cols = cols_cymy, label = F)

# saveRDS(req_subset, file = "results/wLiv_CCl4/seuratRPCAmerge_mus_wLiv_CCl4_Samples_HSC.rds")

req_subset <- readRDS(file = "results/wLiv_CCl4/seuratRPCAmerge_mus_wLiv_CCl4_Samples_HSC.rds")
DimPlot(req_subset, reduction = "umap",label = T,pt.size = 3)

# getting the cy-my annot back to wLiv object
head(req_subset@meta.data)
head(req_subset@meta.data$Cell.types.cymy)
data.combined <- SetIdent(data.combined, value = "Cell.types")
DimPlot(data.combined,label = T)
# myHSC
tmppos <- which(req_subset@meta.data$Cell.types.cymy=="myHSC")
tmpCells <- rownames(req_subset@meta.data)[tmppos]
Idents(object = data.combined, cells = tmpCells) <- "myHSC"
DimPlot(data.combined,label = T)
# cyHSC
tmppos <- which(req_subset@meta.data$Cell.types.cymy=="cyHSC")
tmpCells <- rownames(req_subset@meta.data)[tmppos]
Idents(object = data.combined, cells = tmpCells) <- "cyHSC"
DimPlot(data.combined,label = T)
# writing to a new metadata column
data.combined[["Cell.types.cymy"]] <- Idents(object = data.combined)
DimPlot(data.combined,label = T, group.by = "Cell.types.cymy")

# saveRDS(data.combined, file = "results/wLiv_CCl4/seuratRPCAmerge_mus_wLiv_CCl4_Samples.rds")


#******************************************************************************
#***************************************** pathway analysis
#******************************************************************************
# pathway analysis of hep clusters 10 and 13
library(enrichR)

data.combined <- readRDS(file = "results/wLiv_CCl4/seuratRPCAmerge_mus_wLiv_CCl4_Samples.rds")
DimPlot(data.combined, label = T)
DimPlot(data.combined, reduction = "umap",group.by = "Cell.clusters", label = T)

# reclustering Hepatocytes for progenitors
data.combined <- SetIdent(data.combined, value = "Cell.types")
req_subset <- subset(data.combined, idents = c("Hepatocytes"), invert=F)
req_subset
DimPlot(req_subset,label = T, pt.size = 1)
DimPlot(req_subset,label = T, group.by = "Cell.types")
DimPlot(req_subset,label = T, group.by = "Cell.clusters")
req_subset <- SetIdent(req_subset, value = "Cell.clusters")
DimPlot(req_subset,label = T, pt.size = 1)
levels(req_subset@active.ident)

# get cluster markers
DefaultAssay(req_subset) <- "RNA"
markers <- FindAllMarkers(req_subset,only.pos = T, logfc.threshold = 0.25)
markers <- markers[order(markers$cluster,markers$avg_log2FC,decreasing = T),]
head(markers,20)
# write.csv(markers, file = "results/wLiv_CCl4/markers_mus_wLiv_CCl4_clusters_Hep_only.csv")

markers <- read.csv(file = "results/wLiv_CCl4/markers_mus_wLiv_CCl4_clusters_Hep_only.csv", row.names = 1)
head(markers)
table(markers$cluster)
clustSel <- "10" #"10" "13"
tmppos <- which(markers$cluster==clustSel)
tmp <- markers[tmppos,]
tmp[1:20,]
tmppos <- which(tmp$avg_log2FC>0.5)
tmppos <- tmppos[1:if(length(tmppos)>100) 100 else length(tmppos)] #keep a max of top 100 genes
genelist <- tmp$gene[tmppos]

# pathway analysis
dbs <- listEnrichrDbs()
# dbs
dbs <- c("GO_Biological_Process_2015")
enriched <- enrichr(genelist, dbs)
head(enriched)
# writing to file
df <- enriched[[1]]
head(df)
# write.csv(df, file = paste0("results/wLiv_CCl4/pathway_GO_cluster_",clustSel,".csv"))




#******************************************************************************
#***************************************** working with saved data
#******************************************************************************
data.combined <- readRDS(file = "results/wLiv_CCl4/seuratRPCAmerge_mus_wLiv_CCl4_Samples.rds")
DimPlot(data.combined, label = T)
DimPlot(data.combined, reduction = "umap",group.by = "Cell.types", label = T)
DimPlot(data.combined, reduction = "umap",group.by = "Cell.types.cymy", label = T)
DimPlot(data.combined, reduction = "umap",group.by = "Cell.clusters", label = T)
DimPlot(data.combined,group.by = "orig.ident", label = F,shuffle = T)
DimPlot(data.combined, label = T, split.by = "orig.ident")

DefaultAssay(data.combined) <- "RNA"
DotPlot(data.combined, features = markerGenes,cluster.idents = T, group.by = "Cell.clusters") + RotatedAxis()
DotPlot(data.combined, features = markerGenes, group.by = "Cell.clusters") + RotatedAxis()

# get only the clean clusters
data.combined <- SetIdent(data.combined, value = "Cell.types")
req_subset <- subset(data.combined, idents = c("Mixed"), invert=T)
req_subset
DimPlot(req_subset,label = T, group.by = "ident",pt.size = 1)
DimPlot(req_subset,label = T, group.by = "Cell.types")
DimPlot(req_subset,label = T, group.by = "Cell.clusters")
DimPlot(req_subset,label = T, group.by = "Cell.types.cymy")

# get the cell type annotation
colnames(data.combined@meta.data)
celltypes <- data.combined@meta.data[,c(1,2,3,4,9,10,11)]
head(celltypes)
# write.csv(celltypes, file = "results/wLiv_CCl4/mus_wLiv_CCl4_cellAnnotation.csv")

#************************************************#************************************************#************************************************
#************************************************ get the data required for CPDB
#************************************************#************************************************#************************************************

data_req <- req_subset
DimPlot(data_req,label = T, group.by = "Cell.types")

# removing erythroid cells
data_req <- subset(data_req, idents = c("Erythroid"), invert=T)
data_req
DimPlot(data_req,label = T, group.by = "ident",pt.size = 1)
DimPlot(data_req,label = T, group.by = "Cell.types")

dim(data_req@assays$RNA@counts)
dim(data_req@assays$RNA@data)
dim(data_req@assays$integrated@data)

# getting the normalized counts
rawCounts <- as.matrix(data_req@assays$RNA@counts)
rawCounts[1:5,1:5]
dim(rawCounts)


# getting the cell types
cellLabels <- data.frame(cells=colnames(data_req@assays$RNA@counts), cellType=data_req@meta.data$Cell.types)
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

# Mouse to human gene conversion
homologs_human_mouse <- readRDS(file="D:/pCloud Sync/Columbia/CommonResource/Data/homologs_mouse_human_geneSymbol.rds")
head(homologs_human_mouse)
dim(homologs_human_mouse)
dim(temp_norm_mat)
rownames(temp_norm_mat)[1:10]
#get the mouse genes which are present
temp_pos <- which(rownames(temp_norm_mat) %in% homologs_human_mouse[,"MGI.symbol"])
temp_pos <- temp_pos[!is.na(temp_pos)]
mus2humGenesNormData <- temp_norm_mat[temp_pos,]
mus2humGenesNormData[1:5,1:5]
#match with human genes
temp_pos <- match(rownames(mus2humGenesNormData), homologs_human_mouse[,"MGI.symbol"])
rownames(mus2humGenesNormData)[1:5]
homologs_human_mouse[temp_pos[1:5],]
all(rownames(mus2humGenesNormData)==homologs_human_mouse[temp_pos,"MGI.symbol"])
rownames(mus2humGenesNormData) <- homologs_human_mouse[temp_pos,"HGNC.symbol"]
mus2humGenesNormData[1:5,1:5]
dim(mus2humGenesNormData)
all(cellLabels[,"cells"]==colnames(mus2humGenesNormData))


#**********************************getting the data for CellphoneDB
write.table(t(c("Genes", colnames(mus2humGenesNormData))),file = "wLiv_CCl4/cpdb/data_cellphoneDB_countNorm_mus_hum_mgiJax_wLivCCl4combined.txt", quote = FALSE, sep = "\t", append = FALSE, col.names = FALSE, row.names = FALSE)
write.table(mus2humGenesNormData,file = "wLiv_CCl4/cpdb/data_cellphoneDB_countNorm_mus_hum_mgiJax_wLivCCl4combined.txt", quote = FALSE, sep = "\t", append = TRUE, col.names = FALSE)

# reading back the counts to check the metadata
mus2humGenesNormData <- read.table(file = "wLiv_CCl4/cpdb/data_cellphoneDB_countNorm_mus_hum_mgiJax_wLivCCl4combined.txt", sep = "\t", header = T, row.names = 1)
colnames(mus2humGenesNormData) <- gsub("\\.","-",colnames(mus2humGenesNormData))
mus2humGenesNormData[1:5,1:5]


# #getting the metadata; all cell types
cellLabels_req <- cellLabels 
metaData <- data.frame(Cell=cellLabels_req[,"cells"],cell_type=cellLabels_req[,"cellType"])
all(metaData[,"Cell"]==colnames(mus2humGenesNormData))
head(metaData)
dim(metaData)
table(metaData[,"cell_type"])
# write.table(metaData, file = "wLiv_CCl4/cpdb/data_cellphoneDB_metaData_mus_wLivCCl4combined_cellTypes.txt", sep = '\t', quote = F, row.names = F, col.names = T)

# getting cyHSC and myHSC
# getting the cell types
cellLabels <- data.frame(cells=colnames(data_req@assays$RNA@counts), cellType=data_req@meta.data$Cell.types.cymy, cellClust=data_req@meta.data$Cell.clusters)
head(cellLabels)
dim(cellLabels)
table(cellLabels$cellType)
all(colnames(mus2humGenesNormData)==as.character(cellLabels[,"cells"]))
cellLabels_req <- cellLabels
head(cellLabels_req)
table(cellLabels_req[,"cellType"])
metaData <- data.frame(Cell=cellLabels_req[,"cells"],cell_type=cellLabels_req[,"cellType"])
all(metaData[,"Cell"]==colnames(mus2humGenesNormData))
head(metaData)
dim(metaData)
table(metaData[,"cell_type"])
# write.table(metaData, file = "wLiv_CCl4/cpdb/data_cellphoneDB_metaData_mus_wLivCCl4combined_cellTypes_HSC_subtypes.txt", sep = '\t', quote = F, row.names = F, col.names = T)


# getting cyHSC and myHSC and Hep subclusters
# renaming the Hep
tmp <- as.character(cellLabels_req$cellType)
temp_pos <- which(cellLabels_req$cellType=="Hepatocytes")
table(tmp[temp_pos])
head(cellLabels_req[temp_pos,])
tmp[temp_pos] <- paste0("Hep_C",cellLabels_req$cellClust[temp_pos])
table(tmp[temp_pos])
table(tmp)
cellLabels_req$cellSubType <- tmp
head(cellLabels_req)

table(cellLabels_req[,"cellType"])
table(cellLabels_req[,"cellSubType"])
metaData <- data.frame(Cell=cellLabels_req[,"cells"],cell_type=cellLabels_req[,"cellSubType"])
all(metaData[,"Cell"]==colnames(mus2humGenesNormData))
head(metaData)
dim(metaData)
table(metaData[,"cell_type"])
# write.table(metaData, file = "wLiv_CCl4/cpdb/data_cellphoneDB_metaData_mus_wLivCCl4combined_cellTypes_HSC_Hep_subtypes.txt", sep = '\t', quote = F, row.names = F, col.names = T)
# 

metaData <- read.delim("cpdb/wLiv_CCl4/data_cellphoneDB_metaData_mus_wLivCCl4combined_cellTypes_HSC_Hep_subtypes.txt")
head(metaData)
data_req[["cellSubType"]] <- metaData[,"cell_type"]
DimPlot(data_req,label = T, group.by = "cellSubType",pt.size = 1)

DimPlot(object = data_req, group.by = "cellSubType", label = T,  pt.size = 1)  + labs(title = "Mouse CCl4 whole-liver")
# ggsave(filename = paste0("results/wLiv_CCl4/figs/mus_wLiv_CCl4_cellsubtypes_label.pdf"), width = 8, height = 8)

DimPlot(object = data_req, group.by = "cellSubType", label = F,  pt.size = 1)  + labs(title = "Mouse CCl4 whole-liver")+ NoLegend()
# ggsave(filename = paste0("results/wLiv_CCl4/figs/mus_wLiv_CCl4_cellsubtypes_noLabel.pdf"), width = 8, height = 8)


DimPlot(data_req,label = T,pt.size = 1)
DimPlot(object = data_req, group.by = "cellSubType", label = F,  pt.size = 3)  + labs(title = "Mouse CCl4 whole-liver Hep")+ NoLegend() + xlim(-2,10) +ylim(-10,2)
# ggsave(filename = paste0("results/wLiv_CCl4/figs/mus_wLiv_CCl4_cellsubtypes_Hep_noLabel.pdf"), width = 8, height = 8)

FeaturePlot(object = data_req, features = "Mki67", max.cutoff = "q99", pt.size = 3,order = F) + scale_color_viridis_c() +  xlim(-2,10) +ylim(-10,2)
# ggsave(filename = paste0("results/wLiv_CCl4/figs/mus_wLiv_CCl4_cellsubtypes_Hep_Mki67.pdf"), width = 8, height = 8)


#******************************************************#******************************************************
#******************************************************PLOTS
#******************************************************#******************************************************
library(RColorBrewer)
library("ggsci")
library("scales")

data.combined <- readRDS(file = "results/wLiv_CCl4/seuratRPCAmerge_mus_wLiv_CCl4_Samples.rds")
# get only the clean clusters
data.combined <- SetIdent(data.combined, value = "Cell.types")
req_subset <- subset(data.combined, idents = c("Mixed"), invert=T)
req_subset
DimPlot(req_subset,label = T, group.by = "ident",pt.size = 1)
DimPlot(req_subset,label = T, group.by = "Cell.types")
DimPlot(req_subset,label = T, group.by = "Cell.clusters")
DimPlot(req_subset,label = T, group.by = "Cell.types.cymy")
# DimPlot(req_subset,label = T, group.by = "cellSubType")


# CELL TYPES UMAP
data_req <- req_subset
DimPlot(object = data_req, group.by = "ident", label = T)
levels(data_req@active.ident)


# [1] "Cholangiocytes" "Hepatocytes"    "HSC"            "Mesothelial"    "Endothelial"    "Macrophages"    "Kupffer cell"   "Monocytes"     
# [9] "B cell"         "Plasma cell"    "T_NKT_NK cell"  "Erythroid"      "pDC"            "Neutrophils"   

cols <- c("#9900CCFF","#00468BFF","#1B1919FF","#FFCC00FF","#99991EFF","#FDBF6F","#FF7F00","#6699FFFF","#42B540FF","#B2df8a","#AD002AFF","grey","#58593FFF","#CC0000FF","lightgreen")
DimPlot(object = data_req, group.by = "ident", label = F, cols = cols, pt.size = 2) 

DimPlot(object = data_req, group.by = "Cell.types", label = T, cols = cols, pt.size = 1)  + labs(title = "Mouse CCl4 whole-liver")
ggsave(filename = paste0("results/wLiv_CCl4/figs/mus_wLiv_CCl4_celltypes_label.pdf"))
DimPlot(object = data_req, group.by = "Cell.types", label = F, cols = cols, pt.size = 1)  + labs(title = "Mouse CCl4 whole-liver")+ NoLegend()
ggsave(filename = paste0("results/wLiv_CCl4/figs/mus_wLiv_CCl4_celltypes_noLabel.pdf"))

# with cy-my
# ,"#1B1919FF"
# cols_cymy <- c("#FF9900","#000099")
cols <- c("#FFCCCCFF","#FF00CCFF","#9900CCFF","#00468BFF","#FFCC00FF","#99991EFF","#FDBF6F","#FF7F00","#6699FFFF","#42B540FF","#B2df8a","#AD002AFF","grey","#58593FFF","#CC0000FF","lightgreen")
DimPlot(object = data_req, group.by = "Cell.types.cymy", label = T, cols = cols, pt.size = 1)  + labs(title = "Mouse CCl4 whole-liver")

DimPlot(object = data_req, group.by = "Cell.types.cymy", label = T, cols = cols, pt.size = 1)  + labs(title = "Mouse CCl4 whole-liver")
ggsave(filename = paste0("results/wLiv_CCl4/figs/mus_wLiv_CCl4_celltypes_cymy_label.pdf"))
DimPlot(object = data_req, group.by = "Cell.types.cymy", label = F, cols = cols, pt.size = 1)  + labs(title = "Mouse CCl4 whole-liver")+ NoLegend()
ggsave(filename = paste0("results/wLiv_CCl4/figs/mus_wLiv_CCl4_celltypes_cymy_noLabel.pdf"))


DimPlot(object = data_req, group.by = "orig.ident", label = F, cols = cols, pt.size = 0.5,shuffle =T) + scale_color_aaas() + labs(title = "Mouse CCl4 whole-liver")
ggsave(filename = paste0("results/wLiv_CCl4/figs/mus_wLiv_CCl4_samples.pdf"))


# CELL TYPES DOT PLOT
markerGenes <- c("Alb","Serpina1a",
                 'Epcam', 'Krt7',
                 "Lrat","Rgs5","Col3a1","Col1a1", #,"Acta2","Lox"
                 'Upk3b','Msln',
                 "Kdr","Aqp1",
                 "Lgals3","Ccr2",
                 "Clec4f","Vsig4",
                 "Klrd1","Itgax",
                 "Top2a","Mki67",
                 "Cd79a","Ebf1",
                 "Mzb1","Jchain",
                 "Il7r","Cd3d",
                 "Nkg7",
                 "Hba-a1","Snca",
                 "Siglech","Klk1",
                 "Cxcr2","S100a9")
DotPlot(data_req, features = markerGenes, group.by = "Cell.types") + RotatedAxis() & scale_color_distiller(palette = "RdBu")
ggsave(filename = "results/wLiv_CCl4/figs/mus_wLiv_CCl4_celltypes_dotplot_markers.pdf", width = 10, height = 4)



genes <- c(unlist(strsplit("Ddr1, Met, Has2, Hhip, Cxcl12",split = ", "))) #Hgf, Col1a1, 
# markers required
ptsize <- 1.5
genes <- intersect(genes, rownames(data_req@assays$RNA@data))
FeaturePlot(object = data_req, features = c(genes[1]), max.cutoff = "q99", pt.size = ptsize,order = T) + scale_color_viridis_c()

for(gene in genes[3:4]){
  FeaturePlot(object = data_req, features = c(gene), max.cutoff = "q99", pt.size = ptsize,order = T) + scale_color_viridis_c()
  ggsave(filename = paste0("results/wLiv_CCl4/figs/mus_wLiv_CCl4_gene_",gene,".pdf"))
}

gene <- "Col1a1"
FeaturePlot(object = data_req, features = c(gene), max.cutoff = "q99", pt.size = ptsize,order = F) + scale_color_viridis_c()
ggsave(filename = paste0("results/wLiv_CCl4/figs/mus_wLiv_CCl4_gene_",gene,".pdf"))

gene <- "Hgf"
FeaturePlot(object = data_req, features = c(gene), max.cutoff = "q99", pt.size = ptsize,order = F) + scale_color_viridis_c()
ggsave(filename = paste0("results/wLiv_CCl4/figs/mus_wLiv_CCl4_gene_",gene,".pdf"))

gene <- "Cxcl12"
FeaturePlot(object = data_req, features = c(gene), max.cutoff = "q99", pt.size = ptsize,order = F) + scale_color_viridis_c()
ggsave(filename = paste0("results/wLiv_CCl4/figs/mus_wLiv_CCl4_gene_",gene,".pdf"))


gene <- "Ddr1"
FeaturePlot(object = data_req, features = c(gene),cells = WhichCells(object = data_req, idents = "Hepatocytes"), max.cutoff = "q99", pt.size = 3,order = T) + scale_color_viridis_c() + xlim(-5,12) +ylim(-10,5)
ggsave(filename = paste0("results/wLiv_CCl4/figs/mus_wLiv_CCl4_Hep_gene_",gene,".pdf"), width = 8, height = 8)


# cy-my signature
library(readxl)
dataset <- data.frame(read_excel("results/cymySig_secreated.xlsx", skip = 1))
head(dataset)

cySig <- na.omit(dataset$cySig2)
mySig <- na.omit(dataset$mySgi2)

sig <- cySig
sigName <- "cySig"
callName <- paste0(sigName,"1")
data_req <- AddModuleScore(data_req,features = list(sig), name = sigName)
FeaturePlot(data_req, callName, max.cutoff = "q99") + scale_color_viridis_c(option = "viridis")
FeaturePlot(data_req, callName,pt.size = ptsize,order = T) + scale_color_viridis_c(option = "viridis") +labs(title = "Mus CCl4 wLiver cy signature")
# ggsave(filename = paste0("results/wLiv_CCl4/figs/mus_wLiv_CCl4_signature_cy.pdf"))


sig <- mySig
sigName <- "mySig"
callName <- paste0(sigName,"1")
data_req <- AddModuleScore(data_req,features = list(sig), name = sigName)
FeaturePlot(data_req, callName,pt.size = ptsize,order = T) + scale_color_viridis_c(option = "viridis")+labs(title = "Mus CCl4 wLiver my signature")
# ggsave(filename = paste0("results/wLiv_CCl4/figs/mus_wLiv_CCl4_signature_my.pdf"))


# Hgf violin plot and stat
gene <- "Hgf" 
p <- VlnPlot(data_req, features = gene, idents = c("HSC","Endothelial")) #+geom_boxplot()
df <- p$data
colnames(df)[1] <- "gene"
head(df)
tmppos <- which(df$ident=="HSC")
tmpx <- df$gene[tmppos]
tmpy <- df$gene[-tmppos]
tmp <- wilcox.test(x=tmpx, y=tmpy)
tmp$p.value

VlnPlot(data_req, features = gene, idents = c("HSC","Endothelial")) +
  ggtitle(paste0("wLiv CCl4 ",gene," Mann Whitney p ",scientific(tmp$p.value,digits = 3)))
ggsave(filename = paste0("results/wLiv_CCl4/figs/mus_wLiv_CCl4_HSC_Endo_vln_",gene,".pdf"),width = 5)

