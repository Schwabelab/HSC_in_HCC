# 2021 March 32th created
# this is the Seurat3 analysis of RS042


setwd("D:/pCloud Sync/Columbia/SchwabeLab/Aveline_HSCnHCC/github/HSCinHCC/")

library(Seurat)
library(ggplot2)
library(monocle3)
library(SeuratWrappers)
library(patchwork)
library(magrittr)

markerGenes <- c("Lrat","Myh11","Actg2","Rgs5","Acta2","Lox","Col3a1","Col1a1","Mgp","Eln","Mfap4","Dpt","Gsn","Kdr","Aqp1","Top2a","Mki67",'Upk3b','Msln', 'Epcam', 'Krt7',"Il7r","Cd3d","Trac","Nkg7","Clec4f","Vsig4","Lgals3","Ccr2","Bpgm","Cd79a","Ebf1","Jchain","Siglech","Alb","Serpina1a")


# ************************read the data
data <- Read10X_h5("data/RS042_HSC_TAZ_FPC/filtered_feature_bc_matrix_RS042.h5")
rs042 <- CreateSeuratObject(counts = data, project = "RS042 TAZ_FPC", min.cells = 3, min.features = 200)
rs042
rs042[["percent.mt"]] <- PercentageFeatureSet(rs042, pattern = "^mt-")
rs042[["percent.ribo"]] <- PercentageFeatureSet(rs042, pattern = "^Rpl|^Rps|Mrps|^Mrpl")#

# Visualize QC metrics as a violin plot
VlnPlot(rs042, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# FeatureScatter(rs042, feature1 = "nCount_RNA", feature2 = "percent.mt")
# FeatureScatter(rs042, feature1 = "nFeature_RNA", feature2 = "percent.mt")
# FeatureScatter(rs042, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p1 <- FeatureScatter(rs042, feature1 = "nCount_RNA", feature2 = "percent.mt") + NoLegend()
p2 <- FeatureScatter(rs042, feature1 = "nFeature_RNA", feature2 = "percent.mt") + NoLegend()
p3 <- FeatureScatter(rs042, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# p4 <- FeatureScatter(liv_01_cov, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "percent.mt")
# # keep only the required cells
wrap_plots(list(p1,p2,p3),nrow = 2,ncol = 2)
# keep only the required cells
rs042 <- subset(rs042, subset = nFeature_RNA > 200 & nCount_RNA < 35000 & percent.mt < 10)

VlnPlot(rs042, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p1 <- FeatureScatter(rs042, feature1 = "nCount_RNA", feature2 = "percent.mt") + NoLegend()
p2 <- FeatureScatter(rs042, feature1 = "nFeature_RNA", feature2 = "percent.mt") + NoLegend()
p3 <- FeatureScatter(rs042, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# p4 <- FeatureScatter(liv_01_cov, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "percent.mt")
# # keep only the required cells
wrap_plots(list(p1,p2,p3),nrow = 2,ncol = 2)


rs042 <- NormalizeData(rs042, normalization.method = "LogNormalize", scale.factor = 10000)
rs042 <- FindVariableFeatures(rs042, selection.method = "vst", nfeatures = 2000)
# all.genes <- rownames(rs042)
rs042 <- ScaleData(rs042, features = VariableFeatures(object = rs042), vars.to.regress = c("nCount_RNA", "percent.mt"))
rs042 <- RunPCA(rs042, features = VariableFeatures(object = rs042))
ElbowPlot(rs042, ndims = 50)
ndim=20
rs042 <- FindNeighbors(rs042, dims = 1:ndim)
rs042 <- FindClusters(rs042, resolution = 0.5)
rs042 <- RunUMAP(rs042, reduction = "pca", dims = 1:ndim)
DimPlot(rs042,label = T)
markerGenes <- c("Lrat","Myh11","Actg2","Rgs5","Acta2","Lox","Col3a1","Col1a1","Mgp","Eln","Mfap4","Dpt","Gsn","Kdr","Aqp1","Top2a","Mki67",'Upk3b','Msln', 'Epcam', 'Krt7',"Il7r","Cd3d","Trac","Nkg7","Clec4f","Vsig4","Lgals3","Ccr2","Bpgm","Cd79a","Ebf1","Jchain","Siglech","Alb","Serpina1a")
DotPlot(rs042, features = markerGenes) + RotatedAxis()


# getting marker genes
marker <- FindAllMarkers(rs042,only.pos = T)
# write.csv(marker, file = "results/RS042/markers.csv")

FeaturePlot(object = rs042, features = c("percent.mt"))
FeaturePlot(object = rs042, features = c("percent.ribo"))
FeaturePlot(object = rs042, features = c("nCount_RNA"))
FeaturePlot(object = rs042, features = c("nFeature_RNA"))
VlnPlot(object = rs042, features = c("nCount_RNA"))
VlnPlot(object = rs042, features = c("nFeature_RNA"))

# giving cell types
DimPlot(rs042,label = T)
# rs042[["Cell.clusters"]] <- Idents(object = rs042)
rs042 <- SetIdent(rs042, value = "Cell.clusters")
rs042 <- RenameIdents(object = rs042, 
                              '0' = "HSC",
                              '1' = "HSC",
                              '2' = "HSC",
                              '3' = "HSC",
                              '5' = "HSC",
                              '6' = "HSC",
                              
                              '9' = "Cycling/Mixed",
                              '8' = "VSMC",
                              
                              '15' = "Cholangiocyte",
                              '13' = "PF/Mesothelial",
                              '14' = "Hepatocytes",
                              
                              '4' = "Endothelial",
                              '12' = "Endothelial",
                              
                              '7' = "Immune",
                              '10' = "Low quality",
                              '11' = "Mixed")




rs042[["Cell.types"]] <- Idents(object = rs042)
DimPlot(rs042,label = T, group.by = "Cell.types")
table(rs042@active.ident)

DimPlot(rs042, reduction = "umap",group.by = "Cell.types", label = F)
DimPlot(rs042, reduction = "umap",group.by = "Cell.clusters", label = F)



# saveRDS(rs042, file = "results/RS042/seurat4_RS042.rds")
#****************************************************************saved object
rs042 <- readRDS(file = "results/RS042/seurat4_RS042.rds")
DimPlot(rs042,label = T)
DimPlot(rs042, reduction = "umap",group.by = "Cell.types", label = T)
DimPlot(rs042, reduction = "umap",group.by = "Cell.clusters", label = T)


# Recluster the HSC population 
rs042 <- SetIdent(rs042, value = "Cell.clusters")
rs042_subset <- subset(rs042, idents = c(0:3,5,6))
DimPlot(rs042_subset,label = T, group.by = "ident")
rs042_subset <- FindVariableFeatures(rs042_subset, selection.method = "vst", nfeatures = 2000)
# all.genes <- rownames(rs042_subset)
rs042_subset <- ScaleData(rs042_subset, features = VariableFeatures(object = rs042_subset), vars.to.regress = c("nCount_RNA", "percent.mt"))
rs042_subset <- RunPCA(rs042_subset, features = VariableFeatures(object = rs042_subset))
ElbowPlot(rs042_subset, ndims = 50)
ndim=8
rs042_subset <- FindNeighbors(rs042_subset, dims = 1:ndim)
rs042_subset <- FindClusters(rs042_subset, resolution = 0.5)
rs042_subset <- RunUMAP(rs042_subset, reduction = "pca", dims = 1:ndim)
DimPlot(rs042_subset, reduction = "umap",label = T)
DotPlot(rs042_subset, features = markerGenes) + RotatedAxis()
FeaturePlot(object = rs042_subset, features = c("Col1a1")) & scale_color_viridis_c()

# saveRDS(rs042_subset, file = "results/RS042/seurat4_RS042_HSC.rds")
#****************************************************************saved object
rs042_subset <- readRDS(file = "results/RS042/seurat4_RS042_HSC.rds")
DimPlot(rs042_subset,label = T)


# ************************************************** trajectory analysis
# getting to Monocle
cds <- as.cell_data_set(rs042_subset)
cds <- cluster_cells(cds)#,resolution = 0.05
p1 <- plot_cells(cds, show_trajectory_graph = FALSE)
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
wrap_plots(p1, p2)
cds <- learn_graph(cds)
plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE)
plot_cells(cds, label_groups_by_cluster = F, label_leaves = T, label_branch_points = T,graph_label_size=5)
# assigning root interactively
cds <- order_cells(cds)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)

# saveRDS(cds, file = "results/RS042/monocle_RS042_HSC.rds")
#saved object
cds <- readRDS(file = "results/RS042/monocle_RS042_HSC.rds")
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)

# Set the assay back as 'integrated'
# cds@assays$
# integrated.sub <- as.Seurat(cds, assay = "integrated")
# integrated.sub <- as.Seurat(cds, assay = "counts",project = "cell_data_set")
integrated.sub <- as.Seurat(cds, assay = NULL,project = "cell_data_set")
FeaturePlot(integrated.sub, "monocle3_pseudotime")

rs042_subset@meta.data$monocle3_pseudotime <- integrated.sub@meta.data$monocle3_pseudotime
FeaturePlot(rs042_subset, "monocle3_pseudotime")

# saveRDS(rs042_subset, file = "results/RS042/seurat4_RS042_HSC.rds")


#************************************************************* identifying most and least activated clusters
rs042_subset <- readRDS(file = "results/RS042/seurat4_RS042_HSC.rds")
DimPlot(rs042_subset,label = T)
DimPlot(rs042_subset, reduction = "umap",group.by = "Cell.types", label = T)
DimPlot(rs042_subset, reduction = "umap",group.by = "Cell.clusters", label = T)
FeaturePlot(rs042_subset, "monocle3_pseudotime") + scale_color_viridis_c(option = "plasma")
DotPlot(rs042_subset, features = markerGenes) + RotatedAxis()

FeaturePlot(object = rs042_subset, features = c("percent.mt"))
FeaturePlot(object = rs042_subset, features = c("percent.ribo"))
FeaturePlot(object = rs042_subset, features = c("nCount_RNA"))
FeaturePlot(object = rs042_subset, features = c("nFeature_RNA"))

# myoFib signature
sig <- unlist(strsplit("Acta2, Col1a1, Lox, Timp1",split = ", "))

rs042_subset <- AddModuleScore(rs042_subset,features = list(sig), name = "Myofib")
FeaturePlot(rs042_subset, "Myofib1") + scale_color_viridis_c(option = "viridis")
VlnPlot(rs042_subset, "Myofib1", sort = T, group.by = "Cell.clusters.HSC")

clust_my <- 4
clust_cy <- 2
# compute avg scre for the clusters
tmppos <- which(rs042_subset@meta.data$Cell.clusters.HSC==clust_my)
score_my <- mean(rs042_subset@meta.data$Myofib1[tmppos])
tmppos <- which(rs042_subset@meta.data$Cell.clusters.HSC==clust_cy)
score_cy <- mean(rs042_subset@meta.data$Myofib1[tmppos])
avg_myofibScore_rs042 <- data.frame(score_my, score_cy)
colnames(avg_myofibScore_rs042) <- c(paste0("Clust_",clust_my),paste0("Clust_",clust_cy))
avg_myofibScore_rs042

# get the cy-my DEGs
DimPlot(rs042_subset, reduction = "umap", label = T)
DEGs <- FindMarkers(rs042_subset, ident.1 = c(clust_my), ident.2 = clust_cy)
DEGs <- DEGs[order(DEGs$avg_log2FC, decreasing = T),]
head(DEGs,20)
tail(DEGs,20)
# write.csv(DEGs, file = "results/RS042/DEGs_cymy.csv")

# getting marker genes
marker <- FindAllMarkers(rs042_subset,only.pos = T)
marker <- marker[order(marker$cluster, marker$avg_log2FC, decreasing = T),]
head(marker,20)
tail(marker,20)
# write.csv(marker, file = "results/RS042/markers_HSC.csv")



# ************************************************** cy-my signature analysis
# this part is performed after identifying the cy-my signature using script 'createSignature.R'
library(readxl)
dataset <- data.frame(read_excel("results/cymySig_secreated.xlsx", skip = 1))
head(dataset)

cySig <- na.omit(dataset$cySig2)
mySig <- na.omit(dataset$mySgi2)

sig <- cySig
sigName <- "cySig"
callName <- paste0(sigName,"1")
rs042_subset <- AddModuleScore(rs042_subset,features = list(sig), name = sigName)
FeaturePlot(rs042_subset, callName, max.cutoff = "q99") + scale_color_viridis_c(option = "viridis")
FeaturePlot(rs042_subset, callName,pt.size = ptsize,order = T) + scale_color_viridis_c(option = "viridis")+labs(title = "RS042 TAZFPC cy signature")

sig <- mySig
sigName <- "mySig"
callName <- paste0(sigName,"1")
rs042_subset <- AddModuleScore(rs042_subset,features = list(sig), name = sigName)
FeaturePlot(rs042_subset, callName,pt.size = ptsize,order = T) + scale_color_viridis_c(option = "viridis")+labs(title = "RS042 TAZFPC my signature")

VlnPlot(rs042_subset, features = c("mySig1", "cySig1"), sort = F)

# quantify the cy-my score average in different clusters
df <- data.frame(clust=rs042_subset@active.ident, cy=rs042_subset@meta.data$cySig1, my=rs042_subset@meta.data$mySig1)
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

# rs042_subset[["Cell.clusters.HSC"]] <- Idents(object = rs042_subset)
DimPlot(rs042_subset, reduction = "umap",group.by = "Cell.clusters.HSC", label = T)
rs042_subset <- SetIdent(rs042_subset, value = "Cell.clusters.HSC")
rs042_subset <- RenameIdents(object = rs042_subset, 
                             '2' = "cyHSC",
                             '3' = "cyHSC",
                             
                             '0' = "myHSC",
                             '1' = "myHSC",
                             '4' = "myHSC",
                             '5' = "myHSC",
                             '6' = "myHSC")


rs042_subset[["Cell.types.cymy"]] <- Idents(object = rs042_subset)
DimPlot(rs042_subset,label = T, group.by = "Cell.types.cymy")
table(rs042_subset@active.ident)

cols_cymy <- c("#FF9900","#000099")
DimPlot(rs042_subset, reduction = "umap",group.by = "Cell.types.cymy", pt.size = 2,cols = cols_cymy, label = F)


# saveRDS(rs042_subset, file = "results/RS042/seurat4_RS042_HSC.rds")

#****************************************************************write the cell types annotation to file
# get the cell type annotation
colnames(rs042@meta.data)
celltypes <- rs042@meta.data[,c(1,2,3,4,8,9)]
head(celltypes)

colnames(rs042_subset@meta.data)
tmp <- rs042_subset@meta.data[,c(1,2,3,8,9,10,13,14)]
head(tmp)

tmppos <- match(rownames(rs042_subset@meta.data),rownames(rs042@meta.data))
head(celltypes[tmppos,])
head(tmp)


celltypes$Cell.clusters.HSC <- factor(NA, levels=levels(tmp$Cell.clusters.HSC))
celltypes$monocle3_pseudotime <- as.numeric(NA)
celltypes$Cell.types.cymy <- factor(NA, levels=levels(tmp$Cell.types.cymy))
all(rownames(celltypes[tmppos,])==rownames(tmp))
celltypes[tmppos,c("Cell.clusters.HSC", "monocle3_pseudotime", "Cell.types.cymy")] <- tmp[,c("Cell.clusters.HSC", "monocle3_pseudotime", "Cell.types.cymy")]

head(celltypes)
head(celltypes[tmppos,])

write.csv(celltypes, file = "results/RS042/rs042_TAZ_FPC_cellAnnotation.csv")


#****************************************************************saved object
rs042_subset <- readRDS(file = "results/RS042/seurat4_RS042_HSC.rds")
DimPlot(rs042_subset,label = T, pt.size = 1.5)
DimPlot(rs042_subset, reduction = "umap",group.by = "Cell.types", label = T)
DimPlot(rs042_subset, reduction = "umap",group.by = "Cell.clusters", label = T)
DimPlot(rs042_subset, reduction = "umap",group.by = "Cell.clusters.HSC", label = T)
DimPlot(rs042_subset, reduction = "umap",group.by = "Cell.types.cymy", label = T)
FeaturePlot(rs042_subset, "monocle3_pseudotime") + scale_color_viridis_c(option = "plasma")
DotPlot(rs042_subset, features = markerGenes) + RotatedAxis()



#******************************************************#******************************************************
#******************************************************PLOTS
#******************************************************#******************************************************
library(RColorBrewer)
library("ggsci")
library("scales")

# original data figs
rs042 <- readRDS(file = "results/RS042/seurat4_RS042.rds")
DimPlot(rs042,label = T)
DimPlot(rs042, reduction = "umap",group.by = "Cell.types", label = T)
DimPlot(rs042, reduction = "umap",group.by = "Cell.clusters", label = T)
ggsave(filename = paste0("results/RS042/figs/mus_RS042_TAZFPC_all_clusters_label.pdf"), width = 6, height = 6)
DimPlot(rs042, reduction = "umap",group.by = "Cell.clusters", label = F)
ggsave(filename = paste0("results/RS042/figs/mus_RS042_TAZFPC_all_clusters_noLabel.pdf"), width = 6, height = 6)
DimPlot(rs042, reduction = "umap",group.by = "Cell.types", label = T)
ggsave(filename = paste0("results/RS042/figs/mus_RS042_TAZFPC_all_cellTypes_label.pdf"), width = 6, height = 6)

markerGenes <- c("Lrat","Colec11","Des", #"Col1a1", #"Lox","Col3a1",
                 "Dcn","Pdgfra","Pdgfrb",
                 "Myh11","Actg2", #"Rgs5","Acta2",
                 # "Dpt","Gpx3", #,"Gsn", #"Mgp","Eln","Mfap4"
                 'Upk3b','Msln',
                 "Kdr","Aqp1",
                 # "Il7r","Cd3d","Nkg7", #,"Trac"
                 "Ctss","Ptprc", #"Trem2",
                 # "Clec4f","Vsig4",
                 # "Lgals3","Ccr2","Bpgm",
                 # "Cd79a","Ebf1","Jchain","Siglech",
                 'Epcam', 'Krt7',
                 "Alb","Serpina1a",
                 "Top2a","Mki67")
# 'Epcam', 'Krt7',
# "Alb","Serpina1a")
DotPlot(rs042, features = markerGenes, group.by = "Cell.clusters", cluster.idents = F) + RotatedAxis()+ scale_color_viridis_c() #, cols = "RdBu"
ggsave(filename = paste0("results/RS042/figs/mus_RS042_TAZFPC_all_clusters_dotplot.pdf"))#, width = 6, height = 6)

# CELL TYPES UMAP
data_req <- rs042_subset
DimPlot(object = data_req, group.by = "ident", label = T)
# levels(data_req@active.ident)
ptsize <- 2

DimPlot(object = data_req, group.by = "Cell.clusters.HSC", label = F, pt.size = ptsize)
ggsave(filename = paste0("results/RS042/figs/mus_RS042_TAZFPC_HSC_clusters.pdf"))
DimPlot(object = data_req, group.by = "Cell.clusters.HSC", label = T, pt.size = ptsize)
ggsave(filename = paste0("results/RS042/figs/mus_RS042_TAZFPC_HSC_clusters_label.pdf"))

# cy-my cluster selection
sig <- unlist(strsplit("Acta2, Col1a1, Lox, Timp1",split = ", "))
sigName <- "Myofib"
callName <- paste0(sigName,"1")
data_req <- AddModuleScore(data_req,features = list(sig), name = sigName)
p1 <- DimPlot(object = data_req, group.by = "Cell.clusters.HSC", label = T, pt.size = ptsize)
p2 <- VlnPlot(data_req, group.by = "Cell.clusters.HSC", callName)
wrap_plots(p1,p2)
ggsave(filename = paste0("results/RS042/figs/mus_RS042_TAZFPC_HSC_clusters_cymy_clusters.pdf"),width = 12, height = 6)



genes <- c(unlist(strsplit("Lrat, Col1a1, Acta2, Hgf",split = ", ")))
# markers required
genes <- intersect(genes, rownames(data_req@assays$RNA@data))
FeaturePlot(object = data_req, features = c(genes[1]), max.cutoff = "q99", pt.size = ptsize,order = T) + scale_color_viridis_c()

for(gene in genes[4]){
  FeaturePlot(object = data_req, features = c(gene), max.cutoff = "q99", pt.size = ptsize,order = T) + scale_color_viridis_c()
  ggsave(filename = paste0("results/RS042/figs/mus_RS042_TAZFPC_HSC_gene_",gene,".pdf"))
}



# pseudotime
FeaturePlot(data_req, "monocle3_pseudotime", pt.size = ptsize) + scale_color_viridis_c(option = "plasma") + labs(title = "Mouse RS042 TAZ FPC")
ggsave(filename = paste0("results/RS042/figs/mus_RS042_TAZFPC_HSC_pseudotime.pdf"))
#
#saved monocle object
cds <- readRDS(file = "results/RS042/monocle_RS042_HSC.rds")
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)
ggsave(filename = paste0("results/RS042/figs/mus_RS042_TAZFPC_HSC_pseudotime_monocle.pdf"))

# # cy-my signature
# library(readxl)
# dataset <- data.frame(read_excel("results/cymySig_secreated.xlsx", skip = 1))
# head(dataset)
# 
# cySecrSig <- na.omit(dataset$cySig2_secreted)
# mySecrSig <- na.omit(dataset$mySgi2_secreted)
# 
# sig <- cySecrSig
# sigName <- "cySecrSig"
# callName <- paste0(sigName,"1")
# data_req <- AddModuleScore(data_req,features = list(sig), name = sigName)
# FeaturePlot(data_req, callName, max.cutoff = "q99",pt.size = ptsize,order = T) + scale_color_viridis_c(option = "viridis") +labs(title = "Mouse RS042 TAZ FPC cy-secreted signature")
# ggsave(filename = paste0("results/RS042/figs/mus_RS042_TAZFPC_HSC_signature_cy_secreted.pdf"))
# 
# # VlnPlot(data_req, callName, group.by = "Cell.types", sort = T)
# # VlnPlot(data_req, callName, group.by = "Cell.clusters", sort = T)
# 
# sig <- mySecrSig
# sigName <- "mySecrSig"
# callName <- paste0(sigName,"1")
# data_req <- AddModuleScore(data_req,features = list(sig), name = sigName)
# FeaturePlot(data_req, callName, max.cutoff = "q99",pt.size = ptsize,order = T) + scale_color_viridis_c(option = "viridis") +labs(title = "Mouse RS042 TAZ FPC my-secreted signature")
# ggsave(filename = paste0("results/RS042/figs/mus_RS042_TAZFPC_HSC_signature_my_secreted.pdf"))
# 
# #

# cy-my signature

FeaturePlot(data_req, "cySig1", max.cutoff = "q99",pt.size = ptsize,order = T) + scale_color_viridis_c(option = "viridis") +labs(title = "Mouse RS042 TAZFPC cy signature")
ggsave(filename = paste0("results/RS042/figs/mus_RS042_TAZFPC_HSC_noPF_signature_cy.pdf"))

FeaturePlot(data_req, callName, max.cutoff = "q99",pt.size = ptsize,order = T) + scale_color_viridis_c(option = "viridis") +labs(title = "Mouse RS042 TAZFPC my signature")
ggsave(filename = paste0("results/RS042/figs/mus_RS042_TAZFPC_HSC_noPF_signature_my.pdf"))

cols_cymy <- c("#FF9900","#000099")
DimPlot(rs042_subset, reduction = "umap",group.by = "Cell.types.cymy", pt.size = 2,cols = cols_cymy, label = F)
ggsave(filename = paste0("results/RS042/figs/mus_RS042_TAZFPC_HSC_noPF_signature_clust_cymy.pdf"))


# violin plots
gene <- "Col1a1"
p <- VlnPlot(data_req, features = gene, group.by = "Cell.types.cymy", cols = cols_cymy) #+geom_boxplot()
df <- p$data
colnames(df)[1] <- "gene"
head(df)
tmppos <- which(df$ident=="cyHSC")
tmpx <- df$gene[tmppos]
tmpy <- df$gene[-tmppos]
tmp <- wilcox.test(x=tmpx, y=tmpy)
tmp$p.value

ggplot(data = df, aes(x=ident, y=gene)) +
  geom_violin(aes(fill=ident)) +
  geom_boxplot(width=0.1) +scale_fill_manual(values = cols_cymy) +
  theme_classic()+
  ylab(gene)+
  ggtitle(paste0("RS042 ",gene," Mann Whitney p ",scientific(tmp$p.value,digits = 3)))
ggsave(filename = paste0("results/RS042/figs/mus_RS042_TAZFPC_HSC_noPF_signature_clust_cymy_vln_",gene,".pdf"),width = 5)


gene <- "Hgf" 
p <- VlnPlot(data_req, features = gene, group.by = "Cell.types.cymy", cols = cols_cymy) #+geom_boxplot()
df <- p$data
colnames(df)[1] <- "gene"
head(df)
tmppos <- which(df$ident=="cyHSC")
tmpx <- df$gene[tmppos]
tmpy <- df$gene[-tmppos]
tmp <- wilcox.test(x=tmpx, y=tmpy)
tmp$p.value

ggplot(data = df, aes(x=ident, y=gene)) +
  geom_violin(aes(fill=ident)) +
  geom_boxplot(width=0.1) +scale_fill_manual(values = cols_cymy) +
  theme_classic()+
  ylab(gene)+
  ggtitle(paste0("RS042 ",gene," Mann Whitney p ",scientific(tmp$p.value,digits = 3)))
ggsave(filename = paste0("results/RS042/figs/mus_RS042_TAZFPC_HSC_noPF_signature_clust_cymy_vln_",gene,".pdf"),width = 5)

# get the cy my correlation plots
library(ggpubr)
df <- data.frame(cyHSC=rs042_subset@meta.data$cySig1, myHSC=rs042_subset@meta.data$mySig1, score=rs042_subset@meta.data$mySig1-rs042_subset@meta.data$cySig1)
head(df)
ggscatter(df, x = "cyHSC", y = "myHSC",
          add = "reg.line",                                 # Add regression line
          conf.int = TRUE,                                  # Add confidence interval
          add.params = list(color = "black", fill = "lightgray"),
          color = "score")+
  stat_cor(method = "pearson", label.x = 1, label.y = 1.2)+  # Add correlation coefficient
  gradient_color(c("blue", "yellow", "red"))
ggsave(filename = paste0("results/RS042/figs/mus_RS042_TAZFPC_HSC_noPF_signature_clust_cymy_correln.pdf"))

