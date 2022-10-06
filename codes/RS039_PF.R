# 2021 March 23th created
# this is the Seurat analysis of RS039

setwd("D:/pCloud Sync/Columbia/SchwabeLab/Aveline_HSCnHCC/github/HSCinHCC/")

library(Seurat)
library(ggplot2)
library(monocle3)
library(SeuratWrappers)
library(patchwork)
library(magrittr)
library(enrichR)

markerGenes <- c("Lrat","Myh11","Actg2","Rgs5","Acta2","Lox","Col3a1","Col1a1","Mgp","Eln","Mfap4","Dpt","Gsn","Kdr","Aqp1","Top2a","Mki67",'Upk3b','Msln', 'Epcam', 'Krt7',"Il7r","Cd3d","Trac","Nkg7","Clec4f","Vsig4","Lgals3","Ccr2","Bpgm","Cd79a","Ebf1","Jchain","Siglech","Alb","Serpina1a")

data <- Read10X_h5("data/RS039_HSC_PF/filtered_feature_bc_matrix_RS039.h5")#Aveline
rs039 <- CreateSeuratObject(counts = data, project = "RS039 TAZ_FPC", min.cells = 3, min.features = 200)
rs039
rs039[["percent.mt"]] <- PercentageFeatureSet(rs039, pattern = "^mt-")
rs039[["percent.ribo"]] <- PercentageFeatureSet(rs039, pattern = "^Rpl|^Rps|Mrps|^Mrpl")#

# Visualize QC metrics as a violin plot
VlnPlot(rs039, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(rs039, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(rs039, feature1 = "nFeature_RNA", feature2 = "percent.mt")
FeatureScatter(rs039, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# keep only the required cells
rs039 <- subset(rs039, subset = nFeature_RNA > 200 & nCount_RNA < 35000 & percent.mt < 10)

VlnPlot(rs039, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(rs039, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(rs039, feature1 = "nFeature_RNA", feature2 = "percent.mt")
FeatureScatter(rs039, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")


rs039 <- NormalizeData(rs039, normalization.method = "LogNormalize", scale.factor = 10000)
rs039 <- FindVariableFeatures(rs039, selection.method = "vst", nfeatures = 2000)
# all.genes <- rownames(rs039)
rs039 <- ScaleData(rs039, features = VariableFeatures(object = rs039), vars.to.regress = c("nCount_RNA", "percent.mt"))
rs039 <- RunPCA(rs039, features = VariableFeatures(object = rs039))
ElbowPlot(rs039, ndims = 50)
ndim=17
rs039 <- FindNeighbors(rs039, dims = 1:ndim)
rs039 <- FindClusters(rs039, resolution = 0.5)
rs039 <- RunUMAP(rs039, reduction = "pca", dims = 1:ndim)
DimPlot(rs039,label = T)
DotPlot(rs039, features = markerGenes) + RotatedAxis()

# getting marker genes
marker <- FindAllMarkers(rs039,only.pos = T)
# write.csv(marker, file = "results/RS039/markers.csv")

FeaturePlot(object = rs039, features = c("percent.mt"))
FeaturePlot(object = rs039, features = c("percent.ribo"))
FeaturePlot(object = rs039, features = c("nCount_RNA"))
FeaturePlot(object = rs039, features = c("nFeature_RNA"))
VlnPlot(object = rs039, features = c("nCount_RNA"))
VlnPlot(object = rs039, features = c("nFeature_RNA"))

# giving cell types
DimPlot(rs039,label = T)
rs039[["Cell.clusters"]] <- Idents(object = rs039)
rs039 <- SetIdent(rs039, value = "Cell.clusters")
rs039 <- RenameIdents(object = rs039, 
                      '0' = "HSC",
                      '1' = "HSC",
                      '2' = "HSC",
                      '3' = "HSC",
                      '4' = "HSC",
                      '5' = "HSC",
                      '9' = "HSC",
                      '6' = "Cycling/Mixed",
                      '10' = "Mesothelial",
                      '8' = "Endothelial",
                      '11' = "Immune",
                      '7' = "Mixed/Low-quality")

rs039[["Cell.types"]] <- Idents(object = rs039)
DimPlot(rs039,label = T, group.by = "Cell.types")
table(rs039@meta.data$Cell.types)
table(rs039@active.ident)

DimPlot(rs039, reduction = "umap",group.by = "Cell.types", label = F)
DimPlot(rs039, reduction = "umap",group.by = "Cell.clusters", label = F)



# saveRDS(rs039, file = "results/RS039/seurat4_RS039.rds")
#****************************************************************saved object
rs039 <- readRDS(file = "results/RS039/seurat4_RS039.rds")
DimPlot(rs039,label = T)
DimPlot(rs039, reduction = "umap",group.by = "Cell.types", label = T)
DimPlot(rs039, reduction = "umap",group.by = "Cell.clusters", label = T)

# Recluster the HSC population 
rs039 <- SetIdent(rs039, value = "Cell.clusters")
rs039_subset <- subset(rs039, idents = c(0:5,9))
DimPlot(rs039_subset,label = T, group.by = "ident")
rs039_subset <- FindVariableFeatures(rs039_subset, selection.method = "vst", nfeatures = 2000)
# all.genes <- rownames(rs039_subset)
rs039_subset <- ScaleData(rs039_subset, features = VariableFeatures(object = rs039_subset), vars.to.regress = c("nCount_RNA", "percent.mt"))
rs039_subset <- RunPCA(rs039_subset, features = VariableFeatures(object = rs039_subset))
ElbowPlot(rs039_subset, ndims = 50)
ndim=10
rs039_subset <- FindNeighbors(rs039_subset, dims = 1:ndim)
rs039_subset <- FindClusters(rs039_subset, resolution = 0.5)
rs039_subset <- RunUMAP(rs039_subset, reduction = "pca", dims = 1:ndim)
DimPlot(rs039_subset, reduction = "umap",label = T)
FeaturePlot(object = rs039_subset, features = c("Col1a1")) & scale_color_viridis_c()

# saveRDS(rs039_subset, file = "results/RS039/seurat4_RS039_HSC.rds")
#****************************************************************saved object
rs039_subset <- readRDS(file = "results/RS039/seurat4_RS039_HSC.rds")
DimPlot(rs039_subset,label = T)

# ************************************************** trajectory analysis
# getting to Monocle
cds <- as.cell_data_set(rs039_subset)
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

# saveRDS(cds, file = "results/RS039/monocle_RS039_HSC.rds")
#saved object
cds <- readRDS(file = "results/RS039/monocle_RS039_HSC.rds")
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)

# Set the assay back as 'integrated'
integrated.sub <- as.Seurat(cds, assay = NULL,project = "cell_data_set")
FeaturePlot(integrated.sub, "monocle3_pseudotime")

rs039_subset@meta.data$monocle3_pseudotime <- integrated.sub@meta.data$monocle3_pseudotime
# saveRDS(rs039_subset, file = "results/RS039/seurat4_RS039_HSC.rds")


#************************************************************* identifying most and least activated clusters
rs039_subset <- readRDS(file = "results/RS039/seurat4_RS039_HSC.rds")
DimPlot(rs039_subset,label = T)
DimPlot(rs039_subset, reduction = "umap",group.by = "Cell.types", label = T)
DimPlot(rs039_subset, reduction = "umap",group.by = "Cell.clusters", label = T)
FeaturePlot(rs039_subset, "monocle3_pseudotime") + scale_color_viridis_c(option = "plasma")
DotPlot(rs039_subset, features = markerGenes) + RotatedAxis()

FeaturePlot(object = rs039_subset, features = c("percent.mt"))
FeaturePlot(object = rs039_subset, features = c("percent.ribo"))
FeaturePlot(object = rs039_subset, features = c("nCount_RNA"))
FeaturePlot(object = rs039_subset, features = c("nFeature_RNA"))

# myoFib signature
sig <- unlist(strsplit("Acta2, Col1a1, Lox, Timp1",split = ", "))

rs039_subset <- AddModuleScore(rs039_subset,features = list(sig), name = "Myofib")
FeaturePlot(rs039_subset, "Myofib1") + scale_color_viridis_c(option = "viridis")
VlnPlot(rs039_subset, "Myofib1", sort = T, group.by = "Cell.clusters.HSC")

clust_my <- 5
clust_cy <- 3
# compute avg scre for the clusters
tmppos <- which(rs039_subset@meta.data$Cell.clusters.HSC==clust_my)
score_my <- mean(rs039_subset@meta.data$Myofib1[tmppos])
tmppos <- which(rs039_subset@meta.data$Cell.clusters.HSC==clust_cy)
score_cy <- mean(rs039_subset@meta.data$Myofib1[tmppos])
avg_myofibScore_rs039 <- data.frame(score_my, score_cy)
colnames(avg_myofibScore_rs039) <- c(paste0("Clust_",clust_my),paste0("Clust_",clust_cy))
avg_myofibScore_rs039

# get the cy-my DEGs
DimPlot(rs039_subset, reduction = "umap", label = T)
DEGs <- FindMarkers(rs039_subset, ident.1 = c(clust_my), ident.2 = clust_cy)
DEGs <- DEGs[order(DEGs$avg_log2FC, decreasing = T),]
head(DEGs,20)
tail(DEGs,20)
# write.csv(DEGs, file = "results/RS039/DEGs_cymy.csv")

# getting marker genes
marker <- FindAllMarkers(rs039_subset,only.pos = T)
marker <- marker[order(marker$cluster, marker$avg_log2FC, decreasing = T),]
head(marker,20)
tail(marker,20)
# write.csv(marker, file = "results/RS039/markers_HSC.csv")


# ************************************************** cy-my signature
# this part is performed after identifying the cy-my signature using script 'createSignature.R'
library(readxl)
dataset <- data.frame(read_excel("results/cymySig_secreated.xlsx", skip = 1))
head(dataset)

cySig <- na.omit(dataset$cySig2)
mySig <- na.omit(dataset$mySgi2)

ptsize <- 1

sig <- cySig
sigName <- "cySig"
callName <- paste0(sigName,"1")
rs039_subset <- AddModuleScore(rs039_subset,features = list(sig), name = sigName)
FeaturePlot(rs039_subset, callName,pt.size = ptsize,order = T) + scale_color_viridis_c(option = "viridis") +labs(title = "RS039 PF cy signature")


sig <- mySig
sigName <- "mySig"
callName <- paste0(sigName,"1")
rs039_subset <- AddModuleScore(rs039_subset,features = list(sig), name = sigName)
FeaturePlot(rs039_subset, callName,pt.size = ptsize,order = T) + scale_color_viridis_c(option = "viridis")+labs(title = "RS039 PF my signature")

VlnPlot(rs039_subset, features = c("mySig1", "cySig1"), sort = F)

# quantify the cy-my score median in different clusters
df <- data.frame(clust=rs039_subset@active.ident, cy=rs039_subset@meta.data$cySig1, my=rs039_subset@meta.data$mySig1)
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

# rs039_subset[["Cell.clusters.HSC"]] <- Idents(object = rs039_subset)
# DimPlot(rs039_subset, reduction = "umap",group.by = "Cell.clusters.HSC", label = T)
rs039_subset <- SetIdent(rs039_subset, value = "Cell.clusters.HSC")
rs039_subset <- RenameIdents(object = rs039_subset,
                             '1' = "cyHSC",
                             '3' = "cyHSC",
                             '6' = "cyHSC",
                             
                             '0' = "myHSC",
                             '2' = "myHSC",
                             '4' = "myHSC",
                             '5' = "myHSC",
                             '7' = "myHSC")


rs039_subset[["Cell.types.cymy"]] <- Idents(object = rs039_subset)
DimPlot(rs039_subset,label = T, group.by = "Cell.types.cymy")
table(rs039_subset@active.ident)

cols_cymy <- c("#FF9900","#000099")
DimPlot(rs039_subset, reduction = "umap",group.by = "Cell.types.cymy", pt.size = 2,cols = cols_cymy, label = F)

# saveRDS(rs039_subset, file = "results/RS039/seurat4_RS039_HSC.rds")


#****************************************************************write the cell types annotation to file
# get the cell type annotation
colnames(rs039@meta.data)
celltypes <- rs039@meta.data[,c(1,2,3,4,8,9)]
head(celltypes)

colnames(rs039_subset@meta.data)
tmp <- rs039_subset@meta.data[,c(1,2,3,8,9,10,13,14)]
head(tmp)

tmppos <- match(rownames(rs039_subset@meta.data),rownames(rs039@meta.data))
head(celltypes[tmppos,])
head(tmp)


celltypes$Cell.clusters.HSC <- factor(NA, levels=levels(tmp$Cell.clusters.HSC))
celltypes$monocle3_pseudotime <- as.numeric(NA)
celltypes$Cell.types.cymy <- factor(NA, levels=levels(tmp$Cell.types.cymy))
all(rownames(celltypes[tmppos,])==rownames(tmp))
celltypes[tmppos,c("Cell.clusters.HSC", "monocle3_pseudotime", "Cell.types.cymy")] <- tmp[,c("Cell.clusters.HSC", "monocle3_pseudotime", "Cell.types.cymy")]

head(celltypes)
head(celltypes[tmppos,])

# write.csv(celltypes, file = "results/RS039/rs039_PF_cellAnnotation.csv")


#****************************************************************saved object
rs039_subset <- readRDS(file = "results/RS039/seurat4_RS039_HSC.rds")
DimPlot(rs039_subset,label = T)
DimPlot(rs039_subset, reduction = "umap",group.by = "Cell.types", label = T)
DimPlot(rs039_subset, reduction = "umap",group.by = "Cell.clusters", label = T)
DimPlot(rs039_subset, reduction = "umap",group.by = "Cell.clusters.HSC", label = T)
DimPlot(rs039_subset, reduction = "umap",group.by = "Cell.types.cymy", label = T)
FeaturePlot(rs039_subset, "monocle3_pseudotime") + scale_color_viridis_c(option = "plasma")
DotPlot(rs039_subset, features = markerGenes) + RotatedAxis()



#******************************************************#******************************************************
#******************************************************PLOTS
#******************************************************#******************************************************
library(RColorBrewer)
library("ggsci")
library("scales")
library(patchwork)

# original data figs
rs039 <- readRDS(file = "results/RS039/seurat4_RS039.rds")
DimPlot(rs039,label = T)
DimPlot(rs039, reduction = "umap",group.by = "Cell.types", label = T)
DimPlot(rs039, reduction = "umap",group.by = "Cell.clusters", label = T)
ggsave(filename = paste0("results/RS039/figs/mus_RS039_PF_all_clusters_label.pdf"), width = 6, height = 6)
DimPlot(rs039, reduction = "umap",group.by = "Cell.clusters", label = F)
ggsave(filename = paste0("results/RS039/figs/mus_RS039_PF_all_clusters_noLabel.pdf"), width = 6, height = 6)
DimPlot(rs039, reduction = "umap",group.by = "Cell.types", label = T)
ggsave(filename = paste0("results/RS039/figs/mus_RS039_PF_all_cellTypes_label.pdf"), width = 6, height = 6)

markerGenes <- c("Lrat","Colec11","Des", #"Col1a1", #"Lox","Col3a1",
                 "Dcn","Pdgfra","Pdgfrb",
                 # "Myh11","Actg2", #"Rgs5","Acta2",
                 # "Dpt","Gpx3", #,"Gsn", #"Mgp","Eln","Mfap4"
                 'Upk3b','Msln',
                 "Kdr","Aqp1",
                 # "Il7r","Cd3d","Nkg7", #,"Trac"
                 "Ctss","Ptprc", #"Trem2",
                 # "Clec4f","Vsig4",
                 # "Lgals3","Ccr2","Bpgm",
                 # "Cd79a","Ebf1","Jchain","Siglech",
                 "Top2a","Mki67")
                 # 'Epcam', 'Krt7',
                 # "Alb","Serpina1a")
DotPlot(rs039, features = markerGenes, group.by = "Cell.clusters", cluster.idents = F) + RotatedAxis()+ scale_color_viridis_c() #, cols = "RdBu"
ggsave(filename = paste0("results/RS039/figs/mus_RS039_PF_all_clusters_dotplot.pdf"))#, width = 6, height = 6)




# CELL TYPES UMAP
data_req <- rs039_subset
DimPlot(object = data_req, group.by = "ident", label = T)
DimPlot(rs039_subset, reduction = "umap",group.by = "Cell.clusters.HSC", label = T)
# levels(data_req@active.ident)
ptsize <- 0.5


DimPlot(object = data_req, group.by = "Cell.clusters.HSC", label = F, pt.size = ptsize)
ggsave(filename = paste0("results/RS039/figs/mus_RS039_PF_HSC_clusters.pdf"))
DimPlot(object = data_req, group.by = "Cell.clusters.HSC", label = T, pt.size = ptsize)
ggsave(filename = paste0("results/RS039/figs/mus_RS039_PF_HSC_clusters_label.pdf"))

# cy-my cluster selection
sig <- unlist(strsplit("Acta2, Col1a1, Lox, Timp1",split = ", "))
sigName <- "Myofib"
callName <- paste0(sigName,"1")
data_req <- AddModuleScore(data_req,features = list(sig), name = sigName)
p1 <- DimPlot(object = data_req, group.by = "Cell.clusters.HSC", label = T, pt.size = ptsize)
p2 <- VlnPlot(data_req, group.by = "Cell.clusters.HSC", callName)
wrap_plots(p1,p2)
ggsave(filename = paste0("results/RS039/figs/mus_RS039_PF_HSC_clusters_cymy_clusters.pdf"),width = 12, height = 6)

# 
genes <- c(unlist(strsplit("Lrat, Col1a1, Acta2, Hgf",split = ", ")))
# markers required
genes <- intersect(genes, rownames(data_req@assays$RNA@data))
FeaturePlot(object = data_req, features = c(genes[1]), max.cutoff = "q99", pt.size = ptsize,order = T) + scale_color_viridis_c()

for(gene in genes){
  FeaturePlot(object = data_req, features = c(gene), max.cutoff = "q99", pt.size = ptsize,order = T) + scale_color_viridis_c()
  ggsave(filename = paste0("results/RS039/figs/mus_RS039_PF_HSC_gene_",gene,".pdf"))
}

gene <- "Has2"
FeaturePlot(object = data_req, features = c(gene), max.cutoff = "q99", pt.size = 2,order = T) + scale_color_viridis_c()
ggsave(filename = paste0("results/RS039/figs/mus_RS039_PF_HSC_gene_",gene,".pdf"))

# pseudotime
FeaturePlot(data_req, "monocle3_pseudotime", pt.size = ptsize) + scale_color_viridis_c(option = "plasma") + labs(title = "Mouse RS039 PF")
ggsave(filename = paste0("results/RS039/figs/mus_RS039_PF_HSC_pseudotime.pdf"))

#saved monocle object
cds <- readRDS(file = "results/RS039/monocle_RS039_HSC.rds")
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           cell_size = 0.35, #0.35
           graph_label_size=1.5)
ggsave(filename = paste0("results/RS039/figs/mus_RS039_PF_HSC_pseudotime_monocle.pdf"))


# cy-my signature
cols_cymy <- c("#FF9900","#000099")

FeaturePlot(data_req, "cySig1", max.cutoff = "q99",pt.size = ptsize,order = T) + scale_color_viridis_c(option = "viridis") +labs(title = "Mouse RS039 PF cy signature")
ggsave(filename = paste0("results/RS039/figs/mus_RS039_PF_HSC_signature_cy.pdf"))

FeaturePlot(data_req, callName, max.cutoff = "q99",pt.size = ptsize,order = T) + scale_color_viridis_c(option = "viridis") +labs(title = "Mouse RS039 PF my signature")
ggsave(filename = paste0("results/RS039/figs/mus_RS039_PF_HSC_signature_my.pdf"))


DimPlot(rs039_subset, reduction = "umap",group.by = "Cell.types.cymy", pt.size = ptsize,cols = cols_cymy, label = F)
ggsave(filename = paste0("results/RS039/figs/mus_RS039_PF_HSC_signature_clust_cymy.pdf"))


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
  ggtitle(paste0("RS039 ",gene," Mann Whitney p ",scientific(tmp$p.value,digits = 3)))
ggsave(filename = paste0("results/RS039/figs/mus_RS039_PF_HSC_signature_clust_cymy_vln_",gene,".pdf"),width = 5)


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
  ggtitle(paste0("RS039 ",gene," Mann Whitney p ",scientific(tmp$p.value,digits = 3)))
ggsave(filename = paste0("results/RS039/figs/mus_RS039_PF_HSC_signature_clust_cymy_vln_",gene,".pdf"),width = 5)

gene <- "Has2" 
p <- VlnPlot(data_req, features = gene, group.by = "Cell.types.cymy", cols = cols_cymy) #+geom_boxplot()
df <- p$data
colnames(df)[1] <- "gene"
head(df)
tmppos <- which(df$ident=="cyHSC")
tmpx <- df$gene[tmppos]
tmpy <- df$gene[-tmppos]
tmp <- wilcox.test(x=tmpx, y=tmpy, alternative = "two.sided")
tmp$p.value
ggplot(data = df, aes(x=ident, y=gene,fill=ident)) +
  geom_dotplot(binaxis = "y",dotsize=0.5,stackdir="center",binwidth = 1/20, drop = T, color=NA)+
  scale_fill_manual(values = cols_cymy) +
  theme_classic()+
  ylab(gene)+
  ggtitle(paste0("RS039 ",gene," Mann Whitney (two.sided) p ",scientific(tmp$p.value,digits = 3)))
ggsave(filename = paste0("results/RS039/figs/mus_RS039_PF_HSC_signature_clust_cymy_vln_",gene,".pdf"),width = 5)

tmp <- wilcox.test(x=tmpx, y=tmpy, alternative = "less")
tmp$p.value
ggplot(data = df, aes(x=ident, y=gene,fill=ident)) +
  geom_dotplot(binaxis = "y",dotsize=0.5,stackdir="center",binwidth = 1/20, drop = T, color=NA)+
  scale_fill_manual(values = cols_cymy) +
  theme_classic()+
  ylab(gene)+
  ggtitle(paste0("RS039 ",gene," Mann Whitney (less) p ",scientific(tmp$p.value,digits = 3)))
ggsave(filename = paste0("results/RS039/figs/mus_RS039_PF_HSC_signature_clust_cymy_vln_",gene,"_1.pdf"),width = 5)

# get the cy my correlation plots
library(ggpubr)
df <- data.frame(cyHSC=rs039_subset@meta.data$cySig1, myHSC=rs039_subset@meta.data$mySig1, score=rs039_subset@meta.data$mySig1-rs039_subset@meta.data$cySig1)
head(df)
ggscatter(df, x = "cyHSC", y = "myHSC",
          add = "reg.line",                                 # Add regression line
          conf.int = TRUE,                                  # Add confidence interval
          add.params = list(color = "black", fill = "lightgray"),
          color = "score")+
  stat_cor(method = "pearson", label.x = 0.75, label.y = 1.2)+  # Add correlation coefficient
  gradient_color(c("blue", "yellow", "red"))
ggsave(filename = paste0("results/RS039/figs/mus_RS039_PF_HSC_signature_clust_cymy_correln.pdf"))


