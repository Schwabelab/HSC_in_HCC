# 2021 March 31th created
# this is the Seurat3 analysis of RS043

setwd("D:/pCloud Sync/Columbia/SchwabeLab/Aveline_HSCnHCC/github/HSCinHCC/")

library(Seurat)
library(ggplot2)
library(monocle3)
library(SeuratWrappers)
library(patchwork)
library(magrittr)

markerGenes <- c("Lrat","Myh11","Actg2","Rgs5","Acta2","Lox","Col3a1","Col1a1","Mgp","Eln","Mfap4","Dpt","Gsn","Kdr","Aqp1","Top2a","Mki67",'Upk3b','Msln', 'Epcam', 'Krt7',"Il7r","Cd3d","Trac","Nkg7","Clec4f","Vsig4","Lgals3","Ccr2","Bpgm","Cd79a","Ebf1","Jchain","Siglech","Alb","Serpina1a")

# 
data <- Read10X_h5("data/RS043_HSC_CDAA/filtered_feature_bc_matrix_RS043.h5")
rs043 <- CreateSeuratObject(counts = data, project = "RS043 CDAA", min.cells = 3, min.features = 200)
rs043
rs043[["percent.mt"]] <- PercentageFeatureSet(rs043, pattern = "^mt-")
rs043[["percent.ribo"]] <- PercentageFeatureSet(rs043, pattern = "^Rpl|^Rps|Mrps|^Mrpl")#

# Visualize QC metrics as a violin plot
VlnPlot(rs043, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# FeatureScatter(rs043, feature1 = "nCount_RNA", feature2 = "percent.mt")
# FeatureScatter(rs043, feature1 = "nFeature_RNA", feature2 = "percent.mt")
# FeatureScatter(rs043, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p1 <- FeatureScatter(rs043, feature1 = "nCount_RNA", feature2 = "percent.mt") + NoLegend()
p2 <- FeatureScatter(rs043, feature1 = "nFeature_RNA", feature2 = "percent.mt") + NoLegend()
p3 <- FeatureScatter(rs043, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# p4 <- FeatureScatter(liv_01_cov, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "percent.mt")
# # keep only the required cells
wrap_plots(list(p1,p2,p3),nrow = 2,ncol = 2)
# keep only the required cells
rs043 <- subset(rs043, subset = nFeature_RNA > 200 & nCount_RNA < 35000 & percent.mt < 10)

VlnPlot(rs043, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# FeatureScatter(rs043, feature1 = "nCount_RNA", feature2 = "percent.mt")
# FeatureScatter(rs043, feature1 = "nFeature_RNA", feature2 = "percent.mt")
# FeatureScatter(rs043, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p1 <- FeatureScatter(rs043, feature1 = "nCount_RNA", feature2 = "percent.mt") + NoLegend()
p2 <- FeatureScatter(rs043, feature1 = "nFeature_RNA", feature2 = "percent.mt") + NoLegend()
p3 <- FeatureScatter(rs043, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# p4 <- FeatureScatter(liv_01_cov, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "percent.mt")
# # keep only the required cells
wrap_plots(list(p1,p2,p3),nrow = 2,ncol = 2)


rs043 <- NormalizeData(rs043, normalization.method = "LogNormalize", scale.factor = 10000)
rs043 <- FindVariableFeatures(rs043, selection.method = "vst", nfeatures = 2000)
# all.genes <- rownames(rs043)
rs043 <- ScaleData(rs043, features = VariableFeatures(object = rs043), vars.to.regress = c("nCount_RNA", "percent.mt"))
rs043 <- RunPCA(rs043, features = VariableFeatures(object = rs043))
ElbowPlot(rs043, ndims = 50)
ndim=27
rs043 <- FindNeighbors(rs043, dims = 1:ndim)
rs043 <- FindClusters(rs043, resolution = 0.5)
rs043 <- RunUMAP(rs043, reduction = "pca", dims = 1:ndim)

DimPlot(rs043,label = T)
markerGenes <- c("Lrat","Myh11","Actg2","Rgs5","Acta2","Lox","Col3a1","Col1a1","Mgp","Eln","Mfap4","Dpt","Gsn","Kdr","Aqp1","Top2a","Mki67",'Upk3b','Msln', 'Epcam', 'Krt7',"Il7r","Cd3d","Trac","Nkg7","Clec4f","Vsig4","Lgals3","Ccr2","Bpgm","Cd79a","Ebf1","Jchain","Siglech","Alb","Serpina1a")
DotPlot(rs043, features = markerGenes) + RotatedAxis()

# getting marker genes
marker <- FindAllMarkers(rs043,only.pos = T)
# write.csv(marker, file = "results/rs043/markers.csv")

FeaturePlot(object = rs043, features = c("percent.mt"))
FeaturePlot(object = rs043, features = c("percent.ribo"))
FeaturePlot(object = rs043, features = c("nCount_RNA"))
FeaturePlot(object = rs043, features = c("nFeature_RNA"))
VlnPlot(object = rs043, features = c("nCount_RNA"))
VlnPlot(object = rs043, features = c("nFeature_RNA"))

# giving cell types
DimPlot(rs043,label = T)
# rs043[["Cell.clusters"]] <- Idents(object = rs043)
rs043 <- SetIdent(rs043, value = "Cell.clusters")
rs043 <- RenameIdents(object = rs043, 
                      '0' = "HSC",
                      '1' = "HSC",
                      '2' = "HSC",
                      '3' = "HSC",
                      '4' = "HSC",
                      '7' = "HSC",
                      
                      '8' = "Cycling/Mixed",
                      '5' = "VSMC",
                      
                      '10' = "Cholangiocyte",
                      '12' = "PF/Mesothelial",
                      '11' = "Hepatocytes",
                      
                      '9' = "Endothelial",
                      # '12' = "Endothelial",
                      
                      '6' = "Immune")
rs043[["Cell.types"]] <- Idents(object = rs043)
DimPlot(rs043,label = T, group.by = "Cell.types")
table(rs043@meta.data$Cell.types)
table(rs043@meta.data$Cell.clusters)

DimPlot(rs043, reduction = "umap",group.by = "Cell.types", label = F)
DimPlot(rs043, reduction = "umap",group.by = "Cell.clusters", label = F)


# saveRDS(rs043, file = "results/RS043/seurat4_RS043.rds")
#****************************************************************saved object
rs043 <- readRDS(file = "results/rs043/seurat4_RS043.rds")
DimPlot(rs043,label = T)
DimPlot(rs043, reduction = "umap",group.by = "Cell.types", label = T)
DimPlot(rs043, reduction = "umap",group.by = "Cell.clusters", label = T)
DotPlot(rs043, features = markerGenes) + RotatedAxis()


# Recluster the HSC population 
rs043 <- SetIdent(rs043, value = "Cell.clusters")
rs043_subset <- subset(rs043, idents = c(0:4,7))
DimPlot(rs043_subset,label = T, group.by = "ident")
rs043_subset <- FindVariableFeatures(rs043_subset, selection.method = "vst", nfeatures = 2000)
# all.genes <- rownames(rs043_subset)
rs043_subset <- ScaleData(rs043_subset, features = VariableFeatures(object = rs043_subset), vars.to.regress = c("nCount_RNA", "percent.mt"))
rs043_subset <- RunPCA(rs043_subset, features = VariableFeatures(object = rs043_subset))
ElbowPlot(rs043_subset, ndims = 50)
ndim=14
rs043_subset <- FindNeighbors(rs043_subset, dims = 1:ndim)
rs043_subset <- FindClusters(rs043_subset, resolution = 0.5)
rs043_subset <- RunUMAP(rs043_subset, reduction = "pca", dims = 1:ndim)
DimPlot(rs043_subset, reduction = "umap",label = T)
DotPlot(rs043_subset, features = markerGenes) + RotatedAxis()
FeaturePlot(object = rs043_subset, features = c("Col1a1")) & scale_color_viridis_c()

# saveRDS(rs043_subset, file = "results/rs043/seurat4_RS043_HSC.rds")

#****************************************************************saved object
rs043_subset <- readRDS(file = "results/rs043/seurat4_RS043_HSC.rds")
DimPlot(rs043_subset,label = T)
DimPlot(rs043_subset, reduction = "umap",group.by = "Cell.types", label = T)
DimPlot(rs043_subset, reduction = "umap",group.by = "Cell.clusters", label = T)


# ************************************************** trajectory analysis
# getting to Monocle
cds <- as.cell_data_set(rs043_subset)
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

# saveRDS(cds, file = "results/RS043/monocle_RS043_HSC.rds")
#saved object
cds <- readRDS(file = "results/RS043/monocle_RS043_HSC.rds")
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

rs043_subset@meta.data$monocle3_pseudotime <- integrated.sub@meta.data$monocle3_pseudotime
FeaturePlot(rs043_subset, "monocle3_pseudotime")


#****************************************************************saved object
rs043_subset <- readRDS(file = "results/RS043/seurat4_RS043_HSC.rds")
DimPlot(rs043_subset,label = T)
DimPlot(rs043_subset, reduction = "umap",group.by = "Cell.types", label = T)
DimPlot(rs043_subset, reduction = "umap",group.by = "Cell.clusters", label = T)
FeaturePlot(rs043_subset, "monocle3_pseudotime") + scale_color_viridis_c(option = "plasma")
DotPlot(rs043_subset, features = markerGenes) + RotatedAxis()

FeaturePlot(object = rs043_subset, features = c("percent.mt"))
FeaturePlot(object = rs043_subset, features = c("percent.ribo"))
FeaturePlot(object = rs043_subset, features = c("nCount_RNA"))
FeaturePlot(object = rs043_subset, features = c("nFeature_RNA"))

# myoFib signature
sig <- unlist(strsplit("Acta2, Col1a1, Lox, Timp1",split = ", "))

rs043_subset <- AddModuleScore(rs043_subset,features = list(sig), name = "Myofib")
FeaturePlot(rs043_subset, "Myofib1") + scale_color_viridis_c(option = "viridis")
VlnPlot(rs043_subset, "Myofib1", sort = T, group.by = "Cell.clusters.HSC")

clust_my <- 1
clust_cy <- 3
# compute avg scre for the clusters
tmppos <- which(rs043_subset@meta.data$Cell.clusters.HSC==clust_my)
score_my <- mean(rs043_subset@meta.data$Myofib1[tmppos])
tmppos <- which(rs043_subset@meta.data$Cell.clusters.HSC==clust_cy)
score_cy <- mean(rs043_subset@meta.data$Myofib1[tmppos])
avg_myofibScore_rs043 <- data.frame(score_my, score_cy)
colnames(avg_myofibScore_rs043) <- c(paste0("Clust_",clust_my),paste0("Clust_",clust_cy))
avg_myofibScore_rs043

# get the cy-my DEGs
DimPlot(rs043_subset, reduction = "umap", label = T)
DEGs <- FindMarkers(rs043_subset, ident.1 = c(clust_my), ident.2 = clust_cy)
DEGs <- DEGs[order(DEGs$avg_log2FC, decreasing = T),]
head(DEGs,20)
tail(DEGs,20)
# write.csv(DEGs, file = "results/RS043/DEGs_cymy.csv")

# getting marker genes
marker <- FindAllMarkers(rs043_subset,only.pos = T)
marker <- marker[order(marker$cluster, marker$avg_log2FC, decreasing = T),]
head(marker,20)
tail(marker,20)
# write.csv(marker, file = "results/RS043/markers_HSC.csv")


# ************************************************** cy-my signature analysis
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
rs043_subset <- AddModuleScore(rs043_subset,features = list(sig), name = sigName)
FeaturePlot(rs043_subset, callName,pt.size = ptsize,order = T) + scale_color_viridis_c(option = "viridis") +labs(title = "RS043 CDAA cy signature")


sig <- mySig
sigName <- "mySig"
callName <- paste0(sigName,"1")
rs043_subset <- AddModuleScore(rs043_subset,features = list(sig), name = sigName)
FeaturePlot(rs043_subset, callName,pt.size = ptsize,order = T) + scale_color_viridis_c(option = "viridis")+labs(title = "RS043 CDAA my signature")

VlnPlot(rs043_subset, features = c("mySig1", "cySig1"), sort = F)

# quantify the cy-my score average in different clusters
df <- data.frame(clust=rs043_subset@active.ident, cy=rs043_subset@meta.data$cySig1, my=rs043_subset@meta.data$mySig1)
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

# rs043_subset[["Cell.clusters.HSC"]] <- Idents(object = rs043_subset)
DimPlot(rs043_subset, reduction = "umap",group.by = "Cell.clusters.HSC", label = T)
rs043_subset <- SetIdent(rs043_subset, value = "Cell.clusters.HSC")
rs043_subset <- RenameIdents(object = rs043_subset,
                             '0' = "cyHSC",
                             '2' = "cyHSC",
                             '3' = "cyHSC",
                             '4' = "cyHSC",
                             
                             '1' = "myHSC")


rs043_subset[["Cell.types.cymy"]] <- Idents(object = rs043_subset)
DimPlot(rs043_subset,label = T, group.by = "Cell.types.cymy")
table(rs043_subset@active.ident)

cols_cymy <- c("#FF9900","#000099")
DimPlot(rs043_subset, reduction = "umap",group.by = "Cell.types.cymy", pt.size = 2,cols = cols_cymy, label = F)



# saveRDS(rs043_subset, file = "results/RS043/seurat4_RS043_HSC.rds")

#****************************************************************write the cell types annotation to file
# get the cell type annotation
colnames(rs043@meta.data)
celltypes <- rs043@meta.data[,c(1,2,3,4,8,9)]
head(celltypes)

colnames(rs043_subset@meta.data)
tmp <- rs043_subset@meta.data[,c(1,2,3,8,9,10,13,14)]
head(tmp)

tmppos <- match(rownames(rs043_subset@meta.data),rownames(rs043@meta.data))
head(celltypes[tmppos,])
head(tmp)


celltypes$Cell.clusters.HSC <- factor(NA, levels=levels(tmp$Cell.clusters.HSC))
celltypes$monocle3_pseudotime <- as.numeric(NA)
celltypes$Cell.types.cymy <- factor(NA, levels=levels(tmp$Cell.types.cymy))
all(rownames(celltypes[tmppos,])==rownames(tmp))
celltypes[tmppos,c("Cell.clusters.HSC", "monocle3_pseudotime", "Cell.types.cymy")] <- tmp[,c("Cell.clusters.HSC", "monocle3_pseudotime", "Cell.types.cymy")]

head(celltypes)
head(celltypes[tmppos,])

# write.csv(celltypes, file = "results/RS043/rs043_CDAA_cellAnnotation.csv")



#****************************************************************saved object
rs043_subset <- readRDS(file = "results/RS043/seurat4_RS043_HSC.rds")
DimPlot(rs043_subset,label = T)
DimPlot(rs043_subset, reduction = "umap",group.by = "Cell.types", label = T)
DimPlot(rs043_subset, reduction = "umap",group.by = "Cell.clusters", label = T)
DimPlot(rs043_subset, reduction = "umap",group.by = "Cell.clusters.HSC", label = T)
DimPlot(rs043_subset, reduction = "umap",group.by = "Cell.types.cymy", label = T)
FeaturePlot(rs043_subset, "monocle3_pseudotime") + scale_color_viridis_c(option = "plasma")
DotPlot(rs043_subset, features = markerGenes) + RotatedAxis()



#******************************************************#******************************************************
#******************************************************PLOTS
#******************************************************#******************************************************
library(RColorBrewer)
library("ggsci")
library("scales")

# original data figs
rs043 <- readRDS(file = "results/rs043/seurat4_RS043.rds")
DimPlot(rs043,label = T)
DimPlot(rs043, reduction = "umap",group.by = "Cell.types", label = T)
DimPlot(rs043, reduction = "umap",group.by = "Cell.clusters", label = T)
ggsave(filename = paste0("results/RS043/figs/mus_RS043_CDAA_all_clusters_label.pdf"), width = 6, height = 6)
DimPlot(rs043, reduction = "umap",group.by = "Cell.clusters", label = F)
ggsave(filename = paste0("results/RS043/figs/mus_RS043_CDAA_all_clusters_noLabel.pdf"), width = 6, height = 6)
DimPlot(rs043, reduction = "umap",group.by = "Cell.types", label = T)
ggsave(filename = paste0("results/RS043/figs/mus_RS043_CDAA_all_cellTypes_label.pdf"), width = 6, height = 6)

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
# '
DotPlot(rs043, features = markerGenes, group.by = "Cell.clusters", cluster.idents = F) + RotatedAxis()+ scale_color_viridis_c() #, cols = "RdBu"
ggsave(filename = paste0("results/RS043/figs/mus_RS043_CDAA_all_clusters_dotplot.pdf"))#, width = 6, height = 6)



# CELL TYPES UMAP
rs043_subset <- readRDS(file = "results/RS043/seurat4_RS043_HSC.rds")
DimPlot(rs043_subset,label = T)
data_req <- rs043_subset
DimPlot(object = data_req, group.by = "ident", label = T)
DimPlot(rs043_subset, reduction = "umap",group.by = "Cell.clusters.HSC", label = T)
# levels(data_req@active.ident)
ptsize <- 1.5


DimPlot(object = data_req, group.by = "Cell.clusters.HSC", label = F, pt.size = ptsize)
ggsave(filename = paste0("results/RS043/figs/mus_RS043_CDAA_HSC_clusters.pdf"))
DimPlot(object = data_req, group.by = "Cell.clusters.HSC", label = T, pt.size = ptsize)
ggsave(filename = paste0("results/RS043/figs/mus_RS043_CDAA_HSC_clusters_label.pdf"))

# cy-my cluster selection
sig <- unlist(strsplit("Acta2, Col1a1, Lox, Timp1",split = ", "))
sigName <- "Myofib"
callName <- paste0(sigName,"1")
data_req <- AddModuleScore(data_req,features = list(sig), name = sigName)
p1 <- DimPlot(object = data_req, group.by = "Cell.clusters.HSC", label = T, pt.size = ptsize)
p2 <- VlnPlot(data_req, group.by = "Cell.clusters.HSC", callName)
wrap_plots(p1,p2)
ggsave(filename = paste0("results/RS043/figs/mus_RS043_CDAA_HSC_clusters_cymy_clusters.pdf"),width = 12, height = 6)


# 
genes <- c(unlist(strsplit("Lrat, Col1a1, Acta2, Hgf",split = ", ")))
# markers required
genes <- intersect(genes, rownames(data_req@assays$RNA@data))
FeaturePlot(object = data_req, features = c(genes[1]), max.cutoff = "q99", pt.size = ptsize,order = T) + scale_color_viridis_c()

for(gene in genes){
  FeaturePlot(object = data_req, features = c(gene), max.cutoff = "q99", pt.size = ptsize,order = T) + scale_color_viridis_c()
  ggsave(filename = paste0("results/RS043/figs/mus_RS043_CDAA_HSC_gene_",gene,".pdf"))
}



# pseudotime
FeaturePlot(data_req, "monocle3_pseudotime", pt.size = ptsize) + scale_color_viridis_c(option = "plasma") + labs(title = "Mouse RS043 CDAA")
ggsave(filename = paste0("results/RS043/figs/mus_RS043_CDAA_HSC_pseudotime.pdf"))

#saved monocle object
cds <- readRDS(file = "results/RS043/monocle_RS043_HSC.rds")
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           cell_size = 0.35, #0.35
           graph_label_size=1.5)
ggsave(filename = paste0("results/RS043/figs/mus_RS043_CDAA_HSC_pseudotime_monocle.pdf"))


# cy-my signature

FeaturePlot(data_req, "cySig1", max.cutoff = "q99",pt.size = ptsize,order = T) + scale_color_viridis_c(option = "viridis") +labs(title = "Mouse RS043 CDAA cy signature")
ggsave(filename = paste0("results/RS043/figs/mus_RS043_CDAA_HSC_signature_cy.pdf"))

FeaturePlot(data_req, callName, max.cutoff = "q99",pt.size = ptsize,order = T) + scale_color_viridis_c(option = "viridis") +labs(title = "Mouse RS043 CDAA my signature")
ggsave(filename = paste0("results/RS043/figs/mus_RS043_CDAA_HSC_signature_my.pdf"))

cols_cymy <- c("#FF9900","#000099")
DimPlot(rs043_subset, reduction = "umap",group.by = "Cell.types.cymy", pt.size = ptsize,cols = cols_cymy, label = F)
ggsave(filename = paste0("results/RS043/figs/mus_RS043_CDAA_HSC_signature_clust_cymy.pdf"))


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
  ggtitle(paste0("RS043 ",gene," Mann Whitney p ",scientific(tmp$p.value,digits = 3)))
ggsave(filename = paste0("results/RS043/figs/mus_RS043_CDAA_HSC_signature_clust_cymy_vln_",gene,".pdf"),width = 5)


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
  ggtitle(paste0("RS043 ",gene," Mann Whitney p ",scientific(tmp$p.value,digits = 3)))
ggsave(filename = paste0("results/RS043/figs/mus_RS043_CDAA_HSC_signature_clust_cymy_vln_",gene,".pdf"),width = 5)

# get the cy my correlation plots
library(ggpubr)
df <- data.frame(cyHSC=rs043_subset@meta.data$cySig1, myHSC=rs043_subset@meta.data$mySig1, score=rs043_subset@meta.data$mySig1-rs043_subset@meta.data$cySig1)
head(df)
ggscatter(df, x = "cyHSC", y = "myHSC",
          add = "reg.line",                                 # Add regression line
          conf.int = TRUE,                                  # Add confidence interval
          add.params = list(color = "black", fill = "lightgray"),
          color = "score")+
  stat_cor(method = "pearson", label.x = 1, label.y = 1.)+  # Add correlation coefficient
  gradient_color(c("blue", "yellow", "red"))
ggsave(filename = paste0("results/RS043/figs/mus_RS043_CDAA_HSC_signature_clust_cymy_correln.pdf"))


