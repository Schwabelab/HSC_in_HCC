# 2022 Feb 10 created
# this is the Seurat4 analysis of RS024

setwd("D:/pCloud Sync/Columbia/SchwabeLab/Aveline_HSCnHCC/github/HSCinHCC/")

library(Seurat)
library(ggplot2)
library(monocle3)
library(SeuratWrappers)
library(patchwork)
library(magrittr)
markerGenes <- c("Lrat","Myh11","Actg2","Rgs5","Acta2","Lox","Col3a1","Col1a1","Mgp","Eln","Mfap4","Dpt","Gsn","Kdr","Aqp1","Top2a","Mki67",'Upk3b','Msln', 'Epcam', 'Krt7',"Il7r","Cd3d","Trac","Nkg7","Clec4f","Vsig4","Lgals3","Ccr2","Bpgm","Cd79a","Ebf1","Jchain","Siglech","Alb","Serpina1a")


data <- Read10X_h5("data/RS024_Mdr2KO/RS024_filtered_feature_bc_matrix.h5")
dim(data)
rs024 <- CreateSeuratObject(counts = data, project = "RS024 Mdr2KO", min.cells = 3, min.features = 200)
rs024
rs024[["percent.mt"]] <- PercentageFeatureSet(rs024, pattern = "^mt-")
rs024[["percent.ribo"]] <- PercentageFeatureSet(rs024, pattern = "^Rpl|^Rps|Mrps|^Mrpl")#

# Visualize QC metrics as a violin plot
VlnPlot(rs024, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p1 <- FeatureScatter(rs024, feature1 = "nCount_RNA", feature2 = "percent.mt") + NoLegend()
p2 <- FeatureScatter(rs024, feature1 = "nFeature_RNA", feature2 = "percent.mt") + NoLegend()
p3 <- FeatureScatter(rs024, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# # keep only the required cells
wrap_plots(list(p1,p2,p3),nrow = 2,ncol = 2)
# keep only the required cells



rs024 <- NormalizeData(rs024, normalization.method = "LogNormalize", scale.factor = 10000)
rs024 <- FindVariableFeatures(rs024, selection.method = "vst", nfeatures = 2000)
# all.genes <- rownames(rs024)

DimPlot(rs024,label = T)
DotPlot(rs024, features = markerGenes, cluster.idents = T) + RotatedAxis()
DotPlot(rs024, features = markerGenes) + RotatedAxis()


# getting marker genes
marker <- FindAllMarkers(rs024,only.pos = T)
# write.csv(marker, file = "results/RS024/markers.csv")

FeaturePlot(object = rs024, features = c("percent.mt"))
FeaturePlot(object = rs024, features = c("percent.ribo"))
FeaturePlot(object = rs024, features = c("nCount_RNA"))
FeaturePlot(object = rs024, features = c("nFeature_RNA"))
VlnPlot(object = rs024, features = c("nCount_RNA"))
VlnPlot(object = rs024, features = c("nFeature_RNA"))

# giving cell types
DimPlot(rs024,label = T)
# rs024[["Cell.clusters"]] <- Idents(object = rs024)
rs024 <- SetIdent(rs024, value = "Cell.clusters")
rs024 <- RenameIdents(object = rs024, 
                              '0' = "HSC",
                              '1' = "HSC",
                              '4' = "HSC",
                              '5' = "HSC",
                              '6' = "HSC",
                              '9' = "PF",
                              '12' = "VSMC",
                              '2' = "Cholangiocyte",
                              '3' = "Cholangiocyte",
                              '11' = "Cycling",
                              '10' = "PF/Mesothelial",
                              '13' = "Hepatocytes",
                              '8' = "Endothelial",
                              '7' = "Immune")

rs024[["Cell.types"]] <- Idents(object = rs024)
DimPlot(rs024,label = T, group.by = "Cell.types")
table(rs024@active.ident)

DimPlot(rs024, reduction = "umap",group.by = "Cell.types", label = F)
DimPlot(rs024, reduction = "umap",group.by = "Cell.clusters", label = F)

# saveRDS(rs024, file = "results/RS024/seurat4_RS024.rds")


#****************************************************************saved object for HSC-noPF analysis
rs024 <- readRDS(file = "results/RS024/seurat4_RS024.rds")
DimPlot(rs024,label = T)
DimPlot(rs024, reduction = "umap",group.by = "Cell.types", label = T)
DotPlot(rs024, features = markerGenes) + RotatedAxis()


# Recluster the HSC population 
rs024 <- SetIdent(rs024, value = "Cell.clusters")
rs024_subset <- subset(rs024, idents = c(0,1,4:6))

DimPlot(rs024_subset, reduction = "umap",label = T)
DotPlot(rs024_subset, features = markerGenes) + RotatedAxis()
DimPlot(rs024_subset, reduction = "umap",group.by = "Cell.clusters", label = T)
FeaturePlot(object = rs024_subset, features = c("Col1a1")) & scale_color_viridis_c()
# removing mixed cells

DimPlot(rs024_subset, reduction = "umap",label = T)
DotPlot(rs024_subset, features = markerGenes) + RotatedAxis()
DimPlot(rs024_subset, reduction = "umap",group.by = "Cell.clusters", label = T)
FeaturePlot(object = rs024_subset, features = c("Col1a1")) & scale_color_viridis_c()
VlnPlot(rs024_subset, features = c("Lrat","Hgf","Col1a1"))

# saveRDS(rs024_subset, file = "results/RS024/seurat4_RS024_HSC_noPF.rds")
#****************************************************************saved object
rs024_subset <- readRDS(file = "results/RS024/seurat4_RS024_HSC_noPF.rds")
DimPlot(rs024_subset,label = T)


# get the HGF and COL1a1 GE in clusters
VlnPlot(rs024_subset, features = c("Lrat","Hgf","Col1a1"))
markers <- FindAllMarkers(rs024_subset,features =c("Lrat","Hgf","Col1a1"), min.pct = 0.001,logfc.threshold = 0.0001,return.thresh=1)
markers





# ************************************************** trajectory analysis
# getting to Monocle


#************************************************************* identifying most activated and quiescent clusters
rs024_subset <- readRDS(file = "results/RS024/seurat4_RS024_HSC_noPF.rds")
DimPlot(rs024_subset,label = T)
DimPlot(rs024_subset, reduction = "umap",group.by = "Cell.types", label = T)
DimPlot(rs024_subset, reduction = "umap",group.by = "Cell.clusters", label = T)
FeaturePlot(rs024_subset, "monocle3_pseudotime") + scale_color_viridis_c(option = "plasma")
DotPlot(rs024_subset, features = markerGenes) + RotatedAxis()


# myoFib signature
sig <- unlist(strsplit("Acta2, Col1a1, Lox, Timp1",split = ", "))

rs024_subset <- AddModuleScore(rs024_subset,features = list(sig), name = "Myofib")
FeaturePlot(rs024_subset, "Myofib1") + scale_color_viridis_c(option = "viridis")
VlnPlot(rs024_subset, "Myofib1", sort = T, group.by = "Cell.clusters.HSC")

clust_my <- 2
clust_cy <- 1
# compute avg scre for the clusters
tmppos <- which(rs024_subset@meta.data$Cell.clusters.HSC==clust_my)
score_my <- mean(rs024_subset@meta.data$Myofib1[tmppos])
tmppos <- which(rs024_subset@meta.data$Cell.clusters.HSC==clust_cy)
score_cy <- mean(rs024_subset@meta.data$Myofib1[tmppos])
avg_myofibScore_rs024 <- data.frame(score_my, score_cy)
colnames(avg_myofibScore_rs024) <- c(paste0("Clust_",clust_my),paste0("Clust_",clust_cy))
avg_myofibScore_rs024

# get the cy-my DEGs
DimPlot(rs024_subset, reduction = "umap", label = T)
DEGs <- FindMarkers(rs024_subset, ident.1 = c(clust_my), ident.2 = clust_cy)
DEGs <- DEGs[order(DEGs$avg_log2FC, decreasing = T),]
head(DEGs,20)
tail(DEGs,20)
# write.csv(DEGs, file = "results/RS024/DEGs_cymy_noPF.csv")

# getting marker genes
marker <- FindAllMarkers(rs024_subset,only.pos = T)
marker <- marker[order(marker$cluster, marker$avg_log2FC, decreasing = T),]
head(marker,20)
tail(marker,20)
# write.csv(marker, file = "results/RS024/markers_HSC_noPF.csv")



# ************************************************** cy-my signature analysis
library(readxl)
dataset <- data.frame(read_excel("results/cymySig_secreated.xlsx", skip = 1))
head(dataset)

cySig <- na.omit(dataset$cySig2)
mySig <- na.omit(dataset$mySgi2)

ptsize <- 2
sig <- cySig
sigName <- "cySig"
callName <- paste0(sigName,"1")
rs024_subset <- AddModuleScore(rs024_subset,features = list(sig), name = sigName)
FeaturePlot(rs024_subset, callName,pt.size = ptsize,order = T) + scale_color_viridis_c(option = "viridis") +labs(title = "RS024 Mdr2ko cy signature")


sig <- mySig
sigName <- "mySig"
callName <- paste0(sigName,"1")
rs024_subset <- AddModuleScore(rs024_subset,features = list(sig), name = sigName)
FeaturePlot(rs024_subset, callName,pt.size = ptsize,order = T) + scale_color_viridis_c(option = "viridis")+labs(title = "RS024 Mdr2ko my signature")

VlnPlot(rs024_subset, features = c("mySig1", "cySig1"), sort = F)

# quantify the cy-my score average in different clusters

# rs024_subset[["Cell.clusters.HSC"]] <- Idents(object = rs024_subset)
# DimPlot(rs024_subset, reduction = "umap",group.by = "Cell.clusters.HSC", label = T)
rs024_subset <- SetIdent(rs024_subset, value = "Cell.clusters.HSC")
rs024_subset <- RenameIdents(object = rs024_subset, 
                      '0' = "cyHSC",
                      '1' = "cyHSC",
                      '3' = "cyHSC",
                      '5' = "cyHSC",
                      '6' = "cyHSC",
                      
                      '2' = "myHSC",
                      '4' = "myHSC")


rs024_subset[["Cell.types.cymy"]] <- Idents(object = rs024_subset)
DimPlot(rs024_subset,label = T, group.by = "Cell.types.cymy")
table(rs024_subset@active.ident)

cols_cymy <- c("#FF9900","#000099")
DimPlot(rs024_subset, reduction = "umap",group.by = "Cell.types.cymy", pt.size = 2,cols = cols_cymy, label = F)

# saveRDS(rs024_subset, file = "results/RS024/seurat4_RS024_HSC_noPF.rds")



#****************************************************************saved object
rs024_subset <- readRDS(file = "results/RS024/seurat4_RS024_HSC_noPF.rds")
DimPlot(rs024_subset,label = T)
DimPlot(rs024_subset, reduction = "umap",group.by = "Cell.types", label = T)
DimPlot(rs024_subset, reduction = "umap",group.by = "Cell.clusters", label = T)
DimPlot(rs024_subset, reduction = "umap",group.by = "Cell.clusters.HSC", label = T)
DimPlot(rs024_subset, reduction = "umap",group.by = "Cell.types.cymy", label = T)
FeaturePlot(rs024_subset, "monocle3_pseudotime") + scale_color_viridis_c(option = "plasma")
DotPlot(rs024_subset, features = markerGenes) + RotatedAxis()



#******************************************************#******************************************************
#******************************************************PLOTS
#******************************************************#******************************************************
library(RColorBrewer)
library("ggsci")
library("scales")

# original data figs
rs024 <- readRDS(file = "results/RS024/seurat4_RS024.rds")
DimPlot(rs024,label = T)
DimPlot(rs024, reduction = "umap",group.by = "Cell.types", label = T)
DimPlot(rs024, reduction = "umap",group.by = "Cell.clusters", label = T)
ggsave(filename = paste0("results/RS024/figs/mus_RS024_Mdr2ko_all_clusters_label.pdf"), width = 6, height = 6)
DimPlot(rs024, reduction = "umap",group.by = "Cell.clusters", label = F)
ggsave(filename = paste0("results/RS024/figs/mus_RS024_Mdr2ko_all_clusters_noLabel.pdf"), width = 6, height = 6)

markerGenes <- c("Lrat","Colec11","Des", #"Col1a1", #"Lox","Col3a1",
                 "Dcn","Pdgfra","Pdgfrb",
                 "Myh11","Actg2", #"Rgs5","Acta2",
                 "Dpt","Gpx3", #,"Gsn", #"Mgp","Eln","Mfap4"
                 'Upk3b','Msln',
                 "Kdr","Aqp1",
                 # "Il7r","Cd3d","Nkg7", #,"Trac"
                 "Ctss","Ptprc", #"Trem2",
                 # "Clec4f","Vsig4",
                 # "Lgals3","Ccr2","Bpgm",
                 # "Cd79a","Ebf1","Jchain","Siglech",
                 "Top2a","Mki67",
                 'Epcam', 'Krt7',
                 "Alb","Serpina1a")
DotPlot(rs024, features = markerGenes, group.by = "Cell.clusters", cluster.idents = F) + RotatedAxis()+ scale_color_viridis_c() #, cols = "RdBu"
ggsave(filename = paste0("results/RS024/figs/mus_RS024_Mdr2ko_all_clusters_dotplot.pdf"))#, width = 6, height = 6)

DimPlot(rs024, reduction = "umap",group.by = "Cell.types", label = T)
ggsave(filename = paste0("results/RS024/figs/mus_RS024_Mdr2ko_all_cellTypes_label.pdf"), width = 6, height = 6)

# CELL TYPES UMAP
data_req <- rs024_subset
DimPlot(object = data_req, group.by = "ident", label = T)
DimPlot(rs024_subset, reduction = "umap",group.by = "Cell.clusters.HSC", label = T)
# levels(data_req@active.ident)
ptsize <- 2


DimPlot(object = data_req, group.by = "Cell.clusters.HSC", label = F, pt.size = ptsize)
ggsave(filename = paste0("results/RS024/figs/mus_RS024_Mdr2ko_HSC_noPF_clusters.pdf"))
 

# 
genes <- c(unlist(strsplit("Lrat, Col1a1, Acta2, Hgf",split = ", ")))
# markers required
genes <- intersect(genes, rownames(data_req@assays$RNA@data))
FeaturePlot(object = data_req, features = c(genes[1]), max.cutoff = "q99", pt.size = ptsize,order = T) + scale_color_viridis_c()

for(gene in genes){
  FeaturePlot(object = data_req, features = c(gene), max.cutoff = "q99", pt.size = ptsize,order = T) + scale_color_viridis_c()
  ggsave(filename = paste0("results/RS024/figs/mus_RS024_Mdr2ko_HSC_noPF_gene_",gene,".pdf"))
}



# pseudotime
FeaturePlot(data_req, "monocle3_pseudotime", pt.size = ptsize) + scale_color_viridis_c(option = "plasma") + labs(title = "Mouse RS024 Mdr2ko")
ggsave(filename = paste0("results/RS024/figs/mus_RS024_Mdr2ko_HSC_noPF_pseudotime.pdf"))

#saved monocle object
cds <- readRDS(file = "results/RS024/monocle_RS024_HSC_noPF.rds")
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           cell_size = 0.35, #0.35
           graph_label_size=1.5)
ggsave(filename = paste0("results/RS024/figs/mus_RS024_Mdr2ko_HSC_noPF_pseudotime_monocle.pdf"))


# cy-my signature

FeaturePlot(data_req, "cySig1", max.cutoff = "q99",pt.size = ptsize,order = T) + scale_color_viridis_c(option = "viridis") +labs(title = "Mouse RS024 Mdr2ko cy signature")
ggsave(filename = paste0("results/RS024/figs/mus_RS024_Mdr2ko_HSC_noPF_signature_cy.pdf"))

FeaturePlot(data_req, callName, max.cutoff = "q99",pt.size = ptsize,order = T) + scale_color_viridis_c(option = "viridis") +labs(title = "Mouse RS024 Mdr2ko my signature")
ggsave(filename = paste0("results/RS024/figs/mus_RS024_Mdr2ko_HSC_noPF_signature_my.pdf"))


DimPlot(rs024_subset, reduction = "umap",group.by = "Cell.types.cymy", pt.size = ptsize,cols = cols_cymy, label = F)
ggsave(filename = paste0("results/RS024/figs/mus_RS024_Mdr2ko_HSC_noPF_signature_clust_cymy.pdf"))


# violin plots
gene <- "Col1a1"
p <- VlnPlot(data_req, features = gene, group.by = "Cell.types.cymy", cols = cols_cymy) #+geom_boxplot()


gene <- "Hgf" 
p <- VlnPlot(data_req, features = gene, group.by = "Cell.types.cymy", cols = cols_cymy) #+geom_boxplot()

# get the cy my correlation plots
library(ggpubr)
df <- data.frame(cyHSC=rs024_subset@meta.data$cySig1, myHSC=rs024_subset@meta.data$mySig1, score=rs024_subset@meta.data$mySig1-rs024_subset@meta.data$cySig1)
head(df)
ggscatter(df, x = "cyHSC", y = "myHSC",
          add = "reg.line",                                 # Add regression line
          conf.int = TRUE,                                  # Add confidence interval
          add.params = list(color = "black", fill = "lightgray"),
          color = "score")+
  stat_cor(method = "pearson", label.x = 0.75, label.y = 1.2)+  # Add correlation coefficient
  gradient_color(c("blue", "yellow", "red"))
ggsave(filename = paste0("results/RS024/figs/mus_RS024_Mdr2ko_HSC_noPF_signature_clust_cymy_correln.pdf"))

