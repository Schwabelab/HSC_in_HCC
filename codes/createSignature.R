# 2022 Feb 07 created
# this is the script for creating the cy-my signature

setwd("D:/pCloud Sync/Columbia/SchwabeLab/Aveline_HSCnHCC/github/HSCinHCC/")

library(Seurat)
library(ggplot2)
library(monocle3)
library(SeuratWrappers)
library(patchwork)
library(magrittr)

# common vars
logFC_thresh <- 0.5

# function for integrating the p-values, FC, pct from multiple DEG tests for a set of common genes
df_commonDEGs <- function(commonDEGs,df1,df2,df3,df4){
  # i=commonDEGs[1]
  df_cmnDEGs <- data.frame(matrix(nrow = length(commonDEGs), ncol = 4))
  rownames(df_cmnDEGs) <- commonDEGs
  colnames(df_cmnDEGs) <- c("avg_logFC", "pct.1", "pct.2","p_val_adj")
  head(df_cmnDEGs)
  for (i in commonDEGs) {
    # get the gene from the diff datasets
    temp_pos <- which(df1$gene==i)
    temp_df <- df1[temp_pos,]
    temp_pos <- which(df2$gene==i)
    temp_df <- rbind(temp_df,df2[temp_pos,])
    temp_pos <- which(df3$gene==i)
    temp_df <- rbind(temp_df,df3[temp_pos,])
    temp_pos <- which(df4$gene==i)
    temp_df <- rbind(temp_df,df4[temp_pos,])
    
    # compute median values
    temp_final <- apply(temp_df[,c("avg_log2FC", "pct.1", "pct.2")],2,mean)
    
    # convert to z-scores, stouffer integrate, reconvert to pvalues
    # oneSided_logFC <- any(temp_df$avg_logFC<0) #consider the gene only if the log_FC are all in one direction; but here only positive logFC are considered
    # if(oneSided_logFC){temp <- qnorm(temp_df$p_val_adj,lower.tail = F)}
    temp <- qnorm(temp_df$p_val_adj,lower.tail = F)
    temp[which(temp>40)] <- 40 #saturating off the z-scores to remove Inf; pnorm(40, lower.tail = F) = 0; after z-score >37 the p-vals are zero
    temp <- sum(temp)/sqrt(length(temp)) #stouffer integrate
    temp_final["p_val_adj"] <- pnorm(temp, lower.tail = F)
    # put all the values to a data frame
    df_cmnDEGs[i,] <- temp_final
  }
  df_cmnDEGs <- df_cmnDEGs[order(df_cmnDEGs$avg_logFC, decreasing = T),]
  return(df_cmnDEGs)
}



# ****************************************************************
# creating signature from just cy-my markers
# ****************************************************************

marker024 <- read.csv("results/RS024/DEGs_cymy_noPF.csv", row.names=1)
head(marker024)
tail(marker024)
# read in the myHSC markers
df_my024 <- marker024[which(marker024$avg_log2FC>=logFC_thresh),]
df_my024$gene <- rownames(df_my024)
head(df_my024)
tail(df_my024)
# read in the cyHSC markers
df_cy024 <- marker024[which(marker024$avg_log2FC<=-logFC_thresh),]
df_cy024$avg_log2FC <- df_cy024$avg_log2FC*-1
df_cy024 <- df_cy024[order(df_cy024$avg_log2FC, decreasing = T),]
df_cy024$gene <- rownames(df_cy024)
head(df_cy024)
tail(df_cy024)


marker039 <- read.csv("results/RS039/DEGs_cymy.csv", row.names=1)
head(marker039)
tail(marker039)
# read in the myHSC markers
df_my039 <- marker039[which(marker039$avg_log2FC>=logFC_thresh),]
df_my039$gene <- rownames(df_my039)
head(df_my039)
tail(df_my039)
# read in the cyHSC markers
df_cy039 <- marker039[which(marker039$avg_log2FC<=-logFC_thresh),]
df_cy039$avg_log2FC <- df_cy039$avg_log2FC*-1
df_cy039 <- df_cy039[order(df_cy039$avg_log2FC, decreasing = T),]
df_cy039$gene <- rownames(df_cy039)
head(df_cy039)
tail(df_cy039)

marker042 <- read.csv("results/RS042/DEGs_cymy.csv", row.names=1)
head(marker042)
tail(marker042)
# read in the myHSC markers
df_my042 <- marker042[which(marker042$avg_log2FC>=logFC_thresh),]
df_my042$gene <- rownames(df_my042)
head(df_my042)
tail(df_my042)
# read in the cyHSC markers
df_cy042 <- marker042[which(marker042$avg_log2FC<=-logFC_thresh),]
df_cy042$avg_log2FC <- df_cy042$avg_log2FC*-1
df_cy042 <- df_cy042[order(df_cy042$avg_log2FC, decreasing = T),]
df_cy042$gene <- rownames(df_cy042)
head(df_cy042)
tail(df_cy042)

marker043 <- read.csv("results/RS043/DEGs_cymy.csv", row.names=1)
head(marker043)
tail(marker043)
# read in the myHSC markers
df_my043 <- marker043[which(marker043$avg_log2FC>=logFC_thresh),]
df_my043$gene <- rownames(df_my043)
head(df_my043)
tail(df_my043)
# read in the cyHSC markers
df_cy043 <- marker043[which(marker043$avg_log2FC<=-logFC_thresh),]
df_cy043$avg_log2FC <- df_cy043$avg_log2FC*-1
df_cy043 <- df_cy043[order(df_cy043$avg_log2FC, decreasing = T),]
df_cy043$gene <- rownames(df_cy043)
head(df_cy043)
tail(df_cy043)

# integrate the common DEGs
# myHSC
# method: mySig2_cmy_MDR
commonDEGs_my <- intersect(intersect(intersect(df_my024$gene,df_my039$gene),df_my042$gene),df_my043$gene)
commonDEGs_my
df_commonDEGs_my <- df_commonDEGs(commonDEGs_my,df_my024,df_my039,df_my042,df_my043)
head(df_commonDEGs_my)
tail(df_commonDEGs_my)

# write.csv(df_commonDEGs_my, file = "results/myHSC_mySig2_cmy_MDR.csv")


# cyHSC
# method: cySig2_cmy_MDR
commonDEGs_cy <- intersect(intersect(intersect(df_cy024$gene,df_cy039$gene),df_cy042$gene),df_cy043$gene)
commonDEGs_cy
df_commonDEGs_cy <- df_commonDEGs(commonDEGs_cy,df_cy024,df_cy039,df_cy042,df_cy043)
colnames(df_commonDEGs_cy) <- c("avg_logFC","pct.2", "pct.1", "p_val_adj")
head(df_commonDEGs_cy)
tail(df_commonDEGs_cy)

# write.csv(df_commonDEGs_cy, file = "results/cyHSC_cySig2_cmy_MDR.csv")


#***************************************************************
# pathway analysis of the cy-my signatures
#****************************************************************
library(enrichR)
library(readxl)
dataset <- data.frame(read_excel("results/cymySig_secreated.xlsx", skip = 1))
head(dataset)

cySig <- na.omit(dataset$cySig2)
mySig <- na.omit(dataset$mySgi2)


# pathway analysis
genelist <- as.character(cySig)
genelist
dbs <- listEnrichrDbs()
# dbs
dbs <- c("GO_Biological_Process_2015")
enriched <- enrichr(genelist, dbs)
head(enriched)
# writing to file
df <- enriched[[1]]
head(df)
# write.csv(df, file = paste0("results/pathway_GO_cymySig_cy.csv"))

genelist <- as.character(mySig)
genelist
enriched <- enrichr(genelist, dbs)
head(enriched)
# writing to file
df <- enriched[[1]]
head(df)
# write.csv(df, file = paste0("results/pathway_GO_cymySig_my.csv"))


