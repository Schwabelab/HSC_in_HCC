# 2021 April 05 created 
# this plots the required heatmaps for wLiv CCl4



library(pheatmap)
setwd("D:/pCloud Sync/Columbia/SchwabeLab/Aveline_HSCnHCC/github/HSCinHCC/cpdb/")


heatmaps_plot = function(count_matrix, show_rownames = T, show_colnames = T,scale="none", cluster_cols = T,border_color='white', cluster_rows = T, fontsize_row=11,fontsize_col = 11, 
                         main = '',treeheight_row=0, family='Arial', treeheight_col = 0,col1 = "dodgerblue4", col2 = 'peachpuff', col3 = 'deeppink4', count_filename=count_filename){
  col.heatmap <- colorRampPalette(c(col1,col2,col3 ))( 1000 )
  # pheatmap(count_matrix)
  pheatmap(count_matrix, show_rownames = show_rownames, show_colnames = show_colnames, scale=scale, cluster_cols = cluster_cols,
           border_color=border_color, cluster_rows = cluster_rows, fontsize_row = fontsize_row, fontsize_col = fontsize_col,
           main = main, treeheight_row = treeheight_row, family = family,color = col.heatmap, treeheight_col = treeheight_col, filename = count_filename)
}

heatmaps_plot_angle = function(count_matrix, show_rownames = T, show_colnames = T,scale="none", cluster_cols = T,border_color='white', cluster_rows = T, fontsize_row=11,fontsize_col = 11, 
                         main = '',treeheight_row=0, family='Arial', treeheight_col = 0,col1 = "dodgerblue4", col2 = 'peachpuff', col3 = 'deeppink4', count_filename=count_filename){
  col.heatmap <- colorRampPalette(c(col1,col2,col3 ))( 1000 )
  # pheatmap(count_matrix)
  pheatmap(count_matrix, show_rownames = show_rownames, show_colnames = show_colnames, scale=scale, cluster_cols = cluster_cols,
           border_color=border_color, cluster_rows = cluster_rows, fontsize_row = fontsize_row, fontsize_col = fontsize_col,
           main = main, treeheight_row = treeheight_row, family = family,color = col.heatmap, treeheight_col = treeheight_col, filename = count_filename, angle_col = 45,width = 10, height=10)
}

#*****************************************************************************************
#************************#creating my heatmap- hsc subtypes
#*****************************************************************************************

count_network <- read.delim("wLiv_CCl4/out_newDBAll_HSCsubTypes/count_network.txt", row.names=NULL, stringsAsFactors=FALSE)
head(count_network)
count_matrix = matrix(count_network[,"count"], nrow=length(unique(count_network[,"SOURCE"])), ncol=length(unique(count_network[,"TARGET"])), byrow = T)
rownames(count_matrix)= unique(count_network[,"SOURCE"])
colnames(count_matrix)= unique(count_network[,"TARGET"])
count_matrix

count_filename <- "wLiv_CCl4/out_newDBAll_HSCsubTypes/cpdb_wLivCCl4_heatmap_cymy_Hep.pdf"
hep_pos <- 8
pos <- order(count_matrix[hep_pos,], decreasing = T) #ordering the hep interactions
count_matrix[c(rep(hep_pos,1)),pos]

heatmaps_plot_angle(count_matrix=count_matrix[c(rep(hep_pos,14)),pos], count_filename=count_filename, cluster_rows = F, cluster_cols = F)


#*****************************************************************************************
#************************#creating my heatmap- Hep vs hsc subtypes +  others
#*****************************************************************************************
#*
count_network <- read.delim("wLiv_CCl4/out_newDBAll_HSC_HEPsubTypes/count_network.txt", row.names=NULL, stringsAsFactors=FALSE)
head(count_network)
count_matrix = matrix(count_network[,"count"], nrow=length(unique(count_network[,"SOURCE"])), ncol=length(unique(count_network[,"TARGET"])), byrow = T)
rownames(count_matrix)= unique(count_network[,"SOURCE"])
colnames(count_matrix)= unique(count_network[,"TARGET"])
count_matrix

tmppos <- grep("Hep",rownames(count_matrix))
# tmppos <- c(tmppos,grep("Hep",rownames(count_matrix)))
count_matrix[tmppos,]

count_filename <- "wLiv_CCl4/out_newDBAll_HSC_HEPsubTypes/cpdb_wLivCCl4_heatmap_HepSubtypes_vs_others.pdf"

heatmaps_plot(count_matrix=count_matrix[c(tmppos,tmppos),-tmppos], count_filename=count_filename, cluster_rows = F, cluster_cols = T)
