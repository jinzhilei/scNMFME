###### Data preprocessing
# Normalization, selection of high-variable genes, partition and transposing.

library(Seurat)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(dplyr)
library(reshape)
library(reshape2)

setwd('')

### Read data and match feature names
all_data <- readRDS('OneK1K_data_Seurat.rds')
all_matrix <- all_data@assays[["RNA"]]@data
gene_ensg <- rownames(all_matrix)
gene_symbol_factor <- gsub('_','-',all_data@assays[['RNA']]@meta.features[['feature_name']])
gene_symbol_char <- as.character(gene_symbol_factor)
feature_names <- factor(gene_symbol_char, labels = gene_ensg)
rownames(all_matrix) <- gene_symbol_factor
data_symbol <- CreateSeuratObject(counts = all_matrix)
data_symbol@assays[['RNA']]@meta.features[['feature_name']] <- feature_names
data_symbol[["percent.mt"]] <- PercentageFeatureSet(data_symbol,pattern = "^MT-")
data_symbol[['age']] <- all_data[['age']]
data_symbol[["predicted.celltype.l2"]] <- all_data[["predicted.celltype.l2"]]
data_symbol[["predicted.celltype.l2.score"]] <- all_data[["predicted.celltype.l2.score"]]
data_symbol[["sex"]] <- all_data[["sex"]]
figS1a <- VlnPlot(all_data,features = c("nFeature_RNA","nCount_RNA","percent.mt",ncol=3))

### Normalize
data_symbol <- NormalizeData(data_symbol)

### Select the high-variable genes
data_symbol <- FindVariableFeatures(data_symbol, selection.method = "vst", nfeatures = 3000)
top10 <- head(VariableFeatures(data_symbol), 10)
plot1 <- VariableFeaturePlot(data_symbol,pt.size = 2,cols = c("black","red"))
fig2b <- LabelPoints(plot = plot1, points = top10,repel = FALSE) +
  theme_bw()+ 
  theme(axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.key.size = unit(25,'pt'))

### Save high-variable genes
rn_ensg <- rownames(data_symbol)
rn_symbol <- data_symbol@assays[["RNA"]]@meta.features[["feature_name"]]
df <- data.frame(rn_symbol,rn_ensg)
variable_genes <- df[which(df$rn_ensg %in% VariableFeatures(data_symbol)),]
write.csv(variable_genes, "variable3000.csv", sep = ',')

### Partition
age <- as.data.frame(table(data_symbol@meta.data[['age']]))
for (i in age$Var1) {
  file_path <- paste0('/data_1997/age_',i,'_seurat.rds')
  if (file.exists(file_path) == T) {
    next
  } else {
    age_data <- subset(data_symbol, age == i)
    saveRDS(age_data, file_path)
    norm <- age_data@assays[['RNA']]@data
    rownames(norm) <- age_data@assays[["RNA"]]@meta.features[["feature_name"]]
    norm <- as(norm,'matrix')
    hv <- norm[variable_genes$rn_symbol,]
    spath <- paste0('/hv',i,'.csv')
    write.csv(hv, spath, row.names = T)
    print(paste0('file age ', i, 'has been saved'))
  }
}

###### Other Notes：
# Transpose is performed in MATLAB.