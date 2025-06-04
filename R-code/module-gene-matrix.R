###### Module-gene matrix: H19
# Normalized by column

rm(list=ls())
library(Matrix)
library(pheatmap)

setwd('')
H <- read.table('/output/h19.csv', sep = ',')

### Normalization:
for (i in 1:dim(H)[2]) {
  s <- sum(H[,i])
  if (s==0){
    H[,i] <- H[,i]
  } else {
    H[,i] <- H[,i] / sum(H[,i])
  }
}
rownames(H) <- paste0('Module',1:dim(H)[1])

### Visualize the Module-gene matrix
H_reorder <- H[c(4,2,9,10,6,3,7,8,1,5),]
fig2d <- pheatmap(H_reorder,cluster_cols = T,cluster_rows = F,
                    show_colnames = F, show_rownames = T,
                    fontsize_row = 14,fontsize_col = 14,
                    )

### Select the representative gene set for each module.
H1 <- as(t(H),'matrix')
gene_geq <- matrix(ncol = dim(H1)[2],nrow = 200)
for (m in 1:dim(H1)[2]) {
  lm <- length(which(H1[,m]>=0.95))
  if (lm < 20) {
    ef <- extractFeatures(H1,20L)
    gene_geq[1:20,m] <- rownames(H1)[ef[[m]]]
  } else {
    gene_geq[1:lm,m] <- rownames(H1)[which(H1[,m]>=0.95)]
  }
}
write.table(gene_geq,'/output/genesets_module.csv',
            sep = ',',row.names = F,col.names = F)

### Figure 2 (f) plot
go_type <- read.table('go_10type.csv',sep = ',',header = T,row.names = 1)
count <- colSums(!is.na(go_type))
new <- go_type
for (i in 1:10) {new[i,] <- go_type[i,]/count}
fig2f <- pheatmap(new,cluster_cols = F,cluster_rows = F,
                   show_colnames = T, show_rownames = T,
                   fontsize_row = 14,fontsize_col = 14,
                   display_numbers=T,fontsize_number = 14,
                   angle_col = '90')

###### Other Notes：
# Gene Ontology Enrichment Analysis was conducted on the official website of DAVID
# (https://davidbioinformatics.nih.gov/tools.jsp).