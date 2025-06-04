###### Select age-related features and conduct signature scores.
# PCA, decision tree and signature scores

library(rpart)
library(tibble)
library(bitops)
library(rattle)
library(rpart.plot)
library(RColorBrewer)
library(ggplot2)
library(tidyverse)
library(clusterProfiler)
library(stringr)
library(caret)
library(dplyr)
library(tidymodels)

setwd()

### Replace cell type annotation
all_data <- readRDS('OneK1K_data_Seurat.rds')
ct <- data.frame(age = all_data[['age']], celltype = all_data[["predicted.celltype.l2"]])
ct$Age <- as.numeric(substr(ct$age,1,2))
ct1 <- ct %>%
  mutate(celltype = recode(celltype, "erythrocyte" = "Others", "platelet" = "Others", 
                           "dendritic cell" = "DC", "natural killer cell" = "NK",
                           "CD4-positive, alpha-beta T cell" = "CD4+ T",
                           "CD8-positive, alpha-beta T cell" = "CD8+ T",
                           "plasmacytoid dendritic cell" = "DC",
                           "memory B cell" = "B", "naive B cell" = "B",
                           "gamma-delta T cell" = "Other T", "regulatory T cell" = "Other T",
                           "transitional stage B cell" = "B",
                           "naive thymus-derived CD4-positive, alpha-beta T cell" = "CD4+ T",
                           "naive thymus-derived CD8-positive, alpha-beta T cell" = "CD8+ T",
                           "central memory CD4-positive, alpha-beta T cell" = "CD4+ T",
                           "effector memory CD4-positive, alpha-beta T cell" = "CD4+ T",
                           "central memory CD8-positive, alpha-beta T cell" = "CD8+ T",
                           "effector memory CD8-positive, alpha-beta T cell" = "CD8+ T",
                           "CD4-positive, alpha-beta cytotoxic T cell" = "CD4+ T",
                           "CD16-negative, CD56-bright natural killer cell, human" = "NK",
                           "mucosal invariant T cell" = "Other T", "plasmablast" = "B",
                           "conventional dendritic cell" = "DC", "CD14-positive monocyte" = "Monocyte",
                           "innate lymphoid cell" = "Others", 
                           "CD14-low, CD16-positive monocyte" = "Monocyte",
                           "double negative thymocyte" = "Others",
                           "hematopoietic precursor cell" = "Others",
                           "peripheral blood mononuclear cell" = "Others"))
ctls <- split(ct1,ct1$age) # split by age

### Calculate the average expression
avels <- list()
agevec <- c(19:91,93:97)
for (a in 1:78) {
  age <- agevec[a]
  W <- read.table(paste0('/output/age19-97/w',age,'.csv'), sep = ',')
  colnames(W) <- paste0('module',1:dim(W)[2])
  W$cta <- ctls[[a]][["celltype"]]
  W$cta <- factor(W$cta,levels = c('Others', 'DC', 'NK', 'CD4+ T', 'CD8+ T', 
                                   'Other T', 'B', 'Monocyte'))
  avels[[a]] <- aggregate(.~ cta, data = W, FUN = mean)
}
new_row <- c('DC',rep(0,10))
avels[[78]] <- rbind(avels[[78]][1, ], new_row, avels[[78]][2:7, ])
all_fea_mat <- matrix(ncol=80,nrow=78)
for (i in 1:78) {
  mat <- as.matrix(avels[[i]][,2:11])
  all_fea_mat[i,] <- as.numeric(c(mat))
}
rownames(all_fea_mat) <- paste0('age',c(19:91,93:97))
type <- c('Others', 'DC', 'NK', 'CD4+ T', 'CD8+ T', 'Other T', 'B', 'Monocyte')
colnames(all_fea_mat) <- paste0('m',rep(1:10,each=8),'_',rep(type,10))
saveRDS(all_fea_mat,'featuremat.rds')
feadf <- data.frame(all_fea_mat)
feadf$group <- paste0('g',rep(1:4,c(22,20,20,16)))
feadf$group <- factor(feadf$group,levels = paste0('g',1:4))
coln <- colnames(feadf)
coln <- gsub("CD4\\..", "CD4+", coln)
coln <- gsub("CD8\\..", "CD8+", coln)
coln <- gsub("Other\\.", "Other", coln)
colnames(feadf) <- coln

### PCA
pca_result <- prcomp(feadf[, -81], scale. = T)
pca_scores <- data.frame(pca_result$x)
pca_scores$group <- feadf$group

# fig1-PCA plot
rot12 <- pca_result[["rotation"]][,1:2]
rot <- data.frame(id=rownames(rot12), rot12)
rotlong <- pivot_longer(rot, cols = -id, names_to = "Variable", values_to = "Loadings")
fig1pca <- ggplot(rotlong, aes(x = id, y = Loadings, fill = as.factor(Variable))) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c('#363636', '#C18F4E'))+
  labs(x = "Feature", y = "Loadings", fill = "PC") +
  theme_minimal() +
  theme(axis.title = element_text(size = 14,color = 'black',face = 'bold'),
        axis.text.x = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line.y = element_line(color = 'black',size = 0.5),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.background = element_rect(color = "black", linetype = "solid", size = 0.4),
        legend.position = c(0.1, 0.86))

# Figure 4 (d) plot
rot12 <- pca_result[["rotation"]][,1:2]
rot12p <- as.data.frame(abs(rot12))
rot12p$fea <- rownames(rot12)
rot12p$sign1 <- ifelse(rot12[,'PC1'] >= 0, "positive", "negative")
rot12p$sign2 <- ifelse(rot12[,'PC2'] >= 0, "positive", "negative")
rot12p$PC2 <- -rot12p$PC2
rot1p <- rot12p[,c('PC1','fea','sign1')]
colnames(rot1p) <- c('PC','Feature','Sign')
rot2p <- rot12p[,c('PC2','fea','sign2')]
colnames(rot2p) <- c('PC','Feature','Sign')
rot1p$Group <- rep('PC1',80)
rot2p$Group <- rep('PC2',80)
rotfin <- rbind(rot1p, rot2p)
feat <- c('m1_CD4+T','m5_CD4+T','m8_CD4+T','m10_Monocyte','m9_CD8+T','m3_CD8+T',
          'm10_NK','m3_B','m10_CD8+T')
rotfin2 <- rotfin[rotfin$Feature %in% feat,]
rotpc1_order <- rotfin2[which(rotfin2[,'Group']=='PC1'),]
rotpc1_order <- rotpc1_order[order(-rotpc1_order$PC),]
feaorder <- c(rotpc1_order$Feature)
rotfin2$Feature <- factor(rotfin2$Feature,levels = rev(feaorder))
fig4d <- ggplot(rotfin2, aes(x=PC, y=Feature))+
  geom_col(aes(fill=Sign), width = 0.05)+
  geom_point(aes(size = 0.8, color=Sign))+
  scale_x_continuous(breaks = seq(-0.25,0.25, by=0.1), labels = abs(seq(-0.25,0.25,by=0.1)))+
  geom_vline(xintercept = 0, linetype = "dashed")+
  annotate("text", x = -0.25, y = 'm1_CD4+T', label = "PC2", hjust = 0, size=6) + 
  annotate("text", x = 0.25, y = 'm1_CD4+T', label = "PC1", hjust = 0.5, size=6) +
  guides(size='none')+
  labs(fill='Sign',x='Absolute value of rotation of features')+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(color='black',size = 16,hjust = 0.5),
        axis.line = element_line(color = 'black',size = 0.6),
        axis.text = element_text(size = 12,color = 'black'),
        axis.title = element_text(size = 14,color = 'black',face = 'bold'),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.background = element_rect(color = "black", linetype = "solid", size = 0.4),
        legend.position = c(0.9, 0.25))


# Figure 4 (a) plot
center <- pca_scores %>% group_by(group) %>% summarise(xcenter = mean(PC1), ycenter = mean(PC2))
p <- ggplot(pca_scores, aes(x = PC1, y = PC2,shape=group)) +
  geom_point(size=2, color='darkgrey') +
  scale_shape_manual(values = c(15,16,17,3)) +
  stat_ellipse(aes(color = group,fill='none'), geom = "polygon", level = 0.95, alpha = 0) +
  guides(fill='none')+
  theme_minimal()+
  theme(
    panel.background = element_rect(fill = 'white'),
    plot.background = element_rect(fill = 'white', color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
fig4a <- p + geom_point(data = center, aes(x=xcenter, y=ycenter, color = group), size=3.5)+
  geom_line(data = center, aes(x=xcenter, y=ycenter), color = 'blue', linewidth=0.7, group=1)

### Selection of PC number
ind <- sample(2, nrow(feadf), replace = TRUE, prob = c(0.7,0.3))
train <- pca_scores[ind ==1,]
test <- pca_scores[ind ==2, ]
model1 <- rpart(group~.,data = train)
model2 <- rpart(group~PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data = train)
model3 <- rpart(group~PC1+PC2+PC3,data = train)
model4 <- rpart(group~PC1+PC2,data = train)
x <- subset(test,select=-group)
pred1 <- predict(model1,x,type = 'class')
k <- test[,"group"]
ta1 <- table(pred1,k)
accu1 <- sum(diag(ta1))/sum(ta1)
pred2 <- predict(model2,x,type = 'class')
k <- test[,"group"]
ta2 <- table(pred2,k)
accu2 <- sum(diag(ta2))/sum(ta2)
pred3 <- predict(model3,x,type = 'class')
k <- test[,"group"]
ta3 <- table(pred3,k)
accu3 <- sum(diag(ta3))/sum(ta3)
pred4 <- predict(model4,x,type = 'class')
k <- test[,"group"]
ta4 <- table(pred4,k)
accu4 <- sum(diag(ta4))/sum(ta4)
accu <- c(accu1, accu2, accu3, accu4)
names(accu) <- paste0('model',1:4)
best_model <- max(accu)

### Figure 4 (b)
set.seed(2345)
ct <- rpart.control(xval=5)
tree_model <- rpart(formula = group ~ PC1+PC2, data = pca_scores, method = "class", control = ct)
plotcp(tree_model)
tree_model$cptable
best_cp <- tree_model$cptable[which.min(tree_model$cptable[,"xerror"]), "CP"]
ct <- rpart.control(xval = 5, cp=best_cp)
tree_model <- rpart(formula = group ~ PC1+PC2, data = pca_scores, method = "class", control = ct)
fig4b <- prp(tree_model, main = '', yesno=1, split.yshift = 0.8,
    type = 2, extra = 2, 
    box.palette = "GnBu",  
    faclen = 0,branch = 1,tweak = 1.1, 
    clip.left.labs=F,clip.right.labs=F,border.col=0,leaf.round = 0,
    yshift=1,fallen.leaves = T) 
# calculate the accuracy
x5 <- subset(pca_scores,select=-group)
pred5 <- predict(tree_model,x5,type = 'class')
k <- pca_scores[,"group"]
ta5 <- table(pred5,k)
accu5 <- sum(diag(ta5))/sum(ta5)

### Figure 4 (c) plot
pred <- predict(tree_model,pca_scores)
ind <- apply(pred, 1, function(x) which.max(x))
ind1 <- factor(ind,levels = 1:4)
orind <- factor(rep(1:4,c(22,20,20,16)),levels = 1:4)
pred <- predict(tree_model,pca_scores)
ind <- apply(pred, 1, function(x) which.max(x))
ind1 <- factor(ind,levels = 1:4)
rocdf <- data.frame(obs=orind, pred=ind1, pred)
plot_data <- roc_curve(rocdf, obs, g1:g4)
a1 = c(1,0,0,0); a2=c(0,1,0,0); a3=c(0,0,1,0); a4=c(0,0,0,1)
true_lab <- c(rep(a1,22), rep(a2,20), rep(a3,20), rep(a4,16))
true_lab <- factor(true_lab, levels = c(1,0))
mivec <- as.vector(t(as.matrix(rocdf[,3:6])))
midf <- data.frame(label = true_lab, prop = mivec)
midata <- roc_curve(midf, label, prop)
milevel <- rep('micro',nrow(midata))
midata1 <- cbind(milevel,midata)
colnames(midata1)[1] <- '.level'
thre <- as.vector(midata$.threshold)
t_bi_label <- matrix(nrow = 78, ncol = 4)
t_bi_label[,1] <- rep(c(1,0),c(22, 56))
t_bi_label[,2] <- rep(c(0,1,0),c(22, 20, 36))
t_bi_label[,3] <- rep(c(0,1,0),c(42, 20, 16))
t_bi_label[,4] <- rep(c(0,1),c(62, 16))
senmat <- matrix(nrow = length(thre), ncol = 4)
spemat <- matrix(nrow = length(thre), ncol = 4)
for (j in 1:length(thre)) {
  for (i in 1:4) {
    predl <- factor(ifelse(rocdf[,(i+2)] >= thre[j], 1, 0), levels = c(1,0))
    orind <- factor(t_bi_label[,i],levels = c(1,0))
    conmat <- confusionMatrix(data=predl, reference=orind)
    tab <- conmat[["table"]]
    senmat[j,i] <- tab[1,1]/(sum(tab[,1]))
    spemat[j,i] <- tab[2,2]/(sum(tab[,2]))
  }
}
madf <- data.frame(.level=rep('macro',length(thre)), .threshold=midata$.threshold,
                   specificity=rowMeans(spemat), sensitivity=rowMeans(senmat))
sumdata <- rbind(plot_data, midata1, madf)
auc <- c(length=6)
for (i in 1:4) {
  fpr <- rev(c(t(1-sumdata[which(sumdata$.level==as.character(i)),'specificity'])))
  tpr <- rev(c(t(sumdata[which(sumdata$.level==as.character(i)),'sensitivity'])))
  auc[i] <- sum(diff(fpr) * (tpr[-length(tpr)] + tpr[-1]) / 2)
}
fprma <- rev(c(t(1-sumdata[which(sumdata$.level=='macro'),'specificity'])))
tprma <- rev(c(t(sumdata[which(sumdata$.level=='macro'),'sensitivity'])))
auc[5] <- sum(diff(fprma) * (tprma[-length(tprma)] + tprma[-1]) / 2)
fprmi <- rev(c(t(1-sumdata[which(sumdata$.level=='micro'),'specificity'])))
tprmi <- rev(c(t(sumdata[which(sumdata$.level=='micro'),'sensitivity'])))
auc[6] <- sum(diff(fprmi) * (tprmi[-length(tprmi)] + tprmi[-1]) / 2)
fig4c <- ggplot(sumdata, aes(1-specificity, sensitivity, color=.level, linetype = .level)) +
  geom_line(linewidth = 1.3) +
  geom_abline(linetype = 2) +
  scale_linetype_manual(values = rep(c('solid','dashed'), c(4,2)),
                        guide='none')+
  scale_color_manual(values = c("#F8766D", "#7CAE00", "#00BFC4", "#DB72FB", "#619CFF", "#D39200"),
                     labels = c('g1 (AUC=0.91)', 'g2 (AUC=0.75)',
                                'g3 (AUC=0.98)', 'g4 (AUC=0.95)',
                                'micro (AUC=0.93)','macro (AUC=0.93)'))+
  guides(color = guide_legend(title = NULL)) + 
  labs(x='False Positive Rate', y='True Positive Rate')+
  theme_bw() +
  theme(legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12),
        legend.position = c(0.75, 0.25),
        axis.text = element_text(size = 10,color = 'black'),
        axis.title = element_text(size=12, face = 'bold'),
        axis.ticks.length = unit(1.5,'mm'))

### Figure 4 (e) plot
score1 <- (-0.1661) * feadf$`m1_CD4+T` - 0.1656 * feadf$`m5_CD4+T` - 0.1655 * feadf$`m8_CD4+T`
score2 <- (-0.2289) * feadf$m10_Monocyte -0.2215 * feadf$m3_B -0.2093 * feadf$m10_NK -0.2005 * feadf$`m10_CD8+T` + 0.2204 * feadf$`m9_CD8+T` + 0.2078 * feadf$`m3_CD8+T`
sdf <- data.frame(score1 = score1, score2 = score2)
sdf$group <- paste0('g',rep(1:4,c(22,20,20,16)))
sdf$group <- factor(sdf$group,levels = paste0('g',1:4))
sdf$age <- c(19:91,93:97)
y <- sdf$score1
x <- c(19:91,93:97)
data <- data.frame(x = x, y = y)
data$group <- sdf$group
model_cubic <- lm(y ~ poly(x, 4, raw = TRUE), data = data)
m <- model_cubic
eq <- substitute(italic(y) == a + b~italic(x) ~ c~italic(x)^2 + d~italic(x)^3 ~ e~italic(x)^4*","~~italic(r)^2~"="~r2, 
                 list(a = format(unname(coef(m)[1]), digits = 2),
                      b = format(unname(coef(m)[2]), digits = 2),
                      c = format(unname(coef(m)[3]), digits = 2),
                      d = format(unname(coef(m)[4]), digits = 2),
                      e = format(unname(coef(m)[5]), digits = 2),
                      r2 = format(summary(m)$r.squared, digits = 3)))
labelm1 <- as.character(as.expression(eq))
fig4e1 <- ggplot(data, aes(x = x, y = y)) +
  geom_point(size=2,aes(color=group)) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 4, raw = TRUE), color = "blue")+
  geom_text(x = 50, y = -0.07, label = labelm1, parse = TRUE, size = 4.5)+
  labs(x = "Age", y = 'Score1', title = 'Signature score1')+
  theme_bw()+
  theme(plot.title = element_text(size=14, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 10,color = 'black'),
        axis.title = element_text(size=12),
        axis.ticks.length = unit(1.5,'mm'),
        legend.position = 'none',
        strip.text.x = element_text(size = 12),
        strip.background = element_blank(), strip.placement = 'inside',
        panel.spacing = unit(0.01,'lines'),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())
model_linear <- lm(sdf$score2 ~ c(19:91,93:97), data = sdf)
m <- model_linear
eq <- substitute(italic(y) == a ~ + ~ b~italic(x)*","~~italic(r)^2~"="~r2, 
                 list(a = format(unname(coef(m)[1]), digits = 2),
                      b = format(unname(coef(m)[2]), digits = 2),
                      r2 = format(summary(m)$r.squared, digits = 3)))
labelm2 <- as.character(as.expression(eq))
fig4e2 <- ggplot(sdf, aes(x = age,y = score2)) + 
  geom_point(alpha=1, size = 3, shape = 16, aes(color = group)) +
  geom_smooth(method = "lm", formula = y~x, color = "blue")+
  geom_text(x = 40, y = -0.06, label = labelm2, parse = TRUE, size = 4.5)+
  labs(x = "Age", y = 'Score2', title = 'Signature score2')+
  theme_bw()+
  theme(plot.title = element_text(size=14, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 10,color = 'black'),
        axis.title = element_text(size=12),
        axis.ticks.length = unit(1.5,'mm'),
        legend.position = 'none',
        strip.text.x = element_text(size = 12),
        strip.background = element_blank(), strip.placement = 'inside',
        panel.spacing = unit(0.01,'lines'),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())

 