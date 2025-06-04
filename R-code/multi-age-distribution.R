###### Multi-age distribution analysis
# Set thresholds and plot the distribution figures

library(ggplot2)
library(patchwork)
library(tidyverse)
library(dplyr)
library(reshape2)
library(gghalves)

setwd('')
W <- read.table('/output/w19.csv', sep = ',')

### Figure 3 (a) plot
wls <- list()
c0 <- c(0.05,0.1,0.15,0.15,0.1,0.1,0.2,0.1,0.15,0.2)
for (i in 1:10) {
  dense = data.frame(density(W[,i])[c('x','y')])
  m <- data.frame(value = W[,i])
  wls[[i]] <- ggplot(m,aes(x =value))+
    geom_histogram(aes(y=after_stat(density)), color="grey", alpha=.25, fill="#fffbf0", bins=70)+
    geom_density(color=NA)+
    geom_area(data = subset(dense,x >= 0 & x < c0[i]), aes(x, y, fill = "Part i0"), alpha=.4)+
    geom_area(data = subset(dense,x >= c0[i]), aes(x, y, fill = "Part i1"), alpha=.4)+
    scale_fill_manual("Parts", 
                      breaks = c("Part i0", "Part i1"), 
                      values = c("Part i0"="#4b5cc466", "Part i1"="#16a95166"))+
    labs(title=paste0('module',i))+
    theme_bw()+
    theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
          axis.text = element_text(size = 12),
          axis.title = element_blank(),
          legend.position = 'none',
          axis.ticks.length = unit(2,'mm'),
          panel.background = element_blank())
}
gridExtra::grid.arrange(grobs = wls, ncol = 5)

### Figure 3 (b) plot
mi_name <- vector(length = 10)
n <- 0
for (i in 1:10) {
  name <- paste0('m',i)
  mi_name[i] <- name
}
propdf <- data.frame(mi_name=rep(mi_name,each=2))
a <- 0
for (age in 19:97) {
  fpath <- paste0('/output/age19-97/w',age,'.csv')
  if (file.exists(fpath) == F) {
    next
  } else {
    w <- read.table(fpath,sep = ',')
    a <- a+1
    avec <- vector(length = 2*10)
    for (i in 1:10) {
      mi <- w[,i]
      c0 <- c(0.05,0.1,0.15,0.15,0.1,0.1,0.2,0.1,0.15,0.2)
      i0 <- length(which(mi<c0[i])) #'0' means <c0;'1' means >=c0.
      i1 <- length(which(mi>=c0[i]))
      avec[(2*i-1):(2*i)] <- c(i0,i1)/length(mi)
    }
    propdf <- cbind(propdf,data.frame(avec))
    colnames(propdf)[a+1] <- paste0('age',age)
  }
}
mc <- c(2,3,5,6,8)
# mc <- c(1,4,7,9,10) # for supplementary figure
pls <- list()
k <- 0
for (i in mc) {
  k <- k+1
  m <- data.frame(age = c(19:91,93:97),value0=t(propdf[(2*i-1),-1]),
                  value1=t(propdf[(2*i),-1]))
  colnames(m)[2] <- 'value0'
  colnames(m)[3] <- 'value1'
  p1 <- ggplot(m) + 
    geom_point(aes(x = age,y = value0),color = "#715EA9", size = 1.5, shape = 16) +
    coord_cartesian(ylim = c(0, 1))+
    labs(title = paste0('module',i)) +
    theme_bw()+
    theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
          axis.text = element_text(size = 12),
          axis.title = element_blank(),
          legend.position = 'none',
          axis.ticks.length = unit(2,'mm'),
          panel.background = element_blank())
  pls[[k]] <- p1
  # pls[[2*k-1]] <- p1
  AG = rep(c(paste0('g',1:4)),c(22,20,20,16))
  m <- data.frame(age = c(19:91,93:97), prop = t(propdf[(2*i-1),-1]),
                  Group = AG)
  colnames(m)[2] <- 'prop'
  p2 <- ggplot(m,aes(x=Group,y=prop,color=Group)) + 
    geom_half_boxplot(fill=c("#DCC7E1", "#BBA1CB","#A67EB7","#7D5284"),color='#181717',
                      side = 'r',errorbar.draw=F, width=0.6,linewidth=0.4)+
    geom_half_point_panel(aes(color=Group),side = 'l',shape=16,size=1.2,
                          range_scale = 1.2)+
    scale_color_manual(values=c("#DCC7E1", "#BBA1CB","#A67EB7","#7D5284"))+
    coord_cartesian(ylim = c(0, 1))+
    theme_bw()+
    theme(plot.title = element_blank(),
          axis.text = element_text(size = 12),
          axis.title = element_blank(),
          legend.position = 'none',
          axis.ticks.length = unit(2,'mm'),
          panel.background = element_blank()) 
  pls[[k+5]] <- p2
  # pls[[2*k]] <- p2
}
gridExtra::grid.arrange(grobs = pls, ncol = 5)

### Figure 3 (c) plot
dim1 <- c(1,3,4)
dim2 <- c(2,7,10)
ylim <- c(0.28,0.6,0.72)
pls <- list()
for (n in 1:length(dim1)) {
  a <- dim1[n]
  b <- dim2[n]
  m <- data.frame(modulea = W[,a],moduleb=W[,b])
  part <- vector(length = 6606)
  for (i in 1:6606) {
    if (m[i,1]<c0[a] & m[i,2]<c0[b]) {
      part[i] = '00'
    } else if (m[i,1]<c0[a] & m[i,2]>=c0[b]) {
      part[i] = '01'
    } else if (m[i,1]>=c0[a] & m[i,2]>=c0[b]) {
      part[i] = '11'
    } else { part[i] = '10' }
  }
  m$part <- part
  pls[[n]] <- ggplot(m, aes(x = modulea,y = moduleb, color = part)) + 
    geom_point(alpha=0.8,size = 1, shape = 16) +
    scale_color_manual(values = c("00" = "#159A80", "01" = "#715EA9",
                                  "10" = "#F28147", "11" = "#3371B3"))+
    geom_vline(xintercept = c0[a], linetype = "solid", color = "black", 
               alpha = 0.7,linewidth = 0.7) +
    geom_hline(yintercept = c0[b], linetype = "solid", color = "black", 
               alpha = 0.8,linewidth = 0.7) +
    labs(x = paste0("module",a), y = paste0("module",b))+
    theme_bw()+
    theme(plot.title = element_blank(),
          axis.text = element_text(size = 12,color = 'black'),
          axis.title = element_text(size = 12, face = "bold", hjust = 0.5),
          legend.position = 'none',
          axis.ticks.length = unit(2,'mm'),
          panel.background = element_blank())
}
gridExtra::grid.arrange(grobs = pls, nrow = 3)

### Figure 3 (d) plot
mij_name <- vector(length = 45)
n <- 0
for (i in 1:10) {
  for (j in 1:10) {
    if (i<j) {
      n <- n+1
      name <- paste0('m',i,j)
      mij_name[n] <- name
    }
  }
}
propdf <- data.frame(mij_name=rep(mij_name,each=4))
a <- 0
for (age in 19:97) {
  fpath <- paste0('/output/age19-97/w',age,'.csv')
  if (file.exists(fpath) == F) {
    next
  } else {
    w <- read.table(fpath,sep = ',')
    a <- a+1
    avec <- vector(length = 4*45)
    n <- 0
    for (i in 1:10) {
      for (j in 1:10) {
        if (i<j) {
          mi <- w[,i]
          mj <- w[,j]
          c0 <- c(0.05,0.1,0.15,0.15,0.1,0.1,0.2,0.1,0.15,0.2)
          ij00 <- length(which(mi<c0[i] & mj<c0[j])) #'0' means <c0;'1' means >=c0.
          ij01 <- length(which(mi<c0[i] & mj>=c0[j]))
          ij10 <- length(which(mi>=c0[i] & mj<c0[j]))
          avec[(4*n+1):(4*n+3)] <- c(ij00,ij01,ij10)/length(mi)
          avec[4*(n+1)] <- 1-sum(avec[(4*n+1):(4*n+3)]) #ensure the sum of 4 parts equal to 1.
          n <- n+1
        }
      }
    }
    propdf <- cbind(propdf,data.frame(avec))
    colnames(propdf)[a+1] <- paste0('age',age)
  }
}
# figure 3(d) left panel Module3 & Module7:
mijdf <- distp[117:120,2:79]
mijdf$parts <- c('ij00','ij01','ij10','ij11')
mijdf <- melt(mijdf,id.vars = 'parts')
mijdf$age <- rep(c(19:91,93:97),each=4)
AG <- rep(c(paste0('g',1:4)),c(22,20,20,16))
mijdf$AG <- rep(AG,each=4)
mijdf$parts <- factor(mijdf$parts,levels = c('ij00','ij01','ij10','ij11'))
mijdf$AG <- factor(mijdf$AG,levels = c('g1','g2','g3','g4'))
fig3d_l2 <- ggplot(mijdf,aes(x=age,y=value,fill=parts))+
  geom_bar(alpha=0.75,stat = 'identity',position = 'stack')+
  labs(title = 'module4 and module10',fill='Parts')+
  scale_fill_manual(values = c("ij00" = "#159A80", "ij01" = "#715EA9",
                               "ij10" = "#F28147", "ij11" = "#3371B3"))+
  theme_bw()+
  theme(plot.title = element_text(size=14, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 12,color = 'black'),
        axis.title = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.length = unit(2,'mm'),
        axis.text.x = element_blank(),
        legend.position = 'none',
        strip.text.x = element_text(size = 12),
        strip.background = element_blank(), strip.placement = 'inside',
        panel.spacing = unit(0.01,'lines'),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())+
  facet_grid(~AG,scales = 'free',switch = 'both',margins = F)
# figure 3(d) right panel Module4 & Module10
AG = rep(c(paste0('g',1:4)),c(22,20,20,16))
m <- data.frame(age = c(19:91,93:97),value=t(distp[117,-1]),Group = AG)
colnames(m)[2] <- 'value'
fig3d_r3 <- ggplot(m,aes(x=Group,y=value,color=Group)) + 
  geom_half_boxplot(fill=c("#C6E3C6", "#6AB97D","#3F8D56","#306E3E"),color='#181717',
                    side = 'r',errorbar.draw=F, width=0.6,linewidth=0.4)+
  geom_half_point_panel(aes(color=Group),side = 'l',shape=16,size=1.2,
                        range_scale = 1.2)+
  scale_color_manual(values=c("#C6E3C6", "#6AB97D","#3F8D56","#306E3E"))+
  coord_cartesian(ylim = c(0.65, 1))+
  labs(title = 'module4 and module10')+
  theme_bw()+
  theme(plot.title = element_text(size=14, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 12,color = 'black'),
        axis.title = element_blank(),
        axis.ticks.length = unit(2,'mm'),
        legend.position = 'none',
        panel.background = element_blank())

###### Other Notes：
# The other subplots in Figure 3 are generated by similar plot code and will not be repeated here.