# library packages

library(ggplot2)
library(dplyr)
library(Seurat)
library(patchwork)
library(viridis)
library(magrittr)
library(cowplot)
library(clusterProfiler)
library(org.Hs.eg.db)

# expression of NK-related genes  
UCB_HSPC_iNK_ESC.seu=readRDS("~/Documents/LJH_10X/results/data/UCB_HSPC_iNK_ESCiPSC.rds")

genes=c("SLAMF7","NCR3","NCR2","KLRK1","NCAM1","CD69","KLRD1","KLRC1","CD96","NCR1","TNFSF10","FASLG")

p1 = VlnPlot(UCB_HSPC_iNK_ESC.seu, features = genes, pt.size = 0,ncol = 3) & 
  labs(x="",y="LogNormalized Expression") & 
  theme(plot.title = element_text(face = "italic",size = 10),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_text(size = 8)) &
  scale_fill_manual(values = colors[c(5,10)])#&geom_jitter(size=0.001)
p1
#ggsave(p1,filename = "./LJH/NK-Related-exp.pdf",width = 6,height = 9)                                                                                                                                                                                                                                       panel.grid.major = element_blank(),panel.grid.minor = element_blank())

genes=c("KIR2DL1","KIR2DL3","KIR3DL1","KIR3DL2")

p1 = VlnPlot(UCB_HSPC_iNK_ESC.seu, features = genes, pt.size = 0.001,ncol = 4) & 
  labs(x="",y="LogNormalized Expression") & 
  theme(plot.title = element_text(face = "italic",size = 10),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_text(size = 8)) &
  scale_fill_manual(values = colors[c(5,10)])#&geom_jitter(size=0.001)
p1