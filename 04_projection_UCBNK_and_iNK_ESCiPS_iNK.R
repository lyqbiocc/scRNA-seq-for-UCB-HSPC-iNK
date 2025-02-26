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

#projection of UCB-NK, ESC_iNK, iPSC_iNK, iNK(HSPC_iNK)  
## input data 
UCB_HSPC_iNK_ESC.seu=readRDS("~/Documents/LJH_10X/results/data/UCB_HSPC_iNK_ESCiPSC.rds")
UCB_HSPC_iNK_ESC.seu=subset(UCB_HSPC_iNK_ESC.seu,subset=group!="CAR19_HSPC_iNK")

## transform merge  
### FindIntegrationAnchors（ no SCT normalized）

pb_anchor <- FindIntegrationAnchors(object.list = list(subset(UCB_HSPC_iNK_ESC.seu,ident="UCB_NK"),subset(UCB_HSPC_iNK_ESC.seu,ident=c("HSPC_iNK","ESC_iNK","iPSC_iNK"))), dims = 1:30,k.anchor = 5,k.filter = 150)
#pb_integ <- FindIntegrationAnchors(object.list = list(ILC_induced,ILC_ref.obj2), dims = 1:30,k.anchor = 5,k.filter = 180)
pb_integ <- IntegrateData(anchorset = pb_anchor, dims = 1:30)
DefaultAssay(pb_integ) <- "integrated"
# Run the standard workflow for visualization and clustering
pb_integ <- ScaleData(pb_integ, verbose = FALSE)
pb_integ <- RunPCA(pb_integ, features = VariableFeatures(pb_integ),npcs = 30, verbose = FALSE)
ElbowPlot(pb_integ,ndims = 30)

pb_integ <- RunUMAP(pb_integ, reduction = "pca", dims = 1:30,n.neighbors =100)

pb_integ$group2=factor(pb_integ$group2,levels=c("UCB_NK","ESC_iNK","iPSC_iNK","HSPC_iNK"))
p1 <- DimPlot(pb_integ, reduction = "umap", label = F,group.by = "group",pt.size = 0.0001)+labs(title = "")&theme_bw() &
  theme(
    axis.text.x=element_text(colour="black",size=8,face="bold"),
    axis.text.y=element_text(colour="black",size=8,face="bold"),
    axis.title.x=element_text(size = 8,face="bold"),
    axis.title.y=element_text(size = 8,face="bold"),
    panel.border = element_rect(),
    axis.line = element_line(colour = "black",size=0),
    legend.text=element_text(face="bold", colour="black",size=8),
    legend.title=element_text(face="bold",colour="black",size=8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face="bold.italic",hjust = 0.5,size=10)) &#& NoLegend()
  scale_fill_manual(values = c("#F8766D","#00B0F6","#00BF7D","#E76BF3"))#&geom_jitter(size=0.001)
#p2+p1
p1

#ggsave(p1,filename ="~/Documents/LJH_10X/results/plots/UMAP_NKdasets.pdf",width = 6,height = 5)
#saveRDS(pb_integ,file = "~/Documents/LJH_10X/results/data/UCBNK_ESC_iPS_HSPC_iNK_integ_250219.rds")

