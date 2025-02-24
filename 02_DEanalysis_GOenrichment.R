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

# do DE analysis between iNK (UCH-HSPC-iNK) and (UCB-NK, ESC-iNK, iPSC-iNK, CD19-CAR-iNK)  
## load rds files
UCB_HSPC_iNK.seu=readRDS("~/Documents/LJH_10X/results/data/UCB_NK_HSPC_iNK_aggr.rds")

## UCB-NK vesus iNK (HSPC_iNK)  

UCBNK_iNK=FindMarkers(UCB_HSPC_iNK.seu,ident.1 ="UCB_NK",ident.2 = "HSPC_iNK" ,test.use = "MAST",min.pct = 0.25,logfc.threshold = 0.5,group.by = "group")
dim(UCBNK_iNK)
UCBNK_iNK=UCBNK_iNK[UCBNK_iNK$p_val_adj<0.05,]

## plot heatmap  
UCBNK_iNK$gene=rownames(UCBNK_iNK)
UCBNK_iNK$cluster[UCBNK_iNK$avg_log2FC>0]="UCB_NK"
UCBNK_iNK$cluster[UCBNK_iNK$avg_log2FC<0]="HSPC_iNK"
top20 <- UCBNK_iNK %>% group_by(cluster) %>% top_n(n = 20, wt = abs(avg_log2FC))
top20=arrange(top20, cluster,-abs(avg_log2FC))
bk <- c(seq(0,0.5,by=0.05),seq(0.05,6,by=0.05))

UCB_HSPC_iNK.seu_1=subset(UCB_HSPC_iNK.seu,idents=c("UCB_NK","HSPC_iNK"))
UCB_HSPC_iNK.seu_1$group1=factor(UCB_HSPC_iNK.seu_1$group1,levels=c("UCB_NK","HSPC_iNK"))

p=DoHeatmap(UCB_HSPC_iNK.seu_1, features = top20$gene,label=F,slot = "scale.data",group.by = "group1",group.colors = c("#F8766D","#7CAE00"))+theme(axis.text.y = element_text(face="bold",size=15,colour = "black"))+scale_fill_gradientn(colours = c(colorRampPalette(colors = c("#5252A9","white"))(length(bk)/2),                                                                                                                                                                                                                                                    colorRampPalette(colors = c("white","#DA6161"))(length(bk)/2)))#+NoLegend()
p

## UCB-NK vesus CD19-CAR-iNK (HSPC_CAR19_iNK)  
UCBNK_CAR19=FindMarkers(UCB_HSPC_iNK.seu,ident.1 ="UCB_NK",ident.2 = "HSPC_CAR19_iNK" ,test.use = "MAST",min.pct = 0.25,logfc.threshold = 0.5,group.by = "group")
dim(UCBNK_CAR19)
UCBNK_CAR19=UCBNK_CAR19[UCBNK_CAR19$p_val_adj<0.05,]
## plot heatmap  
UCBNK_CAR19$gene=rownames(UCBNK_CAR19)
UCBNK_CAR19$cluster[UCBNK_CAR19$avg_log2FC>0]="UCB_NK"
UCBNK_CAR19$cluster[UCBNK_CAR19$avg_log2FC<0]="HSPC_CAR19_iNK"
top20 <- UCBNK_CAR19 %>% group_by(cluster) %>% top_n(n = 20, wt = abs(avg_log2FC))
top20=arrange(top20, cluster,-abs(avg_log2FC))
bk <- c(seq(0,0.5,by=0.05),seq(0.05,6,by=0.05))

UCB_HSPC_iNK.seu$group1=UCB_HSPC_iNK.seu$group
UCB_HSPC_iNK.seu_1=subset(UCB_HSPC_iNK.seu,idents=c("UCB_NK","HSPC_CAR19_iNK"))
UCB_HSPC_iNK.seu_1$group1=factor(UCB_HSPC_iNK.seu_1$group1,levels=c("UCB_NK","HSPC_CAR19_iNK"))

p=DoHeatmap(UCB_HSPC_iNK.seu_1, features = top20$gene,label=F,slot = "scale.data",group.by = "group1",group.colors = c("#F8766D","#00BFC4"))+theme(axis.text.y = element_text(face="bold",size=15,colour = "black"))+scale_fill_gradientn(colours = c(colorRampPalette(colors = c("#5252A9","white"))(length(bk)/2),
                                                                                                                                                                                                                                                      colorRampPalette(colors = c("white","#DA6161"))(length(bk)/2)))#+NoLegend()
p
#ggsave(p,filename ="~/Documents/LJH_10X/results/plots/top20_UCBNK_vs_CAR19iNK.pdf",width = 6,height = 8)



## CD19-CAR-iNK (HSPC_CAR19_iNK) vesus iNK (HSPC_iNK)  
CAR19iNK_iNK=FindMarkers(UCB_HSPC_iNK.seu,ident.1 ="HSPC_CAR19_iNK",ident.2 = "HSPC_iNK" ,test.use = "MAST",min.pct = 0.25,logfc.threshold = 0.5,group.by = "group")
dim(CAR19iNK_iNK)
CAR19iNK_iNK=CAR19iNK_iNK[CAR19iNK_iNK$p_val_adj<0.05,]

## plot heatmap
CAR19iNK_iNK$gene=rownames(CAR19iNK_iNK)
CAR19iNK_iNK$cluster[CAR19iNK_iNK$avg_log2FC>0]="HSPC_CAR19_iNK"
CAR19iNK_iNK$cluster[CAR19iNK_iNK$avg_log2FC<0]="HSPC_iNK"
top20 <- CAR19iNK_iNK %>% group_by(cluster) %>% top_n(n = 20, wt = abs(avg_log2FC))
top20=arrange(top20, cluster,-abs(avg_log2FC))
bk <- c(seq(0,0.5,by=0.05),seq(0.05,6,by=0.05))


UCB_HSPC_iNK.seu_1=subset(UCB_HSPC_iNK.seu,idents=c("HSPC_CAR19_iNK","HSPC_iNK"))
UCB_HSPC_iNK.seu_1$group1=factor(UCB_HSPC_iNK.seu_1$group1,levels=c("HSPC_iNK","HSPC_CAR19_iNK"))
#UCB_NK$label_f=factor(UCB_NK$label,levels = c("FCER1G_h UCB_rNK","FCER1G_h UCB_aNK","KLRC2_h UCB_rNK",  "KLRC2_h UCB_aNK","CD3G_h UCB_rNK","C3 UCB_rNK","C5 UCB_rNK","C6 UCB_rNK", "C3 UCB_aNK","C4 UCB_aNK","C5 UCB_aNK" ))


p=DoHeatmap(UCB_HSPC_iNK.seu_1, features = top20$gene,label=F,slot = "scale.data",group.by = "group1",group.colors = c("#7CAE00","#00BFC4"))+theme(axis.text.y = element_text(face="bold",size=15,colour = "black"))+scale_fill_gradientn(colours = c(colorRampPalette(colors = c("#5252A9","white"))(length(bk)/2),                                                                                                                                                                                                                                                     colorRampPalette(colors = c("white","#DA6161"))(length(bk)/2)))#+NoLegend()
p

#ggsave(p,filename ="~/Documents/LJH_10X/results/plots/top20_HSPCiNK_vsCAR19iNK.pdf",width = 6,height = 8)



## ESC_iNK vesus iNK (HSPC_iNK)  
### merge data  
ESC_iNK=readRDS("~/Documents/iNK_comp/ESC_OL.rds")
ESC_iNK$orig.ident="ESC_iNK"

iPS_iNK=readRDS("~/Documents/iNK_comp/iPSC_OL.rds")
iPS_iNK$orig.ident="iPSC_iNK"
UCB_HSPC_iNK_ESC.seu=merge(UCB_HSPC_iNK.seu,ESC_iNK,iPS_iNK)
### scale data to do dimension reduction analysis
### scaled data saved in UCB_HSPC_iNK.seu[["RNA"]]@scale.data
all.genes <- rownames(UCB_HSPC_iNK_ESC.seu)
UCB_HSPC_iNK_ESC.seu <- ScaleData(UCB_HSPC_iNK_ESC.seu, features = all.genes)
saveRDS(UCB_HSPC_iNK_ESC.seu,file = "~/Documents/LJH_10X/results/data/UCB_HSPC_iNK_ESCiPSC.rds")
### DE analysis  
ESCiNK_HSPCiNK=FindMarkers(UCB_HSPC_iNK_ESC.seu,ident.1 ="ESC_iNK",ident.2 = "HSPC_iNK" ,test.use = "MAST",min.pct = 0.25,logfc.threshold = 0.5,group.by = "group")
dim(ESCiNK_HSPCiNK)
ESCiNK_HSPCiNK=ESCiNK_HSPCiNK[ESCiNK_HSPCiNK$p_val_adj<0.05,]

### plot heatmap  
ESCiNK_HSPCiNK$gene=rownames(ESCiNK_HSPCiNK)
ESCiNK_HSPCiNK$cluster[ESCiNK_HSPCiNK$avg_log2FC>0]="ESC_iNK"
ESCiNK_HSPCiNK$cluster[ESCiNK_HSPCiNK$avg_log2FC<0]="HSPC_iNK"
top20 <- ESCiNK_HSPCiNK %>% group_by(cluster) %>% top_n(n = 20, wt = abs(avg_log2FC))
top20=arrange(top20, cluster,-abs(avg_log2FC))
bk <- c(seq(0,0.5,by=0.05),seq(0.05,6,by=0.05))

UCB_HSPC_iNK_ESC.seu$group1=UCB_HSPC_iNK_ESC.seu$group
UCB_HSPC_iNK.seu_1=subset(UCB_HSPC_iNK_ESC.seu,idents=c("ESC_iNK","HSPC_iNK"))
UCB_HSPC_iNK.seu_1$group1=factor(UCB_HSPC_iNK.seu_1$group1,levels=c("ESC_iNK","HSPC_iNK"))

p=DoHeatmap(UCB_HSPC_iNK.seu_1, features = top20$gene,label=F,slot = "scale.data",group.by = "group1",group.colors = c("#C77CFF","#7CAE00"))+theme(axis.text.y = element_text(face="bold",size=15,colour = "black"))+scale_fill_gradientn(colours = c(colorRampPalette(colors = c("#5252A9","white"))(length(bk)/2),
                                                                                                                                                                                                                                                      colorRampPalette(colors = c("white","#DA6161"))(length(bk)/2)))#+NoLegend()
p
#ggsave(p,filename ="~/Documents/LJH_10X/results/plots/top20_ESCiNK_vs_HSPCiNK.pdf",width = 6,height = 8)


##iPSC_iNK vs iNK (HSPC_iNK)   

iPSiNK_HSPCiNK=FindMarkers(UCB_HSPC_iNK_ESC.seu,ident.1 ="iPSC_iNK",ident.2 = "HSPC_iNK" ,test.use = "MAST",min.pct = 0.25,logfc.threshold = 0.5,group.by = "group")
dim(iPSiNK_HSPCiNK)
iPSiNK_HSPCiNK=iPSiNK_HSPCiNK[iPSiNK_HSPCiNK$p_val_adj<0.05,]
### plot heatmap  
iPSiNK_HSPCiNK$gene=rownames(iPSiNK_HSPCiNK)
iPSiNK_HSPCiNK$cluster[iPSiNK_HSPCiNK$avg_log2FC>0]="iPSC_iNK"
iPSiNK_HSPCiNK$cluster[iPSiNK_HSPCiNK$avg_log2FC<0]="HSPC_iNK"
top20 <- iPSiNK_HSPCiNK %>% group_by(cluster) %>% top_n(n = 20, wt = abs(avg_log2FC))
top20=arrange(top20, cluster,-abs(avg_log2FC))
bk <- c(seq(0,0.5,by=0.05),seq(0.05,6,by=0.05))

UCB_HSPC_iNK_ESC.seu$group1=UCB_HSPC_iNK_ESC.seu$group
UCB_HSPC_iNK.seu_1=subset(UCB_HSPC_iNK_ESC.seu,idents=c("iPSC_iNK","HSPC_iNK"))
UCB_HSPC_iNK.seu_1$group1=factor(UCB_HSPC_iNK.seu_1$group1,levels=c("iPSC_iNK","HSPC_iNK"))
p=DoHeatmap(UCB_HSPC_iNK.seu_1, features = top20$gene,label=F,slot = "scale.data",group.by = "group1",group.colors = c("#C77CFF","#7CAE00"))+theme(axis.text.y = element_text(face="bold",size=15,colour = "black"))+scale_fill_gradientn(colours = c(colorRampPalette(colors = c("#5252A9","white"))(length(bk)/2),
                                                                                                                                                                                                                                                      colorRampPalette(colors = c("white","#DA6161"))(length(bk)/2)))#+NoLegend()
p
ggsave(p,filename ="~/Documents/LJH_10X/results/plots/top20_iPSCiNK_vs_HSPCiNK.pdf",width = 6,height = 8)


# do GO enrichment for DEGS  

## GO-enrichment function
# the function to get GO term information as gene list input
GO_inbatch <- function(deg_res,clustern,db){
  # get gene symbol of up-regulation in one cluster
  subset(deg_res, cluster == clustern) %>% # get DEGs information of one cluster
    rownames %>% # get gene symbol
    bitr(.,fromType="SYMBOL",toType="ENTREZID",OrgDb = db) %$%ENTREZID %>%
    enrichGO(gene=.,
             Org = db,
             ont="BP",
             pvalueCutoff=0.05,
             pAdjustMethod="BH",
             readable="TRUE") # GO analysis
}

## GO-enrichment

library(org.Hs.eg.db)
UCBNK_iNK_go <- sapply(unique(UCBNK_iNK$cluster), function(clustern){GO_inbatch(UCBNK_iNK, clustern,"org.Hs.eg.db")}, USE.NAMES = TRUE)
UCBNK_CAR19_go <- sapply(unique(UCBNK_CAR19$cluster), function(clustern){GO_inbatch(UCBNK_CAR19, clustern,"org.Hs.eg.db")}, USE.NAMES = TRUE)
CAR19iNK_iNK_go <- sapply(unique(CAR19iNK_iNK$cluster), function(clustern){GO_inbatch(CAR19iNK_iNK, clustern,"org.Hs.eg.db")}, USE.NAMES = TRUE)
ESCiNK_HSPCiNK_go <- sapply(unique(ESCiNK_HSPCiNK$cluster), function(clustern){GO_inbatch(ESCiNK_HSPCiNK, clustern,"org.Hs.eg.db")}, USE.NAMES = TRUE)
iPSiNK_HSPCiNK_go <- sapply(unique(iPSiNK_HSPCiNK$cluster), function(clustern){GO_inbatch(iPSiNK_HSPCiNK, clustern,"org.Hs.eg.db")}, USE.NAMES = TRUE)

## plot the results  
p1=enrichplot::dotplot(UCBNK_iNK_go[[1]],title =names(UCBNK_iNK_go)[1],showCategory=15,font.size = 20)+
  theme(
    axis.text.x=element_text(colour="black",size=8,face="bold"),
    axis.title.x =element_text(colour="black",size=12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5,size=20),legend.position = "bottom")+guides(size = guide_legend(nrow = 3, byrow = TRUE,order = 1,label.theme = element_text(size=8)),color=guide_colorbar(direction="vertical",barwidth = 0.5,barheight = 3,label.theme = element_text(size=8,angle = 0),order = 2))

p2=enrichplot::dotplot(UCBNK_iNK_go[[2]],title =names(UCBNK_iNK_go)[2],showCategory=15,font.size = 20)+
  theme(
    axis.text.x=element_text(colour="black",size=8,face="bold"),
    axis.title.x =element_text(colour="black",size=12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5,size=20),legend.position = "bottom")+ guides(size = guide_legend(nrow = 3, byrow = TRUE,order = 1,label.theme = element_text(size=8)),color=guide_colorbar(direction="vertical",barwidth = 0.5,barheight = 3,label.theme = element_text(size=8,angle = 0),order = 2))

p=plot_grid(p1,p2,ncol = 2,rel_widths = c(1,1))+theme(
  plot.margin = margin(t = 10,  # 顶部边缘距离
                       r = 20,  # 右边边缘距离
                       b = 10,  # 底部边缘距离
                       l = 20)) # 左边边缘距离

p
#ggsave(p,filename ="~/Documents/LJH_10X/results/plots/Goterms_UCBNK_HSPCiNK.pdf",width = 15,height = 10)

p1=enrichplot::dotplot(UCBNK_CAR19_go[[1]],title =names(UCBNK_CAR19_go)[1],showCategory=15,font.size = 20)+
  theme(
    axis.text.x=element_text(colour="black",size=8,face="bold"),
    axis.title.x =element_text(colour="black",size=12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5,size=20),legend.position = "bottom")+guides(size = guide_legend(nrow = 3, byrow = TRUE,order = 1,label.theme = element_text(size=8)),color=guide_colorbar(direction="vertical",barwidth = 0.5,barheight = 3,label.theme = element_text(size=8,angle = 0),order = 2))

p2=enrichplot::dotplot(UCBNK_CAR19_go[[2]],title =names(UCBNK_CAR19_go)[2],showCategory=15,font.size = 20)+
  theme(
    axis.text.x=element_text(colour="black",size=8,face="bold"),
    axis.title.x =element_text(colour="black",size=12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5,size=20),legend.position = "bottom")+ guides(size = guide_legend(nrow = 3, byrow = TRUE,order = 1,label.theme = element_text(size=8)),color=guide_colorbar(direction="vertical",barwidth = 0.5,barheight = 3,label.theme = element_text(size=8,angle = 0),order = 2))

p=plot_grid(p1,p2,ncol = 2,rel_widths = c(1,1))+theme(
  plot.margin = margin(t = 10,  # 顶部边缘距离
                       r = 20,  # 右边边缘距离
                       b = 10,  # 底部边缘距离
                       l = 20)) # 左边边缘距离

p
#ggsave(p,filename ="~/Documents/LJH_10X/results/plots/Goterms_UCBNK_CAR19iNK.pdf",width = 15,height = 10)


p1=enrichplot::dotplot(CAR19iNK_iNK_go[[1]],title =names(CAR19iNK_iNK_go)[1],showCategory=15,font.size = 20)+
  theme(
    axis.text.x=element_text(colour="black",size=8,face="bold"),
    axis.title.x =element_text(colour="black",size=12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5,size=20),legend.position = "bottom")+guides(size = guide_legend(nrow = 3, byrow = TRUE,order = 1,label.theme = element_text(size=8)),color=guide_colorbar(direction="vertical",barwidth = 0.5,barheight = 3,label.theme = element_text(size=8,angle = 0),order = 2))

p2=enrichplot::dotplot(CAR19iNK_iNK_go[[2]],title =names(CAR19iNK_iNK_go)[2],showCategory=15,font.size = 20)+
  theme(
    axis.text.x=element_text(colour="black",size=8,face="bold"),
    axis.title.x =element_text(colour="black",size=12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5,size=20),legend.position = "bottom")+ guides(size = guide_legend(nrow = 3, byrow = TRUE,order = 1,label.theme = element_text(size=8)),color=guide_colorbar(direction="vertical",barwidth = 0.5,barheight = 3,label.theme = element_text(size=8,angle = 0),order = 2))

p=plot_grid(p1,p2,ncol = 2,rel_widths = c(1,1))+theme(
  plot.margin = margin(t = 10,  # 顶部边缘距离
                       r = 20,  # 右边边缘距离
                       b = 10,  # 底部边缘距离
                       l = 20)) # 左边边缘距离

p
#ggsave(p,filename ="~/Documents/LJH_10X/results/plots/Goterms_CAR19iNK_iNK.pdf",width = 15,height = 10)

p1=enrichplot::dotplot(ESCiNK_HSPCiNK_go[[1]],title =names(ESCiNK_HSPCiNK_go)[1],showCategory=15,font.size = 20)+
  theme(
    axis.text.x=element_text(colour="black",size=8,face="bold"),
    axis.title.x =element_text(colour="black",size=12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5,size=20),legend.position = "bottom")+guides(size = guide_legend(nrow = 3, byrow = TRUE,order = 1,label.theme = element_text(size=8)),color=guide_colorbar(direction="vertical",barwidth = 0.5,barheight = 3,label.theme = element_text(size=8,angle = 0),order = 2))

p2=enrichplot::dotplot(ESCiNK_HSPCiNK_go[[2]],title =names(ESCiNK_HSPCiNK_go)[2],showCategory=15,font.size = 20)+
  theme(
    axis.text.x=element_text(colour="black",size=8,face="bold"),
    axis.title.x =element_text(colour="black",size=12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5,size=20),legend.position = "bottom")+ guides(size = guide_legend(nrow = 3, byrow = TRUE,order = 1,label.theme = element_text(size=8)),color=guide_colorbar(direction="vertical",barwidth = 0.5,barheight = 3,label.theme = element_text(size=8,angle = 0),order = 2))

p=plot_grid(p1,p2,ncol = 2,rel_widths = c(1,1))+theme(
  plot.margin = margin(t = 10,  # 顶部边缘距离
                       r = 20,  # 右边边缘距离
                       b = 10,  # 底部边缘距离
                       l = 20)) # 左边边缘距离

p
#ggsave(p,filename ="~/Documents/LJH_10X/results/plots/Goterms_ESCiNK_HSPCiNK.pdf",width = 15,height = 12)

p1=enrichplot::dotplot(iPSiNK_HSPCiNK_go[[1]],title =names(iPSiNK_HSPCiNK_go)[1],showCategory=15,font.size = 20)+
  theme(
    axis.text.x=element_text(colour="black",size=8,face="bold"),
    axis.title.x =element_text(colour="black",size=12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5,size=20),legend.position = "bottom")+guides(size = guide_legend(nrow = 3, byrow = TRUE,order = 1,label.theme = element_text(size=8)),color=guide_colorbar(direction="vertical",barwidth = 0.5,barheight = 3,label.theme = element_text(size=8,angle = 0),order = 2))

p2=enrichplot::dotplot(iPSiNK_HSPCiNK_go[[2]],title =names(iPSiNK_HSPCiNK_go)[2],showCategory=15,font.size = 20)+
  theme(
    axis.text.x=element_text(colour="black",size=8,face="bold"),
    axis.title.x =element_text(colour="black",size=12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5,size=20),legend.position = "bottom")+ guides(size = guide_legend(nrow = 3, byrow = TRUE,order = 1,label.theme = element_text(size=8)),color=guide_colorbar(direction="vertical",barwidth = 0.5,barheight = 3,label.theme = element_text(size=8,angle = 0),order = 2))

p=plot_grid(p1,p2,ncol = 2,rel_widths = c(1,1))+theme(
  plot.margin = margin(t = 10,  # 顶部边缘距离
                       r = 20,  # 右边边缘距离
                       b = 10,  # 底部边缘距离
                       l = 20)) # 左边边缘距离

p
ggsave(p,filename ="~/Documents/LJH_10X/results/plots/Goterms_iPSiNK_HSPCiNK.pdf",width = 15,height = 12)


### save results and ouput 
save(CAR19iNK_iNK,UCBNK_iNK,UCBNK_CAR19,ESCiNK_HSPCiNK,iPSiNK_HSPCiNK,CAR19iNK_iNK_go,UCBNK_iNK_go,UCBNK_CAR19_go,ESCiNK_HSPCiNK_go,iPSiNK_HSPCiNK_go,file = "~/Documents/LJH_10X/results/data/DEGs_GO_250217.RData")


library(openxlsx)
sheets=list("UCBNK vs HSPCiNK"=UCBNK_iNK,"ESCiNK vs HSPCiNK"=ESCiNK_HSPCiNK,"iPSCiNK vs HSPCiNK"=iPSiNK_HSPCiNK,"HSPCiNK vs CAR19iNK"=CAR19iNK_iNK)
write.xlsx(sheets,"~/Documents/LJH_10X/results/data/DEGs_250217.xlsx")


library(openxlsx)
sheets=list("UCBNK"=UCBNK_iNK_go[[1]]@result,"HSPC_iNK"=UCBNK_iNK_go[[2]]@result)
write.xlsx(sheets,"~/Documents/LJH_10X/results/data/DEGs_GOterms_UCBNK_HSPCiNK.xlsx")

sheets=list("ESC_iNK"=ESCiNK_HSPCiNK_go[[2]]@result,"HSPC_iNK"=ESCiNK_HSPCiNK_go[[1]]@result)
write.xlsx(sheets,"~/Documents/LJH_10X/results/data/DEGs_GOterms_ESCiNK_HSPCiNK.xlsx")

sheets=list("iPSiNK"=iPSiNK_HSPCiNK_go[[2]]@result,"HSPC_iNK"=iPSiNK_HSPCiNK_go[[1]]@result)
write.xlsx(sheets,"~/Documents/LJH_10X/results/data/DEGs_GOterms_iPSCiNK_HSPCiNK.xlsx")

sheets=list("HSPC_CAR19_iNK"=CAR19iNK_iNK_go[[2]]@result,"HSPC_iNK"=CAR19iNK_iNK_go[[1]]@result)
write.xlsx(sheets,"~/Documents/LJH_10X/results/data/DEGs_GOterms_CAR19iNK_HSPCiNK.xlsx")
