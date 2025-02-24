
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


# Projection of UCB-NK, iNK and CD19-CAR-iNK  
## input data  
### read 10x data after aggregation (UCB-NK,iNK, CD19-CAR-iNK)

data.seu <- Read10X(data.dir = "/public/home/lin_yunqing/Documents/LJH_10X/HSPC_iNK_aggr_V7/outs/count/filtered_feature_bc_matrix/")
meta=read.delim("/public/home/lin_yunqing/Documents/LJH_10X/HSPC_iNK_aggr_V7/outs//aggregation.csv",header = T,sep=",")

data_meta=data.frame(bcnum=colnames(data.seu),reshape2::colsplit(colnames(data.seu), "[-]", names = c("bc","num")),group="")
table(data_meta$num)

data_meta[data_meta$num==1,]$group=meta$sample_id[1]
data_meta[data_meta$num==2,]$group=meta$sample_id[2]
data_meta[data_meta$num==3,]$group=meta$sample_id[3]

head(data_meta)
rownames(data_meta)=data_meta$bcnum

UCB_HSPC_iNK.seu <- CreateSeuratObject(counts = data.seu, project = "HSPC_iNK", min.cells = 5, min.features = 100,meta.data  = data_meta)
UCB_HSPC_iNK.seu

### global expression of scRNA-seq 
### The [[ operator can add columns to object metadata. This is a great place to stash QC stats
UCB_HSPC_iNK.seu[["percent.mt"]] <- PercentageFeatureSet(UCB_HSPC_iNK.seu, pattern = "^MT-")
UCB_HSPC_iNK.seu[["percent.ribo"]] <- PercentageFeatureSet(UCB_HSPC_iNK.seu, pattern = "^RP[SL]")
#rb.genes <- rownames(sce)[grep("^RP[SL]",rownames(sce))]
#percent.ribo <- Matrix::colSums(C[rb.genes,])/Matrix::colSums(C)*100
#UCB_HSPC_iNK.seu <- AddMetaData(UCB_HSPC_iNK.seu, percent.ribo, col.name = "percent.ribo")

UCB_HSPC_iNK.seu[["log10GenesPerUMI"]] <- log10(UCB_HSPC_iNK.seu$nFeature_RNA)/log10(UCB_HSPC_iNK.seu$nCount_RNA)
UCB_HSPC_iNK.seu[["UMIsPerGene"]] <- UCB_HSPC_iNK.seu$nCount_RNA/UCB_HSPC_iNK.seu$nFeature_RNA

### erythrocyte-related genes-hemoglobin family
### mouse: HB.genes <- "Hba-a1/Hba-a2/Hbb/Hbd/Hbb-y/Hbb-bh1/hbb-bh2/Hbb-b1/Hbb-b2/hbb-ar/hbb-bt/hbb-bs/lrp5/Hbq1a/Hbq1b/hba/Hba-x" %>% strsplit(.,split = "/") %>% unlist %>% capitalize(.)
HB.genes <- "HBA1/HBA2/HBB/HBD/HBE1/HBG1/HBG2/HBM/HBQ1/HBZ" %>% strsplit(.,split = "/") %>% unlist %>% toupper(.) %>% intersect(.,rownames(UCB_HSPC_iNK.seu))
if(!is.null(HB.genes)){
  percent.hb <- Matrix::colSums(UCB_HSPC_iNK.seu@assays$RNA@counts[HB.genes,])/Matrix::colSums(UCB_HSPC_iNK.seu@assays$RNA@counts)*100
  UCB_HSPC_iNK.seu <- AddMetaData(UCB_HSPC_iNK.seu, percent.hb, col.name = "percent.hb")
}

### Visualize QC metrics as a violin plot
VlnPlot(UCB_HSPC_iNK.seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo"), ncol = 4)

plot1 <- FeatureScatter(UCB_HSPC_iNK.seu, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(UCB_HSPC_iNK.seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2


p=list(VlnPlot(UCB_HSPC_iNK.seu, features = c( "nCount_RNA"),group.by = "group",pt.size = 0)+scale_y_continuous(limits = c(min(UCB_HSPC_iNK.seu$nCount_RNA),100000),breaks = seq(0,100000,10000)) +geom_hline(yintercept = c(20000,5000)),
       VlnPlot(UCB_HSPC_iNK.seu, features = c( "nFeature_RNA"),group.by = "group",pt.size = 0)+scale_y_continuous(limits = c(min(UCB_HSPC_iNK.seu$nFeature_RNA),max(UCB_HSPC_iNK.seu$nFeature_RNA)),breaks = seq(0,8000,500)) +geom_hline(yintercept = c(200,500,2000,1000,3000)),
       VlnPlot(UCB_HSPC_iNK.seu, features = c( "percent.mt"),group.by = "group",pt.size = 0)+scale_y_continuous(limits = c(min(UCB_HSPC_iNK.seu$percent.mt),40),breaks = seq(0,40,5)) +geom_hline(yintercept = c(5,10,20)),
       VlnPlot(UCB_HSPC_iNK.seu, features = c( "percent.ribo"),group.by = "group",pt.size = 0)+scale_y_continuous(limits = c(min(UCB_HSPC_iNK.seu$percent.ribo),max(UCB_HSPC_iNK.seu$percent.ribo)),breaks = seq(0,40,5)) +geom_hline(yintercept = c(5,10,20)))
p1=plot_grid(plotlist=p,ncol = 4)
p1

## Seurat object analysis workflow
### lognormalization  
### do log transformation: log(x*scale,factor/total counts of cell_index(x)+1)
UCB_HSPC_iNK.seu <- subset(UCB_HSPC_iNK.seu, subset = nFeature_RNA>1000 &nFeature_RNA<8000&nCount_RNA<60000  & percent.mt < 10) 
UCB_HSPC_iNK.seu <- NormalizeData(UCB_HSPC_iNK.seu, normalization.method = "LogNormalize", scale.factor = median(UCB_HSPC_iNK.seu$nCount_RNA))
UCB_HSPC_iNK.seu


### identify the highly various genes (HVGs)
### FindVariableFeatures() parameter: 
### nfeatures: the number of top highly expression variable gene needed
### selection.method: the method to get HVGs (vst (default)、mean.var.plot、dispersion)
UCB_HSPC_iNK.seu <- FindVariableFeatures(UCB_HSPC_iNK.seu, selection.method = "vst", nfeatures = 2000)

### scale data to do dimension reduction analysis
### scaled data saved in UCB_HSPC_iNK.seu[["RNA"]]@scale.data
all.genes <- rownames(UCB_HSPC_iNK.seu)
UCB_HSPC_iNK.seu <- ScaleData(UCB_HSPC_iNK.seu, features = all.genes)

saveRDS(UCB_HSPC_iNK.seu,file = "~/Documents/LJH_10X/results/data/UCB_NK_HSPC_iNK_aggr.rds")


### sctransform before IntegrateData
#UCB_NK=readRDS("./UCB_NK_merged.rds")
#CAR19_HSPC_iNK=readRDS(file = "~/Documents/LJH_10X/CAR19_HSPC_iNK/CAR19_UCB_HSPC_iNK.rds")
UCB_HSPC_iNK.seu=readRDS("~/Documents/LJH_10X/results/data/UCB_NK_HSPC_iNK_aggr.rds")

### run SCTransform on each object separately
UCB_HSPC_iNK=SCTransform(subset(UCB_HSPC_iNK.seu,subset=orig.ident=="iNK"), verbose = FALSE)
UCB_aNK=SCTransform(subset(UCB_HSPC_iNK.seu,subset=orig.ident=="UCB-NK"), verbose = FALSE)
CAR19_HSPC_iNK=SCTransform(subset(UCB_HSPC_iNK.seu,subset=orig.ident=="CD19-CAR-iNK"), verbose = FALSE)

### set the maximum allowed size of seurat object list
options(future.globals.maxSize = 5000 * 1024^2)
### select features for downstream integration, and run PrepSCTIntegration, which ensures that all necessary Pearson residuals have been calculated
pb.features <- SelectIntegrationFeatures(object.list = list(CAR19_HSPC_iNK,UCB_HSPC_iNK,UCB_aNK), nfeatures = 2000)
pb.list <- PrepSCTIntegration(object.list = list(CAR19_HSPC_iNK,UCB_HSPC_iNK,UCB_aNK), anchor.features = pb.features,verbose = T)

### identify anchors and integrate the datasets. Commands are identical to the standard workflow, but make sure to set normalization.method = 'SCT'
pb.anchors <- FindIntegrationAnchors(object.list = pb.list, normalization.method = "SCT",anchor.features = pb.features, verbose = FALSE,k.filter = 100)
pb.integrated <- IntegrateData(anchorset = pb.anchors, normalization.method = "SCT", verbose = FALSE)


### dimplot
### visualization, clustering
DefaultAssay(pb.integrated) <- "integrated"
pb.integrated <- RunPCA(pb.integrated, verbose = FALSE)
ElbowPlot(pb.integrated,ndims = 30)

pb.integrated <- RunUMAP(pb.integrated, dims = 1:20,n.neighbors = 30)
plots <- DimPlot(pb.integrated, split.by = c("orig.ident"),label = T)
plots & theme(legend.position = "top") & guides(color = guide_legend(nrow = 3, byrow = TRUE, 
                                                                     override.aes = list(size = 3)))

### UMAP plot  
pb.integrated$group=factor(pb.integrated$orig.ident,levels = c("UCB-NK","iNK","CD19-CAR-iNK"))

p1 <- DimPlot(pb.integrated, reduction = "umap", label = F,group.by = "group",pt.size = 0.0001)+labs(title = "")&theme_bw() &
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
    plot.title = element_text(face="bold.italic",hjust = 0.5,size=10))& NoLegend()

#ggsave(p1, filename = "./umap_projection_0706.pdf",width = 13,height = 3)

## save data
saveRDS(pb.integrated,file = "~/Documents/LJH_10X/CAR19_HSPC_iNK/UCB_aNK_HSPC_CAR19_iNK.rds")

