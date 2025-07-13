
# ------------------------------------------------------------------------------
# title: "10x Single nuclei RNA (SnRNA), identifying key modules that responsible for intrinsic and extrinsic mechanisms of tumor development in tumor micro-environment (TME), from 21 mouse scRNA libraries."
# Author: Nasim Rahmatpour
# Date: "6/3/2022"
# ------------------------------------------------------------------------------


# Project Overview

#This R Markdown document summarizes the analysis of 21 single-cell RNA-seq samples across 3 tumor cell lines, each treated with SOS or Veh. The goal is to understand treatment effects on cell composition, tumor pathway activity, and tumor-specific gene expression patterns.
#The pipeline includes QC, annotation,  DEG analysis, and pathway scoring.

---
  
###Task 1: Initial Processing and Broad Characterization 
library(Seurat)
suppressMessages(require(cowplot))
suppressMessages(require(dplyr))
suppressMessages(library(openxlsx))
suppressMessages(library(ggplot2))
library(ggplot2)
library(RColorBrewer)
library(angrycell)


ids <- c("scGEX_152", "scGEX_153", "scGEX_154", "scGEX_155", "scGEX_156", "scGEX_157", "scGEX_158", "scGEX_159", "scGEX_160", "scGEX_161", "scGEX_162", "scGEX_163", "scGEX_164", "scGEX_165", "scGEX_166", "scGEX_167", "scGEX_168", "scGEX_169", "scGEX_170", "scGEX_171", "scGEX_172")

#make Seurat object for each
for (file in ids){
  BI <- Read10X(data.dir = paste0("./", file))
  BIdata <- CreateSeuratObject(counts = BI, project = file)
  BIdata[["percent.ribo"]] <- PercentageFeatureSet(object=BIdata, pattern= '^Rp[sl]')
  BIdata[["percent.mito"]] <- PercentageFeatureSet(object=BIdata, pattern = '^mt-')
  BIdata[["percent.Hb"]] <- PercentageFeatureSet(object=BIdata, pattern ='^Hb')
  assign(file, BIdata)
}

#make output directory
dir.create(paste0("./", "BI_out"))

#merging the files
BI_combined <- merge(`scGEX_152`, c(`scGEX_153`, `scGEX_154`, `scGEX_155`, `scGEX_156`, `scGEX_157`, `scGEX_158`, `scGEX_159`, `scGEX_160`, `scGEX_161`, `scGEX_162`, `scGEX_163`, `scGEX_164`, `scGEX_165`, `scGEX_166`, `scGEX_167`, `scGEX_168`, `scGEX_169`, `scGEX_170`, `scGEX_171`, `scGEX_172`), 
                     add.cell.ids = c("152", "153","154","155", "156", "157", "158", "159", "160", "161", "162", "163", "164", "165", "166", "167", "168", "169", "170", "171", "172"), project = "BI")

#dimension of combined data
dim(BI_combined)

#add Sample column
BI_combined$Sample <- BI_combined$orig.ident

#function for median values
tf <- function(x){tapply(x, Idents(BI_combined), median)}
pdf(paste0("./BI_out/","combined_violin.pdf"),width = 20,height = 10)
Idents(BI_combined) <- BI_combined$Sample
median.value <- apply(BI_combined@meta.data[,c("nCount_RNA", "nFeature_RNA", "percent.ribo", "percent.mito", "percent.Hb")],2,tf)
write.csv(median.value, paste0("./BI_out/",'Combined_median.csv'))
VlnPlot(BI_combined, features = c("nCount_RNA", "nFeature_RNA", "percent.ribo", "percent.mito", "percent.Hb"),group.by='Sample', pt.size = 0,ncol = 3)
dev.off()

#standard analysis on combined data
BI_combined <- NormalizeData(BI_combined, normalization.method = "LogNormalize", scale.factor = 10000)
BI_combined <- FindVariableFeatures(BI_combined, verbose = FALSE, selection.method = "vst", nfeatures = 3000)
BI_combined <- ScaleData(BI_combined)
BI_combined <- RunPCA(BI_combined, features = VariableFeatures(object = BI_combined))
BI_combined <- FindNeighbors(BI_combined,reduction = "pca",  dims=1:30)
BI_combined <- FindClusters(BI_combined)
BI_combined <- RunUMAP(BI_combined,reduction = "pca", dims=1:30)
BI_combined<- RunTSNE(BI_combined,reduction = "pca", dims=1:30)

#umap
pdf(paste0("./BI_out/","combined_umap.pdf"),width = 20,height = 10)
p1 <- DimPlot(BI_combined, reduction = "umap", group.by = "Sample", label = TRUE, repel = TRUE)
p1
dev.off()

#save
saveRDS(BI_combined, paste0("./BI_out/",'BI_combined.rds'))

#filter the low quality cells
BI_filter <- subset(BI_combined, subset=nCount_RNA > 200 & percent.mito < 20)
dim(BI_filter)
BI_filter$Sample <- BI_filter$orig.ident

tf <- function(x){tapply(x, Idents(BI_filter), median)}
pdf(paste0("./BI_out/","filter_violin.pdf"),width = 20,height = 10)
Idents(BI_filter) <- BI_filter$Sample
median.value <- apply(BI_filter@meta.data[,c("nCount_RNA", "nFeature_RNA", "percent.ribo", "percent.mito", "percent.Hb")],2,tf)
write.csv(median.value, paste0("./BI_out/",'Filter_median.csv'))
VlnPlot(BI_filter, features = c("nCount_RNA", "nFeature_RNA", "percent.ribo", "percent.mito", "percent.Hb"),group.by='Sample', pt.size = 0,ncol = 3)
dev.off()

#standard analysis on filtered data
BI_filter <- NormalizeData(BI_filter, normalization.method = "LogNormalize", scale.factor = 10000)
BI_filter <- FindVariableFeatures(BI_filter, verbose = FALSE, selection.method = "vst", nfeatures = 3000)
BI_filter <- ScaleData(BI_filter)
BI_filter <- RunPCA(BI_filter, features = VariableFeatures(object = BI_filter))
BI_filter <- FindNeighbors(BI_filter,reduction = "pca",  dims=1:30)
BI_filter <- FindClusters(BI_filter, resolution = 0.4)
BI_filter <- RunUMAP(BI_filter,reduction = "pca", dims=1:30)


#umap
pdf(paste0("./BI_out/","filtered_umap.pdf"),width = 20,height = 10)
p1 <- DimPlot(BI_filter, reduction = "umap", group.by = "Sample")
p2 <- DimPlot(BI_filter, reduction = "umap", group.by = "seurat_clusters", label= TRUE, size =20)
p1
p2
dev.off()

#gene markers
markers_RNA <- FindAllMarkers(BI_filter, only.pos = TRUE)
write.xlsx(markers_RNA,file = paste0("./BI_out/clusters_markers.xlsx"),row.names = TRUE)

#top 20
top20.markers <- markers_RNA %>% group_by(cluster) %>% top_n(n=20, wt=avg_log2FC)
top20 <- c()
for(c in unique(top20.markers$cluster)){
  tmp <- top20.markers[top20.markers$cluster==c,]$gene
  top20 <- cbind(top20,tmp)
}
colnames(top20) <- paste0('C_',unique(top20.markers$cluster))
write.csv(top20,'./BI_out/clusters_top20.markers.csv')

#save
saveRDS(BI_filter, paste0("./BI_out/",'BI_filter.rds'))

#make some dot and feature plots of gene markes for broad annotation
pdf(paste0('./BI_out/canonical_markers.pdf'))
library(ggplot2)
library(RColorBrewer)
library(angrycell)
Idents(BI_filter) <- BI_filter@meta.data$seurat_clusters
angrycell::DotPlot2(BI_filter, features = c( "Cd8a", "Cd8b1","Cd4", "Ccr7","Il7r", "Tcf7","Foxp3", "Il2ra", "Klra3", "Klra7", "Klrb1", "Hbb-bs", "Hba-a1", "Jchain","Cd19", "Cd79a","Ms4a1", "Blk","Slamf7", "Mzb1", "Irf7", "S100a9", "S100a8", "Csf3r", "Ms4a2", "Csf1r", "C3ar1", "Cd163", "Ms4a4a", "Kit", "Cpa3","Mki67", "Hist1h2ap", "Lilra4", "Pecam1", "Fdgfrb", "Dcn", "Col1a1", "Top2a", "Ube2c", "Krt8", "Krt18", "C1qa", "C1qb", "H2-Aa", "H2-Ab1", "Trem2", "Cd68", "Ly6c2", "Cxcl8", "Il6", "Cxcr1", "Cxcr2", "Mmp9", "Msln", "Acta2", "Srpx2", "Fgf7"), dot.scale = 8) +
  coord_flip() +
  theme(text = element_text(size=12),
        axis.text.x=element_text(size=12),
        axis.text.y=element_text(size=11),
        axis.title.y = element_text(size =11),
        axis.title.x = element_text(size=12)) 
dev.off()


#feature plot1
pdf(paste0('./BI_out/cannical_feature_plot1.pdf'))
FeaturePlot(BI_filter, features = c("Cd8a","Cd4", "Cd3d","Il7r", "Jchain","Cd79a","Ms4a1", "Hbb-bs", "Hba-a1", "Kit", "S100a9", "S100a8","Csf3r", "Lilra4", "Cd68", "Ly6c2")) &theme(text = element_text(size = 8),
                                                                                                                                                                                     axis.text.x=element_text(size=5),
                                                                                                                                                                                     axis.text.y=element_text(size=5),
                                                                                                                                                                                     axis.title.y = element_text(size =5),
                                                                                                                                                                                     axis.title.x = element_text(size=5))
dev.off() 

#feature plot2
pdf(paste0('./BI_out/cannical_feature_plot2.pdf'))
FeaturePlot(BI.filter_YFP_AddModuleScore, features = c("Dcn", "Col1a1", "Krt8", "Krt18", "Pecam1", "YFP", "Mki67")) &theme(text = element_text(size = 8),
                                                                                                                           axis.text.x=element_text(size=5),
                                                                                                                           axis.text.y=element_text(size=5),
                                                                                                                           axis.title.y = element_text(size =5),
                                                                                                                           axis.title.x = element_text(size=5))
dev.off() 

#YFP feature plot
pdf(paste0('./BI_out/YFP_featureplot.pdf'))
FeaturePlot(BI.filter_YFP_AddModuleScore, features = "YFP") & theme(text = element_text(size = 8),
                                                                    axis.text.x=element_text(size=5),
                                                                    axis.text.y=element_text(size=5),
                                                                    axis.title.y = element_text(size =5),
                                                                    axis.title.x = element_text(size=5))
dev.off()


#broad annotation
Idents(BI_filter) <- BI_filter@meta.data$seurat_clusters
BI_filter<- RenameIdents(BI_filter, `0` = "Epithelial", `1` = "Neutrophil", `2` = "Epithelial",
                         `3` = "Myeloid", `4` = "Epithelial", `5` = "Epithelial", `6` = "T", 
                         `7` = "Fibroblast",`8` = "Proliferating Epithelial", `9` = "Endothelial",
                         `10` = "RBC", `11` = "pDC", `12` = "Fibroblast", `13`= "Myeloid", `14`= "Fibroblast", `15`= "B", `16`="Epithelial", `17`= "Myeloid", `18` ="Fibroblast", `19`= "Mast", `20`= "RBC+Myeloid")
BI_filter@meta.data$MajorCellType <- Idents(BI_filter)

pdf(paste0('./BI_out/annotation_plot.pdf'))
DimPlot(BI_filter, reduction = "umap", label = TRUE, repel = TRUE)
dev.off()


# plots to compare the major cell types by Sample, Cell line and Treatments
#per Sample
pdf(paste0('./BI_out/CellMajorType_per_Sample.pdf'))
library(angrycell)
BI_filter@meta.data$MajorCellType <- Idents(BI_filter)
meta <- BI_filter@meta.data
angrycell::plot_fraction(meta, x = 'orig.ident', fill.bar = 'MajorCellType')+
  # theme(axis.ticks.y=element_blank())
  theme(#axis.text.x=element_blank(),
    #axis.ticks.x=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank())
#xlab("")
dev.off()

#by CellLine
pdf(paste0('./BI_out/CellMajorType_per_CellLine.pdf'))
library(angrycell)
BI_filter@meta.data$MajorCellType <- Idents(BI_filter)
meta <- BI_filter@meta.data
angrycell::plot_fraction(meta, x = 'CellLine', fill.bar = 'MajorCellType')
dev.off()

#by Treatment
pdf(paste0('./BI_out/CellMajorType_per_Treatment.pdf'))
library(angrycell)
BI_filter@meta.data$MajorCellType <- Idents(BI_filter)
meta <- BI_filter@meta.data
angrycell::plot_fraction(meta, x = 'Treatment', fill.bar = 'MajorCellType')
dev.off()

#calculate the percentage and the raw count of Major cell types for each sample, splits the umap for each Sample
meta.df <- table(BI_filter@meta.data$MajorCellType,BI_filter@meta.data$orig.ident)
meta.df
meta.totalper.sample <- apply(meta.df,2,sum)
meta.totalper.sample
meta.ratio <- sweep(meta.df, 2, meta.totalper.sample, `/`)
meta.ratio
meta.pct <- meta.ratio*100
meta.pct 
meta.pctsig <- round(meta.pct, digits = 2)
meta.pctsig
write.xlsx(meta.pctsig,file = paste0("./BI_out/MajorCellType_per_Sample.xlsx"),row.names = TRUE)
write.csv(meta.pctsig,file = paste0("./BI_out/MajorCellType_per_Sample.csv"),row.names = TRUE)

#save table of total counts
meta.df <- table(BI_filter@meta.data$MajorCellType,BI_filter@meta.data$orig.ident)
meta.df
write.csv(meta.df,file = paste0("./BI_out/Total_MajorCellType_per_Sample.csv"),row.names = TRUE)


#umap split for each sample
pdf("./BI_out/split_umap.pdf")
DimPlot(BI_filter, reduction = "umap", split.by = "Sample", group.by = "MajorCellType", ncol = 7, size =5)
dev.off()

###Task 2: Identifying Tumor vs. Normal Cells Across All Cell Types
#this part identifies the Tumor vs. Normal Cells based on YFP gene expression. The fastq reads were mapped to YFP sequence by Kallisto. This parts adds the YFP counts of each cell from each sample, adds to the metadata of seurat.

#read the Kallisto outputs by this function
theme_set(theme_bw())
# Slightly modified from BUSpaRse, just to avoid installing a few dependencies not used here
read_count_output <- function(dir, name='output') {
  dir <- normalizePath(dir, mustWork = TRUE)
  m <- readMM(paste0(dir, "/", name, ".mtx"))
  m <- Matrix::t(m)
  m <- as(m, "dgCMatrix")
  # The matrix read has cells in rows
  ge <- ".genes.txt"
  genes <- readLines(file(paste0(dir, "/", name, ge)))
  barcodes <- readLines(file(paste0(dir, "/", name, ".barcodes.txt")))
  colnames(m) <- barcodes
  rownames(m) <- genes
  return(m)
}

#run the function on all samples and combine them
setwd("/home/nasim/BI")
folders = c('scGEX_152', 'scGEX_153', 'scGEX_154', 'scGEX_155', 'scGEX_156', 'scGEX_157', 'scGEX_158', 'scGEX_159', 'scGEX_160', 'scGEX_161', 'scGEX_162', 'scGEX_163', 'scGEX_164', 'scGEX_165', 'scGEX_166', 'scGEX_167', 'scGEX_168', 'scGEX_169', 'scGEX_170', 'scGEX_171', 'scGEX_172')
tmp_yfp_count = data.frame()
for (dir in folders){
  print(dir)
  a = read_count_output(dir)
  b = as.data.frame(a['YFP',])
  colnames(b) = c('YFP_count')
  b$cell.id= rownames(b)
  b$lib= dir
  tmp_yfp_count = rbind(tmp_yfp_count,b)
}

#add the YFP counts of each cell to meta data of seurat
dim(tmp_yfp_count)
tmp_yfp_count$sample = apply(tmp_yfp_count,1,function(x) substring(x[3],7,9))
tmp_yfp_count$cell.id = paste0(tmp_yfp_count$sample,'_',tmp_yfp_count$cell.id,"-1")
table(tmp_yfp_count$cell.id %in% colnames(BI_filter))

coords=data.frame(BI_filter@reductions$umap@cell.embeddings); dim(coords)
meta_all <- BI_filter@meta.data
meta_all$cell.id <- rownames(meta_all)
meta_all <- cbind(meta_all,coords)
new_meta_all <- merge(meta_all, tmp_yfp_count, by ="cell.id")
#replace the index numbers as row with cell.id which is the first column
rownames(new_meta_all) <- new_meta_all[,1]

#add a column to meta data based on the expression of YFP
new_meta_all$gene.pos.YFP[new_meta_all[,'YFP_count'] > 0] <- 1
new_meta_all$gene.pos.YFP[new_meta_all[,'YFP_count'] ==0] <- 0
new_meta_all$YFP.expression[new_meta_all[,'YFP_count'] > 0] <- "Positive"
new_meta_all$YFP.expression[new_meta_all[,'YFP_count'] == 0] <- "Negative"

#make a new dataframe based on YFP dataframe
YFP.sum = new_meta_all %>%
  group_by(MajorCellType, orig.ident) %>%
  summarize( YFP.positive = round( 100 *  sum(gene.pos.YFP) / n(), digits = 2), 
             YFP.counts = sum(YFP_count)
  )

summary(YFP.sum)

#make some visualizations 
pdf(paste0("./BI_out/YFP_percentage_count1.pdf"), width=20, height=10)
library(RColorBrewer)
library(ggplot2)
library(ggthemes)
ggplot(aes ( x = orig.ident, y = YFP.counts, size = YFP.positive), data = YFP.sum) +
  geom_point(color = 'dark red', aes(alpha = YFP.positive )) +  coord_flip() +
  facet_wrap(scales = "free_x", "MajorCellType" ) + theme_hc()+
  theme(axis.text.x = element_text(face = "bold", size=7),
        axis.text.y = element_text(face = "bold", color = "black", size = 7))
dev.off() 

pdf(paste0("./BI_out/YFP_percentage_count2.pdf"), width=20, height=10)
library(RColorBrewer)
library(ggplot2)
library(ggthemes)
ggplot(aes ( x = MajorCellType, y = YFP.counts, size = YFP.positive), data = YFP.sum) +
  geom_point(color = 'dark red', aes(alpha = YFP.positive )) +  coord_flip() + facet_wrap(scales = "free_x", "orig.ident" ) + theme_hc() +
  theme(axis.text.x = element_text(face = "bold", size=7),
        axis.text.y = element_text(face = "bold", color = "black", size = 7))
dev.off()

pdf(paste0("./BI_out/YFP_percentage_count3.pdf"), width=20, height=10)
library(RColorBrewer)
library(ggplot2)
library(ggthemes)
ggplot(aes ( x = MajorCellType, y = YFP.counts, size = YFP.positive), data = YFP.sum) +
  geom_point(color = 'dark red', aes(alpha = YFP.positive )) + facet_wrap(scales = "free_y", "orig.ident" ) + theme_hc()+
  theme(axis.text.x = element_text(angle= 45, hjust= 1, face="bold", size=7),
        axis.text.y = element_text(face = "bold", color = "black", size = 7))
dev.off() 


YFP.sum$YFP.negative <- ""
YFP.sum$YFP.negative <- apply(YFP.sum, 1, function(x) 100-YFP.sum$YFP.positive)

#sample comparison
pdf(paste0("./BI_out/YFP_count1.pdf"), width=20, height=10)
library(RColorBrewer)
library(ggplot2)
library(ggthemes)
ggplot(new_meta_all, aes ( x = orig.ident, fill = YFP.expression))+
  geom_bar(stat= "count", position = "dodge") +
  ylab("cell counts")+
  #scale_fill_brewer()+
  scale_fill_manual(values=c("#006ddb", "#FF0000")) +
  theme_hc(base_size = 22)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face="bold", size=15),
        axis.text.y = element_text(face = "bold", color = "black", size = 15))
dev.off() 

#each celltype
pdf(paste0("./BI_out/YFP_count2.pdf"), width=20, height=10)
library(RColorBrewer)
library(ggplot2)
library(ggthemes)
ggplot(new_meta_all, aes ( x = orig.ident, fill = YFP.expression))+
  geom_bar(stat= "count", position = "dodge") +
  ylab("cell counts")+
  #scale_fill_brewer()+
  scale_fill_manual(values=c("#006ddb", "#FF0000")) +
  theme_hc(base_size = 15)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face="bold", size=10),
        axis.text.y = element_text(face = "bold", color = "black", size = 10))+
  facet_wrap("MajorCellType")
dev.off() 

#plot on meta_all
library(ggplot2)
library(viridis)
ggplot(aes_string (x = 'UMAP_1', y = 'UMAP_2', color= 'YFP_count', alpha= 'YFP_count'), data= new_meta_all) +
  geom_point(size = 1) +
  scale_colour_viridis(option = 'rocket')+
  facet_wrap("orig.ident")
labs(title = ) +
  theme(legend.position = 'none')

#make new seurat, violin and feature on meta-all
RNA_data=BI_filter@assays$RNA@data; dim(RNA_data)
BI_filter_meta_all <- CreateSeuratObject(counts =RNA_data[, rownames(new_meta_all)], 
                                         project = "YFP", 
                                         assay = "RNA", 
                                         meta.data = new_meta_all)
VlnPlot(BI_filter_meta_all, features = c("YFP_count"), group.by = "orig.ident") + NoLegend()
BI_filter_meta_all@reductions <- BI_filter@reductions
FeaturePlot(BI_filter_meta_all, features = "YFP_count", cols = c("grey", " dark red")) #, split.by= "orig.ident", ncol=7)

BI_filter@meta.data <- cbind(BI_filter@meta.data, new_meta_all)
BI_filter@meta.data <- BI_filter@meta.data[,(14:32)]
write.csv(tmp_yfp_count, "./BI_out/tmp_yfp_count.csv")
write.csv(YFP.sum, "./BI_out/YFP.sum.csv")
saveRDS(BI_filter, "./BI_out/BI.filter_YFP.rds")



#Violin plot of YFP expression grouping by Sample and MajorCellType
pdf("./BI_out/YFP_Violin_UMAP.pdf")
VlnPlot(BI.filter_YFP, features = c("YFP_count"), group.by = "orig.ident") + NoLegend()
VlnPlot(BI.filter_YFP, features = c("YFP_count"), group.by = "MajorCellType") + NoLegend()
FeaturePlot(BI.filter_YFP, features = "YFP_count", cols = c("grey", " dark red"))
dev.off()




#this part runs InfeCNV analysis to identify the epithelial tumor cells, setting B & T cells as the reference
#save the matrix
counts.matrix=BI_filter_YFP@assays$RNA@data; dim(counts.matrix)

cellcounts <- GetAssayData(BI.filter_YFP[['RNA']],slot='counts')
saveRDS(cellcounts,'./BI_out/cellcounts.rds')

#downsampling, split cell indices by identity
idx <- split(meta, f=meta$MajorCellType)
#idx1 <- unsplit(idx, meta$MajorCellType)
cs_keep <- lapply(idx, function(i) {
  print(class(i))
  i=i[1:100,]
})  

meta_downsample=data.frame() 
for (i in cs_keep){
  meta_downsample=rbind(meta_downsample,i)
}

#downsampling, another way
meta_downsample1 <- data.frame()
for(c in unique(meta$MajorCellType)){
  print(meta[meta$MajorCellType==c,][1:100,])
  meta_downsample1 <- rbind(meta_downsample1,meta[meta$MajorCellType==c,][1:100,])
}  
write.csv(meta_downsample,file=paste0("./BI_out/meta_downsample.csv"))
write.csv(meta_downsample1,file=paste0("./BI_out/meta_downsample1.csv"))
#meta_downsample <- read.csv(file = "./BI_out/meta_downsample.csv", row.names = 1)


#make seurat object for smaller meta
BI_filter_YFP_downsample <- CreateSeuratObject(counts = counts.matrix[, rownames(meta_downsample)], 
                                               project = "BI_downsample", 
                                               assay = "RNA", 
                                               meta.data = meta_downsample)

#get the smaller matrix
counts.matrix.small=BI_filter_YFP_downsample@assays$RNA@data; dim(counts.matrix.small)

#save
saveRDS(counts.matrix.small, "./BI_out/counts.matrix.small.rds")

#create infercnv object
infercnv_obj = CreateInfercnvObject(raw_counts_matrix="/home/nasim/BI/BI_out/counts.matrix.small.rds",
                                    annotations_file="/home/nasim/BI/BI_out/annotation_down.txt",
                                    delim="\t",
                                    gene_order_file="/home/nasim/BI/BI_out/gene_pos.txt",
                                    ref_group_names=c("T", "B"))

#run infercnv
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir="infercnv_output_dir",  # dir is auto-created for storing outputs
                             cluster_by_groups=T,   # cluster
                             denoise=T,
                             HMM=T,
)

###Task 3: Differential Expression and Enrichment Analysis in Epithelial/Tumor Cells
#This part is for DGE anlysis based on single cell in each cell line between 2 treatments, making volcano plot, runnign GO enrichment analysis and finding the overlap of up/down regulated gene among 3 cell lines

Rcpp::sourceCpp("/home/chris/project/fast_de/fast_de.cpp")
TC3 <- subset(Epithelial, subset= CellLine=="TC3")
idents <- TC3@meta.data$Treatment
markers_rcpp_TC3 <- t_test_1vA(TC3@assays$RNA@data, idents)

Rcpp::sourceCpp("/home/chris/project/fast_de/fast_de.cpp")
TC4 <- subset(Epithelial, subset= CellLine=="TC4")
idents <- TC4@meta.data$Treatment
markers_rcpp_TC4 <- t_test_1vA(TC4@assays$RNA@data, idents)

Rcpp::sourceCpp("/home/chris/project/fast_de/fast_de.cpp")
TC5 <- subset(Epithelial, subset= CellLine=="TC5")
idents <- TC5@meta.data$Treatment
markers_rcpp_TC5 <- t_test_1vA(TC5@assays$RNA@data, idents)

#subset
markers_rcpp_TC3_subset <- markers_rcpp_TC3[markers_rcpp_TC3$cluster=="SOSMEK",]
markers_rcpp_TC4_subset <- markers_rcpp_TC4[markers_rcpp_TC4$cluster=="SOSMEK",]
markers_rcpp_TC5_subset <- markers_rcpp_TC5[markers_rcpp_TC5$cluster=="SOSMEK",]

#TC3 volcano plot
pdf(paste0("./BI_out/", "Epithelial_TC3_DEGs_volcano.pdf"))
markers_rcpp_TC3_subset$expression <- "not significant"
markers_rcpp_TC3_subset$expression[markers_rcpp_TC3_subset$p_val < 0.05 & markers_rcpp_TC3_subset$avg_log2FC > 0.5] <- "upregulated"
markers_rcpp_TC3_subset$expression[markers_rcpp_TC3_subset$p_val < 0.05 & markers_rcpp_TC3_subset$avg_log2FC < -0.5] <- "downregulated"
top_genes <- markers_rcpp_TC3_subset %>% 
  filter(expression != 'not significant')

top_genes

ggplot(data=markers_rcpp_TC3_subset, aes(x=avg_log2FC, y=-log10(p_val), col=expression)) + 
  geom_point(size=0.5) + 
  #theme_minimal() +
  theme_bw() +
  theme(text=element_text(size=11))+
  #geom_text_repel(data=head(markers_rcpp_TC3_subset, 50), aes(label=gene),nudge_x = .001 ,nudge_y = 0.001, size=3) +
  #geom_text_repel(data=markers_rcpp_TC3_subset, aes(label=gene),nudge_x = 0 ,nudge_y = 0, size=3) +
  #geom_text_repel(data=top_genes, aes(label=gene)) +
  geom_label_repel(data = top_genes, 
                   mapping = aes(avg_log2FC, -log10(p_val), label = gene), size = 2, max.overlaps=Inf)+
  
  theme(axis.text = element_text(face="bold"))+
  theme(axis.line.x = element_line(color="black", size = 0.3),
        axis.line.y = element_line(color="black", size = 0.3))+
  #geom _text is about the size of dots
  #geom_text(position=position_dodge(width=0.9),hjust=1.5, size= 0.1)+
  #scale_color_manual(values=c("red", "blue")) +
  scale_color_manual(values=c("blue", "dark grey", "red"))+
  geom_vline(xintercept=c(-0.5, 0.5), col="black", linetype="dashed") +
  geom_hline(yintercept=-log10(0.01), col="black", linetype="dashed") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 
dev.off()

#TC4 volcano plot
pdf(paste0("./BI_out/", "Epithelial_TC4_DEGs_volcano.pdf"))
markers_rcpp_TC4_subset$expression <- "not significant"
markers_rcpp_TC4_subset$expression[markers_rcpp_TC4_subset$p_val < 0.05 & markers_rcpp_TC4_subset$avg_log2FC > 0.5] <- "upregulated"
markers_rcpp_TC4_subset$expression[markers_rcpp_TC4_subset$p_val < 0.05 & markers_rcpp_TC4_subset$avg_log2FC < -0.5] <- "downregulated"
#top <- 20
#top_genes <- bind_rows(
#markers_rcpp_TC4_subset %>% 
#filter(expression == 'upregulated') %>% 
#arrange(p_val, desc(abs(avg_log2FC))) %>% 
#head(top),
#markers_rcpp_TC4_subset %>% 
#filter(expression == 'downregulated') %>% 
#arrange(p_val, desc(abs(avg_log2FC))) %>% 
#head(top)
#)
#top_genes
top_genes <- markers_rcpp_TC4_subset %>% 
  filter(expression != 'not significant')

top_genes
ggplot(data=markers_rcpp_TC4_subset, aes(x=avg_log2FC, y=-log10(p_val), col=expression)) + 
  geom_point(size=0.5) + 
  #theme_minimal() +
  theme_bw() +
  theme(text=element_text(size=11))+
  #geom_text_repel(data=head(markers_rcpp_TC3_subset, 50), aes(label=gene),nudge_x = .001 ,nudge_y = 0.001, size=3) +
  #geom_text_repel(data=markers_rcpp_TC3_subset, aes(label=gene),nudge_x = 0 ,nudge_y = 0, size=3) +
  #geom_text_repel(data=top_genes, aes(label=gene)) +
  geom_label_repel(data = top_genes, 
                   mapping = aes(avg_log2FC, -log10(p_val), label = gene), size = 2, max.overlaps=Inf)+
  
  theme(axis.text = element_text(face="bold"))+
  theme(axis.line.x = element_line(color="black", size = 0.3),
        axis.line.y = element_line(color="black", size = 0.3))+
  #geom _text is about the size of dots
  #geom_text(position=position_dodge(width=0.9),hjust=1.5, size= 0.1)+
  #scale_color_manual(values=c("red", "blue")) +
  scale_color_manual(values=c("blue", "dark grey", "red"))+
  geom_vline(xintercept=c(-0.5, 0.5), col="black", linetype="dashed") +
  geom_hline(yintercept=-log10(0.1), col="black", linetype="dashed") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 
dev.off()

#TC5 volcano plot
pdf(paste0("./BI_out/", "Epithelial_TC5_DEGs_volcano.pdf"))
markers_rcpp_TC5_subset$expression <- "not significant"
markers_rcpp_TC5_subset$expression[markers_rcpp_TC5_subset$p_val < 0.05 & markers_rcpp_TC5_subset$avg_log2FC > 0.5] <- "upregulated"
markers_rcpp_TC5_subset$expression[markers_rcpp_TC5_subset$p_val < 0.05 &  markers_rcpp_TC5_subset$avg_log2FC < -0.5] <- "downregulated"
#top <- 20
#top_genes <- bind_rows(
#markers_rcpp_TC5_subset %>% 
#filter(expression == 'upregulated') %>% 
#arrange(p_val, desc(abs(avg_log2FC))) %>% 
#head(top),
#markers_rcpp_TC5_subset %>% 
#filter(expression == 'downregulated') %>% 
#arrange(p_val, desc(abs(avg_log2FC))) %>% 
#head(top)
#)
#top_genes
top_genes <- markers_rcpp_TC5_subset %>% 
  filter(expression != 'not significant')
ggplot(data=markers_rcpp_TC5_subset, aes(x=avg_log2FC, y=-log10(p_val), col=expression)) + 
  geom_point(size=0.5) + 
  #theme_minimal() +
  theme_bw() +
  theme(text=element_text(size=11))+
  #geom_text_repel(data=head(markers_rcpp_TC3_subset, 50), aes(label=gene),nudge_x = .001 ,nudge_y = 0.001, size=3) +
  #geom_text_repel(data=markers_rcpp_TC3_subset, aes(label=gene),nudge_x = 0 ,nudge_y = 0, size=3) +
  #geom_text_repel(data=top_genes, aes(label=gene)) +
  geom_label_repel(data = top_genes, 
                   mapping = aes(avg_log2FC, -log10(p_val), label = gene), size = 2, max.overlaps=Inf)+
  
  theme(axis.text = element_text(face="bold"))+
  theme(axis.line.x = element_line(color="black", size = 0.3),
        axis.line.y = element_line(color="black", size = 0.3))+
  #geom _text is about the size of dots
  #geom_text(position=position_dodge(width=0.9),hjust=1.5, size= 0.1)+
  #scale_color_manual(values=c("red", "blue")) +
  scale_color_manual(values=c("blue", "dark grey", "red"))+
  geom_vline(xintercept=c(-0.5, 0.5), col="black", linetype="dashed") +
  geom_hline(yintercept=-log10(0.1), col="black", linetype="dashed") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 
dev.off()

TC3 <- read.csv("./BI_out/Epithelial_TC3_SOSMEK_significant.csv", row.names = 1)
TC3_up <- subset(TC3, expression=="upregulated")
TC3_down <- subset(TC3, expression=="downregulated")
TC4 <- read.csv("./BI_out/Epithelial_TC4_SOSMEK_significant.csv", row.names = 1)
TC4_up <- subset(TC4, expression=="upregulated")
TC4_down <- subset(TC4, expression=="downregulated")
TC5 <- read.csv("./BI_out/Epithelial_TC5_SOSMEK_significant.csv", row.names = 1)
TC5_up <- subset(TC5, expression=="upregulated")
TC5_down <- subset(TC5, expression=="downregulated")


#these libraries needed for GO
library(clusterProfiler)
library(org.Mm.eg.db)
library(org.Hs.eg.db)


#get the GO of all categories together
markerlist_TC3_up <- as.character(TC3_up$gene)
eg1 <- bitr(markerlist_TC3_up , fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb = org.Mm.eg.db)
genelist_TC3_up <- eg1$SYMBOL
Markers_GO_enrich_TC3_up <- enrichGO(genelist_TC3_up, OrgDb = org.Mm.eg.db, keyType="SYMBOL", pAdjustMethod = 'BH', ont= "ALL", pvalueCutoff = 0.5)


markerlist_TC4_up <- as.character(TC4_up$gene)
eg2 <- bitr(markerlist_TC4_up , fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb = org.Mm.eg.db)
genelist_TC4_up <- eg2$SYMBOL
Markers_GO_enrich_TC4_up <- enrichGO(genelist_TC4_up, OrgDb = org.Mm.eg.db, keyType="SYMBOL", pAdjustMethod = 'BH', ont= "ALL", pvalueCutoff = 0.5)


markerlist_TC5_up <- as.character(TC5_up$gene)
eg3 <- bitr(markerlist_TC5_up , fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb = org.Mm.eg.db)
genelist_TC5_up <- eg3$SYMBOL
Markers_GO_enrich_TC5_up <- enrichGO(genelist_TC5_up, OrgDb = org.Mm.eg.db, keyType="SYMBOL", pAdjustMethod = 'BH', ont= "ALL", pvalueCutoff = 0.5)


markerlist_TC3_down <- as.character(TC3_down$gene)
eg1 <- bitr(markerlist_TC3_down , fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb = org.Mm.eg.db)
genelist_TC3_down <- eg1$SYMBOL
Markers_GO_enrich_TC3_down <- enrichGO(genelist_TC3_down, OrgDb = org.Mm.eg.db, keyType="SYMBOL", pAdjustMethod = 'BH', ont= "ALL", pvalueCutoff = 0.5)


markerlist_TC4_down <- as.character(TC4_down$gene)
eg2 <- bitr(markerlist_TC4_down , fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb = org.Mm.eg.db)
genelist_TC4_down <- eg2$SYMBOL
Markers_GO_enrich_TC4_down <- enrichGO(genelist_TC4_down, OrgDb = org.Mm.eg.db, keyType="SYMBOL", pAdjustMethod = 'BH', ont= "ALL", pvalueCutoff = 0.5)


markerlist_TC5_down <- as.character(TC5_down$gene)
eg3 <- bitr(markerlist_TC5_down , fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb = org.Mm.eg.db)
genelist_TC5_down <- eg3$SYMBOL
Markers_GO_enrich_TC5_down <- enrichGO(genelist_TC5_down, OrgDb = org.Mm.eg.db, keyType="SYMBOL", pAdjustMethod = 'BH', ont= "ALL", pvalueCutoff = 0.5)


pdf(paste0("./BI_out/Epithelial_SOSMEK_TC3_TC4_TC5_upregulated_enrichment_allcategories.pdf"),width = 10,height = 10)
enrichplot::dotplot(Markers_GO_enrich_TC3_up, split = "ONTOLOGY",showCategory = 8) + facet_grid(ONTOLOGY~., scale = "free") +ggtitle("TC3 upregulated genes, GO enrichment") +
  theme(text = element_text(size = 9),
        axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.y = element_text(size =9),
        axis.title.x = element_text(size=9))

enrichplot::dotplot(Markers_GO_enrich_TC4_up, split = "ONTOLOGY",showCategory = 8) + facet_grid(ONTOLOGY~., scale = "free") +ggtitle("TC4 upregulated genes, GO enrichment") +
  theme(text = element_text(size = 9),
        axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.y = element_text(size =9),
        axis.title.x = element_text(size=9))

enrichplot::dotplot(Markers_GO_enrich_TC5_up, split = "ONTOLOGY",showCategory = 8) + facet_grid(ONTOLOGY~., scale = "free") +ggtitle("TC5 upregulated genes, GO enrichment") +
  theme(text = element_text(size = 9),
        axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.y = element_text(size =9),
        axis.title.x = element_text(size=9))

dev.off()

pdf("./BI_out/Epithelial_SOSMEK_TC3_TC4_TC5_downregulated_enrichment_allcategories.pdf")
enrichplot::dotplot(Markers_GO_enrich_TC3_down, split = "ONTOLOGY",showCategory = 8) + facet_grid(ONTOLOGY~., scale = "free") +ggtitle("TC3 downregulated genes, GO enrichment")+
  theme(text = element_text(size = 9),
        axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.y = element_text(size =9),
        axis.title.x = element_text(size=9))
enrichplot::dotplot(Markers_GO_enrich_TC4_down, split = "ONTOLOGY",showCategory = 8) + facet_grid(ONTOLOGY~., scale = "free") +ggtitle("TC4 downregulated genes, GO enrichment")+
  theme(text = element_text(size = 9),
        axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.y = element_text(size =9),
        axis.title.x = element_text(size=9))
enrichplot::dotplot(Markers_GO_enrich_TC5_down, split = "ONTOLOGY",showCategory = 8) + facet_grid(ONTOLOGY~., scale = "free") +ggtitle("TC5 downregulated genes, GO enrichment")+
  theme(text = element_text(size = 9),
        axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=7),
        axis.title.y = element_text(size =9),
        axis.title.x = element_text(size=9))
dev.off()

#save the table
Epithelial_TC3_upregulated_GO_enrich_allcategories <- Markers_GO_enrich_TC3_up@result
write.csv(Epithelial_TC3_upregulated_GO_enrich_allcategories, "./BI_out/Epithelial_SOSMEK_TC3_upregulated_GO_enrich_allcategories.csv")
Epithelial_TC4_upregulated_GO_enrich_allcategories <- Markers_GO_enrich_TC4_up@result
write.csv(Epithelial_TC4_upregulated_GO_enrich_allcategories, "./BI_out/Epithelial_SOSMEK_TC4_upregulated_GO_enrich_allcategories.csv")
Epithelial_TC5_upregulated_GO_enrich_allcategories <- Markers_GO_enrich_TC5_up@result
write.csv(Epithelial_TC4_upregulated_GO_enrich_allcategories, "./BI_out/Epithelial_SOSMEK_TC5_upregulated_GO_enrich_allcategories.csv")

#save the table
Epithelial_TC3_downregulated_GO_enrich_allcategories <- Markers_GO_enrich_TC3_down@result
write.csv(Epithelial_TC3_downregulated_GO_enrich_allcategories, "./BI_out/Epithelial_SOSMEK_TC3_downregulated_GO_enrich_allcategories.csv")
Epithelial_TC4_downregulated_GO_enrich_allcategories <- Markers_GO_enrich_TC4_down@result
write.csv(Epithelial_TC4_downregulated_GO_enrich_allcategories, "./BI_out/Epithelial_SOSMEK_TC4_downregulated_GO_enrich_allcategories.csv")
Epithelial_TC5_downregulated_GO_enrich_allcategories <- Markers_GO_enrich_TC5_up@result
write.csv(Epithelial_TC4_downregulated_GO_enrich_allcategories, "./BI_out/Epithelial_SOSMEK_TC5_downregulated_GO_enrich_allcategories.csv")

#for venn
set_TC3_up<- c(TC3_up$gene)
set_TC4_up<- c(TC4_up$gene)
set_TC5_up <- c(TC5_up$gene)
set_TC3_down <- c(TC3_down$gene)
set_TC4_down<- c(TC4_down$gene)
set_TC5_down <- c(TC5_down$gene)

#venn diagram for up
venn.diagram(
  x = list(set_TC3_up, set_TC4_up, set_TC5_up),
  category.names = c("TC3 upregulated" , "TC4 upregulated " , "TC5 upregulated"),
  filename = './BI_out/Epithelial_SOSMEK_upregulated_venn.png',
  output = TRUE ,
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  col=c("#440154ff", '#21908dff', '#924900'),
  fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3), alpha('#924900',0.3)),
  cex = 0.5,
  fontfamily = "sans",
  cat.cex = 0.3,
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  cat.col = c("#440154ff", '#21908dff', '#924900'),
  rotation = 1
)

#venn diagram for down
venn.diagram(
  x = list(set_TC3_down, set_TC4_down, set_TC5_down),
  category.names = c("TC3 downregulated" , "TC4 downregulated " , "TC5 downregulated"),
  filename = './BI_out/Epithelial_SOSMEK_downregulated_venn.png',
  output = TRUE ,
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  col=c("#440154ff", '#21908dff', '#924900'),
  fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3), alpha('#924900',0.3)),
  cex = 0.5,
  fontfamily = "sans",
  cat.cex = 0.3,
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  cat.col = c("#440154ff", '#21908dff', '#924900'),
  rotation = 1
)

#save the common genes
Epithelial_SOSMEK_common_up <- intersect(intersect(set_TC3_up,set_TC4_up),set_TC5_up)
write.csv(Epithelial_SOSMEK_common_up, "./BI_out/Epithelial_SOSMEK_common_up.csv", row.names = F)
Epithelial_SOSMEK_common_down <- intersect(intersect(set_TC3_down,set_TC4_down),set_TC5_down)
write.csv(Epithelial_SOSMEK_common_down, "./BI_out/Epithelial_SOSMEK_common_down.csv", row.names = F)

#save other common
Epithelial_SOSMEK_only_TC3_TC4_common_up <- setdiff(intersect(set_TC3_up, set_TC4_up), set_TC5_up)
Epithelial_SOSMEK_only_TC4_TC5_common_up <- setdiff(intersect(set_TC4_up, set_TC5_up), set_TC3_up)
Epithelial_SOSMEK_only_TC3_TC5_common_up <- setdiff(intersect(set_TC3_up, set_TC5_up), set_TC4_up)

write.csv(Epithelial_SOSMEK_only_TC3_TC4_common_up, "./BI_out/Epithelial_SOSMEK_only_TC3_TC4_common_up.csv", row.names = F)
write.csv(Epithelial_SOSMEK_only_TC4_TC5_common_up, "./BI_out/Epithelial_SOSMEK_only_TC4_TC5_common_up.csv", row.names = F)
write.csv(Epithelial_SOSMEK_only_TC3_TC5_common_up, "./BI_out/Epithelial_SOSMEK_only_TC3_TC5_common_up.csv", row.names = F)


Epithelial_SOSMEK_only_TC3_TC4_common_down <- setdiff(intersect(set_TC3_down, set_TC4_down), set_TC5_down)
Epithelial_SOSMEK_only_TC4_TC5_common_down <- setdiff(intersect(set_TC4_down, set_TC5_down), set_TC3_down)
Epithelial_SOSMEK_only_TC3_TC5_common_down <- setdiff(intersect(set_TC3_down, set_TC5_down), set_TC4_down)


write.csv(Epithelial_SOSMEK_only_TC3_TC4_common_down, "./BI_out/Epithelial_SOSMEK_only_TC3_TC4_common_down.csv", row.names = F)
write.csv(Epithelial_SOSMEK_only_TC4_TC5_common_down, "./BI_out/Epithelial_SOSMEK_only_TC4_TC5_common_down.csv", row.names = F)
write.csv(Epithelial_SOSMEK_only_TC3_TC5_common_down, "./BI_out/Epithelial_SOSMEK_only_TC3_TC5_common_down.csv", row.names = F)

#save
Epithelial_SOSMEK_only_TC3_up <- setdiff(setdiff(set_TC3_up, set_TC4_up), set_TC5_up)
Epithelial_SOSMEK_only_TC4_up <- setdiff(setdiff(set_TC4_up, set_TC5_up), set_TC3_up)
Epithelial_SOSMEK_only_TC5_up <- setdiff(setdiff(set_TC5_up, set_TC3_up), set_TC4_up)

write.csv(Epithelial_SOSMEK_only_TC3_up, "./BI_out/Epithelial_SOSMEK_only_TC3_up.csv", row.names = F)
write.csv(Epithelial_SOSMEK_only_TC4_up, "./BI_out/Epithelial_SOSMEK_only_TC4_up.csv", row.names = F)
write.csv(Epithelial_SOSMEK_only_TC5_up, "./BI_out/Epithelial_SOSMEK_only_TC5_up.csv", row.names = F)


Epithelial_SOSMEK_only_TC3_down <- setdiff(setdiff(set_TC3_down, set_TC4_down), set_TC5_down)
Epithelial_SOSMEK_only_TC4_down <- setdiff(setdiff(set_TC4_down, set_TC5_down), set_TC3_down)
Epithelial_SOSMEK_only_TC5_down <- setdiff(setdiff(set_TC5_down, set_TC3_down), set_TC4_down)


write.csv(Epithelial_SOSMEK_only_TC3_down, "./BI_out/Epithelial_SOSMEK_only_TC3_down.csv", row.names = F)
write.csv(Epithelial_SOSMEK_only_TC4_down, "./BI_out/Epithelial_SOSMEK_only_TC4_down.csv", row.names = F)
write.csv(Epithelial_SOSMEK_only_TC5_down , "./BI_out/Epithelial_SOSMEK_only_TC5_down.csv", row.names = F)




#This part is for DGE anlysis based on pesudo bulk DGE in each cell line between 2 treatments, making heatmap and volcano plots, runnign GO enrichment analysis and finding the overlap of up/down regulated gene among 3 cell lines
Epithelial <- subset(BI.filter_YFP_AddModuleScore, subset= MajorCellType == "Epithelial" |  MajorCellType == "Proliferating Epithelial") 

Epithelial@meta.data$celltype <- "Epithelial"

#extract the counts and metadata to make singleCellExperiment
counts <- Epithelial@assays$RNA@counts
metadata <- Epithelial@meta.data

#vector is a list of atomic values, factor is a list of vectors, set up metadata as desired for aggregation
#metadata$MajorCellType <- factor(Epithelial@active.ident)
#metadata$MajorCellType[metadata$MajorCellType == "Proliferating Epithelial"] <- "Epithelial"
metadata$celltype <- factor(Epithelial@meta.data$celltype)
metadata$Sample <- factor(Epithelial@meta.data$Sample)
metadata$Treatment <- factor(Epithelial@meta.data$Treatment)

#create SingleCellExperiment
sce <- SingleCellExperiment(assays = list(counts = counts), 
                            colData = metadata)

#Identify groups for aggregation of counts, (metadata is the colData in sce)
groups <- colData(sce)[, c("celltype", "Sample")]


assays(sce)
dim(counts(sce))
counts(sce)[1:6, 1:6]
dim(colData(sce))
head(colData(sce))

#determine how many cluster and how many samples do we have and what is the length of them
kids <- purrr::set_names(levels(sce$celltype))

kids
# Total number of clusters
nk <- length(kids)
nk

samples <- purrr::set_names(levels(sce$Sample))
# Total number of samples 
n_samples <- length(samples)
n_samples

#generate sample level metadata
#get the number of cells for each sample
table(sce$Sample)
#only get the numbers
n_cells <- as.numeric(table(sce$Sample))
#determine how to reorder the samples or the row of metadata to match the order of sample names in samples (look above)
m <- match(samples, sce$Sample)

#make the metadata, n_cells is a column in ei metadata now, colnames(ei). In ei maybe the rownames and the n_cells column is important
ei <- data.frame(colData(sce)[m, ], 
                 n_cells, row.names = NULL) %>% 
  select(-"celltype")
ei

#aggregate the matrix per sample (per sample per MajorCellType)
groups <- colData(sce)[, c("celltype", "Sample")]

# Aggregate across cluster-sample groups
pb <- aggregate.Matrix(t(counts(sce)), 
                       groupings = groups, fun = "sum")

class(pb)

dim(pb)

pb[1:6, 1:6]

#split the cell type
splitf <- sapply(stringr::str_split(rownames(pb), 
                                    pattern = "_",  
                                    n = 2), 
                 `[`, 1)

#split the matrix data by celltype and that the genes are row names and the sample names are the column
pb <- split.data.frame(pb, 
                       factor(splitf)) %>%
  lapply(function(u) 
    set_colnames(t(u), stringr::str_extract(rownames(u), "(?<=_)[:alnum:]+[:punct:]+[:alnum:]+")))

class(pb)
# Explore the different components of list
str(pb)

#the counts per sample per major cell type can be checked
options(width = 100)
table(sce$celltype, sce$Sample)

#make sample level metadata, (maybe all samples are not present in each cell type)
get_sample_ids <- function(x){
  pb[[x]] %>%
    colnames()
}
de_samples <- map(1:length(kids), get_sample_ids) %>%
  unlist()


samples_list <- map(1:length(kids), get_sample_ids)
get_cluster_ids <- function(x){
  rep(names(pb)[x], 
      each = length(samples_list[[x]]))
}

de_cluster_ids <- map(1:length(kids), get_cluster_ids) %>%
  unlist()

gg_df <- data.frame(celltype = de_cluster_ids,
                    Sample = de_samples)
gg_df <- left_join(gg_df, ei[, c("Sample", "Treatment")]) 

metadata <- gg_df %>%
  dplyr::select(celltype, Sample, Treatment) 

metadata

clusters <- levels(factor(metadata$celltype))
clusters

clusters[1]

cluster_metadata <- metadata[which(metadata$celltype == clusters[1]), ]
rownames(cluster_metadata) <- cluster_metadata$Sample
head(cluster_metadata)

counts <- pb[[clusters[1]]]
cluster_counts <- data.frame(counts[, which(colnames(counts) %in% rownames(cluster_metadata))])
all(rownames(cluster_metadata) == colnames(cluster_counts))

#subset based on the cell line
TC3 <- cluster_counts[which(colnames(cluster_counts) ==c("scGEX_152", "scGEX_153", "scGEX_154", "scGEX_155", "scGEX_156", "scGEX_157", "scGEX_158"))]
cluster_metadata_TC3 <- cluster_metadata[which(rownames(cluster_metadata) ==c("scGEX_152", "scGEX_153", "scGEX_154", "scGEX_155", "scGEX_156", "scGEX_157", "scGEX_158")),]

TC4 <- cluster_counts[which(colnames(cluster_counts) ==c("scGEX_159", "scGEX_160", "scGEX_161", "scGEX_162", "scGEX_163", "scGEX_164", "scGEX_165"))]
cluster_metadata_TC4 <- cluster_metadata[which(rownames(cluster_metadata) ==c("scGEX_159", "scGEX_160", "scGEX_161", "scGEX_162", "scGEX_163", "scGEX_164", "scGEX_165")),]

TC5 <- cluster_counts[which(colnames(cluster_counts) ==c("scGEX_166", "scGEX_167", "scGEX_168", "scGEX_169", "scGEX_170", "scGEX_171", "scGEX_172"))]
cluster_metadata_TC5 <- cluster_metadata[which(rownames(cluster_metadata) ==c("scGEX_166", "scGEX_167", "scGEX_168", "scGEX_169", "scGEX_170", "scGEX_171", "scGEX_172")),]

#make the DEseq object
dds_TC3 <- DESeqDataSetFromMatrix(TC3, 
                                  colData = cluster_metadata_TC3, 
                                  design = ~ Treatment)

dds_TC4 <- DESeqDataSetFromMatrix(TC4, 
                                  colData = cluster_metadata_TC4, 
                                  design = ~ Treatment)

dds_TC5 <- DESeqDataSetFromMatrix(TC5, 
                                  colData = cluster_metadata_TC5, 
                                  design = ~ Treatment)

#plotting based on Treatmant in each cell line
rld_TC3 <- rlog(dds_TC3, blind=TRUE)
DESeq2::plotPCA(rld_TC3, intgroup = "Treatment")

rld_TC4 <- rlog(dds_TC4, blind=TRUE)
DESeq2::plotPCA(rld_TC4, intgroup = "Treatment")

rld_TC5 <- rlog(dds_TC5, blind=TRUE)
DESeq2::plotPCA(rld_TC5, intgroup = "Treatment")

#plot that shows sample relationship for each cell line
rld_mat_TC3 <- assay(rld_TC3)
rld_cor_TC3 <- cor(rld_mat_TC3)
pheatmap(rld_cor_TC3, annotation = cluster_metadata_TC3[ c("celltype", "Treatment"), drop=F])

rld_mat_TC4 <- assay(rld_TC4)
rld_cor_TC4 <- cor(rld_mat_TC4)
pheatmap(rld_cor_TC4, annotation = cluster_metadata_TC4[ c("celltype", "Treatment"), drop=F])

rld_mat_TC5 <- assay(rld_TC5)
rld_cor_TC5 <- cor(rld_mat_TC5)
pheatmap(rld_cor_TC5, annotation = cluster_metadata_TC5[ c("celltype", "Treatment"), drop=F])

#for DE analysis 
dds_TC3 <- DESeq(dds_TC3)
dds_TC4 <- DESeq(dds_TC4)
dds_TC5 <- DESeq(dds_TC5)

plotDispEsts(dds_TC3)
plotDispEsts(dds_TC4)
plotDispEsts(dds_TC5)
levels(factor(cluster_metadata_TC3$Treatment))[2]
levels(factor(cluster_metadata_TC3$Treatment))[1]
contrast_TC3 <- c("Treatment", levels(factor(cluster_metadata_TC3$Treatment))[2], levels(factor(cluster_metadata_TC3$Treatment))[1])

levels(factor(cluster_metadata_TC4$Treatment))[2]
levels(factor(cluster_metadata_TC4$Treatment))[1]
contrast_TC4 <- c("Treatment", levels(factor(cluster_metadata_TC4$Treatment))[2], levels(factor(cluster_metadata_TC4$Treatment))[1])

levels(factor(cluster_metadata_TC5$Treatment))[2]
levels(factor(cluster_metadata_TC5$Treatment))[1]
contrast_TC5 <- c("Treatment", levels(factor(cluster_metadata_TC5$Treatment))[2], levels(factor(cluster_metadata_TC5$Treatment))[1])


res_TC3 <- results(dds_TC3, 
                   contrast = contrast_TC3,
                   alpha = 0.05)

# res_TC3 <- lfcShrink(dds_TC3,
#                  contrast =  contrast_TC3,
#                  res=res_TC3)
res_TC4 <- results(dds_TC4, 
                   contrast = contrast_TC4,
                   alpha = 0.05)

res_TC5 <- results(dds_TC5, 
                   contrast = contrast_TC5,
                   alpha = 0.05)

res_tbl_TC3 <- res_TC3 %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>%
  as_tibble()

res_tbl_TC4 <- res_TC4 %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>%
  as_tibble()

res_tbl_TC5 <- res_TC5 %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>%
  as_tibble()

#save all the genes
write.csv(res_tbl_TC3,
          paste0("./BI_out/", "TC3_DEseq2_results_", clusters[1], "_", levels(cluster_metadata_TC3$Treatment)[2], "_vs_", levels(cluster_metadata_TC3$Treatment)[1], "_all_genes.csv"),
          quote = FALSE, 
          row.names = FALSE)


write.csv(res_tbl_TC4,
          paste0("./BI_out/", "TC4_DEseq2_results_", clusters[1], "_", levels(cluster_metadata_TC4$Treatment)[2], "_vs_", levels(cluster_metadata_TC4$Treatment)[1], "_all_genes.csv"),
          quote = FALSE, 
          row.names = FALSE)

write.csv(res_tbl_TC5,
          paste0("./BI_out/", "TC5_DEseq2_results_", clusters[1], "_", levels(cluster_metadata_TC5$Treatment)[2], "_vs_", levels(cluster_metadata_TC5$Treatment)[1], "_all_genes.csv"),
          quote = FALSE, 
          row.names = FALSE)

#filter and save
padj_cutoff <- 0.05
sig_res_TC3 <- dplyr::filter(res_tbl_TC3, padj < padj_cutoff) %>%
  dplyr::arrange(padj)
sig_res_TC3
write.csv(sig_res_TC3,
          paste0("./BI_out/", "TC3_DEseq2_results_", clusters[1], "_" , levels(cluster_metadata_TC3$Treatment)[2], "_vs_", levels(cluster_metadata_TC3$Treatment)[1], "_sig_genes.csv"),
          quote = FALSE, 
          row.names = FALSE)


padj_cutoff <- 0.05
sig_res_TC4 <- dplyr::filter(res_tbl_TC4, padj < padj_cutoff) %>%
  dplyr::arrange(padj)
sig_res_TC4
write.csv(sig_res_TC4,
          paste0("./BI_out/", "TC4_DEseq2_results_", clusters[1], "_" , levels(cluster_metadata_TC4$Treatment)[2], "_vs_", levels(cluster_metadata_TC4$Treatment)[1], "_sig_genes.csv"),
          quote = FALSE, 
          row.names = FALSE)


padj_cutoff <- 0.05
sig_res_TC5 <- dplyr::filter(res_tbl_TC5, padj < padj_cutoff) %>%
  dplyr::arrange(padj)
sig_res_TC5
write.csv(sig_res_TC5,
          paste0("./BI_out/", "TC5_DEseq2_results_", clusters[1], "_" , levels(cluster_metadata_TC5$Treatment)[2], "_vs_", levels(cluster_metadata_TC5$Treatment)[1], "_sig_genes.csv"),
          quote = FALSE, 
          row.names = FALSE)

#get the normalized counts
# ggplot of top genes
normalized_counts_TC3 <- counts(dds_TC3, 
                                normalized = TRUE)

# Order results by padj values
top20_sig_genes_TC3 <- sig_res_TC3 %>%
  dplyr::arrange(padj) %>%
  dplyr::pull(gene) %>%
  head(n=20)


top20_sig_norm_TC3 <- data.frame(normalized_counts_TC3) %>%
  rownames_to_column(var = "gene") %>%
  dplyr::filter(gene %in% top20_sig_genes_TC3)

gathered_top20_sig_TC3 <- top20_sig_norm_TC3 %>%
  gather(colnames(top20_sig_norm_TC3)[2:length(colnames(top20_sig_norm_TC3))], key = "samplename", value = "normalized_counts")

gathered_top20_sig_TC3 <- inner_join(ei[, c("Sample", "Treatment" )], gathered_top20_sig_TC3, by = c("Sample" = "samplename"))

# plot using ggplot2
ggplot(gathered_top20_sig_TC3) +
  geom_point(aes(x = gene, 
                 y = normalized_counts, 
                 color = Treatment), 
             position=position_jitter(w=0.1,h=0)) +
  scale_y_log10() +
  xlab("Genes") +
  ylab("log10 Normalized Counts") +
  ggtitle("Top 20 Significant DE Genes in TC3") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  theme(plot.title = element_text(hjust = 0.5))

# ggplot of top genes
normalized_counts_TC4 <- counts(dds_TC4, 
                                normalized = TRUE)

# Order results by padj values
top20_sig_genes_TC4 <- sig_res_TC4 %>%
  dplyr::arrange(padj) %>%
  dplyr::pull(gene) %>%
  head(n=20)


top20_sig_norm_TC4 <- data.frame(normalized_counts_TC4) %>%
  rownames_to_column(var = "gene") %>%
  dplyr::filter(gene %in% top20_sig_genes_TC4)

gathered_top20_sig_TC4 <- top20_sig_norm_TC4 %>%
  gather(colnames(top20_sig_norm_TC4)[2:length(colnames(top20_sig_norm_TC4))], key = "samplename", value = "normalized_counts")

gathered_top20_sig_TC4 <- inner_join(ei[, c("Sample", "Treatment" )], gathered_top20_sig_TC4, by = c("Sample" = "samplename"))

## plot using ggplot2
ggplot(gathered_top20_sig_TC4) +
  geom_point(aes(x = gene, 
                 y = normalized_counts, 
                 color = Treatment), 
             position=position_jitter(w=0.1,h=0)) +
  scale_y_log10() +
  xlab("Genes") +
  ylab("log10 Normalized Counts") +
  ggtitle("Top 20 Significant DE Genes in TC4") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  theme(plot.title = element_text(hjust = 0.5))

## ggplot of top genes
normalized_counts_TC5 <- counts(dds_TC5, 
                                normalized = TRUE)

## Order results by padj values
top20_sig_genes_TC5 <- sig_res_TC5 %>%
  dplyr::arrange(padj) %>%
  dplyr::pull(gene) %>%
  head(n=20)


top20_sig_norm_TC5 <- data.frame(normalized_counts_TC5) %>%
  rownames_to_column(var = "gene") %>%
  dplyr::filter(gene %in% top20_sig_genes_TC5)

gathered_top20_sig_TC5 <- top20_sig_norm_TC5 %>%
  gather(colnames(top20_sig_norm_TC5)[2:length(colnames(top20_sig_norm_TC5))], key = "samplename", value = "normalized_counts")

gathered_top20_sig_TC5 <- inner_join(ei[, c("Sample", "Treatment" )], gathered_top20_sig_TC5, by = c("Sample" = "samplename"))

## plot using ggplot2
ggplot(gathered_top20_sig_TC5) +
  geom_point(aes(x = gene, 
                 y = normalized_counts, 
                 color = Treatment), 
             position=position_jitter(w=0.1,h=0)) +
  scale_y_log10() +
  xlab("Genes") +
  ylab("log10 Normalized Counts") +
  ggtitle("Top 20 Significant DE Genes in TC5") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  theme(plot.title = element_text(hjust = 0.5))


#heat map of significant genes
pdf("./BI_out/TC3_DEseq2_results_Epithellial_Vehicle_vs_SOSMEK_sig_genes_heatmap.pdf")
# Extract normalized counts for only the significant genes
sig_norm_TC3 <- data.frame(normalized_counts_TC3) %>%
  rownames_to_column(var = "gene") %>%
  dplyr::filter(gene %in% sig_res_TC3$gene)

# Set a color palette
heat_colors <- brewer.pal(6, "YlOrRd")

# Run pheatmap using the metadata data frame for the annotation
pheatmap(sig_norm_TC3[ , 2:length(colnames(sig_norm_TC3))], 
         color = heat_colors, 
         cluster_rows = T, 
         show_rownames = F,
         annotation = cluster_metadata_TC3[, c("Treatment", "celltype")],
         #make annotation file here
         #annotation = SampleTable[, "condition"],
         border_color = NA, 
         fontsize = 10, 
         scale = "row", 
         fontsize_row = 10, 
         height = 20) 
dev.off()

pdf("./BI_out/TC4_DEseq2_results_Epithellial_Vehicle_vs_SOSMEK_sig_genes_heatmap.pdf")
# Extract normalized counts for only the significant genes
sig_norm_TC4 <- data.frame(normalized_counts_TC4) %>%
  rownames_to_column(var = "gene") %>%
  dplyr::filter(gene %in% sig_res_TC4$gene)

# Set a color palette
heat_colors <- brewer.pal(6, "YlOrRd")

# Run pheatmap using the metadata data frame for the annotation
pheatmap(sig_norm_TC4[ , 2:length(colnames(sig_norm_TC4))], 
         color = heat_colors, 
         cluster_rows = T, 
         show_rownames = F,
         annotation = cluster_metadata_TC4[, c("Treatment", "celltype")],
         #make annotation file here
         #annotation = SampleTable[, "condition"],
         border_color = NA, 
         fontsize = 10, 
         scale = "row", 
         fontsize_row = 10, 
         height = 20) 
dev.off()

pdf("./BI_out/TC5_DEseq2_results_Epithellial_Vehicle_vs_SOSMEK_sig_genes_heatmap.pdf")
# Extract normalized counts for only the significant genes
sig_norm_TC5 <- data.frame(normalized_counts_TC5) %>%
  rownames_to_column(var = "gene") %>%
  dplyr::filter(gene %in% sig_res_TC5$gene)

# Set a color palette
heat_colors <- brewer.pal(6, "YlOrRd")

# Run pheatmap using the metadata data frame for the annotation
pheatmap(sig_norm_TC5[ , 2:length(colnames(sig_norm_TC5))], 
         color = heat_colors, 
         cluster_rows = T, 
         show_rownames = F,
         annotation = cluster_metadata_TC5[, c("Treatment", "celltype")],
         #make annotation file here
         #annotation = SampleTable[, "condition"],
         border_color = NA, 
         fontsize = 10, 
         scale = "row", 
         fontsize_row = 10, 
         height = 20) 
dev.off()

#volcano plot
res_table_thres_TC3 <- res_tbl_TC3 %>% 
  mutate(threshold = padj < 0.05 & abs(log2FoldChange) >= 0.58)

## Volcano plot
ggplot(res_table_thres_TC3) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold)) +
  ggtitle("Volcano plot of Epithelial significant DEGs from DESeq2 (vehicle vs. SOSMEK)") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  scale_y_continuous(limits = c(0,50)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))

pdf(paste0("./BI_out/", "TC3_DEseq2_results_Epithellial_Vehicle_vs_SOSMEK_volcano.pdf"))
res_tbl_TC3$expression <- "not significant"
res_tbl_TC3$expression[res_tbl_TC3$pvalue < 0.05 & res_tbl_TC3$log2FoldChange > 0.5] <- "upregulated"
res_tbl_TC3$expression[res_tbl_TC3$pvalue < 0.05 &  res_tbl_TC3$log2FoldChange < -0.5] <- "downregulated"

top_genes <- res_tbl_TC3 %>% 
  filter(expression != 'not significant')
ggplot(data=res_tbl_TC3, aes(x=log2FoldChange, y=-log10(pvalue), col=expression)) + 
  geom_point(size=0.5) + 
  #theme_minimal() +
  theme_bw() +
  theme(text=element_text(size=11))+
  #geom_text_repel(data=head(markers_rcpp_TC3_subset, 50), aes(label=gene),nudge_x = .001 ,nudge_y = 0.001, size=3) +
  #geom_text_repel(data=markers_rcpp_TC3_subset, aes(label=gene),nudge_x = 0 ,nudge_y = 0, size=3) +
  #geom_text_repel(data=top_genes, aes(label=gene)) +
  geom_label_repel(data = head(top_genes, n= 20), 
                   mapping = aes(log2FoldChange, -log10(pvalue), label = gene), size = 2, max.overlaps=Inf)+
  
  theme(axis.text = element_text(face="bold"))+
  theme(axis.line.x = element_line(color="black", size = 0.3),
        axis.line.y = element_line(color="black", size = 0.3))+
  #geom _text is about the size of dots
  #geom_text(position=position_dodge(width=0.9),hjust=1.5, size= 0.1)+
  #scale_color_manual(values=c("red", "blue")) +
  scale_color_manual(values=c("blue", "dark grey", "red"))+
  geom_vline(xintercept=c(-0.5, 0.5), col="black", linetype="dashed") +
  geom_hline(yintercept=-log10(0.1), col="black", linetype="dashed") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())


pdf(paste0("./BI_out/", "TC4_DEseq2_results_Epithellial_Vehicle_vs_SOSMEK_volcano.pdf"))
res_tbl_TC4$expression <- "not significant"
res_tbl_TC4$expression[res_tbl_TC4$pvalue < 0.05 & res_tbl_TC4$log2FoldChange > 0.5] <- "upregulated"
res_tbl_TC4$expression[res_tbl_TC4$pvalue < 0.05 &  res_tbl_TC4$log2FoldChange < -0.5] <- "downregulated"

top_genes <- res_tbl_TC4 %>% 
  filter(expression != 'not significant')
ggplot(data=res_tbl_TC4, aes(x=log2FoldChange, y=-log10(pvalue), col=expression)) + 
  geom_point(size=0.5) + 
  #theme_minimal() +
  theme_bw() +
  theme(text=element_text(size=11))+
  #geom_text_repel(data=head(markers_rcpp_TC3_subset, 50), aes(label=gene),nudge_x = .001 ,nudge_y = 0.001, size=3) +
  #geom_text_repel(data=markers_rcpp_TC3_subset, aes(label=gene),nudge_x = 0 ,nudge_y = 0, size=3) +
  #geom_text_repel(data=top_genes, aes(label=gene)) +
  geom_label_repel(data = head(top_genes, n= 20), 
                   mapping = aes(log2FoldChange, -log10(pvalue), label = gene), size = 2, max.overlaps=Inf)+
  
  theme(axis.text = element_text(face="bold"))+
  theme(axis.line.x = element_line(color="black", size = 0.3),
        axis.line.y = element_line(color="black", size = 0.3))+
  #geom _text is about the size of dots
  #geom_text(position=position_dodge(width=0.9),hjust=1.5, size= 0.1)+
  #scale_color_manual(values=c("red", "blue")) +
  scale_color_manual(values=c("blue", "dark grey", "red"))+
  geom_vline(xintercept=c(-0.5, 0.5), col="black", linetype="dashed") +
  geom_hline(yintercept=-log10(0.1), col="black", linetype="dashed") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
dev.off()

pdf(paste0("./BI_out/", "TC5_DEseq2_results_Epithellial_Vehicle_vs_SOSMEK_volcano.pdf"))
res_tbl_TC5$expression <- "not significant"
res_tbl_TC5$expression[res_tbl_TC5$pvalue < 0.05 & res_tbl_TC5$log2FoldChange > 0.5] <- "upregulated"
res_tbl_TC5$expression[res_tbl_TC5$pvalue < 0.05 &  res_tbl_TC5$log2FoldChange < -0.5] <- "downregulated"

top_genes <- res_tbl_TC5 %>% 
  filter(expression != 'not significant')
ggplot(data=res_tbl_TC5, aes(x=log2FoldChange, y=-log10(pvalue), col=expression)) + 
  geom_point(size=0.5) + 
  #theme_minimal() +
  theme_bw() +
  theme(text=element_text(size=11))+
  #geom_text_repel(data=head(markers_rcpp_TC3_subset, 50), aes(label=gene),nudge_x = .001 ,nudge_y = 0.001, size=3) +
  #geom_text_repel(data=markers_rcpp_TC3_subset, aes(label=gene),nudge_x = 0 ,nudge_y = 0, size=3) +
  #geom_text_repel(data=top_genes, aes(label=gene)) +
  geom_label_repel(data = head(top_genes, n= 20), 
                   mapping = aes(log2FoldChange, -log10(pvalue), label = gene), size = 2, max.overlaps=Inf)+
  
  theme(axis.text = element_text(face="bold"))+
  theme(axis.line.x = element_line(color="black", size = 0.3),
        axis.line.y = element_line(color="black", size = 0.3))+
  #geom _text is about the size of dots
  #geom_text(position=position_dodge(width=0.9),hjust=1.5, size= 0.1)+
  #scale_color_manual(values=c("red", "blue")) +
  scale_color_manual(values=c("blue", "dark grey", "red"))+
  geom_vline(xintercept=c(-0.5, 0.5), col="black", linetype="dashed") +
  geom_hline(yintercept=-log10(0.1), col="black", linetype="dashed") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
dev.off()

TC3 <- read.csv("./BI_out/TC3_DEseq2_results_Epithelial_Vehicle_vs_SOSMEK_sig_genes.csv")
TC3_up <- subset(TC3, log2FoldChange > 0)
TC3_down <- subset(TC3, log2FoldChange < 0)
TC4 <- read.csv("./BI_out/TC4_DEseq2_results_Epithelial_Vehicle_vs_SOSMEK_sig_genes.csv")
TC4_up <- subset(TC4, log2FoldChange > 0)
TC4_down <- subset(TC4, log2FoldChange < 0)
TC5 <- read.csv("./BI_out/TC5_DEseq2_results_Epithelial_Vehicle_vs_SOSMEK_sig_genes.csv")
TC5_up <- subset(TC5, log2FoldChange > 0)
TC5_down <- subset(TC5, log2FoldChange < 0)

#for venn

set_TC3_up <- c(TC3_up$gene)
set_TC4_up<- c(TC4_up$gene)
set_TC5_up <- c(TC5_up$gene)
set_TC3_down <- c(TC3_down$gene)
set_TC4_down<- c(TC4_down$gene)
set_TC5_down <- c(TC5_down$gene)

#venn diagram for up
venn.diagram(
  x = list(set_TC3_up, set_TC4_up, set_TC5_up),
  category.names = c("TC3 upregulated" , "TC4 upregulated " , "TC5 upregulated"),
  filename = './BI_out/Epithelial_DEseq2_upregulated_venn.png',
  output = TRUE ,
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  col=c("#440154ff", '#21908dff', '#924900'),
  fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3), alpha('#924900',0.3)),
  cex = 0.5,
  fontfamily = "sans",
  cat.cex = 0.3,
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  cat.col = c("#440154ff", '#21908dff', '#924900'),
  rotation = 1
)

#venn diagram for down
venn.diagram(
  x = list(set_TC3_down, set_TC4_down, set_TC5_down),
  category.names = c("TC3 downregulated" , "TC4 downregulated " , "TC5 downregulated"),
  filename = './BI_out/Epithelial_DEseq2_downregulated_venn.png',
  output = TRUE ,
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  col=c("#440154ff", '#21908dff', '#924900'),
  fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3), alpha('#924900',0.3)),
  cex = 0.5,
  fontfamily = "sans",
  cat.cex = 0.3,
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  cat.col = c("#440154ff", '#21908dff', '#924900'),
  rotation = 1
)

#save the common genes
Epithelial_Desq_common_up <- intersect(intersect(set_TC3_up,set_TC4_up),set_TC5_up)
write.csv(Epithelial_Desq_common_up, "./BI_out/Epithelial_DEseq2_common_up.csv", row.names = F)
Epithelial_Desq_common_down <- intersect(intersect(set_TC3_down,set_TC4_down),set_TC5_down)
write.csv(Epithelial_Desq_common_down, "./BI_out/Epithelial_DEseq2_common_down.csv", row.names = F)

#save other common
Epithelial_Desq_only_TC3_TC4_common_up <- setdiff(intersect(set_TC3_up, set_TC4_up), set_TC5_up)
Epithelial_Desq_only_TC4_TC5_common_up <- setdiff(intersect(set_TC4_up, set_TC5_up), set_TC3_up)
Epithelial_Desq_only_TC3_TC5_common_up <- setdiff(intersect(set_TC3_up, set_TC5_up), set_TC4_up)

write.csv(Epithelial_Desq_only_TC3_TC4_common_up, "./BI_out/Epithelial_DEseq2_only_TC3_TC4_common_up.csv", row.names = F)
write.csv(Epithelial_Desq_only_TC4_TC5_common_up, "./BI_out/Epithelial_DEseq2_only_TC4_TC5_common_up.csv", row.names = F)
write.csv(Epithelial_Desq_only_TC3_TC5_common_up, "./BI_out/Epithelial_DEseq2_only_TC3_TC5_common_up.csv", row.names = F)


Epithelial_Desq_only_TC3_TC4_common_down <- setdiff(intersect(set_TC3_down, set_TC4_down), set_TC5_down)
Epithelial_Desq_only_TC4_TC5_common_down <- setdiff(intersect(set_TC4_down, set_TC5_down), set_TC3_down)
Epithelial_Desq_only_TC3_TC5_common_down <- setdiff(intersect(set_TC3_down, set_TC5_down), set_TC4_down)


write.csv(Epithelial_Desq_only_TC3_TC4_common_down, "./BI_out/Epithelial_DEseq2_only_TC3_TC4_common_down.csv", row.names = F)
write.csv(Epithelial_Desq_only_TC4_TC5_common_down, "./BI_out/Epithelial_DEseq2_only_TC4_TC5_common_down.csv", row.names = F)
write.csv(Epithelial_Desq_only_TC3_TC5_common_down, "./BI_out/Epithelial_DEseq2_only_TC3_TC5_common_down.csv", row.names = F)


#save
Epithelial_Desq_only_TC3_up <- setdiff(setdiff(set_TC3_up, set_TC4_up), set_TC5_up)
Epithelial_Desq_only_TC4_up <- setdiff(setdiff(set_TC4_up, set_TC5_up), set_TC3_up)
Epithelial_Desq_only_TC5_up <- setdiff(setdiff(set_TC5_up, set_TC3_up), set_TC4_up)

write.csv(Epithelial_Desq_only_TC3_up, "./BI_out/Epithelial_DEseq2_only_TC3_up.csv", row.names = F)
write.csv(Epithelial_Desq_only_TC4_up, "./BI_out/Epithelial_DEseq2_only_TC4_up.csv", row.names = F)
write.csv(Epithelial_Desq_only_TC5_up, "./BI_out/Epithelial_DEseq2_only_TC5_up.csv", row.names = F)


Epithelial_Desq_only_TC3_down <- setdiff(setdiff(set_TC3_down, set_TC4_down), set_TC5_down)
Epithelial_Desq_only_TC4_down <- setdiff(setdiff(set_TC4_down, set_TC5_down), set_TC3_down)
Epithelial_Desq_only_TC5_down <- setdiff(setdiff(set_TC5_down, set_TC3_down), set_TC4_down)


write.csv(Epithelial_Desq_only_TC3_down, "./BI_out/Epithelial_DEseq2_only_TC3_down.csv", row.names = F)
write.csv(Epithelial_Desq_only_TC4_down, "./BI_out/Epithelial_DEseq2_only_TC4_down.csv", row.names = F)
write.csv(Epithelial_Desq_only_TC5_down, "./BI_out/Epithelial_DEseq2_only_TC5_down.csv", row.names = F)

library(clusterProfiler)
library(org.Mm.eg.db)
library(org.Hs.eg.db)

markerlist_TC3_up <- as.character(TC3_up$gene)
eg1 <- bitr(markerlist_TC3_up , fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb = org.Mm.eg.db)
genelist_TC3_up <- eg1$SYMBOL
Markers_GO_enrich_TC3_up <- enrichGO(genelist_TC3_up, OrgDb = org.Mm.eg.db, keyType="SYMBOL", pAdjustMethod = 'BH', ont= "ALL", pvalueCutoff = 0.5)


markerlist_TC4_up <- as.character(TC4_up$gene)
eg2 <- bitr(markerlist_TC4_up , fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb = org.Mm.eg.db)
genelist_TC4_up <- eg2$SYMBOL
Markers_GO_enrich_TC4_up <- enrichGO(genelist_TC4_up, OrgDb = org.Mm.eg.db, keyType="SYMBOL", pAdjustMethod = 'BH', ont= "ALL", pvalueCutoff = 0.5)


markerlist_TC5_up <- as.character(TC5_up$gene)
eg3 <- bitr(markerlist_TC5_up , fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb = org.Mm.eg.db)
genelist_TC5_up <- eg3$SYMBOL
Markers_GO_enrich_TC5_up <- enrichGO(genelist_TC5_up, OrgDb = org.Mm.eg.db, keyType="SYMBOL", pAdjustMethod = 'BH', ont= "ALL", pvalueCutoff = 0.5)





markerlist_TC3_down <- as.character(TC3_down$gene)
eg1 <- bitr(markerlist_TC3_down , fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb = org.Mm.eg.db)
genelist_TC3_down <- eg1$SYMBOL
Markers_GO_enrich_TC3_down <- enrichGO(genelist_TC3_down, OrgDb = org.Mm.eg.db, keyType="SYMBOL", pAdjustMethod = 'BH', ont= "ALL", pvalueCutoff = 0.5)


markerlist_TC4_down <- as.character(TC4_down$gene)
eg2 <- bitr(markerlist_TC4_down , fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb = org.Mm.eg.db)
genelist_TC4_down <- eg2$SYMBOL
Markers_GO_enrich_TC4_down <- enrichGO(genelist_TC4_down, OrgDb = org.Mm.eg.db, keyType="SYMBOL", pAdjustMethod = 'BH', ont= "ALL", pvalueCutoff = 0.5)


markerlist_TC5_down <- as.character(TC5_down$gene)
eg3 <- bitr(markerlist_TC5_down , fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb = org.Mm.eg.db)
genelist_TC5_down <- eg3$SYMBOL
Markers_GO_enrich_TC5_down <- enrichGO(genelist_TC5_down, OrgDb = org.Mm.eg.db, keyType="SYMBOL", pAdjustMethod = 'BH', ont= "ALL", pvalueCutoff = 0.5)



pdf(paste0("./BI_out/Epithelial_DEseq2_TC3_TC4_TC5_upregulated_enrichment_allcategories.pdf"),width = 10,height = 10)
enrichplot::dotplot(Markers_GO_enrich_TC3_up, split = "ONTOLOGY",showCategory = 8) + facet_grid(ONTOLOGY~., scale = "free") +ggtitle("TC3 upregulated genes, GO enrichment") +
  theme(text = element_text(size = 9),
        axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.y = element_text(size =9),
        axis.title.x = element_text(size=9))

enrichplot::dotplot(Markers_GO_enrich_TC4_up, split = "ONTOLOGY",showCategory = 8) + facet_grid(ONTOLOGY~., scale = "free") +ggtitle("TC4 upregulated genes, GO enrichment") +
  theme(text = element_text(size = 9),
        axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.y = element_text(size =9),
        axis.title.x = element_text(size=9))

enrichplot::dotplot(Markers_GO_enrich_TC5_up, split = "ONTOLOGY",showCategory = 8) + facet_grid(ONTOLOGY~., scale = "free") +ggtitle("TC5 upregulated genes, GO enrichment") +
  theme(text = element_text(size = 9),
        axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.y = element_text(size =9),
        axis.title.x = element_text(size=9))

dev.off()

pdf("./BI_out/Epithelial_DEseq2_TC3_TC4_TC5_downregulated_enrichment_allcategories.pdf")
enrichplot::dotplot(Markers_GO_enrich_TC3_down, split = "ONTOLOGY",showCategory = 8) + facet_grid(ONTOLOGY~., scale = "free") +ggtitle("TC3 downregulated genes, GO enrichment")+
  theme(text = element_text(size = 9),
        axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.y = element_text(size =9),
        axis.title.x = element_text(size=9))
enrichplot::dotplot(Markers_GO_enrich_TC4_down, split = "ONTOLOGY",showCategory = 8) + facet_grid(ONTOLOGY~., scale = "free") +ggtitle("TC4 downregulated genes, GO enrichment")+
  theme(text = element_text(size = 9),
        axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.y = element_text(size =9),
        axis.title.x = element_text(size=9))
enrichplot::dotplot(Markers_GO_enrich_TC5_down, split = "ONTOLOGY",showCategory = 8) + facet_grid(ONTOLOGY~., scale = "free") +ggtitle("TC5 downregulated genes, GO enrichment")+
  theme(text = element_text(size = 9),
        axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9),
        axis.title.y = element_text(size =9),
        axis.title.x = element_text(size=9))
dev.off()

#save the table
Epithelial_TC3_upregulated_GO_enrich_allcategories <- Markers_GO_enrich_TC3_up@result
write.csv(Epithelial_TC3_upregulated_GO_enrich_allcategories, "./BI_out/Epithelial_DEseq2_TC3_upregulated_GO_enrich_allcategories.csv")
Epithelial_TC4_upregulated_GO_enrich_allcategories <- Markers_GO_enrich_TC4_up@result
write.csv(Epithelial_TC4_upregulated_GO_enrich_allcategories, "./BI_out/Epithelial_DEseq2_TC4_upregulated_GO_enrich_allcategories.csv")
Epithelial_TC5_upregulated_GO_enrich_allcategories <- Markers_GO_enrich_TC5_up@result
write.csv(Epithelial_TC4_upregulated_GO_enrich_allcategories, "./BI_out/Epithelial_DEseq2_TC5_upregulated_GO_enrich_allcategories.csv")

#save the table
Epithelial_TC3_downregulated_GO_enrich_allcategories <- Markers_GO_enrich_TC3_down@result
write.csv(Epithelial_TC3_downregulated_GO_enrich_allcategories, "./BI_out/Epithelial_DEseq2_TC3_downregulated_GO_enrich_allcategories.csv")
Epithelial_TC4_downregulated_GO_enrich_allcategories <- Markers_GO_enrich_TC4_down@result
write.csv(Epithelial_TC4_downregulated_GO_enrich_allcategories, "./BI_out/Epithelial_DEseq2_TC4_downregulated_GO_enrich_allcategories.csv")
Epithelial_TC5_downregulated_GO_enrich_allcategories <- Markers_GO_enrich_TC5_up@result
write.csv(Epithelial_TC4_downregulated_GO_enrich_allcategories, "./BI_out/Epithelial_DEseq2_TC5_downregulated_GO_enrich_allcategories.csv")





###Task 4: Pathway Scoring and Visualization on Epithelial/tumor cells

install.packages("qusage")
library("qusage")
gsets <- qusage::read.gmt("mh.all.v0.3.symbols.gmt")
cell.object <- AddModuleScore(Epithelial, features = gsets)
saveRDS(cell.object, "./BI_out/Epithelial_AddModuleScore.rds")

#Epithelials
colnames(Epithelial_AddModuleScore@meta.data)[21:70] <- names(gsets)

#features <- c(colnames(Epithelial_AddModuleScore@meta.data))[21:70]
features <- names(gsets)


#box plot comaprison
pdf(paste0("./BI_out/","Pathway_score_boxplot_Treatment_Cellline.pdf"),width = 20,height = 10)

features <- names(gsets)

library(RColorBrewer)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library("dplyr")
library(ggpubr)
compare_list = list(
  c('SOSMEK', 'Vehicle')
)
for (i in features){
  print(ggplot(aes ( x = Treatment, y = Epithelial_AddModuleScore@meta.data[,i]), data = Epithelial_AddModuleScore@meta.data) +
          geom_boxplot(inherit.aes = FALSE, outlier.shape =NA, aes(x = Treatment, y = Epithelial_AddModuleScore@meta.data[,i], fill= CellLine),outlier.color="black")+
          geom_point(color = 'darkblue')+
          facet_wrap( 'CellLine', scales = "fixed") + theme_hc() +
          stat_compare_means(comparisons = compare_list,method = 't.test')+
          theme(legend.position="none")+
          ylab(i)+
          theme(text = element_text(size = 16),
                axis.text.x=element_text(size=16),
                axis.text.y=element_text(size=16),
                axis.title.y = element_text(size =16),
                axis.title.x = element_text(size=16)))    
  
}

dev.off()

pdf(paste0("./BI_out/","Pathway_score_boxplot_Cellline_Treatment.pdf"),width = 20,height = 10)

features <- names(gsets)

library(RColorBrewer)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library("dplyr")
library(ggpubr)
compare_list = list(
  c('TC3','TC4'), c('TC3', 'TC5'), c('TC4', 'TC5'))
for (i in features){
  print(ggplot(aes ( x = CellLine, y = Epithelial_AddModuleScore@meta.data[,i]), data = Epithelial_AddModuleScore@meta.data) +
          geom_boxplot(inherit.aes = FALSE, outlier.shape =NA, aes(x = CellLine, y = Epithelial_AddModuleScore@meta.data[,i], fill= Treatment),outlier.color="black")+
          geom_point(color = 'darkblue')+
          facet_wrap( 'Treatment', scales = "free") + theme_hc() +
          stat_compare_means(comparisons = compare_list,method = 't.test')+
          theme(legend.position="none")+
          ylab(i)+
          theme(text = element_text(size = 16),
                axis.text.x=element_text(size=16),
                axis.text.y=element_text(size=16),
                axis.title.y = element_text(size =16),
                axis.title.x = element_text(size=16)))    
  
}

dev.off()

pdf(paste0("./BI_out/","Pathway_score_boxplot_Cellline_Treatment.pdf"),width = 20,height = 10)

features <- names(gsets)

library(RColorBrewer)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library("dplyr")
library(ggpubr)
compare_list = list(
  c('TC3','TC4'), c('TC3', 'TC5'), c('TC4', 'TC5'))
for (i in features){
  print(ggplot(aes ( x = CellLine, y = Epithelial_AddModuleScore@meta.data[,i]), data = Epithelial_AddModuleScore@meta.data) +
          geom_boxplot(inherit.aes = FALSE, outlier.shape =NA, aes(x = CellLine, y = Epithelial_AddModuleScore@meta.data[,i], fill= Treatment),outlier.color="black")+
          geom_point(color = 'darkblue')+
          facet_wrap( 'Treatment', scales = "free") + theme_hc() +
          stat_compare_means(comparisons = compare_list,method = 't.test')+
          theme(legend.position="none")+
          ylab(i)+
          theme(text = element_text(size = 16),
                axis.text.x=element_text(size=16),
                axis.text.y=element_text(size=16),
                axis.title.y = element_text(size =16),
                axis.title.x = element_text(size=16)))    
  
}

dev.off()





###Task 5: COX Pathway Gene Expression in Epithelial cells

rna = Matrix::t(Epithelial@assays$RNA@counts)
meta <- Epithelial@meta.data
gitr = function( mylist, d = rna, m = meta ) {
  cbind(
    m,
    as.matrix(d[, colnames(d) %in% mylist, drop = FALSE])
  )
}
genelist <- c("Ptgs1","Ptgs2","Ptges", "Ptger1", "Ptger2", "Ptger3", "Ptger4")
new_meta = gitr(genelist)

new_meta$expression.Ptgs1[new_meta[,'Ptgs1'] > 0] <- "Positive"
new_meta$expression.Ptgs1[new_meta[,'Ptgs1'] == 0] <- "Negative"
new_meta$expression.Ptgs2[new_meta[,'Ptgs2'] > 0] <- "Positive"
new_meta$expression.Ptgs2[new_meta[,'Ptgs2'] == 0] <- "Negative"
new_meta$expression.Ptges[new_meta[,'Ptges'] > 0] <- "Positive"
new_meta$expression.Ptges[new_meta[,'Ptges'] == 0] <- "Negative"
new_meta$expression.Ptger1[new_meta[,'Ptger1'] > 0] <- "Positive"
new_meta$expression.Ptger1[new_meta[,'Ptger1'] == 0] <- "Negative"
new_meta$expression.Ptger2[new_meta[,'Ptger2'] > 0] <- "Positive"
new_meta$expression.Ptger2[new_meta[,'Ptger2'] == 0] <- "Negative"
new_meta$expression.Ptger3[new_meta[,'Ptger3'] > 0] <- "Positive"
new_meta$expression.Ptger3[new_meta[,'Ptger3'] == 0] <- "Negative"
new_meta$expression.Ptger4[new_meta[,'Ptger4'] > 0] <- "Positive"
new_meta$expression.Ptger4[new_meta[,'Ptger4'] == 0] <- "Negative"
# new_meta$double_expression[new_meta$expression.Ptgs1 =="Positive" & new_meta$expression.Ptges =="Positive"] <- "Ptgs1_Ptges_double_positive"
# new_meta$double_expression[new_meta$expression.Ptgs2 =="Positive" & new_meta$expression.Ptges =="Positive"] <- "Ptgs2_Ptges_double_positive"

new_meta$Ptgs1_Ptges_double_positive <- "Negative"
new_meta$Ptgs1_Ptges_double_positive[new_meta$expression.Ptgs1 =="Positive" & new_meta$expression.Ptges =="Positive"] <- "Positive"
new_meta$Ptgs2_Ptges_double_positive <- "Negative"
new_meta$Ptgs2_Ptges_double_positive[new_meta$expression.Ptgs2 =="Positive" & new_meta$expression.Ptges =="Positive"] <- "Positive"

new_meta$triple_expression <- "Negative"
new_meta$triple_expression[new_meta$expression.Ptgs1 =="Positive" & new_meta$expression.Ptges =="Positive" & new_meta$expression.Ptgs2 =="Positive"] <- "Positive"

new_meta$receptor_expression <- "Negative"
new_meta$receptor_expression[new_meta$expression.Ptger1 =="Positive" | new_meta$expression.Ptger2 == "Positive" | new_meta$expression.Ptger3 == "Positive" | new_meta$expression.Ptger4 == "Positive"] <- "Positive"

genes <- colnames(new_meta[28:38])

pdf(paste0("./BI_out/COX_count_Epithelial.pdf"), width=20, height=10)
library(RColorBrewer)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library("dplyr")
library(ggpubr)
for (i in genes) {
  #play with x axis
  print(ggplot(new_meta, aes ( x = orig.ident, fill = new_meta[,i]))+
          #no subtype determined for this
          #print(ggplot(new_meta, aes ( x = SubCellType, fill = new_meta[,i]))+
          #or dodge
          geom_bar(stat= "count", position = "stack") +
          ylab("cell counts")+
          #scale_fill_brewer()+
          scale_fill_manual(values=c("#006ddb", "#FF0000")) +
          #size =15 for x= orig.ident
          theme_hc(base_size = 25)+
          #size 10 for x = orig.ident
          theme(axis.text.x = element_text(angle = 45, hjust = 1, face="bold", size=15),
                axis.text.y = element_text(face = "bold", color = "black", size = 15))+
          labs(fill=i))
  #facet_wrap("SubCellType"))
  
  
}

dev.off() 

#box for Treatment
tmp_new_meta = as.data.frame.array(table(new_meta$orig.ident, new_meta$double_expression,new_meta$Treatment))
tmp_new_meta_percentage= as.data.frame(apply(tmp_new_meta, 2, function(x) x*100/sum(x)))
tmp_new_meta_percentage$Sample <- rownames(tmp_new_meta_percentage)
tmp_new_meta_percentage = melt(tmp_new_meta_percentage)
colnames(tmp_new_meta_percentage)[2] <- "expressiontype"
tmp_new_meta_percentage$Treatment <- ""
tmp_new_meta_percentage$Treatment <- apply(tmp_new_meta_percentage, 1, function(x) str_split(x[2], "\\." ))
tmp_new_meta_percentage$Treatment <- apply(tmp_new_meta_percentage,1,function(x) x[[4]][[1]][2])
tmp_new_meta_percentage$expressiontype <- apply(tmp_new_meta_percentage, 1, function(x) str_split(x[2], "\\." ))
tmp_new_meta_percentage$expressiontype <- apply(tmp_new_meta_percentage,1,function(x) x[[2]][[1]][1])
tmp_new_meta_percentage <- subset(tmp_new_meta_percentage, expressiontype != "either or both Negative")
tmp_new_meta_percentage$CellLine <- ""
tmp_new_meta_percentage$CellLine[tmp_new_meta_percentage$Sample=="scGEX_152"] <- "TC3"
tmp_new_meta_percentage$CellLine[tmp_new_meta_percentage$Sample=="scGEX_153"] <- "TC3"
tmp_new_meta_percentage$CellLine[tmp_new_meta_percentage$Sample=="scGEX_154"] <- "TC3"
tmp_new_meta_percentage$CellLine[tmp_new_meta_percentage$Sample=="scGEX_155"] <- "TC3"
tmp_new_meta_percentage$CellLine[tmp_new_meta_percentage$Sample=="scGEX_156"] <- "TC3"
tmp_new_meta_percentage$CellLine[tmp_new_meta_percentage$Sample=="scGEX_157"] <- "TC3"
tmp_new_meta_percentage$CellLine[tmp_new_meta_percentage$Sample=="scGEX_158"] <- "TC3"
tmp_new_meta_percentage$CellLine[tmp_new_meta_percentage$Sample=="scGEX_159"] <- "TC4"
tmp_new_meta_percentage$CellLine[tmp_new_meta_percentage$Sample=="scGEX_160"] <- "TC4"
tmp_new_meta_percentage$CellLine[tmp_new_meta_percentage$Sample=="scGEX_161"] <- "TC4"
tmp_new_meta_percentage$CellLine[tmp_new_meta_percentage$Sample=="scGEX_162"] <- "TC4"
tmp_new_meta_percentage$CellLine[tmp_new_meta_percentage$Sample=="scGEX_163"] <- "TC4"
tmp_new_meta_percentage$CellLine[tmp_new_meta_percentage$Sample=="scGEX_164"] <- "TC4"
tmp_new_meta_percentage$CellLine[tmp_new_meta_percentage$Sample=="scGEX_165"] <- "TC4"
tmp_new_meta_percentage$CellLine[tmp_new_meta_percentage$Sample=="scGEX_166"] <- "TC5"
tmp_new_meta_percentage$CellLine[tmp_new_meta_percentage$Sample=="scGEX_167"] <- "TC5"
tmp_new_meta_percentage$CellLine[tmp_new_meta_percentage$Sample=="scGEX_168"] <- "TC5"
tmp_new_meta_percentage$CellLine[tmp_new_meta_percentage$Sample=="scGEX_169"] <- "TC5"
tmp_new_meta_percentage$CellLine[tmp_new_meta_percentage$Sample=="scGEX_170"] <- "TC5"
tmp_new_meta_percentage$CellLine[tmp_new_meta_percentage$Sample=="scGEX_171"] <- "TC5"
tmp_new_meta_percentage$CellLine[tmp_new_meta_percentage$Sample=="scGEX_172"] <- "TC5"
tmp_new_meta_percentage <- subset(tmp_new_meta_percentage, value > 0)

pdf(paste0("./BI_out/","COX_Epithelial_boxplot_Treatment_Cellline.pdf"),width = 20,height = 10)
Ptgs1 <- subset(tmp_new_meta_percentage, subset=expressiontype=="Ptgs1_Ptges_double_positive")
Ptgs2 <- subset(tmp_new_meta_percentage, subset=expressiontype=="Ptgs2_Ptges_double_positive")
library(RColorBrewer)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library("dplyr")
library(ggpubr)

compare_list = list(
  c('SOSMEK', 'Vehicle')
)

p1 <-ggplot(aes ( x = Treatment, y = value), data = Ptgs1) +
  geom_boxplot(inherit.aes = FALSE, outlier.shape =NA, aes(x = Treatment, y = value, fill= Treatment),outlier.color="black")+
  geom_point(color = 'darkblue')+
  ylab("Percentage of Ptgs1 and Ptges double positive cells")+
  facet_wrap( "CellLine" , scales = "fixed") + theme_hc() +
  stat_compare_means(comparisons = compare_list,method = 't.test')+
  theme(legend.position="none")+
  theme(text = element_text(size = 12),
        axis.text.x=element_text(size=10),
        axis.text.y=element_text(size=10),
        axis.title.y = element_text(size =16),
        axis.title.x = element_text(size=16))
print(p1)



p2 <-ggplot(aes ( x = Treatment, y = value), data = Ptgs2) +
  geom_boxplot(inherit.aes = FALSE, outlier.shape =NA, aes(x = Treatment, y = value, fill= Treatment),outlier.color="black")+
  geom_point(color = 'darkblue')+
  ylab("Percentage of Ptgs2 and Ptges double positive cells")+
  facet_wrap( "CellLine" , scales = "fixed") + theme_hc() +
  stat_compare_means(comparisons = compare_list,method = 't.test')+
  theme(legend.position="none")+
  theme(text = element_text(size = 12),
        axis.text.x=element_text(size=10),
        axis.text.y=element_text(size=10),
        axis.title.y = element_text(size =16),
        axis.title.x = element_text(size=16))
print(p2)

dev.off()

#box for cellline
tmp_new_meta = as.data.frame.array(table(new_meta$orig.ident, new_meta$double_expression,new_meta$CellLine))
tmp_new_meta_percentage= as.data.frame(apply(tmp_new_meta, 2, function(x) x*100/sum(x)))
tmp_new_meta_percentage$Sample <- rownames(tmp_new_meta_percentage)
tmp_new_meta_percentage = melt(tmp_new_meta_percentage)
colnames(tmp_new_meta_percentage)[2] <- "expressiontype"
tmp_new_meta_percentage$CellLine <- ""
tmp_new_meta_percentage$CellLine <- apply(tmp_new_meta_percentage, 1, function(x) str_split(x[2], "\\." ))
tmp_new_meta_percentage$CellLine <- apply(tmp_new_meta_percentage,1,function(x) x[[4]][[1]][2])
tmp_new_meta_percentage$expressiontype <- apply(tmp_new_meta_percentage, 1, function(x) str_split(x[2], "\\." ))
tmp_new_meta_percentage$expressiontype <- apply(tmp_new_meta_percentage,1,function(x) x[[2]][[1]][1])
tmp_new_meta_percentage <- subset(tmp_new_meta_percentage, expressiontype != "either or both Negative")
tmp_new_meta_percentage$Treatment <- ""
tmp_new_meta_percentage$Treatment[tmp_new_meta_percentage$Sample=="scGEX_152"] <- "Vehicle"
tmp_new_meta_percentage$Treatment[tmp_new_meta_percentage$Sample=="scGEX_153"] <- "Vehicle"
tmp_new_meta_percentage$Treatment[tmp_new_meta_percentage$Sample=="scGEX_154"] <- "Vehicle"
tmp_new_meta_percentage$Treatment[tmp_new_meta_percentage$Sample=="scGEX_155"] <- "SOSMEK"
tmp_new_meta_percentage$Treatment[tmp_new_meta_percentage$Sample=="scGEX_156"] <- "SOSMEK"
tmp_new_meta_percentage$Treatment[tmp_new_meta_percentage$Sample=="scGEX_157"] <- "SOSMEK"
tmp_new_meta_percentage$Treatment[tmp_new_meta_percentage$Sample=="scGEX_158"] <- "SOSMEK"
tmp_new_meta_percentage$Treatment[tmp_new_meta_percentage$Sample=="scGEX_159"] <- "Vehicle"
tmp_new_meta_percentage$Treatment[tmp_new_meta_percentage$Sample=="scGEX_160"] <- "Vehicle"
tmp_new_meta_percentage$Treatment[tmp_new_meta_percentage$Sample=="scGEX_161"] <- "Vehicle"
tmp_new_meta_percentage$Treatment[tmp_new_meta_percentage$Sample=="scGEX_162"] <- "SOSMEK"
tmp_new_meta_percentage$Treatment[tmp_new_meta_percentage$Sample=="scGEX_163"] <- "SOSMEK"
tmp_new_meta_percentage$Treatment[tmp_new_meta_percentage$Sample=="scGEX_164"] <- "SOSMEK"
tmp_new_meta_percentage$Treatment[tmp_new_meta_percentage$Sample=="scGEX_165"] <- "SOSMEK"
tmp_new_meta_percentage$Treatment[tmp_new_meta_percentage$Sample=="scGEX_166"] <- "Vehicle"
tmp_new_meta_percentage$Treatment[tmp_new_meta_percentage$Sample=="scGEX_167"] <- "Vehicle"
tmp_new_meta_percentage$Treatment[tmp_new_meta_percentage$Sample=="scGEX_168"] <- "Vehicle"
tmp_new_meta_percentage$Treatment[tmp_new_meta_percentage$Sample=="scGEX_169"] <- "SOSMEK"
tmp_new_meta_percentage$Treatment[tmp_new_meta_percentage$Sample=="scGEX_170"] <- "SOSMEK"
tmp_new_meta_percentage$Treatment[tmp_new_meta_percentage$Sample=="scGEX_171"] <- "SOSMEK"
tmp_new_meta_percentage$Treatment[tmp_new_meta_percentage$Sample=="scGEX_172"] <- "SOSMEK"
tmp_new_meta_percentage <- subset(tmp_new_meta_percentage, value > 0)

pdf(paste0("./BI_out/","COX_Epithelial_boxplot_Cellline_Treatment.pdf"),width = 20,height = 10)
Ptgs1 <- subset(tmp_new_meta_percentage, subset=expressiontype=="Ptgs1_Ptges_double_positive")
Ptgs2 <- subset(tmp_new_meta_percentage, subset=expressiontype=="Ptgs2_Ptges_double_positive")

compare_list = list(
  c('TC3','TC4'), c('TC3', 'TC5'), c('TC4', 'TC5'))

p1 <-ggplot(aes ( x = CellLine, y = value), data = Ptgs1) +
  geom_boxplot(inherit.aes = FALSE, outlier.shape =NA, aes(x = CellLine, y = value, fill=CellLine),outlier.color="black")+
  geom_point(color = 'darkblue')+
  ylab("Percentage of Ptgs1 and Ptges double positive cells")+
  facet_wrap( "Treatment", scales = "fixed") + theme_hc() +
  stat_compare_means(comparisons = compare_list,method = 't.test')+
  theme(legend.position="none")+
  theme(text = element_text(size = 12),
        axis.text.x=element_text(size=10),
        axis.text.y=element_text(size=10),
        axis.title.y = element_text(size =16),
        axis.title.x = element_text(size=16))
print(p1)



p2 <-ggplot(aes ( x = CellLine, y = value), data = Ptgs2) +
  geom_boxplot(inherit.aes = FALSE, outlier.shape =NA, aes(x = CellLine, y = value, fill=CellLine), outlier.color="black")+
  geom_point(color = 'darkblue')+
  ylab("Percentage of Ptgs2 and Ptges double positive cells")+
  facet_wrap( "Treatment" , scales = "fixed") + theme_hc() +
  stat_compare_means(comparisons = compare_list,method = 't.test')+
  theme(legend.position="none")+
  theme(text = element_text(size = 12),
        axis.text.x=element_text(size=10),
        axis.text.y=element_text(size=10),
        axis.title.y = element_text(size =16),
        axis.title.x = element_text(size=16))
print(p2)
dev.off()   




