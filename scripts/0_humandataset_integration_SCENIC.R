library(Seurat)
library(data.table)
library(tidyverse)

# realignment information and qUMI normalizations can be found at github.com/bweatherbee/PeriImplantation

Petropoulos <- readRDS("/Petropoulos/Petropoulos_qUMI.RDS")
MZG <- readRDS("/MZG/MZG.RDS")
MZG$percent.mt<-MZG$percent.mito
MZG$species<-"human"
MZG$paper<-"Mole"
DefaultAssay(MZG)<-"RNA"
MZG$percent.mt<-MZG$percent.mito
MZG$percent.mito<-NULL
MZG$nCount_SCENIC<-NULL
MZG$nFeature_SCENIC<-NULL
MZG@assays$SCENIC<-NULL
MZG@assays$dorothea<-NULL
MZG<-NormalizeData(MZG, normalization.method="RC", scale.factor = 10e6)

Yan <- readRDS("/Yan_qUMI.RDS")
Blakely <- readRDS("/Blakely_qUMI.RDS")
Zhou <- readRDS("/Zhou_UMIcounts.rds")
Zhou$Age<-Zhou$timepoint
Zhou$timepoint<-NULL
Zhou<-NormalizeData(Zhou, normalization.method = "RC", scale.factor = 10e6)

Xiang <- readRDS("/Xiang/untrimmed/kallisto_out/gene_matrix_estimated_counts_output/August2021/Xiang_qUMI.RDS")


Blakely_matrix<-Blakely@assays[["RNA_qumi"]]@counts
Blakely_meta<-Blakely@meta.data
MZG_matrix<-MZG@assays[["RNA"]]@counts
MZG_meta<-MZG@meta.data
Petropoulos_matrix<-Petropoulos@assays[["RNA_qumi"]]@counts
Petropoulos_meta<-Petropoulos@meta.data
Xiang_matrix<-Xiang@assays[["RNA_qumi"]]@counts
Xiang_meta<-Xiang@meta.data
Yan_matrix<-Yan@assays[["RNA_qumi"]]@counts
Yan_meta<-Yan@meta.data
Zhou_matrix<-Zhou@assays[["RNA"]]@counts
Zhou_meta<-Zhou@meta.data

common.features<-intersect(rownames(Blakely_matrix), rownames(MZG_matrix))
common.features<-intersect(rownames(Zhou_matrix), common.features)
Blakely_matrix<-Blakely_matrix[common.features,]
MZG_matrix<-MZG_matrix[common.features,]
Petropoulos_matrix<-Petropoulos_matrix[common.features,]
Xiang_matrix<-Xiang_matrix[common.features,]
Yan_matrix<-Yan_matrix[common.features,]
Zhou_matrix<-Zhou_matrix[common.features,]

merged_matrix<-merge(x=Blakely_matrix, y=c(MZG_matrix, Petropoulos_matrix, Xiang_matrix, Yan_matrix, Zhou_matrix),by=0)
rownames(merged_matrix)<-merged_matrix$Row.names
merged_matrix<-merged_matrix[,2:10225]

common.meta<-intersect(colnames(Blakely_meta), colnames(MZG_meta))
Blakely_meta<-Blakely_meta[,common.meta]
MZG_meta<-MZG_meta[,common.meta]
Petropoulos_meta$cell_type<-"TBC"
Petropoulos_meta<-Petropoulos_meta[,common.meta]
Xiang_meta<-Xiang_meta[,common.meta]
Yan_meta$cell_type<-"TBC"
Yan_meta<-Yan_meta[,common.meta]
Zhou_meta$percent.mt<-"TBC"
Zhou_meta$seurat_clusters<-"TBC"
Zhou_meta<-Zhou_meta[,common.meta]

l=list(Blakely_meta, MZG_meta, Petropoulos_meta, Xiang_meta, Yan_meta, Zhou_meta)
merged_meta<-rbindlist(l, use.names=TRUE)
rownames(merged_meta)<-colnames(merged_matrix)

merged<-CreateSeuratObject(counts=merged_matrix, project="Human_merged", assay="UMI_qumi",
                           meta.data=merged_meta)

merged[['percent.mt']]<-PercentageFeatureSet(merged, pattern="^MT")
VlnPlot(merged, features=c('nFeature_UMI_qumi', 'nCount_UMI_qumi', 'percent.mt'))

merged<-NormalizeData(merged, normalization.method = "RC", scale.factor = 10e6)

merged@meta.data<-merged@meta.data %>% 
  mutate(Age = str_replace(Age, "D1", "D01"))
merged@meta.data<-merged@meta.data %>% 
  mutate(Age = str_replace(Age, "D2", "D02"))
merged@meta.data<-merged@meta.data %>% 
  mutate(Age = str_replace(Age, "D3", "D03"))
merged@meta.data<-merged@meta.data %>% 
  mutate(Age = str_replace(Age, "D4", "D04"))
merged@meta.data<-merged@meta.data %>% 
  mutate(Age = str_replace(Age, "D5", "D05"))
merged@meta.data<-merged@meta.data %>% 
  mutate(Age = str_replace(Age, "D6", "D06"))
merged@meta.data<-merged@meta.data %>% 
  mutate(Age = str_replace(Age, "D7", "D07"))
merged@meta.data<-merged@meta.data %>% 
  mutate(Age = str_replace(Age, "D010", "D10"))
merged@meta.data<-merged@meta.data %>% 
  mutate(Age = str_replace(Age, "D011", "D11"))
merged@meta.data<-merged@meta.data %>% 
  mutate(Age = str_replace(Age, "D012", "D12"))
merged@meta.data<-merged@meta.data %>% 
  mutate(Age = str_replace(Age, "D014", "D14"))

#Realized there is a repeat from the end of Petropoulos, removing
cells<-colnames(merged)
removeWords <- function(str, stopwords) {
  x <- unlist(strsplit(str, " "))
  paste(x[!x %in% stopwords])
}
cells<-removeWords(cells, 'E3.53.3438.1')

merged<-subset(merged, cells=cells)

#### cycle and mitocondrial content regressions + Merging ###
DefaultAssay(merged)<-"UMI_qumi"
s.genes<-cc.genes$s.genes
g2m.genes<-cc.genes$g2m.genes
merged<-CellCycleScoring(merged, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
RidgePlot(merged, features=c("PCNA", "TOP2A", "MCM6", "MKI67", ncol=2))
VlnPlot(merged, features="percent.mt")
merged<-ScaleData(merged, vars.to.regress = c("S.Score", "G2M.Score", "percent.mt"))


merged_list<-SplitObject(merged, split.by="paper")
merged_list<-lapply(X=merged_list, FUN=function(x) {
  x<-NormalizeData(x, normalization.method = 'RC', scale.factor = 10e6)
  x<-FindVariableFeatures(x, selection.method = 'vst', nfeatures=2000)
})

features<-SelectIntegrationFeatures(object.list=merged_list)
anchors<-FindIntegrationAnchors(object.list=merged_list, anchor.features=features)
merged<-IntegrateData(anchorset=anchors,k.weight=80)
saveRDS(merged, 'merged.RDS')

DefaultAssay(merged)<-"integrated"
merged<-ScaleData(merged, vars.to.regress=c("S.Score", "G2M.Score", "percent.mt"),features=rownames(merged))
merged<-RunPCA(merged, features=rownames(merged))
merged<-RunUMAP(merged, dims=1:10)
ElbowPlot(merged, ndims = 50)
merged<-FindNeighbors(merged, reduction="pca", dims=1:50)
merged<-FindClusters(merged, resolution = 0.2)
DimPlot(merged)

DimPlot(merged, group.by = 'paper')
DimPlot(merged, group.by='Age')
DimPlot(merged)

#DefaultAssay(merged)<-"UMI_qumi"
FeaturePlot(merged, features=c("SOX2", "NANOG", "POU5F1", "KLF17"), order=TRUE, slot="data", min.cutoff = 'q0', max.cutoff = 'q98')
FeaturePlot(merged, features=c("SOX17", "GATA4", "GATA6", "PDGFRA"), order=TRUE, slot="data", min.cutoff = 'q0', max.cutoff = 'q98')
FeaturePlot(merged, features=c("GATA3", "GATA2","TFAP2C", "CDX2"),order=TRUE, slot="data", min.cutoff = 'q0', max.cutoff = 'q98')
FeaturePlot(merged, features=c("SDC1", "CGB1", "CGA", "ITGB3"), order=TRUE, slot="data", min.cutoff = 'q0', max.cutoff = 'q98')
DimPlot(merged, cells.highlight = (colnames(subset(merged, subset=Age=="D0"))))


new.cluster.ids<-c("TE", "STB", "TE", "TE", "EVT", "EPI", "PreLin", "STB", "PreLin", "HYPO")
names(new.cluster.ids)<-levels(merged)
merged<-RenameIdents(merged, new.cluster.ids)
DimPlot(merged)

DimPlot(merged, cells.highlight = (colnames(subset(merged, idents="PreLin"))))

merged$new_cell_type<-Idents(merged)

metadata<-merged@meta.data


df<-metadata$Age
df<-as.data.frame(df)
df$new_cell_type<-metadata$new_cell_type
rownames(df)<-rownames(metadata)
colnames(df)<-c("Age", "new_cell_type")


df <- within(df, {
  f <- Age == 'D0'
  Age[f] <- 'D0'
  new_cell_type[f] <- 'PreLin'
}) 

df <- within(df, {
  f <- Age == 'D01'
  Age[f] <- 'D01'
  new_cell_type[f] <- 'PreLin'
}) 

df <- within(df, {
  f <- Age == 'D02'
  Age[f] <- 'D02'
  new_cell_type[f] <- 'PreLin'
}) 

df <- within(df, {
  f <- Age == 'D03'
  Age[f] <- 'D03'
  new_cell_type[f] <- 'PreLin'
}) 

df$cell_type<-metadata$cell_type


df <- within(df, {
  f <- Age == 'D06' & new_cell_type=="PreLin" & cell_type!="TBC"
  Age[f] <- 'D06'
  new_cell_type[f] <- cell_type[f]
}) 

df['D14A1B16','new_cell_type']<-"EPI"
df['D14A1B7','new_cell_type']<-"EPI"
df['D14A1S21','new_cell_type']<-"EPI"
df['D14A1S26','new_cell_type']<-"EPI"
df['hv_D12_IVC7_E2_B1_11','new_cell_type']<-"HYPO"
df['hv_D12_IVC7_E1_B1_4','new_cell_type']<-"TE"
df['hv_D12_IVC7_E1_B1_7','new_cell_type']<-"TE"
df['hv_D12_IVC7_E1_B1_10','new_cell_type']<-"TE"
df['hv_D12_IVC7_E1_B1_13','new_cell_type']<-"TE"
df['hv_D12_IVC7_E1_B1_17','new_cell_type']<-"TE"
df['hv_D12_IVC7_E1_B1_21','new_cell_type']<-"TE"
df['hv_D12_IVC7_E1_B1_23','new_cell_type']<-"TE"
df['hv_D12_IVC7_E1_B1_28','new_cell_type']<-"TE"
df['hv_D12_IVC7_E1_B1_43','new_cell_type']<-"TE"
df['embryo3_AAGTCTGTCGGAAATA.1','new_cell_type']<-"EPI"
df['embryo3_CAACCAATCCAAGTAC.1','new_cell_type']<-"EPI"
df['D10_IVC4_E2_B3_81','new_cell_type']<-"HYPO"
df['D10A4B3','new_cell_type']<-"TE"
df['hm_D10_IVC5_E9_B3_57','new_cell_type']<-"TE"
df['embryo7_TCCCGATAGACGCAAC.1','new_cell_type']<-"EPI"
df['D8_8S3','new_cell_type']<-"HYPO"
df['hm_D8_IVC2_E2_B1_91','new_cell_type']<-"HYPO"
df['hm_D8_IVC2_E3_B2_95','new_cell_type']<-"HYPO"
df['D6A2S6','new_cell_type']<-"TE"
df['D6A3S4','new_cell_type']<-"TE"
df['D6N2B28','new_cell_type']<-"TE"
df['hv_D8_IVC2_E1_B1_73','new_cell_type']<-"HYPO"


df <- within(df, {
  f <- new_cell_type=="PreLin" & cell_type=="MIX"
  new_cell_type[f] <- "TE"
})

df <- within(df, {
  f <- Age == 'D14' & new_cell_type=="PreLin" & cell_type!="TBC"
  Age[f] <- 'D14'
  new_cell_type[f] <- cell_type[f]
})

df <- within(df, {
  f <- Age == 'D10' & new_cell_type=="PreLin" & cell_type!="TBC"
  Age[f] <- 'D10'
  new_cell_type[f] <- cell_type[f]
})

df <- within(df, {
  f <- Age == 'D12' & new_cell_type=="PreLin" & cell_type!="TBC"
  Age[f] <- 'D12'
  new_cell_type[f] <- cell_type[f]
})

df <- within(df, {
  f <- Age == 'D09' & new_cell_type=="PreLin" & cell_type!="TBC"
  Age[f] <- 'D09'
  new_cell_type[f] <- cell_type[f]
})

df <- within(df, {
  f <- Age == 'D08' & new_cell_type=="PreLin" & cell_type!="TBC"
  Age[f] <- 'D08'
  new_cell_type[f] <- cell_type[f]
})


df <- within(df, {
  f <- Age == 'D08' & new_cell_type=="EVT"
  Age[f] <- 'D08'
  new_cell_type[f] <- "TE"
})

df <- within(df, {
  f <- Age == 'D07' & new_cell_type=="EVT"
  Age[f] <- 'D07'
  new_cell_type[f] <- "TE"
})

df <- within(df, {
  f <- Age == 'D06' & new_cell_type=="EVT"
  Age[f] <- 'D06'
  new_cell_type[f] <- "TE"
})

df <- within(df, {
  f <- Age == 'D05' & new_cell_type=="EVT"
  Age[f] <- 'D05'
  new_cell_type[f] <- "TE"
})


df <- within(df, {
  f <- Age == 'D08' & new_cell_type=="STB"
  Age[f] <- 'D08'
  new_cell_type[f] <- "TE"
})

df <- within(df, {
  f <- Age == 'D07' & new_cell_type=="STB"
  Age[f] <- 'D07'
  new_cell_type[f] <- "TE"
})

df <- within(df, {
  f <- Age == 'D06' & new_cell_type=="STB"
  Age[f] <- 'D06'
  new_cell_type[f] <- "TE"
})

df <- within(df, {
  f <- Age == 'D05' & new_cell_type=="STB"
  Age[f] <- 'D05'
  new_cell_type[f] <- "TE"
})


df <- within(df, {
  f <- Age == 'D07' & new_cell_type=="PreLin"
  Age[f] <- 'D07'
  new_cell_type[f] <- "TE"
}) 

df <- within(df, {
  f <- Age == 'D06' & new_cell_type=="PreLin"
  Age[f] <- 'D06'
  new_cell_type[f] <- "TE"
}) 

merged$new_cell_type<-df$new_cell_type

DimPlot(merged, group.by = 'new_cell_type')
Idents(merged)<-merged$new_cell_type


levels(merged$new_cell_type)
my_levels<-c("PreLin", "EPI", "HYPO", "TE", "STB", "EVT")
merged$new_cell_type <- factor(x = merged$new_cell_type, levels = my_levels)

merged<-saveRDS(merged, 'merged_integrated.RDS')



library(Seurat)
library(tidyverse)

merged <- readRDS("D:/realigning_for_int/Combined_blast_and_Post/SCENIC/results/merged_integrated_withSCENIC.RDS")
DefaultAssay(merged)<-'SCENIC'

HYPO_regsvTE<-FindMarkers(merged, ident.1="HYPO", ident.2="TE", only.pos = T, logfc.threshold = 0.025, min.pct = 0.95, test.use = 'wilcox')
HYPO_regsvEPI<-FindMarkers(merged, ident.1="HYPO", ident.2="EPI", only.pos = T, logfc.threshold = 0.025, min.pct=0.95, test.use = 'wilcox')
HYPO_regsvTE<-subset(HYPO_regsvTE, subset=p_val_adj<0.05)
HYPO_regsvEPI<-subset(HYPO_regsvEPI, subset=p_val_adj<0.05)
Hypo_regs<-intersect(rownames(HYPO_regsvTE), rownames(HYPO_regsvEPI))

EPI_regsvTE<-FindMarkers(merged, ident.1="EPI", ident.2="TE", only.pos=T, logfc.threshold = 0.025, min.pct = 0.95)
EPI_regsvHYPO<-FindMarkers(merged, ident.1="EPI", ident.2="HYPO", only.pos=T, logfc.threshold = 0.025, min.pct=0.95)
EPI_regsvTE<-subset(EPI_regsvTE, subset=p_val_adj<0.05)
EPI_regsvHYPO<-subset(EPI_regsvHYPO, subset=p_val_adj<0.05)
EPI_regs<-intersect(rownames(EPI_regsvTE), rownames(EPI_regsvHYPO))

TE_regsvEPI<-FindMarkers(merged, ident.1="TE", ident.2="EPI", only.pos=T, logfc.threshold = 0.025, min.pct=0.95)
TE_regsvHYPO<-FindMarkers(merged, ident.1="TE", ident.2="HYPO", only.pos=T, logfc.threshold = 0.025, min.pct=0.95)
TE_regsvEPI<-subset(TE_regsvEPI, subset=p_val_adj<0.05)
TE_regsvHYPO<-subset(TE_regsvHYPO, subset=p_val_adj<0.05)
TE_regs<-intersect(rownames(TE_regsvEPI), rownames(TE_regsvHYPO))

list<-list('EPI'=EPI_regs, 'HYPO'=Hypo_regs, 'TE'=TE_regs)
capture.output(list, file = "SCENIC_pairwise_enriched_regs.txt")
regulonTargetsInfo <- readRDS("C:/Users/baile/Desktop/realigning_for_int/Combined_blast_and_Post/SCENIC/results/int/2.5_regulonTargetsInfo.Rds")
write.csv(regulonTargetsInfo, 'regulonTargetsInfo.csv')
key_genes<-c(EPI_regs, Hypo_regs, TE_regs)
key_genes_2<-gsub("-.*","",key_genes)
subset_target_table<-subset(regulonTargetsInfo, subset=(regulonTargetsInfo$gene %in% key_genes_2 & regulonTargetsInfo$TF %in% key_genes_2))
subset_target_table<-subset(subset_target_table, subset=(subset_target_table$TF != subset_target_table$gene))
subset_target_table$EPI<-subset_target_table$TF%in%gsub("-.*", "", EPI_regs)
subset_target_table$HYPO<-subset_target_table$TF%in%gsub("-.*", "", Hypo_regs)
subset_target_table$TE<-subset_target_table$TF%in%gsub("-.*", "", TE_regs)
write.csv(subset_target_table, 'trimmed_network_table.csv')

#Trimmed network table was used to build network in Cytoscape.