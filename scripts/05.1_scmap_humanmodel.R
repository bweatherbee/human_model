library(Seurat)
library(Signac)
library(dyplr)
library(scmap)
library(biomaRt)
library(SingleCellExperiment)
library(RColorBrewer)

##### Preprocessing Nakamura 2016 data #####
setwd("D:/CynMonkey/Nakamura")
metadata_from_sups <- read_excel("metadata_from_sups.xlsx")
ReadsMapped<-as.data.frame(metadata_from_sups$`SC3seq mapped`)
rownames(ReadsMapped)<-metadata_from_sups$SampleID
colnames(ReadsMapped)<-'reads'

RPM<-read_tsv('SC3seq_Cy_ProcessedData.txt')
RPM<-as.data.frame(RPM)
gene_ids<-as.character(RPM$macFas5_entrez_id)
gene_names<-RPM$macFas5_gene_symbol
RPM<-RPM[,3:423]
common<-intersect(rownames(ReadsMapped),colnames(RPM))

ReadsMapped<-as.data.frame(ReadsMapped[common,])
rownames(ReadsMapped)<-common
colnames(ReadsMapped)<-'reads'

RPM<-as.matrix(RPM)

Counts<-RPM%*%diag(ReadsMapped$reads)
Counts<-Counts/10e6
colnames(Counts)<-common

rownames(RPM)<-as.character(gene_ids)

### subsetting for embryo cells in metadata and changing gene names to hgnc symbols ###
meta<-read.csv('metadata.csv')
meta<-meta[1:390,]
rownames(meta)<-meta$ï..SampleID
emb<-intersect(rownames(meta), colnames(Counts))
Counts<-Counts[,emb]
RPM<-RPM[,emb]

library(biomaRt)
mart <- useDataset("mfascicularis_gene_ensembl", useMart("ensembl"))

symbol <- getBM(filters = "ensembl_gene_id",
                attributes = c("external_gene_name","ensembl_gene_id","hgnc_symbol"),
                values = gene_ids, 
                mart = mart)
symbol<-symbol[!(symbol$hgnc_symbol==""),]

RPM<-RPM[symbol$ensembl_gene_id,]
rownames(RPM)<-make.unique(symbol$hgnc_symbol)

saveRDS(RPM, 'RPM_matrix_hgncsymbols.rds')

##### SCMAP ANALYSIS BASED ON http://bioconductor.org/packages/release/bioc/vignettes/scmap/inst/doc/scmap.html #####

### Creating SCE obj Object for Nakamura data ###
setwd("D:/human_model")
library(SingleCellExperiment)
library(scmap)
library(tidyverse)

rpm_matrix<- readRDS("D:/CynMonkey/Nakamura/RPM_matrix_hgncsymbols.rds")
metadata <- read_csv("D:/CynMonkey/Nakamura/metadata.csv")
metadata<-as.data.frame(metadata[1:390,])
rownames(metadata)<-metadata$SampleID

sce <- SingleCellExperiment(assays = list(normcounts = rpm_matrix), colData = metadata)
logcounts(sce) <- log2(normcounts(sce) + 1)
rowData(sce)$feature_symbol <- rownames(sce)
sce <- sce[!duplicated(rownames(sce)), ]
sce

human_model <- readRDS("D:/human_model/human_model_complete.rds")
DefaultAssay(human_model)<-'RNA'
human_model[['peaks']]<-NULL
human_model[['ATAC']]<-NULL
human_model[['chromvar']]<-NULL
human_model[['GeneActivity']]<-NULL
human_model[['SCT']]<-NULL
human_model<-as.SingleCellExperiment(human_model)
rowData(human_model)$feature_symbol <- rownames(human_model)


sce <- selectFeatures(sce, n_features=500, suppress_plot = FALSE)
table(rowData(sce)$scmap_features)
sce <- indexCluster(sce, cluster_col='cell_type')
head(metadata(sce)$scmap_cluster_index)
heatmap(as.matrix(metadata(sce)$scmap_cluster_index))


#Testing within Nakamura Dataset
scmapCluster_results <- scmapCluster(
  projection = human_model,
  threshold  = 0.5,
  index_list = list(
    nakamura = metadata(sce)$scmap_cluster_index
  )
)

plot<-
  getSankey(
    colData(human_model)$sample_type, 
    scmapCluster_results$scmap_cluster_labs[,'nakamura'],
    plot_height = 400,
    colors=c('#3944BC', '#710193', '#ED7014')
  )

print(plot, file='D:/human_model/scmap_humanmodel/gvis_sankey_nakamura_scmapclusters.html')


scmap_nakamura<-scmapCluster_results$combined_labs
scmap_nakamura<-as.data.frame(scmap_nakamura)
rownames(scmap_nakamura)<-human_model@colData@rownames

saveRDS(scmap_nakamura,'./scmap_humanmodel/scmap_nakamura_results.rds')
saveRDS(scmapCell_results, './scmap_humanmodel/scmap_nakamura_CELL.rds')
saveRDS(scmapCluster_results, './scmap_humanmodel/scmapnakamura_CLUSTER.rds')

##### Running with Yang et al 2021 Cynomolgus dataset (Using annotations from website and count table from Geo) #####
setwd("D:/human_model")
library(SingleCellExperiment)
library(scmap)
library(tidyverse)

counts<-read_tsv('./human_model_analysis/scmap_humanmodel/yang/CM_filtered_SoupX_corrected.tsv')
metadata<-readRDS('./human_model_analysis/scmap_humanmodel/yang/metadata_full.Rds')
metadata<-rbind(metadata[['d10']],metadata[['d12']],metadata[['d14']])
rownames(metadata)<-as.data.frame(metadata$cell)
genes<-counts$...1
counts<-counts[,rownames(metadata)]
rownames(counts)<-genes

library(biomaRt)
mart <- useDataset("mfascicularis_gene_ensembl", useMart("ensembl"))

symbol <- getBM(filters = "ensembl_gene_id",
                attributes = c("external_gene_name","ensembl_gene_id","hgnc_symbol"),
                values = genes, 
                mart = mart)
symbol<-symbol[!(symbol$hgnc_symbol==""),]

counts<-counts[symbol$ensembl_gene_id,]
rownames(counts)<-make.unique(symbol$hgnc_symbol)

Yang <- SingleCellExperiment(assays = list(normcounts = counts), colData = metadata)
logcounts(Yang) <- log2(normcounts(Yang) + 1)
rowData(Yang)$feature_symbol <- rownames(Yang)
Yang <- Yang[!duplicated(rownames(Yang)), ]
Yang

human_model <- readRDS("D:/human_model/human_model_complete.rds")
DefaultAssay(human_model)<-'RNA'
human_model[['peaks']]<-NULL
human_model[['ATAC']]<-NULL
human_model[['chromvar']]<-NULL
human_model[['SCT']]<-NULL
human_model<-as.SingleCellExperiment(human_model)
rowData(human_model)$feature_symbol <- rownames(human_model)

Yang <- selectFeatures(Yang, n_features=500, suppress_plot = FALSE)
table(rowData(Yang)$scmap_features)
Yang <- indexCluster(Yang, cluster_col='lineage')
head(metadata(Yang)$scmap_cluster_index)
heatmap(as.matrix(metadata(Yang)$scmap_cluster_index))


#Testing within Yang Dataset
set.seed(123)
Yang <- indexCell(Yang)
names(metadata(Yang)$scmap_cell_index)
length(metadata(Yang)$scmap_cell_index$subcentroids)
dim(metadata(Yang)$scmap_cell_index$subcentroids[[1]])
metadata(Yang)$scmap_cell_index$subcentroids[[1]][,1:5]
dim(metadata(Yang)$scmap_cell_index$subclusters)
metadata(Yang)$scmap_cell_index$subclusters[1:5,1:5]

scmapCell_results <- scmapCell(
  human_model,
  w = 10,
  list(
    Yang = metadata(Yang)$scmap_cell_index
  )
)

names(scmapCell_results)
scmapCell_results$Yang$cells[,1:3]
scmapCell_results$Yang$similarities[,1:3]

scmapCell_clusters <- scmapCell2Cluster(
  scmapCell_results,
  w=2,
  threshold=0.25,
  list(
    as.character(colData(Yang)$lineage)
  )
)
head(scmapCell_clusters$scmap_cluster_labs)
head(scmapCell_clusters$scmap_cluster_siml)

plot<-
  getSankey(
    colData(human_model)$sample_type, 
    scmapCell_clusters$scmap_cluster_labs[,"Yang"],
    plot_height = 400,
    colors=c('#3944BC', '#710193', '#ED7014')
  )

print(plot, file='D:/human_model/human_model_analysis/scmap_humanmodel/gvis_sankey_yang_scmapCELL.html')

scmapcell<-scmapCell_clusters$combined_labs
scmap_Yang$scmapCELL_Yang<-as.character(scmapcell)

saveRDS(scmapCell_results, './human_model_analysis/scmap_humanmodel/scmap_yang_CELL.rds')

##### Preprocessing Ma 2019 dataset based on deposited count table and metadata#####
library(readxl)

counts<-read_csv("scmap_humanmodel/ma/GSE130114_MF1453.csv")
counts<-as.data.frame(counts)
rownames(counts)<-counts$...1
counts<-counts[,2:1454]
saveRDS(counts, './scmap_humanmodel/ma/counts.rds')

library(biomaRt)
mart <- useDataset("mfascicularis_gene_ensembl", useMart("ensembl"))

symbol <- getBM(filters = "external_gene_name",
                attributes = c("external_gene_name","ensembl_gene_id","hgnc_symbol"),
                values = rownames(counts), 
                mart = mart)
symbol<-symbol[!(symbol$hgnc_symbol==""),]

counts<-counts[symbol$external_gene_name,]
rownames(counts)<-make.unique(symbol$hgnc_symbol)

saveRDS(counts, "D:/human_model/scmap_humanmodel/ma/counts_hgncsymbols.rds")

##### Creating SCE obj Object for MA data and scmap analysis #####
setwd("D:/human_model")
library(SingleCellExperiment)
library(scmap)
library(tidyverse)
meta <- read_csv("scmap_humanmodel/ma/GSE130114_MF1453-meta.csv")
meta <- as.data.frame(meta)
rownames(meta)<-meta$...1
meta <- meta[,2:3]
counts <- readRDS("D:/human_model/scmap_humanmodel/ma/counts_hgncsymbols.rds")

sce <- SingleCellExperiment(assays = list(normcounts = counts), colData = meta)
logcounts(sce) <- log2(normcounts(sce) + 1)
rowData(sce)$feature_symbol <- rownames(sce)
sce <- sce[!duplicated(rownames(sce)), ]
sce

human_model <- readRDS("D:/human_model/human_model_complete.rds")
DefaultAssay(human_model)<-'RNA'
human_model[['peaks']]<-NULL
human_model[['ATAC']]<-NULL
human_model[['chromvar']]<-NULL
human_model[['SCT']]<-NULL
human_model<-as.SingleCellExperiment(human_model)
rowData(human_model)$feature_symbol <- rownames(human_model)


sce <- selectFeatures(sce, n_features=500, suppress_plot = FALSE)
table(rowData(sce)$scmap_features)
sce <- indexCluster(sce, cluster_col='annotation')
head(metadata(sce)$scmap_cluster_index)
heatmap(as.matrix(metadata(sce)$scmap_cluster_index))


#Testing within Nakamura Dataset
scmapCluster_results <- scmapCluster(
  projection = human_model,
  threshold  = 0.5,
  index_list = list(
    nakamura = metadata(sce)$scmap_cluster_index
  )
)

plot<-(
  getSankey(
    colData(human_model)$sample_type, 
    scmapCluster_results$scmap_cluster_labs[,'nakamura'],
    plot_height = 400,
    colors=c('#3944BC', '#710193', '#ED7014')
  )
)

print(plot, file='D:/human_model/scmap_humanmodel/gvis_sankey_ma_scmapcluster.html')


scmap_ma<-scmapCluster_results$combined_labs
scmap_ma<-as.data.frame(scmap_ma)
rownames(scmap_ma)<-human_model@colData@rownames

saveRDS(scmap_ma,'./scmap_humanmodel/scmap_ma_results.rds')
saveRDS(scmapCluster_results, './scmap_humanmodel/scmap_CLUSTER_ma.rds')


##### Creating SCE for Tyser 2021 data and analysis #####
expression_values <- readRDS("D:/human_model/human_model_analysis/scmap_humanmodel/tyser/expression_values.rds")
annot_umap <- readRDS("D:/human_model/human_model_analysis/scmap_humanmodel/tyser/annot_umap.rds")
rownames(annot_umap)<-annot_umap$cell_name
rownames(expression_values)<-annot_umap$cell_name

Tyser <- SingleCellExperiment(assays = list(normcounts = t(expression_values)), colData = annot_umap)
logcounts(Tyser) <- log2(normcounts(Tyser) + 1)
rowData(Tyser)$feature_symbol <- rownames(Tyser)
Tyser <- Tyser[!duplicated(rownames(Tyser)), ]
Tyser

human_model <- readRDS("D:/human_model/human_model_complete.rds")
DefaultAssay(human_model)<-'RNA'
human_model[['peaks']]<-NULL
human_model[['ATAC']]<-NULL
human_model[['chromvar']]<-NULL
human_model[['GeneActivity']]<-NULL
human_model[['SCT']]<-NULL
human_model<-as.SingleCellExperiment(human_model)
rowData(human_model)$feature_symbol <- rownames(human_model)

Tyser <- selectFeatures(Tyser, n_features=500, suppress_plot = FALSE)
table(rowData(Tyser)$scmap_features)
Tyser <- indexCluster(Tyser, cluster_col='sub_cluster')
head(metadata(Tyser)$scmap_cluster_index)
heatmap(as.matrix(metadata(Tyser)$scmap_cluster_index))


#Testing within Tyser Dataset
scmapCluster_results <- scmapCluster(
  projection = human_model,
  threshold  = 0.5,
  index_list = list(
    Tyser = metadata(Tyser)$scmap_cluster_index
  )
)

plot<-
  getSankey(
    colData(human_model)$sample_type, 
    scmapCluster_results$scmap_cluster_labs[,'Tyser'],
    plot_height = 400,
    colors=c('#3944BC', '#710193', '#ED7014')
  )

print(plot, file='D:/human_model/human_model_analysis/scmap_humanmodel/gvis_sankey_Tyser_scmapclusters.html')

saveRDS(scmapCluster_results, './human_model_analysis/scmap_humanmodel/scmap_cluster_results_Tyser.rds')

scmap_Tyser<-scmapCluster_results$combined_labs
scmap_Tyser<-as.data.frame(scmap_Tyser)
rownames(scmap_Tyser)<-human_model@colData@rownames

saveRDS(scmap_Tyser,'./human_model_analysis/scmap_humanmodel/scmap_Tyser_results.rds')

##### Creating SCE for Mole 2021 data and analysis #####
MZG <- readRDS("D:/realigning_for_int/MZG/embryos_integrated.RDS")
Mole <- as.SingleCellExperiment(MZG)
rowData(Mole)$feature_symbol <- rownames(Mole)
Mole <- Mole[!duplicated(rownames(Mole)), ]
Mole

human_model <- readRDS("D:/human_model/human_model_complete.rds")
DefaultAssay(human_model)<-'RNA'
human_model[['peaks']]<-NULL
human_model[['ATAC']]<-NULL
human_model[['chromvar']]<-NULL
human_model[['SCT']]<-NULL
human_model<-as.SingleCellExperiment(human_model)
rowData(human_model)$feature_symbol <- rownames(human_model)

Mole <- selectFeatures(Mole, n_features=500, suppress_plot = FALSE)
table(rowData(Mole)$scmap_features)
Mole <- indexCluster(Mole, cluster_col='cell_type')
head(metadata(Mole)$scmap_cluster_index)
heatmap(as.matrix(metadata(Mole)$scmap_cluster_index))


#Testing within Mole Dataset
set.seed(123)
Mole <- indexCell(Mole)
names(metadata(Mole)$scmap_cell_index)
length(metadata(Mole)$scmap_cell_index$subcentroids)
dim(metadata(Mole)$scmap_cell_index$subcentroids[[1]])
metadata(Mole)$scmap_cell_index$subcentroids[[1]][,1:5]
dim(metadata(Mole)$scmap_cell_index$subclusters)
metadata(Mole)$scmap_cell_index$subclusters[1:5,1:5]

scmapCell_results <- scmapCell(
  human_model,
  w = 10,
  list(
    Mole = metadata(Mole)$scmap_cell_index
  )
)

names(scmapCell_results)
scmapCell_results$Mole$cells[,1:3]
scmapCell_results$Mole$similarities[,1:3]

scmapCell_clusters <- scmapCell2Cluster(
  scmapCell_results,
  w=2,
  threshold=0.25,
  list(
    as.character(colData(Mole)$cell_type)
  )
)
head(scmapCell_clusters$scmap_cluster_labs)
head(scmapCell_clusters$scmap_cluster_siml)

plot<-
  getSankey(
    colData(human_model)$sample_type, 
    scmapCell_clusters$scmap_cluster_labs[,"Mole"],
    plot_height = 400,
    colors=c('#3944BC', '#710193', '#ED7014')
  )

print(plot, file='D:/human_model/human_model_analysis/scmap_humanmodel/gvis_sankey_Mole_scmapCELL.html')

scmapcell<-scmapCell_clusters$combined_labs
scmap_Mole$scmapCELL_Mole<-as.character(scmapcell)

saveRDS(scmap_Mole,'./human_model_analysis/scmap_humanmodel/scmap_Mole_results.rds')

##### Adding SCMAP results as metadata for vis and downstream filtering #####
human_model <- readRDS("D:/human_model/human_model_complete.rds")
scmap_nakamura <- readRDS("D:/human_model/human_model_analysis/scmap_humanmodel/scmap_nakamura_results.rds")
scmap_yang <- readRDS("D:/human_model/human_model_analysis/scmap_humanmodel/scmap_yang_results.rds")
scmap_ma <- readRDS("D:/human_model/human_model_analysis/scmap_humanmodel/scmap_ma_results.rds")
scmap_tyser<-readRDS("D:/human_model/human_model_analysis/scmap_humanmodel/scmap_Tyser_results.rds")
scmap_mole<-readRDS("D:/human_model/human_model_analysis/scmap_humanmodel/scmap_Mole_results.rds")

human_model<-AddMetaData(human_model, metadata=scmap_nakamura)
human_model<-AddMetaData(human_model, metadata=scmap_yang)
human_model<-AddMetaData(human_model, metadata=scmap_ma)
human_model<-AddMetaData(human_model, metadata=scmap_tyser)
human_model<-AddMetaData(human_model, metadata=scmap_mole)

DimPlot(human_model, reduction='wnn.umap', group.by=c('scmap_nakamura', 
                                                      'scmapCELL_Yang', 
                                                      'scmap_ma',
                                                      'scmap_Tyser',
                                                      'scmapCELL_Mole'))

saveRDS(human_model, 'human_model.rds')

DimPlot(human_model, reduction='wnn.umap', order=T, sizes.highlight = 2, 
        cells.highlight=
          rownames(subset(human_model@meta.data,scmap_Tyser=="Primitive Streak")))

DefaultAssay(human_model)<-'RNA'
VlnPlot(human_model, features=c('TFAP2C', 'NANOG', 'SOX17',
                                'NANOS3', 'PRDM1', 'GFP'), group.by='scmap_ma')

#Based on these results, will subcluster the current cluster '1' further and assign names to model for major cell types
DimPlot(human_model, reduction='wnn.umap')
human_model<-FindSubCluster(
  human_model,
  '1',
  'wsnn',
  subcluster.name = "sub.cluster",
  resolution = 0.25,
  algorithm = 2
)

DimPlot(human_model, reduction='wnn.umap', group.by='sub.cluster')

Idents(human_model)<-human_model$sub.cluster
new.cluster.ids<-c('MESO-2', 'L-EPI', 'AM-1', 'HYPO/VE', 'AM-2', 'MESO-1', 'AM-3', 'EXMC')
names(new.cluster.ids)<-levels(human_model)
human_model<-RenameIdents(human_model, new.cluster.ids)
DimPlot(human_model, reduction='wnn.umap')
human_model$cell_assignment<-Idents(human_model)

saveRDS(human_model, 'human_model.rds')