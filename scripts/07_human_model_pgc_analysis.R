library(Seurat)
library(Signac)
library(RColorBrewer)
library(Scillus)
library(tidyverse)

setwd("D://human_model")
human_model <- readRDS("D:/human_model/human_model.rds")
DefaultAssay(human_model)<-'RNA'

#Using PGC genes identified from Jo et al 2022, elife Efficient Differentiation of human PGCs in micopatterns
PGC_genes<-list(c('SOX17','PRDM1','NANOS3','CRHBP','TFAP2C',
             'AC007326.5','LBH','GLIPR2','TMEM64','PCAT14',
             'MPC2','CA2','PLPP1','FAM184A','NTF3','AL353807.5','PRTG','DBI'))

human_model<-AddModuleScore(human_model, features=PGC_genes, name='PGC_Module')


FeaturePlot(human_model, reduction='wnn.umap', features='PGC_Module1', 
            min.cutoff = 0.5, max.cutoff = 0.85, order=T, pt.size = 1,
            cols=c('light grey','dark green'))

metadata<-human_model@meta.data
df<-metadata$cell_assignment
df<-as.data.frame(df)
df$PGC_Module1<-metadata$PGC_Module1
rownames(df)<-rownames(metadata)
colnames(df)<-c("cell_assignment", "PGC_Module1")
df$cell_assignment<-as.character(df$cell_assignment)

#Going with top 2% of scores (0.50 module score)
df <- within(df, {
  f <- PGC_Module1 >0.50
  cell_assignment[f] <- 'PGC'
  PGC_Module1[f] <- PGC_Module1[f]
}) 

human_model$cell_assignment2<-df$cell_assignment

Idents(human_model)<-human_model$cell_assignment2
Idents(human_model)

#ordering assigned cell types for plotting
Idents(human_model)<-factor(x=Idents(human_model),
                            levels=c('L-EPI', 'AM-1', 'AM-2', 'AM-3', 'MESO-1', 'MESO-2', 'EXMC', 'HYPO/VE', 'PGC'))
human_model$cell_assignment2<-Idents(human_model)

Idents(human_model)<-human_model$cell_assignment
Idents(human_model)<-factor(x=Idents(human_model),
                            levels=c('L-EPI', 'AM-1', 'AM-2', 'AM-3', 'MESO-1', 'MESO-2', 'EXMC', 'HYPO/VE'))
human_model$cell_assignment<-Idents(human_model)

Idents(human_model)<-human_model$course_cell_assignment
Idents(human_model)<-factor(x=Idents(human_model),
                            levels=c('L-EPI', 'AME', 'MESO', 'EXMC', 'HYPO/VE'))
human_model$course_cell_assignment<-Idents(human_model)

#plotting key PGCLC gene expression
DimPlot(human_model, reduction='wnn.umap',group.by='cell_assignment2')
DotPlot(human_model, features=c('TFAP2C', 'NANOG', 'SOX17',
                                'PRDM1', 'NANOS3', 'TFAP2A','CRHBP',
                                'AC007326.5','LBH','GLIPR2','TMEM64','PCAT14',
                                'MPC2','CA2','PLPP1','FAM184A','NTF3','AL353807.5','PRTG','DBI'),
        group.by='cell_assignment2', scale.max = 50, col.max=1, col.min = -0.5)

genes<-c('TFAP2C', 'NANOG', 'SOX17',
         'PRDM1', 'NANOS3', 'TFAP2A','CRHBP',
         'AC007326.5','LBH','GLIPR2','TMEM64','PCAT14',
         'MPC2','CA2','PLPP1','FAM184A','NTF3','AL353807.5','PRTG','DBI')
Average_Expression<-as.data.frame(AverageExpression(human_model, assays='RNA', slot='scale.data', group.by='cell_assignment2'))
forplot5<-Average_Expression[genes,]
ComplexHeatmap::Heatmap(forplot5, cluster_rows = T)
saveRDS(human_model, 'human_model.rds')
