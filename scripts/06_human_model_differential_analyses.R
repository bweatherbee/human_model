library(Seurat)
library(Signac)
library(dplyr)
library(tidyverse)
library(Scillus)
library(RColorBrewer)

set.seed(1234)
setwd('D:/human_model/')

#combining mesoderm and amnion clusters for differential expression analysis
human_model <- readRDS("D:/human_model/human_model.rds")
Idents(human_model)
course.clusters<-c('MESO', 'L-EPI', 'AME', 'HYPO/VE', 'AME', 'MESO', 'AME', 'EXMC')
names(course.clusters)<-levels(human_model)
human_model<-RenameIdents(human_model, course.clusters)
DimPlot(human_model, reduction='wnn.umap')
human_model$course_cell_assignment<-Idents(human_model)
saveRDS(human_model, 'human_model.rds')

#ordering levels for future plotting purposes
Idents(human_model)<-human_model$cell_assignment
Idents(human_model)<-factor(x=Idents(human_model),
                                    levels=c('L-EPI', 'AM-1', 'AM-2', 'AM-3', 'MESO-1', 'MESO-2', 'HYPO/VE', 'EXMC'))
human_model$cell_assignment<-Idents(human_model)
Idents(human_model)<-human_model$course_cell_assignment
Idents(human_model)<-factor(x=Idents(human_model),
                            levels=c('L-EPI', 'AME', 'MESO', 'HYPO/VE', 'EXMC'))
human_model$course_cell_assignment<-Idents(human_model)

### RNA-level differential expressed genes ###
DefaultAssay(human_model)<-'RNA'
DEGs<-FindAllMarkers(human_model, only.pos=T, test.use='wilcox', min.pct=0.05)
L_EPI_DEGs<-subset(DEGs, subset=cluster=='L-EPI')
HYPO_VE_DEGs<-subset(DEGs, subset=cluster=='HYPO/VE')
AME_DEGs<-subset(DEGs, subset=cluster=='AME')
Meso_DEGs<-subset(DEGs, subset=cluster=='MESO')
EXMC_DEGs<-subset(DEGs, subset=cluster=='EXMC')

genes<-unique(c(head(L_EPI_DEGs$gene, 20L), head(HYPO_VE_DEGs$gene, 20L), 
                head(AME_DEGs$gene, 20L), head(Meso_DEGs$gene, 20L), head(EXMC_DEGs$gene, 20L)))

DefaultAssay(human_model)<-'RNA'
plot_heatmap(dataset=human_model, markers=genes, sort_var=c('cell_assignment'),
             anno_var=c('cell_assignment', 'sample_type'), anno_colors=list('Set1', 'Dark2'),
             hm_limit=c(-1,0,1))

write.csv(DEGs, './human_model_analysis/differential_analyses/DEGs.csv')


### Analysis on Chromatin Peaks ###

DefaultAssay(human_model)<-'peaks'
peaks_DAPs<-FindAllMarkers(human_model, only.pos = T, test.use='wilcox', min.pct=0.01)
closestfeatures<-ClosestFeature(human_model, peaks_DAPs$gene)
peaks_DAPs$closest_feature<-closestfeatures$gene_name

write.csv(peaks_DAPs, './human_model_analysis/differential_analyses/DAPs.csv')


### Analysis on Motif Activity ###

DefaultAssay(human_model)<-'chromvar'
chromVAR_DAMs<-FindAllMarkers(human_model, only.pos=T, mean.fxn=rowMeans, fc.name='avg_diff')
DefaultAssay(human_model)<-'peaks'
chromVAR_DAMs$motif.name<-ConvertMotifID(human_model, id=chromVAR_DAMs$gene)

L_EPI_DAMS<-subset(chromVAR_DAMs, cluster=='L-EPI')
HYPO_VE_DAMS<-subset(chromVAR_DAMs, cluster=='HYPO/VE')
AME_DAMS<-subset(chromVAR_DAMs, cluster=='AME')
Meso_DAMS<-subset(chromVAR_DAMs, cluster=='MESO')
EXMC_DAMS<-subset(chromVAR_DAMs, cluster=='EXMC')


key_motifs<-unique(c(head(L_EPI_DAMS$gene, 20L), head(HYPO_VE_DAMS$gene, 20L), 
                    head(AME_DAMS$gene, 20L), head(Meso_DAMS$gene, 20L), head(EXMC_DAMS$gene, 20L)))

DefaultAssay(human_model)<-'chromvar'
human_model<-ScaleData(human_model)
plot_heatmap(dataset=human_model, markers=key_motifs, sort_var=c('cell_assignment'),
             anno_var=c('cell_assignment', 'sample_type'), anno_colors=list('Set1', 'Dark2'),
             hm_limit=c(-1,0,1))


write.csv(chromVAR_DAMs, './human_model_analysis/differential_analyses/chromvarDAMs.csv')

### Additional Plots ###

#average per cluster heatmaps
Average_Accessibility<-as.data.frame(AverageExpression(human_model, assays='peaks', slot = 'scale.data'))
forplot<-Average_Accessibility[key_peaks,]
ComplexHeatmap::Heatmap(forplot)

Average_Expression<-as.data.frame(AverageExpression(human_model, assays='RNA', slot = 'scale.data'))
forplot2<-Average_Expression[genes,]
ComplexHeatmap::Heatmap(forplot2)

Average_motif_act<-as.data.frame(AverageExpression(human_model, assays='chromvar', slot = 'scale.data'))
forplot3<-Average_motif_act[key_motifs,]
ComplexHeatmap::Heatmap(forplot3)

### Generate a heatmap and dotplots of selected DEGs and motifs ###

subset_DEGs<-subset(DEGs, subset=SCT_DEGs$avg_log2FC>0.5 & SCT_DEGs$p_val_adj<0.05)
write.csv(subset_DEGs, './human_model_analysis/differential_analyses/DEGs_minlogFC05.csv')

genes<-(c('TDGF1','SOX2','UTF1','NANOG','SFRP2','FGF2','POU5F1','DNMT3B','GDF3','ZIC3',
          'GABRP','ID1','TFAP2A','VTCN1','GRHL1','ISL1','WNT6','TFAP2C','TBX3','MEIS1',
          'MIXL1','MESP1', 'CER1','WNT5A','GAL', 'GATA4', 'FOXH1', 'TBXT','SNAI1', 'EOMES',
          'POSTN','RSPO2','COL6A3','COL6A1','MMP2','HAND1','BMP4','IGF2','FLRT2','LEF1',
          'APOA1','BMP6','CDH2','FN1','FLRT3','HNF1B','HNF4A','FOXA2','SOX17','COL4A1'))
DefaultAssay(human_model)<-'RNA'
plot_heatmap(dataset=human_model, markers=unique(genes), sort_var=c('cell_assignment'),
             anno_var=c('cell_assignment','sample_type'), anno_colors=list('Set1', 'Dark2'),
             hm_limit=c(-1,0,1), row_font_size=5)
forplot4<-Average_Expression[genes,]
ComplexHeatmap::Heatmap(forplot4, cluster_rows = T)
DotPlot(human_model, assay='RNA', features=genes, cols=c('light grey', 'blue'), dot.min=0.05, col.min=0.25)


Motifs<-c('MA0142.1','MA0792.1','MA1115.1','MA1515.1','MA0039.4','MA0139.1','MA1650.1','MA0143.4','MA0162.4','MA0724.1', #L-EPI motifs: Pou5f1::Sox2,POU5F1B, POU5F1, KLF2,KLF4,CTCF,ZBTB14,SOX2,EGR1,VENTX 
          'MA0810.1','MA0815.1','MA1648.1','MA0522.3','MA1476.1','MA0878.2','MA1608.1','MA0708.1','MA0465.2','MA0508.3', # AME enriched motifs:TFAP2A(var.2),TFAP2C(var.3),TCF12(var.2), TCF3, DLX5, CDX1,Isl1,MSX2,CDX2,PRDM1
          'MA0482.2','MA0800.1','MA0513.1','MA0009.2','MA0769.2','MA0479.1','MA0521.1','MA0894.1','MA1491.1','MA0140.2', #Meso enriched motifs: GATA4, EOMES,SMAD2::SMAD3::SMAD4,TBXT,TCF7,FOXH1,Tcf12,HESX1,GLI3,GATA1::TAL1
          'MA1130.1','MA1141.1','MA1123.2','MA0092.1','MA1633.1','MA0682.2','MA0712.2','MA0591.1','MA0762.1','MA1544.1',#EXMC enriched motifs: FOSL2::JUN,FOS::JUND,TWIST1,Hand1:Tcf3,BACH1,PITX1,OTX2,Bach1::Mafk,ETV2,OVOL1
          'MA0153.2','MA0046.2','MA1494.1','MA0047.3','MA0849.1','MA1421.1','MA1104.2','MA0648.1','MA0883.1','MA0830.2') #HYPO/VE enriched motifs: HNF1B,HNF1A,HNF4A(var.2),FOXA2,FOXO6,TCF7L1,GATA6,GSC,Dmbx1,TCF4 
DefaultAssay(human_model)<-'chromvar'
plot_heatmap(dataset=human_model, markers=unique(Motifs), sort_var=c('cell_assignment'),
             anno_var=c('cell_assignment','sample_type'), anno_colors=list('Set1', 'Dark2'),
             hm_limit=c(-1,0,1), row_font_size=5)
forplot5<-Average_motif_act[Motifs,]
ComplexHeatmap::Heatmap(forplot5, cluster_rows = T)
DotPlot(human_model, assay='chromvar', features=Motifs, cols=c('light grey', 'deeppink'), dot.min=0.2, col.min=-0.5,col.max=1)