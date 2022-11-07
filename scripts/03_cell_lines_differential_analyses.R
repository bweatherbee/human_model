library(Seurat)
library(Signac)
library(dplyr)
library(tidyverse)
library(Scillus)
library(RColorBrewer)

set.seed(1234)
setwd('D:/human_model/')

cell_lines <- readRDS("D:/human_model/cell_lines_complete.rds")

### RNA-level differential expressed genes ###
DefaultAssay(cell_lines)<-'RNA'
SCT_DEGs<-FindAllMarkers(cell_lines, only.pos=T, test.use='wilcox', min.pct=0.05)
WT_DEGs<-subset(SCT_DEGs, subset=cluster=='RSeT_wt')
G6_S17_DEGs<-subset(SCT_DEGs, subset=cluster=='g6_s17_RSeT_induced')
G3_AP2y_DEGs<-subset(SCT_DEGs, subset=cluster=='g3_ap2y_RSeT_induced')

genes<-unique(c(head(WT_DEGs$gene, 20L), head(G6_S17_DEGs$gene, 20L), head(G3_AP2y_DEGs$gene, 20L)))

DefaultAssay(cell_lines)<-'RNA'
plot_heatmap(dataset=cell_lines, markers=genes, sort_var=c('sample_type'),
             anno_var=c('sample_type'), anno_colors=list('Set1', 'Dark2'),
             hm_limit=c(-1,0,1))

write.csv(SCT_DEGs, './cell_lines_analysis/differential_analyses/DEGs.csv')


### Analysis on Chromatin Peaks ###

DefaultAssay(cell_lines)<-'peaks'
peaks_DAPs<-FindAllMarkers(cell_lines, only.pos = T, test.use='wilcox', min.pct=0.01)
closestfeatures<-ClosestFeature(cell_lines, peaks_DAPs$gene)
peaks_DAPs$closest_feature<-closestfeatures$gene_name

WT_peaks<-subset(peaks_DAPs, cluster=='RSeT_wt')
WT_enrich.motifs<-FindMotifs(cell_lines, features=WT_peaks$gene)

G6_S17_peaks<-subset(peaks_DAPs, cluster=='g6_s17_RSeT_induced')
G6_S17_enrich.motifs<-FindMotifs(cell_lines, features=G6_S17_peaks$gene)

G3_AP2y_peaks<-subset(peaks_DAPs, cluster=='g3_ap2y_RSeT_induced')
G3_AP2y_enrich.motifs<-FindMotifs(cell_lines, features=G3_AP2y_peaks$gene)

key_peaks<-unique(c(head(WT_peaks$gene, 20L), head(G6_S17_peaks$gene, 20L),head(G3_AP2y_peaks$gene, 20L)))

plot_heatmap(dataset=cell_lines, markers=key_peaks, sort_var=c('sample_type'),
             anno_var=c('sample_type'), anno_colors=list('Set1', 'Dark2'),
             hm_limit=c(-1,0,1))

write.csv(peaks_DAPs, './cell_lines_analysis/differential_analyses/DAPs.csv')
write.csv(WT_enrich.motifs, './cell_lines_analysis/differential_analyses/WT_peaks_enrichedmotifs.csv')
write.csv(G6_S17_enrich.motifs, './cell_lines_analysis/differential_analyses/G6_S17_peaks_enrichedmotifs.csv')
write.csv(G3_AP2y_enrich.motifs, './cell_lines_analysis/differential_analyses/G3_AP2y_peaks_enrichedmotifs.csv')


### Analysis on Motif Activity ###

DefaultAssay(cell_lines)<-'chromvar'
chromVAR_DAMs<-FindAllMarkers(cell_lines, only.pos=T, mean.fxn=rowMeans, fc.name='avg_diff')
DefaultAssay(cell_lines)<-'peaks'
chromVAR_DAMs$motif.name<-ConvertMotifID(cell_lines, id=chromVAR_DAMs$gene)

WT_DAMS<-subset(chromVAR_DAMs, cluster=='RSeT_wt')
G6_S17_DAMS<-subset(chromVAR_DAMs, cluster=='g6_s17_RSeT_induced')
G3_AP2y_DAMS<-subset(chromVAR_DAMs, cluster=='g3_ap2y_RSeT_induced')


key_motifs<-unique(c(head(WT_DAMS$gene, 20L), head(G6_S17_DAMS$gene, 20L),head(G3_AP2y_DAMS$gene, 20L)))
DefaultAssay(cell_lines)<-'chromvar'
cell_lines<-ScaleData(cell_lines)
plot_heatmap(dataset=cell_lines, markers=key_motifs, sort_var=c('sample_type'),
             anno_var=c('sample_type'), anno_colors=list('Set1', 'Dark2'),
             hm_limit=c(-1,0,1))


write.csv(chromVAR_DAMs, './cell_lines_analysis/differential_analyses/chromvarDAMs.csv')


### Average per cell line heatmaps ###
Average_Accessibility<-as.data.frame(AverageExpression(cell_lines, assays='peaks', slot = 'scale.data'))
forplot<-Average_Accessibility[key_peaks,]
ComplexHeatmap::Heatmap(forplot)

Average_Expression<-as.data.frame(AverageExpression(cell_lines, assays='RNA', slot = 'scale.data'))
forplot2<-Average_Expression[genes,]
ComplexHeatmap::Heatmap(forplot2)

Average_motif_act<-as.data.frame(AverageExpression(cell_lines, assays='chromvar', slot = 'scale.data'))
forplot3<-Average_motif_act[key_motifs,]
ComplexHeatmap::Heatmap(forplot3)

### Generate a heatmap of selected DEGs and differentially accessible motifs with log FC >0.5 and adjusted pval <0.05 ###

subset_DEGs<-subset(SCT_DEGs, subset=SCT_DEGs$avg_log2FC>0.5 & SCT_DEGs$p_val_adj<0.05)
write.csv(subset_DEGs, './cell_lines_analysis/differential_analyses/DEGs_minlogFC05.csv')

WT_DEGs<-subset(subset_DEGs, subset=cluster=='RSeT_wt')
G6_S17_DEGs<-subset(subset_DEGs, subset=cluster=='g6_s17_RSeT_induced')
G3_AP2y_DEGs<-subset(subset_DEGs, subset=cluster=='g3_ap2y_RSeT_induced')

genes<-(c('SOX2', 'POU5F1', 'FGF2', 'TDGF1', 'DNMT3B', 'FABP5', 'LIN28A', 'HSPE1', 'CDH1', 'UBP1', 'UTF1', 'PRDM14', 'TERF1', 'BNC2', 'PTMA', 'GRID2', 'ESRG', 'PODXL', 'AASS', 'APOE',
                'GATA6', 'SOX17', 'SOX21', 'COL15A1', 'CDH12', 'EGF', 'RSPO3', 'NEDD4', 'IGF1', 'TET2', 'GPX2', 'COL18A1', 'FZD6', 'NID2', 'TWIST1', 'ETV6', 'CDH8', 'LAMA4', 'PRDM1', 'ADAM19',
                'GATA3', 'TFAP2C', 'KRT19', 'KRT18', 'KRT8', 'ID2', 'BAMBI', 'HES1', 'JUNB', 'MEIS3', 'SMAD6', 'ELF3', 'TGIF1', 'KLF3', 'BMP4', 'SOX11', 'HAND1', 'CDH11', 'FOXP2', 'GFP'))
DefaultAssay(cell_lines)<-'RNA'
plot_heatmap(dataset=cell_lines, markers=unique(genes), sort_var=c('sample_type'),
             anno_var=c('sample_type'), anno_colors=list('Set1', 'Dark2'),
             hm_limit=c(-1,0,1), row_font_size=5)
forplot4<-Average_Expression[genes,]
ComplexHeatmap::Heatmap(forplot4, cluster_rows = T)
DotPlot(cell_lines, assay='RNA', features=genes, cols=c('light grey', 'blue'), dot.min=0.05, col.min=0.25)


Motifs<-c('MA0142.1','MA0792.1','MA0697.1','MA0894.1','MA0712.2','MA1566.1','MA0724.1','MA1115.1','MA0039.4','MA0508.3', #RSeT motifs: POU5F1::SOX2, POU5F1B, ZIC3, HESX1, OTX2, TBX3, VENTX,POU5F1,KLF4, PRDM1
          'MA0078.1','MA1104.2','MA0482.2','MA1152.1','MA1153.1','MA0479.1', 'MA1129.1','MA0465.2','MA1545.1', 'MA0527.1', #G6_S17 induced motifs: Sox17, GATA6, GATA4, SOX15, Smad4, FOXH1, FOSL1::JUN(var.2), CDX2, OVOL2, TCF4
          'MA1121.1','MA0809.2', 'MA0769.2', 'MA0003.4', 'MA1622.1','MA0814.2','MA1111.1','MA1636.1','MA0496.3','MA0037.3') #G3_AP2y induced motifs: TEAD2, TEAD4, TCF7, TFAP2A, Smad2::Smad3, TFAP2C(var.2), NR2F2, CEBPG(var.2), MAFK, GATA3 (note not diff in this pop, likely gata factor redun.)
DefaultAssay(cell_lines)<-'chromvar'
plot_heatmap(dataset=cell_lines, markers=unique(Motifs), sort_var=c('sample_type'),
             anno_var=c('sample_type'), anno_colors=list('Set1', 'Dark2'),
             hm_limit=c(-1,0,1), row_font_size=5)
forplot5<-Average_motif_act[Motifs,]
ComplexHeatmap::Heatmap(forplot5, cluster_rows = T)
DotPlot(cell_lines, assay='chromvar', features=Motifs, cols=c('light grey', 'deeppink'), dot.min=0.1)
