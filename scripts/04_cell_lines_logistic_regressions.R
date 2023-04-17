options(stringsAsFactors = F)
library(Seurat)
library(dplyr)
library(magrittr)
library(ggplot2)
library(cowplot)
library(patchwork)
library(viridis)
library(data.table)
library(RColorBrewer)

setwd("D:/human_model")
source("./scripts/logisticRegression.R")

MZG <- readRDS("D:/realigning_for_int/MZG/MZG.RDS")
DefaultAssay(MZG)<-'RNA'
MZG[['dorothea']]<-NULL
MZG[['SCENIC']]<-NULL
MZG[['integrated']]<-NULL

cell_lines <- readRDS("D:/human_model/cell_lines_complete.rds")
DefaultAssay(cell_lines)<-'RNA'
cell_lines[['SCT']]<-NULL
cell_lines[['peaks']]<-NULL
cell_lines[['ATAC']]<-NULL
cell_lines[['chromvar']]<-NULL
cell_lines[['GeneActivity']]<-NULL
gc()
cell_lines$cell_type<-cell_lines$sample_type

combo=merge(x = MZG, y = cell_lines, 
            add.cell.ids=c("Embryo","Cell_lines"),merge.data = TRUE)

dat=combo@assays$RNA@counts
Idents(combo)<-combo$cell_type
combo$cell_type<-Idents(combo)

select=rownames(combo@meta.data)[grepl("Embryo",rownames(combo@meta.data))]
cells_select=rownames(combo@meta.data)[grepl("Cell_lines",rownames(combo@meta.data))]

training_dat=dat[,select]
classes=combo@meta.data[select,"cell_type"]
test_dat=dat[,cells_select]

fit=trainModel(refMat=training_dat,classes=classes)
saveRDS(fit,"./cell_lines_analysis/log_reg/logistic_regression_MZGembryos_data_fit_allgenes.rds")

preds=predictSimilarity(fit=fit,tgtData=test_dat,classes=combo@meta.data[cells_select,"cell_type"],logits=F)

preds=preds[c('RSeT_wt', 'g6_s17_RSeT_induced', 'g3_ap2y_RSeT_induced'),c("Epiblast","Hypoblast","Cytotrophoblasts")]
note_color=rep('black',9)
note_color[t(preds)>0.5]="white"
pdf("./cell_lines_analysis/log_reg/cells_MZG_embryo_logistic_regression_allgenes.pdf")
gplots::heatmap.2(t(as.matrix(preds)),
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          dendrogram="none",
          notecol = "black",
          col=colorRampPalette(brewer.pal(9,"Blues")[2:9])(100),
          Rowv = F,
          Colv = F,
          scale='none',
          key=F,
          cellnote=round(t(as.matrix(preds)),digits=2),
          cexRow=1.5,
          cexCol=1.5,
          notecex=1.0,
          mar=c(10,6.5))
dev.off()
write.csv(preds, './cell_lines_analysis/log_reg/preds_all_genes.csv')

###On DEGs###
library(readr)
DEGs <- read_csv("cell_lines_analysis/differential_analyses/DEGs_minlogFC05.csv")

dat=combo@assays$RNA@counts[DEGs$gene,]

training_dat=dat[,select]
classes=combo@meta.data[select,"cell_type"]
test_dat=dat[,cells_select]

fit=trainModel(refMat=training_dat,classes=classes)
saveRDS(fit,"./cell_lines_analysis/log_reg/logistic_regression_MZGembryos_data_fit_cell_lineDEGs.rds")

preds=predictSimilarity(fit=fit,tgtData=test_dat,classes=combo@meta.data[cells_select,"cell_type"],logits=F)

preds=preds[c('RSeT_wt', 'g6_s17_RSeT_induced', 'g3_ap2y_RSeT_induced'),c("Epiblast","Hypoblast","Cytotrophoblasts")]
note_color=rep('black',9)
note_color[t(preds)>0.5]="white"
pdf("./cell_lines_analysis/log_reg/cells_MZG_embryo_logistic_regression_cell_line_DEGs.pdf")
gplots::heatmap.2(t(as.matrix(preds)),
                  density.info="none",  # turns off density plot inside color legend
                  trace="none",         # turns off trace lines inside the heat map
                  dendrogram="none",
                  notecol = "black",
                  col=colorRampPalette(brewer.pal(9,"Blues")[2:9])(100),
                  Rowv = F,
                  Colv = F,
                  scale='none',
                  key=F,
                  cellnote=round(t(as.matrix(preds)),digits=2),
                  cexRow=1.5,
                  cexCol=1.5,
                  notecex=1.0,
                  mar=c(10,6.5))
dev.off()
write.csv(preds, './cell_lines_analysis/log_reg/preds_DEGs.csv')


##### Comparisons with Kagawa et al ####
DEGs <- read_csv("cell_lines_analysis/differential_analyses/DEGs_minlogFC05.csv")

counts <- readRDS("D:/human_model/human_model_analysis/log_reg/Kagawa/counts.RDS")
rownames(counts)<-make.unique(gsub(".*-","", counts[,1]))
counts<-counts[,2:2716]
metadata<-as.data.frame(read_tsv('./human_model_analysis/log_reg/Kagawa/metadata.tsv'))
rownames(metadata)<-metadata$sampleid
metadata_0<-subset(metadata, metadata$blastoid.Fig2b.lowres=='0')
metadata_0<-subset(metadata_0, metadata_0$treatment!='PXGL')
metadata_1<-subset(metadata, metadata$blastoid.Fig2b.lowres=='1')
metadata_3<-subset(metadata, metadata$blastoid.Fig2b.lowres=='3')

blastoid_metadata<-rbind(metadata_0, metadata_1, metadata_3)
blastoid_counts<-counts[,rownames(blastoid_metadata)]
blastoid<-CreateSeuratObject(counts=blastoid_counts, meta.data=blastoid_metadata)
blastoid$cell_type<-blastoid$blastoid.Fig2b.lowres
blastoid<-NormalizeData(blastoid)
blastoid<-ScaleData(blastoid)

combo=merge(x = blastoid, y = cell_lines, 
            add.cell.ids=c("Blastoid","Cell_lines"),merge.data = TRUE)

dat=combo@assays$RNA@counts[DEGs$gene,]
Idents(combo)<-combo$cell_type
combo$cell_type<-Idents(combo)

select=rownames(combo@meta.data)[grepl("Blastoid",rownames(combo@meta.data))]
cells_select=rownames(combo@meta.data)[grepl("Cell_lines",rownames(combo@meta.data))]

training_dat=dat[,select]
classes=combo@meta.data[select,"cell_type"]
test_dat=dat[,cells_select]

fit=trainModel(refMat=training_dat,classes=classes)
saveRDS(fit,"./cell_lines_analysis/log_reg/Kagawa/logistic_regression_Kagawa_data_fit_cell_lineDEGs.rds")

preds=predictSimilarity(fit=fit,tgtData=test_dat,classes=combo@meta.data[cells_select,"cell_type"],logits=F)

preds=preds[c('RSeT_wt', 'g6_s17_RSeT_induced', 'g3_ap2y_RSeT_induced'),c("0", "1", "3")]
note_color=rep('black',9)
note_color[t(preds)>0.5]="white"
  pdf("./cell_lines_analysis/log_reg/Kagawa/cells_Kagawa_embryo_logistic_regression_cell_line_DEGs.pdf")
  gplots::heatmap.2(t(as.matrix(preds)),
                    density.info="none",  # turns off density plot inside color legend
                    trace="none",         # turns off trace lines inside the heat map
                    dendrogram="none",
                    notecol = "black",
                    col=colorRampPalette(brewer.pal(9,"Blues")[2:9])(100),
                    Rowv = F,
                    Colv = F,
                    scale='none',
                    key=F,
                    cellnote=round(t(as.matrix(preds)),digits=2),
                    cexRow=1.5,
                    cexCol=1.5,
                    notecex=1.0,
                    mar=c(10,6.5))
  dev.off()
  write.csv(preds, './cell_lines_analysis/log_reg/Kagawa/preds_DEGs.csv')
  
  #### Pham et al ####
  
  count_matrix <- read.csv("human_model_analysis/log_reg/Pham/count_matrix.csv", row.names = 1)
  annotation <- read.csv("human_model_analysis/log_reg/Pham/annotation.csv", sep=';', row.names = 1)
  
  pham<-CreateSeuratObject(counts=count_matrix, meta.data=annotation)
  
  pham$cell_type<-pham$Cell.type
  pham<-NormalizeData(pham)
  pham<-ScaleData(pham)
  
  combo=merge(x = pham, y = cell_lines, 
              add.cell.ids=c("Pham","Cell_lines"),merge.data = TRUE)
  
  dat=combo@assays$RNA@counts[DEGs$gene,]
  Idents(combo)<-combo$cell_type
  combo$cell_type<-Idents(combo)
  
  select=rownames(combo@meta.data)[grepl("Pham",rownames(combo@meta.data))]
  cells_select=rownames(combo@meta.data)[grepl("Cell_lines",rownames(combo@meta.data))]
  
  training_dat=dat[,select]
  classes=combo@meta.data[select,"cell_type"]
  test_dat=dat[,cells_select]
  
  fit=trainModel(refMat=training_dat,classes=classes)
  saveRDS(fit,"./cell_lines_analysis/log_reg/Pham/logistic_regression_Pham_data_fit_cell_lineDEGs.rds")
  
  preds=predictSimilarity(fit=fit,tgtData=test_dat,classes=combo@meta.data[cells_select,"cell_type"],logits=F)
  
  preds=preds[c('RSeT_wt', 'g6_s17_RSeT_induced', 'g3_ap2y_RSeT_induced'),c("naive", "primed", "EXMC", "TSC")]
  note_color=rep('black',9)
  note_color[t(preds)>0.5]="white"
    pdf("./cell_lines_analysis/log_reg/Pham/cells_Pham_embryo_logistic_regression_cell_line_DEGs.pdf")
    gplots::heatmap.2(t(as.matrix(preds)),
                      density.info="none",  # turns off density plot inside color legend
                      trace="none",         # turns off trace lines inside the heat map
                      dendrogram="none",
                      notecol = "black",
                      col=colorRampPalette(brewer.pal(9,"Blues")[2:9])(100),
                      Rowv = F,
                      Colv = F,
                      scale='none',
                      key=F,
                      cellnote=round(t(as.matrix(preds)),digits=2),
                      cexRow=1.5,
                      cexCol=1.5,
                      notecex=1.0,
                      mar=c(10,6.5))
    dev.off()
    write.csv(preds, './cell_lines_analysis/log_reg/Pham/preds_DEGs.csv')
