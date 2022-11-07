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
gc()
cell_lines$cell_type<-cell_lines$sample_type

combo=merge(x = MZG, y = cell_lines, 
            add.cell.ids=c("Embryo","Cell_lines"),merge.data = TRUE)

dat=combo@assays$RNA@counts
Idents(combo)<-combo$cell_type
combo$cell_type<-Idents(combo)

select=rownames(combo@meta.data)[grepl("Embryo",rownames(combo@meta.data))]
cells_select=rownames(combo@meta.data)[grepl("Cell_lines",rownames(combo@meta.data))]

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