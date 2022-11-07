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

source("logisticRegression.R")

human_model <- readRDS("D:/human_model/human_model.rds")
DefaultAssay(human_model)<-'RNA'
human_model[['peaks']]<-NULL
human_model[['ATAC']]<-NULL
human_model[['GeneActivity']]<-NULL
human_model[['chromvar']]<-NULL
human_model[['SCT']]<-NULL

DEGs<-read.csv("./human_model_analysis/differential_analyses/DEGs.csv")
Genes<-unique(DEGs$gene)

##### Log Reg with Nakamura Data #####
setwd("D:/human_model")

rpm_matrix<- readRDS("D:/CynMonkey/Nakamura/RPM_matrix_hgncsymbols.rds")
metadata <- read.csv("D:/CynMonkey/Nakamura/metadata.csv")
metadata<-as.data.frame(metadata[1:390,])
rownames(metadata)<-metadata$ï..SampleID

Nakamura<-CreateSeuratObject(counts=rpm_matrix, meta.data=metadata)
Nakamura$cell_assignment<-Nakamura$cell_type
Nakamura<-NormalizeData(Nakamura)
Nakamura<-ScaleData(Nakamura)


Nakamura_int=merge(x = human_model, y = Nakamura, 
                             add.cell.ids=c("Own","Nakamura"),merge.data = TRUE)
Nakamura_int=subset(Nakamura_int, features = Genes)

dat=Nakamura_int@assays$RNA@counts

select=rownames(Nakamura_int@meta.data)[grepl("Own",rownames(Nakamura_int@meta.data))]
Nakamura_select=rownames(Nakamura_int@meta.data)[grepl("Nakamura",rownames(Nakamura_int@meta.data))]

training_dat=dat[,select]
test_dat=dat[,Nakamura_select]
classes=Nakamura_int@meta.data[select,"cell_assignment"]

fit=trainModel(refMat=training_dat,classes=classes)
saveRDS(fit,"./human_model_analysis/log_reg/Nakamura/fit.rds")

preds=predictSimilarity(fit=fit,tgtData=test_dat,classes=Nakamura_int@meta.data[Nakamura_select,"cell_assignment"],logits=F)

preds=preds[,c("L-EPI", "AM-1", "AM-2", "AM-3", "MESO-1", "MESO-2", "EXMC", "HYPO/VE")]
note_color=rep('black',9)
note_color[t(preds)>0.2]="white"
pdf("./human_model_analysis/log_reg/Nakamura/Nakamura_logreg.pdf")
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
write.csv(preds, './human_model_analysis/log_reg/Nakamura/Nakamura_preds.csv')

##### Log reg with Yang Data #####
library(readr)
counts<-read_tsv('./human_model_analysis/scmap_humanmodel/yang/CM_filtered_SoupX_corrected.tsv')
metadata<-readRDS('./human_model_analysis/scmap_humanmodel/yang/metadata_full.Rds')
metadata<-rbind(metadata[['d10']],metadata[['d12']],metadata[['d14']])
rownames(metadata)<-metadata$cell
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
colnames(counts)<-rownames(metadata)
metadata<-as.data.frame(metadata)
rownames(metadata)<-colnames(counts)

Yang<-CreateSeuratObject(counts=counts, meta.data=metadata)
Yang$cell_assignment<-Yang$lineage
Yang<-NormalizeData(Yang)
Yang<-ScaleData(Yang)


Yang_int=merge(x = human_model, y = Yang, 
                   add.cell.ids=c("Own","Yang"),merge.data = TRUE)
Yang_int=subset(Yang_int, features = Genes)

dat=Yang_int@assays$RNA@counts

select=rownames(Yang_int@meta.data)[grepl("Own",rownames(Yang_int@meta.data))]
Yang_select=rownames(Yang_int@meta.data)[grepl("Yang",rownames(Yang_int@meta.data))]

training_dat=dat[,select]
test_dat=dat[,Yang_select]
classes=Yang_int@meta.data[select,"cell_assignment"]

fit=trainModel(refMat=training_dat,classes=classes)
saveRDS(fit,"./human_model_analysis/log_reg/Yang/fit.rds")

preds=predictSimilarity(fit=fit,tgtData=test_dat,classes=Yang_int@meta.data[Yang_select,"cell_assignment"],logits=F)

preds=preds[,c("L-EPI", "AM-1", "AM-2", "AM-3", "MESO-1", "MESO-2", "EXMC", "HYPO/VE")]
note_color=rep('black',9)
note_color[t(preds)>0.2]="white"
pdf("./human_model_analysis/log_reg/Yang/Yang_logreg.pdf")
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
write.csv(preds, './human_model_analysis/log_reg/Yang/Yang_preds.csv')


##### Log reg with Ma Data #####
counts <- readRDS("D:/human_model/human_model_analysis/scmap_humanmodel/ma/counts_hgncsymbols.rds")
meta <- read_csv("human_model_analysis/scmap_humanmodel/ma/GSE130114_MF1453-meta.csv")
meta <- as.data.frame(meta)
rownames(meta)<-meta$...1
meta <- meta[,2:3]

Ma<-CreateSeuratObject(counts=counts, meta.data=meta)
Ma$cell_assignment<-Ma$annotation
Ma<-NormalizeData(Ma)
Ma<-ScaleData(Ma)


Ma_int=merge(x = huMan_model, y = Ma, 
               add.cell.ids=c("Own","Ma"),merge.data = TRUE)
Ma_int=subset(Ma_int, features = Genes)

dat=Ma_int@assays$RNA@counts

select=rownames(Ma_int@meta.data)[grepl("Own",rownames(Ma_int@meta.data))]
Ma_select=rownames(Ma_int@meta.data)[grepl("Ma",rownames(Ma_int@meta.data))]

training_dat=dat[,select]
test_dat=dat[,Ma_select]
classes=Ma_int@meta.data[select,"cell_assignment"]

fit=trainModel(refMat=training_dat,classes=classes)
saveRDS(fit,"./huMan_model_analysis/log_reg/Ma/fit.rds")

preds=predictSimilarity(fit=fit,tgtData=test_dat,classes=Ma_int@meta.data[Ma_select,"cell_assignment"],logits=F)

preds=preds[c("TE and derivatives", "postE-EPI", "postL-EPI", "E-AM", "L-AM1", "L-AM2", "E-Gast", "L-Gast1", "L-Gast2", "EXMC", "VE/YE", "E-PGC" ),c("L-EPI", "AM-1", "AM-2", "AM-3", "MESO-1", "MESO-2", "EXMC", "HYPO/VE")]
note_color=rep('black',9)
note_color[t(preds)>0.2]="white"
pdf("./huMan_model_analysis/log_reg/Ma/Ma_logreg.pdf")
gplots::heatMap.2(t(as.matrix(preds)),
                  density.info="none",  # turns off density plot inside color legend
                  trace="none",         # turns off trace lines inside the heat Map
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
write.csv(preds, './human_model_analysis/log_reg/Ma/Ma_preds.csv')


##### Log reg with Tyser Data #####

counts <- readRDS("D:/huTysern_model/huTysern_model_analysis/scTyserp_huTysernmodel/tyser/expression_values.rds")
metadata <- readRDS("D:/huTysern_model/huTysern_model_analysis/scTyserp_huTysernmodel/tyser/annot_uTyserp.rds")
rownames(metadata)<-metadata$cell_name
rownames(counts)<-metadata$cell_name

Tyser<-CreateSeuratObject(counts=t(counts), meta.data = metadata)
Tyser$cell_assignment<-Tyser$sub_cluster
Tyser<-NormalizeData(Tyser)
Tyser<-ScaleData(Tyser)


Tyser_int=merge(x = human_model, y = Tyser, 
             add.cell.ids=c("Own","Tyser"),merge.data = TRUE)
Tyser_int=subset(Tyser_int, features = Genes)

dat=Tyser_int@assays$RNA@counts

select=rownames(Tyser_int@meta.data)[grepl("Own",rownames(Tyser_int@meta.data))]
Tyser_select=rownames(Tyser_int@meta.data)[grepl("Tyser",rownames(Tyser_int@meta.data))]

training_dat=dat[,select]
test_dat=dat[,Tyser_select]
classes=Tyser_int@meta.data[select,"cell_assignment"]

fit=trainModel(refMat=training_dat,classes=classes)
saveRDS(fit,"./human_model_analysis/log_reg/Tyser/fit.rds")

preds=predictSimilarity(fit=fit,tgtData=test_dat,classes=Tyser_int@meta.data[Tyser_select,"cell_assignment"],logits=F)

preds=preds[,c("L-EPI", "AM-1", "AM-2", "AM-3", "MESO-1", "MESO-2", "EXMC", "HYPO/VE")]
note_color=rep('black',9)
note_color[t(preds)>0.2]="white"
pdf("./human_model_analysis/log_reg/Tyser/Tyser_logreg.pdf")
gplots::heatmap.2(t(as.matrix(preds)),
                  density.info="none",  # turns off density plot inside color legend
                  trace="none",         # turns off trace lines inside the heat Tyserp
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
                  notecex=0.5,
                  mar=c(20,10))
dev.off()
write.csv(preds, './human_model_analysis/log_reg/Tyser/Tyser_preds.csv')

##### Log reg with Mole Data #####

MZG <- readRDS("D:/realigning_for_int/MZG/embryos_integrated.RDS")
MZG[['integrated']]<-NULL
MZG$cell_assignment<-MZG$cell_type

Mole_int=merge(x = human_model, y = MZG, 
               add.cell.ids=c("Own","Mole"),merge.data = TRUE)
Mole_int=subset(Mole_int, features = Genes)

dat=Mole_int@assays$RNA@counts

select=rownames(Mole_int@meta.data)[grepl("Own",rownames(Mole_int@meta.data))]
Mole_select=rownames(Mole_int@meta.data)[grepl("Mole",rownames(Mole_int@meta.data))]

training_dat=dat[,select]
test_dat=dat[,Mole_select]
classes=Mole_int@meta.data[select,"cell_assignment"]

fit=trainModel(refMat=training_dat,classes=classes)
saveRDS(fit,"./human_model_analysis/log_reg/Mole/fit.rds")

preds=predictSimilarity(fit=fit,tgtData=test_dat,classes=Mole_int@meta.data[Mole_select,"cell_assignment"],logits=F)

preds=preds[,c("L-EPI", "AM-1", "AM-2", "AM-3", "MESO-1", "MESO-2", "EXMC", "HYPO/VE")]
note_color=rep('black',9)
note_color[t(preds)>0.2]="white"
pdf("./human_model_analysis/log_reg/Mole/Mole_logreg.pdf")
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
write.csv(preds, './human_model_analysis/log_reg/Mole/Mole_preds.csv')

