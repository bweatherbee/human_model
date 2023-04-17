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
library(tidyverse)

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

#Train on Embryo#
training_dat=dat[,Nakamura_select]
test_dat=dat[,select]
classes=Nakamura_int@meta.data[Nakamura_select,"cell_assignment"]

fit=trainModel(refMat=training_dat,classes=classes)
saveRDS(fit,"./human_model_analysis/log_reg/Nakamura/embryo_train_fit.rds")

preds=predictSimilarity(fit=fit,tgtData=test_dat,classes=Nakamura_int@meta.data[select,"cell_assignment"],logits=F)

#Can't do similarity to VE/YS because less than 8 cells so excluding
preds=preds[c("L-EPI", "AM-1", "AM-2", "AM-3", "MESO-1", "MESO-2", "EXMC", "HYPO/VE"),c('ICM', 'Pre-EPI', 'PostE-EPI', 'PostL-EPI', 'Gast1', 'Gast2a', 'Gast2b', 'PreE-TE', 'PreL-TE', 'Post-paTE', 'Hypoblast', 'EXMC')]
note_color=rep('black',9)
note_color[t(preds)>0.2]="white"
  pdf("./human_model_analysis/log_reg/Nakamura/Nakamura_logreg_embryotrained.pdf")
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
  write.csv(preds, './human_model_analysis/log_reg/Nakamura/Nakamura_embryotrain_preds.csv')

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

#embryo trained#
training_dat=dat[,Yang_select]
test_dat=dat[,select]
classes=Yang_int@meta.data[Yang_select,"cell_assignment"]

fit=trainModel(refMat=training_dat,classes=classes)
saveRDS(fit,"./human_model_analysis/log_reg/Yang/embryotrained_fit.rds")

preds=predictSimilarity(fit=fit,tgtData=test_dat,classes=Yang_int@meta.data[select,"cell_assignment"],logits=F)

#everything has >0.90 with extraembryonic mesoderm so excluding from plot as potentially a strange population for this geneset
preds=preds[c("L-EPI", "AM-1", "AM-2", "AM-3", "MESO-1", "MESO-2", "EXMC", "HYPO/VE"),c('epiblast', 'transition', 'amnion', 'amnion 1', 'amnion 2', 'mesoderm 1', 'mesoderm 2', 'extraembryonic mesenchyme', 'endoderm', 'trophoblast')]
note_color=rep('black',9)
note_color[t(preds)>0.2]="white"
  pdf("./human_model_analysis/log_reg/Yang/Yang_embryotrained_logreg.pdf")
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
                    mar=c(10,9))
  dev.off()
  write.csv(preds, './human_model_analysis/log_reg/Yang/Yang_embryotrained_preds.csv')


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

#embryo trained#
training_dat=dat[,Ma_select]
test_dat=dat[,select]
classes=Ma_int@meta.data[Ma_select,"cell_assignment"]

fit=trainModel(refMat=training_dat,classes=classes)
saveRDS(fit,"./huMan_model_analysis/log_reg/Ma/embryotrained_fit.rds")

preds=predictSimilarity(fit=fit,tgtData=test_dat,classes=Ma_int@meta.data[select,"cell_assignment"],logits=F)

preds=preds[c("L-EPI", "AM-1", "AM-2", "AM-3", "MESO-1", "MESO-2", "EXMC", "HYPO/VE"),c("TE and derivatives", "postE-EPI", "postL-EPI", "E-AM", "L-AM1", "L-AM2", "E-Gast", "L-Gast1", "L-Gast2", "EXMC", "VE/YE", "E-PGC" )]
note_color=rep('black',9)
note_color[t(preds)>0.2]="white"
  pdf("./huMan_model_analysis/log_reg/Ma/Ma_logreg_embryotrained.pdf")
  gplots::heatmap.2(t(as.matrix(preds)),
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
  write.csv(preds, './human_model_analysis/log_reg/Ma/Ma_embryotrained_preds.csv')


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

#Trained on embryo
training_dat=dat[,Tyser_select]
test_dat=dat[,select]
classes=Tyser_int@meta.data[Tyser_select,"cell_assignment"]

fit=trainModel(refMat=training_dat,classes=classes)
saveRDS(fit,"./human_model_analysis/log_reg/Tyser/embryotrained_fit.rds")

preds=predictSimilarity(fit=fit,tgtData=test_dat,classes=Tyser_int@meta.data[select,"cell_assignment"],logits=F)

preds=preds[c("L-EPI", "AM-1", "AM-2", "AM-3", "MESO-1", "MESO-2", "EXMC", "HYPO/VE"),c('Epiblast', 'Amnion', 'Primitive Streak', 
                                                                                        'Nascent Mesoderm', 'Emergent Mesoderm', 'Advanced Mesoderm', 
                                                                                        'Axial Mesoderm', 'DE(P)', 'DE(NP)', 'NNE', 'Hypoblast', 
                                                                                        'YS Endoderm', 'YS Mesoderm', 'PGC', 'Blood Progenitors', 
                                                                                        'Hemogenic Endothelium','Erythro-Myeloid Progenitors', 
                                                                                        'Erythroblasts', 'Myeloid Progenitors')]
note_color=rep('black',9)
note_color[t(preds)>0.2]="white"
  pdf("./human_model_analysis/log_reg/Tyser/Tyser_logreg_embryotrained.pdf")
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
                    mar=c(10,20))
  dev.off()
  write.csv(preds, './human_model_analysis/log_reg/Tyser/Tyser_preds_embryotrained.csv')

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

#Trained on embryo
training_dat=dat[,Mole_select]
test_dat=dat[,select]
classes=Mole_int@meta.data[Mole_select,"cell_assignment"]

fit=trainModel(refMat=training_dat,classes=classes)
saveRDS(fit,"./human_model_analysis/log_reg/Mole/fit_embryotrained.rds")

preds=predictSimilarity(fit=fit,tgtData=test_dat,classes=Mole_int@meta.data[select,"cell_assignment"],logits=F)

preds=preds[c("L-EPI", "AM-1", "AM-2", "AM-3", "MESO-1", "MESO-2", "EXMC", "HYPO/VE"),]
note_color=rep('black',9)
note_color[t(preds)>0.2]="white"
  pdf("./human_model_analysis/log_reg/Mole/Mole_logreg_embryotrained.pdf")
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
  write.csv(preds, './human_model_analysis/log_reg/Mole/Mole_preds_embryotrained.csv')
  
  ##### Log reg with Kagawa Blastoid Data #####
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
  blastoid$cell_assignment<-blastoid$blastoid.Fig2b.lowres
  blastoid<-NormalizeData(blastoid)
  blastoid<-ScaleData(blastoid)
  
  
  blastoid_int=merge(x = human_model, y = blastoid, 
                     add.cell.ids=c("Own","Kagawa"),merge.data = TRUE)
  blastoid_int=subset(blastoid_int, features = Genes)
  
  dat=blastoid_int@assays$RNA@counts
  
  select=rownames(blastoid_int@meta.data)[grepl("Own",rownames(blastoid_int@meta.data))]
  blastoid_select=rownames(blastoid_int@meta.data)[grepl("Kagawa",rownames(blastoid_int@meta.data))]

  #Trained on blastoids
  training_dat=dat[,blastoid_select]
  test_dat=dat[,select]
  classes=blastoid_int@meta.data[blastoid_select,"cell_assignment"]
  
  fit=trainModel(refMat=training_dat,classes=classes)
  saveRDS(fit,"./human_model_analysis/log_reg/Kagawa/blastoid_trained_fit.rds")
  
  preds=predictSimilarity(fit=fit,tgtData=test_dat,classes=blastoid_int@meta.data[select,"cell_assignment"],logits=F)
  
  preds=preds[c("L-EPI", "AM-1", "AM-2", "AM-3", "MESO-1", "MESO-2", "EXMC", "HYPO/VE"), c('0', '1', '3')]
  note_color=rep('black',9)
  note_color[t(preds)>0.2]="white"
    pdf("./human_model_analysis/log_reg/Kagawa/Kagawa_logreg_blastoid_trained.pdf")
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
    write.csv(preds, './human_model_analysis/log_reg/Kagawa/Kagawa_preds_blastoid_trained.csv')
    
##### PHAM et al #####
  count_matrix <- read.csv("human_model_analysis/log_reg/Pham/count_matrix.csv", row.names = 1)
  annotation <- read.csv("human_model_analysis/log_reg/Pham/annotation.csv", sep=';', row.names = 1)
    
  pham<-CreateSeuratObject(counts=count_matrix, meta.data=annotation)
    
  pham$cell_assignment<-pham$Cell.type
  pham<-NormalizeData(pham)
  pham<-ScaleData(pham)
    
    
  pham_int=merge(x = human_model, y = pham, 
                   add.cell.ids=c("Own","pham"),merge.data = TRUE)
  pham_int=subset(pham_int, features = Genes)
    
  dat=pham_int@assays$RNA@counts
    
  select=rownames(pham_int@meta.data)[grepl("Own",rownames(pham_int@meta.data))]
  pham_select=rownames(pham_int@meta.data)[grepl("pham",rownames(pham_int@meta.data))]

#Trained on Pham data
  training_dat=dat[,pham_select]
  test_dat=dat[,select]
  classes=pham_int@meta.data[pham_select,"cell_assignment"]
  
  fit=trainModel(refMat=training_dat,classes=classes)
  saveRDS(fit,"./human_model_analysis/log_reg/pham/fit_phamtrained.rds")
  
  preds=predictSimilarity(fit=fit,tgtData=test_dat,classes=pham_int@meta.data[select,"cell_assignment"],logits=F)
  
  preds=preds[c("L-EPI", "AM-1", "AM-2", "AM-3", "MESO-1", "MESO-2", "EXMC", "HYPO/VE"), c('naive', 'primed', 'TSC', 'EXMC')]
  note_color=rep('black',9)
  note_color[t(preds)>0.2]="white"
    pdf("./human_model_analysis/log_reg/pham/pham_logreg_phamtrained.pdf")
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
    write.csv(preds, './human_model_analysis/log_reg/pham/pham_preds_phamtrained.csv')
    
##### Zheng et al #####
  PASE <- readRDS("D:/human_model/human_model_analysis/log_reg/Zheng/GSE134571_Posterior48h_H9_Amnion_Merged.rds")
    
  new.cluster.ids<-c("hESCs", "Transwell-AMLC", "AMLC", "hPGCLC", "MeLC1", "MeLC2")
  names(new.cluster.ids)<-levels(PASE)
  PASE<-RenameIdents(PASE, new.cluster.ids)
    
  PASE$cell_assignment<-Idents(PASE)
  PASE<-NormalizeData(PASE)
  PASE<-ScaleData(PASE)
    
    
  PASE_int=merge(x = human_model, y = PASE, 
                   add.cell.ids=c("Own","PASE"),merge.data = TRUE)
  PASE_int=subset(PASE_int, features = Genes)
    
  dat=PASE_int@assays$RNA@counts
    
  select=rownames(PASE_int@meta.data)[grepl("Own",rownames(PASE_int@meta.data))]
  PASE_select=rownames(PASE_int@meta.data)[grepl("PASE",rownames(PASE_int@meta.data))]
  
#Trained on PASE
  training_dat=dat[,PASE_select]
  test_dat=dat[,select]
  classes=PASE_int@meta.data[PASE_select,"cell_assignment"]
  
  fit=trainModel(refMat=training_dat,classes=classes)
  saveRDS(fit,"./human_model_analysis/log_reg/Zheng/fit_pasetrained.rds")
  
  preds=predictSimilarity(fit=fit,tgtData=test_dat,classes=PASE_int@meta.data[select,"cell_assignment"],logits=F)
  
#Excluding transwell-AMLCs since they are not in the 3D PASE Model  
  preds=preds[c("L-EPI", "AM-1", "AM-2", "AM-3", "MESO-1", "MESO-2", "EXMC", "HYPO/VE"),c('hESCs', 'AMLC', 'MeLC1', 'MeLC2', 'hPGCLC')]
  note_color=rep('black',9)
  note_color[t(preds)>0.2]="white"
    pdf("./human_model_analysis/log_reg/Zheng/Zheng_logreg_pasetrained.pdf")
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
    write.csv(preds, './human_model_analysis/log_reg/Zheng/Zheng_preds_pasetrained.csv')
