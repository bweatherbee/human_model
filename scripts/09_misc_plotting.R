library(Seurat)
library(Signac)
library(RColorBrewer)
library(ggplot2)
library(SCpubr)

human_model <- readRDS("D:/human_model/human_model.rds")
cell_lines <- readRDS("D:/human_model/cell_lines_complete.rds")
DefaultAssay(human_model)<-'RNA'

#### QC Plots ####
Idents(human_model)<-human_model$sample_type
VlnPlot(human_model, features=c('nCount_RNA', 'percent.mt', 'nCount_peaks'),
        pt.size=0)
DefaultAssay(human_model)<-'peaks'
FragmentHistogram(human_model)
DefaultAssay(human_model)<-'ATAC'
TSSPlot(human_model)

Idents(cell_lines)<-cell_lines$sample_type
VlnPlot(cell_lines, features=c('nCount_RNA', 'percent.mt', 'nCount_peaks'),
        pt.size=0, cols=c('chartreuse2', 'darkorange', '#9d388b'))
DefaultAssay(cell_lines)<-'peaks'
FragmentHistogram(cell_lines)
DefaultAssay(cell_lines)<-'ATAC'
TSSPlot(cell_lines)


##### Gene Expression Plots ####
#colored dimplot for cell lines Figure 1
DimPlot(cell_lines, reduction='umap.rna', cols=c('chartreuse2', 'darkorange', '#9d388b'))

#Featureplots for cell_lines
DefaultAssay(cell_lines)<-'RNA'
FeaturePlot(cell_lines, features=c('POU5F1', 'NANOG', 'SOX2', 'DPPA5',
                                   'CDX2', 'ISL1', 'TFAP2A', 'VTCN1', 
                                   'GATA6', 'SOX17', 'GATA4', 'CDH2',
                                   'GATA3', 'TFAP2C', 'GATA2', 'GFP',
                                   'TBXT', 'EOMES', 'SNAI1', 'CDH1'), order=T,
            min.cutoff='q20', max.cutoff = 'q98')


#cardinal marker expression for Figure 3 Supplement
FeaturePlot(human_model, features=c('SOX2', 'POU5F1', 'NANOG', 'SFRP2',
                                    'GATA6', 'SOX17', 'HNF4A', 'PDGFRA',
                                    'GATA3', 'TFAP2C', 'GATA2', 'KRT7',
                                    'TFAP2A', 'ISL1', 'GABRP', 'VTCN1',
                                    'COL6A2', 'TBXT', 'EOMES', 'SNAI1'), reduction='wnn.umap',
            order=T, min.cutoff = 'q10')

#EXMC gene expression for Figure 3 supplement
FeaturePlot(human_model, features=c('HAND1', 'TBX20', 'FOXF1'),
            order=T, reduction='wnn.umap', min.cutoff = 'q2', max.cutoff = 'q95')

#Nebulosa plot for PGCLCs
SCpubr::do_NebulosaPlot(human_model, features=c('PRDM1', 'NANOS3', 'SOX17', 'TFAP2C', 'NANOG'), joint=T,
                        return_only_joint = T, plot.axes=F, slot='data')

#Expression patterns in GATA3-AP2y cells marked by GFP expression for Figure 3 Supplement
GFP_pos_ids = names(which(human_model@assays$RNA@counts['GFP',]>=1))
DimPlot(human_model, reduction='wnn.umap', cells.highlight = GFP_pos_ids)

human_model@meta.data$GFP<-ifelse(
  rownames(human_model@meta.data) %in% GFP_pos_ids,
  "GFP_positive", "GFP_negative"
)

Idents(human_model)<-human_model$GFP
VlnPlot(human_model, features=c('GATA3', 'TFAP2C', 'rna_GFP', 'CDX2', 'GATA2', 'TEAD4', 'KRT19',
                                'GATA6', 'GATA4', 'SOX17', 'NANOG', 'SOX2'), pt.size=0)

#ID gene expression plots for Figure 4
FeaturePlot(human_model, features=c('ID1', 'ID2', 'ID3', 'ID4'),
            order=T, reduction='wnn.umap', min.cutoff = 'q5', max.cutoff = 'q99')



#Anterior hypoblast marker expression for Figure 5
FeaturePlot(human_model, features=c('CER1'), order=T, reduction='wnn.umap') + xlim(8, 15)+ylim(4,12)
FeaturePlot(human_model, features=c('LEFTY1'), order=T, reduction='wnn.umap') + xlim(8, 15)+ylim(4,12)


##### Motif Accessibility Plots #####
DefaultAssay(human_model)<-'chromvar'

#SMAD5 and SMAD2/3/4 motifs for Figure 4
FeaturePlot(human_model, features=c('MA1557.1'), cols=c('light grey', 'deep pink'), reduction='wnn.umap',
            min.cutoff = 'q1', max.cutoff = 'q90', order=T) #SMAD5 motif activity

FeaturePlot(human_model, features=c('MA0513.1'), cols=c('light grey', 'deep pink'), reduction='wnn.umap',
            min.cutoff = 'q5', max.cutoff = 'q95', order=T) #SMAD2::SMAD3::SMAD4 motif

##### Dimplots #####
#Colored cell assignment plots for Figure 3 + Figure 3 supplement
DimPlot(human_model, reduction='wnn.umap', 
        group.by='cell_assignment', cols=c('gold', 'cyan', 'deepskyblue', 'blue',
                                           'green', 'darkgreen', 'orangered', 'blueviolet'))

DimPlot(human_model, reduction='wnn.umap', 
        group.by='cell_assignment2', cols=c('gold', 'cyan', 'deepskyblue', 'blue',
                                            'green', 'darkgreen', 'orangered', 'blueviolet',
                                           'deeppink'))

DimPlot(human_model, reduction='wnn.umap', 
        group.by='course_cell_assignment', cols=c('gold', 'blue', 'darkgreen', 'orangered', 'blueviolet'))


DimPlot(human_model, reduction='wnn.umap', 
        group.by='scmap_Tyser', cols=c('antiquewhite3', 'aquamarine3', 'bisque2', 'blue', 'blueviolet', 
                                       'brown', 'burlywood4', 'cadetblue4', 'chartreuse4', 'coral2',
                                       'cornflowerblue', 'cyan', 'darkgrey', 'darkgreen', 'darkred'), pt.size=2)

DimPlot(human_model, reduction='wnn.umap', 
        group.by='scmap_ma', cols=c('aquamarine1', 'hotpink', 'firebrick', 
                                    'darkcyan', 'blue', 'magenta', 'coral2',
                                    'gold', 'chartreuse', 'darkgrey', 'darkgreen'))

DimPlot(human_model, reduction='wnn.umap', 
        group.by='scmap_nakamura', cols=c('darkred', 'chartreuse', 'darkgoldenrod4', 'darkolivegreen4', 'darkgreen',
                                          'darkorchid', 'bisque4', 'darkgrey', 'darkgreen'))


DimPlot(human_model, reduction='wnn.umap', 
        group.by='scmapCELL_Yang', cols=c('aquamarine', 'cadetblue4', 'blue', 'darkgreen',
                                          'bisque4', 'darkred', 'coral4', 'chartreuse4', 'darkolivegreen3',
                                          'coral', 'darkorchid', 'darkgrey'))


DimPlot(human_model, reduction='wnn.umap', 
        group.by='scmapCELL_Mole', cols=c('darkorchid', 'bisque4', 'darkgreen', 'darkgoldenrod1', 'darkgrey'))

DimPlot(human_model, reduction='wnn.umap',
        group.by='scmap_Kagawa', cols=c('gold', 'cyan', 'magenta', 'darkgrey'))

DimPlot(human_model, reduction='wnn.umap',
        group.by='scmapCELL_PASE', cols=c('darkorchid1', 'goldenrod1', 'magenta', 'darkolivegreen3', 'darkgreen', 'purple4', 'darkgrey'))

DimPlot(human_model, reduction='wnn.umap',
        group.by='scmapCELL_Pham', cols=c('firebrick', 'goldenrod1', 'cyan', 'darkgrey'))


##### Alluvial Plots ####

do_AlluvialPlot(human_model, "sample_type","cell_assignment", 
                colors.use=c('gold', 'cyan', 'deepskyblue', 'blue',
                             'green', 'darkgreen', 'orangered', 'blueviolet'),
                )

##### Violin Plots #####
DefaultAssay(human_model)<-'RNA'
VlnPlot(human_model, pt.size=0, features=c('BMP4', 'BMP2', 'BMP6', 'BMP7', 'BMP5'), split.by = 'sample_type')

##### switchde/latent time plots for Figure 4 supplement #####
library(switchde)
library(SingleCellExperiment)
library(tidyverse)
library(ggplot2)
library(SeuratWrappers)

multivelo_obs <- read_csv("human_model_analysis/multivelo/adata_results/obs.csv")
rownames(multivelo_obs)<-multivelo_obs$seurat_id
latent_time<-as.data.frame(multivelo_obs$latent_time)
rownames(latent_time)<-rownames(multivelo_obs)
human_model<-AddMetaData(human_model, metadata=latent_time, col.name='latent_time')
#FeaturePlot(human_model, reduction='wnn.umap', features='latent_time', cols=c('blue', 'red'))
gc()
Idents(human_model)<-human_model$cell_assignment

##### gene expression #####

### Amnion ###
AME<-subset(human_model, idents=c('L-EPI', 'AM-1', 'AM-2', 'AM-3'))
DefaultAssay(AME)<-"RNA"
FeaturePlot(AME,reduction='wnn.umap', 'latent_time', cols=c('blue', 'red'))
AME_expr<-AME@assays[["RNA"]]@data
AME_expr<-as.matrix(AME_expr)
AME_filtered <- AME_expr[rowMeans(AME_expr) > 0.1 & rowMeans(AME_expr > 0) > 0.2,]
AME_RNA_out<-switchde(AME_filtered, pseudotime=AME$latent_time)
saveRDS(AME_RNA_out, "human_model_analysis/switchde/AME_RNA_out.RDS")

pars <- extract_pars(AME_RNA_out, "ID1")
switchplot(AME_filtered["ID1",], AME$latent_time, pars) + theme_classic() +
  geom_point(aes(color=factor(AME$cell_assignment2)), alpha=1, ) + scale_color_manual(values=c('gold','cyan', 'deepskyblue', 'blue', 'magenta')) +
  labs(x="latent_time", y="Expression", color="cell_type", title = "ID1")

pars <- extract_pars(AME_RNA_out, "ID2")
switchplot(AME_filtered["ID2",], AME$latent_time, pars) + theme_classic() +
  geom_point(aes(color=factor(AME$cell_assignment2)), alpha=1, ) + scale_color_manual(values=c('gold','cyan', 'deepskyblue', 'blue', 'magenta')) +
  labs(x="latent_time", y="Expression", color="cell_type", title = "ID2")

pars <- extract_pars(AME_RNA_out, "ID3")
switchplot(AME_filtered["ID3",], AME$latent_time, pars) + theme_classic() +
  geom_point(aes(color=factor(AME$cell_assignment2)), alpha=1, ) + scale_color_manual(values=c('gold','cyan', 'deepskyblue', 'blue', 'magenta')) +
  labs(x="latent_time", y="Expression", color="cell_type", title = "ID3")

pars <- extract_pars(AME_RNA_out, "ID4")
switchplot(AME_filtered["ID4",], AME$latent_time, pars) + theme_classic() +
  geom_point(aes(color=factor(AME$cell_assignment2)), alpha=1, ) + scale_color_manual(values=c('gold','cyan', 'deepskyblue', 'blue' , 'magenta')) +
  labs(x="latent_time", y="Expression", color="cell_type", title = "ID4")

#Chromvar
DefaultAssay(AME)<-"chromvar"
AME_expr<-AME@assays[["chromvar"]]@data
AME_expr<-as.matrix(AME_expr)
AME_filtered <- AME_expr
AME_chromvar_out<-switchde(AME_filtered, pseudotime=AME$latent_time)
saveRDS(AME_chromvar_out, "human_model_analysis/switchde/AME_chromvar_out.RDS")

pars <- extract_pars(AME_chromvar_out, "MA1557.1")
switchplot(AME_filtered["MA1557.1",], AME$latent_time, pars) + theme_classic() +
  geom_point(aes(color=factor(AME$cell_assignment2)), alpha=1, ) + scale_color_manual(values=c('gold','cyan', 'deepskyblue', 'blue', 'magenta')) +
  labs(x="latent_time", y="Relative Accessibility", color="cell_type", title = "SMAD5 Motif")

pars <- extract_pars(AME_chromvar_out, "MA0513.1")
switchplot(AME_filtered["MA0513.1",], AME$latent_time, pars) + theme_classic() +
  geom_point(aes(color=factor(AME$cell_assignment2)), alpha=1, ) + scale_color_manual(values=c('gold','cyan', 'deepskyblue', 'blue', 'magenta')) +
  labs(x="latent_time", y="Relative Accessibility", color="cell_type", title = "SMAD2::3::4 Motif")

### EXMC ###
EXMC<-subset(human_model, idents=c('L-EPI', 'MESO-1', 'MESO-2', 'EXMC'))
DefaultAssay(EXMC)<-"RNA"
FeaturePlot(EXMC,reduction='wnn.umap', 'latent_time', cols=c('blue', 'red'))
EXMC_expr<-EXMC@assays[["RNA"]]@data
EXMC_expr<-as.matrix(EXMC_expr)
EXMC_filtered <- EXMC_expr[rowMeans(EXMC_expr) > 0.1 & rowMeans(EXMC_expr > 0) > 0.1,]
EXMC_RNA_out<-switchde(EXMC_filtered, pseudotime=EXMC$latent_time)
saveRDS(EXMC_RNA_out, "human_model_analysis/switchde/EXMC_RNA_out.RDS")

pars <- extract_pars(EXMC_RNA_out, "ID1")
switchplot(EXMC_filtered["ID1",], EXMC$latent_time, pars) + theme_classic() +
  geom_point(aes(color=factor(EXMC$cell_assignment2)), alpha=1, ) + scale_color_manual(values=c('gold','green', 'darkgreen', 'orangered', 'magenta')) +
  labs(x="latent_time", y="Expression", color="cell_type", title = "ID1")

pars <- extract_pars(EXMC_RNA_out, "ID2")
switchplot(EXMC_filtered["ID2",], EXMC$latent_time, pars) + theme_classic() +
  geom_point(aes(color=factor(EXMC$cell_assignment2)), alpha=1, ) + scale_color_manual(values=c('gold','green', 'darkgreen', 'orangered', 'magenta')) +
  labs(x="latent_time", y="Expression", color="cell_type", title = "ID2")

pars <- extract_pars(EXMC_RNA_out, "ID3")
switchplot(EXMC_filtered["ID3",], EXMC$latent_time, pars) + theme_classic() +
  geom_point(aes(color=factor(EXMC$cell_assignment2)), alpha=1, ) + scale_color_manual(values=c('gold','green', 'darkgreen', 'orangered', 'magenta')) +
  labs(x="latent_time", y="Expression", color="cell_type", title = "ID3")

pars <- extract_pars(EXMC_RNA_out, "ID4")
switchplot(EXMC_filtered["ID4",], EXMC$latent_time, pars) + theme_classic() +
  geom_point(aes(color=factor(EXMC$cell_assignment2)), alpha=1, ) + scale_color_manual(values=c('gold','green', 'darkgreen', 'orangered', 'magenta')) +
  labs(x="latent_time", y="Expression", color="cell_type", title = "ID4")

#Chromvar
DefaultAssay(EXMC)<-"chromvar"
EXMC_expr<-EXMC@assays[["chromvar"]]@data
EXMC_expr<-as.matrix(EXMC_expr)
EXMC_filtered <- EXMC_expr
EXMC_chromvar_out<-switchde(EXMC_filtered, pseudotime=EXMC$latent_time)
saveRDS(EXMC_chromvar_out, "human_model_analysis/switchde/EXMC_chromvar_out.RDS")

pars <- extract_pars(HYPO_chromvar_out, "MA1557.1")
switchplot(EXMC_filtered["MA1557.1",], EXMC$latent_time, pars) + theme_classic() +
  geom_point(aes(color=factor(EXMC$cell_assignment2)), alpha=1, ) + scale_color_manual(values=c('gold','green', 'darkgreen', 'orangered', 'magenta')) +
  labs(x="latent_time", y="Relative Accessibility", color="cell_type", title = "SMAD5 Motif")

pars <- extract_pars(EXMC_chromvar_out, "MA0513.1")
switchplot(EXMC_filtered["MA0513.1",], EXMC$latent_time, pars) + theme_classic() +
  geom_point(aes(color=factor(EXMC$cell_assignment2)), alpha=1, ) + scale_color_manual(values=c('gold','green', 'darkgreen', 'orangered', 'magenta')) +
  labs(x="latent_time", y="Relative Accessibility", color="cell_type", title = "SMAD2::3::4 Motif")


###### Generating input files for cellphoneDB #####

matrix<-as.matrix(human_model@assays[["RNA"]]@data)
meta<-cbind(human_model@meta.data[,'cell_assignment', drop=F])
meta2<-cbind(human_model@meta.data[,'cell_assignment2', drop=F])
write.table(matrix, './human_model_analysis/cellphonedb/matrix.txt', sep='\t', quote=F)
write.table(meta, './human_model_analysis/cellphonedb/meta.txt', sep='\t', quote=F)
write.table(meta, './human_model_analysis/cellphonedb/meta2.txt', sep='\t', quote=F)

matrix<-as.matrix(cell_lines@assays[["RNA"]]@data)
meta<-cbind(cell_lines@meta.data[,'sample_type', drop=F])
write.table(matrix, './cell_lines_analysis/cellphonedb/matrix.txt', sep='\t', quote=F)
write.table(meta, './cell_lines_analysis/cellphonedb/meta.txt', sep='\t', quote=F)


#In text editor, add column name "cell" to meta file
#in ubuntu window, install cellphonedb and R4.0
#in ubuntu window, change to directory with counts & run
#cellphonedb method statistical_analysis meta_file.txt matrix_file.txt --iterations=10 --threads=2 --counts-data=gene_name
#cellphonedb plot dot_plot --means-path=./out/means.txt --pvalues-path=./out/pvalues.txt --rows rows_to_plot.txt