library(Seurat)
library(Signac)
library(RColorBrewer)
library(ggplot2)

human_model <- readRDS("D:/human_model/human_model.rds")
cell_lines <- readRDS("D:/human_model/cell_lines_complete.rds")
DefaultAssay(human_model)<-'RNA'

##### Gene Expression Plots ####
#colored dimplot for cell lines Figure 1
DimPlot(cell_lines, reduction='umap.rna', cols=c('chartreuse2', 'darkorange', '#9d388b'))

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
                                       'cornflowerblue', 'cyan', 'darkgrey', 'darkgreen', 'darkred'))

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