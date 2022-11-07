#largely following Joint RNA and ATAC: 10x multiomic vignette from Sajita Lab#
#https://satijalab.org/signac/articles/pbmc_multiomic.html#

library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(dplyr)
library(tidyverse)

set.seed(1234)
setwd('D:/human_model/')

### PEAK CALLING ###
#did this in linux because macs2 not supported in windows

### In Ubuntu terminal ###

#R
#reticulate::conda_create('reticulate')
#library(reticulate)
#conda_list()
#use_condaenv('reticulate')
#py_install('macs2', pip=T)
#library(Seurat)
#library(Signac)
#merged<-readRDS('/mnt/d/human_model/intermediary_objects/merge_seurat_noQC.rds')
#DefaultAssay(merged)<-'ATAC'
#merged@assays[["ATAC"]]@fragments[[1]]@path<-'/mnt/d/human_model/raw_outs/exp48_d4/atac_fragments.tsv.gz'
#merged@assays[["ATAC"]]@fragments[[2]]@path<-'/mnt/d/human_model/raw_outs/exp48_d6/atac_fragments.tsv.gz'
#merged@assays[["ATAC"]]@fragments[[3]]@path<-'/mnt/d/human_model/raw_outs/exp48_d8/atac_fragments.tsv.gz'
#merged@assays[["ATAC"]]@fragments[[4]]@path<-'/mnt/d/human_model/raw_outs/shef6/atac_fragments.tsv.gz'
#merged@assays[["ATAC"]]@fragments[[5]]@path<-'/mnt/d/human_model/raw_outs/g6_s17/atac_fragments.tsv.gz'
#merged@assays[["ATAC"]]@fragments[[6]]@path<-'/mnt/d/human_model/raw_outs/g3_ap2y/atac_fragments.tsv.gz'
#peaks<-CallPeaks(merged, macs2.path="/home/bailey/miniconda3/envs/reticulate/bin/macs2")
#saveRDS(peaks, '/mnt/d/human_model/intermediary_objects/macs2_peaks_GRangesobj.rds')


peaks <- readRDS("D:/human_model/intermediary_objects/macs2_peaks_GRangesobj.rds")
merged <- readRDS("D:/human_model/intermediary_objects/merge_seurat_noQC.rds")
peaks<-keepStandardChromosomes(peaks, pruning.mode='coarse')
peaks<-subsetByOverlaps(x=peaks, ranges=blacklist_hg38_unified, invert=T)
DefaultAssay(merged)<-'ATAC'
annotation<-GetGRangesFromEnsDb(ensdb=EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- "UCSC"
genome(annotation)<-"hg38"

macs2_counts<-FeatureMatrix(fragments=Fragments(merged), features=peaks, cells=colnames(merged))
merged[["peaks"]]<-CreateChromatinAssay(counts=macs2_counts, fragments=merged@assays[["ATAC"]]@fragments, annotation=annotation)


### QC and Normalizations ###
merged<-NucleosomeSignal(merged)
merged<-TSSEnrichment(merged, fast=F)
VlnPlot(merged, features=c('nCount_RNA', 'nCount_ATAC', 'TSS.enrichment', 'nucleosome_signal'),
        ncol=4, pt.size=0)
summary(merged$TSS.enrichment)
summary(merged$nCount_RNA)
summary(merged$nCount_ATAC)
saveRDS(merged, './intermediary_objects/merged_withpeaks_unfiltered.rds')



###Weight Nearest Neighbor Analysis###
#based on 10x multiome vignette in WNNA in Seurat
#https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis.html
merged_filtered<-subset(merged, subset=nCount_RNA>500 & nCount_ATAC>500 & percent.mt<20 & TSS.enrichment>1)
DefaultAssay(merged_filtered)<-'RNA'

### Doublet calling ###
library(scDblFinder)
library(data.table)
set.seed(1234)
sce<-as.SingleCellExperiment(merged_filtered, assay = 'RNA')
sce<-scDblFinder(sce, samples='sample_type')
doublet_table<-as.data.frame(sce$scDblFinder.class)
table(truth=sce$sample_type, call=sce$scDblFinder.class)
merged_filtered$doublet<-sce$scDblFinder.class
merged_filtered<-subset(merged_filtered, subset=doublet=='singlet')
saveRDS(merged_filtered, './intermediary_objects/merged_filtered.rds')

##### POST-IMPLANTATION EMBRYO MODELS #####
### Subsetting to only structures ###
Idents(merged_filtered)<-merged_filtered$sample_type
human_model<-subset(merged_filtered, idents=c('double_structure_day4', 'double_structure_day6', 'double_structure_day8'))

#Cell Cycle Scoring
s.genes<-cc.genes$s.genes
g2m.genes<-cc.genes$g2m.genes
human_model<-CellCycleScoring(human_model, s.features = s.genes, g2m.features = g2m.genes, set.ident=F)

#Removing mitochondrial, cell cycle and long non-coding RNA genes from assay to clean up downstream analysis 
human_model@assays[["RNA"]]@counts<-human_model@assays[["RNA"]]@counts[-(which(str_detect(rownames(human_model), "^HIST|^MT-|^TOP|^CDK|^CCN|^CDC|^CCDC|^MKI|^MALAT|^NUSAP|^SMC|^CENP|^UBE|^SGO|^ASPM|^PLK|^KPN|^RP|^PTTG|^SNHG|^CK|^AUR^|^BUB|^KIF|^KCNQ|^SMO|^HMG|^S100|^LINC|^ATP"))),]
human_model@assays[["RNA"]]@data<-human_model@assays[["RNA"]]@data[-(which(str_detect(rownames(human_model), "^HIST|^MT-|^TOP|^CDK|^CCN|^CDC|^CCDC|^MKI|^MALAT|^NUSAP|^SMC|^CENP|^UBE|^SGO|^ASPM|^PLK|^KPN|^RP|^PTTG|^SNHG|^CK|^AUR^|^BUB|^KIF|^KCNQ|^SMO|^HMG|^S100|^LINC|^ATP"))),]

#Using SCTransform Normalization with no batch correction followed by RNA-based UMAP
human_model<-SCTransform(human_model, vars.to.regress=c('percent.mt', 'S.Score', 'G2M.Score'), verbose=F)
DefaultAssay(human_model)<-'SCT'
human_model<-RunPCA(human_model)
ElbowPlot(human_model)
human_model<-RunUMAP(human_model, dims=1:15, reduction.name='umap.rna', reduction.key='rnaUMAP_')
DimPlot(human_model, reduction='umap.rna', group.by='sample_type')
FeaturePlot(human_model, features=c('POU5F1', 'NANOG', 'SOX2', 'BCL11A',
                                    'CDX2', 'ISL1', 'TFAP2A', 'VTCN1',
                                    'GATA6', 'HNF4A', 'SOX17', 'CDH2',
                                    'GATA3', 'TFAP2C', 'GATA2', 'GFP',
                                    'TBXT', 'EOMES', 'SNAI1', 'CDH1'), order=T, min.cutoff = 'q25', max.cutoff = 'q90')


#Now projecting UMAP from peaks ATAC assay
DefaultAssay(human_model)<-'peaks'
human_model<-FindTopFeatures(human_model, min.cutoff='q0')
human_model<-RunTFIDF(human_model)
human_model<-RunSVD(human_model)
human_model<-RunUMAP(human_model, reduction='lsi', dims=2:15, reduction.name='umap.atac', reduction.key='atacUMAP_')
DimPlot(human_model, reduction='umap.atac', group.by='sample_type')


### Joint UMAP Visualization and Clustering ###
human_model<-FindMultiModalNeighbors(human_model,
                                         reduction.list=c('pca', 'lsi'),
                                         dims.list=list(1:10, 2:10),
                                         modality.weight.name = c('SCT.weight', 'peaks.weight'),
                                         verbose=T)

human_model <- RunUMAP(human_model, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_", dims.list = list(1:10, 2:10))
human_model <- FindClusters(human_model, graph.name = "wsnn", algorithm = 3, verbose = FALSE, resolution = 0.1)
DimPlot(human_model, reduction='wnn.umap')
DimPlot(human_model, group.by = 'sample_type', reduction='wnn.umap')
DefaultAssay(human_model)<-'SCT'
FeaturePlot(human_model, features=c('POU5F1', 'NANOG', 'SOX2', 'BCL11A',
                                    'CDX2', 'ISL1', 'TFAP2A', 'VTCN1',
                                    'GATA6', 'HNF4A', 'SOX17', 'CDH2',
                                    'GATA3', 'TFAP2C', 'GATA2', 'GFP',
                                    'TBXT', 'EOMES', 'SNAI1', 'CDH1'), reduction='wnn.umap', order=T, min.cutoff = 'q10', max.cutoff = 'q95')

DimPlot(human_model, reduction='umap.rna', group.by = 'sample_type') + DimPlot(human_model, reduction='umap.atac', group.by='sample_type') + DimPlot(human_model, reduction='wnn.umap', group.by='sample_type')


DimPlot(human_model, reduction='umap.rna') + DimPlot(human_model, reduction='umap.atac') + DimPlot(human_model, reduction='wnn.umap')
DimPlot(human_model, reduction='wnn.umap')
DimPlot(human_model, reduction='wnn.umap', group.by = 'sample_type')
saveRDS(human_model, './intermediary_objects/human_model.rds')

### Linking Peaks to Genes ###

DefaultAssay(human_model) <- "peaks"

# first compute the GC content for each peak
human_model <- RegionStats(human_model, genome = BSgenome.Hsapiens.UCSC.hg38)

# link peaks to genes - this can take a long time
human_model <- LinkPeaks(
  object = human_model,
  peak.assay = "peaks",
  expression.assay = "SCT",
  genes.use = NULL
)

saveRDS(human_model, './intermediary_objects/human_model_linkedpeaks.rds')

### QC Plots (based on pmbc signac vignette) ###
DefaultAssay(human_model)<-'ATAC'
TSSPlot(human_model)
FragmentHistogram(human_model)


### Adding Motif information to object based on Signac Motif Vignette ###
DefaultAssay(human_model)<-'peaks'

library(JASPAR2020)
library(TFBSTools)

#get a list of motif position frequencies matrices from JASPAR database
pfm<-getMatrixSet(x=JASPAR2020, opts=list(collection='CORE', tax_group='vertebrates', all_version=F))

#add motif info
human_model<-AddMotifs(human_model, genome=BSgenome.Hsapiens.UCSC.hg38, pfm)
rm(pfm)
gc()

#computing motif accessibility scores using chromVAR per each single cell (creates new assay)
human_model<-RunChromVAR(human_model, genome=BSgenome.Hsapiens.UCSC.hg38)
saveRDS(human_model, 'human_model_linkedpeaks.rds')


### Scaling Data and saving ###
DefaultAssay(human_model)<-'RNA'
human_model<-ScaleData(human_model, vars.to.regress=c('nCount_RNA', 'percent.mt', 'S.Score', 'G2M.Score'))
human_model<-NormalizeData(human_model)
DefaultAssay(human_model)<-'peaks'
human_model<-ScaleData(human_model)
saveRDS(human_model, 'human_model_complete.rds')



##### Cell Lines Only ####
Idents(merged_filtered)<-merged_filtered$sample_type
cell_lines<-subset(merged_filtered, idents=c('RSeT_wt', 'g6_s17_RSeT_induced', 'g3_ap2y_RSeT_induced'))


#Cell Cycle Scoring
s.genes<-cc.genes$s.genes
g2m.genes<-cc.genes$g2m.genes
cell_lines<-CellCycleScoring(cell_lines, s.features = s.genes, g2m.features = g2m.genes, set.ident=F)

#Removing mitochondrial, cell cycle and long non-coding RNA genes 
cell_lines@assays[["RNA"]]@counts<-cell_lines@assays[["RNA"]]@counts[-(which(str_detect(rownames(cell_lines), "^HIST|^MT-|^TOP|^CDK|^CCN|^CDC|^CCDC|^MKI|^MALAT|^NUSAP|^SMC|^CENP|^UBE|^SGO|^ASPM|^PLK|^KPN|^RP|^PTTG|^SNHG|^CK|^AUR^|^BUB|^KIF|^KCNQ|^SMO|^HMG|^S100|^LINC|^ATP"))),]
cell_lines@assays[["RNA"]]@data<-cell_lines@assays[["RNA"]]@data[-(which(str_detect(rownames(cell_lines), "^HIST|^MT-|^TOP|^CDK|^CCN|^CDC|^CCDC|^MKI|^MALAT|^NUSAP|^SMC|^CENP|^UBE|^SGO|^ASPM|^PLK|^KPN|^RP|^PTTG|^SNHG|^CK|^AUR^|^BUB|^KIF|^KCNQ|^SMO|^HMG|^S100|^LINC|^ATP"))),]

#Using SCTransform Normalization with no batch correction
cell_lines<-SCTransform(cell_lines, vars.to.regress=c('percent.mt', 'S.Score', 'G2M.Score'), verbose=F)
DefaultAssay(cell_lines)<-'SCT'
cell_lines<-RunPCA(cell_lines)
ElbowPlot(cell_lines)
cell_lines<-RunUMAP(cell_lines, dims=1:5, reduction.name='umap.rna', reduction.key='rnaUMAP_')
DimPlot(cell_lines, reduction='umap.rna', group.by='sample_type')
DefaultAssay(cell_lines)<-'RNA'
FeaturePlot(cell_lines, features=c('POU5F1', 'NANOG', 'SOX2', 'DPPA5',
                                    'CDX2', 'ISL1', 'TFAP2A', 'VTCN1',
                                    'GATA6', 'HNF4A', 'SOX17', 'CDH2',
                                    'GATA3', 'TFAP2C', 'GATA2', 'GFP',
                                    'TBXT', 'EOMES', 'SNAI1', 'CDH1'), order=T, min.cutoff = 'q25', max.cutoff = 'q90')


#Now projecting UMAP from peaks ATAC assay
DefaultAssay(cell_lines)<-'peaks'
cell_lines<-FindTopFeatures(cell_lines, min.cutoff='q0')
cell_lines<-RunTFIDF(cell_lines)
cell_lines<-RunSVD(cell_lines)
cell_lines<-RunUMAP(cell_lines, reduction='lsi', dims=2:5, reduction.name='umap.atac', reduction.key='atacUMAP_')
DimPlot(cell_lines, reduction='umap.atac', group.by='sample_type')


### Joint UMAP Visualization and Clustering ###
cell_lines<-FindMultiModalNeighbors(cell_lines,
                                     reduction.list=c('pca', 'lsi'),
                                     dims.list=list(1:5, 2:5),
                                     modality.weight.name = c('SCT.weight', 'peaks.weight'),
                                     verbose=T)

cell_lines <- RunUMAP(cell_lines, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_", dims.list = list(1:5, 2:5))
DimPlot(cell_lines, reduction='wnn.umap')
DimPlot(cell_lines, group.by = 'sample_type', reduction='wnn.umap')
Idents(cell_lines)<-cell_lines$sample_type
DefaultAssay(cell_lines)<-'RNA'
FeaturePlot(cell_lines, features=c('POU5F1', 'NANOG', 'SOX2', 'DPPA5',
                                    'CDX2', 'ISL1', 'TFAP2A', 'VTCN1',
                                    'GATA6', 'HNF4A', 'SOX17', 'CDH2',
                                    'GATA3', 'TFAP2C', 'GATA2', 'GFP'), reduction='wnn.umap', order=T, min.cutoff = 'q10', max.cutoff = 'q95')

DimPlot(cell_lines, reduction='umap.rna', group.by = 'sample_type') + DimPlot(cell_lines, reduction='umap.atac', group.by='sample_type') + DimPlot(cell_lines, reduction='wnn.umap', group.by='sample_type')


saveRDS(cell_lines, './intermediary_objects/cell_lines.rds')

### Linking Peaks to Genes ###

DefaultAssay(cell_lines) <- "peaks"

# first compute the GC content for each peak
cell_lines <- RegionStats(cell_lines, genome = BSgenome.Hsapiens.UCSC.hg38)

# link peaks to genes
cell_lines <- LinkPeaks(
  object = cell_lines,
  peak.assay = "peaks",
  expression.assay = "SCT",
  genes.use = NULL
)

saveRDS(cell_lines, './intermediary_objects/cell_lines_linkedpeaks.rds')

### QC Plots (based on pmbc signac vignette) ###
DefaultAssay(cell_lines)<-'ATAC'
TSSPlot(cell_lines)
FragmentHistogram(cell_lines)

### Adding Motif information to object based on Signac Motif Vignette ###
DefaultAssay(cell_lines)<-'peaks'

library(JASPAR2020)
library(TFBSTools)

#get a list of motif position frequencies matrices from JASPAR database
pfm<-getMatrixSet(x=JASPAR2020, opts=list(collection='CORE', tax_group='vertebrates', all_version=F))

#add motif info
cell_lines<-AddMotifs(cell_lines, genome=BSgenome.Hsapiens.UCSC.hg38, pfm)
rm(pfm)
gc()

#computing motif accessibility scores using chromVAR per each single cell (creates new assay)
cell_lines<-RunChromVAR(cell_lines, genome=BSgenome.Hsapiens.UCSC.hg38)
saveRDS(cell_lines, 'cell_lines_linkedpeaks.rds')

### Scaling data and saving ###
DefaultAssay(cell_lines)<-'RNA'
cell_lines<-ScaleData(cell_lines, vars.to.regress=c('nCount_RNA', 'percent.mt', 'S.Score', 'G2M.Score'))
cell_lines<-NormalizeData(cell_lines)
DefaultAssay(cell_lines)<-'peaks'
cell_lines<-ScaleData(cell_lines)
saveRDS(cell_lines, 'cell_lines_complete.rds')
