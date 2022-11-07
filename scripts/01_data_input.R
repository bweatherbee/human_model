#Following the vignette at: https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis.html#wnn-analysis-of-10x-multiome-rna-atac-1

library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(ggplot2)

setwd('D://human_model')

# the 10x hdf5 file contains both data types
inputdata_10x_d4 <- Read10X_h5('raw_outs/exp48_d4/filtered_feature_bc_matrix.h5')
inputdata_10x_d6 <- Read10X_h5('raw_outs/exp48_d6/filtered_feature_bc_matrix.h5')
inputdata_10x_d8 <- Read10X_h5('raw_outs/exp48_d8/filtered_feature_bc_matrix.h5')

inputdata_10x_s6 <- Read10X_h5('raw_outs/shef6/filtered_feature_bc_matrix.h5')
inputdata_10x_g6_s17 <- Read10X_h5('raw_outs/g6_s17/filtered_feature_bc_matrix.h5')
inputdata_10x_g3_ap2y <- Read10X_h5('raw_outs/g3_ap2y/filtered_feature_bc_matrix.h5')

# extract RNA and ATAC data 
rna_counts_d4 <- inputdata_10x_d4$`Gene Expression`
atac_counts_d4 <- inputdata_10x_d4$Peaks

rna_counts_d6 <- inputdata_10x_d6$`Gene Expression`
atac_counts_d6 <- inputdata_10x_d6$Peaks

rna_counts_d8 <- inputdata_10x_d8$`Gene Expression`
atac_counts_d8 <- inputdata_10x_d8$Peaks

rna_counts_shef6 <- inputdata_10x_s6$`Gene Expression`
atac_counts_shef6 <- inputdata_10x_s6$Peaks

rna_counts_g6_s17 <- inputdata_10x_g6_s17$`Gene Expression`
atac_counts_g6_s17 <- inputdata_10x_g6_s17$Peaks

rna_counts_g3_ap2y <- inputdata_10x_g3_ap2y$`Gene Expression`
atac_counts_g3_ap2y <- inputdata_10x_g3_ap2y$Peaks

# Create Seurat object - one for each sample. RNA first, plus calculation of mitocondrial content.
d4_seurat <- CreateSeuratObject(counts = rna_counts_d4)
d4_seurat[["percent.mt"]] <- PercentageFeatureSet(d4_seurat, pattern = "^MT-")

d6_seurat <- CreateSeuratObject(counts = rna_counts_d6)
d6_seurat[["percent.mt"]] <- PercentageFeatureSet(d6_seurat, pattern = "^MT-")

d8_seurat <- CreateSeuratObject(counts = rna_counts_d8)
d8_seurat[["percent.mt"]] <- PercentageFeatureSet(d8_seurat, pattern = "^MT-")

s6_seurat <- CreateSeuratObject(counts = rna_counts_shef6)
s6_seurat[["percent.mt"]] <- PercentageFeatureSet(s6_seurat, pattern = "^MT-")

g6_s17_seurat <- CreateSeuratObject(counts = rna_counts_g6_s17)
g6_s17_seurat[["percent.mt"]] <- PercentageFeatureSet(g6_s17_seurat, pattern = "^MT-")

g3_ap2y_seurat <- CreateSeuratObject(counts = rna_counts_g3_ap2y)
g3_ap2y_seurat[["percent.mt"]] <- PercentageFeatureSet(g3_ap2y_seurat, pattern = "^MT-")

# Now add in the ATAC-seq data - one for each time point.
# we'll only use peaks in standard chromosomes
# day 4 structures from exp 48
grange_counts_d4 <- StringToGRanges(rownames(atac_counts_d4), sep = c(":", "-"))
grange_use_d4 <- seqnames(grange_counts_d4) %in% standardChromosomes(grange_counts_d4)
atac_counts_d4 <- atac_counts_d4[as.vector(grange_use_d4), ]
annotations_d4 <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations_d4) <- 'UCSC'
genome(annotations_d4) <- "hg38"

frag_file_d4 <- "raw_outs/exp48_d4/atac_fragments.tsv.gz"
chrom_assay_d4 <- CreateChromatinAssay(
  counts = atac_counts_d4,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = frag_file_d4,
  min.cells = 10,
  annotation = annotations_d4
)
d4_seurat[["ATAC"]] <- chrom_assay_d4

# d6 structures from exp 48
grange_counts_d6 <- StringToGRanges(rownames(atac_counts_d6), sep = c(":", "-"))
grange_use_d6 <- seqnames(grange_counts_d6) %in% standardChromosomes(grange_counts_d6)
atac_counts_d6 <- atac_counts_d6[as.vector(grange_use_d6), ]
annotations_d6 <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations_d6) <- 'UCSC'
genome(annotations_d6) <- "hg38"

frag_file_d6 <- "raw_outs/exp48_d6/atac_fragments.tsv.gz"
chrom_assay_d6 <- CreateChromatinAssay(
  counts = atac_counts_d6,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = frag_file_d6,
  min.cells = 10,
  annotation = annotations_d6
)
d6_seurat[["ATAC"]] <- chrom_assay_d6

# day 8 structures from exp 48
grange_counts_d8 <- StringToGRanges(rownames(atac_counts_d8), sep = c(":", "-"))
grange_use_d8 <- seqnames(grange_counts_d8) %in% standardChromosomes(grange_counts_d8)
atac_counts_d8 <- atac_counts_d8[as.vector(grange_use_d8), ]
annotations_d8 <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations_d8) <- 'UCSC'
genome(annotations_d8) <- "hg38"

frag_file_d8 <- "raw_outs/exp48_d8/atac_fragments.tsv.gz"
chrom_assay_d8 <- CreateChromatinAssay(
  counts = atac_counts_d8,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = frag_file_d8,
  min.cells = 10,
  annotation = annotations_d8
)
d8_seurat[["ATAC"]] <- chrom_assay_d8

# shef6 wt RSeT hESCs

grange_counts_s6 <- StringToGRanges(rownames(atac_counts_shef6), sep = c(":", "-"))
grange_use_s6 <- seqnames(grange_counts_s6) %in% standardChromosomes(grange_counts_s6)
atac_counts_s6 <- atac_counts_shef6[as.vector(grange_use_s6), ]
annotations_s6 <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations_s6) <- 'UCSC'
genome(annotations_s6) <- "hg38"

frag_file_s6 <- "raw_outs/shef6/atac_fragments.tsv.gz"
chrom_assay_s6 <- CreateChromatinAssay(
  counts = atac_counts_shef6,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = frag_file_s6,
  min.cells = 10,
  annotation = annotations_s6
)
s6_seurat[["ATAC"]] <- chrom_assay_s6

# 3 day induced GATA6-SOX17 overexpression in RSeT hESCs

grange_counts_g6_s17 <- StringToGRanges(rownames(atac_counts_g6_s17), sep = c(":", "-"))
grange_use_g6_s17 <- seqnames(grange_counts_g6_s17) %in% standardChromosomes(grange_counts_g6_s17)
atac_counts_g6_s17 <- atac_counts_g6_s17[as.vector(grange_use_g6_s17), ]
annotations_g6_s17 <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations_g6_s17) <- 'UCSC'
genome(annotations_g6_s17) <- "hg38"

frag_file_g6_s17 <- "raw_outs/g6_s17/atac_fragments.tsv.gz"
chrom_assay_g6_s17 <- CreateChromatinAssay(
  counts = atac_counts_g6_s17,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = frag_file_g6_s17,
  min.cells = 10,
  annotation = annotations_g6_s17
)
g6_s17_seurat[["ATAC"]] <- chrom_assay_g6_s17


#g3_ap2y 3day overexpression in RSeT hESCs
grange_counts_g3_ap2y <- StringToGRanges(rownames(atac_counts_g3_ap2y), sep = c(":", "-"))
grange_use_g3_ap2y <- seqnames(grange_counts_g3_ap2y) %in% standardChromosomes(grange_counts_g3_ap2y)
atac_counts_g3_ap2y <- atac_counts_g3_ap2y[as.vector(grange_use_g3_ap2y), ]
annotations_g3_ap2y <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations_g3_ap2y) <- 'UCSC'
genome(annotations_g3_ap2y) <- "hg38"

frag_file_g3_ap2y <- "raw_outs/g3_ap2y/atac_fragments.tsv.gz"
chrom_assay_g3_ap2y <- CreateChromatinAssay(
  counts = atac_counts_g3_ap2y,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = frag_file_g3_ap2y,
  min.cells = 10,
  annotation = annotations_g3_ap2y
)
g3_ap2y_seurat[["ATAC"]] <- chrom_assay_g3_ap2y


# merge unnormalised seurat objects
# add sample metadata
d4_seurat[['sample_type']] <- 'double_structure_day4'
d6_seurat[['sample_type']] <- 'double_structure_day6'
d8_seurat[['sample_type']] <- 'double_structure_day8'
s6_seurat[['sample_type']] <- 'RSeT_wt'
g6_s17_seurat[['sample_type']] <- 'g6_s17_RSeT_induced'
g3_ap2y_seurat[['sample_type']] <- 'g3_ap2y_RSeT_induced'

# save unfiltered, time-separated objects
saveRDS(d4_seurat, file = "./intermediary_objects/d4_seurat_noQC.rds")
saveRDS(d6_seurat, file = "./intermediary_objects/d6_seurat_noQC.rds")
saveRDS(d8_seurat, file = "./intermediary_objects/d8_seurat_noQC.rds")
saveRDS(s6_seurat, file = "./intermediary_objects/s6_seurat_noQC.rds")
saveRDS(g6_s17_seurat, file = "./intermediary_objects/g6_s17_seurat_noQC.rds")
saveRDS(g3_ap2y_seurat, file = "./intermediary_objects/g3_ap2y_seurat_noQC.rds")


# merge and add cell id prefixes.
merge_seurat <- merge(d4_seurat, y=c(d6_seurat, d8_seurat, s6_seurat, g6_s17_seurat, g3_ap2y_seurat), add.cell.ids = c('d4_double', 'd6_double', 'd8_double', 'WT', 'g6_s17', 'g3_ap2y'), project = 'human_model_multiome')


saveRDS(merge_seurat, file = "./intermediary_objects/merge_seurat_noQC.rds")





