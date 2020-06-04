library(Seurat)
library(Signac)
library(GenomeInfoDb)
library(tidyverse)

################################################################################
# Step 01: Load 10X aggregated peaks by cells matrix into seurat
################################################################################

peaks <- Read10X_h5("cellranger_aggr_dir/outs/filtered_peak_bc_matrix.h5")
cellranger_meta <-  read.csv("cellranger_aggr/outs/singlecell.csv", header=T, row.names=1)
NucSeq.atac <-  CreateSeuratObject(
  counts=peaks,
  assay="peaks",
  meta.data=cellranger_meta,
  min.cells=1
)
fragment.path <- "cellranger_aggr/outs/fragments.tsv.gz"
NucSeq.atac <- SetFragments(NucSeq.atac, file = fragment.path)

# binarize peaks
NucSeq.atac@assays$peaks@counts@x[NucSeq.atac@assays$peaks@counts@x > 0] <- 1

# add sample level metadata (Diagnosis, Age, Sex, etc)
my.data <- read.csv(file = "data/CellSampleID.csv", header=TRUE, row.names=1)
all.equal(rownames(my.data), colnames(NucSeq.atac))
NucSeq.atac <- AddMetaData(NucSeq.atac, metadata=my.data)

################################################################################
# Step 02: Quality Control
################################################################################

NucSeq.atac$pct_reads_in_peaks <- NucSeq.atac$peak_region_fragments / NucSeq.atac$nCount_peaks * 100
NucSeq.atac$blacklist_ratio <- NucSeq.atac$blacklist_region_fragments / NucSeq.atac$peak_region_fragments

# for each sample, compute nucleosome banding
# note: download chromosome sizes file from UCSC https://hgdownload-test.gi.ucsc.edu/goldenPath/hg38/bigZips/
chrom.sizes = read.csv('hg38.chrom.sizes', sep='\t', header=F)
colnames(chrom.sizes) <- c('chr', 'size')
NucSeq.atac$nucleosome_signal <- NA
NucSeq.atac$nucleosome_group <- NA
samples <- unique(NucSeq.atac$Sample.ID)
for(i in 1:length(samples)){
  print(samples[i])
  temp <- NucleosomeSignal(
    subset(NucSeq.atac, Sample.ID == samples[i]),
    region=paste0(chrom.sizes[1,][[1]],"-1-",chrom.sizes[1,][[2]])
  )
  temp$nucleosome_group <- ifelse(temp$nucleosome_signal > 10, 'NS > 10', 'NS < 10')
  NucSeq.atac$nucleosome_signal <- ifelse(NucSeq.atac$Sample.ID == samples[i], temp$nucleosome_signal, NucSeq.atac$nucleosome_signal)
  NucSeq.atac$nucleosome_group <- ifelse(NucSeq.atac$Sample.ID == samples[i], temp$nucleosome_group, NucSeq.atac$nucleosome_group)
}
NucSeq.atac$pass_qc <- ifelse(NucSeq.atac$peak_region_fragments > 300 & NucSeq.atac$peak_region_fragments < 10000 & NucSeq.atac$pct_reads_in_peaks > 15 & NucSeq.atac$blacklist_ratio < 0.01 & NucSeq.atac$nucleosome_signal < 10, TRUE, FALSE)
NucSeq.atac <- NucSeq.atac[,NucSeq.atac$pass_qc]

################################################################################
# Step 03: Primary Processing
################################################################################

# option 1: Process with Seurat / Signac (not used)
NucSeq.atac <- RunTFIDF(NucSeq.atac)
NucSeq.atac <- RunSVD(
  object = NucSeq.atac,
  assay = 'peaks',
  reduction.key = 'LSI_',
  reduction.name = 'lsi'
)
NucSeq.atac <- RunUMAP(NucSeq.atac, reduction='lsi', dims=1:30)

# option 2: Process with monocle3
peak_matrix <- GetAssayData(NucSeq.atac, slot='counts')
peaks <- data.frame(as.character(rownames(peak_matrix)))
rownames(peaks) <- rownames(peak_matrix)
peaks <- as.data.frame(cbind(peaks,peaks))
colnames(peaks) <- c("GeneSymbol", "gene_short_name")
NucSeq.atac_cds <- new_cell_data_set(
  peak_matrix,
  cell_metadata=NucSeq.atac@meta.data,
  gene_metadata=genes
)
NucSeq.atac_cds <- detect_genes(NucSeq.atac_cds)
NucSeq.atac_cds <- estimate_size_factors(NucSeq.atac_cds)
NucSeq.atac_cds <- preprocess_cds(NucSeq.atac_cds, method = "LSI")
NucSeq.atac_cds <- align_cds(NucSeq.atac_cds, preprocess_method='LSI', alignment_group = "Batch")
NucSeq.atac_cds <- reduce_dimension(NucSeq.atac_cds, reduction_method = 'UMAP', preprocess_method = "LSI")
NucSeq.atac_cds <- cluster_cells(NucSeq.atac_cds)
monocle_umap <- NucSeq.atac_cds@reducedDims[["UMAP"]][rownames(NucSeq.atac_cds@reducedDims[["UMAP"]]) %in% colnames(NucSeq.atac),]
colnames(monocle_umap) <- c('UMAP_1', 'UMAP_2')
all.equal(rownames(monocle_umap), colnames(NucSeq.atac))
NucSeq.atac@reductions$umap@cell.embeddings <- monocle_umap
NucSeq.atac$monocle_clusters_umap <- clusters(NucSeq.atac_cds, reduction_method='UMAP')
monocle_aligned <- NucSeq.atac_cds@reducedDims[["Aligned"]]
rownames(monocle_aligned) <- colnames(NucSeq.atac_cds)
colnames(monocle_aligned) <- paste0('LSI_', seq(1:ncol(monocle_aligned)))
all.equal(rownames(monocle_aligned), colnames(NucSeq.atac))
NucSeq.atac@reductions$lsi <- CreateDimReducObject(
  embeddings=monocle_aligned,
  key="LSI_",
  assay="peaks"
)

################################################################################
# Step 04: Construct gene activity matrix
################################################################################


library(EnsDb.Hsapiens.v86)
gene.coords <- genes(EnsDb.Hsapiens.v86, filter = ~ gene_biotype == "protein_coding")
seqlevelsStyle(gene.coords) <- 'UCSC'
genebody.coords <- keepStandardChromosomes(gene.coords, pruning.mode = 'coarse')
genebodyandpromoter.coords <- Extend(x = gene.coords, upstream = 2000, downstream = 0)
gene.activities <- FeatureMatrix(
  fragments = fragment.path,
  features = genebodyandpromoter.coords,
  cells = colnames(NucSeq.atac),
  chunk = 10
)
gene.key <- genebodyandpromoter.coords$gene_name
names(gene.key) <- GRangesToString(grange = genebodyandpromoter.coords)
rownames(gene.activities) <- gene.key[rownames(gene.activities)]
NucSeq.atac[['RNA']] <- CreateAssayObject(counts = gene.activities)
NucSeq.atac <- NormalizeData(
  object = NucSeq.atac,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(NucSeq.atac$nCount_RNA)
)
DefaultAssay(NucSeq.atac) <- 'RNA'


################################################################################
# Step 04: Construct TF activity matrix (warning: takes a lot of time/memory!!!)
################################################################################

library(JASPAR2018)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)

pfm <- getMatrixSet(x = JASPAR2018,opts = list(species = 9606, all_versions = FALSE))
motif.matrix <- CreateMotifMatrix(features = StringToGRanges(rownames(NucSeq.atac), sep = c(":", "-")),pwm = pfm,genome = 'hg38',sep = c(":", "-"))
motif <- CreateMotifObject(data = motif.matrix,pwm = pfm)
NucSeq.atac[['peaks']] <- AddMotifObject(object = NucSeq.atac[['peaks']],motif.object = motif)
NucSeq.atac <- RunChromVAR(object = NucSeq.atac,genome = BSgenome.Hsapiens.UCSC.hg38)

################################################################################
# Step 05: integrate with snRNA-seq data
################################################################################

NucSeq.rna <- readRDS('data/snRNA.rds')
NucSeq.rna$tech <- 'rna'; NucSeq.atac$tech <- 'atac';
DefaultAssay(NucSeq.atac) <- 'RNA'


# compute anchors between RNA and ATAC
transfer.anchors <- FindTransferAnchors(
 reference=NucSeq.rna,
 query=NucSeq.atac,
 features=VariableFeatures(NucSeq.rna),
 reference.assay="RNA",
 query.assay="RNA",
 reduction="cca",
 verbose=T,
 dims=1:40
)
celltype.predictions <- TransferData(
  anchorset=transfer.anchors,
  refdata=NucSeq.rna$Cell.Type,
  weight.reduction=NucSeq.atac[["lsi"]]
)
NucSeq.atac <- AddMetaData(NucSeq.atac, celltype.predictions)
NucSeq.atac$predicted.id <- factor(NucSeq.atac$predicted.id, levels = levels(as.factor(NucSeq.rna$Cell.Type)))
Idents(NucSeq.atac) <- NucSeq.atac$monocle_clusters_umap_Cell.Type
Idents(NucSeq.rna) <- NucSeq.rna$Cell.Type
Idents(NucSeq.atac) <- paste0(Idents(NucSeq.atac), 'atac')
Idents(NucSeq.rna) <- paste0(Idents(NucSeq.rna), 'rna')
genes.use <- VariableFeatures(NucSeq.rna)
refdata <- GetAssayData(NucSeq.rna, assay = "RNA", slot = "data")[genes.use, ]
imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = NucSeq.atac[["lsi"]])
NucSeq.atac[["RNA"]] <- imputation
NucSeq.atac <- RenameCells(NucSeq.atac, add.cell.id='atac')
NucSeq.rna <- RenameCells(NucSeq.rna, add.cell.id='rna')
NucSeq.coembed <- merge(x = NucSeq.rna, y = NucSeq.atac)
cells_to_keep <- c(colnames(NucSeq.rna), colnames(NucSeq.atac)[NucSeq.atac$prediction.score.max >= 0.5])
NucSeq.coembed <- NucSeq.coembed[,cells_to_keep]
NucSeq.coembed <- ScaleData(NucSeq.coembed, features = genes.use, do.scale = FALSE)
NucSeq.coembed <- RunPCA(NucSeq.coembed, features = rownames(NucSeq.coembed), verbose = FALSE)
NucSeq.coembed <- RunUMAP(NucSeq.coembed, dims = 1:30)
expr_matrix <- GetAssayData(NucSeq.coembed, slot='data', assay='RNA')
genes <- data.frame(as.character(rownames(expr_matrix)))
rownames(genes) <- rownames(expr_matrix)
genes <- as.data.frame(cbind(genes,genes))
colnames(genes) <- c("GeneSymbol", "gene_short_name")
NucSeq_cds <- new_cell_data_set(
  expr_matrix,
  cell_metadata=NucSeq.coembed@meta.data,
  gene_metadata=genes
)
NucSeq_cds@reducedDims[['PCA']] <- NucSeq.coembed@reductions$pca@cell.embeddings
NucSeq_cds <- align_cds(NucSeq_cds, preprocess_method='PCA', alignment_group = "Batch")
NucSeq_cds <- reduce_dimension(NucSeq_cds, reduction_method = 'UMAP', preprocess_method = "Aligned")
NucSeq_cds <- cluster_cells(NucSeq_cds, reduction_method='UMAP')
