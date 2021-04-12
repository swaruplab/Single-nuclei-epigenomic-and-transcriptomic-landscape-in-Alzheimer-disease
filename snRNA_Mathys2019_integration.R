library(tidyverse)
library(Seurat)
library(cowplot)
theme_use(theme_cowplot())

NucSeq <- readRDS(file="/dfs7/dfs3/swaruplab/smorabit/analysis/cross-disorder/data/seurat_objects/swarup_AD_2019_unprocessed.rds")

rosmap_metadata <- read.csv(file='/dfs3/swaruplab/smorabit/analysis/AD_NucSeq_2019/batch_correction/liger/update/joint_analysis/data/ROSMAP_Clinical_2019-05_v3.csv')
tsai_metadata <- read.csv('/dfs3/swaruplab/smorabit/analysis/AD_NucSeq_2019/batch_correction/liger/update/joint_analysis/data/TsaiMetadata.txt', sep='\t')
intersect(tsai_metadata$projid, rosmap_metadata$projid)

# load tsai seurat object:
load('~/swaruplab/vswarup/AD2019/IntegrateTsai/NucSeq.Tsai.Scaled.rda')

# add metadata to Tsai object
for(meta in names(tsai_metadata)){
  NucSeq.Tsai[[meta]] <- tsai_metadata[[meta]]
}
Idents(NucSeq.Tsai) <- NucSeq.Tsai$Subcluster

for(meta in names(select(rosmap_metadata, -c(projid)))){
  NucSeq.Tsai[[meta]] <- rosmap_metadata[[meta]][match(NucSeq.Tsai$projid, rosmap_metadata$projid)]
}

# make Diagnosis column based on CERAD score:
NucSeq.Tsai$Diagnosis <- ifelse(NucSeq.Tsai$ceradsc %in% c(1,2), 'AD', 'Control')

# process their data!

NucSeq.Tsai <- NormalizeData(NucSeq.Tsai)
NucSeq.Tsai <- FindVariableFeatures(NucSeq.Tsai, nfeatures=4500)
NucSeq.Tsai <- ScaleData(NucSeq.Tsai, features=VariableFeatures(NucSeq.Tsai))
NucSeq.Tsai <- RunPCA(NucSeq.Tsai, features=VariableFeatures(NucSeq.Tsai), dims=1:100)
NucSeq.Tsai <- RunUMAP(NucSeq.Tsai, reduction = "pca", dims = 1:30)

pdf('figures/tsai_umap.pdf', width=7, height=6)
DimPlot(NucSeq.Tsai, group.by='broad.cell.type')
DimPlot(NucSeq.Tsai, group.by='Diagnosis')
DimPlot(NucSeq.Tsai, group.by='individualID')
dev.off()

library(liger)

SeuratList <- list(
  Swarup = GetAssayData(NucSeq, slot="counts"),
  Tsai = GetAssayData(NucSeq.Tsai, slot="counts")
)
colnames(SeuratList$Swarup) <- paste0('Swarup-', colnames(SeuratList$Swarup))
colnames(SeuratList$Tsai) <- paste0('Tsai-', colnames(SeuratList$Tsai))

a.NucSeq <- createLiger(SeuratList)

# normalize, select variable genes, scale
a.NucSeq <- normalize(a.NucSeq)
pdf("../joint_analysis/figures/liger_variable_genes.pdf", width=8, height=8)
a.NucSeq <- selectGenes(a.NucSeq, var.thresh =c(0.4, 0.1), do.plot=T)
dev.off()
a.NucSeq <- scaleNotCenter(a.NucSeq)

a.NucSeq <- optimizeALS(a.NucSeq, k=30)
a.NucSeq <- quantileAlignSNF(a.NucSeq, resolution = 1.0, small.clust.thresh = 20)
save(a.NucSeq, file="../joint_analysis/data/NucSeq_joint_liger.rds")

# plot TSNE
a.NucSeq <- runTSNE(a.NucSeq, method="fftRtsne", fitsne.path="/data/users/smorabit/bin/software/FIt-SNE")
pdf("../joint_analysis/figures/liger_tsne.pdf", width=10, height=8)
plotByDatasetAndCluster(a.NucSeq, return.plots=T)
dev.off()

a.NucSeq@tsne.coords <- NucSeq@reductions$tsne@cell.embeddings
# plot gene loadings for each factor
pdf("../joint_analysis/figures/liger_geneLoadings.pdf", width=10, height=10)
plotGeneLoadings(a.NucSeq, num.genes=50)
dev.off()

NucSeq.joint <- customLigerToSeurat(a.NucSeq)
NucSeq.joint <- RenameCells(NucSeq.joint, new.names=do.call('rbind', strsplit(colnames(NucSeq.joint), "_"))[,2])

# set up metadata ##############################################################
temp <- do.call('rbind', strsplit(colnames(NucSeq.joint), "-"))
NucSeq.joint$dataset <- temp[,1]
NucSeq.joint$barcode <- paste0(temp[,2], '-', temp[,3])

# are cells in the same order as they are in their dataset of origin?
all.equal(subset(NucSeq.joint@meta.data, dataset=='Swarup') %>% .$barcode, colnames(NucSeq))
all.equal(subset(NucSeq.joint@meta.data, dataset=='Tsai') %>% .$barcode, colnames(NucSeq.Tsai))

# take Diagnosis and cell type from dataset of origin
NucSeq.joint$Diagnosis <- c(NucSeq$Diagnosis, NucSeq.Tsai$Diagnosis)
NucSeq.joint$celltype <- c(NucSeq$Cell.Type, as.character(NucSeq.Tsai$broad.cell.type))
NucSeq.joint$SampleID <- c(NucSeq$SampleID, NucSeq.Tsai$individualID)

NucSeq.joint <- NormalizeData(NucSeq.joint)
NucSeq.joint <- ScaleData(NucSeq.joint, features=VariableFeatures(NucSeq.joint))
NucSeq.joint <- RunUMAP(NucSeq.joint, reduction = "inmf", dims = 1:dim(NucSeq.joint[["inmf"]])[2])

NucSeq.joint$UMAP_1 <- NucSeq.joint@reductions$umap@cell.embeddings[,1]
NucSeq.joint$UMAP_2 <- NucSeq.joint@reductions$umap@cell.embeddings[,2]


saveRDS(NucSeq.joint, file='../joint_analysis/data/NucSeq_joint_inmf_integration.rds')

pdf('figures/joint_umap_inmf.pdf', width=12, height=10)
DimPlot(NucSeq.joint, group.by='dataset')
DimPlot(NucSeq.joint, group.by='Diagnosis')
DimPlot(NucSeq.joint, group.by='SampleID')
DimPlot(NucSeq.joint, group.by='celltype')
dev.off()

pdf('figures/joint_umap_inmf_split.pdf', width=12, height=5)
DimPlot(NucSeq.joint, group.by='Diagnosis', split.by='dataset')
DimPlot(NucSeq.joint, group.by='SampleID', split.by='dataset')
DimPlot(NucSeq.joint, group.by='celltype', split.by='dataset')
dev.off()

png('figures/pngs/joint_umap_inmf_dataset.png', width=6, height=5, res=500, units='in')
DimPlot(NucSeq.joint, group.by='dataset')
dev.off()

png('figures/pngs/joint_umap_inmf_celltype.png', width=6, height=5, res=500, units='in')
DimPlot(NucSeq.joint, group.by='celltype')
dev.off()

png('figures/pngs/joint_umap_inmf_split.png', width=12, height=5, res=500, units='in')
DimPlot(NucSeq.joint, group.by='celltype', split.by='dataset')
dev.off()
