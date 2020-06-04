library(dplyr)
library(Seurat)
library(ggplot2)
library(liger)
library(Matrix)

################################################################################
# Step 01: Load 10X aggregated peaks by cells matrix into seurat
################################################################################

NucSeq.data <- Read10X(data.dir = "cellranger_aggr_dir/outs/filtered_feature_bc_matrix")
NucSeq <- CreateSeuratObject(
  counts = NucSeq.data,
  min.cells = 3,
  min.features = 200
)

# add sample level metadata (Diagnosis, Age, Sex, etc)
my.data <- read.table(file = "/dfs3/swaruplab/smorabit/analysis/AD_NucSeq_2019/data/CellSampleID.tsv", header=TRUE, stringsAsFactors=FALSE)
rownames(my.data) <- my.data$Barcode10X_Orig
my.data <- my.data[rownames(my.data) %in% colnames(NucSeq),]
NucSeq <- AddMetaData(NucSeq, metadata=my.data)
NucSeq$BRcode_sample <- paste0(do.call('rbind', strsplit(NucSeq$Barcode10X_Orig, "-"))[,1], "-", do.call("rbind", strsplit(as.character(NucSeq$Sample.ID), "-"))[,2])

################################################################################
# Step 02: Quality Control
################################################################################

NucSeq[["percent.mt"]] <- PercentageFeatureSet(object = NucSeq, pattern = "^MT-")
NucSeq <- subset(x = NucSeq, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 10)
NucSeq <- NucSeq[!grepl("^MT-", rownames(NucSeq)),]

################################################################################
# Step 03: integrative Non-negative Matrix Factorization (iNMF)
################################################################################

SeuratList <- list(
  b1 = GetAssayData(subset(NucSeq, Batch == 1), slot="counts"),
  b2 = GetAssayData(subset(NucSeq, Batch == 2), slot="counts"),
  b3 = GetAssayData(subset(NucSeq, Batch == 3), slot="counts")
)

a.NucSeq <- createLiger(SeuratList)
a.NucSeq <- normalize(a.NucSeq)
a.NucSeq <- selectGenes(a.NucSeq, var.thresh =0.4, do.plot=T)
a.NucSeq <- scaleNotCenter(a.NucSeq)
a.NucSeq <- optimizeALS(a.NucSeq, k=30)
a.NucSeq <- quantileAlignSNF(a.NucSeq, resolution = 1.0, small.clust.thresh = 20)

MergeSparseDataAll <- function(datalist, library.names = NULL) {
  col_offset <- 0
  allGenes <- unique(unlist(lapply(datalist, rownames)))
  allCells <- c()
  for (i in 1:length(datalist)) {
    curr <- datalist[[i]]
    curr_s <- summary(curr)
  curr_s[, 2] <- curr_s[, 2] + col_offset

    if (!is.null(library.names)) {
      cellnames <- paste0(library.names[i], "_", colnames(curr))
    } else {
      cellnames <- colnames(curr)
    }
    allCells <- c(allCells, cellnames)
    idx <- match(rownames(curr), allGenes)
    newgenescurr <- idx[curr_s[, 1]]
    curr_s[, 1] <- newgenescurr
    if (!exists("full_mat")) {
      full_mat <- curr_s
    } else {
      full_mat <- rbind(full_mat, curr_s)
    }
    col_offset <- length(allCells)
  }
  M <- sparseMatrix(
    i = full_mat[, 1],
    j = full_mat[, 2],
    x = full_mat[, 3],
    dims = c(
      length(allGenes),
      length(allCells)
    ),
    dimnames = list(
      allGenes,
      allCells
    )
  )
  return(M)
}

customLigerToSeurat <- function(liger_object){
  raw.data <- MergeSparseDataAll(liger_object@raw.data, names(liger_object@H))
  scale.data <- do.call(rbind, liger_object@scale.data)
  rownames(scale.data) <- colnames(raw.data)
  var.genes <- liger_object@var.genes
  var.genes <- gsub("_", replacement = "-", var.genes)
  # inmf.obj <- new(Class = "DimReduc", feature.loadings = t(liger_object@W),
  #     cell.embeddings = liger_object@H.norm, key = "iNMF_")
  inmf.obj <- CreateDimReducObject(
    loadings=t(liger_object@W),
    embeddings=liger_object@H.norm,
    key="iNMF_",
    assay="RNA"
  )
  rownames(inmf.obj@feature.loadings) <- var.genes
  rownames(inmf.obj@cell.embeddings) <- rownames(scale.data)
  new.seurat <- CreateSeuratObject(raw.data)
  new.seurat@assays$RNA@var.features <- var.genes
  new.seurat <- SetAssayData(new.seurat, slot = "scale.data",
             t(scale.data), assay = "RNA")
  new.seurat@reductions$inmf <- inmf.obj
  return(new.seurat)
}


NucSeq <- customLigerToSeurat(a.NucSeq)
my.data <- read.table(file = "data/CellSampleID.tsv", header=TRUE, stringsAsFactors=F)
for(meta in names(my.data)){
  print(meta)
  NucSeq@meta.data[[meta]] <- my.data[[meta]]
}

################################################################################
# Step 04: Primary Processing
################################################################################

NucSeq <- NormalizeData(NucSeq)
NucSeq <- ScaleData(NucSeq, features=rownames(NucSeq))
NucSeq <- RunPCA(NucSeq, dims=1:100)
NucSeq <- RunUMAP(NucSeq, reduction = "inmf", dims = 1:dim(NucSeq[["inmf"]])[2])
NucSeq <- RunTSNE(NucSeq, reduction = "inmf", dims = 1:dim(NucSeq[["inmf"]])[2])
NucSeq <- FindNeighbors(NucSeq, reduction = "inmf", dims = 1:dim(NucSeq[["inmf"]])[2], nn.eps=0.5)
NucSeq <- FindClusters(NucSeq, resolution = 0.90, n.start=10)

################################################################################
# Step 05: Re-processing Mathys et al. 2019
################################################################################

rosmap_metadata <- read.csv(file='data/ROSMAP_Clinical_2019-05_v3.csv')
tsai_metadata <- read.csv('data/TsaiMetadata.txt', sep='\t')
mathys_metadata <- read.csv('data/metaData.Tsai.merged_Final.csv', sep=',', stringsAsFactors=F)
intersect(tsai_metadata$projid, rosmap_metadata$projid)
load('data/NucSeq.Tsai.Scaled.rda')
for(meta in names(tsai_metadata)){
  NucSeq.Tsai[[meta]] <- tsai_metadata[[meta]]
}
Idents(NucSeq.Tsai) <- NucSeq.Tsai$Subcluster
for(meta in names(select(rosmap_metadata, -c(projid)))){
  NucSeq.Tsai[[meta]] <- rosmap_metadata[[meta]][match(NucSeq.Tsai$projid, rosmap_metadata$projid)]
}
NucSeq.Tsai$pathology.group <- mathys_metadata$pathology.group
NucSeq.Tsai$Diagnosis <- ifelse(NucSeq.Tsai$ceradsc %in% c(1,2), 'AD', 'Control')
NucSeq.Tsai <- NormalizeData(NucSeq.Tsai)
NucSeq.Tsai <- FindVariableFeatures(NucSeq.Tsai, nfeatures=4500)
NucSeq.Tsai <- ScaleData(NucSeq.Tsai, features=VariableFeatures(NucSeq.Tsai))
NucSeq.Tsai <- RunPCA(NucSeq.Tsai, features=VariableFeatures(NucSeq.Tsai), dims=1:100)
NucSeq.Tsai <- RunUMAP(NucSeq.Tsai, reduction = "pca", dims = 1:30)

################################################################################
# Step 06: Joint analysis of Mathys et al & UCI snRNA-seq data
################################################################################

SeuratList <- list(
  uci = GetAssayData(NucSeq, slot="counts"),
  mathys = GetAssayData(NucSeq.Tsai, slot="counts")
)
a.NucSeq <- createLiger(SeuratList)
a.NucSeq <- normalize(a.NucSeq)
a.NucSeq <- optimizeALS(a.NucSeq, k=30)
a.NucSeq <- quantileAlignSNF(a.NucSeq, resolution = 1.0, small.clust.thresh = 20)
NucSeq.joint <- customLigerToSeurat(a.NucSeq)
NucSeq.joint <- RenameCells(NucSeq.joint, new.names=do.call('rbind', strsplit(colnames(NucSeq.joint), "_"))[,2])
temp <- do.call('rbind', strsplit(colnames(NucSeq.joint), "-"))
NucSeq.joint$dataset <- temp[,1]
NucSeq.joint$barcode <- paste0(temp[,2], '-', temp[,3])
NucSeq.joint <- NormalizeData(NucSeq.joint)
NucSeq.joint <- ScaleData(NucSeq.joint, features=VariableFeatures(NucSeq.joint))
NucSeq.joint <- RunUMAP(NucSeq.joint, reduction = "inmf", dims = 1:dim(NucSeq.joint[["inmf"]])[2])

################################################################################
# Step 07: metacell aggregation
################################################################################

library(cicero)
seurat_list <- list()
k = 50
celltypes <- unique(NucSeq.joint$Cell.Type)
celltypes <- celltypes[celltypes != 'PER.END']
for(cur_celltype in celltypes){
  condition_list <- list()
  for(condition in unique(NucSeq.joint$Diagnosis)){
    print(paste(cur_celltype, condition))
    cur_seurat <- subset(NucSeq.joint, Cell.Type == cur_celltype & Diagnosis == condition)
    expr_matrix <- GetAssayData(cur_seurat, slot='data')
    genes <- data.frame(as.character(rownames(expr_matrix)))
    rownames(genes) <- rownames(expr_matrix)
    genes <- as.data.frame(cbind(genes,genes))
    colnames(genes) <- c("GeneSymbol", "gene_short_name")
    cds <- new_cell_data_set(
      expr_matrix,
      cell_metadata=cur_seurat@meta.data,
      gene_metadata=genes
    )
    cds@reducedDims[['UMAP']] <- cur_seurat@reductions$umap@cell.embeddings
    umap_coords <- reducedDims(cds)$UMAP
    metacell_cds <- make_cicero_cds(cds, reduced_coordinates=umap_coords, k=k, size_factor_normalize=FALSE)
    metacell_seurat <- CreateSeuratObject(
      counts = exprs(metacell_cds) / k
    metacell_seurat$Cell.Type <- cur_celltype
    metacell_seurat$Diagnosis <- condition
    metacell_seurat <- RenameCells(metacell_seurat, new.names=paste0(cur_celltype, '_', condition, '_', seq(1:ncol(metacell_seurat))))
    condition_list[[condition]] <- metacell_seurat
  }
  seurat_list[[cur_celltype]] <- merge(condition_list[[1]], y=condition_list[2:length(condition_list)])
}
metacell_seurat <- merge(x=seurat_list[[1]], y=seurat_list[2:length(seurat_list)])
metacell_seurat <- FindVariableFeatures(metacell_seurat, nfeatures=3000)
metacell_seurat <- ScaleData(metacell_seurat, features = VariableFeatures(metacell_seurat))
metacell_seurat <- RunPCA(metacell_seurat, features=VariableFeatures(metacell_seurat))
metacell_seurat <- RunUMAP(metacell_seurat, reduction='pca', dims=1:25)
