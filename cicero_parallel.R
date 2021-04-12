
# script to run cicero in parallel for all clusters:
library(optparse)

option_list = list(
  make_option(
    c('-f', '--file'), type='character', default=NULL,
    help='Cell Dataset .rds file', metavar='character'
  ),
  make_option(
    c('-o', '--outdir'), type='character', default='./',
    help='Directory to place output files', metavar='character'
  ),
  make_option(
    c('-c', '--clusters'), type='character', default=NULL,
    help='name of the cluster subset of the CDS to process cicero'
  ),
  make_option(
    c('-n', '--num'), type='numeric', default=NULL,
    help='SGE task array number goes here, selects which cluster to process on this task'
  ),
  make_option(
    c('-g', '--genome'), type='character', default=NULL,
    help='Genome sizes file'
  )
)

# parse arguments
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

################################################################################

print(opt)

print('Loading libraries...')
library(cicero);
library(monocle3);
library(Seurat);
library(Signac);
print('Finished loading libraries.')

# load Seurat object:
print(paste('Loading Seurat object from file', opt$file))
NucSeq.atac <- readRDS(file=opt$file)
# NucSeq.atac <- readRDS('/dfs3/swaruplab/smorabit/analysis/AD_NucSeq_2019/atac_analysis/all_data/ArchR3/data/NucSeq_archrPeaks_only_seurat.rds')
print('Finished loading seurat object')

# get a list of clusters
clusters <- unique(as.character(NucSeq.atac@meta.data[[opt$clusters]]))
clusters <- clusters[order(clusters)]

print(clusters)

# load genome sizes:
human.hg38.genome <- read.csv(file=opt$genome, sep='\t', header=FALSE)
human.hg38.genome <- read.csv("/dfs3/swaruplab/smorabit/resources/hg38.chrom.sizes", sep='\t', header=FALSE)

# select cluster based on SGE task ID
cluster <- clusters[opt$num]
print(paste('Processing on group', cluster))

# subset cell dataset object by current cluster
print('Subsetting Seurat object...')
cur_seurat <- NucSeq.atac[,NucSeq.atac@meta.data[[opt$clusters]] == cluster]
print(paste('CDS dimensions:', dim(cur_seurat)))

# # are there peaks that are not accessible in this cell type?
# test <- rowSums(expr_matrix)
# sum(test == 0)

# create CDS from Seurat object:
expr_matrix <- GetAssayData(cur_seurat, slot='counts', assay='peaks')
genes <- data.frame(as.character(rownames(expr_matrix)))
rownames(genes) <- rownames(expr_matrix)
genes <- as.data.frame(cbind(genes,genes))
colnames(genes) <- c("GeneSymbol", "gene_short_name")
cur_cds <- new_cell_data_set(
  expr_matrix,
  cell_metadata=cur_seurat@meta.data,
  gene_metadata=genes
)

# # add UMAP from original analysis
# NucSeq <- readRDS('~/swaruplab/smorabit/analysis/AD_NucSeq_2019/atac_analysis/all_data/celltype-analysis/data/NucSeq_processed_activity_qc_batch_correct.rds')
# NucSeq <- NucSeq[,colnames(NucSeq)%in%colnames(NucSeq.atac)]
#
# NucSeq.atac@reductions$umap <- CreateDimReducObject(
#   embeddings = NucSeq@reductions$umap@cell.embeddings,
#   key= "UMAP_",
#   assay="peaks"
# )
#
# saveRDS(NucSeq.atac, '/dfs3/swaruplab/smorabit/analysis/AD_NucSeq_2019/atac_analysis/all_data/ArchR3/data/NucSeq_archrPeaks_only_seurat.rds')
#
#
#
# all.equal(colnames(NucSeq), colnames(NucSeq.atac))


cur_cds@reducedDims[['UMAP']] <- cur_seurat[['umap']]@cell.embeddings

# subset by Diagnosis
print('Subsetting by Diagnosis...')
cur_cds_AD <- cur_cds[, rownames(pData(cur_cds))[pData(cur_cds)$Diagnosis == 'AD']]
cur_cds_control <- cur_cds[, rownames(pData(cur_cds))[pData(cur_cds)$Diagnosis == 'Control']]
print(paste('AD CDS dimensions:', dim(cur_cds_AD)))
print(paste('Control CDS dimensions:', dim(cur_cds_control)))


# make cicero cds
print('Making cicero CDS...')
cur_cicero_cds <- make_cicero_cds(cur_cds, reduced_coordinates = reducedDims(cur_cds)$UMAP)
cur_cicero_cds_AD <- make_cicero_cds(cur_cds_AD, reduced_coordinates = reducedDims(cur_cds_AD)$UMAP)
cur_cicero_cds_control <- make_cicero_cds(cur_cds_control, reduced_coordinates = reducedDims(cur_cds_control)$UMAP)
print('Done making cicero CDS.')

# compute co-accessibility:
print('Computing co-accessibility')
conns <- run_cicero(cur_cicero_cds, human.hg38.genome)
conns_AD <- run_cicero(cur_cicero_cds_AD, human.hg38.genome)
conns_control <- run_cicero(cur_cicero_cds_control, human.hg38.genome)
print('Done computing co-accessibility')

# find co-accessibility networks
print('Compute CCANs...')
CCAN_assigns <- generate_ccans(conns)
CCAN_AD_assigns <- generate_ccans(conns_AD)
CCAN_control_assigns <- generate_ccans(conns_control)
print('Done computing CCANs.')

# save all the connections
print('Saving results')
save(conns, conns_AD, conns_control, file=paste0(opt$outdir,'/',cluster,'_cicero_connections.rda'))

# save CCANs
save(CCAN_assigns, CCAN_AD_assigns, CCAN_control_assigns, file=paste0(opt$outdir,'/',cluster,'_CCANs.rda'))

# save cicero CDS:
saveRDS(cur_cicero_cds, file=paste0(opt$outdir,'/',cluster,'_cicero_cds.rds'))
saveRDS(cur_cicero_cds_AD, file=paste0(opt$outdir,'/',cluster,'_cicero_cds_AD.rds'))
saveRDS(cur_cicero_cds_control, file=paste0(opt$outdir,'/',cluster,'_cicero_cds_control.rds'))
print('Done!')
