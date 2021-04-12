library(optparse)
library(Seurat)
library(Signac)
library(tidyverse)
library(ArchR)
library(future.apply)
library(ggpubr)
library(reshape2)

option_list = list(
  make_option(
    c('-c', '--cicerodir'), type='character', default='./',
    help='directory with cicero output files', metavar='character'
  ),
  make_option(
    c('-o', '--outdir'), type='character', default=NULL,
    help='output directory'
  ),
  make_option(
    c('-i', '--index', type='character', default=NULL,
    help='index corresponding to task ID, values should start at 1 bc R is 1-indexed!!')
  ),
  make_option(
    c('-r', '--rna', type='character', default=NULL,
    help='path to Seurat snRNA-seq object')
  ),
  make_option(
    c('-a', '--atac', type='character', default=NULL,
    help='path to Seurat snATAC-seq object')
  ),
  make_option(
    c('-g', '--group', type='character', default="Cell.Type",
    help='name of clusters/groups, ie Cell.Type')
  ),
  make_option(
    c('-p', '--proj', type='character', default=NULL,
    help='Path to ArchR project')
  )
)

# parse arguments
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# set args
cicero_data_dir <- opt$cicerodir
output_dir <- opt$outdir
celltype_index <- as.numeric(opt$index)
rna_path <- opt$rna
atac_path <- opt$atac
group_name <- opt$group
proj_path <- opt$proj

################################################################################
# Load data
################################################################################
print('Loading data...')

# loar Seurat objs
NucSeq.rna <- readRDS(rna_path)
NucSeq.atac <- readRDS(atac_path)

# load ArchR project
proj <- loadArchRProject(path = proj_path)
proj@peakSet$site_name <- rownames(NucSeq.atac)

# get list of unique groups, and get current group
celltypes <- as.character(unique(NucSeq.rna@meta.data[[group_name]]))
cur_celltype <- celltypes[celltype_index]

# load MG cicero results:
load(paste0(cicero_data_dir, cur_celltype, "_CCANs.rda"))
load(paste0(cicero_data_dir, cur_celltype, "_cicero_connections.rda"))

# subset seurat objects by microglia:
cur_rna <- NucSeq.rna[,NucSeq.rna@meta.data[[group_name]] == cur_celltype]
cur_atac <- NucSeq.atac[,NucSeq.atac@meta.data[[group_name]] == cur_celltype]

# subset atac by samples that are also in RNA:
cur_atac <- cur_atac[,cur_atac$Sample.ID %in% cur_rna$Sample.ID]
cur_atac$Sample.ID <- factor(cur_atac$Sample.ID, levels=unique(cur_atac$Sample.ID))
cur_rna$Sample.ID <- factor(cur_rna$Sample.ID, levels=unique(cur_atac$Sample.ID))

# set idents to sampleIDs:
Idents(cur_atac) <- cur_atac$Sample.ID
Idents(cur_rna) <- cur_rna$Sample.ID

# delete stuff to save memory
rm(NucSeq.atac, NucSeq.rna); gc();

################################################################################
# Set-up before looping
################################################################################
print('Computing average expression & accessibility...')

# compute average exp/acc
average_acc_all <- AverageExpression(cur_atac)
average_exp_all <- AverageExpression(cur_rna)

# ccan cutoff from output file data/cicero-cells.o4463757.4
ccan_cutoff <- 0.01

# add peakType and nearestGene to conns:
tmp1 <- proj@peakSet[na.omit(match(conns$Peak1, proj@peakSet$site_name), c('peakType', 'nearestGene'))]
tmp2 <- proj@peakSet[na.omit(match(conns$Peak2, proj@peakSet$site_name), c('peakType', 'nearestGene'))]

# add columns to conns
conns$Peak1_type <- tmp1$peakType
conns$Peak1_nearestGene <- tmp1$nearestGene
conns$Peak2_type <- tmp2$peakType
conns$Peak2_nearestGene <- tmp2$nearestGene
conns <- conns[!is.na(conns$Peak1_nearestGene),]

# get only links that are connected to a gene promoter region
conns <- subset(conns, Peak1_type == 'Promoter')

# split conns by target gene:
conns_list <- group_split(conns, Peak1_nearestGene)
names(conns_list) <- sort(unique(conns$Peak1_nearestGene))

genes <- names(conns_list)[names(conns_list) %in% rownames(cur_rna)]

################################################################################
# loop over all genes to compute correlation between accessibility & expression
################################################################################

df <- data.frame()
corr_list <- lapply(genes, function(cur_gene){

  # subset by cur_gene
  cur_conns <- conns_list[[cur_gene]]
  cur_conns <- cur_conns[!(cur_conns$Peak2 %in% cur_conns$Peak1),]
  cur_conns <- subset(cur_conns, coaccess >= ccan_cutoff)

  # skip this gene if there are no co-accessible connections
  if(nrow(cur_conns) == 0){return(data.frame())}

  # get average exp and acc for this gene and peaks that are co-accessible
  # average_acc <- average_acc_all$peaks[as.character(unique(cur_conns$Peak2)),]
  average_acc <- average_acc_all$peaks[as.character(cur_conns$Peak2),]
  average_exp <- average_exp_all$RNA[cur_gene,]

  # correlation between expression and accessibility:
  cor_mat <- apply(average_acc, 1, function(x){
    correlation <- cor.test(as.numeric(average_exp), as.numeric(x), method='pearson')
    data.frame("pval"=correlation$p.value, "pcc"=correlation$estimate)
  })

  # collapse individual correlation dfs, add info
  cor_df <- Reduce(rbind, cor_mat)
  # cor_df$Peak <- rownames(average_acc)
  # rownames(cor_df) <- cor_df$Peak
  # cor_df$target_gene <- cur_gene
  # cor_df

  cur_conns$pcc <- cor_df$pcc
  cur_conns$pval <- cor_df$pval
  cur_conns

})

# combine into one df and remove incomplete entries:
df <- Reduce(rbind, corr_list)
df <- df[!is.na(df$pcc),]

# add column for peakType:
# df$PeakType <- proj@peakSet$peakType[na.omit(match(df$Peak, proj@peakSet$site_name))]
#
# # compute FDR:
# df$FDR <- p.adjust(df$pval, method='fdr')

# write df to file:
write.table(df, file=paste0(output_dir, cur_celltype, '_peak_gene_correlation.csv'), sep=',', quote=FALSE, row.names=FALSE)
