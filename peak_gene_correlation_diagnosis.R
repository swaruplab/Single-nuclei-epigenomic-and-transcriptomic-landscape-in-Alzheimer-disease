library(optparse)
library(Seurat)
library(Signac)
library(tidyverse)
library(ArchR)
library(future.apply)
library(ggpubr)
library(reshape2)
library(tictoc)
library(patchwork)
library(ggridges)
library(RColorBrewer)
library(GenomicRanges)



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

# set args manually for testing:
# cicero_data_dir <- "/dfs3/swaruplab/smorabit/analysis/AD_NucSeq_2019/atac_analysis/all_data/cicero/update/data/"
# output_dir <- "/dfs3/swaruplab/smorabit/analysis/AD_NucSeq_2019/atac_analysis/all_data/cicero/update/data/"
# celltype_index <- 1
# rna_path <- "/dfs3/swaruplab/smorabit/analysis/AD_NucSeq_2019/batch_correction/liger/update/celltype-analysis/data/NucSeq_batch_correct_seurat.rds"
# atac_path <- "/dfs3/swaruplab/smorabit/analysis/AD_NucSeq_2019/atac_analysis/all_data/ArchR3/data/NucSeq_archrPeaks_only_seurat.rds"
# group_name <- "Cell.Type"
# proj_path <- "/dfs3/swaruplab/smorabit/analysis/AD_NucSeq_2019/atac_analysis/all_data/ArchR3/all_samples"
#

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

# split objects by control & AD
cur_atac_AD <- subset(cur_atac, Diagnosis == 'AD')
cur_atac_Control <- subset(cur_atac, Diagnosis == 'Control')
cur_rna_AD <- subset(cur_rna, Diagnosis == 'AD')
cur_rna_Control <- subset(cur_rna, Diagnosis == 'Control')

# set idents to sampleIDs:
Idents(cur_atac_AD) <- cur_atac_AD$Sample.ID
Idents(cur_atac_Control) <- cur_atac_Control$Sample.ID
Idents(cur_rna_AD) <- cur_rna_AD$Sample.ID
Idents(cur_rna_Control) <- cur_rna_Control$Sample.ID


# delete stuff to save memory
rm(NucSeq.atac, NucSeq.rna); gc();
rm(cur_atac, cur_rna)

################################################################################
# Set-up before looping
################################################################################
print('Computing average expression & accessibility...')

# compute average exp/acc
average_acc_AD <- AverageExpression(cur_atac_AD)
average_acc_Control <- AverageExpression(cur_atac_Control)

average_exp_AD <- AverageExpression(cur_rna_AD)
average_exp_Control <- AverageExpression(cur_rna_Control)

# ccan cutoff from output file data/cicero-cells.o4463757.4
ccan_cutoff <- 0.01

# add peakType and nearestGene to conns_AD:
tmp1 <- proj@peakSet[na.omit(match(conns_AD$Peak1, proj@peakSet$site_name), c('peakType', 'nearestGene'))]
tmp2 <- proj@peakSet[na.omit(match(conns_AD$Peak2, proj@peakSet$site_name), c('peakType', 'nearestGene'))]

# add columns to conns_AD
conns_AD$Peak1_type <- tmp1$peakType
conns_AD$Peak1_nearestGene <- tmp1$nearestGene
conns_AD$Peak2_type <- tmp2$peakType
conns_AD$Peak2_nearestGene <- tmp2$nearestGene
conns_AD <- conns_AD[!is.na(conns_AD$Peak1_nearestGene),]

tmp1 <- proj@peakSet[na.omit(match(conns_control$Peak1, proj@peakSet$site_name), c('peakType', 'nearestGene'))]
tmp2 <- proj@peakSet[na.omit(match(conns_control$Peak2, proj@peakSet$site_name), c('peakType', 'nearestGene'))]

# add columns to conns
conns_control$Peak1_type <- tmp1$peakType
conns_control$Peak1_nearestGene <- tmp1$nearestGene
conns_control$Peak2_type <- tmp2$peakType
conns_control$Peak2_nearestGene <- tmp2$nearestGene
conns_control <- conns_control[!is.na(conns_control$Peak1_nearestGene),]

# get only links that are connected to a gene promoter region
conns_AD <- subset(conns_AD, Peak1_type == 'Promoter')
conns_control <- subset(conns_control, Peak1_type == 'Promoter')

# split conns by target gene:
conns_list_AD <- group_split(conns_AD, Peak1_nearestGene)
names(conns_list_AD) <- sort(unique(conns_AD$Peak1_nearestGene))
conns_list_control <- group_split(conns_control, Peak1_nearestGene)
names(conns_list_control) <- sort(unique(conns_control$Peak1_nearestGene))


################################################################################
# loop for AD
################################################################################
genes <- names(conns_list_AD)[names(conns_list_AD) %in% rownames(cur_rna_AD)]


df <- data.frame()
corr_list <- lapply(genes, function(cur_gene){
  print(cur_gene)
  # subset by cur_gene
  cur_conns_AD <- conns_list_AD[[cur_gene]]
  cur_conns_AD <- cur_conns_AD[!(cur_conns_AD$Peak2 %in% cur_conns_AD$Peak1),]
  cur_conns_AD <- subset(cur_conns_AD, coaccess >= ccan_cutoff)

  # skip this gene if there are no co-accessible connections
  if(nrow(cur_conns_AD) == 0){return(data.frame())}

  # get average exp and acc for this gene and peaks that are co-accessible
  # average_acc <- average_acc_AD$peaks[as.character(unique(cur_conns_AD$Peak2)),]
  average_acc <- average_acc_AD$peaks[as.character(cur_conns_AD$Peak2),]
  average_exp <- average_exp_AD$RNA[cur_gene,]

  # correlation between expression and accessibility:
  cor_mat <- apply(average_acc, 1, function(x){
    correlation <- cor.test(as.numeric(average_exp), as.numeric(x), method='pearson')
    data.frame("pval"=correlation$p.value, "pcc"=correlation$estimate)
  })

  # collapse individual correlation dfs, add info
  cor_df <- Reduce(rbind, cor_mat)
  # cor_df$Peak <- rownames(average_acc)
  # rownames(cor_df) <- cor_df$Peak
  # cor_df$Peak <- sub('[.].*', "", cor_df$Peak)
  # cor_df$target_gene <- cur_gene

  cur_conns_AD$pcc <- cor_df$pcc
  cur_conns_AD$pval <- cor_df$pval
  cur_conns_AD
})

# combine into one df and remove incomplete entries:
df <- Reduce(rbind, corr_list)
df <- df[!is.na(df$pcc),]

# # add column for peakType:
# df$PeakType <- proj@peakSet$peakType[na.omit(match(df$Peak, proj@peakSet$site_name))]
#
# # compute FDR:
# df$FDR <- p.adjust(df$pval, method='fdr')

# write df to file:
write.table(df, file=paste0(output_dir, cur_celltype, '_peak_gene_correlation_AD.csv'), sep=',', quote=FALSE, row.names=FALSE)


################################################################################
# loop for Control
################################################################################
genes <- names(conns_list_control)[names(conns_list_control) %in% rownames(cur_rna_Control)]


df <- data.frame()
corr_list <- lapply(genes, function(cur_gene){
  print(cur_gene)
  # subset by cur_gene
  cur_conns_control <- conns_list_control[[cur_gene]]
  cur_conns_control <- cur_conns_control[!(cur_conns_control$Peak2 %in% cur_conns_control$Peak1),]
  cur_conns_control <- subset(cur_conns_control, coaccess >= ccan_cutoff)

  # skip this gene if there are no co-accessible connections
  if(nrow(cur_conns_control) == 0){return(data.frame())}

  # get average exp and acc for this gene and peaks that are co-accessible
  # average_acc <- average_acc_Control$peaks[as.character(unique(cur_conns_control$Peak2)),]
  average_acc <- average_acc_Control$peaks[as.character(cur_conns_control$Peak2),]
  average_exp <- average_exp_Control$RNA[cur_gene,]

  # correlation between expression and accessibility:
  cor_mat <- apply(average_acc, 1, function(x){
    correlation <- cor.test(as.numeric(average_exp), as.numeric(x), method='pearson')
    data.frame("pval"=correlation$p.value, "pcc"=correlation$estimate)
  })

  # collapse individual correlation dfs, add info
  cor_df <- Reduce(rbind, cor_mat)
  # cor_df$Peak <- rownames(average_acc)
  # rownames(cor_df) <- cor_df$Peak
  # cor_df$Peak <- sub('[.].*', "", cor_df$Peak)
  # cor_df$target_gene <- cur_gene

  cur_conns_control$pcc <- cor_df$pcc
  cur_conns_control$pval <- cor_df$pval
  cur_conns_control
})

# combine into one df and remove incomplete entries:
df <- Reduce(rbind, corr_list)
df <- df[!is.na(df$pcc),]

# # add column for peakType:
# df$PeakType <- proj@peakSet$peakType[na.omit(match(df$Peak, proj@peakSet$site_name))]
#
# # compute FDR:
# df$FDR <- p.adjust(df$pval, method='fdr')

# write df to file:
write.table(df, file=paste0(output_dir, cur_celltype, '_peak_gene_correlation_Control.csv'), sep=',', quote=FALSE, row.names=FALSE)
