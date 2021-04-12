library(tidyverse)
library(Seurat)
library(ggplot2)
library(Matrix)
library(ggrepel)
library(grid)
library(gridExtra)
library(ggpubr)
library(reshape2)
library(RColorBrewer)

# load seurat object
NucSeq <- readRDS("celltype-analysis/data/NucSeq_batch_correct_seurat.rds")

################################################################################
# Cell type marker DEGs:
################################################################################

Idents(NucSeq) <- NucSeq$Cell.Type
celltype.markers <- FindAllMarkers(object = NucSeq, only.pos = FALSE, min.pct = 0.10, test.use = "MAST",logfc.threshold = 0.25,do.print =TRUE, verbose=TRUE, assay='RNA')
save(celltype.markers,file="data/celltype_markers.rda")
write.csv(celltype.markers, file='data/celltype_markers.csv', quote=F)

################################################################################
# Control vs AD in each cell type
################################################################################

celltypes <- unique(NucSeq$Cell.Type)
combined_df <- data.frame()
for(i in 1:length(celltypes)){

  cur_celltype <- celltypes[[i]]
  print(cur_celltype)

  cur_seurat_obj <- subset(NucSeq, Cell.Type == cur_celltype)
  Idents(cur_seurat_obj) <- cur_seurat_obj$Diagnosis
  print(unique(Idents(cur_seurat_obj)))

  markers <- FindMarkers(object = cur_seurat_obj, ident.1='AD', ident.2='Control', only.pos = FALSE, min.pct = 0.10, test.use = "MAST",logfc.threshold = 0.10,do.print =TRUE, verbose=TRUE, assay='RNA')
  markers$diff <- markers$pct.1 - markers$pct.2
  markers$celltype <- cur_celltype

  combined_df <- rbind(combined_df, markers)
}
celltype.diagnosis.markers <- combined_df

save(celltype.diagnosis.markers,file="data/celltype_diagnosis_markers.rda")
write.csv(celltype.diagnosis.markers, file='data/celltype_diagnosis_markers.csv', quote=F)

################################################################################
# DEGs between different clusters within the same cell type:
################################################################################

celltypes <- unique(NucSeq$Cell.Type)
combined_df <- data.frame()
for(i in 1:length(celltypes)){

  cur_celltype <- celltypes[[i]]
  print(cur_celltype)

  cur_seurat_obj <- subset(NucSeq, Cell.Type == cur_celltype)
  Idents(cur_seurat_obj) <- cur_seurat_obj$monocle_clusters_umap_ID
  print(unique(Idents(cur_seurat_obj)))

  # compute marker genes
  markers <- FindAllMarkers(object = cur_seurat_obj, only.pos = FALSE, min.pct = 0.10, test.use = "MAST",logfc.threshold = 0.25,do.print =TRUE, verbose=TRUE, assay='RNA')
  markers$diff <- markers$pct.1 - markers$pct.2

  combined_df <- rbind(combined_df, markers)
}
cluster.markers <- combined_df

save(cluster.markers,file="data/cluster_markers.rda")
write.csv(cluster.markers, file='data/cluster_markers.csv', quote=F)

################################################################################
# control vs AD in each cluster:
################################################################################

clusters <- unique(NucSeq$monocle_clusters_umap_ID)
combined_df <- data.frame()
for(i in 1:length(clusters)){

  cur_cluster <- clusters[[i]]
  print(cur_cluster)

  cur_seurat_obj <- subset(NucSeq, monocle_clusters_umap_ID == cur_cluster)
  Idents(cur_seurat_obj) <- cur_seurat_obj$Diagnosis
  print(unique(Idents(cur_seurat_obj)))

  markers <- FindAllMarkers(object = cur_seurat_obj, only.pos = FALSE, min.pct = 0.10, test.use = "MAST",logfc.threshold = 0.10,do.print =TRUE, verbose=TRUE, assay='RNA')
  markers$diff <- markers$pct.1 - markers$pct.2
  markers$celltype <- cur_cluster

  combined_df <- rbind(combined_df, markers)
}
cluster.diagnosis.markers <- combined_df

cluster.diagnosis.markers[grepl('HSP', cluster.diagnosis.markers$gene),] %>% subset(cluster=='AD' & avg_logFC > 0)

save(cluster.diagnosis.markers,file="data/cluster_diagnosis_markers.rda")
write.csv(cluster.diagnosis.markers, file='data/cluster_diagnosis_markers.csv', quote=F)

save(celltype.markers, celltype.diagnosis.markers, cluster.markers, cluster.diagnosis.markers, file='data/all_DEGs.rda')

################################################################################
# Dot Plot
################################################################################

# reverse idents
Idents(NucSeq) <- factor(
  as.character(Idents(NucSeq)),
  levels = rev(levels(Idents(NucSeq)))
)

deg_list <- c(
  "SLC6A11",
  "TNC",
  "CHI3L1",
  "SLC15A2",
  "LINC00507",
  "IL1RAPL2",
  "RIT2",
  "CHRM2",
  "RAB3A",
  "CALB2",
  "OPRD1",
  "KLF5",
  "SV2C",
  "TLR2",
  "CX3CR1", #OR P2RY12
  "ETS1",
  "ADAMTS18",
  "CNDP1",
  "CA2",
  "OPALIN",
  "CHRM5",
  "ENPP6",
  "SLC5A11",
  "CRYAB",
  "KLHL1",
  "PMP2",
  "WDR86",
  "NOTCH3",
  "FLT1"
)

pdf('figures/DEG_dotplot.pdf', width=9, height=5, useDingbats=F)
DotPlot(NucSeq, features=deg_list, col.min=0.0, col.max=1.0, dot.min=0.20) + RotatedAxis() + scale_color_viridis(direction=-1) + xlab('') + ylab('')
dev.off()


################################################################################
# volcano plots
################################################################################

load('data/all_DEGs.rda')


## Volcano for cell type & cluster DEGs:
cur_degs <- dplyr::rename(cluster.markers, c(group=cluster))
color_scheme <- unlist(color_scheme_snRNA_clusters)
name = 'clusters'; w=18; h=12; ncols=7
cur_degs$group <- factor(as.character(cur_degs$group), levels=unique(cur_degs$group)[order(as.character(unique(cur_degs$group)))])

cur_degs <- dplyr::rename(celltype.markers, c(group=cluster))
color_scheme <- unlist(color_scheme_snRNA_celltype)
name = 'celltype'; w=18; h=3; ncols=7
cur_degs$group <- factor(as.character(cur_degs$group), levels=unique(cur_degs$group)[order(as.character(unique(cur_degs$group)))])


# get genes to annotate:
cur_degs <- Reduce(rbind, lapply(unique(cur_degs$group), function(x){
  cur <- subset(cur_degs, group == x)
  cur$anno <- ifelse(cur$avg_logFC > cur$avg_logFC[rev(order(cur$avg_logFC))][6], cur$gene, '')
  #cur$anno <- ifelse(cur$avg_logFC < cur$avg_logFC[order(cur$avg_logFC)][6], cur$gene, cur$anno)
  cur
}))


p <- ggplot(cur_degs, aes(x=avg_logFC, y=-log10(p_val_adj), color=group)) +
   geom_point(alpha=0.5) +
   scale_color_manual(values=color_scheme) +
   geom_text_repel(
    aes(label=anno), color='black') +
   # xlim(-1*max(abs(DARs_df$avg_logFC))-0.1, max(abs(DARs_df$avg_logFC))+0.1) +
   theme(panel.grid.major = element_line(colour = "lightgrey")) +

pdf(paste0('figures/volcano_', name, '_DEGs.pdf'), width=w, height=h, useDingbats=FALSE)
p + facet_wrap( ~ group, scales='free', ncol=ncol)
dev.off()


## Volcano for diagnosis DEGs:

cur_degs <- dplyr::rename(cluster.diagnosis.markers, c(group=celltype)) %>% subset(cluster=='AD')
name = 'diagnosis_clusters'; w=18; h=12; ncols=7
cur_degs$group <- factor(as.character(cur_degs$group), levels=unique(cur_degs$group)[order(as.character(unique(cur_degs$group)))])

cur_degs <- dplyr::rename(celltype.diagnosis.markers, c(group=celltype)) %>% subset(cluster=='AD')
name = 'diagnosis_celltype'; w=18; h=3; ncols=7
cur_degs$group <- factor(as.character(cur_degs$group), levels=unique(cur_degs$group)[order(as.character(unique(cur_degs$group)))])


# get genes to annotate:
cur_degs <- Reduce(rbind, lapply(unique(cur_degs$group), function(x){
  cur <- subset(cur_degs, group == x)
  cur$anno <- ifelse(cur$avg_logFC > cur$avg_logFC[rev(order(cur$avg_logFC))][6], cur$gene, '')
  cur$anno <- ifelse(cur$avg_logFC < cur$avg_logFC[order(cur$avg_logFC)][6], cur$gene, cur$anno)
  cur$color <- ifelse(cur$avg_logFC > 0, 'AD', 'Control')
  cur
}))


p <- ggplot(cur_degs, aes(x=avg_logFC, y=-log10(p_val_adj), color=color)) +
   geom_point(alpha=0.5) +
   geom_text_repel(
    aes(label=anno), color='black') +
   # xlim(-1*max(abs(DARs_df$avg_logFC))-0.1, max(abs(DARs_df$avg_logFC))+0.1) +
   theme(panel.grid.major = element_line(colour = "lightgrey"))


pdf(paste0('figures/volcano_', name, '_DEGs.pdf'), width=w, height=h, useDingbats=FALSE)
p + facet_wrap( ~ group, scales='free', ncol=ncol)
dev.off()

################################################################################
# GO term analysis
################################################################################

library(enrichR)
load('data/all_DEGs.rda')

# helper fucntion
wrapText <- function(x, len) {
    sapply(x, function(y) paste(strwrap(y, len), collapse = "\n"), USE.NAMES = FALSE)
}

dbs <- c('GO_Biological_Process_2018','GO_Cellular_Component_2018',
         'GO_Molecular_Function_2018','TRANSFAC_and_JASPAR_PWMs',
         'KEGG_2019_Human','LINCS_L1000_Chem_Pert_up',
         'LINCS_L1000_Chem_Pert_down', 'ENCODE_TF_ChIP-seq_2015',
         'TF_Perturbations_Followed_by_Expression', 'ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X',
         'TF-LOF_Expression_from_GEO', 'Transcription_Factor_PPIs')

dbs <- c('GO_Biological_Process_2018','GO_Cellular_Component_2018',
         'GO_Molecular_Function_2018')

my_theme <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(color="black"), text = element_text(size=12))

marker.genes <- cluster.diagnosis.markers
n_genes <- 200

collapsed_output <- data.frame()
for(cur in as.character(unique(marker.genes$celltype))){
  print(cur)
  # select genes
  cur_genes_AD <- marker.genes %>%
    subset(celltype == cur & cluster == 'AD' & avg_logFC >= 0) %>%
    top_n(n_genes, -p_val_adj) %>% .$gene

  cur_genes_control <- marker.genes %>%
    subset(celltype == cur & cluster == 'Control' & avg_logFC >= 0) %>%
    top_n(n_genes, -p_val_adj) %>% .$gene

  # run enrichR on different gene sets:
  cur_result_AD <- enrichr(cur_genes_AD, dbs)
  cur_result_control <- enrichr(cur_genes_control, dbs)

  # collapse results into one dataframe
  for(db in dbs){
    cur_result_AD[[db]]$cluster <- cur
    cur_result_control[[db]]$cluster <- cur
    cur_result_AD[[db]]$db <- db
    cur_result_control[[db]]$db <- db
    cur_result_AD[[db]]$Diagnosis <- 'AD'
    cur_result_control[[db]]$Diagnosis <- 'Control'
    collapsed_output <- rbind(collapsed_output, cur_result_AD[[db]])
    collapsed_output <- rbind(collapsed_output, cur_result_control[[db]])
  }
}

collapsed_output %>%
  write.csv(file='data/cluster_diagnosis_DEGs_GO_terms.csv')


colfunc <- colorRampPalette((brewer.pal(9, 'GnBu' )[3:9]))
cur_terms_AD <- read.csv('data/ProcessedGOTerms_UpAD.csv')
cur_terms_control <- read.csv('data/ProcessedGOTerms_DownAD.csv')
names(cur_terms_AD)[1] <- 'Term'; names(cur_terms_control)[1] <- 'Term';
cur_terms_AD$cluster <- factor(cur_terms_AD$cluster, levels=cur_terms_AD$cluster)
cur_terms_AD$Term <- factor(cur_terms_AD$Term, levels=unique(cur_terms_AD$Term))
cur_terms_AD$logp <- -1*log(as.numeric(cur_terms_AD$P.value))
cur_terms_AD$wrap <- wrapText(cur_terms_AD$Term, 45)
cur_terms_AD$wrap <- factor(cur_terms_AD$wrap, levels=rev(unique(cur_terms_AD$wrap)))

pdf('figures/go_terms_up_AD.pdf', width=10, height=14, useDingbats=FALSE)
g <- ggplot(data=cur_terms_AD, aes(cluster, wrap)) +
  geom_point(aes(col=logp, size=Odds.Ratio)) +
  scale_size(range=c(3,10)) +
  scale_color_gradientn(colors=colfunc(256)) +
  theme(
    axis.text.x = element_text(angle=90, hjust=1),
    axis.line = element_blank(),
    panel.grid.major = element_line(color='gray', size=0.25),
    panel.border = element_rect(color='black', fill=NA, size=1)
  )
print(g)
dev.off()

# control
cur_terms_control$cluster <- factor(cur_terms_control$cluster, levels=cur_terms_control$cluster)
cur_terms_control$Term <- factor(cur_terms_control$Term, levels=rev(unique(cur_terms_control$Term)))
cur_terms_control$logp <- -1*log(as.numeric(cur_terms_control$P.value))
cur_terms_control$wrap <- wrapText(cur_terms_control$Term, 45)
cur_terms_control$wrap <- factor(cur_terms_control$wrap, levels=rev(unique(cur_terms_control$wrap)))

pdf('figures/go_terms_down_AD.pdf', width=10, height=14, useDingbats=FALSE)
g <- ggplot(data=cur_terms_control, aes(cluster, wrap)) +
  geom_point(aes(col=logp, size=Odds.Ratio)) +
  scale_size(range=c(3,10)) +
  scale_color_gradientn(colors=colfunc(256)) +
  theme(
    axis.text.x = element_text(angle=90, hjust=1),
    axis.line = element_blank(),
    panel.grid.major = element_line(color='gray', size=0.25),
    panel.border = element_rect(color='black', fill=NA, size=1)
  )
print(g)
dev.off()
