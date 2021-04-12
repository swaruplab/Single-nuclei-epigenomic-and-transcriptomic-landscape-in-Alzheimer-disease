#!/bin/bash
#$ -N cicero-cells
#$ -q bio
#$ -pe openmp 4
#$ -t 1-9

SEURAT="/dfs3/swaruplab/smorabit/analysis/AD_NucSeq_2019/atac_analysis/all_data/ArchR3/data/NucSeq_archrPeaks_only_seurat.rds"
OUTDIR="/dfs3/swaruplab/smorabit/analysis/AD_NucSeq_2019/atac_analysis/all_data/cicero/update/data"
CLUST="monocle_clusters_umap_Cell.Type"
GENOME="/dfs3/swaruplab/smorabit/resources/hg38.chrom.sizes"

conda activate r-sc
Rscript --vanilla ~/swaruplab/smorabit/pipelines/sn-atac-seq/bin/cicero_parallel_new.R \
    -f $SEURAT \
    -o $OUTDIR \
    -c $CLUST \
    -g $GENOME \
    -n $SGE_TASK_ID
