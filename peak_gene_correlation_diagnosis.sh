#!/bin/bash
#$ -N cCREs_diagnosis
#$ -q bio
#$ -pe openmp 4
#$ -t 1-7

# call R script:
conda activate r-sc
Rscript --vanilla ~/swaruplab/smorabit/pipelines/sn-atac-seq/bin/peak_gene_correlation_diagnosis.R \
  -c "/dfs3/swaruplab/smorabit/analysis/AD_NucSeq_2019/atac_analysis/all_data/cicero/update/data/" \
  -o "/dfs3/swaruplab/smorabit/analysis/AD_NucSeq_2019/atac_analysis/all_data/cicero/update/data/" \
  -i $SGE_TASK_ID \
  -r "/dfs3/swaruplab/smorabit/analysis/AD_NucSeq_2019/batch_correction/liger/update/celltype-analysis/data/NucSeq_batch_correct_seurat.rds" \
  -a "/dfs3/swaruplab/smorabit/analysis/AD_NucSeq_2019/atac_analysis/all_data/ArchR3/data/NucSeq_archrPeaks_only_seurat.rds" \
  -g "Cell.Type" \
  -p "/dfs3/swaruplab/smorabit/analysis/AD_NucSeq_2019/atac_analysis/all_data/ArchR3/all_samples"
