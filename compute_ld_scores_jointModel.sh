#!/bin/bash
#$ -N LD-Scores-joint
#$ -q bio
#$ -pe openmp 4
#$ -t 1-22

# script to compute annot file and ld scores given a bed file as input

# how to run the script:
# IN_DIR="/dfs3/swaruplab/smorabit/analysis/AD_NucSeq_2019/atac_analysis/all_data/ldsc/joint/annot/"
# PLINK_DIR="/dfs3/swaruplab/smorabit/analysis/AD_NucSeq_2019/atac_analysis/all_data/ldsc/plink_files/"
# qsub /dfs3/swaruplab/smorabit/pipelines/sn-atac-seq/bin/compute_ld_scores_join.sh \
#   $IN_DIR \
#   $PLINK_DIR
conda activate ldsc

IN_DIR=$1
PLINK_DIR=$2

PLINK_FILES=($(ls $PLINK_DIR))
PLINK_NAME=$(basename ${PLINK_FILES[1]} .bim | cut -d '.' -f 1-3)
NAME="joint"

# compute LD scores for each annotation
ANNOT_FILES=($(ls $IN_DIR*.annot.gz))
let index="$SGE_TASK_ID - 1"
annot=${ANNOT_FILES[$index]}

name=$(basename $annot .annot.gz)
echo $name
chrom=$(echo $name | cut -d '.' -f 2)
echo $chrom
~/bin/software/ldsc/ldsc.py \
  --l2 \
  --bfile $PLINK_DIR/$PLINK_NAME.$chrom \
  --ld-wind-cm 1 \
  --annot $annot \
  --out $IN_DIR/$name
