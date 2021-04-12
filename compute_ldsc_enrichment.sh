#!/bin/bash
#$ -N LDSC-enrich
#$ -q bio
#$ -pe openmp 1
#$ -t 1-10


conda activate ldsc

SUMSTATS_DIR=$1
ANNOTATION=$2
CONTROL=$3
BASELINE=$4
OUTDIR=$5
WEIGHTS=$6
FRQ=$7

# select sumstats
SUMSTATS_FILES=($(ls $SUMSTATS_DIR*.sumstats.gz))
let index="$SGE_TASK_ID - 1"
sumstats=${SUMSTATS_FILES[$index]}

NAME=$(basename $sumstats .sumstats.gz)

~/bin/software/ldsc/ldsc.py \
  --h2 $sumstats \
  --w-ld-chr $WEIGHTS \
  --overlap-annot \
  --frqfile-chr $FRQ \
  --out $OUTDIR/$NAME \
  --ref-ld-chr $ANNOTATION,$BASELINE,$CONTROL \
  --print-coefficients
