#!/bin/bash

# example script to train NCD using the kingdom fragment method
set -e

# Input paths
TRAIN_DATA="data/train"
TAXA_DATA="data/metadata"

# Parameters
FRAGMENT_SIZE=126
RANK="superkingdom"
N=100        # total samples/superkingdom
OUTPUT_DATASET="data/test_dataset"

# Run NCD training pipeline
python NCD.py \
  --train_data "$TRAIN_DATA" \
  --taxa_data "$TAXA_DATA" \
  --subsample_mode taxa \
  --fragment_size "$FRAGMENT_SIZE" \
  --rank "$RANK" \
  --N "$N" \
  --save_dataset "$OUTPUT_DATASET"

