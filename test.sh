#!/bin/bash

#example script to run testing on trained dataset
set -e

# Input paths
TEST_DATA="data/test"
PRECOMPUTED_TRAIN="data/test_dataset.npz"

# Parameters
NUM_THREADS=8
OUTPUT="ncd_kingdom_output.csv"

# Run NCD test pipeline
python NCD.py \
  --train_data "$PRECOMPUTED_TRAIN" \
  --test_data "$TEST_DATA" \
  --precomputed_lengths \
  --num_threads "$NUM_THREADS" \
  --output "$OUTPUT"

