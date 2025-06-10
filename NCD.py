from gziplength import GzipLengthCalc
from concurrent.futures import ProcessPoolExecutor
from Bio import SeqIO

import pandas as pd
import gzip
import argparse
import numpy as np
import os
import random
import time
import sys
import csv
import concurrent.futures
import glob
import pickle

from dataloader import Fasta_DataLoader
from utils import TaxaLevel

clen = lambda data: len(gzip.compress(data))
NCD = lambda c1, c2, c12: (c12 - min(c1, c2)) / max(c1, c2)

def calculate_ncd(test_data, train_data, precomputed_lengths):
    n_test = len(test_data)
    n_train = len(train_data)

    D = np.zeros((n_test, n_train))

    count = 0

    for i, t1 in enumerate(test_data):
        g = GzipLengthCalc(t1)
        l1 = g.length1
        for j, t2 in enumerate(train_data):
            l2 = precomputed_lengths[j]
            l12 = g.length2(t2)
            D[i, j] = NCD(l1, l2, l12)
            count += 1
    return D


def calculate_ncd_nmin(test_data, train_data, precomputed_lengths, k):
    """
    calculates Normalized Compression Distance (NCD) between test and train data,
    returning only the top k results for each test sequence.
    This is useful for reducing memory usage when calculating NCD on large datasets.
    """
    n_test = len(test_data)

    D = np.ones((n_test, k), dtype=float)
    top_args = np.ones((n_test, k), dtype=int)

    for i, t1 in enumerate(test_data):
        g = GzipLengthCalc(t1)
        l1 = g.length1
        for j, t2 in enumerate(train_data):
            l2 = precomputed_lengths[j]
            l12 = g.length2(t2)
            dist = NCD(l1, l2, l12)
            if dist < max(D[i]):
                idx = np.argmax(D[i]) 
                D[i, idx] = dist
                top_args[i, idx] = j
    D = np.hstack([top_args, D])
    return D


def train_NCD(args):

 # validate arguments
    if args.subsample_mode is not None and args.subsample_mode not in {"genome", "taxa"}:
        raise ValueError("Invalid --subsample_mode. Must be 'genome' or 'taxa' if specified.")
    if args.subsample_mode == "taxa" and args.taxa_data is None:
        raise ValueError("When --subsample_mode is 'taxa', --taxa_data must be provided.")
    else:
        # Load taxa data if subsampling by taxa
        if args.taxa_data:
            args.taxa_data = TaxaLevel(args.taxa_data)
            args.taxa_data.parse_files()

    # Load training data and calculate lengths
    train_data = Fasta_DataLoader(args.train_data)
    train_data.load_data()
    if args.subsample_mode:
        train_data.subsample(mode=args.subsample_mode, taxa_data=args.taxa_data, fragment_size=args.fragment_size, S=args.S, rank=args.rank, N=args.N)

    # encode and get lengths
    train_data.data['clen'] = [clen(seq.encode('utf8')) for seq in train_data.data['seq']]
    
    if args.save_dataset:
        np.savez(f"{args.save_dataset}.npz", data=train_data.data)


    return train_data


def main():
    parser = argparse.ArgumentParser(description="Calculate Normalized Compression Distance (NCD) between test and train data.")
    # general parameters
    parser.add_argument("--test_data", type=str, default=None, help="Path to the test data directory. If None, testing will be skipped. This is useful for isolating training")
    parser.add_argument("--train_data", type=str, required=True, help="Path to the train data. If raw fasta data, this should be a directory. Otherwise this should be the path to a precomputed NCD matrix (npz file).")
    parser.add_argument("--output", type=str, default="NCD_output", help="Path to save the output NCD matrix.")
    parser.add_argument("--precomputed_lengths", action="store_true", help="Flag to indicate if training data lengths are already precomputed. If True, train_data should be a precomputed NCD matrix (npz file). If False, the script will compute lengths from the fasta files.")
    # training parameters
    parser.add_argument("--subsample_mode", type=str, default=None, help="Subsampling mode: 'genome' or 'taxa'.")
    parser.add_argument("--fragment_size", type=int, default=126, help="Size of fragments for subsampling.")
    parser.add_argument("--S", type=int, default=1000, help="Number of samples per taxa for subsampling.")
    parser.add_argument("--rank", type=str, choices=["superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species"], default="species", help="Taxonomic rank for subsampling.")
    parser.add_argument("--N", type=int, default=10, help="Number of samples per class in taxa for subsampling.")
    parser.add_argument("--taxa_data", type=str, default=None, help="Path to folder containing nodes.dmp and names.dmp. This is required if taxa-level sampling is used.") 
    parser.add_argument("--save_dataset",type=str, default=None , help="Path to save the dataset after training. If None, the dataset is not saved.")
    # testing parameters
    parser.add_argument("--num_threads", type=int, default=4, help="Number of threads to use for parallel processing.")
    parser.add_argument("--k_min", type=int, default=None, help="Number of top NCD results to save for each test sequence. If None, all results are saved.")
    args = parser.parse_args()

    # train NCD if precomputed lengths are not provided
    if not args.precomputed_lengths:
        print("[INFO] Precomputed lengths not provided. Training NCD...")
        train_data = train_NCD(args)
    else:  
        print("[INFO] Precomputed lengths provided. Loading training data...")
        train_data = Fasta_DataLoader(None)
        train_data.data = np.load(f"{args.train_data}", allow_pickle=False)["data"]

    if args.test_data is None:
        print("[INFO] No --test_data provided. Skipping NCD calculation.")
        return


    # load test data
    test_data = Fasta_DataLoader(args.test_data)
    test_data.load_data()

    # in the case of k_min returns (for reduced memory usage when calculating NCD on large datasets)
    if args.k_min is not None:
        if args.k_min <= 0:
            raise ValueError("k_min must be a positive integer.")
        if args.k_min > len(train_data.data['seq']):
            raise ValueError("k_min cannot be greater than the number of training sequences.")
    # calculate NCD with multithreading
    encoded_train_data = [seq.encode('utf8') for seq in train_data.data['seq']]
    with ProcessPoolExecutor(max_workers=args.num_threads) as executor:
        print(f"[INFO] Calculating NCD with {args.num_threads} threads...")
        futures = []
        for i in range(0, len(test_data.data['seq']), args.num_threads):
            batch = test_data.data['seq'][i:i + args.num_threads]
            batch = [seq.encode('utf8') for seq in batch]  # encode test sequences
            if args.k_min is not None:
                futures.append(executor.submit(calculate_ncd_nmin, batch, encoded_train_data, train_data.data['clen'], args.k_min))
            else:
                futures.append(executor.submit(calculate_ncd, batch, encoded_train_data, train_data.data['clen']))
        results = [future.result() for future in concurrent.futures.as_completed(futures)]

    print("[INFO] NCD calculation completed.")
    # save results
    result_matrix = np.vstack(results)
    test_labels = test_data.data['label']
    train_labels = train_data.data['label']
    rows = []

    for i, row in enumerate(result_matrix):
        if args.k_min is not None:
            k = args.k_min
            indices = row[:k].astype(int)
            distances = row[k:]
            selected_train_labels = train_labels[indices]
        else:
            # Sort the full NCD distances and align train labels
            sort_idx = np.argsort(row)
            distances = row[sort_idx]
            selected_train_labels = train_labels[sort_idx]

        rows.append({
            "test_label": test_labels[i],
            "NCD": ';'.join(f"{d:.6f}" for d in distances),
            "train_labels": ';'.join(selected_train_labels)
        })

    df = pd.DataFrame(rows)
    df.to_csv(args.output, index=False)
    print(f"[INFO] Results saved to {args.output}")



if __name__ == "__main__":
    main()

