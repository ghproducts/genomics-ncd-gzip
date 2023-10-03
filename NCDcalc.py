from gziplength import GzipLengthCalc
import gzip
import argparse
import pickle
from tqdm import tqdm
import numpy as np


clen = lambda data: len(gzip.compress(data))
NCD = lambda c1, c2, c12: (c12 - min(c1, c2)) / max(c1, c2)

def calculate_ncd(test_data, train_data, precomputed_lengths):
    n_test = len(test_data)
    n_train = len(train_data)

    D = np.zeros((n_test, n_train))

    for i, t1 in enumerate(test_data):
        g = GzipLengthCalc(t1)
        l1 = g.length1

        for j, t2 in enumerate(train_data):
            l2 = precomputed_lengths[j]
            l12 = g.length2(t2)
            D[i, j] = NCD(l1, l2, l12)

    return D

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--dataset', required=True, help="Path of .pkl file for dataset")

    args = parser.parse_args()

    # Load test and train datasets
    ds = pickle.load(open(args.dataset, "rb"))

    # Convert strings to bytes
    train_data = [t.encode('utf8') for t in ds['train_data']]
    test_data = [t.encode('utf8') for t in ds['test_data']]

    train_lengths = []
    for j, t2 in enumerate(tqdm(train_data)):
        train_lengths.append(clen(t2))

    # Calculate NCD matrix
    ncd_matrix = calculate_ncd(test_data, train_data, train_lengths)

    # Save the NCD matrix to an output file
    output_file = f'ncd_matrix.npy'
    np.save(output_file, ncd_matrix)
    print(f"NCD matrix saved to {output_file}")

if __name__ == "__main__":
    main()
