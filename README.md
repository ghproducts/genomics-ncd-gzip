
# Genomic Similarity via Normalized Compression Distance (NCD)

This project implements a Python pipeline to compute pairwise distances between genomic sequences using **Normalized Compression Distance (NCD)** — a method rooted in Kolmogorov complexity and approximated via gzip compression. It supports both genome- and taxonomy-level subsampling, multithreaded distance computation, and top-K filtering for scalability. 

---

## Features

- Enabled usage with .fasta format for input sequences
- Subsampling by **genome length** or **taxonomic hierarchy** (e.g., kingdom, phylum, species)
- Optional saving/loading of precomputed compressed lengths for fast reuse
- Parallelized NCD calculation with `ProcessPoolExecutor`
- Export of results as `.csv` file, including distances and matched labels
- Support for **top-K filtering** (`--k_min`) to reduce output memory footprint

---

## Project Structure

```
.
├── NCD.py                 # Main script for training/testing
├── dataloader.py          # FASTA loader and fragment subsampling logic
├── gziplength.py          # Gzip-based compression length utility
├── utils.py               # Taxonomic parsing, sequence validation
├── train.sh               # Example training script
├── test.sh                # Example testing script
└── data/
    ├── train/             # Directory of training FASTA files
    ├── test/              # Directory of test FASTA files
    └── metadata/          # Contains `nodes.dmp` and optionally `names.dmp`
```

---

## Requirements

- Python 3.8+
- Biopython
- NumPy
- pandas

Install dependencies with:

```bash
pip install -r requirements.txt
```

---

## Example Usage

### Training (with sampling and saving dataset)

```bash
./train.sh
```

This will:
- Load and fragment genomic FASTA files in `data/train`
- Optionally subsample by kingdom (using `--subsample_mode taxa`)
- Save the resulting dataset as `data/test_dataset.npz`

### Testing (using precomputed training dataset)

```bash
./test.sh
```

This will:
- Load test sequences from `data/test`
- Load precomputed training fragments
- Compute NCD distances
- Save a CSV result with test label, sorted distances, and matched training labels

---

## Output Format

```
test_label,NCD,train_labels
NC_000001,0.123;0.234;0.245;...,Escherichia;Klebsiella;Salmonella;...
...
```

If `--k_min` is used, only the top-k smallest NCD values are saved for each test sequence.

---

## Taxonomy Metadata Format
The taxonomy higherarchy is based on the NCBI taxanomic tree structure. This is required if any taxa-level sampling is used.

The `--taxa_data` folder must contain:
- `nodes.dmp`: Taxonomic tree structure (from NCBI)
- `names.dmp`: (optional) TaxID-to-name mapping

Used for sampling genomes by taxonomic rank (e.g., 100 fragments per species, genus, etc.)

---

## Notes

- NCD is symmetric and bounded between 0 and 1
- In this case, we utilize 126bp fragments for training and testing. We have found that 
- For large datasets, use `--k_min` to avoid memory and storage issues. The NCD matrix can become large as the number of compared sequences increases.


