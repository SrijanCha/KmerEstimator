# KmerEstimator

## Overview

The `KmerEstimator.py` script is designed to estimate the optimal k-mer length for a given sequence dataset. The script reads sequences from a FASTA file, generates k-mers of varying lengths, and computes the optimal k-mer length based on the distribution of k-mer counts.

## Features

- **K-mer Generation**: Generate k-mers of varying lengths from input sequences.
- **Optimal K-mer Length Estimation**: Calculate the optimal k-mer length based on k-mer count distributions.
- **Histogram Plotting**: Visualize k-mer count distributions as histograms.
- **Random Reads Generation**: Generate random reads for testing purposes.

To be done:
- **Error Handling**: Account for sequencing errors and other noise in the data.
- **Computational Efficiency**: Balance between capturing meaningful patterns and computational resource usage.

## Requirements

- Python 3.6 or higher
- Biopython
- NumPy
- Matplotlib

## Installation

1. Clone the repository or download the `KmerEstimator.py` file.
2. Ensure you have the required libraries installed:
   ```bash
   pip install biopython numpy matplotlib

## Usage
   ```bash
   python KmerEstimator.py <input_file> <min_k> <max_k>
```
- `<input_file>`: Path to the input file containing sequences.
- `<min_k>`: Minimum k-mer length to consider (should set default).
- `<max_k>`: Maximum k-mer length to consider (should set default).

## Output 

The script generates an output file containing the optimal k-mer length along with relevant statistics and analysis for each k-mer length considered.
