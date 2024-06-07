# KmerEstimator

## Overview

The `KmerEstimator.py` script is designed to estimate the optimal k-mer length for a given sequence dataset. The script reads sequences from a FASTQ file, generates k-mers of varying lengths, and computes the optimal k-mer length based on the distribution of k-mer counts.

## Features

- **K-mer Generation**: Generate k-mers of varying lengths from input sequences using the compute_kmer_dict function.
- **K-mer Count Distribution**: Calculate the k-mer count distribution using the compute_kmer_dist function.
- **Optimal K-mer Length Estimation**: Determine the optimal k-mer length based on the k-mer count distributions using the optimal_kmer_length function.
- **Histogram Plotting**: Visualize k-mer count distributions as histograms using the get_kmer_histogram and kmer_num_histogram functions.
- **Memory Usage Monitoring**: Monitor memory usage during k-mer analysis using the memory_usage function.
- **Run K-mer Analysis**: Run a comprehensive k-mer analysis over a range of k values and refine the optimal k-mer length using the run_kmer_analysis function.

## Requirements

- Python 3.6 or higher
- NumPy
- Matplotlib

## Installation

1. Clone the repository or download the `KmerEstimator.py` file.
2. Ensure you have the required libraries installed:
   ```bash
   pip install biopython numpy matplotlib

## Usage
   ```bash
   python KmerEstimator.py <input_file> -o <output_directory>
```
- `<input_file>`: Path to the input file containing sequences.
- `-o`, `--output`: (Optional) Directory to save the histograms. Default is `"histograms"`.

## Complete Usage
   ```bash
   python KmerEstimator.py [fastq_file] -o [output_directory]
```
- `[fastq_file]`: (Optional) Path to the input FASTQ file containing sequences.
- `-o`, `--output`: (Optional) Directory to save the histograms. Default is "histograms".

## Example
   ```bash
   python KmerEstimator.py seq.fq -o histogram_dir
```
In this example:

- `test_2.fq` is the input FASTQ file.
- `new_histogram` is the directory where the histograms will be saved.

## Output 

The script generates an output containing:

- Optimal k-mer length.
- Relevant statistics and analysis for each k-mer length considered.
- Histogram images saved in the specified output directory.

## Credits


