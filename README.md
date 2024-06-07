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

Installation requires the `numpy` and `matplotlib` libraries to be installed. You can install these with `pip`:
   ```bash
   pip install numpy matplotlib
```
Once required libraries are installed, you can install `kmerestimator` with the following command:
   ```bash
   python setup.py install
```
If the install was successful, typing `kmerestimator --help` should show a useful message.

## Basic Usage
The basic usage of `kmerestimator` is:
   ```bash
   kmerestimator [fastq_file] -o [output_directory]
```
- `[fastq_file]`: (Optional) Path to the input FASTQ file containing sequences.
- `-o`, `--output`: (Optional) Directory to save the histograms. Default is "histograms".

To run `kmerestimator` on a small test example (using files in this repo):
   ```bash
   kmerestimator test_2.fq -o new_histogram
```
- `test_2.fq` is the input FASTQ file.
- `new_histogram` is the directory where the histograms will be saved.

## kmerestimator options
The only required input to `kmerestimator` is a FASTQ file. Users may additionally specify the options below:
- `-o directory`, `--output directory`: Output histograms to directory. By default, the histograms are written to `histograms` directory. 

## Credits


