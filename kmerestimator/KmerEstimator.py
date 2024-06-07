import os
import numpy as np
from typing import List, Tuple
import matplotlib.pyplot as plt
from collections import Counter
import sys
import random


def compute_kmer_dict(reads: list[str], k: int) -> Counter:
    """
    Compute k-mer counts from a list of reads.
    
    Parameters
    ----------
    reads : List[str]
        List of sequence reads.
    k : int 
        Length of k-mers.
    
    Returns
    -------
    Counter
        Dictionary of k-mer counts.
    """
    kmer_counts = Counter()
    for read in reads:
        for i in range(len(read) - k + 1):
            kmer = read[i:i+k]
            if kmer in kmer_counts:
                kmer_counts[kmer] += 1
            else:
                kmer_counts[kmer] = 1
    return kmer_counts

def compute_kmer_dist(kmer_counts: Counter) -> Counter:
    """
    Compute k-mer count distribution from k-mer counts.
    
    Parameters
    ----------
    kmer_counts : Counter
        Dictionary of k-mer counts.
        
    Returns
    -------
    Counter
        Dictionary of k-mer count distribution.
    """
    distribution = Counter(kmer_counts.values())

    return distribution

def get_kmer_histogram(distribution: Counter, k: int, max_xval:int = None, max_yval: int = None, save_dir: str = None) -> None:
    """ Plot a histogram of the kmer counts

    Parameters
    ----------
    distribution : Counter
        Dictionary of k-mer count distribution
    k : int
        Length of k-mers.
    max_xval : int, optional
        Maximum x-axis value for the plot.
    max-yval : int, optional
        Maximum y-axis value for the plot.

    """
    fig, ax = plt.subplots()
    ax.bar(distribution.keys(), distribution.values())
    ax.set_xlabel("Count")
    ax.set_ylabel("Number of kmers")
    if max_xval is not None:
        ax.set_xlim(left=0, right=max_xval)
    if max_yval is not None:
        ax.set_ylim(bottom=0, top=max_yval)
    ax.set_yscale('log')
    
    plt.title(f"Histogram of k-mer counts for k={k}")
    if save_dir:
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
        filename = os.path.join(save_dir, f"k_{k}_distribution.png")
        plt.savefig(filename)
    else:
        plt.show()

    
def read_fastq(filename: str) -> List[str]:
    """
    Read a FASTQ file and return the sequences as a list of strings.

    Parameters
    ----------
    filename : str
        Path to the FASTQ file.

    Returns
    -------
    List[str]
        List of sequences from the FASTQ file.
    """
    sequences = []
    with open(filename, 'r') as file:
        line_count = 0
        for line in file:
            line_count += 1
            if line_count % 4 == 2:
                sequences.append(line.strip())
    return sequences


def optimal_kmer_length(kmer_counts_list: List[Counter]) -> int:
    """
    Return optimal kmer length based on list of kmer counts.

    Parameters
    ----------
    kmer_counts_list: List[Counter]
        List of counts for each kmer length

    Returns
    -------
    int
        Optimal kmer length
    """
    # int set to 1/20 of kmer_counts_list (arbitrary value)
    param = round(len(kmer_counts_list) / 20)
    max_unique_kmers = 0
    kmer_length = 0
    check_decrease = param * 3 + 1
    
    def is_decreasing(index):
        """Check if the counts are decreasing from a given index."""
        return (len(kmer_counts_list[index]) < len(kmer_counts_list[index - param]) and
                len(kmer_counts_list[index - param]) < len(kmer_counts_list[index - param * 2]))
    
    for index in range(len(kmer_counts_list)):
        # ignore kmers that are too small
        if index < param:
            continue
        
        unique_kmers = len(kmer_counts_list[index])

        if max_unique_kmers < unique_kmers:
            max_unique_kmers = unique_kmers
            kmer_length = index
        
        # end search early if kmer count is consistently decreasing
        if index >= check_decrease and is_decreasing(index):
            break
        
        if index == check_decrease:
            check_decrease += param * 3
    return kmer_length


def kmer_num_histogram(kmer_counts_list: List[Counter], max_xval: int, save_dir: str = None):
    """
    Create histogram of kmer size vs number of kmers. This should be concave and show a clear global maximum for best k to be accurate.
    Parameters
    ----------
    kmer_counts_list : List[Counter]
        List of counts for each kmer-length
    max_xval: 
        Max k-mer size given by user
    """
    # int set to 1/20 of kmer_counts_list (arbitrary value)
    
    M_for_k = dict()
    
    for index in range(len(kmer_counts_list)):
        M_for_k[index] = len(kmer_counts_list[index])
    
    fig, ax = plt.subplots()
    
    ax.bar(M_for_k.keys(), M_for_k.values())
    ax.set_xlabel("k-mer Size")
    ax.set_ylabel("Number of Unique k-mers")
    ax.set_xlim(left=0, right=max_xval)
    
    plt.title(f"Histogram of Unique Kmers for values of k")
    if save_dir:
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
        filename = os.path.join(save_dir, f"all_k_distribution.png")
        plt.savefig(filename)
    else:
        plt.show()
    
    

def run_kmer_analysis(sequences: List[str]) -> Tuple[int, List[Counter]]:
    """
    Run k-mer analysis over a range of k values and determine the optimal k-mer length.
    
    Parameters
    ----------
    sequences : List[str]
        List of sequence reads.
    min_k : int
        Minimum k-mer size to consider.
    max_k : int
        Maximum k-mer size to consider.
    
    Returns
    -------
    Tuple[int, List[Counter]]
        Optimal k-mer length and list of k-mer counts.
    """
    kmer_counts_list = []
    save_dir = "histograms"
    k_values = list(range(11, 132, 10))
    min_seq_len = len(sequences[0])

    for k in k_values:
        if k > min_seq_len:
            continue
        print(k)
        combined_kmer_counts = Counter()
        for sequence in sequences:
            kmer_counts = compute_kmer_dict([sequence], k)
            combined_kmer_counts.update(kmer_counts)
        kmer_counts_list.append(combined_kmer_counts)
        distribution = compute_kmer_dist(combined_kmer_counts)
        get_kmer_histogram(distribution, k, save_dir=save_dir)

    optimal_k_initial = optimal_kmer_length(kmer_counts_list) + 11

    # Refine the search in the range [optimal_k_initial - 10, optimal_k_initial + 10]
    refined_k_values = list(range(max(11, optimal_k_initial - 10), min(131, optimal_k_initial + 10) + 1))
    kmer_counts_list_refined = []

    for k in refined_k_values:
        if k > min_seq_len:
            continue
        print(k)
        combined_kmer_counts = Counter()
        for sequence in sequences:
            kmer_counts = compute_kmer_dict([sequence], k)
            combined_kmer_counts.update(kmer_counts)
        kmer_counts_list_refined.append(combined_kmer_counts)
        distribution = compute_kmer_dist(combined_kmer_counts)
        get_kmer_histogram(distribution, k, save_dir=save_dir)

    optimal_k_final = optimal_kmer_length(kmer_counts_list_refined) + refined_k_values[0]

    kmer_num_histogram(kmer_counts_list_refined, max_xval=refined_k_values[-1], save_dir=save_dir)

    return optimal_k_final, kmer_counts_list_refined


def main(fastq_file: str):
    sequences = read_fastq(fastq_file)
    print(f"read {len(sequences)} sequences from {fastq_file}")
    if not sequences:
        print("No sequences found. Please check the FASTQ file.")
        return

    optimal_k, kmer_counts_list = run_kmer_analysis(sequences)
    print(f"Optimal k-mer length: {optimal_k}")
    # for k, kmer_counts in zip(range(min_k, max_k + 1), kmer_counts_list):
    #     print(f"\nk-mer counts for k={k}:")
    #     for kmer, count in kmer_counts.items():
    #         print(f"{kmer}: {count}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python KmerEstimator.py <fastq_file>")
        sys.exit(1)
    
    fastq_file = sys.argv[1]

    main(fastq_file)
