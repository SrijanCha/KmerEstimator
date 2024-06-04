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
    # kmer_counts = compute_kmer_dict(reads, k)
    # distribution = dict()

    # for key in kmer_counts.keys():
    #     val = kmer_counts[key]

    #     if (val in distribution.keys()):
    #         distribution[val] += 1
    #     else:
    #         distribution[val] = 1
    return distribution

def get_kmer_histogram(distribution: Counter, k: int, max_xval:int = None, max_yval: int = None) -> None:
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
    plt.title(f"Histogram of k-mer counts for k={k}")
    plt.show()
    # if plot: 
    #     fig = plt.figure()
    #     ax = fig.add_subplot(111)
    #     ax.bar(distribution.keys(), distribution.values())
    #     ax.set_xlabel("Count")
    #     ax.set_ylabel("Number of kmers")
    #     if max_xval is not None: 
    #         ax.set_xlim(left=0, right=max_xval)
    #     if max_yval is not None: 
    #         ax.set_ylim(bottom=0, top=max_yval)
    #     plt.title(f"Histogram of k-mer counts for k={k}")
    #     plt.show()
    # return distribution
    
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
    current_sequence = []
    with open(filename, 'r') as file:
        line_count = 0
        for line in file:
            line_count += 1
            if line_count % 4 == 2:
                sequences.append(line.strip())
        # for line in file:
        #     if line.startswith('>'):
        #         if current_sequence:
        #             sequences.append(''.join(current_sequence))
        #             current_sequence = []
        #     else:
        #         current_sequence.append(line.strip())
        # if current_sequence:
        #     sequences.append(''.join(current_sequence))
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


def random_reads_generator(n, l):
    alphabet = ["A", "C", "T", "G"]
    reads = [''.join(random.choices(alphabet, k=l)) for _ in range(n)]
    return reads

def run_kmer_analysis(sequences: List[str], min_k: int, max_k: int) -> Tuple[int, List[Counter]]:
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
    for k in range(min_k, max_k + 1):
        combined_kmer_counts = Counter()
        for sequence in sequences:
            kmer_counts = compute_kmer_dict([sequence], k)
            combined_kmer_counts.update(kmer_counts)
        kmer_counts_list.append(combined_kmer_counts)
        distribution = compute_kmer_dist(combined_kmer_counts)
        get_kmer_histogram(distribution, k)

        # Debug print to check kmer counts for each k
        # print(f"k-mer counts for k={k}: {combined_kmer_counts}")

    optimal_k = optimal_kmer_length(kmer_counts_list) + min_k
    return optimal_k, kmer_counts_list

def main(fastq_file: str, min_k: int, max_k: int):
    sequences = read_fastq(fastq_file)
    print(f"read {len(sequences)} sequences from {fastq_file}")
    if not sequences:
        print("No sequences found. Please check the FASTQ file.")
        return
    
    # Print the first few sequences
    # print("First few sequences:")
    # for seq in sequences[:5]:
    #     print(seq)

    optimal_k, kmer_counts_list = run_kmer_analysis(sequences, min_k, max_k)
    print(f"Optimal k-mer length: {optimal_k}")
    # for k, kmer_counts in zip(range(min_k, max_k + 1), kmer_counts_list):
    #     print(f"\nk-mer counts for k={k}:")
    #     for kmer, count in kmer_counts.items():
    #         print(f"{kmer}: {count}")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python KmerEstimator.py <fastq_file> <min_k> <max_k>")
        sys.exit(1)
    
    fastq_file = sys.argv[1]
    min_k = int(sys.argv[2])
    max_k = int(sys.argv[3])

    main(fastq_file, min_k, max_k)
