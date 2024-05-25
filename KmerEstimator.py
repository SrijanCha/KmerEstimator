import numpy as np
from typing import List
import matplotlib.pyplot as plt
from collections import Counter
import sys


def compute_kmer_dict(reads: list[str], k: int) -> dict[str, int]:
    kmer_counts = dict()
    for read in reads:
        for i in range(len(read) - k + 1):
            kmer = read[i:i+k]
            if kmer in kmer_counts:
                kmer_counts[kmer] += 1
            else:
                kmer_counts[kmer] = 1
    return kmer_counts

def compute_kmer_dist(reads: list[str], k: int) -> dict[int, int]:
    kmer_counts = compute_kmer_dict(reads, k)
    distribution = dict()

    for key in kmer_counts.keys():
        val = kmer_counts[key]

        if (val in distribution.keys()):
            distribution[val] += 1
        else:
            distribution[val] = 1
    return distribution

def get_kmer_histogram(distribution: dict[int, int], k: int, max_xval:int = None, max_yval: int = None, plot: bool = True) -> dict[int, int]:
    """ Plot a histogram of the kmer counts

    Parameters
    ----------
    distribution dictionary : dict[int, int]
        Path to the distribution of k-mers
    Returns
    -------
    Histogram of the distribution of kmer counts
    """
    if plot: 
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.bar(distribution.keys(), distribution.values())
        ax.set_xlabel("Count")
        ax.set_ylabel("Number of kmers")
        if max_xval is not None: 
            ax.set_xlim(left=0, right=max_xval)
        if max_yval is not None: 
            ax.set_ylim(bottom=0, top=max_yval)
        plt.title(f"Histogram of k-mer counts for k={k}")
        plt.show()
    return distribution
    
def read_fasta(filename: str) -> List[str]:
    """
    Read a FASTA file and return the sequences as a list of strings.

    Parameters
    ----------
    filename : str
        Path to the FASTA file.

    Returns
    -------
    List[str]
        List of sequences from the FASTA file.
    """
    sequences = []
    current_sequence = []
    with open(filename, 'r') as file:
        for line in file:
            if line.startswith('>'):
                if current_sequence:
                    sequences.append(''.join(current_sequence))
                    current_sequence = []
            else:
                current_sequence.append(line.strip())
        if current_sequence:
            sequences.append(''.join(current_sequence))
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
    maxvalue = 0
    kmer_length = 0
    check_decrease = param * 3 + 1
    
    def is_decreasing(index):
        """Check if the counts are decreasing from a given index."""
        return (kmer_counts_list[index] < kmer_counts_list[index - param] and
                kmer_counts_list[index - param] < kmer_counts_list[index - param * 2])
    
    for index in range(len(kmer_counts_list)):
        # ignore kmers that are too small
        if index < param:
            continue
        
        if maxvalue < kmer_counts_list[index]:
            maxvalue = kmer_counts_list[index]
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

def main(fasta_file: str, min_k: int, max_k: int):
    sequences = read_fasta(fasta_file)
    kmer_counts_list = []

    for k in range(min_k, max_k + 1):
        combined_kmer_counts = Counter()
        for sequence in sequences:
            kmer_counts = compute_kmer_dict(sequence, k)
            combined_kmer_counts.update(kmer_counts)
        kmer_counts_list.append(combined_kmer_counts)
        distribution = compute_kmer_dist(sequences, k)
        get_kmer_histogram(distribution, k)

    optimal_k = optimal_kmer_length(kmer_counts_list) + min_k
    print(f"Optimal k-mer length: {optimal_k}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python Counts.py <fasta_file> <min_k> <max_k>")
        sys.exit(1)
    
    fasta_file = sys.argv[1]
    min_k = int(sys.argv[2])
    max_k = int(sys.argv[3])

    main(fasta_file, min_k, max_k)

