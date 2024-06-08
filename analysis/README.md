# Analysis
KmerEstimator and KmerGenie were used on the trimmed [SRR19392958](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR19392958&display=metadata) FASTQ file. The commands to run each program are as follows:
```
KmerEstimator (path_to_fastq)
kmergenie (path_to_fastq)
```
## Results
The results may slightly vary. In this run, KmerEstimator found an optimal k-mer length of 32 while KmerGenie found the optimal k-mer length to be 67. The output kmerhistograms are in the folders labeled KmerEstimator_output and KmerGenie_output. The output for KmerEstimator is only the kmer distribution histograms. KmerGenie has a histograms_report.html that provides an overview of the kmer-length vs abundance and optimal k-mer length it finds.
