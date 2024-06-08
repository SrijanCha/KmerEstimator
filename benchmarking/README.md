# Benchmarking KmerEstimator against KmerGenie
For benchmarking, a version of the FASTQ file for [SRR19392958](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR19392958&display=metadata) trimmed using [sickle](https://github.com/najoshi/sickle), was used.
## Running KmerEstimator
KmerEstimator outputs the memory usage and runtime. These can be found by running 
```KmerEstimator (path_to_fastq)```
The peak memory usage can be found using the [memusg tool](https://github.com/jhclark/memusg) made by jhclark. The time can also be measured through the "time" command.
Combined, the command to find runtime, memory usage, and peak memory usage is similar to the following:
```
memusg KmerEstimator (path_to_fastq)
time KmerEstimator (path_to_fastq)
```
## Running KmerGenie
To run [KmerGenie](http://kmergenie.bx.psu.edu/) on the same dataset, the command is as follows:
```kmergenie (path_to_fastq)```
The time and memusg commands can be used to find the runtime and peak memory usage.
```
memusg kmergenie (path_to_fastq)
time kmergenie (path_to_fastq)
```
These values can then be compared for the purpose of benchmarking.
## Results
KmerEstimator was found to have an average runtime that doubled that of KmerGenie. This is expected to worsen as the dataset grows in size. However, the peak memory usage of KmerEstimator was much less than KmerGenie, being about 10x less.
