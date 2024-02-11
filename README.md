# RNAseqQualCheck
This repository contains script to assess the quality of RNAseq data sets.


## rRNA check ##
This script analyzes all reads in a FASTQ file. If multiple k-mers match an rRNA sequence, the read is classified as rRNA derived. This allows the comparison of different FASTQ files regarding the contamination of rRNA in different samples.

```
Usage
python3 rRNA_check.py --fastq <FILE> --rRNA <FILE> --out <DIR>

Mandatory:
--fastq    STR   FASTQ input file
--rRNA     STR   rRNA FASTA input file
--out      STR   Output folder

Optional:
--kmer     STR    Kmer size [21]
--cutoff   STR    Cutoff for read classification [3]
```

`--fastq` specifies a FASTQ input file that is analyzed. Each read is checked for the presence of rRNA sequence k-mers.

`--rRNA` specifies a FASTA input file that contains rRNA sequences of the species of interest. These sequences are processed to extract k-mers.

`--out` specifies an output folder. This folder will be created if it does not exist already.

`--kmer` specifies the size of k-mers used to screen reads for rRNA similarity. Default value is 21.

`--cutoff` specifies the number of k-mers that need to be detected in a read in order to classify it as rRNA read. Default is 3.


Example:

![rRNA content](https://github.com/bpucker/RNAseqQualCheck/blob/main/rRNA_content_example.png?raw=true)


## RNA-seq coverage across transcripts analysis ##
This script analyzes the coverage of RNA-seq reads across transcripts. The results allow conclusions about the quality of the processed sample.




![RNA-seq coverage across transcripts](https://github.com/bpucker/RNAseqQualCheck/blob/main/RNAseq_coverage_across_transcripts.png?raw=true)





