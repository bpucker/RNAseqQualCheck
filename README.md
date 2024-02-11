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


```
Usage1
python3 read_distr_checker.py --bam <FILE> --gff <FILE> --out <DIR>

Usage2
python3 read_distr_checker.py --cov <FILE> --gff <FILE> --out <DIR>

Mandatory:
--bam        STR   BAM input file    |  --cov   STR    Coverage input file
--gff        STR   GFF input file
--out        STR   Output folder

Optional:
--sample     STR    Sample name
--samtools   STR    Path to samtools [samtools]
--bedtools   STR    Path to bedtools [bedtools]
--chunks     INT    Number of chunks [100]
--minexpcut  INT    Minimal coverage [100]
```

`--bam` specifies a BAM input file. This will be converted into a coverage file (COV) to analyze the distribution of reads across transcripts. This argument can also be used to provide a comma-separated list of files for automatic processing of large batches of files.
`--cov` specifies a COV input file. This file is the basis to analyze the distribution of reads across transcripts. This argument can also be used to provide a comma-separated list of files for automatic processing of large batches of files.
`--gff` specifies a GFF input file. This file is processed to identify the positions of exons that form an intron. 
`--out` specifies an output folder. If this folder does not exist, it will be created.

`--sample` specifies sample names. This argument can also be used to provide a comma-separated list of sample names for automatic processing of large batches of files. The length of this list should match the number of BAM/COV files provided.
`--samtools` specifies the full path to samtools. Default: samtools.
`--bedtools` specifies the full path to genomeCoverageBed. Default: genomeCoverageBed.
`--chunks` specifies the number of chunks to create for each transcript when assessing the coverage. Default: 100. Only transcripts with at least this length will be considered for the analysis.
`--minexpcut` specifies the minimal total number of sequenced basis. Only transcripts with at least this number of sequenced bases is considered for the analysis.


![RNA-seq coverage across transcripts](https://github.com/bpucker/RNAseqQualCheck/blob/main/RNAseq_coverage_across_transcripts.png?raw=true)



## Dependencies ##
[Python3](https://www.python.org/downloads/) with standard modules is required to run this analysis. Plotting is based on matplotlib and seaborn.
```
sudo apt update
sudo apt install python3
```

[Samtools](http://www.htslib.org/) is needed to handle BAM files.
```
sudo apt update
sudo apt install samtools
```


[Bedtools](https://bedtools.readthedocs.io/en/latest/) is needed to calculate coverage values.
```
sudo apt update
sudo apt install bedtools
```


