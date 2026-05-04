# Mitogenes2tree

<img width="1024" height="1024" alt="image" src="https://github.com/user-attachments/assets/ff8c5cb7-1123-4dad-9aeb-d0b33f3954d9" />
Image credit: Laura Jimenez Burney


# Getting started

## Requirements:
- reference genome 
- mitos bed file
- bam file or fastq reads that are either paired end or merged (museum)

## Installing

```bash
git clone X
cd X
conda env create -f environment.yml
```
## Running

A detailed description of how to generate the input files is provided here.

To run:

```bash
bash mitogenes2tree.sh -r <reference> -a <annotation> -i <input_path> -o <output_prefix> -f <sample_file> [-s 10] [-d 20] [-b 0.2]"

Parameters:
-r    Mitochondrial reference genome (required)
-a    Annotation bed file (required)
-i    Path to bam files (required)
-o    Output prefix (required)
-f    File containing sample names (required) and depths (optional)
-s    Minimum depth to keep sample (optional)
-d    Minimum mean depth for genotype call (optional)  
-b    Maximum allele balance allowed for genotype call  (optional)
-h    Show help message
```

## Filters

There are three optional filters that can be performed:
- If a depth file is provided, samples with a depth lower than a threshold can be exlcuded (-s)
- Genotype calls can be masked (set to N) if the allele depth is below a certain threshold (-d)
- If the sample is heterozygous at a position, genotypes can be set to N if allele balance is above a threshold (-b). i.e. setting -b 0.2 will mask any sites where more than 20% of the reads support an alternate allele.

# Output files

## To do:

Conda environment


Clear files from previous runs
