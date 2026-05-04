# Mitogenes2tree

<img width="1024" height="1024" alt="image" src="https://github.com/user-attachments/assets/ff8c5cb7-1123-4dad-9aeb-d0b33f3954d9" />
Image credit: Laura Jimenez Burney


# Description

Produce a tree of concatenated mitochondrial DNA genes from samples within the same species.

# Getting started

## Requirements:
- reference genome 
- mitos bed file
- bam file or fastq reads that are either paired end or merged (museum)

## Dependencies

The main script only requires the following tools:

- bcftools
- bedtools
- samtools
- iqtree

However, the bonus scripts may require additional software, which are listed as a comment at the top.

Alternatively, a conda environment containing all required packages is attached.

## Installing

```bash
git clone https://github.com/mmerkin/Mitogenes2tree.git
cd Mitogenes2tree
conda env create -f environment.yml
conda activate mitogenes2tree
```
## Running

A detailed description of how to generate the input files is provided [here](Instructions/Input_files.md)

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

- Tree in newick format: "_tree.treefile"
- Concatenated mitogene sequences: "_all_sequences.fa"
- Consenus sequences folder containing the fasta gene sequences ("_genes_consensus.fasta"), the concatenated gene sequences ("_full_mt_consensus.fasta") and the bed files of positions masked during filtering ("_mask.bed")
- Tree folder containing the extra outputs of iqtree
- vcf files containing raw variants and with indels removed
- A list of samples with a high number of Ns ("_high_missing.tsv")
- log file ("logfile_mitogenes2tree.txt")
