#!/bin/bash

# Usage check
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <fasta_file> <sample_list_file> <output_prefix>"
    exit 1
fi

# Input variables
msa="$1"
sample_list="$2"
output_prefix="$3"

# Output file names
subset_msa="${output_prefix}_subset.fasta"
output_vcf="${output_prefix}_snp_only.vcf"
output_fasta="${output_prefix}_snp_only.fasta"

awk -v samples_file="$sample_list" '
BEGIN {
    while ((getline line < samples_file) > 0) {
        samples[line] = 1
    }
    close(samples_file)
}
{
    if (substr($0,1,1)==">") {
        header=$0
        seq_name = substr($0,2)
        print_flag = (seq_name in samples)
        if (print_flag) { print header }
    } else {
        if (print_flag) { print $0 }
    }
}' "$msa" > "$subset_msa"

snp-sites -v "$subset_msa" > "$output_vcf" 
snp-sites "$subset_msa" > "$output_fasta"
