#!/bin/bash
# Usage: bash script.sh bedfile output *.bam

set -euo pipefail

if [[ $# -lt 3 ]]; then
    echo "Usage: bash $0 bedfile output *.bam" >&2
    exit 1
fi

bedfile=$1
output=$2
shift 2

{
    echo -e "Sample_id\tMean_depth"

    for i in "$@"; do
        filetag=$(basename "$i")
        filetag=${filetag%.realn.bam}

        depth=$(samtools depth -a -Q 20 -q 30 -b "$bedfile" "$i" \
            | awk '{sum+=$3} END { if (NR>0) printf "%.4f\n", sum/NR; else print 0 }')

        echo -e "${filetag}\t${depth}"
    done
} | tee "$output"
