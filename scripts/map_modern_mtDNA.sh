#!/bin/bash
# Usage: bash script.sh reference read_path output_path

# Variables to change:
threads=32
remove_temp=false

set -euo pipefail

if [[ $# -ne 3 ]]; then
    echo "Usage: $0 reference.fa read_path output_path" >&2
    exit 1
fi

reference="$1"
read_path="$2"
output_path="$3"

# Create output directory
mkdir -p "$output_path"
temp="${output_path}/temp"

# Send errors to log file

exec 3>&1 1>"map_modern_mtDNA_logfile.txt" 2>&1
trap "echo 'ERROR: An error occurred during execution, check map_modern_mtDNA_logfile.txt for details.' >&3" ERR
trap '{ set +x; } 2>/dev/null; echo -n "[$(date -Is)]  "; set -x' DEBUG

# Process the bam files

for file in $read_path/*R1*; do

# Set variables for next task
sample_name="${f%%R1*}"
sample_name="${sample_name%[_\.]}"
filepath="${read_path}/${sample_name}"
echo -e "\nMoving to sample $sample_name" | tee /dev/fd/3 # tee displays the echo output in the terminal now that stdout is sent to a log file
mkdir -p "$temp/$filetag"
mkdir -p "$output_path/$filetag"

R1="$file"
R2="${file/R1/R2}"
RG="@RG\tID:${sample_name}\tSM:${sample_name}\tPL:ILLUMINA"
bwa-mem2 mem -t $threads-R "$RG" $reference "$R1" "$R2" | \
samtools view -@ $threads -b -F 4 -q 20 -o | \
samtools sort -n -@ $threads -o "$temp/${sample_name}.sort.n.bam"

samtools fixmate -m -@ $threads "$temp/${sample_name}.sort.n.bam" "$temp/${sample_name}.fixmate.bam"

samtools sort -@ $threads -o "$temp/${sample_name}.sort.p.bam" "$temp/${sample_name}.fixmate.bam"

samtools markdup -r -@ $threads "$temp/${sample_name}.sort.p.bam" "$temp/${sample_name}.dedup.bam"

samtools index "$temp/${sample_name}.dedup.bam"

bamUtil clipOverlap --in "$temp/$filetag/$filetag.rmdup.bam" --out "$temp/$filetag/$filetag.noclip.bam"

samtools index "$temp/$filetag/$filetag.noclip.bam"

if $remove_temp; then
rm -r $temp
fi

done

echo -e "\nAll samples processed!" | tee /dev/fd/3
