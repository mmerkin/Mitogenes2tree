#!/bin/bash
# Usage: bash script.sh reference.fa read_path output_path

set -euo pipefail

# Variables to change

threads=8 # Number of threads to use per sample
jobs=4 # Number of samples to run in parallel
remove_temp=false # Whether to delete intermediate files true|false

gatk=false # Whether to realign around indels with GATK true|false
GATK38=GATK38 # Name of conda environment containing GATK3.8

# Check variables are set and environment is activated

if [[ $# -ne 3 ]]; then
    echo "Usage: bash $0 reference.fa read_path output_path" >&2
    exit 1
fi

reference="$1"
read_path="$2"
output_path="$3"

mkdir -p "$output_path/bams"
mkdir -p "$output_path/logs"
temp="${output_path}/temp"


missing=()
for cmd in bwa-mem2 samtools bam parallel; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        missing+=("$cmd")
    fi
done

if ((${#missing[@]} > 0)); then
    echo "There are missing commands. Have you activated the conda environment?" >&3
    exit 1
fi


# Map and filter reads function

process_sample() {
    file="$1"

    sample_name="${file##*/}"
    sample_name="${sample_name%%R1*}"
    sample_name="${sample_name%[_\.]}"

    log="$output_path/logs/${sample_name}.log" # Create log file

    {
        echo "[$(date -Is)] START $sample_name"

        sample_temp="$temp/$sample_name"
        mkdir -p "$sample_temp"

        R1="$file"
        R2="${file/R1/R2}"
        RG="@RG\tID:${sample_name}\tSM:${sample_name}\tPL:ILLUMINA"

        # Map reads and apply basic filters
        bwa-mem2 mem -t "$threads" -R "$RG" "$reference" "$R1" "$R2" | \
        samtools view -@ "$threads" -b -F 4 -q 20 | \
        samtools sort -n -@ "$threads" -o "$sample_temp/${sample_name}.sort.n.bam"

        echo "Mapped reads successfully"

        # Remove duplicate reads
        samtools fixmate -m -@ "$threads" \
            "$sample_temp/${sample_name}.sort.n.bam" \
            "$sample_temp/${sample_name}.fixmate.bam"

        samtools sort -@ "$threads" \
            -o "$sample_temp/${sample_name}.sort.p.bam" \
            "$sample_temp/${sample_name}.fixmate.bam"

        samtools markdup -r -@ "$threads" \
            "$sample_temp/${sample_name}.sort.p.bam" \
            "$sample_temp/${sample_name}.dedup.bam"

        samtools index "$sample_temp/${sample_name}.dedup.bam"

        echo "Removed PCR duplicates successfully"

        # If GATK is not being used, the pipeline stops after clipping overlaps
        if $gatk; then
            clip_bam="$sample_temp/${sample_name}.noclip.bam"
        else
            clip_bam="$output_path/bams/${sample_name}.bam"
        fi

        bam clipOverlap \
            --in "$sample_temp/${sample_name}.dedup.bam" \
            --out "$clip_bam"

        samtools index "$clip_bam"

        echo "Clipped overlaps successfully"

        # GATK realignment
        if $gatk; then
            intervals="$sample_temp/${sample_name}.intervals"

            conda run -n "$GATK38" gatk3 \
                -T RealignerTargetCreator \
                -R "$reference" \
                -I "$clip_bam" \
                -o "$intervals"

            conda run -n "$GATK38" gatk3 \
                -T IndelRealigner \
                -R "$reference" \
                -I "$clip_bam" \
                -targetIntervals "$intervals" \
                -o "$output_path/bams/${sample_name}.bam"

            echo "Realigned reads successfully"
        fi

        if $remove_temp; then
            rm -rf "$sample_temp"
        fi

        echo "[$(date -Is)] DONE $sample_name"

    } >"$log" 2>&1
}

export -f process_sample
export threads reference output_path temp gatk GATK38 remove_temp

# Run job
find "$read_path" -name "*R1*" | \
parallel --bar -j "$jobs" process_sample {} \
    --joblog "$output_path/parallel_joblog.txt"

echo "All samples processed!"
