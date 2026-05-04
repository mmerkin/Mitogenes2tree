#!/bin/bash
# Usage: bash script.sh reference read_path output_path

# Variables to change:
threads=32
remove_temp=false

gatk=false
conda=~/miniconda3/bin/
GATK38=GATK38
mitogenes2tree=mitogenes2tree


set -euo pipefail

if [[ $# -ne 3 ]]; then
    echo "Usage: bash $0 reference.fa read_path output_path" >&2
    exit 1
fi

reference="$1"
read_path="$2"
output_path="$3"

# Create output directory
mkdir -p "$output_path/bams"
temp="${output_path}/temp"

# Send errors to log file

logfile="$output_path/modern_mapping_log.txt"

exec 3>&1
exec >"$logfile" 2>&1

missing=()

for cmd in bwa-mem2 bam; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        missing+=("$cmd")
    fi
done

if ((${#missing[@]} > 0)); then
    echo "There are missing commands. Have you activated the conda environment?" >&3
    exit 1
fi

trap 'echo "ERROR: An error occurred during execution, check '"$logfile"' for details." >&3' ERR
export PS4='[$(date -Is)] + '
trap 'set -x' DEBUG

# Command to output to terminal only

terminal() {
    echo "[$(date -Is)] $*" >&3
}

# Command to output to both terminal and log

all() {
    "$@" 2>&1 | tee /dev/fd/3
}


# Clear screen

printf '\n%.0s' {1..50} >&3
all echo -e "Mapping and filtering reads\n\n"


# Process the bam files


for file in $read_path/*R1*; do

# Set variables for next task

sample_name="${file##*/}"
sample_name="${sample_name%%R1*}"
sample_name="${sample_name%[_\.]}"
all echo -e "\nMoving to sample $sample_name"
mkdir -p "$temp/$sample_name"

R1="$file"
R2="${file/R1/R2}"
RG="@RG\tID:${sample_name}\tSM:${sample_name}\tPL:ILLUMINA"
bwa-mem2 mem -t $threads -R "$RG" $reference "$R1" "$R2" | \
samtools view -@ $threads -b -F 4 -q 20 | \
samtools sort -n -@ $threads -o "$temp/$sample_name/${sample_name}.sort.n.bam"

all echo -e "Mapped reads"


samtools fixmate -m -@ $threads "$temp/$sample_name/${sample_name}.sort.n.bam" "$temp/$sample_name/${sample_name}.fixmate.bam"

samtools sort -@ $threads -o "$temp/$sample_name/${sample_name}.sort.p.bam" "$temp/$sample_name/${sample_name}.fixmate.bam"

samtools markdup -r -@ $threads "$temp/$sample_name/${sample_name}.sort.p.bam" "$temp/$sample_name/${sample_name}.dedup.bam"

samtools index "$temp/$sample_name/${sample_name}.dedup.bam"

all echo -e "Removed PCR duplicates"

if $gatk; then
final_path="$temp/$sample_name/${sample_name}.noclip.bam"
else
final_path="$output_path/bams/${sample_name}.bam"
fi

bam clipOverlap --in "$temp/$sample_name/${sample_name}.dedup.bam" --out $final_path

samtools index "$final_path"

all echo -e "Clipped overlaps\n"

if $gatk; then

source $conda/activate $GATK38

gatk3 \
-T RealignerTargetCreator \
-R $REF \
-o "$temp/$sample_name/${sample_name}.intervals" \
-I "$temp/$sample_name/${sample_name}.noclip.bam"

gatk3 \
-T IndelRealigner \
-R $REF \
-targetIntervals "$temp/$sample_name/${sample_name}.intervals" \
-I "$temp/$sample_name/${sample_name}.noclip.bam" \
-o "$output_path/bams/${sample_name}.bam"

all echo -e "Realigned reads\n"

conda deactivate

source $conda/activate $mitogenes2tree

fi


if $remove_temp; then
rm -r $temp
fi

done

all echo -e "\nAll samples processed!"
