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
for cmd in bwa samtools bam parallel dedup realignsamfile; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        missing+=("$cmd")
    fi
done

if ((${#missing[@]} > 0)); then
    echo "There are missing commands. Have you activated the conda environment?" >&2
    exit 1
fi

circularised_reference="${reference%.*}_500.${reference##*.}"

if [[ -f "$circularised_reference" ]]; then
    echo "Circularised reference already exists."
else
    echo "Circularising reference"
    header=$(head -n1 "$reference" | sed 's/^>//')
    circulargenerator -e 500 -i "$reference" -s "$header"
    bwa index -a bwtsw "$circularised_reference"
fi

# Map and filter reads function

process_sample() {
    file="$1"

    sample_name="${file##*/}"
    sample_name="${sample_name%%.fastp.merged.fastq.gz}"

    log="$output_path/logs/${sample_name}.log" # Create log file

    {
        echo "[$(date -Is)] START $sample_name"

        sample_temp="$temp/$sample_name"
        mkdir -p "$sample_temp"

        RG="@RG\tID:${sample_name}\tSM:${sample_name}\tPL:ILLUMINA"

        # Map reads and apply basic filters
        bwa aln -l 16500 -n 0.01 -o 0 -t $threads "$circularised_reference" "$file" > "$sample_temp/${sample_name}.reads.sai"
        bwa samse -r "$RG" "$circularised_reference" "$sample_temp/${sample_name}.reads.sai" "$file" | samtools view -bS - > "$sample_temp/${sample_name}.raw.bam"
        realignsamfile -e 500 -i "$sample_temp/${sample_name}.raw.bam" -r "$reference"
        
        echo "Mapped reads successfully"

        samtools view -@ $threads -b -F 4 -q 20 "$sample_temp/${sample_name}.raw_realigned.bam" -o "$sample_temp/${sample_name}.filtered.bam"

mkdir -p "$sample_temp"/dedup
dedup -i "$sample_temp"/"${sample_name}".filtered.bam -m -u -o "$sample_temp"/dedup


        # If GATK is not being used, the pipeline stops after clipping overlaps
        if $gatk; then
            sort_bam="$sample_temp/${sample_name}.sorted.bam"
        else
            sort_bam="$output_path/bams/${sample_name}.bam"
        fi

        samtools sort -@ $threads "$sample_temp/dedup/${sample_name}.filtered_rmdup.bam" -o "$sort_bam" 2> /dev/null
        samtools index "$sort_bam"

        echo "Sorted bam successfully"

        # GATK realignment
        if $gatk; then
            intervals="$sample_temp/${sample_name}.intervals"

            conda run -n "$GATK38" gatk3 \
                -T RealignerTargetCreator \
                -R "$reference" \
                -I "$sort_bam" \
                -o "$intervals"

            conda run -n "$GATK38" gatk3 \
                -T IndelRealigner \
                -R "$reference" \
                -I "$sort_bam" \
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
export threads reference circularised_reference output_path temp gatk GATK38 remove_temp

# Run job
find "$read_path" -name "*fastq*" | \
parallel --bar -j "$jobs" process_sample {} \
    --joblog "$output_path/parallel_joblog.txt"

echo "All samples processed!"
jobs=4
remove_temp=false

gatk=true
GATK38=GATK38

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
    echo "There are missing commands. Have you activated the conda environment?" >&2
    exit 1
fi

circularised_reference="${reference%.*}_500.${reference##*.}"

if [[ -f "$circularised_reference" ]]; then
    echo "Circularised reference already exists."
else
    echo "Circularising reference"
    header=$(head -n1 "$reference" | sed 's/^>//')
    circulargenerator -e 500 -i "$reference" -s "$header"
    bwa index -a bwtsw "$circularised_reference"
fi

# Map and filter reads function

process_sample() {
    file="$1"

    sample_name="${file##*/}"
    sample_name="${sample_name%%.fastp.merged.fastq.gz}"

    log="$output_path/logs/${sample_name}.log" # Create log file

    {
        echo "[$(date -Is)] START $sample_name"

        sample_temp="$temp/$sample_name"
        mkdir -p "$sample_temp"

        RG="@RG\tID:${sample_name}\tSM:${sample_name}\tPL:ILLUMINA"

        # Map reads and apply basic filters
        bwa aln -l 16500 -n 0.01 -o 0 -t $threads "$circularised_reference" "$file" > "$sample_temp/${sample_name}.reads.sai"
        bwa samse -r "$RG" "$circularised_reference" "$sample_temp/${sample_name}.reads.sai" "$file" | samtools view -bS - > "$sample_temp/${sample_name}.raw.bam"
        realignsamfile -e 500 -i "$sample_temp/${sample_name}.raw.bam" -r "$reference"
        
        echo "Mapped reads successfully"

        samtools view -@ $threads -b -F 4 -q 20 "$sample_temp/${sample_name}.raw_realigned.bam" -o "$sample_temp/${sample_name}.filtered.bam"

mkdir -p "$sample_temp"/dedup
dedup -i "$sample_temp"/"${sample_name}".filtered.bam -m -u -o "$sample_temp"/dedup


        # If GATK is not being used, the pipeline stops after clipping overlaps
        if $gatk; then
            sort_bam="$sample_temp/${sample_name}.sorted.bam"
        else
            sort_bam="$output_path/bams/${sample_name}.bam"
        fi

        samtools sort -@ $threads "$sample_temp/dedup/${sample_name}.filtered_rmdup.bam" -o "$sort_bam" 2> /dev/null
        samtools index "$sort_bam"

        echo "Sorted bam successfully"

        # GATK realignment
        if $gatk; then
            intervals="$sample_temp/${sample_name}.intervals"

            conda run -n "$GATK38" gatk3 \
                -T RealignerTargetCreator \
                -R "$reference" \
                -I "$sort_bam" \
                -o "$intervals"

            conda run -n "$GATK38" gatk3 \
                -T IndelRealigner \
                -R "$reference" \
                -I "$sort_bam" \
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
export threads reference circularised_reference output_path temp gatk GATK38 remove_temp

# Run job
find "$read_path" -name "*fastq*" | \
parallel --bar -j "$jobs" process_sample {} \
    --joblog "$output_path/parallel_joblog.txt"

echo "All samples processed!"
