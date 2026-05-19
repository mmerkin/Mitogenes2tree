#!/bin/bash
# Haplophylo v1.2
# Usage: bash script.sh -r reference.fa -i /path/to/bams -o output_prefix -f depth.tsv
set -euo pipefail

# Set variables

log_file="logfile_haplophylo.txt"
> "$log_file"
BQ=20
MQ=20
N_THRESHOLD=20
reference=""
annotation=""
input_path=""
output_prefix=""
sample_file=""
min_depth=""
min_allele_depth=""
allele_balance=""


# Send EVERYTHING (stdout + stderr) to log file
exec >> "$log_file" 2>&1

# Functions

show_help() {
  echo "Usage: $0 -r <reference> -a <annotation> -i <input_path> -o <output_prefix> -f <sample_file> [-s 10] [-d 20] [-b 0.2]"
  echo
  echo "Options:"
  echo "  -r    Mitochondrial reference genome (required)"
  echo "  -a    Annotation bed file (required)"
  echo "  -i    Path to bam files (required)"
  echo "  -o    Output prefix (required)"
  echo "  -f    File containing sample names (required) and depths (optional)"
  echo "  -s    Minimum depth to keep sample (optional)"  
  echo "  -d    Minimum mean depth for genotype call (optional)"  
  echo "  -b    Maximum allele balance allowed for genotype call (optional)"
  echo "  -h    Show help message"
}

log() {
    local msg="$1"
    local dest="${2:-both}"

    case "$dest" in
        terminal)
            echo -e "$msg" > /dev/tty
            ;;
        log)
            echo -e "$msg" >> "$log_file"
            ;;
        both)
            echo -e "$msg" > /dev/tty
            echo -e "$msg" >> "$log_file"
            ;;
        *)
            echo "Invalid log destination: $dest" >&2
            ;;
    esac
}

print_banner() {
clear > /dev/tty
cat << 'EOF' > /dev/tty
          _   _    _    ____  _     ___  ____  _   ___   ___     ___  
         | | | |  / \  |  _ \| |   / _ \|  _ \| | | \ \ / / |   / _ \ 
         | |_| | / _ \ | |_) | |  | | | | |_) | |_| |\ V /| |  | | | |
         |  _  |/ ___ \|  __/| |__| |_| |  __/|  _  | | | | |__| |_| |
         |_| |_/_/   \_\_|   |_____\___/|_|   |_| |_| |_| |_____\___/ 
                                                                  
EOF
}


# Arguments

while getopts "r:a:i:o:f:s:d:b:h" opt; do
  case $opt in
    r) reference="$OPTARG" ;;
    a) annotation="$OPTARG" ;;
    i) input_path="$OPTARG" ;;
    o) output_prefix="$OPTARG" ;;
    f) sample_file="$OPTARG" ;;
    s) min_depth="$OPTARG" ;;
    d) min_allele_depth="$OPTARG" ;;
    b) allele_balance="$OPTARG" ;;
    h)
      show_help > /dev/tty
      exit 0
      ;;
    \?)
      log "Invalid option: -$OPTARG"
      show_help > /dev/tty
      exit 1
      ;;
    :)
      log "Option -$OPTARG requires an argument"
      show_help > /dev/tty
      exit 1
      ;;
  esac
done

trap 'log "Command failed: $BASH_COMMAND" terminal' ERR


if [[ -z "$reference" || -z "$annotation" || -z "$input_path" || -z "$output_prefix" || -z "$sample_file" ]]; then
  log "ERROR: Missing required arguments"
  show_help > /dev/tty
  exit 1
fi

missing=()
for cmd in bcftools bedtools samtools iqtree snp-sites; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        missing+=("$cmd")
    fi
done

if ((${#missing[@]} > 0)); then
    log "There are missing commands. Have you activated the conda environment?"
    exit 1
fi


run_cmd() {
    local cmd="$*"
    echo "$cmd" >> "$log_file"
    eval "$cmd"
}


# Start

print_banner

log "Script started at $(date)" log
log "Starting run of Mitogenes2tree\nPlease report any errors as an issue on github\n"
log "Reference: $reference"
log "Annotation file: $annotation"
log "Input path: $input_path"
log "Output prefix: $output_prefix"
log "Sample file: $sample_file"
[[ -n "$min_depth" ]] && log "Minimum sample depth: $min_depth"
[[ -n "$min_allele_depth" ]] && log "Minimum allele depth: $min_allele_depth"
[[ -n "$allele_balance" ]] && log "Maximum allele balance: $allele_balance"
log "\nBegin analysis"

# Find samples

mkdir -p ${output_prefix}_consensus


cols=$(awk 'NR>1 {print NF; exit}' "$sample_file")

if [[ -n "$min_depth" && "$cols" -lt 2 ]]; then
  log "ERROR: The sample file only contains one column. A 2-column tsv file is required to filter by sample depth with the -s option"
  exit 1
fi

if [[ -z "$min_depth" ]]; then
  mapfile -t SAMPLES < <(awk 'NR>1 {print $1}' "$sample_file")
else
  mapfile -t SAMPLES < <(awk -v min_dp="$min_depth" 'NR>1 && $2 > min_dp {print $1}' "$sample_file")
fi

retained=${#SAMPLES[@]}
total=$(awk 'NR>1' "$sample_file" | wc -l)

log "Using $retained samples out of $total"

for sample in "${SAMPLES[@]}"; do
    sample_base=$(basename "$sample" .bam)
    printf "sample_base='%s'\n" "$sample_base"
done

# Call variants


BAMS=""
> readgroup_map.txt
for sample in "${SAMPLES[@]}"; do
    bam="$input_path/$sample"
    if [ -f "$bam" ]; then
        BAMS="$BAMS $bam"
    fi
    sample_base=$(basename "$sample" .bam)
    RG=$(samtools view -H $bam | grep '^@RG'| sed -n 's/.*SM:\([^[:space:]]*\).*/\1/p')
    echo -e "$RG $sample_base" >> readgroup_map.txt
done

if [ -z "$BAMS" ]; then
    log "ERROR: No bam files available for variant calling"
    exit 1
fi

run_cmd bcftools mpileup \
-f "$reference" \
-q $MQ -Q $BQ \
-G readgroup_map.txt \
-a AD,DP \
-Ou \
$BAMS | \
bcftools call \
--ploidy 1 \
-mv \
-Oz \
-o "variants.raw.vcf.gz"
run_cmd bcftools index variants.raw.vcf.gz

# Remove indels
run_cmd bcftools view -v snps -Oz -o variants.snps.vcf.gz variants.raw.vcf.gz
run_cmd bcftools index variants.snps.vcf.gz

snp_count=$(bcftools view -H variants.snps.vcf.gz | wc -l)
log "Called $snp_count SNPs"

# Make gene counts output file

run_cmd bedtools intersect -a "$annotation" -b variants.snps.vcf.gz -c > ${output_prefix}_gene_variant_counts.tsv

# Mask sequences with low depth or high allele imbalance

total_masked=0
processed_samples=0

for sample in "${SAMPLES[@]}"; do
BAM="$input_path/$sample"
sample_base=$(basename "$sample" .bam)
MASK="${output_prefix}_consensus/${sample_base}_mask.bed"
> "$MASK"

if [[ -n "${min_allele_depth:-}" ]]; then
samtools depth -aa "$BAM" | awk -v min_dp="$min_allele_depth" '
$3 < min_dp {
print $1 "\t" ($2-1) "\t" $2
}' >> "$MASK"
fi

if [[ -n "${allele_balance:-}" ]]; then
bcftools query -f '%CHROM\t%POS[\t%GT][\t%AD]\n' -s "$sample_base" "variants.snps.vcf.gz" |
awk -v thresh=$allele_balance -F'\t' '{
gt = $3
split($4, ad, ",")
for(i=1;i<=length(ad);i++){
if(ad[i]=="") ad[i]=0
}
total = 0
for(i=1;i<=length(ad);i++) total += ad[i]
sum_other = 0
for(i=1;i<=length(ad);i++){
if(i-1 != gt) sum_other += ad[i]  # i-1 because GT is 0-based
}
if(total>0){
ab = sum_other / total
if(ab > thresh){
print $1"\t"($2-1)"\t"$2  # BED: 0-based start, 1-based end
}
}
}' >> "$MASK"
fi  
  
if [[ -s "$MASK" ]]; then
masked_sites=$(wc -l < "$MASK")
else
masked_sites=0
fi
echo "Sample $sample_base masked sites: $masked_sites"
total_masked=$((total_masked + masked_sites))
processed_samples=$((processed_samples + 1))
sort -k1,1 -k2,2n -o "$MASK" "$MASK"
done


if [[ "$processed_samples" -gt 0 ]]; then
avg_masked=$((total_masked / processed_samples))
log "Average masked sites per sample: $avg_masked"
else
log "No samples processed"
fi


# Create consensus sequence and extract gene sequences

for sample in "${SAMPLES[@]}"; do
sample_base=$(basename "$sample" .bam)
bcftools consensus \
-f "$reference" \
-s "$sample_base" \
-m "${output_prefix}_consensus/${sample_base}_mask.bed" \
variants.snps.vcf.gz \
> "${output_prefix}_consensus/${sample_base}_full_mt_consensus.fasta"

bedtools getfasta -fi "${output_prefix}_consensus/${sample_base}_full_mt_consensus.fasta" \
-bed "$annotation" -nameOnly -fo - | \
awk -v sample="$sample_base" 'BEGIN{RS=">"; ORS=""} NR>1 {
split($0, lines, "\n")
gene = lines[1]
seq = ""
for(i=2;i<=length(lines);i++) seq = seq lines[i]
print ">" sample "|" gene "\n" seq "\n"
}' > "${output_prefix}_consensus/${sample_base}_genes_consensus.fasta"
done


# Create parition file and concatenated gene file


concat_sequence="${output_prefix}_all_sequences.fa"
partition_file="${output_prefix}_partitions.nex"
> "$concat_sequence"
> "$partition_file"
pos=1

first_sample="${SAMPLES[0]}"  # Use first sample to get gene lengths
first_sample=$(basename "$first_sample" .bam)

# Read gene names from BED

gene_list=()
while read -r chrom start end gene; do
gene_list+=("$gene")
done < "$annotation"

# Initialize per-sample concatenated strings

declare -A concat_seqs
for sample in "${SAMPLES[@]}"; do
sample_base=$(basename "$sample" .bam)
concat_seqs[$sample_base]=""
done

# Loop through genes
for gene in "${gene_list[@]}"; do
seq=$(awk -v sample="$first_sample" -v gene="$gene" '
BEGIN {found=0; seq=""}
/^>/ {
gsub(/^ +| +$/,"",$0)
if ($0 ~ "^>" sample "\\|" gene "$") {found=1} else {found=0}
next
}
found {seq=seq $0}
END {print seq}' "${output_prefix}_consensus/${first_sample}_genes_consensus.fasta")
seq_len=$(echo -n "$seq" | tr -d '\n' | wc -c)

if [ "$seq_len" -eq 0 ]; then
log "Warning: gene $gene sequence length is 0 in sample $first_sample"
exit 1
fi

end_pos=$((pos + seq_len - 1))
echo "DNA, $gene = $pos-$end_pos" >> "$partition_file"

# Append sequences for all samples
for sample in "${SAMPLES[@]}"; do
sample_base=$(basename "$sample" .bam)
gene_seq=$(awk -v sample_base="$sample_base" -v gene="$gene" '
BEGIN {found=0; seq=""}
/^>/ {
gsub(/^ +| +$/,"",$0)
if ($0 ~ "^>" sample_base "\\|" gene "$") {found=1} else {found=0}
next
}
found {seq=seq $0}
END {print seq}' "${output_prefix}_consensus/${sample_base}_genes_consensus.fasta")
concat_seqs[$sample_base]="${concat_seqs[$sample_base]}$gene_seq"
done
pos=$((end_pos + 1))
done

# Write concatenated fasta
for sample in "${SAMPLES[@]}"; do
sample_base=$(basename "$sample" .bam)
echo ">$sample_base" >> "$concat_sequence"
echo "${concat_seqs[$sample_base]}" >> "$concat_sequence"
done

# Final sanity check

ALN_LEN=$(awk '/^[^>]/ {print length($0); exit}' "$concat_sequence")
LAST_PART=$(tail -n1 "$partition_file" | awk -F"=" '{print $2}' | tr -d ' ' | awk -F"-" '{print $2}')

if [ "$ALN_LEN" -ne "$LAST_PART" ]; then
    log "ERROR: alignment length ($ALN_LEN) does not match last partition end ($LAST_PART). Please report this bug."
    exit 1
fi

GENOME_LEN=$(awk '!/^>/ {total += length($0)} END {print total}' "$reference")
PERCENTAGE=$(awk -v u="$ALN_LEN" -v g="$GENOME_LEN" 'BEGIN {printf "%.2f", (u/g)*100}')

log "\nAlignment length: $ALN_LEN (${PERCENTAGE}% of genome used)"

# snp sites only

run_cmd snp-sites "$concat_sequence" > ${output_prefix}_variants_only_alignment.fa

# Check Ns

OUTFILE="${output_prefix}_high_missing.tsv"

: > "$OUTFILE"

read AVG_N COUNT <<EOF
$(awk -v thresh="$N_THRESHOLD" -v out="$OUTFILE" '
BEGIN {
  OFS="\t"
  print "sample_id", "N_count" > out
}

/^>/ {
  if (seqs > 0) {
    total += n
    if (n > thresh) {
      print id, n >> out
      count++
    }
  }
  id = substr($0,2)
  seqs++
  n = 0
  next
}

{
  n += gsub(/[Nn]/, "&")
}

END {
  total += n
  if (n > thresh) {
    print id, n >> out
    count++
  }

  if (seqs > 0)
    print total / seqs, count
  else
    print 0, 0
}
' "$concat_sequence")
EOF

log "Average number of Ns: $AVG_N"

if (( COUNT > 0 )); then
  log "WARNING: $COUNT samples found with more than $N_THRESHOLD missing sites."
  log "See $OUTFILE for details"
fi

# IQTree

mkdir -p ${output_prefix}_tree

run_cmd iqtree3 -s "$concat_sequence" -p "$partition_file" -m MFP -bb 1000 -T AUTO -pre ${output_prefix}_tree/${output_prefix}_tree

cp ${output_prefix}_tree/${output_prefix}_tree.treefile ${output_prefix}_tree.treefile

log "Run finished successfully. Find output tree in ${output_prefix}_tree.treefile"

log "Script finished at $(date)" log

# End of script
