#!/bin/bash
# Usage: bash script.sh species reference depth_file

species=$1 # e.g. B1_Cupido_minimus
reference=$2 # e.g. Cupido_minimus_MIQ.fa
depth_file=$3 # e.g. B1_all_mtDNA_depths.tsv

# Frequently changed variables

min_depth=3
min_mean_depth=10
allele_balance=0.2


# File paths

mito_dir=~/Velocity/mitochondria/$species/gene_sequences
REF=~/Velocity/mitochondria/reference/reordered/$reference
BED=~/Velocity/mitochondria/reference/annotation/$species.bed
depth_path=~/Velocity/mitochondria/$species/qc/depths/$depth_file
logfile=$mito_dir/${species}_call_mtDNA_gene_sequences.log
consensus_dir=$mito_dir/consensus_genes
concat_sequence=$mito_dir/${species}_mtDNA_all_genes_aligned.fa
partition_file=$mito_dir/${species}_mtDNA_all_genes_partitions.nex

# Conda

conda=~/miniconda3/bin/
mtDNA=mtDNA
source $conda/activate $mtDNA

# Determine whether variables have remained unset

if [[ -z $species ]] || [[ -z $reference ]] || [[ -z $depth_file ]]; then
echo "Usage: bash script.sh species reference depth_file"
exit 1
fi

# Send errors to log file

exec 3>&1 1>"$logfile" 2>&1
trap "echo 'ERROR: An error occurred during execution, check $logfile for details.' >&3" ERR
trap '{ set +x; } 2>/dev/null; echo -n "[$(date -Is)]  "; set -x' DEBUG
set -e

printf "\n\n\n" | tee /dev/fd/3
printf "Starting analysis...\n" | tee /dev/fd/3

# CREATE OUTPUT DIRECTORIES

mkdir -p "$consensus_dir"

# Select samples with high enough depth

mapfile -t SAMPLES < <(awk -v min_dp="$min_mean_depth" 'NR>1 && $2 > min_dp {print $1}' "$depth_path")

retained=${#SAMPLES[@]}
total=$(awk 'NR>1' "$depth_path" | wc -l)

printf "Using %d samples out of %d\n" "$retained" "$total" | tee /dev/fd/3

# Call variants

printf "Calling variants...\n" | tee /dev/fd/3

BAMS=""
for sample in "${SAMPLES[@]}"; do
    bam="$mito_dir/bams/${sample}.realn.bam"
    if [ -f "$bam" ]; then
        BAMS="$BAMS $bam"
    fi
done

if [ -z "$BAMS" ]; then
    printf "Error: No BAMs available for variant calling\n" | tee /dev/fd/3
    exit 1
fi

bcftools mpileup \
-f "$REF" \
-q 20 -Q 20 \
-a AD,DP \
-Ou \
$BAMS | \
bcftools call \
--ploidy 1 \
-mv \
-Oz \
-o "$mito_dir/${species}_mtDNA.raw.vcf.gz"
bcftools index $mito_dir/${species}_mtDNA.raw.vcf.gz

echo "Filtering VCF to SNPs only..."
bcftools view -v snps -Oz -o $mito_dir/${species}_mtDNA.snps.vcf.gz $mito_dir/${species}_mtDNA.raw.vcf.gz
bcftools index $mito_dir/${species}_mtDNA.snps.vcf.gz

snp_count=$(bcftools view -H -v snps $mito_dir/${species}_mtDNA.snps.vcf.gz | wc -l)
printf "Total SNPs: %d\n" "$snp_count" | tee /dev/fd/3

# Mask sequences with low depth or high allele imbalance

printf "Masking regions...\n" | tee /dev/fd/3

for sample in "${SAMPLES[@]}"; do
BAM="$mito_dir/bams/${sample}.realn.bam"
MASK="${consensus_dir}/${sample}_mask.bed"
> "$MASK"
samtools depth -aa "$BAM" | awk -v min_dp=$min_depth '{if($3<min_dp) print $1"\t"$2-1"\t"$2}' >> "$MASK"

#bcftools query -f '%CHROM\t%POS[\t%AD]\n' -s "$sample" $mito_dir/${species}_mtDNA.snps.vcf.gz | \
#awk -v thresh=$allele_balance -F'\t' '{
#split($3, ad, ",")
#ad_ref = ad[1]
#ad_alt = ad[2]
#total = ad_ref + ad_alt
#if(total>0){
#minor = (ad_ref < ad_alt) ? ad_ref : ad_alt
#ab = minor / total
#if(ab > thresh){
#print $1"\t"($2-1)"\t"$2 # subtract 1 as bed file is 0-based
#}
#}
#}' >> "$MASK"



bcftools query -f '%CHROM\t%POS[\t%GT][\t%AD]\n' -s "$sample" "$mito_dir/${species}_mtDNA.snps.vcf.gz" | \
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

sort -k1,1 -k2,2n -o "$MASK" "$MASK"
done


# Create consensus sequence and extract gene sequences

printf "Creating per-gene consensus sequences...\n" | tee /dev/fd/3

for sample in "${SAMPLES[@]}"; do
bcftools consensus \
-f "$REF" \
-s "$sample" \
-m "${consensus_dir}/${sample}_mask.bed" \
$mito_dir/${species}_mtDNA.snps.vcf.gz \
> "${consensus_dir}/${sample}_full_mt_consensus.fasta"

bedtools getfasta -fi "${consensus_dir}/${sample}_full_mt_consensus.fasta" \
-bed "$BED" -nameOnly -fo - | \
awk -v sample="$sample" 'BEGIN{RS=">"; ORS=""} NR>1 {
split($0, lines, "\n")
gene = lines[1]
seq = ""
for(i=2;i<=length(lines);i++) seq = seq lines[i]
print ">" sample "|" gene "\n" seq "\n"
}' > "${consensus_dir}/${sample}_genes_consensus.fasta"
done


# CONCATENATE PER-GENE SEQUENCES AND CREATE PARTITION FILE

printf "Concatenating sequences and building partition file...\n" | tee /dev/fd/3
> "$concat_sequence"
> "$partition_file"
pos=1

first_sample="${SAMPLES[0]}"  # Use first sample to get gene lengths

# Read gene names from BED

gene_list=()
while read -r chrom start end gene; do
gene_list+=("$gene")
done < "$BED"

# Initialize per-sample concatenated strings

declare -A concat_seqs
for sample in "${SAMPLES[@]}"; do
concat_seqs[$sample]=""
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
END {print seq}' "${consensus_dir}/${first_sample}_genes_consensus.fasta")
seq_len=$(echo -n "$seq" | tr -d '\n' | wc -c)

if [ "$seq_len" -eq 0 ]; then
echo "Warning: gene $gene sequence length is 0 in sample $first_sample"
exit 1
fi

end_pos=$((pos + seq_len - 1))
echo "DNA, $gene = $pos-$end_pos" >> "$partition_file"

# Append sequences for all samples
for sample in "${SAMPLES[@]}"; do
gene_seq=$(awk -v sample="$sample" -v gene="$gene" '
BEGIN {found=0; seq=""}
/^>/ {
gsub(/^ +| +$/,"",$0)
if ($0 ~ "^>" sample "\\|" gene "$") {found=1} else {found=0}
next
}
found {seq=seq $0}
END {print seq}' "${consensus_dir}/${sample}_genes_consensus.fasta")
concat_seqs[$sample]="${concat_seqs[$sample]}$gene_seq"
done
pos=$((end_pos + 1))
done

# Write concatenated fasta
for sample in "${SAMPLES[@]}"; do
echo ">$sample" >> "$concat_sequence"
echo "${concat_seqs[$sample]}" >> "$concat_sequence"
done

printf "Gene sequence calling completed.\n\n\n" | tee /dev/fd/3

# Final sanity check

ALN_LEN=$(grep -v '^>' "$concat_sequence" | head -n1 | tr -d '\n' | wc -c)
LAST_PART=$(tail -n1 "$partition_file" | awk -F"=" '{print $2}' | tr -d ' ' | awk -F"-" '{print $2}')

printf "Alignment length: $ALN_LEN\n" | tee /dev/fd/3
printf "Last partition end: $LAST_PART\n" | tee /dev/fd/3
