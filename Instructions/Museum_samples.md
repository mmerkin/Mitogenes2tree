## Museum samples

This section is a WIP. Briefly, museum samples should be collapsed to form a single merged set of reads and then mapped with bwa aln instead of mem2.

It is recommended to first collapse paired end sequences with fastp:

```bash
input_path=
output_path=

mkdir -p $output_path/qc/discard
mkdir -p $output_path/qc/reports
mkdir -p $output_path/merged

for r1 in "$input_path"/*R1*.fastq.gz; do
filetag=$(basename "$r1" .veladapt.clean_R1.fastq.gz)
echo "Collapsing $filetag"
fastp -m -A \
-i "$input_path/${filetag}.veladapt.clean_R1.fastq.gz" \
-I "$input_path/${filetag}.veladapt.clean_R2.fastq.gz" \
-o "$output_path/qc/discard/${filetag}.fastp.discard_R1.fastq.gz" \
-O "$output_path/qc/discard/${filetag}.fastp.discard_R2.fastq.gz" \
--merged_out "$output_path/merged/${filetag}.fastp.merged.fastq.gz" \
-j "$output_path/qc/reports/${filetag}_fastp_report.json" \
-h "$output_path/qc/reports/${filetag}_fastp_report.html" \
2> "$output_path/qc/reports/${filetag}_fastp_report.txt"
done
```
The script Map_museum_reads.sh can then be run in a similar manner to map_reads_parallel.sh
