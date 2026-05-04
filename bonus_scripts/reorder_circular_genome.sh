#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 genome.fasta gene.fasta" >&2
    exit 1
fi

GENOME="$1"
GENE="$2"

genome_seq=$(grep -v "^>" "$GENOME" | tr -d '\n')
genome_len=${#genome_seq}

if [[ ! -f "${GENOME}.nhr" ]]; then
    makeblastdb -in "$GENOME" -dbtype nucl >/dev/null
fi

# Perform blast search

BLAST_HIT=$(blastn -query "$GENE" -db "$GENOME" \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
    | sort -nrk12,12 | head -n1)

if [[ -z "$BLAST_HIT" ]]; then
    echo "ERROR: Gene did not align to the genome." >&2
    exit 1
fi

sstart=$(echo "$BLAST_HIT" | awk '{print $9}')
send=$(echo "$BLAST_HIT"   | awk '{print $10}')

if (( send < sstart )); then
    strand="-"
    genome_seq=$(echo "$genome_seq" | tr "ACGTacgt" "TGCAtgca" | rev)
    sstart0=$(( genome_len - sstart ))
else
    strand="+"
    sstart0=$(( sstart - 1 ))
fi

prefix=${genome_seq:0:$sstart0}
suffix=${genome_seq:$sstart0}
rotated="${suffix}${prefix}"

echo ">reordered_genome_forward"
echo "$rotated" | fold -w 60
