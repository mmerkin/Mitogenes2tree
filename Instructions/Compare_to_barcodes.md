The gene sequence outputs can be used to compare to existing datasets.

For example, there is a large dataset containing a 658bp barcoding region from over 20,000 butterfly samples. 

This is a walkthrough of how to make a full tree containing additional butterfly samples that you have obtained.

## Convert dataset to fasta

The dataset can only be downloaded in csv and json formats. Since the csv file has inconsistent formatting that leads to column errors, I used the json file to extract the fasta nucleotide individuals. I then filtered this to only retain samples that had exactly 658bp and renamed any samples with duplicate names.

```bash
cat all_COI_samples_metadata.json | awk '
{
    species=""; family=""; subfamily=""; country=""; sample=""; seq="";

    if (match($0, /"species":"([^"]+)"/, a)) species=a[1];
    if (match($0, /"family":"([^"]+)"/, b)) family=b[1];
    if (match($0, /"subfamily":"([^"]+)"/, c)) subfamily=c[1];
    if (match($0, /"country_iso":"([^"]+)"/, d)) country=d[1];
    if (match($0, /"sampleid":"([^"]+)"/, e)) sample=e[1];
    if (match($0, /"nuc":"([^"]+)"/, f)) seq=f[1];

    if (species && family && subfamily && country && sample && seq) {
        printf(">%s | %s | %s | %s | %s\n%s\n",
               species, family, subfamily, country, sample, seq);
    }
}
' | awk '
{
  if ($0 ~ /^>/) {
    sub(/^>/, "", $0)

    # split ONLY on first " | "
    n = index($0, " | ")

    header = (n ? substr($0, 1, n-1) : $0)
    rest   = (n ? substr($0, n+3)   : "")

    split(header, name, /[[:space:]]+/)

    printf(">%s_%s", name[1], name[2])

    if (rest != "")
      print " | " rest
    else
      print ""
  }
  else {
    print
  }
}
' > all_COI_samples.fa


seqkit seq -m 658 -M 658 all_COI_samples.fa > COI_samples_658bp.fa

seqkit rename COI_samples_658bp.fa > renamed_COI_samples_658bp.fa
```

## Extract barcode from mitogenes2tree output

The COI gene can be extracted from each sample using seqkit. After some experimentation, I found out that it is typically between 39-696, but a few species have small indels that changes this range. Perform an msa using e.g. clustalomega from one sample with a few of the 658bp barcode sequences, then find the range of the overlapping region.

```bash
fasta_files=consensus_genes/*genes_consensus.fasta
output=species_COI_barcodes.fa

seqkit grep -r -p "cox1" $fasta_files \
| seqkit subseq -r 39:696 \
| sed 's/N/-/g' > "$output"
```

Finally, use iqtree to make a tree

```bash
iqtree -s all_COI_barcodes.fa -T AUTO
```
