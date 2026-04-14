# Create input files for Mitogenes2tree

Mitogenes2tree requires the following three inputs:
- [Mitochondrial genome reference fasta file](https://github.com/mmerkin/Mitogenes2tree/blob/Update1/Instructions/1.Input_files.md#generate-relinearised-mitochondrial-reference-sequence)
- [Bed file with locations of each gene]()
- [Bam file of sequencing reads mapped to the reference]()

- [Optionally, a file of read depths can be provided for filtering]()


## Generate relinearised mitochondrial reference sequence

If a reference genome currently exists for the species of interest or a close relative, the mitochondrial sequence of the fasta file should be downloaded. Alternatively, a program such as [getOrganelle](https://github.com/Kinggerm/GetOrganelle) or [NOVOPlasty](https://github.com/ndierckx/NOVOPlasty) can be used to generate a reference sequence from short read data.

Since mitochondria possess circular genomes, the start point of this sequence will be random and could even be in the middle of a gene. As such, it is recommended to first relinearise the reference genome to start at a more appropriate position. The reference genome should be annotated to find the locations of all genes, visualised to choose an appropriate starting location and finally relinearised.

### Annotation with mitos2

If an annotation already exists on ncbi (identified by the page listing the locations of each gene), the file should be downloaded in genbank format. Otherwise, the fasta sequence can be annotated easily using mitos2 on the [galaxy webserver](https://usegalaxy.eu/). 

Create a galaxy account, upload the fasta file and run the tool MITOS2 to generate SEQ and nucleotide FASTA files, ensuring to also change the genetic code if necessary:

<img width="1173" height="381" alt="image" src="https://github.com/user-attachments/assets/4984a625-002c-4d09-a8db-76fb3428e00d" />

### Create genbank file for visualisation

The SEQ file produced by mitos2 can be converted to genbank file for visualisation. Download this file and rename it to ${species}_annotation.tbl
There are sometimes invalid gene codes, such as anticodons in brackets or tRNA-SER2/LEU2, so run this command to rename them:
```bash
sed -E 's/\([^)]*\)//g; s/\b(LEU|SER)[12]\b/\1/g' "${species}_annotation.tbl" > "${species}_mtDNA.tbl"
```
If you encounter any errors related such as "Qualifier had bad value", edit the tbl file to rename the affected gene(s).

An sbt file containing the header information is also required. This can be generated from the [ncbi submission portal](https://submit.ncbi.nlm.nih.gov/genbank/template/submission/) or [this example file](https://github.com/mmerkin/Mitogenes2tree/blob/Update1/Example_files/Example_template.sbt) can be used instead.

```bash
table2asn -i ${species}_mtDNA.fa -t genbank_template.sbt -f ${species}_mtDNA.tbl
asn2gb -i ${species}_mtDNA.sqn -o ${species}_mtDNA.gb
```

Finally, run the script update_genbank.awk to produce the genbank file

```bash
./update_genbank.awk $genus $species ${species}_mtDNA.gb > ${species}_mtDNA_updated.gb
```

### Visualisation

The genbank file can be visualised in R with chloroplot, which will come installed in the conda environment:

```bash
R
library(chloroplot)
table <- PlotTab("Aglais_urticae_mtDNA_updated.gb", T)
PlotMitGenome(table, organelle_type = F)
q()
```

Alterantively the [Chloroplot](https://irscope.shinyapps.io/Chloroplot/) or [OGDraw](https://chlorobox.mpimp-golm.mpg.de/OGDraw.html) webservers can be used instead.

Here is an example chloroplot output with the mitochondrial sequence of the butterfly *Papillio machaon*

<img width="1016" height="855" alt="image" src="https://github.com/user-attachments/assets/53243c0a-ffc5-4b74-921e-263ba7d35b8a" />

There is a large gap on the right that corresponds to a large non-coding region (D-loop). Traditionally, mitochondrial genomes are ordered such that the D-loop sequence is at the end, so the tRNA-M gene appears to be a good starting point.

### Relinearise sequence

Once a new starting gene has been chosen, a new fasta file should be created using the gene sequence from the mitos2 nucleotide fasta output. Then, the script reorder_mtDNA.sh can be run
```bash
bash reorder_mtDNA.sh ${species}_mtDNA.fa ${gene}.fa
```

## Generate bed file

If the reference genome has been reordered, mitos2 needs to be run again on the new reference genome. However, a bed file should be creating instead of seq and nucleotide fasta. 

The bed files should then be modified to extract just the first 4 columns and remove any tRNA anticodon sequences in brackets:
```bash
awk -F'\t' '{
for(i=1;i<=4;i++) gsub(/\([^)]*\)/,"",$i);
print $1,$2,$3,$4
}' OFS='\t' "$bed" > "$output"
```

## Generate bam file from fastq files

Reads should be mapped to the mitochondrial reference sequence. For modern samples, the reference genome should first be indexed:

```bash
bwa-mem2 index ${species}_reordered_mtDNA.fa
samtools faidx ${species}_reordered_mtDNA.fa
picard CreateSequenceDictionary -R $i -O ${species}_reordered_mtDNA.dict
```
Mapping can then be performed using the script map_modern_mem.sh. 

However, museum samples require many additional steps of preprocessing, which are explained elsewhere.

## Generate depth file (optional)

A depth file can also be generated to filter low depth samples with the script create_depth_file.sh.

