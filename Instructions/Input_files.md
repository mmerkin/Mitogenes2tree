# Create input files for Mitogenes2tree

Mitogenes2tree requires the following three inputs:
- [Mitochondrial genome reference fasta file](https://github.com/mmerkin/Mitogenes2tree/blob/main/Instructions/Input_files.md#generate-relinearised-mitochondrial-reference-sequence)
- [Bed file with locations of each gene](https://github.com/mmerkin/Mitogenes2tree/blob/main/Instructions/Input_files.md#generate-bed-file)
- [Bam file of sequencing reads mapped to the reference](https://github.com/mmerkin/Mitogenes2tree/blob/main/Instructions/Input_files.md#generate-bam-file-from-fastq-files)
- [A list of sample names (required) and depths (optional)](https://github.com/mmerkin/Mitogenes2tree/blob/main/Instructions/Input_files.md#generate-sample-file)


## Generate relinearised mitochondrial reference sequence

If a reference genome currently exists for the species of interest or a close relative, the mitochondrial sequence of the fasta file should be downloaded. Alternatively, a program such as [getOrganelle](https://github.com/Kinggerm/GetOrganelle) or [NOVOPlasty](https://github.com/ndierckx/NOVOPlasty) can be used to generate a reference sequence from short read data.

Since mitochondria possess circular genomes, the start point of this sequence will be random and could even be in the middle of a gene. As such, it is recommended to first relinearise the reference genome to start at a more appropriate position. Many taxa also have standard start regions for mitogenomes, so relinearisation may aid with literature comparisons.


### Annotation with mitos2

If an annotation already exists on ncbi (identified by the page listing the locations of each gene), the file should be downloaded in genbank format.

Otherwise, the reference fasta sequence can be annotated using mitos2 on the [galaxy webserver](https://usegalaxy.eu/). 

Create a galaxy account, upload the fasta file and run the tool MITOS2 to generate SEQ and nucleotide FASTA files, ensuring to also change the genetic code if necessary:

<img width="1173" height="381" alt="image" src="https://github.com/user-attachments/assets/4984a625-002c-4d09-a8db-76fb3428e00d" />

The SEQ file can then be used for visualisation if desired. Instructions to do so are provided [here](https://github.com/mmerkin/Mitogenes2tree/blob/main/Visualisation/Visualisation_env.yml)

### Relinearise sequence

Once a new starting gene has been chosen, a fasta file should be created using the gene sequence from the mitos2 nucleotide fasta output. Then, the script reorder_circular_genome.sh can be run.

```bash
bash bonus_scripts/reorder_circular_gemome.sh ${species}_mtDNA.fa ${gene}.fa > ${output}.fa
```

## Generate bed file

Mitos2 on the [galaxy webserver](https://usegalaxy.eu/) can be used to make a bed file of the gene locations. If the reference genome was reordered, mitos2 needs to be run again on the new reference genome as the positions of each gene will have changed.

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
```

Mapping can then be performed using the script Map_reads_parallel.sh.

The reads must all be found in the same folder and contain R1 and R2 in the forward and reverse read names.

```bash
bash bonus_scripts/Map_reads_parallel.sh ${species}_reordered_mtDNA.fa /path/to/reads /path/to/output/folder
```

### GATK

Optionally, GATK can be used to realign around indels. To enable this, create a new conda environment containing gatk38 and replace the variables at the top of the script. A dictionary will also need to be created with picard
```bash
picard CreateSequenceDictionary -R $i -O ${species}_reordered_mtDNA.dict
```

### Museum samples

Museum samples also require many additional steps of preprocessing. which are explained [here](https://github.com/mmerkin/Mitogenes2tree/blob/main/Instructions/Museum_samples.md))

## Generate sample file

A basic sample file can be generated with ls:

```bash
ls /path/to/bams/*.bam > samples.txt
```

However if you wish to filter samples by average read depth, a tsv should be supplied instead containing depths in the second column. Such a file can be generated with the script create_depth_file.sh.

```bash
bash bonus_scripts/create_depth_file.sh $bedfile $output /path/to/bams/*.bam
```


