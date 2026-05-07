# Visualise circular genomes

When deciding which gene to use as a starting location, it may be valuable to visualise the mitochondrial genome. However, many existing tools require a genbank file to do so.

A separate conda environment is provided for visualisation, so install and activate this with

```bash
conda env create -f Visualisation/conda env create -f Visualisation/Visualisation_env.yml
conda activate Chloroplot_R
```

## Create genbank file for visualisation

If the annotation is already provided on ncbi, skip this step. Otherwise, the SEQ file produced by mitos2 can be converted to genbank file for visualisation. Download this file and rename it to ${species}_annotation.tbl

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
./bonus_scripts/update_genbank.awk $genus $species ${species}_mtDNA.gb > ${species}_mtDNA_updated.gb
```

## Visualisation

The genbank file can be visualised in R with chloroplot, which will come installed in the conda environment:

```bash
R
library(chloroplot)
table <- PlotTab("species_mtDNA_updated.gb", T)
PlotMitGenome(table, organelle_type = F)
q()
```

Alterantively the [Chloroplot](https://irscope.shinyapps.io/Chloroplot/) or [OGDraw](https://chlorobox.mpimp-golm.mpg.de/OGDraw.html) webservers can be used instead.

Here is an example chloroplot output with the mitochondrial sequence of the butterfly *Papillio machaon*

<img width="1016" height="855" alt="image" src="https://github.com/user-attachments/assets/53243c0a-ffc5-4b74-921e-263ba7d35b8a" />

There is a large gap on the right that corresponds to a large non-coding region (D-loop). Mitochondrial genomes are often ordered such that the D-loop sequence is at the end, so the tRNA-M gene appears to be a good starting point.
