## Installation

- You shall first install anaconda or minconda, then run

```bash
conda env create -f rna-pipe-env.yaml --verbose
```

- If you want to run the statistical part, several R packages are required

```{R}
# R should already been installed in environment rna-pipe-env
# run this in R environment
install.packages("argparse")
install.packages("aod")
install.packages("BiocManager")
BiocManager::install(c("limma", "edgeR"))
```

## Overview of used software and packages

- trim_galore: for removing adapter and low qulity reads
- bowtie2: for removing reads derived from unwanted sequences, and circular RNA mapping
- STAR: for genome mapping
- samtools
- bedtools
- featureCounts in subread
- rmats
- kraken2
- snakemake: for pipeline management
- bedGraphToBigWig

- python packages
  - numpy
  - pandas: manipulate tables
  - HTSeq: for handling paired end reads in bam file
  - tqdm: show a progress bar in for loop
  - pyBigWig: for parsing bigwig

- R packages
  - argparse: for parsing arguments
  - limma: for differential analysis of counts
  - edger: for differential analysis of counts
  - aod: for beta-binomial glm modeling

