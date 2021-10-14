## Installation

- You shall first install anaconda or minconda, then run

```bash
conda env create -f rna-pipe-env.yaml --verbose
```

- If you want to run the statistical part, several R packages are required

```{R}
# R should already been installed in environment rna-pipe-env
# type "conda activate rna-pipe-env" to activate the conda environment
# type "R" to enter R prompt
# run this in R interactive environment
install.packages("argparse")
install.packages("aod")
install.packages("BiocManager")
BiocManager::install(c("limma", "edgeR"))
```

## Overview of used software and packages

- [trim_galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/): for removing adapter and low qulity reads
- [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml): for removing reads derived from unwanted sequences, and circular RNA mapping
- [STAR](https://github.com/alexdobin/STAR): for genome mapping
- [samtools](http://www.htslib.org/)
- [bedtools](https://bedtools.readthedocs.io/en/latest/)
- featureCounts in [subread](http://subread.sourceforge.net/) package
- [rmats](http://rnaseq-mats.sourceforge.net/): for splicing analysis
- [kraken2](https://ccb.jhu.edu/software/kraken2/)
- [snakemake](https://snakemake.readthedocs.io/en/stable/): for pipeline management
- UCSC's bedGraphToBigWig in <http://hgdownload.soe.ucsc.edu/admin/exe/>

- python packages
  - [numpy](https://numpy.org/)
  - [pandas](https://pandas.pydata.org/): manipulate tables
  - [HTSeq](https://htseq.readthedocs.io/en/master/): for handling paired end reads in bam file
  - [tqdm](https://tqdm.github.io/): show a progress bar in for loop
  - [pyBigWig](https://github.com/deeptools/pyBigWig): for parsing bigwig

- R packages
  - [argparse](https://cran.r-project.org/web/packages/argparse/index.html): for parsing arguments
  - [limma](https://bioconductor.org/packages/release/bioc/html/limma.html): for differential analysis of counts
  - [edger](https://bioconductor.org/packages/release/bioc/html/edgeR.html): for differential analysis of counts
  - [aod](https://cran.r-project.org/web/packages/aod/index.html): for beta-binomial glm modeling

