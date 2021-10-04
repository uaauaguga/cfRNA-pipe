## Overview of the requirements

- The versions software tested in this pipeline are shown in the bracket
- trim_galore: for removing adapter and low qulity reads
- bowtie2: for removing reads derived from unwanted sequences, and circular RNA mapping
- STAR: for genome mapping
- samtools
- bedtools
- featureCounts
- rmats
- kraken2
- snakemake: for pipeline management
- bedGraphToBigWig
  ```bash
  conda install -c bioconda ucsc-bedgraphtobigwig
   ```

- python packages
  - numpy
  - pandas
  - HTSeq
  - tqdm: show a progress bar in for loop
  - pyBigWig: for parsing bigwig

- R packages
  - limma
  - edger
  - aod
  - argparse


## Installation


