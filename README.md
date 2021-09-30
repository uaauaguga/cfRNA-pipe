## RNA-pipe
### Introduction
- This repo implement a simple, efficient, and highly modular RNA-seq analysis pipeline
- It was a general purposed pipeling, although it was originally designed for analyzing cell free RNA in human plasma, it could be migrated to analysis any paired-end long RNA sequencing data
- Take advantage of snakemake to automatically solve for inter dependency between intermediate files
- We provide the following features
  - Gene and circRNA expression analysis, using the traditional mapping and counting schema
  - Analysis of A to I editing events, count editing events at kown RNA editing sites
  - Coverage pattern at 3' UTR: a much more efficient implementation for method descript in [dapar](https://github.com/ZhengXia/dapars)
  - Alternative splicing analysis: a wrapper for [rmats](http://rnaseq-mats.sourceforge.net/), provide a more user friendly interface
  - Analysis of unmapped reads: perform metagenomic classification with [kraken2](https://ccb.jhu.edu/software/kraken2/)
  - All of functionalities have separated a quantification part (generate matrix for different features) and a statistical analysis part (multiple statistical testing for difference between cases / controls or different treatment, allow adjust for covariates and batch effects)
  - Several useful scripts for data visualization
