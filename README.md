## RNA-pipe

### Introduction
- This repo implements a simple, efficient, and highly modular RNA-seq analysis pipeline
- It was a general purposed pipeling. It was originally designed for analyzing cell free RNA in human plasma, but in principle it could be applied to analysis any paired-end long RNA sequencing data
- It takes advantage of workflow management tool [snakemake](https://snakemake.readthedocs.io/en/stable/) to automatically solve for inter-dependency between intermediate files, and execute analysis steps in parallel

- We provide the following features

  - Quantification of gene expression profile and several post transcriptional regulation events
    - Gene and circRNA expression analysis: using the traditional mapping and counting schema
    - Analysis of A to I editing events: count editing events at known RNA editing sites efficiently
    - Coverage pattern at 3' UTR: a much more faster implementation for the method descript in [dapar](https://github.com/ZhengXia/dapars)
    - Alternative splicing analysis: a wrapper for [rmats](http://rnaseq-mats.sourceforge.net/), provide a interface that is more friendly to user and parallel computing environment 
    - Analysis of unmapped reads: perform metagenomic classification with [kraken2](https://ccb.jhu.edu/software/kraken2/)

  - Statistical analysis
    - The quantification part generate two types of data 
      - Counts data (for gene, circRNA expression profile and microbial abundance): linear model in edger and limma
      - Relative proportion of counts data for two isoform (for splicing, 3' UTR coverage pattern and RNA editing analysis): statistcal testing based on beta-binomial generalized linear model (GLM), quasi-binomial GLM, or bionomial GLM.
    - Allow adjust for covariates and batch effects. 

  - Several useful scripts for data visualization
