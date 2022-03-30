## Usage

- Enter the conda environment `rna-pipe-env` before running the pipeline

```bash
conda activate rna-pipe-env
```

### Quantification part

- This part takes paired end fastq files as input, and produce outputs specified in the configure file

- The workflow is wrapped in a single Snakefile 

- A configure file fully control the behavior of this workflow

- In the configure file, the following parameters are required
  - `input_dir`: directory contains raw paired end fastq files. These files should be gzip compressed, and named as `{sample_id}_1.fastq.gz` and `{sample_id}_2.fastq.gz`
  - `sample_ids`: path of a text file contains ids of samples to run. One sample id per line
  - `output_dir`: directory to save outputs

- Other available parameters have a defult setting in `config/default.yaml`
  - The default settings will be overwriten if the parameters are set in configure file
  - Some parameters control which steps to run, some control behaviors of the programs to run
  - See comments in `config/default.yaml` for detail

```bash
snakemake --jobs 8 --configfile config/test_pe.yaml
```

- If you want to know which steps will be run, instead of actually running these steps, you can use the "dry run" option of snakemake

```bash
# -np is short for --dry-run (-n) and --printshellcmds (-p) 
snakemake --jobs 8 --configfile config/test_pe.yaml -np
```

- If you use a LSF cluster for parallel computation

```bash
# use bsub to submit job 
# --latency-wait 100: allow 100s of file system delay
snakemake --jobs 32 --configfile config/test_pe.yaml --cluster "bsub -R span[hosts=1] -q queue_name -n {threads}" --latency-wait 100 
```

### Statistical analysis part

- This part provides two scripts for differential analysis

- `scripts/differential-expression-analysis.R` is used for modeling counts data

```bash
  scripts/differential-expression-analysis.R --matrix count.matrix.txt --label-field label [--covariate-fields batch] --normalize TMM --output diff.table.txt --metadata metadata.txt --case-label T --control-label N --test edger-glmlrt
```

- `scripts/differential-proportion-analysis.R` is used for modeling relative abundance of two counts, each for abundance of one isoform

```bash
  scripts/differential-proportion-analysis.R --matrix-1 counts_1.txt --matrix-2 counts_2.txt --metadata metadata.txt --label-field label [--covariate-fields batch] --case-label T --control-label N --output diff.table.txt --cores 8
```


### Visualization part

- PCA & MDS plot

- Volcano plot

- For visualization of RNA editing and splicing events, you shall load bam files (sorted by genome coordinate, present in `{output_dir}/{sample_id}/genome.sorted.bam` ) with IGV. 

- [rmats2sashimiplot](https://github.com/Xinglab/rmats2sashimiplot) is also a great tool for visualze alternative splicing events. It can be simply installed by typing `conda install -c bioconda rmats2sashimiplot`. See their documentation at <https://github.com/Xinglab/rmats2sashimiplot> for detailed usage.

- For visualization coverage pattern at 3' UTR, you can simply load the `bigwig` file (present in `{output_dir}/bigwig`) with [IGV](https://software.broadinstitute.org/software/igv/)


- We recommand [Krona](https://github.com/marbl/Krona) for visualze taxonomy composition of reads unmapped to human genome. 
  - Installation and preparing taxo database
```bash
conda install -c bioconda krona # install krona into current conda environment
ktUpdateTaxonomy.sh # prepare taxonomy database. this script is shipped with krona
```

  - Import kraken2 report to (files in `{output_dir}/microbe/report` ) Krona's html. You can open this file with any web browser for further exploration 

```bash
ktImportTaxonomy -m 3 -t 5 {output_dir}/microbe/report/{sample_id}.txt -o krona.html
```

- Visualize the taxonomy composition is useful under the following circumstances

  - You have a large unmapped rate, but don't know what's going wrong. Many potential reasons could lead to such problems, and a krona plot is helpful for debugging. 
    - For low-input RNA sequencing, (cell free RNA for example), microbial contaminations are prevalent. If in krona plot, the majority of unmapped reads were assigned to a few microbial species, this suggests a considerable fraction of in your RNA-seq reads were derived from microbe contaminations.
    - If the majority of unmapped reads were assigned to human genome, the large unmapped rate may be attributed to poor data quality (kraken2 is less stringent than STAR aligner), failed adapter trimming, or out-of-paired fastq file (reads id in `{sample_id}_1.fastq.gz` and `{sample_id}_2.fastq.gz` should have same order), etc.

  - Sometimes you may interest in microbial reads in your sequencing data. For example, several studies suggest some tumors contain clinically relevant living bacteria.


