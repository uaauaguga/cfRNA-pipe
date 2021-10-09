## Usage

### Quantification part

- This part takes paired end fastq files as input, and produce specified in configure file
- The workflow is wrapped in a single Snakefile 
- A configure file fully control the behavior of this workflow
- In the configure file, the following parameter is required
  - `input_dir`: directory contains raw paired end fastq files. These files should be gzip compressed, and named as {sample_id}_1.fastq.gz and {sample_id}_2.fastq.gz
  - `sample_ids`: path of a text file contains ids of samples to run. One sample id per line
  - `output_dir`: directory to save outputs

- Other available parameters have a defult setting in `config/default.yaml`
  - The default settings will be overwriten if the parameters are setted in configure file
  - Some parameters control which steps to run
  - Some parameters are controls behaviors of programs to run
  - See comments in `config/default.yaml` for detail

```bash
snakemake --jobs 100 --configfile your.configuration.yaml
```

- If you want to known which steps will be run instead of actually running these steps, you can use the "dry run" option of snakemake

```bash
# -np is short for --dry-run (-n) and --printshellcmds (-p) 
snakemake --jobs 100 --configfile your.configuration.yaml -np
```

### Statistical analysis part

- This part provides two scripts for differential analysis
  - `scripts/differential-expression-analysis.R` is used for modeling counts data
  ```bash
  scripts/differential-expression-analysis.R --matrix count.matrix.txt --label-field label --covariate-fields batch --normalize TMM --output diff.table.txt --metadata metadata.txt --case-label T --control-label N --test edger-glmlrt
  ```
  - `scripts/differential-proportion-analysis.R` is used for modeling relative abundance of two counts
  ```bash
  scripts/differential-proportion-analysis.R --matrix-1 counts_1.txt --matrix-2 counts_2.txt --metadata metadata.txt --label-field label --covariate-fields batch --case-label T --control-label N --output diff.table.txt --cores 8
  ```

