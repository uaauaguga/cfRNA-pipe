#!/bin/bash
snakemake --jobs 100 --configfile config/test_pe.yaml --cluster "bsub -R span[hosts=1] -q TEST-1U -n {threads}" --latency-wait 100 --keep-going
