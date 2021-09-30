#!/bin/bash
STAR --runMode genomeGenerate --genomeDir reference/star-index/hg38 --genomeFastaFiles reference/fasta/hg38.fa --sjdbGTFfile reference/gtf/gencode.v38.annotation.gtf --runThreadN 10
