#!/bin/bash
for s in spikein circRNA rRNA UniVec;do
  bowtie2-build reference/fasta/${s}.fa reference/bowtie2-index/${s} 
done

bowtie2-build reference/fasta/spikein.fa,reference/fasta/rRNA.fa,reference/fasta/UniVec.fa reference/bowtie2-index/unwanted-with-spikein
bowtie2-build reference/fasta/rRNA.fa,reference/fasta/UniVec.fa reference/bowtie2-index/unwanted
