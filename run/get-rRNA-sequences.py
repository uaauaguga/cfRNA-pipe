#!/usr/bin/env python
import argparse
from tqdm import tqdm
import gzip
import pyfaidx

def parseAttr(s):
    info = {}
    s = s.strip()
    for data in s.split(";"):
        data = data.strip()
        key,value = data.split("=")
        info[key] = value
    return info



def main():
    fasta = pyfaidx.Fasta("reference/fasta/hg38.fa")
    fin = gzip.open("reference/source/GRCh38_latest_genomic.gff.gz")
    RefSeq2UCSC = {}
    UCSC2gencode = {}
    with open("reference/source/GRCh38_RefSeq2UCSC.txt") as f:
        for line in f:
            refseq, ucsc = line.strip().split("\t")
            RefSeq2UCSC[refseq] = ucsc
    with open("reference/source/GRCh38_UCSC2gencode.txt") as f:
        for line in f:
            ucsc, gencode = line.strip().split("\t")
            UCSC2gencode[ucsc] = gencode
    f = open("reference/fasta/rRNA.fa","w")
    for line in tqdm(fin):
        line = line.decode()
        if line.startswith("#"):
            continue
        fields = line.strip().split("\t")
        if fields[2] != "rRNA":
            continue 
        attrs = parseAttr(fields[8])
        if "Name" not in attrs or "product" not in attrs:
            continue
        name = attrs["Name"] + " "  + attrs["product"].replace("%"," ")
        chrom, start, end, strand = fields[0], int(fields[3]) - 1, int(fields[4]), fields[6]
        if chrom not in RefSeq2UCSC:
            continue
        if RefSeq2UCSC[chrom] not in UCSC2gencode:
            continue
        chrom = UCSC2gencode[RefSeq2UCSC[chrom]]
        print(f"{chrom}\t{start}\t{end}\t{name}\t.\t{strand}") 
        if strand == "+":
            sequence = str(fasta[chrom][start:end])
        else:
            sequence = str(fasta[chrom][start:end].reverse.complement)
        print(f">{name}",file=f)
        print(sequence,file=f)
    f.close()
        


if __name__ == "__main__":
    main()
