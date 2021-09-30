#!/usr/bin/env python
import os
import sys
import argparse

def extraction_3putr(tx_bed_path,utr3p_bed_path):    
    fout = open(utr3p_bed_path,'w')
    scanned_3UTRs = set()
    num_saved = 0
    for line in open(tx_bed_path,'r'):
        fields = line.strip('\n').split('\t')
        chrom, start, end, tx_id, strand = fields[0], int(fields[1]), int(fields[2]), fields[3], fields[5]
        if strand == "+":
            utr_end = end
            # use the last exon in transcript as model 3p utr
            # in bed12 format, field 11 is blocksize, field 12 is start position relative to tx start position
            utr_start = start+int(fields[11].strip(',').split(',')[-1])
        elif strand == "-":
            utr_start = start
            utr_end   = start + int(fields[10].split(',')[0]) 
        utr_id = "|".join([chrom,tx_id,str(utr_start),str(utr_end),strand])          
        if utr_id not in scanned_3UTRs:
            record = [chrom, str(utr_start), str(utr_end), utr_id, '0', strand]
            fout.writelines('\t'.join(record) + '\n')
            scanned_3UTRs.add(utr_id)
            num_saved += 1
    fout.close()   
    print("Total extracted 3' UTR: " + str(num_saved))

def main():
    parser = argparse.ArgumentParser(description="Extract 3'UTR annotation from gene model in bed12 format")
    parser.add_argument('--bed','-b',type=str,required=True,help="Input gene model in bed12 format, one line corresponds to a transcript")
    parser.add_argument('--utr','-u',type=str,help="3' UTR annotation in bed format",required=True)
    args = parser.parse_args()
    print("Extract 3' UTRs ...")
    extraction_3putr(args.bed,args.utr)
    print("Finished.")

if __name__ == '__main__':
    main()
    
