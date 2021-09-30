#!/usr/bin/env python
import argparse
from collections import OrderedDict
def main():
    parser = argparse.ArgumentParser(description="Select a single UTR region fro each transcript")
    parser.add_argument('--input','-i',type=str,help="Input 3' UTR annotation",required=True)
    parser.add_argument('--output','-o',type=str,help="Output 3' utr annotation",required=True)
    args = parser.parse_args()
    utr_dict = OrderedDict()
    with open(args.input) as f:
        for line in f:
            fields = line.strip().split("\t")
            chrom, start, end, utr_id, score, strand = fields
            if utr_id not in utr_dict: 
                utr_dict[utr_id] = (chrom, int(start), int(end), utr_id, score, strand)
            else:
                start, end = int(start), int(end) 
                if (strand == "+" and start < utr_dict[utr_id][1]) or (strand == "-" and end > utr_dict[utr_id][2]):
                    utr_dict[utr_id] = (chrom, start, end, utr_id, score, strand)
                else:
                    continue
    with open(args.output,"w") as fout:
        for utr_id in utr_dict:
            chrom, start, end, utr_id, score, strand = utr_dict[utr_id]
            start, end = str(start), str(end)
            fout.write("\t".join([chrom, start, end, utr_id, score, strand])+"\n")
        
                    
if __name__ == "__main__":
    main()
