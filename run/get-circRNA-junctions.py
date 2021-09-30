#!/usr/bin/env python
import gzip
from tqdm import tqdm

def load_fasta(path):
    f = gzip.open(path) 
    sequences = {}
    for line in f:
        line = line.decode().strip()
        if line.startswith(">"):
            seq_id = line.replace(">","").split("|")[0] 
            sequences[seq_id] = ""
        else:
            sequences[seq_id] += line
    return sequences


def main():
    anchor_length = 300
    min_length = 100
    print("Load spliced sequences ...")
    sequences = load_fasta("reference/source/human_hg19_circRNAs_putative_spliced_sequence.fa.gz")
    print("Done .")
     
    fout = open("reference/fasta/circRNA.fa","w")
    print("Write junction sequences to reference/fasta/circRNA.fa")
    for seq_id in tqdm(sequences):
        seq = sequences[seq_id]
        if len(seq) < min_length:
            continue
        s = min(len(seq), anchor_length)
        fout.write(f'>{seq_id}\n')
        fout.write(seq[-s:] + seq[:s])
        fout.write('\n')
    fout.close()
    print("All done .")



if __name__ == "__main__":
    main()

