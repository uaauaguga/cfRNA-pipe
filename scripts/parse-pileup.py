#!/usr/bin/env python
import re
import argparse
def main():
    parser = argparse.ArgumentParser(description="Extract Editing Info From pileup results")
    parser.add_argument("--input","-i",help="input mpileup",required=True)
    parser.add_argument("--output","-o",help="output",required=True)
    args = parser.parse_args()

    editing = {"A":"G","T":"C"}

    fout = open(args.output,"w")
    print("editing_id","chrom", "pos", "ref", "n_ref", "n_edit", "n_skip", "n_others", sep="\t", file=fout)
    with open(args.input) as fin:
        for line in fin:
            chrom, pos, ref, depth, match, qual = line.strip().split("\t")
            ref = ref.upper()
            if ref not in editing:
                continue
            if "+" in match or "-" in match:
                # remove insertions and deletions in match field
                match2 = ""
                for i, m in enumerate(re.split(r"[+-]",match)):
                    if i == 0:
                        match2 = m
                    else:
                        offset = re.match(r"\d+",m).group(0)
                        m = m[len(offset)+int(offset):]
                        match2 += m
                match = match2
            match = match.upper()
            n_ref = match.count(".") + match.count(",")
            n_skip = match.count("<") + match.count(">")
            n_edit = match.count(editing[ref])
            if n_edit == 0:
                continue
            print(f"{chrom}|{pos}|{ref}|{editing[ref]}",chrom, pos, ref, n_ref, n_edit, n_skip, int(depth) - n_ref - n_edit - n_skip, sep="\t", file=fout)
    fout.close()


if __name__ == "__main__":
    main()

