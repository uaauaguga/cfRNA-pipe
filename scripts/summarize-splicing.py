#!/usr/bin/env python
import pandas as pd
import numpy as np
import argparse
import sys

import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(name)s] [%(levelname)s] %(message)s')
logger = logging.getLogger("summarize tables")


def se2mat(se):
    lines=se.apply(lambda x:np.array(x.strip().split(",")))#.astype(datatype))
    return np.vstack(lines.values)


def main():
    parser=argparse.ArgumentParser(description='Extract Data From rMATs result')
    parser.add_argument("--indir", "-i", type=str, required=True, help="Input rMATs result table for parsing")
    parser.add_argument("--inclusion", "-inc", type=str, required=True, help="Counts of inclusion / retaintion events")
    parser.add_argument("--exclusion", "-exc", type=str, required=True, help="Counts of exclusion / skipping events")
    parser.add_argument("--method", "-m", type=str, default="JC", help="Count method, either JC or JCEC", choices=["JC","JCEC"])
    parser.add_argument("--sample-ids", "-s", type=str, required=True, help="File contains sample ids")
    args = parser.parse_args()
    
    sample_ids = open(args.sample_ids).read().strip().split("\n")
    id_in_common = ['GeneID', 'geneSymbol', 'chr', 'strand']
    up_down_esee=['upstreamES','upstreamEE', 'downstreamES', 'downstreamEE']
    MXE = id_in_common + ['1stExonStart_0base','1stExonEnd', '2ndExonStart_0base', '2ndExonEnd'] + up_down_esee
    SE = id_in_common + ['exonStart_0base', 'exonEnd'] + up_down_esee
    RI = id_in_common + ['riExonStart_0base','riExonEnd'] + up_down_esee
    A3SS = id_in_common + ['longExonStart_0base','longExonEnd', 'shortES', 'shortEE', 'flankingES', 'flankingEE']
    A5SS = A3SS
    columns_dict={"MXE":MXE,"SE":SE,"RI":RI,"A3SS":A3SS,"A5SS":A5SS}

    logger.info("Load input data ...")     
    inclusions = []
    exclusions = [] 
    for event_type in columns_dict:
        path = f"{args.indir}/{event_type}.MATS.{args.method}.txt"
        logger.info(f"Load {path} ...")
        df = pd.read_csv(path,sep="\t")
        id_columns = columns_dict[event_type]
        junction_ids = event_type + "|" + df.loc[:,id_columns].apply(lambda x:"|".join([str(each) for each in x]),axis=1).values
        inclusion = se2mat(df.loc[:,'IJC_SAMPLE_1'])
        inclusion = pd.DataFrame(index=junction_ids,columns=sample_ids,data=inclusion)
        exclusion = se2mat(df.loc[:,'SJC_SAMPLE_1'])
        exclusion = pd.DataFrame(index=junction_ids,columns=sample_ids,data=exclusion)
        inclusions.append(inclusion)
        exclusions.append(exclusion)
    inclusions = pd.concat(inclusions, axis=0)
    exclusions = pd.concat(exclusions, axis=0)
    logger.info(f"Save inclusion counts to {args.inclusion} ...")
    inclusions.to_csv(args.inclusion, sep="\t")
    logger.info(f"Save inclusion counts to {args.exclusion} ...")
    exclusions.to_csv(args.exclusion, sep="\t")
    logger.info("All done .")


if __name__ == "__main__":
    main()
