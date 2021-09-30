#!/usr/bin/env python
import numpy as np
import argparse
import pyBigWig
import pandas as pd
from tqdm import tqdm
from multiprocessing import Pooimport loggilogging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(name)s] %(message)s')
logger = logging.getLogger("APA Analysis")

def rle(inarray):
    """ 
    refer to https://stackoverflow.com/questions/1066758/find-length-of-sequences-of-identical-values-in-a-numpy-array-run-length-encoding 
    run length encoding. Partial credit to R rle function. 
    Multi datatype arrays catered for including non Numpy
    returns: (values, run length) 
    """
    ia = np.asarray(inarray)                # force numpy
    n = len(ia)
    if n == 0: 
        return (None, None)
    else:
        y = ia[1:] != ia[:-1]               # pairwise unequal (string safe)
        i = np.append(np.where(y), n - 1)   # must include last element posi
        z = np.diff(np.append(-1, i))       # run lengths
        #p = np.cumsum(np.append(0, z))[:-1] # positions
    return(ia[i], z)


def rld(v, l):
    return np.repeat(v, l).astype(int)
    

def load_coverage(bigwig_path):
    coverages = {}
    bw = pyBigWig.open(bigwig_path)
    i = 0
    for chrom, start, end, utr_id, score, strand in intervals:
        start, end = int(start), int(end)
        if strand == "-":
            start += int((end - start)*0.2)
        else:
            end -= int((end - start)*0.2)
        depths = np.nan_to_num(np.array(bw.values(chrom, start, end)),0).astype(int)
        if strand  == "-":
            depths = depths[::-1]
        # use run length encoding to reduce memory overhead
        coverages[utr_id] = rle(depths)
    libsize = bw.header()['sumData']
    return coverages, libsize
        

def load_bed(path,min_length=300):
    intervals = []
    n_too_short = 0
    n_total = 0
    n_duplicated = 0
    existed = set()
    with open(path) as f:
        for line in f:
            n_total += 1
            fields = line.strip().split("\t")
            chrom, start, end, utr_id, score, strand = fields
            if f"{chrom}|{start}|{end}|{strand}" in existed:
                n_duplicated += 1
                continue
            existed.add(f"{chrom}|{start}|{end}|{strand}")
            start, end = int(start), int(end)
            if end - start < min_length:
                n_too_short += 1
                continue
            intervals.append((chrom, start, end, utr_id, score, strand))
    logger.info(f"Among {n_total} provided 3' UTRs")
    logger.info(f"{n_too_short} are shorter than {min_length}")
    logger.info(f"{n_duplicated} are deuplicated")
    logger.info(f"{n_total - n_too_short - n_duplicated} are used")
    return intervals


def get_mse(tx_coverage,position):
    """
    Input : n by L matrix
    n is number of samples, L is length of current 3'UTR
    """
    shorter_mean = tx_coverage[:,:position].mean(axis=1).reshape(-1,1)
    longer_mean = tx_coverage[:,position:].mean(axis=1).reshape(-1,1)
    shorter_mean = shorter_mean - longer_mean  
    shorter_mean[shorter_mean<0] = 0
    shorter_diff = tx_coverage[:,:position] - shorter_mean - longer_mean
    longer_diff = tx_coverage[:,position:] - longer_mean
    diff = np.concatenate([shorter_diff,longer_diff],axis=1)
    return np.mean(diff**2)


def infer_apa_sites(utr_id,libsizes, min_coverage, min_sample):
    coverages = []
    for i in range(len(coverages_list)):
        v, l = coverages_list[i][utr_id]
        coverage = rld(v,l)
        coverages.append(coverage.reshape(1,-1))
    coverages = np.concatenate(coverages,axis=0)
    if (coverages.mean(axis=1) > min_coverage).sum() < min_sample:
        return None
    coverages_rpm = 1000000*coverages/np.array(libsizes).reshape(-1,1)
    mses = []
    positions = []
    for position in range(200,int(coverages_rpm.shape[1]*0.9)):
        mse = get_mse(coverages_rpm,position)
        mses.append(mse)
        positions.append(position)
    min_position = positions[np.argmin(mses)]
    shorter_abundance = coverages[:,:min_position].mean(axis=1)
    longer_abundance = coverages[:,min_position:].mean(axis=1)
    shorter_abundance = shorter_abundance - longer_abundance
    shorter_abundance[shorter_abundance<0] = 0
    PDUI = np.empty(shorter_abundance.shape[0])
    PDUI.fill(np.nan)
    abundance = shorter_abundance+longer_abundance
    zero_mask = abundance==0
    PDUI[~zero_mask] = longer_abundance[~zero_mask]/abundance[~zero_mask]
    chrom, tx_id, start, end, strand = utr_id.split("|")
    if strand == "-":
        infered_site = int(end) - min_position - 1
    else:
        infered_site = int(start) + min_position
    return utr_id, PDUI, infered_site

    
def main():
    parser = argparse.ArgumentParser(description='Infer poly A sites and isoform abundance')
    parser.add_argument('--bed','-b',type=str,required=True,help="3' UTR annotations")
    parser.add_argument('--bigwig-paths','-bws',type=str,help="A text file contains path of bigwig file to process",required=True)
    parser.add_argument('--min-length','-ml',type=int,help="Min length of 3' UTR to consider",default=300)
    parser.add_argument('--min-coverage','-mc',type=int,help="Min coverage required for 3' UTR",default=10)
    parser.add_argument('--min-sample','-ms',type=int,help="To keep a transcript, More than this number should meet the min coverage",default=3)
    parser.add_argument('--njobs','-j',type=int,help="Number of process to run",default=1)
    parser.add_argument('--PDUI',help="Output PDUI matrix",required=True)
    args = parser.parse_args()   

    global coverages_list
    global intervals

    coverages_list = []
    libsizes = []
    sample_ids = []

    logger.info("Load bed files ...")
    intervals = load_bed(args.bed,args.min_length)

    paths = open(args.bigwig_paths).read().strip().split("\n")
    pool = Pool(args.njobs)
    workers = []
    logger.info(f"Load coverages with {args.njobs} processes, it would take a while ...")
    for path in paths:
        sample_id = ".".join(path.split("/")[-1].split(".")[:-1])
        sample_ids.append(sample_id)
        workers.append(pool.apply_async(func=load_coverage,args=(path,)))
    for worker in workers:
        coverages, libsize = worker.get()
        coverages_list.append(coverages)
        libsizes.append(libsize)
        logger.info(f"{path} loaded .")

    # Analysis for different transcripts are independent
    logger.info(f"Infer APA sites for 3' UTRs with {args.njobs} processes...")
    pool = Pool(args.njobs)
    workers = [] 
    for utr_id in coverages_list[0].keys():
        workers.append(pool.apply_async(func=infer_apa_sites,args=(utr_id,libsizes, args.min_coverage, args.min_sample)))
    results = []
    for worker in tqdm(workers):
        result = worker.get()
        if result is None:
            continue
        results.append(result)
    logger.info(f"{len(results)} in {len(workers)} UTRs passed the coverage filter .")
    PDUIs = []
    utr_ids_passed_filter = []
    infered_sites = []

    logger.info("Save PDUIs and infered ployA sites ...")
    for utr_id, PDUI, infered_site in results:
        PDUIs.append(PDUI.reshape(1,-1))
        infered_sites.append(infered_site)
        utr_ids_passed_filter.append(utr_id)
    PDUIs = np.concatenate(PDUIs,axis=0)
    PDUIs = np.round(PDUIs,3)
    PDUIs = pd.DataFrame(data=PDUIs,index=utr_ids_passed_filter,columns=sample_ids)    
    PDUIs["infered-polyA-sites"] = np.array(infered_sites).reshape(-1,1)
    PDUIs.to_csv(args.PDUI, sep="\t")         
    logger.info("All done .")
if __name__ == "__main__":
    main()
