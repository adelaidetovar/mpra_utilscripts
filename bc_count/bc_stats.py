import argparse
import pandas as pd
import numpy as np
import os
import glob

def main():
    opts = argparse.ArgumentParser()

    opts.add_argument("--in_dir", required=True)
    opts.add_argument("--out_tab", required=True)

    o = opts.parse_args()

    files = glob.glob(os.path.join(o.in_dir, "*.bc_cluster.txt"))
    out_tab = o.out_tab

    if not files:
        print("Check that you've specified the correct directory for barcode clusters and that all files are present.")
        return

    # initialize
    topthresh = np.zeros(len(files))
    midthresh = np.zeros(len(files))
    lowthresh = np.zeros(len(files))
    maxbc = np.zeros(len(files))
    topfrac = np.zeros(len(files))
    midfrac = np.zeros(len(files))
    lowfrac = np.zeros(len(files))

    # process for each sample
    for i, file in enumerate(files):
        t = pd.read_csv(file, sep='\t', header=None, names=['barcode', 'counts', 'cluster'], na_values=['NA'], dtype={'barcode': str, 'counts': float, 'cluster': str})
        t = t.sort_values(by='counts', ascending=False).dropna(subset=['counts'])
        cfvs = np.cumsum(t['counts']) / np.sum(t['counts'])

        topthresh[i] = np.min(np.where(cfvs >= 0.95))
        midthresh[i] = np.min(np.where(cfvs >= 0.9))
        lowthresh[i] = np.min(np.where(cfvs >= 0.8))
        maxbc[i] = len(cfvs)
        topfrac[i] = (topthresh[i] / maxbc[i])*100
        midfrac[i] = (midthresh[i] / maxbc[i])*100
        lowfrac[i] = (lowthresh[i] / maxbc[i])*100

    # make results table
    results = pd.DataFrame({
        'Total_number_bc': maxbc,
        'BC_at_80pct_of_reads': lowthresh,
        'BC_at_90pct_of_reads': midthresh,
        'BC_at_95pct_of_reads': topthresh,
        'Pct_of_BC_at_80pct_of_reads': lowfrac,
        'Pct_of_BC_at_90pct_of_reads': midfrac,
        'Pct_of_BC_at_95pct_of_reads': topfrac
    })

    # get sample names from input
    results['Sample'] = results['files'].apply(lambda x: x.split('/')[-1].replace('.bc_cluster.txt', ''))

    # write to file
    results.to_csv(o.out_tab, sep='\t', index=False)


if __name__ == '__main__':
    main()
