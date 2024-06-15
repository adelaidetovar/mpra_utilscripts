# Author: Jacob Kitzman, kitzmanj

import sys
import argparse
from collections import OrderedDict

import scipy.misc as sm

import pandas as pd

if __name__=='__main__':

    # process and filter consensus barcode filter list from 

    opts = argparse.ArgumentParser()
   
    opts.add_argument('--in_histo',dest='in_histo')
    opts.add_argument('--in_has_header',default=False,action='store_true',dest='in_has_header')
    opts.add_argument('--fxn_pass', default="lambda l:True", dest='fxn_pass')
    opts.add_argument('--out_histo',dest='out_histo')
    opts.add_argument('--out_summary',dest='out_summary')

    o = opts.parse_args()

    fxn_pass = eval(o.fxn_pass)

    if not o.in_has_header:
        tbl_in = pd.read_csv( o.in_histo, header=None, sep='\t' )
    else:
        tbl_in = pd.read_csv( o.in_histo, sep='\t' )

    if not o.in_has_header:
        tbl_in.columns = [ 'seq', 'count', 'tags' ]

    tbl_in['filtered'] = [ 'pass' if fxn_pass( l ) else 'FAIL' for _,l in tbl_in.iterrows() ]

    tbl_in = tbl_in[ ['seq','filtered','count','tags'] ]

    nItems = tbl_in.shape[0]
    nItemsPass = (tbl_in.filtered == 'pass').sum()
    totalCounts = tbl_in['count'].sum()
    countsPass = (tbl_in.loc[tbl_in.filtered == 'pass', 'count']).sum()

    outSummary = pd.DataFrame(
        OrderedDict(
        [ ('histo',[o.in_histo]),
           ('n_bcs',[nItems]),
           ('n_bcs_pass',[nItemsPass]),
           ('frac_bcs_pass',[nItemsPass/float(nItems) if nItemsPass>0 else 0]),
           ('n_reads',[totalCounts]),
           ('n_reads_pass',[countsPass]),
           ('frac_reads_pass',[countsPass/float(totalCounts) if totalCounts>0 else 0])
        ]))

    outSummary.to_csv(o.out_summary, sep='\t', index=False)

    tbl_in.to_csv( o.out_histo, sep='\t', index=False )
