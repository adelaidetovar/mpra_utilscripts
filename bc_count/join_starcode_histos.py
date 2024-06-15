# Author: Jacob Kitzman, kitzmanj

import sys
import socket
import os
import os.path
import argparse
import re

import pysam

import numpy as np

#from jkutils import *

import pandas as pd

import Bio.Seq

if __name__=='__main__':
    
    opts = argparse.ArgumentParser()

    opts.add_argument('--inSampleKey',dest='inSampleKey')
    opts.add_argument('--colSampname',default=None,dest='colSampname')
    opts.add_argument('--colHistoPath',default=None,dest='colHistoPath')

    opts.add_argument('--fxnInputPassFilter',default=None,dest='fxnInputPassFilter')

    opts.add_argument('--fxnTransformBarcodeSeq',default=None,dest='fxnTransformBarcodeSeq')

    opts.add_argument('--dropNbarcs',default=False, action='store_true',dest='dropNbarcs')

    opts.add_argument('--columnNormalize',default=False, action='store_true',dest='columnNormalize')

    opts.add_argument('--outHisto',dest='outHisto')

    o = opts.parse_args()

    if o.fxnTransformBarcodeSeq is not None:
        fnbc=eval(o.fxnTransformBarcodeSeq)
    else:
        fnbc=None

    sampkey = pd.read_csv( o.inSampleKey, sep='\t' )

    assert( sampkey.shape[0] == len(sampkey[o.colSampname].unique()) ), 'should not be any repeated samples.'
    
    sBcs = set()

    hist_joint = pd.DataFrame( {'barcode':[] } )

    if o.fxnInputPassFilter is not None:
        fxnInputPassFilter=eval(o.fxnInputPassFilter)
    else:
        fxnInputPassFilter=None

    rctr=0

    for _, rsamp in sampkey.iterrows():
        hist_cur = pd.read_csv( rsamp[o.colHistoPath], sep='\t', header=None )
        samp_name = rsamp[o.colSampname]

        print('processing %d/%d'%(rctr,sampkey.shape[0]),samp_name,' ',rsamp[o.colHistoPath])
        rctr+=1

        hist_cur.columns = ['barcode', samp_name, 'consituent_barcodes' ]

        # discard the consituent barcodes for now
        hist_cur = hist_cur[ ['barcode',samp_name] ]

        oldN,oldCt=hist_cur.shape[0],hist_cur[samp_name].sum()
        print('   %d lines, %d counts'%(oldN,oldCt))

        # discard N-containing barcodes
        if o.dropNbarcs:
            hist_cur.drop( hist_cur.index[hist_cur['barcode'].str.upper().str.contains('N')], inplace=True )
            filtN,filtCt=hist_cur.shape[0],hist_cur[samp_name].sum()
            print('  dropN --> %d lines (%.2f%%), %d counts (%.2f%%)'%( filtN,100.*filtN/float(oldN) if oldN>0 else 0, filtCt, 100.*filtCt/float(oldCt) if oldCt>0 else 0 ))
            oldN,oldCt=filtN,filtCt

        if fxnInputPassFilter:
            hist_cur=hist_cur.loc[ fxnInputPassFilter(hist_cur) ].copy()
            filtN,filtCt=hist_cur.shape[0],hist_cur[samp_name].sum()
            print('  filt  --> %d lines (%.2f%%), %d counts (%.2f%%)'%( filtN,100.*filtN/float(oldN) if oldN>0 else 0, filtCt, 100.*filtCt/float(oldCt) if oldCt>0 else 0 ))

        hist_cur = hist_cur.sort_values(by='barcode')
        
        if hist_cur.shape[0]==0: continue

        hist_joint = pd.merge( hist_joint, hist_cur, how='outer', on='barcode' )
        hist_joint = hist_joint.sort_values(by='barcode')
        print('joined -> %d lines'%( hist_joint.shape[0] ))

        sys.stdout.flush()

    mean_counts = hist_joint[hist_joint.columns[1:]].sum(axis=1)
    hist_joint = hist_joint.iloc[ np.argsort(-mean_counts) ]

    for samp_name in sampkey[o.colSampname]:
        hist_joint.loc[ pd.isnull(hist_joint[samp_name]), samp_name ] = 0.

    if o.columnNormalize:
        for samp_name in sampkey[o.colSampname]:
            if hist_joint[ samp_name ].sum() > 0:
                hist_joint[ samp_name ] = hist_joint[ samp_name ].astype('f')/hist_joint[ samp_name ].sum()


    if fnbc is not None:
        hist_joint['barcode'] = [ fnbc(x) for x in hist_joint['barcode'] ]

    hist_joint.to_csv( o.outHisto, sep='\t', index=False )
