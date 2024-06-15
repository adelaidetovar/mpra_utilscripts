# Author: Jacob Kitzman, kitzmanj

import sys
import socket
import os
import os.path
from optparse import OptionParser
from collections import defaultdict
import tempfile
import time

import re

import pysam

import numpy as np

import pandas as pd

#from jkutils import *

import scipy.stats as stats

def dotProd( jh, lsamps, suffix='_sum_hits' ):
    Nsamps=len(lsamps)
    mtx = jh[ [ '%s%s'%(s,suffix) for s in lsamps ] ]
    mtx.columns = [ s for s in lsamps ]
    mtxfrac = mtx.copy()
    for s in mtxfrac:
        mtxfrac[s] = np.array( mtxfrac[s], 'f') / sum(mtx[s])
    mtxout = np.zeros( (Nsamps,Nsamps), dtype=np.float32 )
    for isamp in range(Nsamps):
        for jsamp in range(Nsamps):
            mtxout[ isamp, jsamp ] = np.sum(mtxfrac[ lsamps[isamp] ] * mtxfrac[ lsamps[jsamp] ])
    tblout = pd.DataFrame(mtxout, columns=lsamps)
    tblout['sample']=lsamps
    tblout=tblout[ ['sample']+lsamps ]
    return tblout

def spearmanr( jh, lsamps, suffix='_sum_hits' ):
    Nsamps=len(lsamps)
    mtx = jh[ [ '%s%s'%(s,suffix) for s in lsamps ] ]
    mtx.columns = [ s for s in lsamps ]
    mtxfrac = mtx.copy()
    for s in mtxfrac:
        mtxfrac[s] = np.array( mtxfrac[s], 'f') / sum(mtx[s])
    mtxout = np.zeros( (Nsamps,Nsamps), dtype=np.float32 )
    for isamp in range(Nsamps):
        for jsamp in range(Nsamps):
            mtxout[ isamp, jsamp ], _ = stats.spearmanr( mtxfrac[ lsamps[isamp] ], mtxfrac[ lsamps[jsamp] ] )
    tblout = pd.DataFrame(mtxout, columns=lsamps)
    tblout['sample']=lsamps
    tblout=tblout[ ['sample']+lsamps ]
    return tblout

def pearsonr( jh, lsamps, suffix='_sum_hits' ):
    Nsamps=len(lsamps)
    mtx = jh[ [ '%s%s'%(s,suffix) for s in lsamps ] ]
    mtx.columns = [ s for s in lsamps ]
    mtxfrac = mtx.copy()
    for s in mtxfrac:
        mtxfrac[s] = np.array( mtxfrac[s], 'f') / sum(mtx[s])
    mtxout = np.zeros( (Nsamps,Nsamps), dtype=np.float32 )
    for isamp in range(Nsamps):
        for jsamp in range(Nsamps):
            if np.isclose(mtxfrac[ lsamps[isamp] ],0).all() and np.isclose(mtxfrac[ lsamps[jsamp] ],0).all():
                mtxout[isamp,jsamp]=0
            else:
                mtxout[ isamp, jsamp ], _ = stats.pearsonr( mtxfrac[ lsamps[isamp] ], mtxfrac[ lsamps[jsamp] ] )
    tblout = pd.DataFrame(mtxout, columns=lsamps)
    tblout['sample']=lsamps
    tblout=tblout[ ['sample']+lsamps ]
    return tblout

def pctOverlap( jh, lsamps, suffix='_sum_hits' ):
    Nsamps=len(lsamps)
    mtx = jh[ [ '%s%s'%(s,suffix) for s in lsamps ] ]
    mtx.columns = [ s for s in lsamps ]
    mtxout = np.zeros( (Nsamps,Nsamps), dtype=np.float32 )
    for isamp in range(Nsamps):
        for jsamp in range(Nsamps):
            mtxout[ isamp, jsamp ] =  (mtx[ lsamps[isamp] ][  mtx[ lsamps[jsamp] ]>0 ] > 0).mean()
    tblout = pd.DataFrame(mtxout, columns=lsamps)
    tblout['sample']=lsamps
    tblout=tblout[ ['sample']+lsamps ]
    return tblout

def pctOverlapSym( jh, lsamps, suffix='_sum_hits' ):
    Nsamps=len(lsamps)
    mtx = jh[ [ '%s%s'%(s,suffix) for s in lsamps ] ]
    mtx.columns = [ s for s in lsamps ]
    mtxout = np.zeros( (Nsamps,Nsamps), dtype=np.float32 )
    for isamp in range(Nsamps):
        for jsamp in range(Nsamps):
            mtxout[ isamp, jsamp ] = ((mtx[ lsamps[isamp] ]>0) & (mtx[ lsamps[jsamp] ]>0)).sum() / float( ((mtx[ lsamps[isamp] ]>0) | (mtx[ lsamps[jsamp] ]>0)).sum() )

    tblout = pd.DataFrame(mtxout, columns=lsamps)
    tblout['sample']=lsamps
    tblout=tblout[ ['sample']+lsamps ]
    return tblout

# i,j --> sum fraction in sample j over top elements making X% of sample i.
# def overlapBeneathCumulSum( jh, lsamps, frac=0.90, suffix='_sum_hits' ):
#     Nsamps=len(lsamps)
#     mtx = jh[ [ '%s%s'%(s,suffix) for s in lsamps ] ]
#     mtx.columns = [ s for s in lsamps ]
#     mtxout = np.zeros( (Nsamps,Nsamps), dtype=np.float32 )
#     for isamp in range(Nsamps):
#         ttlSampI = mtx[ lsamps[isamp] ].sum()
#         ipermSortByI = np.argsort( -(mtx[lsamps[isamp]]) )

#         # print ttlSampI
#         # print mtx[lsamps[isamp]][ipermSortByI].cumsum()

#         ieltmin = np.argmax(mtx[lsamps[isamp]][ipermSortByI].cumsum()/float(ttlSampI) >= frac )

#         for jsamp in range(Nsamps):
#             ttlSampJ = mtx[ lsamps[jsamp] ].sum()
#             ttlBelowThreshJ = mtx[ lsamps[jsamp] ][ ipermSortByI ][ :ieltmin ]

#             mtxout[ isamp, jsamp ] = ttlBelowThreshJ.sum() / float(ttlSampJ)

#     tblout = pd.DataFrame(mtxout, columns=lsamps)
#     tblout['sample']=lsamps
#     tblout=tblout[ ['sample']+lsamps ]
#     return tblout

def overlapBeneathCumulSum( jh, lsamps, frac=0.90, suffix='_sum_hits' ):
    Nsamps=len(lsamps)
    mtx = np.array( jh[ [ '%s%s'%(s,suffix) for s in lsamps ] ] )
    # mtx.columns = [ s for s in lsamps ]
    mtxout = np.zeros( (Nsamps,Nsamps), dtype=np.float32 )
    for isamp in range(Nsamps):
        ttlSampI = mtx[ :, isamp ].sum()
        ipermSortByI = np.argsort( -(mtx[:, isamp]) )

        # print ttlSampI
        # print mtx[lsamps[isamp]][ipermSortByI].cumsum()

        ieltmin = np.argmax(mtx[:,isamp][ipermSortByI].cumsum()/float(ttlSampI) >= frac )

        for jsamp in range(Nsamps):
            ttlSampJ = mtx[ :, jsamp ].sum()
            ttlBelowThreshJ = mtx[  ipermSortByI, jsamp ][ :ieltmin ]

            mtxout[ isamp, jsamp ] = ttlBelowThreshJ.sum() / float(ttlSampJ)

    tblout = pd.DataFrame(mtxout, columns=lsamps)
    tblout['sample']=lsamps
    tblout=tblout[ ['sample']+lsamps ]
    return tblout

if __name__=='__main__':
    
    opts = OptionParser()    

    opts.add_option('','--inKey',dest='inKey')
    opts.add_option('','--colLibname',dest='colLibname')
    opts.add_option('','--inJointHisto',dest='inJointHisto')
    opts.add_option('','--outBase',dest='outBase')

    opts.add_option('','--suffix',default='',dest='suffix')

    (o, args) = opts.parse_args()
    
    sk = pd.read_table( o.inKey )

    jh = pd.read_table( o.inJointHisto )

    lLibs = list( sk[o.colLibname] ) 

    mtxOverlapBeneathCumulSum = overlapBeneathCumulSum( jh, lLibs, 0.90, o.suffix )
    mtxOverlapBeneathCumulSum.to_csv( '%s.cumulSumOfColUnderRowCumulsum90.txt'%o.outBase, index=False, sep='\t' )

    mtxPearsonr = pearsonr( jh, lLibs, o.suffix )
    mtxPearsonr.to_csv( '%s.pearsonr.txt'%o.outBase, index=False, sep='\t' )

    mtxSpearmanr = spearmanr( jh, lLibs, o.suffix )
    mtxSpearmanr.to_csv( '%s.spearmanr.txt'%o.outBase, index=False, sep='\t' )

    mtxDotProd = dotProd( jh, lLibs, o.suffix )
    mtxDotProd.to_csv( '%s.freqDotProd.txt'%o.outBase, index=False, sep='\t' )

    mtxPctOverlap = pctOverlap( jh, lLibs,  o.suffix)
    mtxPctOverlap.to_csv( '%s.pctovl.txt'%o.outBase, index=False, sep='\t' )

    mtxPctOverlapSym = pctOverlapSym( jh, lLibs, o.suffix)
    mtxPctOverlapSym.to_csv( '%s.pctovlsym.txt'%o.outBase, index=False, sep='\t' )

    
    
    
    




