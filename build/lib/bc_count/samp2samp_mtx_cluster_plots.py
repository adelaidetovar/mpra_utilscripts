# Author: Jacob Kitzman, kitzmanj

import sys
import os
import os.path
from optparse import OptionParser
from collections import defaultdict

import os.path as op 

import pandas as pd
import seaborn as sns
import pylab as plt

import scipy.stats as stats


# we can cluster on a number of different measures
lMeas = ['pearsonr','spearmanr','cumulSumOfColUnderRowCumulsum90','freqDotProd','pctovl','pctovlsym']
mMeasDesc = {'pearsonr':'pearson correlation',
             'spearmanr':'rank correlation',
             'cumulSumOfColUnderRowCumulsum90': 'Col sum of barcodes within row top 90%ile (column comes from row)',
             'freqDotProd':'Frequency dot product',
             'pctovl':'Percent overlap',
             'pctovlsym':'Percent overlap symmetric'             
            }


def make_heatmap( 
    fnmtx,
    fnplot,
    subsetRows = None,
    subsetCols = None,
    cluster = False,
    title = None,
    figsize=(12,8) ):

    tbls2s = pd.read_table( fnmtx )
    tbls2s = tbls2s.set_index('sample')

    if subsetRows is not None:
        if subsetCols is not None:
            tbls2s = tbls2s.loc[ subsetRows, subsetCols ]
        else:
            tbls2s = tbls2s.loc[ subsetRows ]
    else:
        if subsetCols is not None:
            tbls2s = tbls2s[ subsetCols ]


    if not cluster:    
        plt.figure(figsize=figsize)
        hm = sns.heatmap( tbls2s )
        _=hm.set_title( title if title is not None else fnmtx )

        plt.tight_layout()
        _=plt.savefig(fnplot)
    else:
        clgrid = sns.clustermap( tbls2s, figsize=figsize  )
        for tl in clgrid.ax_heatmap.yaxis.get_ticklabels():
            tl.set_rotation(0.)

        _=clgrid.fig.suptitle( title if title is not None else fnmtx )

        clgrid.savefig(fnplot)



def heatmap_samp2samp_set(
    fnmtxbase,
    fnplotbase,
    cluster=False,
    extension='png',
    **kwargs ):

    for meas in lMeas:
        if op.exists('%s.%s.txt'%( fnmtxbase, meas )):
            make_heatmap( fnmtx='%s.%s.txt'%( fnmtxbase, meas ),
                        fnplot='%s.%s.%s'%( fnplotbase, meas, extension ),
                        cluster=cluster,
                        title=mMeasDesc[meas],
                        **kwargs )
        else:
            print('warning: could not find file %s.%s.txt'%( fnmtxbase, meas ))





def main():    
    opts = OptionParser()    

    opts.add_option('','--inSamp2sampMtxBase',dest='inSamp2sampMtxBase')
    opts.add_option('','--outPlotBase',dest='outPlotBase')
    opts.add_option('','--cluster',default=False,action='store_true',dest='cluster')
    opts.add_option('','--figSize',default='(12,8)',dest='figSize')

    opts.add_option('','--subsetToRowList',default=None,dest='subsetToRowList')
    opts.add_option('','--subsetToColList',default=None,dest='subsetToColList')

    (o, args) = opts.parse_args()

    subsetRows = None if o.subsetToRowList is None else [ x.strip() for x in open(o.subsetToRowList,'r').readlines() ]
    subsetCols = None if o.subsetToColList is None else [ x.strip() for x in open(o.subsetToColList,'r').readlines() ]
    
    heatmap_samp2samp_set( o.inSamp2sampMtxBase, 
                           o.outPlotBase,
                           cluster=o.cluster,
                           figsize=eval(o.figSize),
                           subsetRows=subsetRows,
                           subsetCols=subsetCols )



if __name__ == '__main__':
    main()
