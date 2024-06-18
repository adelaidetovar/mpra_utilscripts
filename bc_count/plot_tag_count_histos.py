# Author: Jacob Kitzman, kitzmanj

import sys
import argparse
from collections import defaultdict

import numpy as np

# import matplotlib as mpl
# import matplotlib.cm as cm
import matplotlib.pyplot as plt
# import matplotlib.gridspec as gridspec

# import scipy.misc as sm

import pandas as pd

def main():
    # count histograms are expected to be tab delimited:  tag count [other columns,e.g., constituent tags]
    # no header lines

    opts = argparse.ArgumentParser()
   
    opts.add_argument('--linCountHisto',dest='linCountHisto')
    opts.add_argument('--lNames',dest='lNames')
    
    opts.add_argument('--outPlot',dest='outPlot')
    
    opts.add_argument('--fxnApplyPreselection',default="lambda x:np.ones(x.shape[0],'bool')",dest='fxnApplyPreselection')

    opts.add_argument('--normalizeY',default=False,action='store_true',dest='normalizeY')    

    opts.add_argument('--x_linear',default=True,action='store_false',dest='x_log')
    opts.add_argument('--y_linear',default=True,action='store_false',dest='y_log')

    opts.add_argument('--annotMaxY',default=False,action='store_true',dest='annotMaxY')

    opts.add_argument('--ylabFormatterExpr',default="%.1f",dest='ylabFormatterExpr')

    o = opts.parse_args()

    fxnApplyPreselection = eval(o.fxnApplyPreselection)
        
    fig = plt.figure( figsize=(8,8) )
    ax = fig.add_subplot( 111 )

    vRange=(np.inf,-np.inf)

    if o.y_log:
        plt.yscale('log')
    if o.x_log:
        plt.xscale('log')

    ax2 = ax.twinx()

    rkMax=0

    for it in range(len(o.linCountHisto.split(','))):
        tblIn = pd.read_csv(o.linCountHisto.split(',')[it],sep='\t',header=None)
        vIn = tblIn[ 1 ].values

        vIn = vIn[ fxnApplyPreselection(vIn) ]

        vIn = np.sort( vIn )[::-1]
        rkIn = np.arange( 1, vIn.shape[0]+1 )        

        rkMax = max(rkMax, rkIn.shape[0] )

        if o.normalizeY:
            ax.step( rkIn, vIn/float(vIn.sum()), label=o.lNames.split(',')[it] )
        else:
            ax.step( rkIn, vIn, label=o.lNames.split(',')[it] )

        vRange=(min(vRange[0],min(vIn)), max(vRange[1],max(vIn)) )

        ax2.step( rkIn, 100.*vIn.cumsum()/float(vIn.sum()),
                    color='red',
                           linestyle='-' )

        if o.annotMaxY:
            ylv = o.ylabFormatterExpr%vIn[0]

            ax.annotate('max Y=%s'%ylv, 
                        xy=(1, vIn[0]),  xycoords='data',
                        xytext=(0.1, 1.0), textcoords='axes fraction',
                        arrowprops=dict(facecolor='black', shrink=0.05),
                        horizontalalignment='left', verticalalignment='top',
                        )


    ax.set_xlim( 0.9, 1.01*rkMax )
    ax.set_ylim(0.1*vRange[0],10.*vRange[1])
    ax2.set_ylim(0,100)
    
    if o.x_log:
        ax2.set_xscale('log')

    ax.set_xlabel('Rank among barcodes')
    ax2.set_ylabel('Cumulative percentage of reads')
    ax.set_ylabel('Number of hits from barcode')

    for tl in ax2.get_yticklabels():
        tl.set_color('r')

    ax.legend()

    plt.savefig( o.outPlot )



if __name__ = '__main__':
    main()
