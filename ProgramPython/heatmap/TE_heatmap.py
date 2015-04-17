#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import re
import os
import argparse
import matplotlib as mp
mp.use('pdf')
import matplotlib.pyplot as plt
import pandas as pd
from urllib2 import urlopen
import numpy as np
import pylab as pl
import matplotlib.colors as colors
import random
from matplotlib import gridspec

def usage():
    test="name"
    message='''
python TE_heatmap.py --input sly1_strain.list

    '''
    print message

def set_ticks_X(ax, ylen, ylab):
    fig = plt.gcf()
    #fig.set_size_inches(8, 8)

    # turn off the frame
    ax.set_frame_on(False)

    # put the major ticks at the middle of each cell
    #ax.set_yticks(np.arange(qtl.shape[0]) + 0.5, minor=False)
    ax.set_xticks(np.arange(ylen) + 0.5, minor=False)

    # want a more natural, table-like display
    ax.invert_yaxis()
    ax.xaxis.tick_top()

    # Set the labels
    ax.set_xticklabels(ylab, minor=False)
    ax.set_yticklabels([])

    # rotate the
    plt.xticks(rotation=90)

    ax.grid(False)

    # Turn off all the ticks
    ax = plt.gca()

    for t in ax.xaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False
    for t in ax.yaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False

    return ax

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

def set_ticks_XY(ax, xlen, ylen, xlab, ylab):
    fig = plt.gcf()
    #fig.set_size_inches(8, 11)

    # turn off the frame
    ax.set_frame_on(False)

    # put the major ticks at the middle of each cell
    ax.set_yticks(np.arange(xlen) + 0.5, minor=False)
    ax.set_xticks(np.arange(ylen) + 0.5, minor=False)

    # want a more natural, table-like display
    ax.invert_yaxis()
    ax.xaxis.tick_top()

    # Set the labels
    ax.set_xticklabels(xlab, minor=False)
    ax.set_yticklabels(ylab, minor=False)

    # rotate the
    plt.xticks(rotation=60)

    ax.grid(False)

    # Turn off all the ticks
    ax = plt.gca()

    for t in ax.xaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False
    for t in ax.yaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False

    return ax

def colorbar(heatmap,  cbaxes, ticks, ticke):
    # Add a colorbar below the heatmap
    # add_axes refer to [left, bottom, width, height], where the coordinates are just fractions that go from 0 to 1 of the plotting area.
    #cbaxes = fig.add_axes([0.15, 0.15, 0.2, 0.02]) 
    cb = plt.colorbar(heatmap, cax = cbaxes, orientation='horizontal', ticks=np.linspace(ticks, ticke, 3))
    cb.ax.set_xticklabels(['NaN', 'Absence', 'Presence']) 
    #cb = pl.colorbar(heatmap, orientation='horizontal', shrink=0.5825,
    #                 fraction=0.02, pad=-0.035,ticks=np.linspace(-1, 1, 5),
    #                 use_gridspec=True)
    #cb.set_label("$\mathrm{Pearson's}\ r$")
    #return fig

def heatmap(gene_list, project):
    # Create a figure.
    figsize=(18,8)
    fig = pl.figure(figsize=figsize)
    pl.subplots_adjust(bottom=0.2)
    gs = gridspec.GridSpec(1, 1, width_ratios=[1])

    # Read Dataframe from file
    qtl = pd.read_table(gene_list, index_col=0)

    #use truncated color table from 0.2 to 0.8
    cmap = plt.cm.get_cmap("Blues")
    new_cmap = truncate_colormap(cmap, 0.2, 0.8) 

    # Draw heatmap for variations
    #ax = fig.add_subplot(131, frame_on=False)
    ax0 = plt.subplot(gs[0])
    heatmap = ax0.pcolor(qtl, cmap=cmap, vmin=-1, vmax=1, alpha=0.8)
    # Set names for x,y_ticks
    ax0 = set_ticks_XY(ax0, qtl.shape[0], qtl.shape[1], qtl.columns, qtl.index)
    
    ticks0 = -1
    ticke0 = 1 
    cbaxes0 = fig.add_axes([0.15, 0.15, 0.2, 0.02])
    colorbar(heatmap, cbaxes0, ticks0, ticke0)
    
    # Save file
    fig.savefig('%s.pdf' %(project), bbox_inches='tight')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)

    if not args.output:
        args.output = 'sly1_Nc.heatmap'

    heatmap(args.input, args.output)

if __name__ == '__main__':
    main()

