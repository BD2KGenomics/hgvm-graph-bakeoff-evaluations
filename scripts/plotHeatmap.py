#!/usr/bin/env python
"""
Create a heatmap in PDF format from a tsv file.  Expects labels in 0th columns and rows. Rows or columns containing None values are omitted.
"""

# note to self: copied over from https://github.com/glennhickey/teHmm/blob/master/parameterAnalysis.py


#Copyright (C) 2014 by Glenn Hickey
#
#Released under the MIT license, see LICENSE.txt
#!/usr/bin/env python

import os
import sys
import logging
import argparse
import numpy as np
import scipy
import scipy.spatial
import scipy.cluster
import matplotlib
import math
matplotlib.use('Agg')
import matplotlib.backends.backend_pdf as pltBack
import matplotlib.lines as lines
import matplotlib.patches as patches
import matplotlib.pylab  as pylab
import matplotlib.pyplot as plt
import matplotlib.markers
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, LogLocator
from matplotlib.ticker import LogFormatter
import matplotlib.mlab as mlab
from matplotlib.mlab import PCA
import matplotlib.cm as cm
from matplotlib.colors import LogNorm, NoNorm, Normalize
import scipy.cluster.hierarchy as sch

from computeVariantsDistances import read_tsv

def parse_args(args):
    parser = argparse.ArgumentParser(description=__doc__, 
        formatter_class=argparse.RawDescriptionHelpFormatter)
        
    # General options
    parser.add_argument("dist_mat", type=str,
                        help="distance matrix tsv file")
    parser.add_argument("out_file", type=str,
                        help="path to output pdf for heatmap")
    parser.add_argument("--log_scale", action="store_true", default=False,
                        help="use log scaling (linear by default)")
    parser.add_argument("--vmax", type=float, default=None,
                        help="maximum value of colorbar (automatic by default)")
    parser.add_argument("--skip", type=str, default="",
                        help="comma-separated list of keywords"
                        " such that input row or column will ignored if its header contains one")
                            
    args = args[1:]

    return parser.parse_args(args)


def hierarchicalCluster(points, normalizeDistances=False, linkageMethod='average'):
    """ Perform simple hiearchical clustering using euclidian distance"""
    assert points is not None and len(points) > 0
    distanceMatrix = scipy.spatial.distance.pdist(points, "minkowski")
    if normalizeDistances is True:
        distanceMatrix /= math.pow(len(points[0]), 0.5)
    hc = scipy.cluster.hierarchy.linkage(distanceMatrix, method=linkageMethod)
    return hc

def rankHierarchies(hcList, rankStat = "branch_length"):
    """ produce a ranking of hierarchical clusterings using one of a variety
    of statistics. Each element of hcList is the output of
    hierarchichalCluster()
    Return list of integer indexes in increasing order"""

    inputRanking = [x for x in xrange(len(hcList))]
        
    if rankStat == "branch_length":
        totalLengths = [np.sum([x[2] for x in hc]) for hc in hcList]
        sortedRanking = sorted(zip(totalLengths, inputRanking), reverse=True)
        ranks = list(zip(*sortedRanking)[1])
        return ranks
    else:
        raise RuntimeError("Rankstat %s not recognized" % rankStat)

def plotHierarchicalClusters(hcList, titles, leafNames, outFile):
    """ print out a bunch of dendrograms to a PDF file.  Each element of
    hcList is the output of hierarchichalCluster()"""
    cols = 4
    rows = int(np.ceil(float(len(hcList)) / float(cols)))
    width=15
    height= 10

    pdf = pltBack.PdfPages(outFile)
    fig = plt.figure(figsize=(width, height))
    plt.clf()
    for i, hc in enumerate(hcList):
        # +1 below is to prevent 1st element from being put last
        # (ie sublot seems to behave as if indices are 1-base)
        ax = plt.subplot(rows, cols, (i + 1) % len(hcList))
        dgram = scipy.cluster.hierarchy.dendrogram(
            hc, color_threshold=100000, labels=leafNames, show_leaf_counts=False)
#            p=6,
#            truncate_mode='lastp')
        plt.title(titles[i])
        plt.setp(plt.xticks()[1], rotation=-90, fontsize=10)
#        ax.set_ylim((0.,0.5))
    fig.tight_layout()
    fig.savefig(pdf, format = 'pdf')
    pdf.close()

def pcaFlatten(points, outDim = 2):
    """ flatten points to given dimensionality using PCA """
    assert outDim == 2
    
    # will get LinAlgError: SVD did not converge exception if all points
    # lie on some plane (ie all values equal for some dimension so we
    # have to check for that first
    dims = []
    for dim in xrange(len(points[0])):
        vals = set()
        for point in points:
            vals.add(point[dim])
        if len(vals) > 1:
            dims.append(dim)
    assert len(dims) > 0
    cleanPoints = np.array([[point[i] for i in dims] for point in points])
    assert len(cleanPoints) > 0
    
    pca = PCA(cleanPoints)
    return pca.Y, np.sum(pca.fracs[:2])

colorList = ['#1f77b4', # dark blue
             '#aec7e8', # light blue
            '#ff7f0e', # bright orange
            '#ffbb78', # light orange
            '#4B4C5E', # dark slate gray
            '#9edae5', # light blue 
            '#7F80AB', # purple-ish slate blue
            '#c7c7c7', # light gray
            '#9467bd', # dark purple
            '#c5b0d5', # light purple
            '#d62728', # dark red
            '#ff9896', # light red
                 ]
def plotPoints2d(distList, titles, stateNames, outFile, xRange=None,
                 yRange=None, ptSize=100, xLabel=None, yLabel=None, cols=2,
                 width=10, rowHeight=5, singleLegendPerRow=False,
                 markerList=["o", "s", "^"], rgbs = None):
    """ plot some points to a pdf file.  Some other marker sets are
                 matplotlib.markers.MarkerStyle.filled_markers
                 matplotlib.markers.MarkerStyle.markers
    """

    rows = int(np.ceil(float(len(distList)) / float(cols)))
    height=rowHeight * rows
    alpha = 0.7

    # pallettes are here : cm.datad.keys()
    if rgbs is None:
        rgbs = [cm.gist_rainbow_r(float(i) / float(len(stateNames)))
                for i in xrange(len(stateNames))]
        for i in xrange(len(rgbs)):
            rgbs[i] = list(rgbs[i])
            rgbs[i][3] = alpha
    
    pdf = pltBack.PdfPages(outFile)
    fig = plt.figure(figsize=(width, height))
    plt.clf()
    for i,  dist in enumerate(distList):
        # +1 below is to prevent 1st element from being put last
        # (ie sublot seems to behave as if indices are 1-base)
        ax = plt.subplot(rows, cols, (i + 1) % len(distList))
        plotList = []
        for j in xrange(len(dist)):
            plotList.append(plt.scatter(dist[j][0], dist[j][1],
                                        c=rgbs[j % len(rgbs)],
                                        s=ptSize,
                                        marker=markerList[j%len(markerList)]))
        #plt.axis('equal')
        plt.grid(True)
        plt.title(titles[i])
        if i % cols == 0 or not singleLegendPerRow:
            # write legend
            plt.legend(plotList, stateNames, 
            scatterpoints=1,
            loc='upper left',
            ncol=3,
            fontsize=8)
        if yRange is not None:
            assert xrange is not None
            ax.set_ylim(yRange)
            ax.set_xlim(xRange)
        if xLabel is not None:
            plt.xlabel(xLabel)
        if yLabel is not None:
            plt.ylabel(yLabel)

    fig.tight_layout()
    fig.savefig(pdf, format = 'pdf')
    pdf.close()

def plotHeatMap(inputArray, rowNames, colNames, outFile, leftTree = False, topTree = False,
                xLabelPosition=None, yLabelPosition=None, aspect='auto', logNorm = False,
                vmax=None):
    """ from here
    http://stackoverflow.com/questions/2455761/reordering-matrix-elements-to-reflect-column-and-row-clustering-in-naiive-python
    """

    # make sure array is a numpy
    array = np.array(inputArray, dtype=np.float)
    
    width=10.6
    height= 9
    sX = -0.17
    sY = 0.15
    pdf = pltBack.PdfPages(outFile)
    fig = plt.figure(figsize=(width, height))
    #print array
    #print rowNames
    #print colNames

    # Compute and plot dendrogram.
    if leftTree is True:
        #axdendro = fig.add_axes([0.09,0.1,0.2,0.6], frame_on=False)
        axdendro = fig.add_axes([0.26+sX,0.1+sY,0.05,0.6], frame_on=False)
        #Y = sch.linkage(array, method='centroid')
        Y = hierarchicalCluster(array, linkageMethod = 'single')
        Z = sch.dendrogram(Y, orientation='right', color_threshold=100000)
        axdendro.set_xticks([])
        axdendro.set_yticks([])
        
    # Compute and plot second dendrogram.
    if topTree is True:
        #ax2 = fig.add_axes([0.3,0.71,0.6,0.2], frame_on=False)
        ax2 = fig.add_axes([0.3+sX,0.69+sY,0.6,0.05], frame_on=False)
        Y = hierarchicalCluster(array.T, linkageMethod = 'single')        
        Z2 = sch.dendrogram(Y, color_threshold=100000)
        ax2.set_xticks([])
        ax2.set_yticks([])

    # plot matrix as heatmap
    #axmatrix = fig.add_axes([0.3,0.1,0.6,0.8])
    axmatrix = fig.add_axes([0.3+sX,0.1+sY,0.6,0.6])
    if leftTree is True:
        idx1 = Z['leaves']
        array = array[idx1,:]
        rowNames = list(np.array(rowNames)[idx1])
    if topTree is True:
        idx2 = Z2['leaves']
        array = array[:,idx2]
        colNames = list(np.array(colNames)[idx2])
    # picture of built-in colormaps http://wiki.scipy.org/Cookbook/Matplotlib/Show_colormaps
    if logNorm is True:
        norm = LogNorm(vmax=vmax)
    else:
        norm = Normalize(vmin=0.0, vmax=vmax)
    im = axmatrix.matshow(array,
                          #interpolation='nearest',
                          aspect=aspect,
                          #origin='lower',
                          cmap=plt.cm.YlGnBu, norm=norm)

    #axmatrix.set_xticks([])
    #axmatrix.set_yticks([])
    
    axmatrix.set_xticks([i for i in xrange(-1, len(colNames))])
    axmatrix.set_yticks([i for i in xrange(-1, len(rowNames))])

    if leftTree is True or yLabelPosition == 'right':
        axmatrix.yaxis.set_label_position('right')
        axmatrix.yaxis.tick_right()
    if topTree is True or xLabelPosition == "bottom":
        axmatrix.xaxis.set_label_position('bottom')
        pylab.xticks(rotation=-90)
        axmatrix.xaxis.tick_bottom()


    axmatrix.set_xticklabels(['']+colNames)
    axmatrix.set_yticklabels(['']+rowNames)

    plt.setp(plt.xticks()[1], rotation=90)

    # Plot colorbar.
    #axcolor = fig.add_axes([0.91,0.1,0.02,0.6])
    axcolor = fig.add_axes([1.09+sX,0.1+sY,0.02,0.6])
    pylab.colorbar(im, cax=axcolor)
    
    fig.savefig(pdf, format = 'pdf')
    pdf.close()
    plt.close()

def remove_skips(mat, col_names, row_names, skips):
    """ remove rows and columns if their name contains something in list
    """
    dead_rows = []
    dead_cols = []
    # scan rows and columns, remembering indexes to delete
    for skip in skips:
        if len(skip) > 0:
            for i in range(len(row_names)):
                if skip in row_names[i]:
                    dead_rows.append(i)
            for j in range(len(col_names)):
                if skip in col_names[j]:
                    dead_cols.append(j)

    # go ahead and delete the indexes, remembering to decrement by offset
    # in dead list to account for shrinking..
    for x in range(len(dead_rows)):
        del mat[dead_rows[x] - x]
        del row_names[dead_rows[x] - x]

    for x in range(len(dead_cols)):
        for i in range(len(row_names)):
            del mat[i][dead_cols[x] - x]
        del col_names[dead_rows[x] - x]
        
    return mat, col_names, row_names

def main(args):

    options = parse_args(args)
    assert os.path.splitext(options.out_file)[1] == ".pdf"
    options.skip = options.skip.split(",")

    mat, col_names, row_names, row_label = read_tsv(options.dist_mat)
    mat, col_names, row_names = remove_skips(mat, col_names, row_names, options.skip)

    plotHeatMap(mat, row_names, col_names,
                options.out_file,
                leftTree=True,
                topTree=True,
                logNorm=options.log_scale,
                vmax=options.vmax)

if __name__ == "__main__":
    sys.exit(main(sys.argv))
