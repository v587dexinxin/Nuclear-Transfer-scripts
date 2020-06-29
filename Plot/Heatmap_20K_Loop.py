# -*- coding: utf-8 -*-
"""
Created on Wed Feb 13 15:24:17 2019

@author: xxli
"""

from __future__ import division
import numpy as np
import pyBigWig
#from tadlib.calfea.analyze import getmatrix
import matplotlib
# Use a non-interactive backend
matplotlib.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
import os
import sys
from scipy import sparse
from tadlib.calfea import analyze
from palettable.colorbrewer.qualitative import Dark2_8
hexcolors = Dark2_8.hex_colors
#--------------------------------------------------------------------------
## Matplotlib Settings
matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
# Our Own Color Map
my_cmap = LinearSegmentedColormap.from_list('interaction',
                                            ['#FFFFFF','#CD0000'])
my_cmap.set_bad('#2672a1')
                
#-----------------------------------------------------Functions---------------------------------------------------------
def binbias(up, down):
    bias = 0
    zeromask = (up != 0) & (down != 0)
    if  zeromask.sum() <= 1:
        return bias
        
    upmean = up.mean(); downmean = down.mean()
    SD_1 = np.sum((up - upmean) ** 2) / (up.size * (up.size - 1))
    SD_2 = np.sum((down - downmean) ** 2) / (down.size * (down.size - 1))
    SD_pool = np.sqrt(SD_1 + SD_2)
    if SD_pool != 0:
        bias = (upmean - downmean) / SD_pool

    return bias
    
    
def CalDI(matrix):
    DIs = []
    n = 0
    shape = matrix.shape
    for j in matrix[:shape[0]]:
        if n <20:
            DIs.append(0)
        elif n >= 20 and n <= int(shape[0])-21:
            if len(j[j!=0])/int(shape[0]) < 0.05:
                bias = 0
            else:
                up = j[n-20:n] 
                down = j[n+1:n+21]
                bias = binbias(up,down) 
            DIs.append(bias)
        else:
            DIs.append(0)
        n += 1
    return DIs


def caxis_H(ax):
    """
    Axis Control for HeatMaps.
    """
    ax.yaxis.set_ticks_position('right')
    ax.xaxis.set_ticks_position('bottom')
    ax.tick_params(axis = 'both', labelsize = 12, length = 5, pad = 7)

def caxis_S(ax, color):
    """
    Axis Control for signal plots.
    """
    for spine in ['left', 'top']:
        ax.spines[spine].set_visible(False)
    ax.tick_params(axis = 'both', bottom = True, top = False, left = False,
                   right = False, labelbottom = False, labeltop = False,
                   labelleft = False, labelright = False)
    ax.spines['bottom'].set_lw(1.5)
    ax.spines['bottom'].set_color(color)
    ax.spines['bottom'].set_alpha(0.9)
    ax.spines['bottom'].set_linestyle('dotted')
    
    ax.spines['right'].set_lw(0.5)
    ax.spines['right'].set_color('#b3b3b3')
    ax.spines['right'].set_alpha(0.9)

def properU(pos):
    """
    Express a genomic position in a proper unit (KB, MB, or both).
    
    """
    i_part = int(pos) // 1000000 # Integer Part
    d_part = (int(pos) % 1000000) // 1000 # Decimal Part
    
    if (i_part > 0) and (d_part > 0):
        return ''.join([str(i_part), 'M', str(d_part), 'K'])
    elif (i_part == 0):
        return ''.join([str(d_part), 'K'])
    else:
        return ''.join([str(i_part), 'M'])
    
def getmatrix(inter,l_bin,r_bin):
    inter_matrix = np.zeros((r_bin - l_bin, r_bin - l_bin),dtype = float )
    mask = (inter['bin1'] >= l_bin) & (inter['bin1'] < r_bin) & \
           (inter['bin2'] >= l_bin) & (inter['bin2'] < r_bin)
    inter_extract = inter[mask]
    for i in inter_extract:
        if i['bin1'] != i['bin2']:
            inter_matrix[i['bin1'] - l_bin][i['bin2'] - l_bin] += i['IF']
            inter_matrix[i['bin2'] - l_bin][i['bin1'] - l_bin] += i['IF']
        else:
            inter_matrix[i['bin1'] - l_bin][i['bin2'] - l_bin] += i['IF']
    return inter_matrix

#--------------------------------------------------Files----------------------------------------------------------------

OutFolder = '/public/home/xxli/data/BDF1_New/HiC/Plot/heatmap_loop_20K_compare_0.1_ww4'
LoopFolder = '/public/home/xxli/data/BDF1_New/HiC/Loop/Raw_20K_0.1_ww4'

chrom = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','X']
cells = {'CCS':0 , 'NT5':1 , 'NT6':2 , 'fESC':3}
cell = ['CCS' , 'NT5' , 'NT6' , 'fESC']
R = 20000
res = '20K'
interval = 400
size = (12, 12)
Left = 0.2 ; HB = 0.2 ; width = 0.6 ; HH = 0.6

tad_type = np.dtype({'names':['chr' , 'start' , 'end'],
                      'formats':['S4' , np.int , np.int]})
 
for c in cell:
    OutFil = c + '_20K_6.pdf'
    pp = PdfPages(os.path.join(OutFolder , OutFil))
    HiCFolder = '/public/home/xxli/data/BDF1_New/HiC/HapHiC_workspace/' + c + '_workspace/' + c + '_20K'
    HiCFil = 'Correct_Merged_Reps_Local_Chromosome_Matrix.npz'
    HiCSource = os.path.join(HiCFolder , HiCFil)
    HiCData = np.load(HiCSource)
    LoopFil = 'Cluster_' + c + '_loops20K_6.txt'
    LoopFil_raw = c + '_loops_20K.txt'
    LoopSource = os.path.join(LoopFolder , LoopFil)
    LoopSource_raw = os.path.join(LoopFolder , LoopFil_raw)
    LoopData = np.loadtxt(LoopSource , skiprows = 1 , dtype = tad_type , usecols = (0 , 1 , 2))
    LoopData_raw = np.loadtxt(LoopSource_raw , skiprows = 1 , dtype = tad_type , usecols = (0 , 1 , 2))
    for g in chrom:
        lib = HiCData['Matrix'][()][g]
        loops = LoopData[LoopData['chr'] == g]
        loops_raw = LoopData_raw[LoopData_raw['chr'] == g]
        for i in range(len(lib) // interval +1):
            startHiC = i*interval
            endHiC = (i+1)*interval
            matrix = lib[startHiC:endHiC , startHiC:endHiC]
            nonzero = matrix[np.nonzero(matrix)]
            if nonzero.shape[0] == 0:
                break
            vmax = np.percentile(nonzero, 95)
            fig = plt.figure(figsize = size)
            ax = fig.add_axes([Left  , HB , width , HH])
            sc = ax.imshow(matrix, cmap = my_cmap, aspect = 'auto', interpolation = 'none',
                           extent = (0, len(matrix), 0, len(matrix)), vmax = vmax, origin = 'lower')
            cxlim = ax.get_xlim()
            cylim = ax.get_ylim()
            ## Ticks and Labels
            ticks = list(np.linspace(0 , interval , 5).astype(float))
            pos = [((startHiC + t) * R) for t in ticks]
            labels = [properU(p) for p in pos]
            ax.set_xticks(ticks)
            ax.set_xticklabels(labels)
            ax.set_xlim(cxlim)
            ax.set_ylim(cylim)
            ax.set_yticks(ticks)
            ax.set_yticklabels(labels , rotation = 'vertical')
            ax.set_xlabel('chr' + g , fontsize=15,labelpad=10)
            
            ## loop data
            mask = (loops['start'] >= startHiC * R) & (loops['end'] < endHiC * R)
            extract = loops[mask]
            mask_raw = (loops_raw['start'] >= startHiC * R) & (loops_raw['end'] <= endHiC * R)
            extract_raw = loops_raw[mask_raw] 
            for p in extract:
                x = p['start'] // R - startHiC
                y = p['end'] // R - startHiC
                ax.scatter(y ,x , color = '', edgecolors = 'b', s = 8)
            for p in extract_raw:
                x = p['end'] // R - startHiC
                y = p['start'] // R - startHiC
                ax.scatter(y + 0.5,x + 0.5 , marker = 's' , color = 'b', edgecolors = '', s = 2)
                
            ax.set_xlim(cxlim)
            ax.set_ylim(cylim)                    
            caxis_H(ax)
            ## Colorbar
            ax = fig.add_axes([Left + width + 0.05, HB, 0.05, HH])
            fig.colorbar(sc, cax = ax)
            pp.savefig(fig)
            plt.close(fig)
    pp.close()