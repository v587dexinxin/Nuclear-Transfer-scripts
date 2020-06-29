# -*- coding: utf-8 -*-
"""
Created on Wed Nov 07 10:08:08 2018

@author: xxli
"""

from __future__ import division
from scipy import sparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LinearSegmentedColormap
import cPickle
import sys
import os


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
    DIs = {}
    for i in chrom:
        n = 0
        DIs[i] =[]
        shape = matrix[i].shape
        for j in matrix[i][:shape[0]]:
            if n <20:
                DIs[i].append(0)
            elif n >= 20 and n <= int(shape[0])-21:
                if len(j[j!=0])/int(shape[0]) < 0.05:
                    bias = 0
                else:
                    up = j[n-20:n] 
                    down = j[n+1:n+21]
                    bias = binbias(up,down) 
                DIs[i].append(bias)
            else:
                DIs[i].append(0)
            n += 1
    return DIs


def caxis_H(ax):
    """
    Axis Control for HeatMaps.
    """
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    ax.tick_params(axis = 'both', labelsize = 12, length = 5, pad = 7)

def caxis_S(ax, color):
    """
    Axis Control for signal plots.
    """
    for spine in ['right', 'top']:
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
    

#--------------------------------------------------Files----------------------------------------------------------------


TadFolder = 'D:\\Workspace_New\\data\\HiC\\HiTAD\\respective\\stable_600K\\bottom_domain'
OutFolder = 'D:\\Workspace_New\\Plot\\Figure1\\HiC-Heatmap\\heatmap_DI_bottom_domain_sliding_window_respective_reps_40K'

chrom = ['1' , '2' , '3' , '4' , '5' , '6' , '7' , '8' , '9' , '10' , '11' , '12' , '13' , '14' , '15' , '16' , '17' , '18' , '19' ,'X']
cells = {'CCS':0 , 'NT5':1 , 'NT6':2 , 'fESC':3}
cell = ['CCS','NT5','NT6']
R = 40000
res = '40K'
interval = 250
size = (12, 12)
Left = 0.15 ; HB = 0.2 ; width = 0.6 ; HH = 0.6

tad_type = np.dtype({'names':['chr' , 'start' , 'end'],
                     'formats':['S4' , np.int , np.int]})
di_type = np.dtype({'names':['chr' , 'DI'],
                     'formats':['S4' , np.float]})
 
for c in cell:
    HiCFolder = 'D:\\Workspace_New\\data\\HiC\\Matrix\\Matrix_40K\\' + c + '_40K'
    HiCFil = 'Correct_Merged_Reps_Local_Chromosome_Matrix.npz'
    HiCSource = os.path.join(HiCFolder , HiCFil)
    HiCData = np.load(HiCSource)
    TadFil = c + '_40K_allreps.txt'
#    DIFil = c + '_40K_DI.txt'
    TadSource = os.path.join(TadFolder , TadFil)
#    DISource = os.path.join(TadFolder , DIFil)
    TadData = np.loadtxt(TadSource , usecols = (0 , 1 , 2) ,dtype = tad_type)
    DIData = CalDI(HiCData['Matrix'][()])
    for g in chrom:
        OutFil = c + '_chr' + g + '_40K_Heatmap_DI_domain.pdf'
        pp = PdfPages(os.path.join(OutFolder , OutFil))
        lib = HiCData['Matrix'][()][g]
        di = DIData[g]
        tads = TadData[TadData['chr'] == g]
        for i in range(len(lib) // interval +1):
            startHiC = i*interval
            endHiC = (i+1)*interval
            matrix = lib[startHiC:endHiC , startHiC:endHiC]
            nonzero = matrix[np.nonzero(matrix)]
            if nonzero.shape[0] == 0:
                break
            vmax = np.percentile(nonzero, 95)
            sig = np.array(di[startHiC : endHiC])
            if len(matrix) < interval:
                matrix0 = lib[startHiC:endHiC , startHiC:endHiC]
                matrix = np.zeros((interval , interval))
                nonzero = matrix0[np.nonzero(matrix0)]
                matrix[np.nonzero(matrix0)] = nonzero
                sig0 = np.array(di[startHiC : endHiC])
                sig = np.zeros(interval)
                nonzero = sig0[np.nonzero(sig0)]
                sig[np.nonzero(sig0)] = nonzero
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
            
            caxis_H(ax)
            
            ## Domain Boundaries
            mask = (tads['end'] > startHiC * R ) & (tads['start'] < endHiC * R )
            extract = tads[mask]
            pairs = [(bi['start']//R - startHiC, bi['end']//R - startHiC)
                     for bi in extract]
            for corner in pairs:
                if (corner[0] <= 0):
                    ax.plot([0, corner[1]], [corner[1], corner[1]], color = '#b3b3b3',
                            linewidth = 1)
                    ax.plot([corner[1], corner[1]], [0, corner[1]], color = '#b3b3b3',
                            linewidth = 1)
                elif (corner[1] >= interval):
                    ax.plot([corner[0], corner[0]], [corner[0], interval],
                            color = '#b3b3b3', linewidth = 1)
                    ax.plot([corner[0], interval], [corner[0], corner[0]],
                            color = '#b3b3b3', linewidth = 1)
                else:
                    ax.plot([corner[0], corner[0]], [corner[0], corner[1]],
                            color = '#b3b3b3', linewidth = 1)
                    ax.plot([corner[0], corner[1]], [corner[0], corner[0]],
                            color = '#b3b3b3', linewidth = 1)
                    ax.plot([corner[0], corner[1]], [corner[1], corner[1]],
                            color = '#b3b3b3', linewidth = 1)
                    ax.plot([corner[1], corner[1]], [corner[0], corner[1]],
                            color = '#b3b3b3', linewidth = 1)
            ax.set_xlim(cxlim)
            ax.set_ylim(cylim)                    
            caxis_H(ax)
            ## Colorbar
            ax = fig.add_axes([Left + 0.45 , HB - 0.15 , 0.15 , 0.035])
            fig.colorbar(sc,cax = ax, orientation='horizontal')

            ##Signal Tracks
            ax = fig.add_axes([Left + width, HB, 0.1, HH])
            ax.fill_betweenx(np.arange(sig.size) , sig , facecolor = 'blue' , edgecolor = 'none' )
            ax.set_xlabel('DI',fontsize=15)
            ax.set_ylim((cxlim[0] , cxlim[1] - 1))
            caxis_S(ax, 'blue')
            pp.savefig(fig)
            plt.close(fig)
        pp.close()
            
            
        
    
    
    

