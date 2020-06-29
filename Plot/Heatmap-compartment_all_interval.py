# -*- coding: utf-8 -*-
"""
Created on Tue Nov 06 16:21:10 2018

@author: xxli
"""

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
# Use a non-interactive backend
matplotlib.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
import os
from matplotlib.colors import LinearSegmentedColormap
import cPickle

# Our Own Color Map
my_cmap = LinearSegmentedColormap.from_list('interaction',
                                            ['#FFFFFF','#CD0000'])
my_cmap.set_bad('#2672a1')

#---------------------------------------------------Functions-------------------------------------------------------------
def caxis_H(ax):
    """
    Axis Control for HeatMaps.
    """
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    ax.tick_params(axis = 'both', bottom = False, top = False, left = False,
                   right = False, labelbottom = False, labeltop = False,
                   labelleft = False, labelright = False)
def caxis_colorbar(ax):
    """
    Axis Control for HeatMaps.
    """
    ax.tick_params(axis = 'both', bottom = True, top = False, left = False,
                   right = False, labelbottom = True, labeltop = False,
                   labelleft = False, labelright = False , labelsize = 25)
    
    
def caxis_S_vertical(ax, color):
    """
    Axis Control for signal plots.
    """
    for spine in ['right', 'top']:
        ax.spines[spine].set_visible(False)
    ax.tick_params(axis = 'y', bottom = True, top = False, left = False,
                   right = False, labelbottom = True, labeltop = False,
                   labelleft = False, labelright = False)
    ax.spines['bottom'].set_lw(1.5)
    ax.spines['bottom'].set_color(color)
    ax.spines['bottom'].set_alpha(0.9)
    ax.spines['bottom'].set_linestyle('dotted')
def caxis_S_horizontal(ax, color):
    """
    Axis Control for PCA plots.
    """
    for spine in ['right', 'top']:
        ax.spines[spine].set_visible(False)
    ax.tick_params(axis = 'both', bottom = False, top = False, left = True,
                   right = False, labelbottom = False, labeltop = False,
                   labelleft = True, labelright = False , labelsize = 23)
    ax.spines['left'].set_lw(1.5)
    ax.spines['left'].set_color(color)
    ax.spines['left'].set_alpha(0)
    ax.spines['left'].set_linestyle('dotted')

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
    
def UpdateDI(DI):
    """
    """
    New_DI = []
    New_index = []

    for index in range(len(DI) - 1):
        if DI[index] * DI[index + 1] < 0:
            New_DI.append(DI[index])
            New_DI.append(0)
            New_index.append(index)
            New_index.append(index + 0.5)
        else:
            New_DI.append(DI[index])
            New_index.append(index)
    
    return np.array(New_index), np.array(New_DI)

    
#-----------------------------------------------Files-------------------------------------------------------------------------    
cells = {'CCS':0 , 'NT5':1 , 'NT6':2 , 'fESC':3}
cell = ['CCS','NT5','NT6','fESC']
chrom=['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','X']
R = 200000
res = '200K'

sig_type = np.dtype({'names':['chr','score'],
                      'formats':['S4',np.float]})
sig_type_1 = np.dtype({'names':['chr','start' , 'end' , 'score'],
                      'formats':['S4',np.int , np.int , np.float]})

PCFolder = 'D:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new'
Outfolder = 'D:\\Workspace_New\\High_Quality_Figures\\Fig2\\heatmap_200K'

interval = ('2' , 105500000 , 110000000)
    
#----------------------------------------------Plot---------------------------------------------------------------------------
size = (12, 12)
Left = 0.19 ; HB = 0.17 ; width = 0.7 ; HH = 0.7

for c in cell:
    HiCFolder = 'D:\\Workspace_New\\data\\HiC\\Matrix\\Matrix_200K\\' + c + '_200K'
    HiCFil = 'Correct_Merged_Reps_Local_Chromosome_Matrix.npz'
    HiCSource = os.path.join(HiCFolder , HiCFil)
    HiCData = np.load(HiCSource)
    PCFil = c + '_compartment_200K.txt'
    PCSource = os.path.join(PCFolder , PCFil)
    PCData = np.loadtxt(PCSource , dtype = sig_type)

    OutFil = c + '_' + res + '_Heatmap_compartment_3.pdf'
    pp = PdfPages(os.path.join(Outfolder , OutFil))
    g = chrom[0]
    for g in chrom:
        matrix = HiCData['Matrix'][()][g]
        sig = PCData[PCData['chr'] == g]['score']
#        atac = atacData[atacData['chr'] == g]['score']
        nonzero = matrix[np.nonzero(matrix)]
        vmax = np.percentile(nonzero, 95)
        ## Heatmap Plotting
        fig = plt.figure(figsize = size)
        ax = fig.add_axes([Left  , HB , width , HH])
        sc = ax.imshow(matrix, cmap = my_cmap, aspect = 'auto', interpolation = 'none',
                       extent = (0, len(matrix), 0, len(matrix)), vmax = vmax, origin = 'lower')
        cxlim = ax.get_xlim()
        cylim = ax.get_ylim()
        ## Ticks and Labels
#        ticks = [0 , int(len(matrix) * 0.75)]
#        pos = [t * R for t in ticks]
#        labels = [properU(p) for p in pos]
#        ax.set_xticks(ticks)
#        ax.set_xticklabels(labels , fontsize=5)
#        ax.set_yticks(ticks)
#        ax.set_yticklabels(labels , fontsize=5 , rotation = 'horizontal')
#        ax.set_xlabel('chr' + g , fontsize=20,labelpad=30)

        ## Selected interval
        if g == interval[0]:
            ax.plot([interval[1]//R, interval[1]//R], [interval[1]//R, interval[2]//R] , color = '#7093DB' , linestyle = '--' , linewidth = 5)
            ax.plot([interval[1]//R, interval[2]//R], [interval[1]//R, interval[1]//R] , color = '#7093DB' , linestyle = '--' , linewidth = 5)
            ax.plot([interval[1]//R, interval[2]//R], [interval[2]//R, interval[2]//R] , color = '#7093DB' , linestyle = '--' , linewidth = 5)
            ax.plot([interval[2]//R, interval[2]//R], [interval[1]//R, interval[2]//R] , color = '#7093DB' , linestyle = '--' , linewidth = 5)
        else:
            pass
                        
#b3b3b3
        ax.set_xlim(cxlim)
        ax.set_ylim(cylim)                    
        caxis_H(ax)
        ## Colorbar
        ax = fig.add_axes([Left + 0.6 , HB - 0.12 , 0.1 , 0.035])
        cbar = fig.colorbar(sc,cax = ax, orientation='horizontal')
        cbar.set_ticks([0 , int(vmax)])
        cbar.set_ticklabels([0 , int(vmax)] )
        caxis_colorbar(ax)
        
        ##PCA Tracks
        PCA_index , PCA = UpdateDI(sig)
        ax = fig.add_axes([Left, HB + width , width , 0.1])
        ax.fill_between(PCA_index , PCA , where = PCA >= 0 , facecolor = '#E47833' , edgecolor = 'none' )
        ax.fill_between(PCA_index , PCA , where = PCA <= 0 , facecolor = '#7093DB' , edgecolor = 'none' )
        ax.set_xlim(0 , PCA_index.max())
        ytick = [round(sig.min() * 0.6667 , 2) , 0.00 , round(sig.max() * 0.6667, 2)]
        ax.set_yticks(ytick)
        ax.set_ylim((sig.min() * 1.1 , sig.max()* 1.1))
#        ax.set_ylabel('PC1',fontsize=40,rotation = 'horizontal' , labelpad = 45)
        
        
        caxis_S_horizontal(ax, 'black')
        pp.savefig(fig)
        plt.close(fig)
    pp.close()
    
