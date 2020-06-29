# -*- coding: utf-8 -*-
"""
Created on Fri Nov 08 21:04:01 2019

@author: han-luo
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

cell = ['CCS' , 'NT5' , 'NT6' , 'fESC']
chrom=['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','X']
R = 200000
res = '200K'

sig_type = np.dtype({'names':['chr','score'],
                      'formats':['S4',np.float]})
sig_type_1 = np.dtype({'names':['chr','start' , 'end' , 'score'],
                      'formats':['S4',np.int , np.int , np.float]})

PCFolder = 'H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\'
Outfolder = 'H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\Plot\\'


#----------------------------------------------Plot---------------------------------------------------------------------------
size = (12, 12)
Left = 0.19 ; HB = 0.17 ; width = 0.7 ; HH = 0.7

for c in cell:
    HiCFolder = 'H:\\Workspace_New\\data\\HiC\\Matrix\\Matrix_200K\\' + c + '_200K\\'
    HiCFil = 'Correct_Merged_Reps_Local_Chromosome_Matrix.npz'
    HiCSource = os.path.join(HiCFolder , HiCFil)
    HiCData = np.load(HiCSource , allow_pickle = True)
    HiCData = HiCData['Matrix'][()]
    PCFil = c + '_compartment_200K.txt'
    PCSource = os.path.join(PCFolder , PCFil)
    PCData = np.loadtxt(PCSource , dtype = sig_type)

    OutFil = c + '_' + res + '_Heatmap_compartment.pdf'
    pp = PdfPages(os.path.join(Outfolder , OutFil))
    g = chrom[0]
    for g in chrom:
        matrix = HiCData[g]
        matrix[np.isnan(matrix)] = 0
        sig = PCData[PCData['chr'] == g]['score']
#        atac = atacData[atacData['chr'] == g]['score']
        nonzero = matrix[np.nonzero(matrix)]
        vmax = np.percentile(nonzero, 95)
        vmin = nonzero.min()
        ## Heatmap Plotting
        fig = plt.figure(figsize = size)
        ax = fig.add_axes([Left  , HB , width , HH])
        sc = ax.imshow(matrix, cmap = my_cmap, aspect = 'auto', interpolation = 'none',
                       extent = (0, len(matrix), 0, len(matrix)), vmax = vmax, vmin = vmin , origin = 'lower')
        cxlim = ax.get_xlim()
        cylim = ax.get_ylim()
        ## Ticks and Labels
        ticks = list(np.linspace(0 , len(matrix) , 5).astype(float))
        pos = [t * R for t in ticks]
        labels = [properU(p) for p in pos]
        ax.set_xticks(ticks)
        ax.set_xticklabels(labels , fontsize=10)
        ax.set_yticks(ticks)
        ax.set_yticklabels(labels , fontsize=10 , rotation = 'horizontal')
        ax.set_xlabel('chr' + g , fontsize=20,labelpad=30)

                        
#b3b3b3
        ax.set_xlim(cxlim)
        ax.set_ylim(cylim)                    
        ## Colorbar
        ax = fig.add_axes([Left + 0.6 , HB - 0.12 , 0.1 , 0.035])
        cbar = fig.colorbar(sc,cax = ax, orientation='horizontal')
        cbar.set_ticks([vmin , vmax])

        
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