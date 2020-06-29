# -*- coding: utf-8 -*-
"""
Created on Wed Nov 21 21:38:06 2018

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
    ax.tick_params(axis = 'both', labelsize = 12, length = 5, pad = 7)

def caxis_S(ax, color):
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
def caxis_PCA(ax, color):
    """
    Axis Control for PCA plots.
    """
    for spine in ['right', 'top']:
        ax.spines[spine].set_visible(False)
    ax.tick_params(axis = 'both', bottom = True, top = False, left = False,
                   right = False, labelbottom = True, labeltop = False,
                   labelleft = False, labelright = False)
    ax.spines['left'].set_lw(1.5)
    ax.spines['left'].set_color(color)
    ax.spines['left'].set_alpha(0.9)
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
    
#-----------------------------------------------Files-------------------------------------------------------------------------    
cells = {'CCS':0 , 'NT5':1 , 'NT6':2 , 'fESC':3}
cell = ['CCS' , 'NT5' , 'NT6' , 'fESC']
chrom=['1' , '2' , '3' , '4' , '5', '6' ,'7' ,'8' , '9' , '10' ,'11' , '12' ,'13' ,'14' ,'15' ,'16' ,'17' ,'18' ,'19' ,'X']
R = 200000
res = '200K'

sig_type = np.dtype({'names':['chr','score'],
                      'formats':['S4',np.float]})
sig_type_1 = np.dtype({'names':['chr','start' , 'end' , 'score'],
                      'formats':['S4',np.int , np.int , np.float]})

PCFolder = 'D:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new'
ATACFolder = 'D:\\Workspace_New\\data\ATAC\\signal\\normalization\\bedgraph_1K'
Outfolder = 'D:\\Workspace_New\\Plot\\Figure1\\HiC-Heatmap\\heapmap_200K'
#atacFolder = 'G:\\data\\nuclear_transfer\\sigFolder\\ATAC_coverage_1K'
    
#----------------------------------------------Plot---------------------------------------------------------------------------
size = (12, 12)
Left = 0.05 ; HB = 0.25 ; width = 0.5 ; HH = 0.5

for c in cell:
    HiCFolder = 'D:\\Workspace_New\\data\\HiC\\Matrix\\Matrix_200K\\' + c + '_200K'
    HiCFil = 'Correct_Merged_Reps_Local_Chromosome_Matrix.npz'
    HiCSource = os.path.join(HiCFolder , HiCFil)
    HiCData = np.load(HiCSource)
    PCFil = c + '_compartment_200K.txt'
    PCFil_1 = c + '_R1_compartment_200K.txt'
    PCFil_2 = c + '_R2_compartment_200K.txt'
    PCSource = os.path.join(PCFolder , PCFil)
    PCSource_1 = os.path.join(PCFolder , PCFil_1)
    PCSource_2 = os.path.join(PCFolder , PCFil_2)
    PCData = np.loadtxt(PCSource , dtype = sig_type)
    PCData_1 = np.loadtxt(PCSource_1 , dtype = sig_type)
    PCData_2 = np.loadtxt(PCSource_2 , dtype = sig_type)
    ATACFil = c + '_ATAC_1K.bedgraph'
    ATACSource = os.path.join(ATACFolder , ATACFil)
    ATACData = np.loadtxt(ATACSource , dtype=sig_type_1)
    OutFil = c + '_' + res + '_Heatmap_compartment_atac_new.pdf'
    pp = PdfPages(os.path.join(Outfolder , OutFil))
    g = chrom[0]
    for g in chrom:
        matrix = HiCData['Matrix'][()][g]
        sig = PCData[PCData['chr'] == g]['score']
        sig_1 = PCData_1[PCData_1['chr'] == g]['score']
        sig_2 = PCData_2[PCData_2['chr'] == g]['score']
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
        ticks = list(np.linspace(0, len(matrix), 5).astype(int))
        pos = [t * R for t in ticks]
        labels = [properU(p) for p in pos[:4]]        
        ax.set_xticks(ticks)
        ax.set_xticklabels(labels , fontsize=5)
        ax.set_yticks(ticks)
        ax.set_yticklabels(labels , fontsize=5 , rotation = 'horizontal')
        ax.set_xlabel('chr' + g , fontsize=20,labelpad=10)
        ax.set_xlim(cxlim)
        ax.set_ylim(cylim)
        caxis_H(ax)
        ## Colorbar
        ax = fig.add_axes([Left + 0.5 , HB - 0.1 , 0.1 , 0.035])
        cbar = fig.colorbar(sc,cax = ax, orientation='horizontal')
        cbar.set_ticks([0 , int(vmax/2) , int(vmax)])
       
        ##PCA Tracks
        pc = {'R1':sig_1 , 'R2':sig_2 , 'all':sig}
        step = 0
        for s in ['R1' , 'R2' ,'all']:
            ax = fig.add_axes([Left + width + step, HB , 0.1 , HH])
            ax.fill_betweenx(np.arange(pc[s].size) , pc[s] , facecolor = 'green' , edgecolor = 'none' )
            xtick = [round(sig.min() + 0.005, 3) , 0 , round(sig.max() - 0.02 , 3)]
            ax.set_xticks(xtick)
            ax.set_xlabel(s,fontsize=15)
            ax.set_xlim((sig.min() , sig.max()))
            ax.set_ylim((cxlim[0] , cxlim[1] - 1))
            caxis_PCA(ax, 'green')
            step += 0.1
        ##signal Tracks
        atac = ATACData[ATACData['chr'] == g]
        length = int(len(matrix) * R/1000)
        new = np.zeros(length)
        for i in atac:
            start = i['start'] // 1000
            end = i['end'] // 1000
            score = i['score']
            n = end - start
            for j in np.arange(n):
                new[start + j] = score
            
        color = 'blue'
        vm = new.max()
        ax = fig.add_axes([Left + width + step, HB, 0.1 , HH])
        ax.fill_betweenx(np.arange(len(new)) , new , facecolor = color , edgecolor = 'none' )
        ax.set_xlabel('ATAC',fontsize=15)
        ax.set_xlim((0 , vm))
        ax.set_ylim((0 , len(new)))
        caxis_S(ax, color)
        pp.savefig(fig)
        plt.close(fig)
    pp.close()