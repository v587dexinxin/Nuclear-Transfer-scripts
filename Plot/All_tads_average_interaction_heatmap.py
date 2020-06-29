# -*- coding: utf-8 -*-
"""
Created on Sat Mar 02 10:21:34 2019

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
                

def caxis_H(ax):
    """
    Axis Control for HeatMaps.
    """
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    ax.tick_params(axis = 'both', bottom = True, top = False, left = True,
                   right = False, labelbottom = True, labeltop = False,
                   labelleft = True, labelright = False , length = 5 ,labelsize = 40  )

def caxis_colorbar(ax):
    """
    Axis Control for HeatMaps.
    """
    ax.tick_params(axis = 'both', bottom = True, top = False, left = False,
                   right = False, labelbottom = True, labeltop = False,
                   labelleft = False, labelright = False , labelsize = 25)
    
                
size = (12, 12)
Left = 0.2 ; HB = 0.2 ; width = 0.6 ; HH = 0.6

tad_type = np.dtype({'names':['chr' , 'start' , 'end'],
                     'formats':['S4' , np.int , np.int]})

TadFolder = 'D:\\Workspace_New\\data\\HiC\\HiTAD\\respective\\stable_600K\\bottom_domain'
cells = ['CCS','NT5','NT6','fESC']


CCS_Lib = np.load('D:/Workspace_New/data/HiC/Matrix/Matrix_40K/CCS_40K/Correct_Merged_Reps_Local_Chromosome_Matrix.npz')
CCS_Lib = CCS_Lib['Matrix'][()]
fESC_Lib = np.load('D:/Workspace_New/data/HiC/Matrix/Matrix_40K/fESC_40K/Correct_Merged_Reps_Local_Chromosome_Matrix.npz')
fESC_Lib = fESC_Lib['Matrix'][()]
NT5_Lib = np.load('D:/Workspace_New/data/HiC/Matrix/Matrix_40K/NT5_40K/Correct_Merged_Reps_Local_Chromosome_Matrix.npz')
NT5_Lib = NT5_Lib['Matrix'][()]
NT6_Lib = np.load('D:/Workspace_New/data/HiC/Matrix/Matrix_40K/NT6_40K/Correct_Merged_Reps_Local_Chromosome_Matrix.npz')
NT6_Lib = NT6_Lib['Matrix'][()]
HiC_Data = {'CCS':CCS_Lib,
            'fESC':fESC_Lib,
            'NT5':NT5_Lib,
            'NT6':NT6_Lib}

tads = np.loadtxt(os.path.join(TadFolder , 'CCS_standard_domain.txt') , dtype = tad_type)

for c in cells:
    pp = PdfPages('D:\\Workspace_New\\High_Quality_Figures\\Fig3\\' + c + '_tads_average_interaction.pdf')
    matrix = np.zeros((50, 50))
    for i in tads:
        g = i['chr']
        s = (i['start'] - 500000) // 40000 
        e = (i['end'] + 500000) // 40000
        lib = HiC_Data[c][g]
        matrix += lib[s:e , s:e]
    matrix = matrix / (len(tads))
    for i in range(50):
        for j in range(50):
            if i == j:
                matrix[i][j] = 0
            else:
                pass
    nonzero = matrix[np.nonzero(matrix)]
    vmax = np.percentile(nonzero, 80)
    fig = plt.figure(figsize = size)
    ax = fig.add_axes([Left  , HB , width , HH])
    sc = ax.imshow(matrix, cmap = my_cmap, aspect = 'auto', interpolation = 'none',
                       extent = (0, len(matrix), 0, len(matrix)), vmax = vmax, origin = 'lower')
    ticks = [0 , 12.5 , 27.5 , 50]
    labels = ['0M','12.5M','27.5M','50M']
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels)
    ax.set_yticks(ticks)
    ax.set_yticklabels(labels, rotation = 'horizontal')
    ax.set_xlabel(c , fontsize = 50 )
    caxis_H(ax)
    ##-----------------------colorbar--------------------------------------------------------------------------
    ax = fig.add_axes([Left + 0.5 , HB - 0.11 , 0.1 , 0.035])
    cbar = fig.colorbar(sc,cax = ax, orientation='horizontal')
    cbar.set_ticks([0 , int(vmax)])
    caxis_colorbar(ax)
    
    pp.savefig(fig)
    pp.close()
    
    
    