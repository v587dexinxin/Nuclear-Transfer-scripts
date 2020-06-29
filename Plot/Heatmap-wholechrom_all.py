# -*- coding: utf-8 -*-
"""
Created on Mon Nov 19 16:15:49 2018

@author: xxli
"""

from __future__ import division
import numpy as np
#from tadlib.calfea.analyze import getmatrix
import matplotlib
# Use a non-interactive backend
matplotlib.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
import os
import sys
#from tadlib.calfea import analyze

#--------------------------------------------------------------------------
## Matplotlib Settings
matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
# Our Own Color Map
my_cmap = LinearSegmentedColormap.from_list('interaction',
                                            ['#FFFFFF' ,'#CD0000'])
my_cmap.set_bad('#D3D3D3')
#2672a1                
#'#FFFFFF','#CD0000'                
cells = {'CCS':0 , 'NT5':1 , 'NT6':2 , 'fESC':3}
cell = ['CCS' , 'NT5' , 'NT6' , 'fESC']
chrom=['1' , '2' , '3' , '4' , '5' , '6' , '7' , '8' , '9' , '10' , '11' , '12' , '13' , '14' , '15' , '16' , '17' , '18' , '19' , 'X']
chrom_1={'1':0 , '2':1 , '3':2 , '4':3 , '5':4 , '6':5 , '7':6 , '8':7 , '9':8 , '10':9 , '11':10 , '12':11 , '13':12 , '14':13 , '15':14 , '16':15 , '17':16 , '18':17 , '19':18 , 'X':19}
HiCFolder = 'D:\\Workspace_New\\data\\HiC\\Matrix\\Matrix_1M'
OutFolder = 'D:\\Workspace_New\\Plot\\HiC-Heatmap\\Selected\\Selected_2\\heatmap_1M'
size = (12, 12)   
Left = 0.25 ; HB = 0.2 ; width = 0.6 ; HH = 0.6 
R = 1000000

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
    
    
    
for c in cell:
    OutFil = c + '_whole_1M_all.pdf'
    pp = PdfPages(os.path.join(OutFolder , OutFil))
    wholeFil = c + '_Correct_Merged_Reps_Whole_Genome_Matrix.npz'
    wholeSource = os.path.join(HiCFolder , wholeFil)
    wholeData = np.load(wholeSource)
    matrix = wholeData['Matrix'][()]
#    m = np.zeros((20 , 20))
#    for i in chrom:
#        i_start = wholeData['Bins'][()][i][0]
#        i_end = wholeData['Bins'][()][i][1] + 1
#        for j in chrom:
#            j_start = wholeData['Bins'][()][j][0]
#            j_end = wholeData['Bins'][()][j][1] + 1
#            inter = matrix[i_start:i_end , j_start:j_end]
#            s = inter.sum() / inter.size
#            m[chrom_1[i]][chrom_1[j]] = s
    Bool = np.zeros_like(matrix, dtype = bool)
    ticks = []
    for g in chrom:
        start = wholeData['Bins'][()][g][0]
        end = wholeData['Bins'][()][g][1] + 1
        Bool[start:end , start:end] = 1
        matrix = np.ma.array(matrix, mask = Bool)
        ticks.append((start + end) / 2)
    pos = [(t * R) for t in ticks]
    labels = [properU(p) for p in pos]
    nonzero = matrix[np.nonzero(matrix)]
    vmax = np.percentile(nonzero, 97)
    fig = plt.figure(figsize = size)
    ax = fig.add_axes([Left  , HB , width , HH])
    sc = ax.imshow(matrix, cmap = my_cmap, aspect = 'auto', interpolation = 'none', vmax = vmax,
                   extent = (0, matrix.shape[1], 0, matrix.shape[0]), origin = 'lower')
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    ax.set_xticklabels(chrom , fontsize=15)
    ax.set_yticklabels(chrom , fontsize=15)

    ax.set_xlim(0, matrix.shape[1])
    ax.set_ylim(0, matrix.shape[0])
    ax.set_title(c , fontsize=25)
    ## Colorbar
    ax = fig.add_axes([Left + 0.45 , HB - 0.1 , 0.15 , 0.035])
    fig.colorbar(sc,cax = ax, orientation='horizontal')
    pp.savefig(fig)
    pp.close()