# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 13:29:51 2020

@author: han-luo
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
                                            ['teal' , '#FFFFFF' ,'firebrick'])
my_cmap.set_bad('#D3D3D3')
#2672a1                
#'#FFFFFF','#CD0000'                
cells = {'CCS':0 , 'NT5':1 , 'NT6':2 , 'fESC':3}
cell = ['CCS' , 'NT5' , 'NT6' , 'fESC']
chrom=['1' , '2' , '3' , '4' , '5' , '6' , '7' , '8' , '9' , '10' , '11' , '12' , '13' , '14' , '15' , '16' , '17' , '18' , '19' , 'X']
chrom_1={'1':0 , '2':1 , '3':2 , '4':3 , '5':4 , '6':5 , '7':6 , '8':7 , '9':8 , '10':9 , '11':10 , '12':11 , '13':12 , '14':13 , '15':14 , '16':15 , '17':16 , '18':17 , '19':18 , 'X':19}
OutFolder = 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\S9_figs'
size = (12, 12)   
Left = 0.2 ; HB = 0.2 ; width = 0.6 ; HH = 0.6 
R = 200000

def properU(pos):
    """
    Express a genomic position in a proper unit (KB, MB, or both).
    
    """
    i_part = int(pos) // 1000000 # Integer Part
    d_part = (int(pos) % 1000000) // 1000 # Decimal Part
    
    if (i_part > 0) and (d_part > 0):
        return ''.join([str(i_part), 'M', str(d_part), 'K'])
    elif (i_part == 0):
        return ''.join([str(d_part), 'M'])
    else:
        return ''.join([str(i_part), 'M'])
    
    


NT5 = np.load('H:\\Workspace_New\\data\\HiC\\Matrix\\Matrix_200K\\NT5_200K\\Correct_Merged_Reps_Local_Chromosome_Matrix.npz' , allow_pickle=True)
NT5 = NT5['Matrix'][()]['7']
NT6 = np.load('H:\\Workspace_New\\data\\HiC\\Matrix\\Matrix_200K\\NT6_200K\\Correct_Merged_Reps_Local_Chromosome_Matrix.npz' , allow_pickle=True)
NT6 = NT6['Matrix'][()]['7'] 
fESC = np.load('H:\\Workspace_New\\data\\HiC\\Matrix\\Matrix_200K\\fESC_200K\\Correct_Merged_Reps_Local_Chromosome_Matrix.npz' , allow_pickle=True)
fESC = fESC['Matrix'][()]['7'] 


for i in range(len(NT5)):
    NT5[i , i] = 0
    NT6[i , i] = 0
    fESC[i , i] = 0

NT = (NT5 / NT5.sum() + NT6 / NT6.sum()) / 2 
fESC = fESC / fESC.sum()
matrix = fESC - NT


NT = (NT5 + NT6) / 2 
matrix = fESC - NT



    
OutFil = 'Chr7_fESC-NT5NT6_heatmap_1.pdf'
pp = PdfPages(os.path.join(OutFolder , OutFil))
ticksx = [0 , 133 , 134 , 728]
ticksy = [0 , 133 , 134 , 728]
posx = [(t * R) for t in ticksx]
posy = [(t * R) for t in ticksy]
labelsx = [properU(posx[0]) , '' , '' , properU(posx[3])]
labelsy = [properU(posy[0]) , '' , '' , properU(posy[3])]
nonzero = matrix[np.nonzero(matrix)]
vmax = np.percentile(nonzero, 97)
vmin = np.percentile(nonzero, 3)
fig = plt.figure(figsize = size)
ax = fig.add_axes([Left  , HB , width , HH])
sc = ax.imshow(matrix, cmap = my_cmap, aspect = 'auto', interpolation = 'none', vmax = vmax, vmin = vmin ,
               extent = (0, matrix.shape[1], 0, matrix.shape[0]), origin = 'lower')
ax.set_xticks(ticksx)
ax.set_yticks(ticksy)
ax.set_xticklabels(labelsx , fontsize=15)
ax.set_yticklabels(labelsy , fontsize=15)
ax.set_xlabel('Chr9' , fontsize=20)
ax.set_ylabel('Chr2' , fontsize=20)
ax.set_xlim(0, matrix.shape[1])
ax.set_ylim(0, matrix.shape[0])
ax.set_title('fESC_NT_difference heatmap' , fontsize=25)
## Colorbar
ax = fig.add_axes([Left , HB - 0.1 , 0.6 , 0.035])
cbar = fig.colorbar(sc,cax = ax, orientation='horizontal')
cbar.set_ticks([vmin , 0 , vmax])
pp.savefig(fig)
pp.close()